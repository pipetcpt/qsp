################################################################################
## Atrial Fibrillation QSP — Interactive Shiny Dashboard (v2.0)
## ============================================================
## Tabs:
##   1. 환자 프로파일  (Patient Profile) — CHA2DS2-VASc, HAS-BLED, demographics
##   2. 약동학 PK      (Pharmacokinetics) — Amiodarone, Apixaban, Metoprolol,
##                     Diltiazem, Dronedarone, Flecainide, Rivaroxaban, Warfarin
##   3. 부정맥 PD      (AF Dynamics)     — AF burden, ERP, HR, QTc, reentry
##   4. 항응고 & 뇌졸중 (Anticoagulation) — FXa, Thrombin, INR proxy, stroke risk
##   5. 시나리오 비교   (Scenario Comparison) — 6 strategies side-by-side
##   6. 바이오마커      (Biomarkers)      — NT-proBNP, CRP, D-dimer, Troponin
##
## Dependencies: shiny, shinydashboard, mrgsolve, ggplot2, dplyr, tidyr,
##               DT, plotly, shinycssloaders
################################################################################

suppressPackageStartupMessages({
  library(shiny)
  library(shinydashboard)
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(mrgsolve)
})

# ==============================================================================
# mrgsolve MODEL — AF QSP v2 (rate control / rhythm control / anticoagulation)
# ==============================================================================
af_model_code <- '
$PARAM
// ---- Amiodarone (Class III, multi-channel) ----
ka_AMIO=0.06  CL_AMIO=3.5   V1_AMIO=40    V2_AMIO=4200
k12_AMIO=0.02 k21_AMIO=0.002 F_AMIO=0.46
IC50_AMIO_ERP=0.5  Emax_AMIO_ERP=0.35
IC50_AMIO_HR=1.0   Emax_AMIO_HR=0.45
IC50_AMIO_QTc=0.8  Emax_AMIO_QTc=60.0

// ---- Dronedarone (Class III, amiodarone analogue) ----
ka_DRON=1.5   CL_DRON=130   V1_DRON=1400  F_DRON=0.70
IC50_DRON_ERP=0.08 Emax_DRON_ERP=0.22
IC50_DRON_HR=0.10  Emax_DRON_HR=0.30

// ---- Flecainide (Class IC) ----
ka_FLEC=1.0   CL_FLEC=8.5   V1_FLEC=600   F_FLEC=0.95
IC50_FLEC_ERP=0.3  Emax_FLEC_ERP=0.28
IC50_FLEC_HR=0.5   Emax_FLEC_HR=0.10

// ---- Apixaban (FXa inhibitor) ----
ka_APIX=1.2   CL_APIX=3.3   V1_APIX=21    F_APIX=0.50
IC50_APIX_FXa=0.08 Emax_APIX_FXa=0.95

// ---- Rivaroxaban (FXa inhibitor) ----
ka_RIVA=1.5   CL_RIVA=10.0  V1_RIVA=50    F_RIVA=0.80
IC50_RIVA_FXa=0.12 Emax_RIVA_FXa=0.90

// ---- Warfarin (Vitamin K antagonist — INR proxy via indirect PD) ----
ka_WARF=0.80  CL_WARF=0.20  V1_WARF=10.0  F_WARF=0.99
Kin_VKA=0.02  Kout_VKA=0.02 IC50_WARF_VKA=0.5 Imax_WARF_VKA=0.85

// ---- Metoprolol (beta-1 blocker) ----
ka_METRO=1.5  CL_METRO=65   V1_METRO=290  F_METRO=0.40
IC50_METRO_HR=25   Emax_METRO_HR=0.40

// ---- Diltiazem (CCB) ----
ka_DILT=1.3   CL_DILT=45    V1_DILT=220   F_DILT=0.40
IC50_DILT_HR=100   Emax_DILT_HR=0.38

// ---- Digoxin (cardiac glycoside) ----
ka_DIGO=0.60  CL_DIGO=8.0   V1_DIGO=500   F_DIGO=0.70
IC50_DIGO_HR=0.6   Emax_DIGO_HR=0.30

// ---- AF / EP disease model ----
AF0=0.60     ERP0=180   kfib=0.001  kfib_ERP=20  kAF_remod=0.005
HR0_AF=140   QTc0=400   kQTc_ERP=0.5
AngII0=1.0   ROS0=1.0   NE0=1.0     IL6_0=1.0    SMAD_base=0.5
kAngII_fib=0.003  kROS_fib=0.002  kNE_decay=0.1

// ---- Coagulation model ----
FXa0=1.0     Thrombin0=1.0  kStroke_base=0.035  kStroke_Thr=0.04
kThr_FXa=0.8  INR0=1.0

// ---- Biomarkers ----
kBNP_base=200.0   kBNP_AF=500.0
kDDimer_base=0.4  kDDimer_AF=1.8   kDDimer_Thr=0.8
kTropI_base=0.02  kTropI_rate=0.005

$CMT
// PK compartments — rate control
GI_METRO C1_METRO
GI_DILT  C1_DILT
GI_DIGO  C1_DIGO
// PK compartments — rhythm control
GI_AMIO  C1_AMIO  C2_AMIO
GI_DRON  C1_DRON
GI_FLEC  C1_FLEC
// PK compartments — anticoagulation
GI_APIX  C1_APIX
GI_RIVA  C1_RIVA
GI_WARF  C1_WARF  VKA_effect
// PD state variables
AF_BURDEN ERP Fibrosis QTc HR_AF
AngII ROS NE IL6 SMAD23
FXa Thrombin STROKE_RISK
NT_proBNP_st CRP_st Ddimer_st TropI_st

$ODE
// ---- PK ODEs ----
// Metoprolol
dxdt_GI_METRO = -ka_METRO * GI_METRO;
dxdt_C1_METRO =  F_METRO * ka_METRO * GI_METRO - (CL_METRO/V1_METRO) * C1_METRO;

// Diltiazem
dxdt_GI_DILT  = -ka_DILT * GI_DILT;
dxdt_C1_DILT  =  F_DILT * ka_DILT * GI_DILT  - (CL_DILT/V1_DILT) * C1_DILT;

// Digoxin
dxdt_GI_DIGO  = -ka_DIGO * GI_DIGO;
dxdt_C1_DIGO  =  F_DIGO * ka_DIGO * GI_DIGO  - (CL_DIGO/V1_DIGO) * C1_DIGO;

// Amiodarone (2-compartment)
dxdt_GI_AMIO  = -ka_AMIO * GI_AMIO;
dxdt_C1_AMIO  =  F_AMIO * ka_AMIO * GI_AMIO - (CL_AMIO/V1_AMIO) * C1_AMIO
                 - k12_AMIO * C1_AMIO + k21_AMIO * C2_AMIO;
dxdt_C2_AMIO  =  k12_AMIO * C1_AMIO - k21_AMIO * C2_AMIO;

// Dronedarone
dxdt_GI_DRON  = -ka_DRON * GI_DRON;
dxdt_C1_DRON  =  F_DRON * ka_DRON * GI_DRON  - (CL_DRON/V1_DRON) * C1_DRON;

// Flecainide
dxdt_GI_FLEC  = -ka_FLEC * GI_FLEC;
dxdt_C1_FLEC  =  F_FLEC * ka_FLEC * GI_FLEC  - (CL_FLEC/V1_FLEC) * C1_FLEC;

// Apixaban
dxdt_GI_APIX  = -ka_APIX * GI_APIX;
dxdt_C1_APIX  =  F_APIX * ka_APIX * GI_APIX  - (CL_APIX/V1_APIX) * C1_APIX;

// Rivaroxaban
dxdt_GI_RIVA  = -ka_RIVA * GI_RIVA;
dxdt_C1_RIVA  =  F_RIVA * ka_RIVA * GI_RIVA  - (CL_RIVA/V1_RIVA) * C1_RIVA;

// Warfarin
dxdt_GI_WARF  = -ka_WARF * GI_WARF;
dxdt_C1_WARF  =  F_WARF * ka_WARF * GI_WARF  - (CL_WARF/V1_WARF) * C1_WARF;
double Cp_WARF = C1_WARF / V1_WARF;
double Inh_WARF = Imax_WARF_VKA * Cp_WARF / (IC50_WARF_VKA + Cp_WARF);
dxdt_VKA_effect = Kin_VKA * (1.0 - Inh_WARF) - Kout_VKA * VKA_effect;

// ---- Derived concentrations ----
double Cp_AMIO  = C1_AMIO  / V1_AMIO;       // ug/mL
double Cp_DRON  = C1_DRON  / V1_DRON;       // ng/mL (already in ng range)
double Cp_FLEC  = C1_FLEC  / V1_FLEC;       // ug/mL
double Cp_APIX  = C1_APIX  / V1_APIX;       // ug/mL
double Cp_RIVA  = C1_RIVA  / V1_RIVA;       // ug/mL
double Cp_METRO = C1_METRO / V1_METRO;      // mg/L -> convert to ng/mL below
double Cp_DILT  = C1_DILT  / V1_DILT;       // mg/L
double Cp_DIGO  = C1_DIGO  / V1_DIGO;       // mg/L -> ng/mL proxy

double Cp_METRO_ng = Cp_METRO * 1000.0;
double Cp_DILT_ng  = Cp_DILT  * 1000.0;
double Cp_DIGO_ng  = Cp_DIGO  * 1000.0;

// ---- PD: Heart Rate effects (additive with ceiling) ----
double Eff_METRO_HR = Emax_METRO_HR * Cp_METRO_ng / (IC50_METRO_HR + Cp_METRO_ng);
double Eff_DILT_HR  = Emax_DILT_HR  * Cp_DILT_ng  / (IC50_DILT_HR  + Cp_DILT_ng);
double Eff_DIGO_HR  = Emax_DIGO_HR  * Cp_DIGO_ng  / (IC50_DIGO_HR  + Cp_DIGO_ng);
double Eff_AMIO_HR  = Emax_AMIO_HR  * Cp_AMIO     / (IC50_AMIO_HR  + Cp_AMIO);
double Eff_DRON_HR  = Emax_DRON_HR  * Cp_DRON     / (IC50_DRON_HR  + Cp_DRON);

double HR_red = Eff_METRO_HR + Eff_DILT_HR + Eff_DIGO_HR
                + Eff_AMIO_HR + Eff_DRON_HR;
if (HR_red > 0.82) HR_red = 0.82;

// ---- PD: ERP prolongation ----
double Eff_AMIO_ERP = Emax_AMIO_ERP * Cp_AMIO / (IC50_AMIO_ERP + Cp_AMIO);
double Eff_DRON_ERP = Emax_DRON_ERP * Cp_DRON / (IC50_DRON_ERP + Cp_DRON);
double Eff_FLEC_ERP = Emax_FLEC_ERP * Cp_FLEC / (IC50_FLEC_ERP + Cp_FLEC);
double Total_ERP_eff = Eff_AMIO_ERP + Eff_DRON_ERP + Eff_FLEC_ERP;
if (Total_ERP_eff > 0.60) Total_ERP_eff = 0.60;

// ---- PD: FXa inhibition ----
double Eff_APIX_FXa = Emax_APIX_FXa * Cp_APIX / (IC50_APIX_FXa + Cp_APIX);
double Eff_RIVA_FXa = Emax_RIVA_FXa * Cp_RIVA / (IC50_RIVA_FXa + Cp_RIVA);
double Eff_WARF_FXa = VKA_effect / (VKA_effect + 0.5) * 0.70;
double Total_FXa_inh = Eff_APIX_FXa + Eff_RIVA_FXa + Eff_WARF_FXa;
if (Total_FXa_inh > 0.95) Total_FXa_inh = 0.95;

// ---- QTc prolongation ----
double Eff_AMIO_QTc = Emax_AMIO_QTc * Cp_AMIO / (IC50_AMIO_QTc + Cp_AMIO);
double QTc_target = QTc0 + Eff_AMIO_QTc + 10.0 * Eff_DRON_ERP - kQTc_ERP * (ERP0 - ERP);
dxdt_QTc = 0.02 * (QTc_target - QTc);

// ---- ERP dynamics ----
double dERP_remodel = -kAF_remod * AF_BURDEN * ERP;
double dERP_drug    = Total_ERP_eff * ERP0;
double dERP_fib     = -kfib_ERP * Fibrosis / 24.0;
dxdt_ERP = (dERP_remodel + dERP_drug / 24.0 + dERP_fib)
           - 0.0001 * (ERP - (ERP0 - kfib_ERP * Fibrosis));

// ---- AF burden ----
double ERP_eff_AF  = 1.0 / (1.0 + exp((ERP - 200.0) / 20.0));
double Fib_eff_AF  = 1.5 * Fibrosis;
double kAF_in  = 0.003 * ERP_eff_AF * (1.0 + Fib_eff_AF);
double kAF_out = 0.002 * (1.0 - ERP_eff_AF);
dxdt_AF_BURDEN = kAF_in * (1.0 - AF_BURDEN) - kAF_out * AF_BURDEN
                 + 0.0001 * NE * AF_BURDEN * (1.0 - AF_BURDEN);
if (AF_BURDEN < 0.001 && dxdt_AF_BURDEN < 0) dxdt_AF_BURDEN = 0;
if (AF_BURDEN > 0.999 && dxdt_AF_BURDEN > 0) dxdt_AF_BURDEN = 0;

// ---- Fibrosis ----
double kFib_in = kfib * (AngII * kAngII_fib + ROS * kROS_fib + SMAD23 * 0.002);
dxdt_Fibrosis = kFib_in * (1.0 - Fibrosis) - 0.00005 * Fibrosis;
if (Fibrosis < 0.001 && dxdt_Fibrosis < 0) dxdt_Fibrosis = 0;
if (Fibrosis > 0.999 && dxdt_Fibrosis > 0) dxdt_Fibrosis = 0;

// ---- HR ----
double HR_target = HR0_AF * (1.0 - HR_red) * (1.0 + 0.3 * NE);
dxdt_HR_AF = 0.05 * (HR_target - HR_AF);

// ---- Neurohormonal ----
dxdt_AngII  = 0.02 * AF_BURDEN * (2.0 - AngII)  - 0.05 * (AngII  - AngII0);
dxdt_ROS    = 0.015 * AF_BURDEN * AngII           - 0.04 * (ROS    - ROS0);
dxdt_SMAD23 = 0.03  * AngII * (1.5 - SMAD23)     - 0.02 * (SMAD23 - SMAD_base);
dxdt_IL6    = 0.01  * AF_BURDEN * (3.0 - IL6)    - 0.03 * (IL6    - IL6_0);
dxdt_NE     = 0.01  * AF_BURDEN * (2.0 - NE)     - kNE_decay * (NE - NE0);

// ---- Coagulation ----
double FXa_prod = 0.1 * AF_BURDEN * (1.0 + 0.5 * Thrombin);
double FXa_inh  = Total_FXa_inh * FXa;
dxdt_FXa = FXa_prod - 0.15 * FXa - FXa_inh + 0.05 * (FXa0 - FXa);
dxdt_Thrombin = kThr_FXa * FXa * AF_BURDEN - 0.2 * Thrombin
                + 0.02 * (Thrombin0 - Thrombin);
double stroke_rate = kStroke_base * Thrombin * AF_BURDEN * 100.0
                     + kStroke_Thr * (Thrombin - 1.0);
if (stroke_rate < 0) stroke_rate = 0;
dxdt_STROKE_RISK = 0.01 * (stroke_rate - STROKE_RISK);

// ---- Biomarker states ----
double bnp_target   = kBNP_base + kBNP_AF * AF_BURDEN * (1.0 + 0.5 * Fibrosis);
dxdt_NT_proBNP_st   = 0.005 * (bnp_target - NT_proBNP_st);

double crp_target   = IL6 * 3.5;
dxdt_CRP_st         = 0.008 * (crp_target - CRP_st);

double ddimer_target = kDDimer_base + kDDimer_AF * AF_BURDEN + kDDimer_Thr * (Thrombin - 1.0);
dxdt_Ddimer_st       = 0.01 * (ddimer_target - Ddimer_st);

double tropI_target  = kTropI_base + kTropI_rate * HR_AF * AF_BURDEN;
dxdt_TropI_st        = 0.003 * (tropI_target - TropI_st);

$TABLE
// concentrations for output
double Cp_AMIO_out   = C1_AMIO / V1_AMIO;
double Cp_DRON_out   = C1_DRON / V1_DRON * 1000.0;  // ng/mL
double Cp_FLEC_out   = C1_FLEC / V1_FLEC * 1000.0;  // ng/mL
double Cp_APIX_ng    = C1_APIX / V1_APIX * 1000.0;  // ng/mL
double Cp_RIVA_ng    = C1_RIVA / V1_RIVA * 1000.0;  // ng/mL
double Cp_WARF_out   = C1_WARF / V1_WARF;
double INR_proxy     = 1.0 / (VKA_effect + 0.001) * 0.5 + 1.0; // simplified
double Cp_METRO_ng_out = C1_METRO / V1_METRO * 1000.0;
double Cp_DILT_ng_out  = C1_DILT  / V1_DILT  * 1000.0;
double Cp_DIGO_ng_out  = C1_DIGO  / V1_DIGO  * 1000.0;

double AntiXa_pct    = (Eff_APIX_FXa + Eff_RIVA_FXa) * 100.0;
double LA_diam       = 38.0 + 12.0 * Fibrosis;

$CAPTURE AF_BURDEN ERP QTc HR_AF STROKE_RISK Fibrosis
         Cp_AMIO_out Cp_DRON_out Cp_FLEC_out
         Cp_APIX_ng Cp_RIVA_ng Cp_WARF_out INR_proxy
         Cp_METRO_ng_out Cp_DILT_ng_out Cp_DIGO_ng_out
         FXa Thrombin AntiXa_pct LA_diam AngII ROS NE IL6
         NT_proBNP_st CRP_st Ddimer_st TropI_st
'

cat("Compiling AF QSP v2 mrgsolve model...\n")
af_mod <- tryCatch(
  mcode("AF_QSP_v2", af_model_code, quiet = TRUE),
  error = function(e) { message("mrgsolve compile error: ", e$message); NULL }
)

# ==============================================================================
# HELPER: CHA2DS2-VASc and HAS-BLED
# ==============================================================================
calc_chads2 <- function(age, htn, dm, hf, prior_stroke, pad, female) {
  score <- 0L
  if (prior_stroke) score <- score + 2L
  if (age >= 75)    score <- score + 2L
  else if (age >= 65) score <- score + 1L
  if (htn)    score <- score + 1L
  if (dm)     score <- score + 1L
  if (hf)     score <- score + 1L
  if (pad)    score <- score + 1L
  if (female) score <- score + 1L
  score
}

chads_stroke_pct <- function(score) {
  tbl <- c(0.2, 0.6, 1.5, 2.8, 4.0, 5.3, 6.6, 7.9, 9.6, 11.2, 12.5)
  tbl[min(max(score + 1L, 1L), length(tbl))]
}

calc_hasbled <- function(age, htn_unc, renal_abn, liver_abn,
                         prior_bleed, labile_inr, drugs_alc) {
  s <- 0L
  if (htn_unc)    s <- s + 1L
  if (renal_abn)  s <- s + 1L
  if (liver_abn)  s <- s + 1L
  if (age > 65)   s <- s + 1L
  if (prior_bleed)s <- s + 1L
  if (labile_inr) s <- s + 1L
  if (drugs_alc)  s <- s + 1L
  s
}

# ==============================================================================
# HELPER: run_sim — builds event table and calls mrgsim
# ==============================================================================
run_sim <- function(
    af_type      = "persistent",
    rate_drug    = "none",    dose_rate   = 0,
    rhythm_drug  = "none",    dose_rhythm = 0,
    anticoag     = "none",    dose_ac     = 0,
    t_days       = 365,
    dt_h         = 12,
    age          = 65,
    chads_score  = 3,
    htn          = TRUE,
    hf           = FALSE,
    dm           = FALSE,
    la_mm        = 44
) {
  if (is.null(af_mod)) return(NULL)

  af_init <- switch(af_type,
    paroxysmal   = 0.20, persistent = 0.60, longstanding = 0.85, 0.60)
  fib_init <- switch(af_type,
    paroxysmal   = 0.08, persistent = 0.22, longstanding = 0.42, 0.22)

  angII_init <- 1.0 + 0.35 * htn + 0.20 * dm + 0.45 * hf
  erp_init   <- max(150, 180 - (max(age - 50, 0)) * 0.35)
  stroke_init <- max(0.1, chads_score * 1.2 * af_init)
  la_fib_init <- max(0, (la_mm - 38) / 12)
  fib_init    <- max(fib_init, la_fib_init)

  iv <- c(
    GI_METRO=0, C1_METRO=0, GI_DILT=0, C1_DILT=0, GI_DIGO=0, C1_DIGO=0,
    GI_AMIO=0, C1_AMIO=0, C2_AMIO=0, GI_DRON=0, C1_DRON=0, GI_FLEC=0, C1_FLEC=0,
    GI_APIX=0, C1_APIX=0, GI_RIVA=0, C1_RIVA=0, GI_WARF=0, C1_WARF=0,
    VKA_effect=1.0,
    AF_BURDEN=af_init, ERP=erp_init, Fibrosis=fib_init,
    QTc=405, HR_AF=138, AngII=angII_init, ROS=1.2,
    FXa=1.0, Thrombin=1.0, STROKE_RISK=stroke_init,
    NE=1.1, IL6=1.3, SMAD23=0.55,
    NT_proBNP_st = 200 + 500 * af_init,
    CRP_st = 1.3 * 3.5,
    Ddimer_st = 0.4 + 1.8 * af_init,
    TropI_st  = 0.02 + 0.005 * 138 * af_init
  )

  ev_list <- list()

  # Rate control drug
  if (rate_drug != "none" && dose_rate > 0) {
    cmt_map <- c(metro=1, dilt=3, digo=5)
    cmt_n   <- switch(rate_drug, metro=1L, dilt=3L, digo=5L, 1L)
    freq_h  <- switch(rate_drug, metro=12, dilt=8, digo=24, 12)
    times   <- seq(0, t_days * 24 - 1, by = freq_h)
    ev_list[["rate"]] <- data.frame(time=times, cmt=cmt_n, amt=dose_rate, evid=1)
  }

  # Rhythm control drug
  if (rhythm_drug != "none" && dose_rhythm > 0) {
    if (rhythm_drug == "amio") {
      load <- seq(0, 28*24-1, by=12)
      ev_list[["amio_load"]] <- data.frame(time=load, cmt=7L, amt=400, evid=1)
      maint <- seq(28*24, t_days*24-1, by=24)
      ev_list[["amio_maint"]] <- data.frame(time=maint, cmt=7L, amt=dose_rhythm, evid=1)
    } else {
      cmt_n  <- switch(rhythm_drug, dron=10L, flec=12L, 10L)
      freq_h <- switch(rhythm_drug, dron=12, flec=12, 12)
      times  <- seq(0, t_days*24-1, by=freq_h)
      ev_list[["rhythm"]] <- data.frame(time=times, cmt=cmt_n, amt=dose_rhythm, evid=1)
    }
  }

  # Anticoagulation
  if (anticoag != "none" && dose_ac > 0) {
    if (anticoag == "apix") {
      times <- seq(0, t_days*24-1, by=12)
      ev_list[["ac"]] <- data.frame(time=times, cmt=14L, amt=dose_ac, evid=1)
    } else if (anticoag == "riva") {
      times <- seq(0, t_days*24-1, by=24)
      ev_list[["ac"]] <- data.frame(time=times, cmt=16L, amt=dose_ac, evid=1)
    } else if (anticoag == "warf") {
      times <- seq(0, t_days*24-1, by=24)
      ev_list[["ac"]] <- data.frame(time=times, cmt=18L, amt=dose_ac, evid=1)
    }
  }

  ev_df <- if (length(ev_list) > 0) {
    do.call(rbind, ev_list)
  } else {
    data.frame(time=0, cmt=1L, amt=0, evid=0L)
  }
  ev_df <- ev_df[order(ev_df$time), ]

  tryCatch({
    out <- mrgsim(
      af_mod,
      idata  = data.frame(ID=1),
      events = ev_df,
      init   = iv,
      tgrid  = tgrid(0, t_days * 24, dt_h),
      output = "df"
    )
    df <- as.data.frame(out)
    df$time_days <- df$time / 24
    df
  }, error = function(e) { message("Sim error: ", e$message); NULL })
}

# ==============================================================================
# CSS / Theme helpers
# ==============================================================================
score_badge <- function(score, low, high, low_col="#27AE60", mid_col="#F39C12", high_col="#E74C3C") {
  col <- if (score <= low) low_col else if (score <= high) mid_col else high_col
  tags$span(
    style = paste0("background:", col,
                   ";color:white;padding:4px 10px;border-radius:12px;",
                   "font-size:18px;font-weight:bold;"),
    score
  )
}

# ==============================================================================
# UI
# ==============================================================================
ui <- dashboardPage(
  skin = "blue",

  dashboardHeader(
    title = span(icon("heartbeat"), " AF QSP v2.0"),
    titleWidth = 260
  ),

  dashboardSidebar(
    width = 260,
    sidebarMenu(
      id = "main_tabs",
      menuItem("환자 프로파일",         tabName = "tab_patient",   icon = icon("user-md")),
      menuItem("약동학 PK",             tabName = "tab_pk",        icon = icon("flask")),
      menuItem("부정맥 PD",             tabName = "tab_ep",        icon = icon("wave-square")),
      menuItem("항응고 & 뇌졸중",       tabName = "tab_anticoag",  icon = icon("tint")),
      menuItem("시나리오 비교",         tabName = "tab_scenarios", icon = icon("chart-bar")),
      menuItem("바이오마커",            tabName = "tab_biomarker", icon = icon("microscope"))
    ),
    hr(),
    div(style = "padding:10px 15px; color:#aaa; font-size:11px;",
      strong("AF QSP Dashboard v2.0"), br(),
      "mrgsolve · ggplot2 · Shiny", br(),
      "Calibrated: AFFIRM · RACE II ·", br(),
      "ARISTOTLE · ROCKET-AF · RE-LY"
    )
  ),

  dashboardBody(
    tags$head(tags$style(HTML("
      .content-wrapper, .right-side { background-color: #F0F4F8; }
      .box { border-radius: 8px; box-shadow: 0 2px 8px rgba(0,0,0,0.08); }
      .score-card {
        background: white; border-radius: 10px; padding: 14px 10px;
        margin: 6px 2px; box-shadow: 0 2px 6px rgba(0,0,0,0.10);
        text-align: center;
      }
      .sc-val { font-size: 26px; font-weight: 700; color: #2C3E50; }
      .sc-lbl { font-size: 11px; color: #7F8C8D; margin-top: 3px; }
      .pill-green  { background:#27AE60; color:white; border-radius:12px;
                     padding:3px 10px; font-size:12px; font-weight:bold; }
      .pill-yellow { background:#F39C12; color:white; border-radius:12px;
                     padding:3px 10px; font-size:12px; font-weight:bold; }
      .pill-red    { background:#E74C3C; color:white; border-radius:12px;
                     padding:3px 10px; font-size:12px; font-weight:bold; }
      .info-callout {
        background:#EBF5FB; border-left:4px solid #2980B9;
        padding:8px 12px; border-radius:0 6px 6px 0;
        font-size:12px; margin-top:8px;
      }
      .drug-header { font-weight:700; color:#1A5276; margin-top:10px; }
      table.dataTable tbody tr:hover { background-color: #EBF5FB !important; }
      .btn-sim { margin-top:8px; }
    ")))
    ,

    tabItems(

      # ====================================================================
      # TAB 1: 환자 프로파일 (Patient Profile)
      # ====================================================================
      tabItem(tabName = "tab_patient",
        fluidRow(
          # --- Demographics ---
          box(
            title = tagList(icon("user"), " 환자 기본정보 (Demographics)"),
            width = 4, status = "primary", solidHeader = TRUE,
            sliderInput("age", "나이 (Age, years)", 18, 90, 68, 1),
            radioButtons("sex", "성별 (Sex)",
              choices = c("남성 (Male)" = "male", "여성 (Female)" = "female"),
              selected = "male", inline = TRUE),
            radioButtons("af_type", "심방세동 유형 (AF Type)",
              choices = c(
                "발작성 (Paroxysmal)"     = "paroxysmal",
                "지속성 (Persistent)"    = "persistent",
                "장기지속 (Long-standing)" = "longstanding"
              ), selected = "persistent"),
            sliderInput("lvef", "좌심실 박출률 LVEF (%)", 15, 75, 55, 1),
            sliderInput("la_diam_inp", "좌심방 직경 LA (mm)", 30, 65, 44, 1),
            sliderInput("egfr", "eGFR (mL/min/1.73m²)", 10, 130, 72, 1),
            sliderInput("sbp", "수축기 혈압 SBP (mmHg)", 80, 200, 130, 5)
          ),

          # --- Comorbidities ---
          box(
            title = tagList(icon("stethoscope"), " 동반질환 & 위험인자 (Comorbidities)"),
            width = 4, status = "warning", solidHeader = TRUE,
            h5(strong("CHA₂DS₂-VASc 인자:")),
            checkboxInput("htn",          "고혈압 (Hypertension)",           TRUE),
            checkboxInput("dm",           "당뇨병 (Diabetes Mellitus)",       FALSE),
            checkboxInput("hf",           "심부전 (Heart Failure)",           FALSE),
            checkboxInput("prior_stroke", "뇌졸중/TIA 기왕력 (Prior Stroke/TIA)", FALSE),
            checkboxInput("pad",          "혈관질환 (PAD/Vascular Disease)",  FALSE),
            hr(),
            h5(strong("HAS-BLED 인자:")),
            checkboxInput("hbled_htn_unc",   "조절안된 고혈압 SBP>160",         FALSE),
            checkboxInput("hbled_renal",     "신기능 이상 (eGFR<60 or dialysis)", FALSE),
            checkboxInput("hbled_liver",     "간기능 이상 (Liver disease)",      FALSE),
            checkboxInput("hbled_bleed",     "출혈력 (Prior bleeding)",          FALSE),
            checkboxInput("hbled_labile",    "불안정 INR (Labile INR >TTR<60%)", FALSE),
            checkboxInput("hbled_drugs_alc", "약물/알코올 (Drugs/Alcohol ≥8u/wk)",FALSE)
          ),

          # --- Risk scores ---
          box(
            title = tagList(icon("calculator"), " 위험도 점수 (Risk Scores)"),
            width = 4, status = "success", solidHeader = TRUE,
            h4("CHA₂DS₂-VASc"),
            div(class = "score-card",
              div(class = "sc-val", textOutput("chads_num", inline = TRUE)),
              div(class = "sc-lbl", "점수"),
              br(),
              div(class = "sc-val", style = "font-size:18px; color:#E74C3C;",
                  textOutput("chads_risk_pct", inline = TRUE)),
              div(class = "sc-lbl", "연간 뇌졸중 위험도 (%/yr)")
            ),
            uiOutput("chads_badge"),
            hr(),
            h4("HAS-BLED"),
            div(class = "score-card",
              div(class = "sc-val", textOutput("hasbled_num", inline = TRUE)),
              div(class = "sc-lbl", "점수 (≥3 = 고위험)"),
              uiOutput("hasbled_badge")
            ),
            hr(),
            h4("항응고 권고 (Anticoagulation Recommendation)"),
            uiOutput("anticoag_rec"),
            hr(),
            h4("신기능별 용량 조정 (Renal Dosing)"),
            uiOutput("renal_guidance")
          )
        ),

        # --- CHA2DS2-VASc detailed breakdown table ---
        fluidRow(
          box(
            title = "CHA₂DS₂-VASc 세부 점수 구성 (Score Breakdown)",
            width = 12, status = "info", collapsible = TRUE, collapsed = FALSE,
            tableOutput("chads_breakdown_table")
          )
        )
      ),

      # ====================================================================
      # TAB 2: 약동학 PK
      # ====================================================================
      tabItem(tabName = "tab_pk",
        fluidRow(
          # --- Drug selection panel ---
          box(
            title = tagList(icon("pills"), " 투여 약물 설정 (Drug Dosing)"),
            width = 3, status = "info", solidHeader = TRUE,

            h5(class = "drug-header", "심박수 조절 (Rate Control)"),
            selectInput("rate_drug", "약물 선택:",
              choices = c("없음 (None)"      = "none",
                          "메토프롤롤 (Metoprolol)" = "metro",
                          "딜티아젬 (Diltiazem)"    = "dilt",
                          "디곡신 (Digoxin)"        = "digo"),
              selected = "metro"),
            conditionalPanel("input.rate_drug != 'none'",
              uiOutput("rate_dose_ui")
            ),

            h5(class = "drug-header", "심율동 전환 (Rhythm Control)"),
            selectInput("rhythm_drug", "약물 선택:",
              choices = c("없음 (None)"          = "none",
                          "아미오다론 (Amiodarone)"   = "amio",
                          "드로네다론 (Dronedarone)"  = "dron",
                          "플레카이나이드 (Flecainide)" = "flec"),
              selected = "none"),
            conditionalPanel("input.rhythm_drug != 'none'",
              uiOutput("rhythm_dose_ui")
            ),

            h5(class = "drug-header", "항응고제 (Anticoagulation)"),
            selectInput("anticoag_drug", "약물 선택:",
              choices = c("없음 (None)"       = "none",
                          "아픽사반 (Apixaban)"   = "apix",
                          "리바록사반 (Rivaroxaban)" = "riva",
                          "와파린 (Warfarin)"     = "warf"),
              selected = "apix"),
            conditionalPanel("input.anticoag_drug != 'none'",
              uiOutput("ac_dose_ui")
            ),

            hr(),
            sliderInput("pk_tmax", "시뮬레이션 기간 (days)", 14, 365, 180, 14),
            actionButton("run_pk", "PK 시뮬레이션 실행",
              icon = icon("play"), class = "btn-primary btn-block btn-sim")
          ),

          # --- PK plots ---
          box(
            title = tagList(icon("chart-line"), " 약물 농도-시간 곡선 (Concentration-Time)"),
            width = 9, status = "primary", solidHeader = TRUE,
            tabsetPanel(
              tabPanel("심율동 조절 (Rhythm)",
                plotOutput("pk_rhythm_plot", height = "320px"),
                div(class = "info-callout",
                  strong("아미오다론:"), " t½β≈50일 (지방 축적), 치료농도 1–2.5 µg/mL | ",
                  strong("드로네다론:"), " t½≈24h, 치료농도 50–150 ng/mL | ",
                  strong("플레카이나이드:"), " t½≈20h, 치료농도 200–1000 ng/mL"
                )
              ),
              tabPanel("심박수 조절 (Rate)",
                plotOutput("pk_rate_plot", height = "320px"),
                div(class = "info-callout",
                  strong("메토프롤롤:"), " t½≈3.5h, 치료농도 20–100 ng/mL | ",
                  strong("딜티아젬:"), " t½≈3h, 치료농도 50–200 ng/mL | ",
                  strong("디곡신:"), " t½≈36h, 치료농도 0.8–2.0 ng/mL"
                )
              ),
              tabPanel("항응고제 (Anticoagulation)",
                plotOutput("pk_ac_plot", height = "320px"),
                div(class = "info-callout",
                  strong("아픽사반:"), " t½≈12h, 치료농도 50–200 ng/mL | ",
                  strong("리바록사반:"), " t½≈9h, 치료농도 50–250 ng/mL | ",
                  strong("와파린:"), " t½≈40h, 목표 INR 2.0–3.0"
                )
              ),
              tabPanel("INR / Anti-Xa",
                plotOutput("pk_inr_antifxa_plot", height = "320px")
              ),
              tabPanel("PK 요약표 (Summary)",
                tableOutput("pk_summary_table")
              )
            )
          )
        )
      ),

      # ====================================================================
      # TAB 3: 부정맥 PD (AF Dynamics)
      # ====================================================================
      tabItem(tabName = "tab_ep",
        fluidRow(
          box(
            title = "PD 시뮬레이션 설정", width = 3, status = "info", solidHeader = TRUE,
            helpText("Tab 1 · 2의 환자·약물 설정을 자동 적용합니다"),
            sliderInput("ep_tmax", "기간 (days)", 30, 730, 365, 30),
            hr(),
            div(class = "info-callout",
              strong("AF Burden 분류:"), br(),
              "< 25%: 발작성 (Paroxysmal)", br(),
              "25–75%: 지속성 (Persistent)", br(),
              "> 75%: 영구성 (Permanent)"
            ),
            hr(),
            actionButton("run_ep", "EP 시뮬레이션 실행",
              icon = icon("play"), class = "btn-primary btn-block btn-sim")
          ),
          box(
            title = tagList(icon("wave-square"), " 심장전기생리학 출력 (Electrophysiology PD)"),
            width = 9, status = "primary", solidHeader = TRUE,
            tabsetPanel(
              tabPanel("AF 부담 (%)",
                plotOutput("ep_af_plot", height = "290px")
              ),
              tabPanel("유효불응기 ERP (ms)",
                plotOutput("ep_erp_plot", height = "290px")
              ),
              tabPanel("심박수 조절 (bpm)",
                plotOutput("ep_hr_plot", height = "290px")
              ),
              tabPanel("QTc 안전성 (ms)",
                plotOutput("ep_qtc_plot", height = "290px")
              ),
              tabPanel("섬유화 진행 (%)",
                plotOutput("ep_fib_plot", height = "290px")
              ),
              tabPanel("ERP-AF 위상도 (Phase)",
                plotOutput("ep_phase_plot", height = "290px")
              )
            )
          )
        )
      ),

      # ====================================================================
      # TAB 4: 항응고 & 뇌졸중
      # ====================================================================
      tabItem(tabName = "tab_anticoag",
        fluidRow(
          box(
            title = "항응고 시뮬레이션 설정", width = 3, status = "info", solidHeader = TRUE,
            helpText("Tab 1 · 2의 설정을 자동 적용합니다"),
            sliderInput("ac_tmax", "기간 (days)", 30, 730, 365, 30),
            hr(),
            div(class = "info-callout",
              strong("ARISTOTLE (Apixaban):"), br(),
              "뇌졸중/SE: 1.27 vs 1.60%/yr", br(),
              "RRR 21%, ARR 0.33%/yr, NNT 303", br(), br(),
              strong("ROCKET-AF (Rivaroxaban):"), br(),
              "뇌졸중/SE: 1.70 vs 2.15%/yr", br(),
              "RRR 21%, NNT 222"
            ),
            hr(),
            actionButton("run_ac", "항응고 시뮬레이션 실행",
              icon = icon("play"), class = "btn-primary btn-block btn-sim")
          ),
          box(
            title = tagList(icon("tint"), " 항응고 & 혈전색전증 위험도"),
            width = 9, status = "danger", solidHeader = TRUE,
            tabsetPanel(
              tabPanel("FXa & Thrombin",
                fluidRow(
                  column(6, plotOutput("ac_fxa_plot",     height = "270px")),
                  column(6, plotOutput("ac_thrombin_plot",height = "270px"))
                )
              ),
              tabPanel("뇌졸중 위험 & Anti-Xa",
                fluidRow(
                  column(6, plotOutput("ac_stroke_plot",  height = "270px")),
                  column(6, plotOutput("ac_antifxa_plot", height = "270px"))
                )
              ),
              tabPanel("INR (와파린)",
                plotOutput("ac_inr_plot", height = "300px")
              ),
              tabPanel("NNT 계산기",
                br(),
                uiOutput("nnt_output")
              )
            )
          )
        )
      ),

      # ====================================================================
      # TAB 5: 시나리오 비교 (Scenario Comparison)
      # ====================================================================
      tabItem(tabName = "tab_scenarios",
        fluidRow(
          box(
            title = "시나리오 비교 설정", width = 3, status = "info", solidHeader = TRUE,
            helpText("6개 표준 치료 시나리오를 일괄 비교합니다"),
            radioButtons("hl_scen", "강조 시나리오:",
              choices = c(
                "무치료 (No Tx)"           = "S1",
                "메토프롤롤 (Metro)"         = "S2",
                "아미오다론 (Amio)"          = "S3",
                "아픽사반 (Apix)"            = "S4",
                "Metro + Apix"             = "S5",
                "Amio + Apix (표준)"        = "S6"
              ), selected = "S6"),
            sliderInput("scen_tmax", "기간 (days)", 90, 730, 365, 90),
            hr(),
            div(class = "info-callout",
              "각 시나리오는 현재 Tab 1의 환자 설정을 기반으로 합니다."
            ),
            hr(),
            actionButton("run_all", "전체 시나리오 실행",
              icon = icon("play-circle"), class = "btn-success btn-block btn-sim"),
            br(),
            downloadButton("dl_results", "결과 CSV 다운로드",
              class = "btn-default btn-block")
          ),
          box(
            title = tagList(icon("chart-bar"), " 치료 시나리오 비교 (Treatment Comparison)"),
            width = 9, status = "primary", solidHeader = TRUE,
            tabsetPanel(
              tabPanel("AF 부담 시계열",
                plotOutput("scen_af_plot", height = "300px")
              ),
              tabPanel("뇌졸중 위험 시계열",
                plotOutput("scen_stroke_plot", height = "300px")
              ),
              tabPanel("최종 요약표 (Summary Table)",
                div(style = "overflow-x:auto;",
                  tableOutput("scen_summary_tbl")
                )
              ),
              tabPanel("막대 비교 (Bar Chart)",
                fluidRow(
                  column(6, plotOutput("scen_bar_af",     height = "280px")),
                  column(6, plotOutput("scen_bar_stroke", height = "280px"))
                )
              ),
              tabPanel("레이더 (Radar)",
                plotOutput("scen_radar_plot", height = "350px")
              )
            )
          )
        )
      ),

      # ====================================================================
      # TAB 6: 바이오마커 (Biomarkers)
      # ====================================================================
      tabItem(tabName = "tab_biomarker",
        fluidRow(
          box(
            title = "바이오마커 설정", width = 3, status = "info", solidHeader = TRUE,
            helpText("Tab 1 · 2의 설정을 자동 적용합니다"),
            sliderInput("bio_tmax", "기간 (days)", 30, 730, 365, 30),
            hr(),
            div(class = "info-callout",
              strong("정상 범위:"), br(),
              "NT-proBNP < 125 pg/mL", br(),
              "CRP < 3 mg/L", br(),
              "D-dimer < 0.5 µg/mL FEU", br(),
              "Troponin I < 0.04 ng/mL", br(),
              "LA diameter < 40 mm", br(),
              "심방 섬유화 < 15%"
            ),
            hr(),
            actionButton("run_bio", "바이오마커 시뮬레이션 실행",
              icon = icon("play"), class = "btn-primary btn-block btn-sim"),
            br(),
            downloadButton("dl_bio_report", "보고서 CSV 다운로드",
              class = "btn-default btn-block")
          ),
          box(
            title = tagList(icon("microscope"), " 바이오마커 대시보드 (Biomarker Dashboard)"),
            width = 9, status = "success", solidHeader = TRUE,
            # KPI strip
            fluidRow(
              column(3, div(class="score-card",
                div(class="sc-val", textOutput("kpi_bnp")),
                div(class="sc-lbl", "NT-proBNP (pg/mL)"),
                uiOutput("kpi_bnp_pill")
              )),
              column(3, div(class="score-card",
                div(class="sc-val", textOutput("kpi_crp")),
                div(class="sc-lbl", "CRP (mg/L)"),
                uiOutput("kpi_crp_pill")
              )),
              column(3, div(class="score-card",
                div(class="sc-val", textOutput("kpi_ddimer")),
                div(class="sc-lbl", "D-dimer (µg/mL)"),
                uiOutput("kpi_ddimer_pill")
              )),
              column(3, div(class="score-card",
                div(class="sc-val", textOutput("kpi_trop")),
                div(class="sc-lbl", "Troponin I (ng/mL)"),
                uiOutput("kpi_trop_pill")
              ))
            ),
            hr(),
            tabsetPanel(
              tabPanel("NT-proBNP & CRP",
                plotOutput("bio_bnp_crp",    height = "270px")
              ),
              tabPanel("D-dimer & Troponin",
                plotOutput("bio_ddimer_trop", height = "270px")
              ),
              tabPanel("구조적 리모델링 (LA / Fibrosis)",
                plotOutput("bio_struct",      height = "270px")
              ),
              tabPanel("신경호르몬 (AngII / ROS / NE)",
                plotOutput("bio_neuro",       height = "270px")
              ),
              tabPanel("신호등 요약 (Traffic Light)",
                plotOutput("bio_traffic",     height = "300px")
              )
            )
          )
        )
      )

    ) # end tabItems
  )   # end dashboardBody
)     # end dashboardPage

# ==============================================================================
# SERVER
# ==============================================================================
server <- function(input, output, session) {

  # ---------- Dynamic dose UIs ----------
  output$rate_dose_ui <- renderUI({
    switch(input$rate_drug,
      metro = sliderInput("dose_rate", "메토프롤롤 용량 (mg BID)", 12.5, 200, 50, 12.5),
      dilt  = sliderInput("dose_rate", "딜티아젬 용량 (mg TID)",   30,   360, 120, 30),
      digo  = sliderInput("dose_rate", "디곡신 용량 (mg QD)",      0.0625, 0.375, 0.125, 0.0625),
      NULL
    )
  })
  output$rhythm_dose_ui <- renderUI({
    switch(input$rhythm_drug,
      amio = sliderInput("dose_rhythm", "아미오다론 유지용량 (mg/day)", 50, 400, 200, 50),
      dron = sliderInput("dose_rhythm", "드로네다론 용량 (mg BID)",    200, 800, 400, 200),
      flec = sliderInput("dose_rhythm", "플레카이나이드 용량 (mg BID)", 50, 200, 100, 50),
      NULL
    )
  })
  output$ac_dose_ui <- renderUI({
    switch(input$anticoag_drug,
      apix = sliderInput("dose_ac", "아픽사반 용량 (mg BID)", 2.5, 10, 5, 2.5),
      riva = sliderInput("dose_ac", "리바록사반 용량 (mg QD)", 10, 20, 20, 5),
      warf = sliderInput("dose_ac", "와파린 용량 (mg QD)",     1,  10, 5,  0.5),
      NULL
    )
  })

  # ---------- Helpers to safely get dose ----------
  get_dose_rate   <- reactive({ if (is.null(input$dose_rate))   0 else input$dose_rate })
  get_dose_rhythm <- reactive({ if (is.null(input$dose_rhythm)) 0 else input$dose_rhythm })
  get_dose_ac     <- reactive({ if (is.null(input$dose_ac))     0 else input$dose_ac })

  # ---------- CHA2DS2-VASc ----------
  chads_val <- reactive({
    calc_chads2(input$age, input$htn, input$dm, input$hf,
                input$prior_stroke, input$pad, input$sex == "female")
  })
  hasbled_val <- reactive({
    calc_hasbled(input$age, input$hbled_htn_unc, input$hbled_renal,
                 input$hbled_liver, input$hbled_bleed,
                 input$hbled_labile, input$hbled_drugs_alc)
  })

  output$chads_num      <- renderText({ chads_val() })
  output$chads_risk_pct <- renderText({ paste0(chads_stroke_pct(chads_val()), "%") })
  output$hasbled_num    <- renderText({ hasbled_val() })

  output$chads_badge <- renderUI({
    sc <- chads_val()
    if (sc >= 2)      tags$span(class="pill-red",    "항응고 적응증 (Indicated)")
    else if (sc == 1) tags$span(class="pill-yellow",  "고려 (Consider)")
    else              tags$span(class="pill-green",   "저위험 (Low Risk)")
  })
  output$hasbled_badge <- renderUI({
    s <- hasbled_val()
    if (s >= 3)       tags$span(class="pill-red",    "고위험 출혈 (High Bleeding)")
    else if (s == 2)  tags$span(class="pill-yellow",  "중위험 (Moderate)")
    else              tags$span(class="pill-green",   "저위험 (Low)")
  })

  output$anticoag_rec <- renderUI({
    sc <- chads_val()
    if (sc >= 2)
      div(class="info-callout",
        icon("check-circle", style="color:#27AE60"),
        strong(" 항응고 치료 권고 (Anticoagulation Recommended)"), br(),
        "DOAC 우선 권고: Apixaban 5mg BID (ARISTOTLE 연구)", br(),
        "INR 관리 불안정시: 와파린 (목표 INR 2–3)"
      )
    else if (sc == 1)
      div(class="info-callout",
        icon("exclamation-triangle", style="color:#F39C12"),
        strong(" 개별 위험-이익 평가 후 고려"), br(),
        "남성 CHA₂DS₂-VASc=1 또는 여성 =2: 항응고 고려 가능"
      )
    else
      div(class="info-callout",
        icon("info-circle", style="color:#2980B9"),
        strong(" 항응고 불필요 (No anticoagulation)"), br(),
        "남성 CHA₂DS₂-VASc=0 / 여성=1: 항혈소판제도 불필요"
      )
  })

  output$renal_guidance <- renderUI({
    eg <- input$egfr
    txt <- if (eg < 15)
      list("red", "eGFR <15: DOAC 금기. 와파린 또는 투석 전 중단 고려")
    else if (eg < 30)
      list("orange", "eGFR 15-29: 다비가트란 금기. 아픽사반 2.5mg BID 우선")
    else if (eg < 50)
      list("#B8860B", "eGFR 30-49: 80세 이상 또는 체중 60kg 이하시 아픽사반 감량 고려")
    else
      list("#27AE60", "eGFR ≥50: 표준용량 적용 가능")
    div(style=paste0("color:", txt[[1]], "; font-weight:bold;"), txt[[2]])
  })

  output$chads_breakdown_table <- renderTable({
    data.frame(
      "위험인자" = c("뇌졸중/TIA 기왕력", "나이 ≥75세",  "나이 65-74세",
                    "고혈압", "당뇨병", "심부전", "혈관질환", "여성"),
      "점수"    = c(2, 2, 1, 1, 1, 1, 1, 1),
      "해당여부" = c(
        ifelse(input$prior_stroke, "해당 ✓", "해당 없음"),
        ifelse(input$age >= 75, "해당 ✓", "해당 없음"),
        ifelse(input$age >= 65 & input$age < 75, "해당 ✓", "해당 없음"),
        ifelse(input$htn, "해당 ✓", "해당 없음"),
        ifelse(input$dm,  "해당 ✓", "해당 없음"),
        ifelse(input$hf,  "해당 ✓", "해당 없음"),
        ifelse(input$pad, "해당 ✓", "해당 없음"),
        ifelse(input$sex == "female", "해당 ✓", "해당 없음")
      ),
      "기여 점수" = c(
        ifelse(input$prior_stroke, 2, 0),
        ifelse(input$age >= 75, 2, 0),
        ifelse(input$age >= 65 & input$age < 75, 1, 0),
        ifelse(input$htn, 1, 0),
        ifelse(input$dm,  1, 0),
        ifelse(input$hf,  1, 0),
        ifelse(input$pad, 1, 0),
        ifelse(input$sex == "female", 1, 0)
      ),
      check.names = FALSE
    )
  }, striped = TRUE, bordered = TRUE, spacing = "s")

  # ============================================================
  # TAB 2 PK — simulation
  # ============================================================
  pk_dat <- eventReactive(input$run_pk, {
    run_sim(
      af_type     = input$af_type,
      rate_drug   = input$rate_drug,   dose_rate   = get_dose_rate(),
      rhythm_drug = input$rhythm_drug, dose_rhythm = get_dose_rhythm(),
      anticoag    = input$anticoag_drug, dose_ac   = get_dose_ac(),
      t_days      = input$pk_tmax,
      age         = input$age,
      chads_score = chads_val(),
      htn         = input$htn, hf = input$hf, dm = input$dm,
      la_mm       = input$la_diam_inp
    )
  }, ignoreNULL = FALSE)

  output$pk_rhythm_plot <- renderPlot({
    df <- pk_dat(); if (is.null(df)) return(NULL)
    rhythm <- input$rhythm_drug
    if (rhythm == "amio") {
      ggplot(df, aes(x=time_days, y=Cp_AMIO_out)) +
        annotate("rect", xmin=-Inf,xmax=Inf, ymin=1.0,ymax=2.5,
                 alpha=0.2, fill="#27AE60") +
        geom_line(color="#1A5276", linewidth=1.3) +
        annotate("text", x=input$pk_tmax*0.75, y=1.75,
                 label="치료 농도 1–2.5 µg/mL", size=3.5, color="darkgreen") +
        labs(title="아미오다론 혈장 농도 (Amiodarone Cp)",
             x="시간 (일)", y="농도 (µg/mL)") +
        theme_bw(base_size=12) + theme(panel.grid.minor=element_blank())
    } else if (rhythm == "dron") {
      ggplot(df, aes(x=time_days, y=Cp_DRON_out)) +
        annotate("rect", xmin=-Inf,xmax=Inf, ymin=50,ymax=150,
                 alpha=0.2, fill="#2ECC71") +
        geom_line(color="#117A65", linewidth=1.3) +
        labs(title="드로네다론 혈장 농도 (Dronedarone Cp)",
             x="시간 (일)", y="농도 (ng/mL)") +
        theme_bw(base_size=12) + theme(panel.grid.minor=element_blank())
    } else if (rhythm == "flec") {
      ggplot(df, aes(x=time_days, y=Cp_FLEC_out)) +
        annotate("rect", xmin=-Inf,xmax=Inf, ymin=200,ymax=1000,
                 alpha=0.2, fill="#82E0AA") +
        geom_line(color="#1D8348", linewidth=1.3) +
        labs(title="플레카이나이드 혈장 농도 (Flecainide Cp)",
             x="시간 (일)", y="농도 (ng/mL)") +
        theme_bw(base_size=12) + theme(panel.grid.minor=element_blank())
    } else {
      ggplot() + annotate("text",x=0.5,y=0.5,label="심율동 조절 약물을 선택하세요",size=5) +
        theme_void()
    }
  })

  output$pk_rate_plot <- renderPlot({
    df <- pk_dat(); if (is.null(df)) return(NULL)
    rd <- input$rate_drug
    if (rd == "metro") {
      ggplot(df, aes(x=time_days, y=Cp_METRO_ng_out)) +
        annotate("rect", xmin=-Inf,xmax=Inf, ymin=20,ymax=100,
                 alpha=0.2, fill="#F39C12") +
        geom_line(color="#B9770E", linewidth=1.3) +
        labs(title="메토프롤롤 혈장 농도 (Metoprolol Cp)",
             x="시간 (일)", y="농도 (ng/mL)") +
        theme_bw(base_size=12) + theme(panel.grid.minor=element_blank())
    } else if (rd == "dilt") {
      ggplot(df, aes(x=time_days, y=Cp_DILT_ng_out)) +
        annotate("rect", xmin=-Inf,xmax=Inf, ymin=50,ymax=200,
                 alpha=0.2, fill="#F1948A") +
        geom_line(color="#922B21", linewidth=1.3) +
        labs(title="딜티아젬 혈장 농도 (Diltiazem Cp)",
             x="시간 (일)", y="농도 (ng/mL)") +
        theme_bw(base_size=12) + theme(panel.grid.minor=element_blank())
    } else if (rd == "digo") {
      ggplot(df, aes(x=time_days, y=Cp_DIGO_ng_out)) +
        annotate("rect", xmin=-Inf,xmax=Inf, ymin=0.8,ymax=2.0,
                 alpha=0.2, fill="#A9CCE3") +
        geom_line(color="#1F618D", linewidth=1.3) +
        labs(title="디곡신 혈장 농도 (Digoxin Cp)",
             x="시간 (일)", y="농도 (ng/mL)") +
        theme_bw(base_size=12) + theme(panel.grid.minor=element_blank())
    } else {
      ggplot() + annotate("text",x=0.5,y=0.5,label="심박수 조절 약물을 선택하세요",size=5) +
        theme_void()
    }
  })

  output$pk_ac_plot <- renderPlot({
    df <- pk_dat(); if (is.null(df)) return(NULL)
    ac <- input$anticoag_drug
    if (ac == "apix") {
      ggplot(df, aes(x=time_days, y=Cp_APIX_ng)) +
        annotate("rect", xmin=-Inf,xmax=Inf, ymin=50,ymax=200,
                 alpha=0.2, fill="#8E44AD") +
        geom_line(color="#6C3483", linewidth=1.3) +
        labs(title="아픽사반 혈장 농도 (Apixaban Cp)",
             x="시간 (일)", y="농도 (ng/mL)") +
        theme_bw(base_size=12) + theme(panel.grid.minor=element_blank())
    } else if (ac == "riva") {
      ggplot(df, aes(x=time_days, y=Cp_RIVA_ng)) +
        annotate("rect", xmin=-Inf,xmax=Inf, ymin=50,ymax=250,
                 alpha=0.2, fill="#5DADE2") +
        geom_line(color="#1A5276", linewidth=1.3) +
        labs(title="리바록사반 혈장 농도 (Rivaroxaban Cp)",
             x="시간 (일)", y="농도 (ng/mL)") +
        theme_bw(base_size=12) + theme(panel.grid.minor=element_blank())
    } else if (ac == "warf") {
      ggplot(df, aes(x=time_days, y=Cp_WARF_out)) +
        geom_line(color="#E74C3C", linewidth=1.3) +
        labs(title="와파린 혈장 농도 (Warfarin Cp)",
             x="시간 (일)", y="농도 (µg/mL)") +
        theme_bw(base_size=12) + theme(panel.grid.minor=element_blank())
    } else {
      ggplot() + annotate("text",x=0.5,y=0.5,label="항응고제를 선택하세요",size=5) +
        theme_void()
    }
  })

  output$pk_inr_antifxa_plot <- renderPlot({
    df <- pk_dat(); if (is.null(df)) return(NULL)
    ac <- input$anticoag_drug
    if (ac == "warf") {
      ggplot(df, aes(x=time_days, y=INR_proxy)) +
        annotate("rect", xmin=-Inf,xmax=Inf, ymin=2.0,ymax=3.0,
                 alpha=0.25, fill="#27AE60") +
        annotate("text", x=input$pk_tmax*0.75, y=2.5,
                 label="목표 INR 2.0–3.0", size=3.5, color="darkgreen") +
        geom_hline(yintercept=c(2,3), linetype="dashed", color="#27AE60", linewidth=0.7) +
        geom_line(color="#C0392B", linewidth=1.3) +
        scale_y_continuous(limits=c(0.8, 5)) +
        labs(title="INR 추정치 (와파린, Warfarin — indirect PD model)",
             x="시간 (일)", y="INR") +
        theme_bw(base_size=12) + theme(panel.grid.minor=element_blank())
    } else if (ac %in% c("apix","riva")) {
      ggplot(df, aes(x=time_days, y=AntiXa_pct)) +
        annotate("rect", xmin=-Inf,xmax=Inf, ymin=50,ymax=100,
                 alpha=0.15, fill="#8E44AD") +
        geom_line(color="#6C3483", linewidth=1.3) +
        labs(title="Anti-FXa 활성도 억제율 (%)",
             x="시간 (일)", y="FXa 억제율 (%)") +
        theme_bw(base_size=12) + theme(panel.grid.minor=element_blank())
    } else {
      ggplot() + annotate("text",x=0.5,y=0.5,label="항응고제를 선택하세요",size=5)+theme_void()
    }
  })

  output$pk_summary_table <- renderTable({
    df <- pk_dat(); if (is.null(df)) return(NULL)
    last <- tail(df, 1)
    data.frame(
      "약물"    = c("아미오다론","드로네다론","플레카이나이드",
                    "아픽사반","리바록사반","와파린",
                    "메토프롤롤","딜티아젬","디곡신"),
      "최종 Cp" = c(
        round(last$Cp_AMIO_out,   3),
        round(last$Cp_DRON_out,   1),
        round(last$Cp_FLEC_out,   1),
        round(last$Cp_APIX_ng,    1),
        round(last$Cp_RIVA_ng,    1),
        round(last$Cp_WARF_out,   3),
        round(last$Cp_METRO_ng_out,1),
        round(last$Cp_DILT_ng_out, 1),
        round(last$Cp_DIGO_ng_out, 3)
      ),
      "단위" = c("µg/mL","ng/mL","ng/mL","ng/mL","ng/mL",
                 "µg/mL","ng/mL","ng/mL","ng/mL"),
      "치료농도" = c("1.0–2.5","50–150","200–1000",
                    "50–200","50–250","—(INR 2–3)",
                    "20–100","50–200","0.8–2.0"),
      check.names = FALSE
    )
  }, striped=TRUE, bordered=TRUE, spacing="s")

  # ============================================================
  # TAB 3 EP — simulation
  # ============================================================
  ep_dat <- eventReactive(input$run_ep, {
    run_sim(
      af_type     = input$af_type,
      rate_drug   = input$rate_drug,   dose_rate   = get_dose_rate(),
      rhythm_drug = input$rhythm_drug, dose_rhythm = get_dose_rhythm(),
      anticoag    = input$anticoag_drug, dose_ac   = get_dose_ac(),
      t_days      = input$ep_tmax,
      age         = input$age,
      chads_score = chads_val(),
      htn = input$htn, hf = input$hf, dm = input$dm,
      la_mm = input$la_diam_inp
    )
  }, ignoreNULL = FALSE)

  output$ep_af_plot <- renderPlot({
    df <- ep_dat(); if (is.null(df)) return(NULL)
    ggplot(df, aes(x=time_days, y=AF_BURDEN*100)) +
      annotate("rect",xmin=-Inf,xmax=Inf,ymin=75,ymax=100,alpha=0.12,fill="#E74C3C")+
      annotate("rect",xmin=-Inf,xmax=Inf,ymin=25,ymax=75, alpha=0.12,fill="#F39C12")+
      annotate("rect",xmin=-Inf,xmax=Inf,ymin=0, ymax=25, alpha=0.12,fill="#27AE60")+
      annotate("text",x=5,y=88,label="영구성 (Permanent)",  size=3.2,hjust=0,color="#922B21")+
      annotate("text",x=5,y=50,label="지속성 (Persistent)", size=3.2,hjust=0,color="#D35400")+
      annotate("text",x=5,y=12,label="발작성 (Paroxysmal)", size=3.2,hjust=0,color="#1E8449")+
      geom_line(color="#1A5276", linewidth=1.4) +
      scale_y_continuous(limits=c(0,100)) +
      labs(title="심방세동 부담 추이 (AF Burden Over Time)",
           x="시간 (일)", y="AF 부담 (%)") +
      theme_bw(base_size=12) + theme(panel.grid.minor=element_blank())
  })

  output$ep_erp_plot <- renderPlot({
    df <- ep_dat(); if (is.null(df)) return(NULL)
    ggplot(df, aes(x=time_days, y=ERP)) +
      geom_hline(yintercept=200, linetype="dashed", color="#2980B9", linewidth=0.9)+
      annotate("text", x=input$ep_tmax*0.65, y=203,
               label="재진입 임계값 (Reentry threshold, 200ms)",
               size=3.2, color="#2980B9") +
      geom_line(color="#1E8449", linewidth=1.4) +
      labs(title="심방 유효불응기 (Atrial ERP)",
           x="시간 (일)", y="ERP (ms)") +
      theme_bw(base_size=12) + theme(panel.grid.minor=element_blank())
  })

  output$ep_hr_plot <- renderPlot({
    df <- ep_dat(); if (is.null(df)) return(NULL)
    ggplot(df, aes(x=time_days, y=HR_AF)) +
      annotate("rect",xmin=-Inf,xmax=Inf,ymin=60,ymax=110,alpha=0.15,fill="#27AE60")+
      geom_line(color="#E74C3C", linewidth=1.4) +
      annotate("text", x=5, y=85, label="목표 심박수 60–110 bpm",
               size=3.5, hjust=0, color="#1E8449") +
      labs(title="심방세동 중 심실박동수 (Ventricular Rate in AF)",
           x="시간 (일)", y="심박수 (bpm)") +
      theme_bw(base_size=12) + theme(panel.grid.minor=element_blank())
  })

  output$ep_qtc_plot <- renderPlot({
    df <- ep_dat(); if (is.null(df)) return(NULL)
    ggplot(df, aes(x=time_days, y=QTc)) +
      annotate("rect",xmin=-Inf,xmax=Inf,ymin=470,ymax=Inf,  alpha=0.15,fill="#E74C3C")+
      annotate("rect",xmin=-Inf,xmax=Inf,ymin=440,ymax=470,  alpha=0.10,fill="#F39C12")+
      annotate("text",x=5,y=477,label="QTc >470ms: TdP 고위험", size=3.2,hjust=0,color="red")+
      annotate("text",x=5,y=453,label="QTc 440–470ms: 모니터링", size=3.2,hjust=0,color="orange")+
      geom_line(color="#8E44AD", linewidth=1.4) +
      labs(title="QTc 간격 안전성 모니터링 (QTc Safety)",
           x="시간 (일)", y="QTc (ms)") +
      theme_bw(base_size=12) + theme(panel.grid.minor=element_blank())
  })

  output$ep_fib_plot <- renderPlot({
    df <- ep_dat(); if (is.null(df)) return(NULL)
    ggplot(df, aes(x=time_days, y=Fibrosis*100)) +
      annotate("rect",xmin=-Inf,xmax=Inf,ymin=35,ymax=Inf,  alpha=0.12,fill="#E74C3C")+
      annotate("rect",xmin=-Inf,xmax=Inf,ymin=15,ymax=35,   alpha=0.10,fill="#F39C12")+
      annotate("text",x=5,y=38,label="중증 섬유화 (Severe)",  size=3.2,hjust=0,color="red")+
      annotate("text",x=5,y=25,label="중등도 (Moderate)", size=3.2,hjust=0,color="orange")+
      geom_line(color="#784212", linewidth=1.4) +
      labs(title="심방 섬유화 진행 (Atrial Fibrosis)",
           x="시간 (일)", y="섬유화 (%)") +
      theme_bw(base_size=12) + theme(panel.grid.minor=element_blank())
  })

  output$ep_phase_plot <- renderPlot({
    df <- ep_dat(); if (is.null(df)) return(NULL)
    df_sub <- df[seq(1, nrow(df), by=max(1,floor(nrow(df)/400))), ]
    ggplot(df_sub, aes(x=ERP, y=AF_BURDEN*100, color=time_days)) +
      geom_point(size=1.2, alpha=0.75) +
      scale_color_gradient(low="#E74C3C", high="#2980B9", name="Day") +
      geom_vline(xintercept=200, linetype="dashed", color="gray50") +
      labs(title="ERP–AF 위상도 (Phase Portrait)",
           x="심방 ERP (ms)", y="AF 부담 (%)") +
      theme_bw(base_size=12) + theme(panel.grid.minor=element_blank())
  })

  # ============================================================
  # TAB 4: Anticoagulation
  # ============================================================
  ac_dat <- eventReactive(input$run_ac, {
    run_sim(
      af_type     = input$af_type,
      rate_drug   = input$rate_drug,   dose_rate   = get_dose_rate(),
      rhythm_drug = input$rhythm_drug, dose_rhythm = get_dose_rhythm(),
      anticoag    = input$anticoag_drug, dose_ac   = get_dose_ac(),
      t_days      = input$ac_tmax,
      age         = input$age,
      chads_score = chads_val(),
      htn = input$htn, hf = input$hf, dm = input$dm,
      la_mm = input$la_diam_inp
    )
  }, ignoreNULL = FALSE)

  output$ac_fxa_plot <- renderPlot({
    df <- ac_dat(); if (is.null(df)) return(NULL)
    ggplot(df, aes(x=time_days, y=FXa)) +
      geom_hline(yintercept=1, linetype="dashed", color="gray60") +
      geom_line(color="#C0392B", linewidth=1.3) +
      labs(title="FXa 활성도 (Factor Xa Activity)",
           x="시간 (일)", y="FXa (상대단위)") +
      theme_bw(base_size=11) + theme(panel.grid.minor=element_blank())
  })

  output$ac_thrombin_plot <- renderPlot({
    df <- ac_dat(); if (is.null(df)) return(NULL)
    ggplot(df, aes(x=time_days, y=Thrombin)) +
      geom_hline(yintercept=1, linetype="dashed", color="gray60") +
      geom_line(color="#922B21", linewidth=1.3) +
      labs(title="트롬빈 활성도 (Thrombin Activity)",
           x="시간 (일)", y="Thrombin (상대단위)") +
      theme_bw(base_size=11) + theme(panel.grid.minor=element_blank())
  })

  output$ac_stroke_plot <- renderPlot({
    df <- ac_dat(); if (is.null(df)) return(NULL)
    ggplot(df, aes(x=time_days, y=STROKE_RISK)) +
      annotate("rect",xmin=-Inf,xmax=Inf,ymin=4,ymax=Inf,  alpha=0.12,fill="#E74C3C")+
      annotate("rect",xmin=-Inf,xmax=Inf,ymin=2,ymax=4,    alpha=0.10,fill="#F39C12")+
      geom_line(color="#1A5276", linewidth=1.3) +
      labs(title="연간 뇌졸중 위험도 (Annual Stroke Risk)",
           x="시간 (일)", y="뇌졸중 위험 (%/년)") +
      theme_bw(base_size=11) + theme(panel.grid.minor=element_blank())
  })

  output$ac_antifxa_plot <- renderPlot({
    df <- ac_dat(); if (is.null(df)) return(NULL)
    ggplot(df, aes(x=time_days, y=AntiXa_pct)) +
      annotate("rect",xmin=-Inf,xmax=Inf,ymin=50,ymax=100,alpha=0.15,fill="#8E44AD")+
      geom_line(color="#6C3483", linewidth=1.3) +
      labs(title="Anti-FXa 억제율 (DOAC)",
           x="시간 (일)", y="FXa 억제율 (%)") +
      theme_bw(base_size=11) + theme(panel.grid.minor=element_blank())
  })

  output$ac_inr_plot <- renderPlot({
    df <- ac_dat(); if (is.null(df)) return(NULL)
    ggplot(df, aes(x=time_days, y=INR_proxy)) +
      annotate("rect",xmin=-Inf,xmax=Inf,ymin=2,ymax=3,alpha=0.25,fill="#27AE60")+
      geom_hline(yintercept=c(2,3), linetype="dashed", color="#27AE60", linewidth=0.7)+
      geom_line(color="#E74C3C", linewidth=1.3) +
      scale_y_continuous(limits=c(0.5,6)) +
      labs(title="INR 추정치 (와파린 간접 PD 모델)",
           x="시간 (일)", y="INR") +
      theme_bw(base_size=12) + theme(panel.grid.minor=element_blank())
  })

  output$nnt_output <- renderUI({
    df <- ac_dat(); if (is.null(df)) return(NULL)
    last <- tail(df, 1)
    base_risk <- last$STROKE_RISK[1] / 100
    # ARISTOTLE apixaban RRR=21%, ROCKET-AF rivaroxaban RRR=21%, warfarin (ref)
    ac  <- input$anticoag_drug
    rrr <- switch(ac, apix=0.21, riva=0.21, warf=0.00, 0.21)
    arr <- base_risk * rrr
    nnt <- if (arr > 0) round(1 / arr) else "N/A"
    trial_ref <- switch(ac,
      apix = "ARISTOTLE (Granger et al. NEJM 2011;365:981)",
      riva = "ROCKET-AF (Patel et al. NEJM 2011;365:883)",
      warf = "와파린은 기준 치료 (Reference arm; RRR vs no anticoag ~64%)",
      "ARISTOTLE / ROCKET-AF"
    )
    tagList(
      div(class="info-callout",
        strong("NNT 계산기 (Number Needed to Treat)"), br(), br(),
        tags$table(style="width:100%;",
          tags$tr(
            tags$td(style="padding:4px;", "기저 뇌졸중 위험 (Baseline stroke risk):"),
            tags$td(style="padding:4px; font-weight:bold;",
                    paste0(round(base_risk*100,2), "%/year"))
          ),
          tags$tr(
            tags$td(style="padding:4px;", "상대위험감소 (RRR):"),
            tags$td(style="padding:4px; font-weight:bold;", paste0(rrr*100, "%"))
          ),
          tags$tr(
            tags$td(style="padding:4px;", "절대위험감소 (ARR):"),
            tags$td(style="padding:4px; font-weight:bold;",
                    paste0(round(arr*100,3), "%/year"))
          ),
          tags$tr(
            tags$td(style="padding:4px;", "1년 NNT:"),
            tags$td(style="padding:4px; font-weight:bold; color:#E74C3C;",
                    paste0(nnt, "명"))
          )
        ),
        em(paste0("근거: ", trial_ref))
      )
    )
  })

  # ============================================================
  # TAB 5: Scenario Comparison
  # ============================================================
  scen_colors <- c(
    "무치료"           = "#E41A1C",
    "메토프롤롤"        = "#377EB8",
    "아미오다론"        = "#4DAF4A",
    "아픽사반"          = "#984EA3",
    "Metro+Apix"      = "#FF7F00",
    "Amio+Apix (표준)" = "#A65628"
  )
  sc_defs <- list(
    S1 = list(label="무치료",           rate="none",  dr=0,   rhythm="none", rr=0,   ac="none", ar=0),
    S2 = list(label="메토프롤롤",        rate="metro", dr=50,  rhythm="none", rr=0,   ac="none", ar=0),
    S3 = list(label="아미오다론",        rate="none",  dr=0,   rhythm="amio", rr=200, ac="none", ar=0),
    S4 = list(label="아픽사반",          rate="none",  dr=0,   rhythm="none", rr=0,   ac="apix", ar=5),
    S5 = list(label="Metro+Apix",       rate="metro", dr=50,  rhythm="none", rr=0,   ac="apix", ar=5),
    S6 = list(label="Amio+Apix (표준)", rate="none",  dr=0,   rhythm="amio", rr=200, ac="apix", ar=5)
  )

  scen_dat <- eventReactive(input$run_all, {
    withProgress(message="6개 시나리오 시뮬레이션 중...", value=0, {
      out <- lapply(names(sc_defs), function(sid) {
        sc <- sc_defs[[sid]]
        incProgress(1/6, detail=sc$label)
        df <- run_sim(
          af_type     = input$af_type,
          rate_drug   = sc$rate,   dose_rate   = sc$dr,
          rhythm_drug = sc$rhythm, dose_rhythm = sc$rr,
          anticoag    = sc$ac,     dose_ac     = sc$ar,
          t_days      = input$scen_tmax,
          age         = input$age,
          chads_score = chads_val(),
          htn = input$htn, hf = input$hf, dm = input$dm,
          la_mm = input$la_diam_inp
        )
        if (!is.null(df)) { df$scen_id <- sid; df$label <- sc$label }
        df
      })
      bind_rows(out)
    })
  })

  output$scen_af_plot <- renderPlot({
    df <- scen_dat(); if (is.null(df)) return(NULL)
    df$label <- factor(df$label, levels=names(scen_colors))
    ggplot(df, aes(x=time_days, y=AF_BURDEN*100, color=label,
                   linewidth=ifelse(scen_id==input$hl_scen, 2.0, 0.8))) +
      geom_line(alpha=0.88) + scale_linewidth_identity() +
      scale_color_manual(values=scen_colors, name="시나리오") +
      labs(title="AF 부담 추이 — 시나리오 비교",
           x="시간 (일)", y="AF 부담 (%)") +
      theme_bw(base_size=12) +
      theme(legend.position="bottom", panel.grid.minor=element_blank())
  })

  output$scen_stroke_plot <- renderPlot({
    df <- scen_dat(); if (is.null(df)) return(NULL)
    df$label <- factor(df$label, levels=names(scen_colors))
    ggplot(df, aes(x=time_days, y=STROKE_RISK, color=label,
                   linewidth=ifelse(scen_id==input$hl_scen, 2.0, 0.8))) +
      geom_line(alpha=0.88) + scale_linewidth_identity() +
      scale_color_manual(values=scen_colors, name="시나리오") +
      labs(title="연간 뇌졸중 위험도 — 시나리오 비교",
           x="시간 (일)", y="뇌졸중 위험 (%/년)") +
      theme_bw(base_size=12) +
      theme(legend.position="bottom", panel.grid.minor=element_blank())
  })

  output$scen_summary_tbl <- renderTable({
    df <- scen_dat(); if (is.null(df)) return(NULL)
    df %>%
      filter(time_days > max(time_days, na.rm=TRUE) - 1) %>%
      group_by("시나리오"=label) %>%
      summarise(
        "AF 부담 (%)"      = round(mean(AF_BURDEN*100, na.rm=TRUE), 1),
        "ERP (ms)"         = round(mean(ERP,            na.rm=TRUE), 1),
        "심박수 (bpm)"      = round(mean(HR_AF,          na.rm=TRUE), 1),
        "QTc (ms)"         = round(mean(QTc,             na.rm=TRUE), 1),
        "뇌졸중 위험 (%/yr)"= round(mean(STROKE_RISK,    na.rm=TRUE), 2),
        "섬유화 (%)"        = round(mean(Fibrosis*100,   na.rm=TRUE), 1),
        "NT-proBNP (pg/mL)"= round(mean(NT_proBNP_st,   na.rm=TRUE), 0),
        .groups="drop"
      )
  }, striped=TRUE, bordered=TRUE, hover=TRUE, spacing="s")

  output$scen_bar_af <- renderPlot({
    df <- scen_dat(); if (is.null(df)) return(NULL)
    s <- df %>%
      filter(time_days > max(time_days,na.rm=TRUE)-1) %>%
      group_by(label) %>%
      summarise(v=round(mean(AF_BURDEN*100,na.rm=TRUE),1), .groups="drop")
    s$label <- factor(s$label, levels=names(scen_colors))
    ggplot(s, aes(x=label, y=v, fill=label)) +
      geom_col(alpha=0.88) +
      geom_text(aes(label=paste0(v,"%")), vjust=-0.3, size=3.2) +
      scale_fill_manual(values=scen_colors) +
      labs(title="최종 AF 부담", x=NULL, y="AF 부담 (%)") +
      theme_bw(base_size=11) +
      theme(axis.text.x=element_text(angle=30,hjust=1),
            legend.position="none", panel.grid.minor=element_blank())
  })

  output$scen_bar_stroke <- renderPlot({
    df <- scen_dat(); if (is.null(df)) return(NULL)
    s <- df %>%
      filter(time_days > max(time_days,na.rm=TRUE)-1) %>%
      group_by(label) %>%
      summarise(v=round(mean(STROKE_RISK,na.rm=TRUE),2), .groups="drop")
    s$label <- factor(s$label, levels=names(scen_colors))
    ggplot(s, aes(x=label, y=v, fill=label)) +
      geom_col(alpha=0.88) +
      geom_text(aes(label=paste0(v,"%")), vjust=-0.3, size=3.2) +
      scale_fill_manual(values=scen_colors) +
      labs(title="최종 뇌졸중 위험", x=NULL, y="뇌졸중 위험 (%/년)") +
      theme_bw(base_size=11) +
      theme(axis.text.x=element_text(angle=30,hjust=1),
            legend.position="none", panel.grid.minor=element_blank())
  })

  output$scen_radar_plot <- renderPlot({
    df <- scen_dat(); if (is.null(df)) return(NULL)
    sumdf <- df %>%
      filter(time_days > max(time_days,na.rm=TRUE)-1) %>%
      group_by(label) %>%
      summarise(
        AF   = mean(AF_BURDEN*100, na.rm=TRUE),
        HR   = mean(HR_AF,         na.rm=TRUE),
        QTc  = mean(QTc,           na.rm=TRUE),
        Stk  = mean(STROKE_RISK,   na.rm=TRUE),
        Fib  = mean(Fibrosis*100,  na.rm=TRUE),
        BNP  = mean(NT_proBNP_st,  na.rm=TRUE) / 10,
        .groups="drop"
      )
    # Normalize to 0–1 per metric (higher = worse) for radar-like bar chart
    normalize <- function(x) (x - min(x)) / (max(x) - min(x) + 1e-9)
    sumdf_n <- sumdf %>%
      mutate(across(AF:BNP, normalize)) %>%
      pivot_longer(-label, names_to="metric", values_to="score")

    metric_labels <- c(AF="AF부담", HR="심박수", QTc="QTc",
                       Stk="뇌졸중위험", Fib="섬유화", BNP="NT-proBNP")
    sumdf_n$metric_kr <- metric_labels[sumdf_n$metric]
    sumdf_n$label     <- factor(sumdf_n$label, levels=names(scen_colors))

    ggplot(sumdf_n, aes(x=metric_kr, y=score, fill=label)) +
      geom_col(position="dodge", alpha=0.85) +
      scale_fill_manual(values=scen_colors, name="시나리오") +
      scale_y_continuous(labels=scales::percent_format()) +
      labs(title="정규화 지표 비교 (낮을수록 좋음)",
           x=NULL, y="정규화 점수 (0=최고, 1=최저)") +
      theme_bw(base_size=12) +
      theme(legend.position="bottom", panel.grid.minor=element_blank())
  })

  output$dl_results <- downloadHandler(
    filename = function() paste0("AF_QSP_scenarios_", Sys.Date(), ".csv"),
    content  = function(f) {
      df <- scen_dat()
      if (!is.null(df)) write.csv(df, f, row.names=FALSE)
    }
  )

  # ============================================================
  # TAB 6: Biomarkers
  # ============================================================
  bio_dat <- eventReactive(input$run_bio, {
    run_sim(
      af_type     = input$af_type,
      rate_drug   = input$rate_drug,   dose_rate   = get_dose_rate(),
      rhythm_drug = input$rhythm_drug, dose_rhythm = get_dose_rhythm(),
      anticoag    = input$anticoag_drug, dose_ac   = get_dose_ac(),
      t_days      = input$bio_tmax,
      age         = input$age,
      chads_score = chads_val(),
      htn = input$htn, hf = input$hf, dm = input$dm,
      la_mm = input$la_diam_inp
    )
  }, ignoreNULL = FALSE)

  bio_last <- reactive({
    df <- bio_dat(); if (is.null(df)) return(NULL); tail(df, 1)
  })

  pill_status <- function(val, g_max, y_max) {
    if (is.null(val)) return(NULL)
    if      (val <= g_max) tags$span(class="pill-green",  "정상 (Normal)")
    else if (val <= y_max) tags$span(class="pill-yellow", "경계 (Borderline)")
    else                   tags$span(class="pill-red",    "이상 (Abnormal)")
  }

  output$kpi_bnp    <- renderText({ if(!is.null(bio_last())) round(bio_last()$NT_proBNP_st,0) else "—" })
  output$kpi_crp    <- renderText({ if(!is.null(bio_last())) round(bio_last()$CRP_st,2)       else "—" })
  output$kpi_ddimer <- renderText({ if(!is.null(bio_last())) round(bio_last()$Ddimer_st,2)    else "—" })
  output$kpi_trop   <- renderText({ if(!is.null(bio_last())) round(bio_last()$TropI_st,3)     else "—" })

  output$kpi_bnp_pill    <- renderUI({ if(!is.null(bio_last())) pill_status(bio_last()$NT_proBNP_st, 125, 500)  })
  output$kpi_crp_pill    <- renderUI({ if(!is.null(bio_last())) pill_status(bio_last()$CRP_st,       3,   10)   })
  output$kpi_ddimer_pill <- renderUI({ if(!is.null(bio_last())) pill_status(bio_last()$Ddimer_st,    0.5, 1.5)  })
  output$kpi_trop_pill   <- renderUI({ if(!is.null(bio_last())) pill_status(bio_last()$TropI_st,     0.04, 0.1) })

  output$bio_bnp_crp <- renderPlot({
    df <- bio_dat(); if (is.null(df)) return(NULL)
    df_l <- df %>% select(time_days, NT_proBNP_st, CRP_st) %>%
      pivot_longer(-time_days, names_to="marker", values_to="value")
    lbl <- c(NT_proBNP_st="NT-proBNP (pg/mL)", CRP_st="CRP (mg/L)")
    df_l$marker_kr <- lbl[df_l$marker]
    ggplot(df_l, aes(x=time_days, y=value, color=marker)) +
      geom_line(linewidth=1.3) +
      facet_wrap(~marker_kr, scales="free_y") +
      scale_color_manual(values=c(NT_proBNP_st="#E74C3C", CRP_st="#E67E22")) +
      labs(title="심장 및 염증 바이오마커",
           x="시간 (일)", y="농도") +
      theme_bw(base_size=11) +
      theme(legend.position="none", panel.grid.minor=element_blank())
  })

  output$bio_ddimer_trop <- renderPlot({
    df <- bio_dat(); if (is.null(df)) return(NULL)
    df_l <- df %>% select(time_days, Ddimer_st, TropI_st) %>%
      pivot_longer(-time_days, names_to="marker", values_to="value")
    lbl <- c(Ddimer_st="D-dimer (µg/mL FEU)", TropI_st="Troponin I (ng/mL)")
    df_l$marker_kr <- lbl[df_l$marker]
    ggplot(df_l, aes(x=time_days, y=value, color=marker)) +
      geom_line(linewidth=1.3) +
      facet_wrap(~marker_kr, scales="free_y") +
      scale_color_manual(values=c(Ddimer_st="#8E44AD", TropI_st="#C0392B")) +
      labs(title="혈전 및 심근 손상 바이오마커",
           x="시간 (일)", y="농도") +
      theme_bw(base_size=11) +
      theme(legend.position="none", panel.grid.minor=element_blank())
  })

  output$bio_struct <- renderPlot({
    df <- bio_dat(); if (is.null(df)) return(NULL)
    df_l <- df %>%
      mutate(Fib_pct=Fibrosis*100, LA=LA_diam) %>%
      select(time_days, Fib_pct, LA) %>%
      pivot_longer(-time_days, names_to="marker", values_to="value")
    lbl <- c(Fib_pct="심방 섬유화 (%)", LA="좌심방 직경 LA (mm)")
    df_l$marker_kr <- lbl[df_l$marker]
    ggplot(df_l, aes(x=time_days, y=value, color=marker)) +
      geom_line(linewidth=1.3) +
      facet_wrap(~marker_kr, scales="free_y") +
      scale_color_manual(values=c(Fib_pct="#784212", LA="#1F618D")) +
      labs(title="구조적 리모델링 지표",
           x="시간 (일)", y="값") +
      theme_bw(base_size=11) +
      theme(legend.position="none", panel.grid.minor=element_blank())
  })

  output$bio_neuro <- renderPlot({
    df <- bio_dat(); if (is.null(df)) return(NULL)
    df_l <- df %>% select(time_days, AngII, ROS, NE) %>%
      pivot_longer(-time_days, names_to="marker", values_to="value")
    ggplot(df_l, aes(x=time_days, y=value, color=marker)) +
      geom_line(linewidth=1.2) +
      geom_hline(yintercept=1, linetype="dashed", color="gray50") +
      scale_color_manual(
        values=c(AngII="#E74C3C", ROS="#8E44AD", NE="#F39C12"),
        labels=c(AngII="안지오텐신 II", ROS="활성산소 (ROS)", NE="노르에피네프린")
      ) +
      labs(title="신경호르몬 매개인자 (Neurohormonal Mediators)",
           x="시간 (일)", y="상대 단위 (정상=1.0)", color=NULL) +
      theme_bw(base_size=11) +
      theme(legend.position="bottom", panel.grid.minor=element_blank())
  })

  output$bio_traffic <- renderPlot({
    df <- bio_dat(); if (is.null(df)) return(NULL)
    last <- tail(df, 1)
    metrics <- data.frame(
      metric    = c("AF 부담 (%)", "심박수 (bpm)", "QTc (ms)",
                    "뇌졸중 위험 (%/yr)", "섬유화 (%)", "NT-proBNP (pg/mL)",
                    "CRP (mg/L)", "D-dimer (µg/mL)", "Troponin I (ng/mL)"),
      value     = c(
        round(last$AF_BURDEN*100,  1),
        round(last$HR_AF,          0),
        round(last$QTc,            0),
        round(last$STROKE_RISK,    2),
        round(last$Fibrosis*100,   1),
        round(last$NT_proBNP_st,   0),
        round(last$CRP_st,         2),
        round(last$Ddimer_st,      2),
        round(last$TropI_st,       3)
      ),
      green_max  = c(25,  80,  440, 2.0, 15,  125,  3,   0.5, 0.04),
      yellow_max = c(75, 110,  470, 4.0, 35,  500,  10,  1.5, 0.10)
    )
    metrics$status <- with(metrics,
      ifelse(value <= green_max, "정상",
             ifelse(value <= yellow_max, "경계", "이상")))
    metrics$status <- factor(metrics$status, levels=c("정상","경계","이상"))
    metrics$label_txt <- paste0(metrics$metric, "\n", metrics$value)

    ggplot(metrics, aes(x=reorder(metric, as.numeric(status)),
                        y=0.5, fill=status)) +
      geom_tile(height=0.9, color="white", linewidth=1.5) +
      geom_text(aes(label=label_txt), size=3.4, fontface="bold") +
      scale_fill_manual(
        values=c("정상"="#27AE60","경계"="#F39C12","이상"="#E74C3C"),
        name="상태"
      ) +
      coord_flip() +
      labs(title="바이오마커 신호등 현황 (시뮬레이션 종료 시점)",
           x=NULL, y=NULL) +
      theme_minimal(base_size=12) +
      theme(axis.text=element_blank(), axis.ticks=element_blank(),
            panel.grid=element_blank(), legend.position="right")
  })

  output$dl_bio_report <- downloadHandler(
    filename = function() paste0("AF_QSP_biomarker_", Sys.Date(), ".csv"),
    content  = function(f) {
      df <- bio_dat()
      if (!is.null(df)) {
        write.csv(
          df %>% select(time_days, AF_BURDEN, ERP, QTc, HR_AF, STROKE_RISK,
                        Fibrosis, LA_diam,
                        NT_proBNP_st, CRP_st, Ddimer_st, TropI_st,
                        AngII, ROS, NE, IL6, FXa, Thrombin,
                        Cp_AMIO_out, Cp_DRON_out, Cp_FLEC_out,
                        Cp_APIX_ng, Cp_RIVA_ng, Cp_WARF_out, INR_proxy,
                        Cp_METRO_ng_out, Cp_DILT_ng_out, Cp_DIGO_ng_out),
          f, row.names=FALSE
        )
      }
    }
  )

} # end server

# Launch
shinyApp(ui = ui, server = server)
