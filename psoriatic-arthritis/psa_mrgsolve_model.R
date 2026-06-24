## ============================================================
## Psoriatic Arthritis (PsA) – mrgsolve QSP Model
## ============================================================
## 22 ODE compartments covering:
##   Drug PK:  Adalimumab (2-cmt SC+TMDD) · Ixekizumab (SC+TMDD) ·
##             Guselkumab (SC+TMDD) · Upadacitinib (oral 1-cmt) ·
##             Apremilast (oral 1-cmt)
##   Disease PD: IL-17A · TNF-α · IL-23 · Th17 cells ·
##               RANKL/OPG · Skin (PASI) · Joint (DAPSA) ·
##               Bone erosion · Calprotectin · CRP
## Clinical trial calibration notes:
##   Adalimumab  – ADEPT (Ann Rheum Dis 2005): ACR20 57% at wk12
##   Ixekizumab  – SPIRIT-P1/P2 (Ann Rheum Dis 2017): ACR20 62-58%
##   Guselkumab  – DISCOVER-1/2 (Lancet 2020): ACR20 59-64%
##   Upadacitinib– SELECT-PsA 1/2 (Ann Rheum Dis 2021): ACR20 71-71%
##   Apremilast  – PALACE 1-3 (Arthritis Rheum 2014): ACR20 38-40%
## ============================================================

library(mrgsolve)
library(dplyr)
library(ggplot2)
library(tidyr)

## --------------------------------------------------------
## Model code block
## --------------------------------------------------------
psa_code <- '
$PROB
  Psoriatic Arthritis QSP Model
  Compartments: Drug PK (ADA/IXE/GUS/UPA/APR) + Disease PD
  Date: 2026-06-17

$PARAM
  // ---- Patient parameters ----
  BWT     = 85      // body weight (kg)
  BSA     = 2.0     // body surface area (m2)

  // ---- Adalimumab (ADA) PK – 2-cmt SC + TMDD ----
  ka_ADA  = 0.0065  // absorption (1/h) SC, t1/2_abs ~107h
  CL_ADA  = 0.0117  // clearance (L/h)   t1/2 ~14 days
  V1_ADA  = 2.75    // central volume (L)
  V2_ADA  = 1.76    // peripheral volume (L)
  Q_ADA   = 0.0096  // intercompartmental CL (L/h)
  F_ADA   = 0.64    // SC bioavailability
  kon_ADA = 0.096   // ADA-TNFa binding on-rate (1/(nM·h))
  koff_ADA= 0.00048 // ADA-TNFa binding off-rate (1/h)
  kdeg_ADA= 0.0042  // ADA-TNFa complex deg (1/h)
  dose_ADA= 40      // mg SC Q2W

  // ---- Ixekizumab (IXE) PK – 1-cmt SC + TMDD ----
  ka_IXE  = 0.0055  // (1/h)
  CL_IXE  = 0.0119  // (L/h)
  V1_IXE  = 3.8     // (L)
  F_IXE   = 0.60    // SC bioavailability
  kon_IXE = 0.21    // IXE-IL17A binding (1/(nM·h))
  koff_IXE= 0.00021 // (1/h)  Kd ~1 pM
  kdeg_IXE= 0.0052  // (1/h)
  dose_IXE= 80      // mg SC: 160mg wk0 → 80mg Q2W (induction), Q4W (maint)

  // ---- Guselkumab (GUS) PK – 1-cmt SC + TMDD ----
  ka_GUS  = 0.0072  // (1/h)
  CL_GUS  = 0.0082  // (L/h)  t1/2 ~17.5 days
  V1_GUS  = 3.2     // (L)
  F_GUS   = 0.49    // SC bioavailability
  kon_GUS = 0.30    // GUS-IL23p19 (1/(nM·h))
  koff_GUS= 0.00018 // (1/h)  Kd ~0.6 pM
  kdeg_GUS= 0.0038  // (1/h)
  dose_GUS= 100     // mg SC Q8W (after wk0+wk4 loading)

  // ---- Upadacitinib (UPA) PK – 1-cmt oral ----
  ka_UPA  = 0.72    // (1/h) rapid oral absorption
  CL_UPA  = 38.4    // (L/h)
  V1_UPA  = 220     // (L)
  F_UPA   = 0.76
  dose_UPA= 15      // mg QD

  // ---- Apremilast (APR) PK – 1-cmt oral ----
  ka_APR  = 0.65    // (1/h)
  CL_APR  = 10.1    // (L/h)
  V1_APR  = 84      // (L)
  F_APR   = 0.73
  dose_APR= 30      // mg BID

  // ---- Disease PD parameters ----
  // IL-17A dynamics
  kin_IL17  = 0.018   // production rate (nM/h) calibrated to SPIRIT-P1
  kout_IL17 = 0.014   // elimination rate (1/h)
  IL17_base = 1.286   // baseline nM (~8 pg/mL)

  // TNF-α dynamics
  kin_TNF   = 0.022
  kout_TNF  = 0.017
  TNF_base  = 1.294   // baseline nM (~22 pg/mL)

  // IL-23 dynamics
  kin_IL23  = 0.015
  kout_IL23 = 0.012
  IL23_base = 1.25

  // Th17 cell dynamics (fold-change)
  kin_Th17  = 0.008
  kout_Th17 = 0.006
  Th17_base = 1.0

  // RANKL/OPG ratio (bone)
  kin_RANKL = 0.010
  kout_RANKL= 0.008
  RANKL_base= 1.25

  // PASI (skin score 0-72, continuous proxy)
  kin_PASI  = 0.004
  kout_PASI = 0.003
  PASI_base = 18.0

  // DAPSA (joint score 0-164, continuous proxy)
  kin_DAPSA = 0.005
  kout_DAPSA= 0.0038
  DAPSA_base= 28.0

  // CRP (mg/L)
  kin_CRP   = 0.030
  kout_CRP  = 0.025
  CRP_base  = 1.2     // ratio, baseline = 20 mg/L

  // Calprotectin (S100A8/9, µg/mL proxy)
  kin_S100  = 0.012
  kout_S100 = 0.009
  S100_base = 1.3

  // Bone erosion (mTSS proxy, slow dynamics)
  kin_MTSS  = 0.00015
  kout_MTSS = 0.00005
  MTSS_base = 0.0

  // Drug effect parameters
  // Emax/EC50 for each drug on relevant targets
  Emax_ADA_TNF   = 0.92
  EC50_ADA_TNF   = 0.18    // µg/mL
  Emax_IXE_IL17  = 0.96
  EC50_IXE_IL17  = 0.12    // µg/mL
  Emax_GUS_IL23  = 0.94
  EC50_GUS_IL23  = 0.08    // µg/mL
  Emax_UPA_JAK   = 0.88    // JAK1/TYK2 inhibition
  EC50_UPA_JAK   = 0.082   // µg/mL  IC50 JAK1 ~45nM
  Emax_APR_PDE4  = 0.78
  EC50_APR_PDE4  = 0.54    // µg/mL

  // Hill coefficients
  n_hill = 1.5

  // Feedback strengths (cytokine cross-regulation)
  fb_IL17_TNF  = 0.35  // IL-17 amplifies TNF
  fb_TNF_IL17  = 0.28
  fb_IL23_Th17 = 0.55  // IL-23 drives Th17
  fb_Th17_IL17 = 0.70

$CMT
  // ADA PK
  ADA_DEPOT ADA_C1 ADA_C2 ADA_RC
  // IXE PK
  IXE_DEPOT IXE_C1 IXE_RC
  // GUS PK
  GUS_DEPOT GUS_C1 GUS_RC
  // UPA PK
  UPA_GI UPA_C1
  // APR PK
  APR_GI APR_C1
  // Disease PD
  IL17 TNFa IL23 TH17 RANKL CRP_ratio
  S100 PASI DAPSA MTSS

$MAIN
  double ADA_conc  = ADA_C1  / V1_ADA;   // µg/mL (=mg/L, MW ~148kDa)
  double IXE_conc  = IXE_C1  / V1_IXE;
  double GUS_conc  = GUS_C1  / V1_GUS;
  double UPA_conc  = UPA_C1  / V1_UPA * 1000.0;  // µg/mL from mg/L
  double APR_conc  = APR_C1  / V1_APR * 1000.0;

  // Emax drug effects (inhibitory)
  double Imax_ADA  = Emax_ADA_TNF  * pow(ADA_conc,  n_hill) /
                     (pow(EC50_ADA_TNF,  n_hill) + pow(ADA_conc,  n_hill));
  double Imax_IXE  = Emax_IXE_IL17 * pow(IXE_conc,  n_hill) /
                     (pow(EC50_IXE_IL17, n_hill) + pow(IXE_conc,  n_hill));
  double Imax_GUS  = Emax_GUS_IL23 * pow(GUS_conc,  n_hill) /
                     (pow(EC50_GUS_IL23, n_hill) + pow(GUS_conc,  n_hill));
  double Imax_UPA  = Emax_UPA_JAK  * pow(UPA_conc,  n_hill) /
                     (pow(EC50_UPA_JAK,  n_hill) + pow(UPA_conc,  n_hill));
  double Imax_APR  = Emax_APR_PDE4 * pow(APR_conc,  n_hill) /
                     (pow(EC50_APR_PDE4, n_hill) + pow(APR_conc,  n_hill));

  // Combined JAKi effect on multiple cytokines
  double Imax_JAK_IL17 = Imax_UPA * 0.75;
  double Imax_JAK_TNF  = Imax_UPA * 0.60;
  double Imax_JAK_IL23 = Imax_UPA * 0.80;

  // Combined PDE4i effect
  double Imax_PDE4_TNF  = Imax_APR * 0.55;
  double Imax_PDE4_IL17 = Imax_APR * 0.45;

  // Initialization of baseline at T=0
  if(NEWIND <= 1) {
    IL17_0  = IL17_base;
    TNFa_0  = TNF_base;
    IL23_0  = IL23_base;
    TH17_0  = Th17_base;
    RANKL_0 = RANKL_base;
    CRP_ratio_0 = CRP_base;
    S100_0  = S100_base;
    PASI_0  = PASI_base;
    DAPSA_0 = DAPSA_base;
    MTSS_0  = MTSS_base;
  }

$INIT
  ADA_DEPOT = 0, ADA_C1 = 0, ADA_C2 = 0, ADA_RC = 0
  IXE_DEPOT = 0, IXE_C1 = 0, IXE_RC = 0
  GUS_DEPOT = 0, GUS_C1 = 0, GUS_RC = 0
  UPA_GI    = 0, UPA_C1 = 0
  APR_GI    = 0, APR_C1 = 0
  IL17  = 1.286, TNFa  = 1.294, IL23  = 1.25
  TH17  = 1.0,   RANKL = 1.25,  CRP_ratio = 1.2
  S100  = 1.3,   PASI  = 18.0,  DAPSA = 28.0,  MTSS  = 0.0

$ODE
  // ============================================================
  // ADA PK
  // ============================================================
  double ka_eff_ADA = ka_ADA * F_ADA;
  dxdt_ADA_DEPOT = -ka_eff_ADA * ADA_DEPOT;
  dxdt_ADA_C1    =  ka_eff_ADA * ADA_DEPOT
                   - (CL_ADA / V1_ADA) * ADA_C1
                   - (Q_ADA  / V1_ADA) * ADA_C1
                   + (Q_ADA  / V2_ADA) * ADA_C2
                   - kon_ADA * ADA_conc * TNFa + koff_ADA * ADA_RC;
  dxdt_ADA_C2    =  (Q_ADA  / V1_ADA) * ADA_C1
                   - (Q_ADA  / V2_ADA) * ADA_C2;
  dxdt_ADA_RC    =  kon_ADA * ADA_conc * TNFa
                   - koff_ADA * ADA_RC
                   - kdeg_ADA * ADA_RC;

  // ============================================================
  // IXE PK
  // ============================================================
  double ka_eff_IXE = ka_IXE * F_IXE;
  dxdt_IXE_DEPOT = -ka_eff_IXE * IXE_DEPOT;
  dxdt_IXE_C1    =  ka_eff_IXE * IXE_DEPOT
                   - (CL_IXE / V1_IXE) * IXE_C1
                   - kon_IXE * IXE_conc * IL17 + koff_IXE * IXE_RC;
  dxdt_IXE_RC    =  kon_IXE * IXE_conc * IL17
                   - koff_IXE * IXE_RC
                   - kdeg_IXE * IXE_RC;

  // ============================================================
  // GUS PK
  // ============================================================
  double ka_eff_GUS = ka_GUS * F_GUS;
  dxdt_GUS_DEPOT = -ka_eff_GUS * GUS_DEPOT;
  dxdt_GUS_C1    =  ka_eff_GUS * GUS_DEPOT
                   - (CL_GUS / V1_GUS) * GUS_C1
                   - kon_GUS * GUS_conc * IL23 + koff_GUS * GUS_RC;
  dxdt_GUS_RC    =  kon_GUS * GUS_conc * IL23
                   - koff_GUS * GUS_RC
                   - kdeg_GUS * GUS_RC;

  // ============================================================
  // UPA PK (oral, 1-cmt)
  // ============================================================
  double ka_eff_UPA = ka_UPA * F_UPA;
  dxdt_UPA_GI =  -ka_eff_UPA * UPA_GI;
  dxdt_UPA_C1 =   ka_eff_UPA * UPA_GI - (CL_UPA / V1_UPA) * UPA_C1;

  // ============================================================
  // APR PK (oral, 1-cmt)
  // ============================================================
  double ka_eff_APR = ka_APR * F_APR;
  dxdt_APR_GI =  -ka_eff_APR * APR_GI;
  dxdt_APR_C1 =   ka_eff_APR * APR_GI - (CL_APR / V1_APR) * APR_C1;

  // ============================================================
  // Disease PD ODEs
  // ============================================================

  // Cytokine cross-talk multipliers
  double TNF_drive  = 1.0 + fb_TNF_IL17 * (TNFa / TNF_base  - 1.0);
  double IL23_drive = 1.0 + fb_IL23_Th17* (IL23  / IL23_base - 1.0);
  double Th17_drive = 1.0 + fb_Th17_IL17* (TH17  / Th17_base - 1.0);

  // Total inhibition of IL-17A
  double Itot_IL17 = 1.0 - (1.0 - (1.0 - Imax_IXE))  *
                            (1.0 - (Imax_JAK_IL17))    *
                            (1.0 - (Imax_PDE4_IL17));

  // IL-17A: produced by Th17/ILC3, IL-23 driven, inhibited by IXE/JAKi/PDE4i
  dxdt_IL17 = kin_IL17 * TNF_drive * IL23_drive * Th17_drive * (1.0 - Itot_IL17)
             - kout_IL17 * IL17;

  // TNF-α: produced by macrophages/Th1, inhibited by ADA/JAKi/PDE4i
  double Itot_TNF = 1.0 - (1.0 - Imax_ADA) *
                           (1.0 - Imax_JAK_TNF) *
                           (1.0 - Imax_PDE4_TNF);
  double IL17_drive_TNF = 1.0 + fb_IL17_TNF * (IL17 / IL17_base - 1.0);
  dxdt_TNFa = kin_TNF * IL17_drive_TNF * (1.0 - Itot_TNF)
             - kout_TNF * TNFa;

  // IL-23: produced by DCs/macrophages, inhibited by GUS/JAKi
  double Itot_IL23 = 1.0 - (1.0 - Imax_GUS) * (1.0 - Imax_JAK_IL23);
  dxdt_IL23 = kin_IL23 * (1.0 - Itot_IL23)
             - kout_IL23 * IL23;

  // Th17 cells: driven by IL-23, suppressed by drug-mediated IL reduction
  double Itot_Th17 = Itot_IL23 * 0.6 + Itot_IL17 * 0.3;
  dxdt_TH17 = kin_Th17 * IL23_drive * (1.0 - Itot_Th17)
             - kout_Th17 * TH17;

  // RANKL (bone resorption signal)
  double IL17_drive_RANKL = IL17 / IL17_base;
  double TNF_drive_RANKL  = TNFa / TNF_base;
  double Drug_RANKL_inh   = 0.5 * Itot_IL17 + 0.4 * Itot_TNF;
  dxdt_RANKL = kin_RANKL * IL17_drive_RANKL * TNF_drive_RANKL * (1.0 - Drug_RANKL_inh)
              - kout_RANKL * RANKL;

  // CRP (log-linear surrogate of inflammation)
  double Infl_drive_CRP = (IL17 / IL17_base + TNFa / TNF_base + IL23 / IL23_base) / 3.0;
  double Drug_CRP_inh   = 0.4 * Itot_IL17 + 0.4 * Itot_TNF + 0.2 * Itot_IL23;
  dxdt_CRP_ratio = kin_CRP * Infl_drive_CRP * (1.0 - Drug_CRP_inh)
                  - kout_CRP * CRP_ratio;

  // Calprotectin (S100A8/9 – neutrophil/macrophage marker)
  double Drug_S100_inh = 0.45 * Itot_IL17 + 0.35 * Itot_TNF + 0.20 * Itot_IL23;
  dxdt_S100 = kin_S100 * (IL17 / IL17_base) * (TNFa / TNF_base) * (1.0 - Drug_S100_inh)
             - kout_S100 * S100;

  // PASI (skin score – driven by IL-17A primarily)
  double Drug_PASI_inh = 0.65 * Itot_IL17 + 0.20 * Itot_IL23 + 0.10 * Itot_TNF + 0.05 * Itot_Th17;
  dxdt_PASI = kin_PASI * (IL17 / IL17_base) * (TNFa / TNF_base) * (1.0 - Drug_PASI_inh)
             - kout_PASI * PASI;

  // DAPSA (joint disease activity – IL-17, TNF, IL-23 driven)
  double Drug_DAPSA_inh = 0.45 * Itot_TNF + 0.35 * Itot_IL17 + 0.15 * Itot_IL23 + 0.05 * Itot_Th17;
  dxdt_DAPSA = kin_DAPSA * (TNFa / TNF_base) * (IL17 / IL17_base) * (1.0 - Drug_DAPSA_inh)
              - kout_DAPSA * DAPSA;

  // mTSS (bone/structural damage – slow accumulation, RANKL driven)
  dxdt_MTSS = kin_MTSS * (RANKL / RANKL_base) * (1.0 - 0.7 * Drug_RANKL_inh);

$TABLE
  double CONC_ADA  = ADA_C1  / V1_ADA;
  double CONC_IXE  = IXE_C1  / V1_IXE;
  double CONC_GUS  = GUS_C1  / V1_GUS;
  double CONC_UPA  = (UPA_C1 / V1_UPA) * 1000.0;
  double CONC_APR  = (APR_C1 / V1_APR) * 1000.0;

  // CRP in mg/L  (baseline assumed 20 mg/L)
  double CRP_mgL   = CRP_ratio * 20.0;

  // PASI75/90/100 response flag (binary)
  double PASI_pct_change = (PASI_base > 0) ? (PASI_base - PASI) / PASI_base * 100.0 : 0.0;
  int PASI75 = (PASI_pct_change >= 75.0) ? 1 : 0;
  int PASI90 = (PASI_pct_change >= 90.0) ? 1 : 0;
  int PASI100= (PASI_pct_change >= 99.0) ? 1 : 0;

  // DAPSA response categories
  int DAPSA_REM = (DAPSA <= 4.0)  ? 1 : 0;
  int DAPSA_LDA = (DAPSA <= 14.0) ? 1 : 0;
  int DAPSA_MDA_flag = (DAPSA <= 14.0 && PASI_pct_change >= 75.0) ? 1 : 0;

  // Estimated ACR20 (logistic proxy based on DAPSA change)
  double DAPSA_pct_change = (DAPSA_base > 0) ? (DAPSA_base - DAPSA) / DAPSA_base * 100.0 : 0.0;
  int ACR20 = (DAPSA_pct_change >= 20.0) ? 1 : 0;
  int ACR50 = (DAPSA_pct_change >= 50.0) ? 1 : 0;
  int ACR70 = (DAPSA_pct_change >= 70.0) ? 1 : 0;

  // Calprotectin in µg/mL (baseline ~3.5 µg/mL)
  double Calprotectin = S100 * 3.5;

  capture CONC_ADA CONC_IXE CONC_GUS CONC_UPA CONC_APR
  capture IL17 TNFa IL23 TH17 RANKL CRP_mgL Calprotectin
  capture PASI DAPSA MTSS
  capture PASI_pct_change DAPSA_pct_change
  capture PASI75 PASI90 PASI100
  capture ACR20 ACR50 ACR70 DAPSA_REM DAPSA_LDA DAPSA_MDA_flag
'

## --------------------------------------------------------
## Compile model
## --------------------------------------------------------
mod <- mread("psa", tempdir(), psa_code)

cat("Model compiled successfully\n")
cat("Compartments:", length(init(mod)), "\n")

## --------------------------------------------------------
## Helper: create dosing event records
## --------------------------------------------------------
make_dosing <- function(scenario, weeks = 52) {
  hours <- weeks * 7 * 24
  ev_list <- list()

  if ("ADA" %in% scenario) {
    # Adalimumab 40mg SC Q2W
    times_ADA <- seq(0, hours, by = 14*24)
    ev_list$ADA <- ev(cmt = 1, amt = 40, time = times_ADA)  # ADA_DEPOT
  }
  if ("IXE" %in% scenario) {
    # Ixekizumab: 160mg wk0, 80mg Q2W wk2-16, Q4W thereafter
    times_induction <- seq(0, 16*7*24, by = 2*7*24)
    times_maint     <- seq(20*7*24, hours, by = 4*7*24)
    amt_IXE <- c(160, rep(80, length(times_induction)-1))
    ev_IXE  <- rbind(
      ev(cmt = 5, amt = amt_IXE, time = times_induction),
      ev(cmt = 5, amt = 80,      time = times_maint)
    )
    ev_list$IXE <- ev_IXE
  }
  if ("GUS" %in% scenario) {
    # Guselkumab: 100mg SC wk0, wk4, then Q8W
    times_GUS <- c(0, 4*7*24, seq(12*7*24, hours, by = 8*7*24))
    ev_list$GUS <- ev(cmt = 8, amt = 100, time = times_GUS)  # GUS_DEPOT
  }
  if ("UPA" %in% scenario) {
    # Upadacitinib 15mg QD
    times_UPA <- seq(0, hours, by = 24)
    ev_list$UPA <- ev(cmt = 11, amt = 15, time = times_UPA)   # UPA_GI
  }
  if ("APR" %in% scenario) {
    # Apremilast 30mg BID
    times_APR <- seq(0, hours, by = 12)
    ev_list$APR <- ev(cmt = 13, amt = 30, time = times_APR)   # APR_GI
  }
  if (length(ev_list) == 0) return(NULL)
  do.call(c, ev_list)
}

## --------------------------------------------------------
## Scenario definitions (5 treatment + 1 placebo)
## --------------------------------------------------------
scenarios <- list(
  "Placebo"        = c(),
  "Adalimumab"     = c("ADA"),
  "Ixekizumab"     = c("IXE"),
  "Guselkumab"     = c("GUS"),
  "Upadacitinib"   = c("UPA"),
  "Apremilast"     = c("APR")
)

sim_weeks  <- 52
sim_hours  <- sim_weeks * 7 * 24
obs_times  <- seq(0, sim_hours, by = 24)   # daily output

## --------------------------------------------------------
## Run simulations
## --------------------------------------------------------
cat("Running", length(scenarios), "scenarios ...\n")

results <- lapply(names(scenarios), function(sc) {
  evnt <- if (length(scenarios[[sc]]) == 0) {
    ev(cmt = 1, amt = 0, time = 0)  # null event for placebo
  } else {
    make_dosing(scenarios[[sc]], weeks = sim_weeks)
  }
  if (is.null(evnt)) {
    evnt <- ev(cmt = 1, amt = 0, time = 0)
  }

  out <- mod %>%
    ev(evnt) %>%
    mrgsim(end = sim_hours, delta = 24, obsonly = TRUE) %>%
    as.data.frame()

  out$Scenario <- sc
  out$Week     <- out$time / (7 * 24)
  out
})

df_all <- bind_rows(results)
cat("Simulation complete. Rows:", nrow(df_all), "\n")

## --------------------------------------------------------
## Summary table: key endpoints at weeks 12, 24, 52
## --------------------------------------------------------
summary_table <- df_all %>%
  filter(Week %in% c(12, 24, 52)) %>%
  group_by(Scenario, Week) %>%
  summarise(
    PASI_mean        = round(mean(PASI), 1),
    PASI75_pct       = round(mean(PASI75) * 100, 1),
    PASI90_pct       = round(mean(PASI90) * 100, 1),
    DAPSA_mean       = round(mean(DAPSA), 1),
    ACR20_pct        = round(mean(ACR20) * 100, 1),
    ACR50_pct        = round(mean(ACR50) * 100, 1),
    DAPSA_REM_pct    = round(mean(DAPSA_REM) * 100, 1),
    CRP_mgL_mean     = round(mean(CRP_mgL), 1),
    IL17_mean        = round(mean(IL17), 3),
    MTSS_mean        = round(mean(MTSS), 4),
    .groups = "drop"
  )

cat("\n=== Key Endpoints Summary ===\n")
print(as.data.frame(summary_table))

## --------------------------------------------------------
## Plot 1: PASI over time
## --------------------------------------------------------
p1 <- ggplot(df_all, aes(x = Week, y = PASI, color = Scenario)) +
  geom_line(size = 1.1) +
  geom_hline(yintercept = c(0.25*18, 0.1*18), linetype = "dashed", alpha = 0.5) +
  annotate("text", x = 50, y = 0.25*18 + 0.5, label = "PASI75", size = 3) +
  annotate("text", x = 50, y = 0.1*18  + 0.5, label = "PASI90", size = 3) +
  scale_color_brewer(palette = "Set1") +
  labs(title = "PsA QSP — PASI over 52 Weeks",
       x = "Week", y = "PASI Score",
       color = "Treatment") +
  theme_bw(base_size = 12)

## --------------------------------------------------------
## Plot 2: DAPSA (joint) over time
## --------------------------------------------------------
p2 <- ggplot(df_all, aes(x = Week, y = DAPSA, color = Scenario)) +
  geom_line(size = 1.1) +
  geom_hline(yintercept = c(4, 14), linetype = "dashed", alpha = 0.5) +
  annotate("text", x = 50, y = 4.5,  label = "Remission (≤4)", size = 3) +
  annotate("text", x = 50, y = 14.5, label = "LDA (≤14)",      size = 3) +
  scale_color_brewer(palette = "Set1") +
  labs(title = "PsA QSP — DAPSA over 52 Weeks",
       x = "Week", y = "DAPSA Score",
       color = "Treatment") +
  theme_bw(base_size = 12)

## --------------------------------------------------------
## Plot 3: IL-17A serum concentration
## --------------------------------------------------------
p3 <- ggplot(df_all, aes(x = Week, y = IL17, color = Scenario)) +
  geom_line(size = 1.1) +
  scale_color_brewer(palette = "Set1") +
  labs(title = "PsA QSP — IL-17A Dynamics",
       x = "Week", y = "IL-17A (relative nM)",
       color = "Treatment") +
  theme_bw(base_size = 12)

## --------------------------------------------------------
## Plot 4: CRP over time
## --------------------------------------------------------
p4 <- ggplot(df_all, aes(x = Week, y = CRP_mgL, color = Scenario)) +
  geom_line(size = 1.1) +
  geom_hline(yintercept = 5, linetype = "dashed", color = "gray40") +
  annotate("text", x = 48, y = 6, label = "Normal (<5 mg/L)", size = 3) +
  scale_color_brewer(palette = "Set1") +
  labs(title = "PsA QSP — CRP over 52 Weeks",
       x = "Week", y = "CRP (mg/L)",
       color = "Treatment") +
  theme_bw(base_size = 12)

## --------------------------------------------------------
## Plot 5: Bone erosion (mTSS) accumulation
## --------------------------------------------------------
p5 <- ggplot(df_all, aes(x = Week, y = MTSS, color = Scenario)) +
  geom_line(size = 1.1) +
  scale_color_brewer(palette = "Set1") +
  labs(title = "PsA QSP — Structural Damage (mTSS proxy) over 52 Weeks",
       x = "Week", y = "mTSS Proxy (relative units)",
       color = "Treatment") +
  theme_bw(base_size = 12)

## --------------------------------------------------------
## Plot 6: ACR20/50/70 bar chart at week 24
## --------------------------------------------------------
df_w24 <- df_all %>%
  filter(abs(Week - 24) < 0.5) %>%
  group_by(Scenario) %>%
  summarise(
    ACR20 = mean(ACR20) * 100,
    ACR50 = mean(ACR50) * 100,
    ACR70 = mean(ACR70) * 100,
    .groups = "drop"
  ) %>%
  pivot_longer(cols = c(ACR20, ACR50, ACR70), names_to = "Response", values_to = "Pct")

p6 <- ggplot(df_w24, aes(x = Scenario, y = Pct, fill = Response)) +
  geom_col(position = "dodge") +
  scale_fill_manual(values = c(ACR20 = "#4CAF50", ACR50 = "#2196F3", ACR70 = "#FF9800")) +
  labs(title = "PsA QSP — ACR Response Rates at Week 24",
       x = "Treatment", y = "Response Rate (%)",
       fill = "Criterion") +
  theme_bw(base_size = 12) +
  theme(axis.text.x = element_text(angle = 30, hjust = 1))

## --------------------------------------------------------
## Plot 7: Drug PK concentrations (biologics)
## --------------------------------------------------------
pk_biologic <- df_all %>%
  filter(Scenario %in% c("Adalimumab", "Ixekizumab", "Guselkumab")) %>%
  select(Week, Scenario, CONC_ADA, CONC_IXE, CONC_GUS) %>%
  pivot_longer(cols = c(CONC_ADA, CONC_IXE, CONC_GUS),
               names_to = "Drug", values_to = "Conc") %>%
  mutate(Conc_use = case_when(
    Scenario == "Adalimumab"   & Drug == "CONC_ADA" ~ Conc,
    Scenario == "Ixekizumab"   & Drug == "CONC_IXE" ~ Conc,
    Scenario == "Guselkumab"   & Drug == "CONC_GUS" ~ Conc,
    TRUE ~ NA_real_
  )) %>%
  filter(!is.na(Conc_use))

p7 <- ggplot(pk_biologic, aes(x = Week, y = Conc_use, color = Scenario)) +
  geom_line(size = 1.1) +
  scale_color_brewer(palette = "Set2") +
  labs(title = "PsA QSP — Biologic Drug Concentrations (µg/mL)",
       x = "Week", y = "Drug Concentration (µg/mL)",
       color = "Drug") +
  theme_bw(base_size = 12)

## --------------------------------------------------------
## Print plots
## --------------------------------------------------------
cat("Generating plots ...\n")
print(p1); print(p2); print(p3); print(p4)
print(p5); print(p6); print(p7)

cat("\n=== Simulation complete ===\n")
cat("Total observations:", nrow(df_all), "\n")
cat("Scenarios:", paste(names(scenarios), collapse = ", "), "\n")
