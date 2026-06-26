# ============================================================
# Sarcoidosis QSP Simulator — Shiny App
# File: sarc_shiny_app.R
# Description: Interactive Shiny dashboard for sarcoidosis
#   QSP simulation using mrgsolve ODE model
# ============================================================

library(shiny)
library(mrgsolve)
library(ggplot2)
library(plotly)
library(DT)
library(dplyr)
library(tidyr)

# ============================================================
# mrgsolve Model Definition (inline via mcode)
# ============================================================
sarc_model_code <- '
$PROB Sarcoidosis QSP Model
  Mechanistic PKPD model for sarcoidosis pharmacotherapy
  Compartments: Pred PK, MTX PK/PD, Granuloma, Cytokines,
                Biomarkers (ACE, sIL-2R, calcitriol, Ca), FVC

$PARAM
  // ---- Drug doses (mg/d or mg/wk) ----
  PRED_DOSE  = 0     // Prednisone daily dose (mg/d)
  MTX_DOSE   = 0     // MTX weekly dose (mg/wk)
  AZA_DOSE   = 0     // Azathioprine daily dose (mg/d)
  IFX_ON     = 0     // Infliximab flag (0/1)
  IFX_DOSE   = 0     // Infliximab dose (mg/kg)
  ADA_ON     = 0     // Adalimumab flag (0/1)
  HCQ_DOSE   = 0     // HCQ daily dose (mg/d)
  BW         = 70    // Body weight (kg)
  TAPER_WK   = 999   // Week to start prednisone taper
  TAPER_RATE = 0.85  // Fraction per 4 wks

  // ---- Prednisolone PK ----
  KA_PRED  = 1.5     // Absorption rate (1/h)
  CL_PRED  = 14.0    // Clearance (L/h)
  V1_PRED  = 42.0    // Central volume (L)
  V2_PRED  = 90.0    // Peripheral volume (L)
  Q_PRED   = 6.5     // Inter-compartmental CL (L/h)
  F_PRED   = 0.82    // Bioavailability

  // ---- MTX PK ----
  KA_MTX   = 0.80    // Absorption rate (1/h)
  CL_MTX   = 5.2     // Clearance (L/h)
  V_MTX    = 55.0    // Volume of distribution (L)
  KPG      = 0.012   // MTX → polyglutamate conversion (1/h)
  KOUT_PG  = 0.0008  // Polyglutamate elimination (1/h)

  // ---- Macrophage & Immune Compartments ----
  KIN_MAC  = 0.005   // Macrophage influx base (1/h)
  KOUT_MAC = 0.002   // Macrophage efflux (1/h)
  MAC0     = 2.5     // Macrophage baseline

  KIN_TH1  = 0.004
  KOUT_TH1 = 0.003
  TH10     = 1.5

  KIN_TREG = 0.003
  KOUT_TREG= 0.002
  TREG0    = 0.8

  // ---- Cytokines ----
  KPROD_TNF = 0.08
  KDEG_TNF  = 0.50
  TNF0      = 1.0

  KPROD_IFNG = 0.06
  KDEG_IFNG  = 0.45
  IFNG0      = 1.0

  KPROD_IL12 = 0.05
  KDEG_IL12  = 0.40
  IL120      = 1.0

  // ---- Granuloma dynamics ----
  KFORM_GRAN = 0.010   // Granuloma formation rate
  KRES_GRAN  = 0.001   // Spontaneous resolution rate
  GRAN0      = 1.0     // Baseline granuloma burden (normalized)

  // ---- Fibrosis ----
  KFIB_ON   = 0.0005
  KFIB_OFF  = 0.0002
  FIBR0     = 0.05

  // ---- Biomarker parameters ----
  // ACE: produced by macrophages in granuloma
  KPROD_ACE = 0.60
  KDEG_ACE  = 0.025
  ACE0_NORM = 65.0    // Baseline serum ACE (U/L)

  // Calcitriol (1,25-VitD3): produced by granuloma macrophages
  KPROD_CALIT = 0.04
  KDEG_CALIT  = 0.08
  CALIT0      = 35.0   // pg/mL normal

  // Serum calcium: driven by calcitriol
  KPROD_CA = 0.030
  KDEG_CA  = 0.012
  CA0      = 9.5       // mg/dL normal

  // sIL-2R: produced by activated T cells
  KPROD_SIL2R = 0.12
  KDEG_SIL2R  = 0.050
  SIL2R0      = 900.0  // U/mL normal upper ~1000

  // FVC %pred: decreases with fibrosis
  FVC_BASE = 80.0
  FVC_FIBR = 15.0   // max FVC loss from max fibrosis
  KFVC_REC = 0.0003

  // ---- Drug PD effect parameters ----
  // Prednisone: inhibits macrophage activation, TNF, IL-12
  IC50_PRED_MAC  = 0.08   // µg/mL
  IC50_PRED_TNF  = 0.05
  EC50_PRED_TREG = 0.12   // Treg induction

  // MTX: inhibits proliferation via polyglutamate
  IC50_PG_TH1  = 0.20   // µM poly-glutamate
  IC50_PG_TNF  = 0.15

  // AZA: inhibits T-cell proliferation (simplified as effect on TH1)
  IC50_AZA_TH1 = 0.30   // µg/mL equivalent

  // Infliximab: TNF neutralization
  IC50_IFX_TNF = 1.5   // µg/mL effective (simplified)

  // Adalimumab: TNF neutralization
  IC50_ADA_TNF = 1.2

  // HCQ: reduces macrophage activation
  IC50_HCQ_MAC = 0.80   // µg/mL

$INIT
  PRED_GUT  = 0
  PRED_C    = 0
  PREDL_C   = 0
  PREDL_P   = 0
  MTX_GUT   = 0
  MTX_C     = 0
  MTX_POLY  = 0
  MAC_ACT   = 2.5
  TH1       = 1.5
  TREG      = 0.8
  TNF_sim   = 1.0
  IFNG_sim  = 1.0
  IL12_sim  = 1.0
  GRAN      = 1.0
  FIBR      = 0.05
  ACE_BM    = 65.0
  CALIT     = 35.0
  SERUM_CA  = 9.5
  SIL2R     = 900.0
  FVC_P     = 80.0

$ODE
  // ---- Dosing variables ----
  // Prednisone: convert mg/d to ug/mL input into gut
  // Doses are handled via event records; here we just model PK
  double PRED_CONC = PRED_C / V1_PRED;       // ug/mL (approx)
  double MTX_CONC  = MTX_C / V_MTX;          // ug/mL
  double PG_CONC   = MTX_POLY;               // nmol/L (simplified)
  double AZA_CONC  = AZA_DOSE / 100.0;       // normalized 0-2
  double HCQ_CONC  = HCQ_DOSE / 400.0;       // normalized 0-1

  // Infliximab / Adalimumab (simplified: flat effect during dosing)
  double IFX_EFF   = IFX_ON * (IFX_DOSE * BW / 70.0) / 5.0;   // 0-1 scale
  double ADA_EFF   = ADA_ON * 0.70;

  // ---- Inhibitory Hill functions ----
  double INH_MAC_PRED = (PRED_CONC > 0) ? PRED_CONC / (PRED_CONC + IC50_PRED_MAC) : 0.0;
  double INH_MAC_HCQ  = HCQ_CONC / (HCQ_CONC + IC50_HCQ_MAC);
  double INH_MAC_TOTAL= 1.0 - (1.0 - INH_MAC_PRED) * (1.0 - INH_MAC_HCQ);

  double INH_TH1_PG   = (PG_CONC > 0) ? PG_CONC / (PG_CONC + IC50_PG_TH1) : 0.0;
  double INH_TH1_AZA  = AZA_CONC / (AZA_CONC + IC50_AZA_TH1);
  double INH_TH1_TOTAL= 1.0 - (1.0 - INH_TH1_PG) * (1.0 - INH_TH1_AZA);

  double INH_TNF_PRED = (PRED_CONC > 0) ? PRED_CONC / (PRED_CONC + IC50_PRED_TNF) : 0.0;
  double INH_TNF_PG   = (PG_CONC > 0) ? PG_CONC / (PG_CONC + IC50_PG_TNF) : 0.0;
  double INH_TNF_IFX  = IFX_EFF / (IFX_EFF + IC50_IFX_TNF / 5.0 + 1e-9);
  double INH_TNF_ADA  = ADA_EFF * 0.70;
  double INH_TNF_TOTAL= 1.0 - (1.0 - INH_TNF_PRED) * (1.0 - INH_TNF_PG) *
                               (1.0 - INH_TNF_IFX) * (1.0 - INH_TNF_ADA);

  double STIM_TREG_PRED = (PRED_CONC > 0) ? PRED_CONC / (PRED_CONC + EC50_PRED_TREG) : 0.0;

  // Taper: reduce PRED_DOSE linearly (handled externally via events)
  // Here PRED_DOSE is set externally per time step

  // ---- Prednisolone PK (2-compartment oral) ----
  dxdt_PRED_GUT = -KA_PRED * PRED_GUT;
  dxdt_PRED_C   =  KA_PRED * PRED_GUT * F_PRED
                   - (CL_PRED + Q_PRED) / V1_PRED * PRED_C
                   + Q_PRED / V2_PRED * PREDL_P;
  dxdt_PREDL_C  =  Q_PRED / V1_PRED * PRED_C - Q_PRED / V2_PRED * PREDL_C;
  dxdt_PREDL_P  =  Q_PRED / V1_PRED * PRED_C - Q_PRED / V2_PRED * PREDL_P;

  // ---- MTX PK + polyglutamate ----
  dxdt_MTX_GUT  = -KA_MTX * MTX_GUT;
  dxdt_MTX_C    =  KA_MTX * MTX_GUT - CL_MTX / V_MTX * MTX_C - KPG * MTX_C;
  dxdt_MTX_POLY =  KPG * MTX_C - KOUT_PG * MTX_POLY;

  // ---- Macrophage activation ----
  // Positive feedback: TNF activates macrophages; Pred, HCQ inhibit
  double MAC_STIM = TNF_sim * IFNG_sim;
  dxdt_MAC_ACT = KIN_MAC * MAC_STIM * (1.0 - INH_MAC_TOTAL) +
                 KOUT_MAC * MAC0 - KOUT_MAC * MAC_ACT;

  // ---- Th1 cells ----
  dxdt_TH1 = KIN_TH1 * IL12_sim * (1.0 - INH_TH1_TOTAL) +
              KOUT_TH1 * TH10 - KOUT_TH1 * TH1 -
              0.05 * TREG * TH1;

  // ---- Treg cells ----
  dxdt_TREG = KIN_TREG * (1.0 + STIM_TREG_PRED) +
              KOUT_TREG * TREG0 - KOUT_TREG * TREG;

  // ---- Cytokines ----
  dxdt_TNF_sim  = KPROD_TNF * MAC_ACT * (1.0 - INH_TNF_TOTAL) -
                  KDEG_TNF * TNF_sim;
  dxdt_IFNG_sim = KPROD_IFNG * TH1 - KDEG_IFNG * IFNG_sim;
  dxdt_IL12_sim = KPROD_IL12 * MAC_ACT - KDEG_IL12 * IL12_sim -
                  0.03 * TREG * IL12_sim;

  // ---- Granuloma burden ----
  double GRAN_FORM = KFORM_GRAN * MAC_ACT * TNF_sim * IFNG_sim;
  double GRAN_RES  = (KRES_GRAN + 0.008 * INH_MAC_TOTAL + 0.005 * INH_TNF_TOTAL) * GRAN;
  dxdt_GRAN = GRAN_FORM - GRAN_RES;

  // ---- Fibrosis (cumulative, partially reversible) ----
  double FIBR_DRIVE = (GRAN > 1.0) ? KFIB_ON * (GRAN - 1.0) : 0.0;
  dxdt_FIBR = FIBR_DRIVE - KFIB_OFF * FIBR;

  // ---- Biomarkers ----
  // ACE: produced proportional to granuloma burden
  dxdt_ACE_BM  = KPROD_ACE * GRAN - KDEG_ACE * ACE_BM;

  // Calcitriol: produced by activated macrophages in granuloma
  dxdt_CALIT   = KPROD_CALIT * GRAN * MAC_ACT - KDEG_CALIT * CALIT;

  // Serum calcium: driven by calcitriol
  dxdt_SERUM_CA = KPROD_CA * CALIT / 35.0 - KDEG_CA * SERUM_CA;

  // sIL-2R: produced by TH1 + MAC_ACT
  dxdt_SIL2R   = KPROD_SIL2R * (TH1 + MAC_ACT) - KDEG_SIL2R * SIL2R;

  // FVC %pred: decreases with fibrosis, partially recovers with treatment
  double FVC_TARGET = FVC_BASE - FVC_FIBR * FIBR;
  dxdt_FVC_P = KFVC_REC * (FVC_TARGET - FVC_P);

$TABLE
  double PRED_CP  = PRED_C / V1_PRED * 1000.0;   // ng/mL
  double MTX_CP   = MTX_C  / V_MTX;              // ug/mL
  double GRAN_OUT = GRAN;
  double ACE_OUT  = ACE_BM;
  double CALIT_OUT= CALIT;
  double CA_OUT   = SERUM_CA;
  double SIL2R_OUT= SIL2R;
  double FVC_OUT  = FVC_P;
  double FIBR_OUT = FIBR;
  double DLCO_OUT = FVC_P * (1.0 - 0.20 * FIBR);   // rough estimate
  // mMRC dyspnea score 0-4 driven by FVC decline
  double mMRC_OUT = pmax(0.0, pmin(4.0, (80.0 - FVC_P) / 10.0));

$CAPTURE PRED_CP MTX_CP GRAN_OUT ACE_OUT CALIT_OUT CA_OUT SIL2R_OUT FVC_OUT FIBR_OUT DLCO_OUT mMRC_OUT MTX_POLY
'

# Build the model once at startup
sarc_mod <- mcode("sarcoidosis_qsp", sarc_model_code, quiet = TRUE)

# ============================================================
# Helper: build event table from UI inputs
# ============================================================
build_events <- function(pred_dose_init, mtx_dose, aza_dose,
                         ifx_on, ifx_dose, ada_on, hcq_dose,
                         bw, taper_on, taper_wk_start, sim_weeks) {

  sim_hours <- sim_weeks * 168

  # Prednisone: daily dosing (every 24h)
  ev_list <- list()

  if (pred_dose_init > 0) {
    # Convert mg → ug for PRED_GUT compartment (amt in ug)
    pred_amt <- pred_dose_init * 1000   # ug
    # Add doses for the full simulation; taper handled by reducing amt
    n_days <- sim_weeks * 7

    if (taper_on && taper_wk_start <= sim_weeks) {
      # Pre-taper doses
      taper_start_d <- taper_wk_start * 7
      pred_days_pre  <- seq(0, taper_start_d - 1, by = 1) * 24
      pred_days_post <- seq(taper_start_d, n_days - 1, by = 1) * 24
      # Exponential taper after taper_wk_start
      wks_since_taper <- (seq(taper_start_d, n_days - 1) - taper_start_d) %/% 28
      taper_amts      <- pred_amt * (0.80 ^ wks_since_taper)
      taper_amts      <- pmax(taper_amts, 5000)  # floor 5 mg

      ev_pre <- ev(amt = pred_amt, cmt = "PRED_GUT", time = pred_days_pre)
      ev_post <- ev(amt = taper_amts, cmt = "PRED_GUT", time = pred_days_post)
      ev_list <- c(ev_list, list(ev_pre, ev_post))
    } else {
      pred_times <- seq(0, (n_days - 1)) * 24
      ev_pred    <- ev(amt = pred_amt, cmt = "PRED_GUT", time = pred_times)
      ev_list    <- c(ev_list, list(ev_pred))
    }
  }

  # MTX: weekly dosing (convert mg → ug)
  if (mtx_dose > 0) {
    mtx_times <- seq(0, sim_weeks - 1) * 168
    ev_mtx    <- ev(amt = mtx_dose * 1000, cmt = "MTX_GUT", time = mtx_times)
    ev_list   <- c(ev_list, list(ev_mtx))
  }

  if (length(ev_list) == 0) {
    return(ev(amt = 0, cmt = 1, time = 0))  # dummy
  }

  do.call(c, ev_list)
}

# ============================================================
# Helper: parse drug selections into numeric doses
# ============================================================
parse_drug_inputs <- function(drug1, drug2, biologic, antimalarial, bw) {
  pred_dose <- switch(drug1,
    "Prednisone 40mg/d" = 40,
    "Prednisone 20mg/d" = 20,
    "Prednisone 10mg/d" = 10,
    0)

  mtx_dose <- 0; aza_dose <- 0
  if (grepl("Methotrexate 10mg", drug2))  mtx_dose <- 10
  if (grepl("Methotrexate 15mg", drug2))  mtx_dose <- 15
  if (grepl("Azathioprine 100",  drug2))  aza_dose <- 100
  if (grepl("Azathioprine 150",  drug2))  aza_dose <- 150

  ifx_on <- 0; ifx_dose <- 0; ada_on <- 0
  if (grepl("Infliximab 3",  biologic)) { ifx_on <- 1; ifx_dose <- 3 }
  if (grepl("Infliximab 5",  biologic)) { ifx_on <- 1; ifx_dose <- 5 }
  if (grepl("Adalimumab",    biologic)) { ada_on <- 1 }

  hcq_dose <- 0
  if (grepl("200", antimalarial)) hcq_dose <- 200
  if (grepl("400", antimalarial)) hcq_dose <- 400

  list(pred = pred_dose, mtx = mtx_dose, aza = aza_dose,
       ifx_on = ifx_on, ifx_dose = ifx_dose, ada_on = ada_on,
       hcq = hcq_dose)
}

# ============================================================
# Helper: run simulation
# ============================================================
run_sim <- function(mod, pred_dose, mtx_dose, aza_dose,
                    ifx_on, ifx_dose, ada_on, hcq_dose,
                    bw, taper_on, taper_wk, sim_weeks,
                    ace_baseline = 65, fvc_baseline = 78,
                    disease_stage = "Stage II",
                    disease_activity = "Moderate") {

  # Scale initial conditions by baseline inputs
  ace_scale  <- ace_baseline / 65.0
  fvc_scale  <- fvc_baseline / 80.0
  gran_init  <- switch(disease_stage,
                  "Stage I"   = 1.2,
                  "Stage II"  = 1.8,
                  "Stage III" = 2.5,
                  "Stage IV"  = 3.2,
                  1.8)
  act_mult   <- switch(disease_activity,
                  "Mild"     = 0.8,
                  "Moderate" = 1.0,
                  "Severe"   = 1.5,
                  1.0)
  gran_init <- gran_init * act_mult

  idata <- data.frame(
    PRED_DOSE  = pred_dose,
    MTX_DOSE   = mtx_dose,
    AZA_DOSE   = aza_dose,
    IFX_ON     = ifx_on,
    IFX_DOSE   = ifx_dose,
    ADA_ON     = ada_on,
    HCQ_DOSE   = hcq_dose,
    BW         = bw,
    TAPER_WK   = ifelse(taper_on, taper_wk, 9999),
    GRAN       = gran_init,
    ACE_BM     = ace_baseline,
    FVC_P      = fvc_baseline,
    CALIT      = 35 * gran_init / 1.5,
    SERUM_CA   = 9.5 + 0.3 * (gran_init - 1.0),
    SIL2R      = 900 * gran_init / 1.5,
    MAC_ACT    = 2.5 * act_mult,
    TH1        = 1.5 * act_mult,
    FIBR       = 0.05 * gran_init
  )

  ev_obj <- build_events(pred_dose, mtx_dose, aza_dose,
                         ifx_on, ifx_dose, ada_on, hcq_dose,
                         bw, taper_on, taper_wk, sim_weeks)

  out <- tryCatch({
    mod %>%
      param(PRED_DOSE  = pred_dose,
            MTX_DOSE   = mtx_dose,
            AZA_DOSE   = aza_dose,
            IFX_ON     = ifx_on,
            IFX_DOSE   = ifx_dose,
            ADA_ON     = ada_on,
            HCQ_DOSE   = hcq_dose,
            BW         = bw,
            TAPER_WK   = ifelse(taper_on, taper_wk, 9999)) %>%
      init(GRAN     = gran_init,
           ACE_BM   = ace_baseline,
           FVC_P    = fvc_baseline,
           CALIT    = 35 * gran_init / 1.5,
           SERUM_CA = 9.5 + 0.3 * (gran_init - 1.0),
           SIL2R    = 900 * gran_init / 1.5,
           MAC_ACT  = 2.5 * act_mult,
           TH1      = 1.5 * act_mult,
           FIBR     = 0.05 * gran_init) %>%
      mrgsim(ev_obj,
             end    = sim_weeks * 168,
             delta  = 24,
             obsonly = TRUE) %>%
      as.data.frame()
  }, error = function(e) NULL)

  if (!is.null(out)) {
    out$week <- out$time / 168
  }
  out
}

# ============================================================
# Scenario definitions
# ============================================================
scenarios_def <- list(
  "No treatment"  = list(pred=0,  mtx=0,  aza=0, ifx_on=0, ifx_dose=0, ada_on=0, hcq=0, taper=FALSE, taper_wk=99),
  "Pred 40mg"     = list(pred=40, mtx=0,  aza=0, ifx_on=0, ifx_dose=0, ada_on=0, hcq=0, taper=FALSE, taper_wk=99),
  "Pred+MTX"      = list(pred=20, mtx=15, aza=0, ifx_on=0, ifx_dose=0, ada_on=0, hcq=0, taper=FALSE, taper_wk=99),
  "Infliximab"    = list(pred=10, mtx=0,  aza=0, ifx_on=1, ifx_dose=5, ada_on=0, hcq=0, taper=FALSE, taper_wk=99),
  "HCQ"           = list(pred=0,  mtx=0,  aza=0, ifx_on=0, ifx_dose=0, ada_on=0, hcq=400, taper=FALSE, taper_wk=99),
  "Pred taper"    = list(pred=40, mtx=0,  aza=0, ifx_on=0, ifx_dose=0, ada_on=0, hcq=0, taper=TRUE,  taper_wk=8)
)

# ============================================================
# UI
# ============================================================
ui <- fluidPage(
  tags$head(
    tags$style(HTML("
      body { font-family: 'Helvetica Neue', Arial, sans-serif; background: #f8f9fa; }
      .navbar { background-color: #2c3e50 !important; }
      .app-title { color: #2c3e50; font-weight: bold; font-size: 26px;
                   padding: 15px 0 5px 0; }
      .app-subtitle { color: #7f8c8d; font-size: 14px; margin-bottom: 15px; }
      .tab-content { background: white; padding: 20px;
                     border-radius: 0 0 8px 8px; border: 1px solid #dee2e6;
                     border-top: none; }
      .well { background: #ecf0f1; border: none; }
      .param-box { background: #eaf4fb; border-left: 4px solid #2980b9;
                   padding: 10px; margin-bottom: 10px; border-radius: 4px; }
      .ref-band { color: #888; font-size: 12px; }
    "))
  ),

  div(class = "app-title", "Sarcoidosis QSP Simulator"),
  div(class = "app-subtitle",
      "Quantitative Systems Pharmacology model for sarcoidosis treatment simulation"),

  tabsetPanel(
    id = "main_tabs",

    # ----------------------------------------------------------
    # TAB 1: Patient Profile
    # ----------------------------------------------------------
    tabPanel("Patient Profile",
      fluidRow(
        column(4,
          wellPanel(
            h4("Patient Characteristics"),
            sliderInput("bw", "Body weight (kg)", 50, 120, 70, step = 5),
            radioButtons("disease_stage", "Scadding Stage",
              choices = c("Stage I", "Stage II", "Stage III", "Stage IV"),
              selected = "Stage II", inline = FALSE),
            radioButtons("disease_activity", "Disease Activity",
              choices = c("Mild", "Moderate", "Severe"),
              selected = "Moderate", inline = TRUE),
            hr(),
            h4("Organ Involvement"),
            checkboxGroupInput("organs", NULL,
              choices  = c("Pulmonary", "Cardiac", "CNS", "Hepatic", "Ocular", "Cutaneous"),
              selected = c("Pulmonary"))
          )
        ),
        column(4,
          wellPanel(
            h4("Baseline Biomarkers"),
            numericInput("ace_base", "Baseline serum ACE (U/L)", 65, 10, 300, 1),
            numericInput("fvc_base", "Baseline FVC %predicted", 78, 20, 120, 1),
            hr(),
            h4("Simulation Settings"),
            sliderInput("sim_weeks", "Simulation duration (weeks)",
                        26, 156, 104, step = 4)
          )
        ),
        column(4,
          h4("Patient Summary"),
          DTOutput("patient_summary_tbl")
        )
      )
    ),  # end TAB 1

    # ----------------------------------------------------------
    # TAB 2: Treatment Selection
    # ----------------------------------------------------------
    tabPanel("Treatment Selection",
      fluidRow(
        column(5,
          wellPanel(
            h4("Drug Selection"),
            selectInput("drug1", "Corticosteroid",
              choices  = c("None", "Prednisone 40mg/d", "Prednisone 20mg/d",
                           "Prednisone 10mg/d"),
              selected = "Prednisone 40mg/d"),
            selectInput("drug2", "Steroid-sparing Agent",
              choices  = c("None", "Methotrexate 10mg/wk", "Methotrexate 15mg/wk",
                           "Azathioprine 100mg/d", "Azathioprine 150mg/d"),
              selected = "None"),
            selectInput("biologic", "Biologic Therapy",
              choices  = c("None", "Infliximab 3mg/kg", "Infliximab 5mg/kg",
                           "Adalimumab 40mg q2wk"),
              selected = "None"),
            selectInput("antimalarial", "Antimalarial",
              choices  = c("None", "Hydroxychloroquine 200mg/d",
                           "Hydroxychloroquine 400mg/d"),
              selected = "None"),
            hr(),
            h4("Steroid Taper Strategy"),
            checkboxInput("taper_on", "Start prednisone taper", FALSE),
            sliderInput("taper_wk", "Begin taper at week:", 4, 24, 8, step = 2),
            sliderInput("steroid_dur", "Duration of initial steroid (weeks):", 4, 24, 12, step = 2)
          )
        ),
        column(7,
          h4("Treatment Schedule"),
          DTOutput("tx_schedule_tbl"),
          br(),
          div(class = "param-box",
            h5("Selected Regimen PK Parameters"),
            verbatimTextOutput("pk_params_text")
          )
        )
      )
    ),  # end TAB 2

    # ----------------------------------------------------------
    # TAB 3: PK Profiles
    # ----------------------------------------------------------
    tabPanel("PK Profiles",
      fluidRow(
        column(6,
          h4("Prednisolone Plasma Concentration (first 2 weeks)"),
          plotlyOutput("pk_pred_plot", height = "320px")
        ),
        column(6,
          h4("MTX Plasma & Polyglutamate Concentrations"),
          plotlyOutput("pk_mtx_plot", height = "320px")
        )
      ),
      fluidRow(
        column(12,
          div(class = "param-box",
            h5("PK Parameter Summary"),
            DTOutput("pk_param_tbl")
          )
        )
      )
    ),  # end TAB 3

    # ----------------------------------------------------------
    # TAB 4: Disease PD — Biomarkers
    # ----------------------------------------------------------
    tabPanel("Disease PD — Biomarkers",
      fluidRow(
        column(6,
          h4("Granuloma Burden Over Time"),
          plotlyOutput("gran_plot", height = "300px")
        ),
        column(6,
          h4("Serum ACE Over Time"),
          plotlyOutput("ace_plot", height = "300px")
        )
      ),
      fluidRow(
        column(6,
          h4("Serum sIL-2R Over Time"),
          plotlyOutput("sil2r_plot", height = "300px")
        ),
        column(6,
          h4("Calcitriol (1,25-(OH)₂D₃) Over Time"),
          plotlyOutput("calit_plot", height = "300px")
        )
      ),
      fluidRow(
        column(12,
          h4("Serum Calcium Over Time"),
          plotlyOutput("ca_plot", height = "280px")
        )
      )
    ),  # end TAB 4

    # ----------------------------------------------------------
    # TAB 5: Pulmonary Function
    # ----------------------------------------------------------
    tabPanel("Pulmonary Function",
      fluidRow(
        column(6,
          h4("FVC %predicted Over Time"),
          plotlyOutput("fvc_plot", height = "310px")
        ),
        column(6,
          h4("DLCO %predicted Over Time"),
          plotlyOutput("dlco_plot", height = "310px")
        )
      ),
      fluidRow(
        column(6,
          h4("Dyspnea Score (mMRC) Over Time"),
          plotlyOutput("mmrc_plot", height = "300px")
        ),
        column(6,
          h4("Stage Distribution at Selected Time Points"),
          selectInput("stage_timepoints", "Select time point (week):",
            choices = c(12, 26, 52, 78, 104), selected = 52),
          plotlyOutput("stage_bar_plot", height = "280px")
        )
      )
    ),  # end TAB 5

    # ----------------------------------------------------------
    # TAB 6: Scenario Comparison
    # ----------------------------------------------------------
    tabPanel("Scenario Comparison",
      fluidRow(
        column(4,
          wellPanel(
            h4("Select Scenarios to Compare"),
            checkboxGroupInput("scenarios", NULL,
              choices  = names(scenarios_def),
              selected = c("No treatment", "Pred 40mg", "Pred+MTX"))
          )
        ),
        column(8,
          h4("Granuloma Burden — All Scenarios"),
          plotlyOutput("scenario_gran_plot", height = "280px"),
          h4("Serum ACE — All Scenarios"),
          plotlyOutput("scenario_ace_plot", height = "280px")
        )
      ),
      fluidRow(
        column(12,
          h4("Scenario Comparison Table (Final Week)"),
          DTOutput("scenario_tbl"),
          br(),
          h4("Treatment Response Assessment"),
          DTOutput("response_tbl")
        )
      )
    ),  # end TAB 6

    # ----------------------------------------------------------
    # TAB 7: Biomarker Correlation
    # ----------------------------------------------------------
    tabPanel("Biomarker Correlation",
      fluidRow(
        column(4,
          h4("Granuloma Burden vs Serum ACE"),
          plotlyOutput("corr_gran_ace", height = "280px")
        ),
        column(4,
          h4("Granuloma Burden vs FVC %pred"),
          plotlyOutput("corr_gran_fvc", height = "280px")
        ),
        column(4,
          h4("Calcitriol vs Serum Calcium"),
          plotlyOutput("corr_calit_ca", height = "280px")
        )
      ),
      fluidRow(
        column(12,
          h4("Biomarker Correlation Table"),
          DTOutput("corr_table")
        )
      )
    )   # end TAB 7
  )   # end tabsetPanel
)   # end fluidPage

# ============================================================
# Server
# ============================================================
server <- function(input, output, session) {

  # ---- Reactive: run main simulation ----
  sim_result <- reactive({
    drugs <- parse_drug_inputs(input$drug1, input$drug2,
                               input$biologic, input$antimalarial,
                               input$bw)
    tryCatch(
      run_sim(sarc_mod,
              pred_dose        = drugs$pred,
              mtx_dose         = drugs$mtx,
              aza_dose         = drugs$aza,
              ifx_on           = drugs$ifx_on,
              ifx_dose         = drugs$ifx_dose,
              ada_on           = drugs$ada_on,
              hcq_dose         = drugs$hcq,
              bw               = input$bw,
              taper_on         = input$taper_on,
              taper_wk         = input$taper_wk,
              sim_weeks        = input$sim_weeks,
              ace_baseline     = input$ace_base,
              fvc_baseline     = input$fvc_base,
              disease_stage    = input$disease_stage,
              disease_activity = input$disease_activity),
      error = function(e) { message("Simulation error: ", e$message); NULL }
    )
  })

  # ---- TAB 1: Patient summary table ----
  output$patient_summary_tbl <- renderDT({
    stage_gran <- c("Stage I"=1.2,"Stage II"=1.8,"Stage III"=2.5,"Stage IV"=3.2)
    act_mult   <- c("Mild"=0.8,"Moderate"=1.0,"Severe"=1.5)
    gran_est   <- stage_gran[input$disease_stage] * act_mult[input$disease_activity]

    df <- data.frame(
      Parameter = c("Body Weight", "Scadding Stage", "Disease Activity",
                    "Organ Involvement", "Baseline ACE", "Baseline FVC",
                    "Estimated Granuloma Burden", "Simulation Duration"),
      Value     = c(paste0(input$bw, " kg"),
                    input$disease_stage,
                    input$disease_activity,
                    paste(input$organs, collapse=", "),
                    paste0(input$ace_base, " U/L"),
                    paste0(input$fvc_base, " %pred"),
                    sprintf("%.2f (normalized)", gran_est),
                    paste0(input$sim_weeks, " weeks"))
    )
    datatable(df, rownames=FALSE, options=list(dom='t', paging=FALSE),
              class='compact stripe')
  })

  # ---- TAB 2: Treatment schedule table ----
  output$tx_schedule_tbl <- renderDT({
    drugs <- parse_drug_inputs(input$drug1, input$drug2,
                               input$biologic, input$antimalarial,
                               input$bw)
    df <- data.frame(
      Drug       = c("Prednisone", "Steroid-sparing", "Biologic", "Antimalarial"),
      Regimen    = c(ifelse(drugs$pred>0, paste0(drugs$pred,"mg/d"), "None"),
                     ifelse(drugs$mtx>0, paste0("MTX ",drugs$mtx,"mg/wk"),
                            ifelse(drugs$aza>0, paste0("AZA ",drugs$aza,"mg/d"), "None")),
                     ifelse(drugs$ifx_on==1, paste0("Infliximab ",drugs$ifx_dose,"mg/kg q8wk"),
                            ifelse(drugs$ada_on==1, "Adalimumab 40mg q2wk","None")),
                     ifelse(drugs$hcq>0, paste0("HCQ ",drugs$hcq,"mg/d"), "None")),
      `Taper`    = c(ifelse(input$taper_on, paste0("Start wk ",input$taper_wk),"No"), "", "", ""),
      `Duration` = c(paste0(input$steroid_dur," wks"), paste0(input$sim_weeks," wks"),
                     paste0(input$sim_weeks," wks"), paste0(input$sim_weeks," wks"))
    )
    datatable(df, rownames=FALSE, options=list(dom='t', paging=FALSE),
              class='compact stripe')
  })

  output$pk_params_text <- renderPrint({
    drugs <- parse_drug_inputs(input$drug1, input$drug2,
                               input$biologic, input$antimalarial,
                               input$bw)
    cat("=== Prednisolone PK ===\n")
    cat(sprintf("  Dose: %g mg/d  |  F: 82%%  |  Cmax ~%.0f ng/mL  |  t½ ~3 h\n",
                drugs$pred,
                drugs$pred * 1000 * 0.82 / 42 * 0.85))
    cat("  Vd: 42 L  |  CL: 14 L/h  |  2-compartment model\n\n")
    cat("=== Methotrexate PK ===\n")
    cat(sprintf("  Dose: %g mg/wk  |  F: 70%%  |  Cmax ~%.2f µg/mL  |  t½ ~6 h\n",
                drugs$mtx,
                if (drugs$mtx>0) drugs$mtx*1000*0.7/55/1000 else 0))
    cat("  Vd: 55 L  |  CL: 5.2 L/h  |  Polyglutamate t½ ~weeks\n\n")
    if (drugs$hcq > 0) {
      cat(sprintf("=== HCQ PK ===\n  Dose: %g mg/d  |  t½ ~50 d  |  Vd: ~800 L/kg\n", drugs$hcq))
    }
  })

  # ---- TAB 3: PK Plots ----
  output$pk_pred_plot <- renderPlotly({
    out <- sim_result()
    if (is.null(out)) return(plot_ly() %>% layout(title="Simulation failed"))
    df_pk <- out %>% filter(week <= 2) %>%
      select(week, PRED_CP) %>%
      mutate(day = week * 7)

    plot_ly(df_pk, x = ~day, y = ~PRED_CP, type='scatter', mode='lines',
            line=list(color='#e74c3c', width=2)) %>%
      layout(xaxis=list(title="Day"),
             yaxis=list(title="Prednisolone (ng/mL)"),
             hovermode='x unified',
             margin=list(t=20))
  })

  output$pk_mtx_plot <- renderPlotly({
    out <- sim_result()
    if (is.null(out)) return(plot_ly() %>% layout(title="Simulation failed"))
    df_pk <- out %>% filter(week <= 12) %>%
      mutate(day = week * 7)

    p <- plot_ly() %>%
      add_trace(data=df_pk, x=~day, y=~MTX_CP, name="MTX plasma (µg/mL)",
                type='scatter', mode='lines',
                line=list(color='#3498db', width=2)) %>%
      add_trace(data=df_pk, x=~day, y=~MTX_POLY/1000, name="MTX-PG (scaled)",
                type='scatter', mode='lines',
                line=list(color='#2ecc71', width=2, dash='dash')) %>%
      layout(xaxis=list(title="Day"),
             yaxis=list(title="Concentration"),
             legend=list(orientation='h'),
             hovermode='x unified',
             margin=list(t=20))
    p
  })

  output$pk_param_tbl <- renderDT({
    df <- data.frame(
      Drug=c("Prednisolone","Methotrexate","Azathioprine","Infliximab","HCQ"),
      `Bioavailability`=c("82%","70%","~60%","100% (IV)","~75%"),
      `t½`=c("2-4 h","6-9 h","1-2 h (thiopurine)", "9-12 d","~50 d"),
      `Vd`=c("42 L","55 L","0.8 L/kg","3-6 L","800 L/kg"),
      `CL`=c("14 L/h","5.2 L/h","renal","proteolysis","hepatic"),
      `Mechanism`=c("GR agonist","DHFR/AICAR inhibition",
                    "Purine antimetabolite","Anti-TNF-α mAb",
                    "Lysosome alkalinization"),
      check.names=FALSE
    )
    datatable(df, rownames=FALSE, options=list(dom='t', paging=FALSE),
              class='compact stripe')
  })

  # ---- TAB 4: Biomarker Plots ----
  output$gran_plot <- renderPlotly({
    out <- sim_result()
    if (is.null(out)) return(plot_ly())
    plot_ly(out, x=~week, y=~GRAN_OUT, type='scatter', mode='lines',
            fill='tozeroy', fillcolor='rgba(231,76,60,0.15)',
            line=list(color='#e74c3c', width=2.5)) %>%
      layout(xaxis=list(title="Week"),
             yaxis=list(title="Granuloma Burden (normalized)", rangemode='tozero'),
             shapes=list(list(type='line', x0=0, x1=input$sim_weeks,
                              y0=1.0, y1=1.0,
                              line=list(color='gray', width=1, dash='dot'))),
             annotations=list(list(x=5, y=1.05, text="Normal",
                                   showarrow=FALSE, font=list(color='gray'))),
             margin=list(t=20))
  })

  output$ace_plot <- renderPlotly({
    out <- sim_result()
    if (is.null(out)) return(plot_ly())
    plot_ly(out, x=~week, y=~ACE_OUT, type='scatter', mode='lines',
            line=list(color='#8e44ad', width=2.5)) %>%
      layout(xaxis=list(title="Week"),
             yaxis=list(title="Serum ACE (U/L)"),
             shapes=list(
               list(type='rect', x0=0, x1=input$sim_weeks, y0=18, y1=67,
                    fillcolor='rgba(46,204,113,0.12)',
                    line=list(color='rgba(46,204,113,0.4)', width=1)),
               list(type='line', x0=0, x1=input$sim_weeks, y0=67, y1=67,
                    line=list(color='#e67e22', width=1.5, dash='dash'))
             ),
             annotations=list(list(x=5, y=70, text="ULN 67 U/L",
                                   showarrow=FALSE,
                                   font=list(color='#e67e22', size=11))),
             margin=list(t=20))
  })

  output$sil2r_plot <- renderPlotly({
    out <- sim_result()
    if (is.null(out)) return(plot_ly())
    plot_ly(out, x=~week, y=~SIL2R_OUT, type='scatter', mode='lines',
            line=list(color='#2980b9', width=2.5)) %>%
      layout(xaxis=list(title="Week"),
             yaxis=list(title="sIL-2R (U/mL)"),
             shapes=list(
               list(type='line', x0=0, x1=input$sim_weeks, y0=1000, y1=1000,
                    line=list(color='#e67e22', width=1.5, dash='dash'))
             ),
             annotations=list(list(x=5, y=1040, text="ULN ~1000",
                                   showarrow=FALSE,
                                   font=list(color='#e67e22', size=11))),
             margin=list(t=20))
  })

  output$calit_plot <- renderPlotly({
    out <- sim_result()
    if (is.null(out)) return(plot_ly())
    plot_ly(out, x=~week, y=~CALIT_OUT, type='scatter', mode='lines',
            line=list(color='#f39c12', width=2.5)) %>%
      layout(xaxis=list(title="Week"),
             yaxis=list(title="1,25-(OH)₂D₃ (pg/mL)"),
             shapes=list(
               list(type='rect', x0=0, x1=input$sim_weeks, y0=18, y1=72,
                    fillcolor='rgba(46,204,113,0.12)',
                    line=list(color='rgba(46,204,113,0.4)', width=1))
             ),
             margin=list(t=20))
  })

  output$ca_plot <- renderPlotly({
    out <- sim_result()
    if (is.null(out)) return(plot_ly())
    plot_ly(out, x=~week, y=~CA_OUT, type='scatter', mode='lines',
            fill='tozeroy', fillcolor='rgba(241,196,15,0.10)',
            line=list(color='#d4ac0d', width=2.5)) %>%
      layout(xaxis=list(title="Week"),
             yaxis=list(title="Serum Calcium (mg/dL)", range=c(8, 13)),
             shapes=list(
               list(type='line', x0=0, x1=input$sim_weeks, y0=10.5, y1=10.5,
                    line=list(color='#e74c3c', width=2, dash='dash'))
             ),
             annotations=list(list(x=5, y=10.7, text="Hypercalcemia threshold 10.5",
                                   showarrow=FALSE,
                                   font=list(color='#e74c3c', size=11))),
             margin=list(t=20))
  })

  # ---- TAB 5: Pulmonary Function ----
  output$fvc_plot <- renderPlotly({
    out <- sim_result()
    if (is.null(out)) return(plot_ly())
    fvc_threshold <- input$fvc_base - 10

    plot_ly(out, x=~week, y=~FVC_OUT, type='scatter', mode='lines',
            line=list(color='#27ae60', width=2.5)) %>%
      layout(xaxis=list(title="Week"),
             yaxis=list(title="FVC %predicted"),
             shapes=list(
               list(type='line', x0=0, x1=input$sim_weeks,
                    y0=fvc_threshold, y1=fvc_threshold,
                    line=list(color='#e74c3c', width=1.5, dash='dash'))
             ),
             annotations=list(
               list(x=5, y=fvc_threshold-1.5,
                    text=sprintf("Sig. decline threshold: %.0f%%", fvc_threshold),
                    showarrow=FALSE, font=list(color='#e74c3c', size=11))
             ),
             margin=list(t=20))
  })

  output$dlco_plot <- renderPlotly({
    out <- sim_result()
    if (is.null(out)) return(plot_ly())
    plot_ly(out, x=~week, y=~DLCO_OUT, type='scatter', mode='lines',
            line=list(color='#1abc9c', width=2.5)) %>%
      layout(xaxis=list(title="Week"),
             yaxis=list(title="DLCO %predicted (estimated)"),
             margin=list(t=20))
  })

  output$mmrc_plot <- renderPlotly({
    out <- sim_result()
    if (is.null(out)) return(plot_ly())
    plot_ly(out, x=~week, y=~mMRC_OUT, type='scatter', mode='lines',
            line=list(color='#8e44ad', width=2.5)) %>%
      layout(xaxis=list(title="Week"),
             yaxis=list(title="mMRC Dyspnea Score", range=c(-0.2, 4.2),
                        tickvals=0:4, ticktext=c("0-None","1-Mild","2-Moderate","3-Severe","4-Very Severe")),
             margin=list(t=20))
  })

  output$stage_bar_plot <- renderPlotly({
    out <- sim_result()
    if (is.null(out)) return(plot_ly())
    tp_wk  <- as.numeric(input$stage_timepoints)
    tp_row <- out[which.min(abs(out$week - tp_wk)), ]

    # Assign Scadding stage from FVC
    fvc_now  <- tp_row$FVC_OUT
    gran_now <- tp_row$GRAN_OUT

    stages <- c("Stage I","Stage II","Stage III","Stage IV")
    # Very rough probabilistic assignment
    if (gran_now < 1.2 && fvc_now > 85)
      probs <- c(0.70,0.25,0.04,0.01)
    else if (gran_now < 1.8 && fvc_now > 75)
      probs <- c(0.20,0.60,0.18,0.02)
    else if (gran_now < 2.5 && fvc_now > 60)
      probs <- c(0.05,0.25,0.55,0.15)
    else
      probs <- c(0.02,0.10,0.33,0.55)

    df_stage <- data.frame(Stage=stages, Probability=probs)
    plot_ly(df_stage, x=~Stage, y=~Probability,
            type='bar', marker=list(color=c('#2ecc71','#f1c40f','#e67e22','#e74c3c'))) %>%
      layout(xaxis=list(title="Scadding Stage"),
             yaxis=list(title=sprintf("Estimated Distribution (wk %g)", tp_wk),
                        range=c(0,1), tickformat='.0%'),
             margin=list(t=20))
  })

  # ---- TAB 6: Scenario Comparison ----
  all_scenario_results <- reactive({
    selected <- input$scenarios
    if (length(selected) == 0) return(list())

    results <- list()
    for (sc_name in selected) {
      sc <- scenarios_def[[sc_name]]
      res <- tryCatch(
        run_sim(sarc_mod,
                pred_dose        = sc$pred,
                mtx_dose         = sc$mtx,
                aza_dose         = sc$aza,
                ifx_on           = sc$ifx_on,
                ifx_dose         = sc$ifx_dose,
                ada_on           = sc$ada_on,
                hcq_dose         = sc$hcq,
                bw               = input$bw,
                taper_on         = sc$taper,
                taper_wk         = sc$taper_wk,
                sim_weeks        = input$sim_weeks,
                ace_baseline     = input$ace_base,
                fvc_baseline     = input$fvc_base,
                disease_stage    = input$disease_stage,
                disease_activity = input$disease_activity),
        error = function(e) NULL
      )
      if (!is.null(res)) {
        res$Scenario <- sc_name
        results[[sc_name]] <- res
      }
    }
    results
  })

  scenario_colors <- c(
    "No treatment" = "#e74c3c",
    "Pred 40mg"    = "#3498db",
    "Pred+MTX"     = "#2ecc71",
    "Infliximab"   = "#9b59b6",
    "HCQ"          = "#f39c12",
    "Pred taper"   = "#1abc9c"
  )

  output$scenario_gran_plot <- renderPlotly({
    res_list <- all_scenario_results()
    if (length(res_list) == 0) return(plot_ly())

    p <- plot_ly()
    for (sc_name in names(res_list)) {
      df_sc <- res_list[[sc_name]]
      col   <- scenario_colors[sc_name]
      if (is.na(col)) col <- "#95a5a6"
      p <- p %>% add_trace(data=df_sc, x=~week, y=~GRAN_OUT,
                           name=sc_name, type='scatter', mode='lines',
                           line=list(color=col, width=2))
    }
    p %>% layout(xaxis=list(title="Week"),
                 yaxis=list(title="Granuloma Burden (normalized)"),
                 legend=list(orientation='h', y=-0.2),
                 hovermode='x unified', margin=list(t=20))
  })

  output$scenario_ace_plot <- renderPlotly({
    res_list <- all_scenario_results()
    if (length(res_list) == 0) return(plot_ly())

    p <- plot_ly()
    for (sc_name in names(res_list)) {
      df_sc <- res_list[[sc_name]]
      col   <- scenario_colors[sc_name]
      if (is.na(col)) col <- "#95a5a6"
      p <- p %>% add_trace(data=df_sc, x=~week, y=~ACE_OUT,
                           name=sc_name, type='scatter', mode='lines',
                           line=list(color=col, width=2))
    }
    p %>% layout(xaxis=list(title="Week"),
                 yaxis=list(title="Serum ACE (U/L)"),
                 shapes=list(
                   list(type='line', x0=0, x1=input$sim_weeks, y0=67, y1=67,
                        line=list(color='#e67e22', width=1, dash='dash'))
                 ),
                 legend=list(orientation='h', y=-0.2),
                 hovermode='x unified', margin=list(t=20))
  })

  output$scenario_tbl <- renderDT({
    res_list <- all_scenario_results()
    if (length(res_list) == 0) return(datatable(data.frame()))

    rows <- lapply(names(res_list), function(sc_name) {
      df_sc <- res_list[[sc_name]]
      last  <- tail(df_sc, 1)
      data.frame(
        Scenario       = sc_name,
        `Granuloma`    = sprintf("%.2f", last$GRAN_OUT),
        `ACE (U/L)`    = sprintf("%.1f", last$ACE_OUT),
        `sIL-2R (U/mL)`= sprintf("%.0f", last$SIL2R_OUT),
        `FVC (%pred)`  = sprintf("%.1f", last$FVC_OUT),
        `DLCO (%pred)` = sprintf("%.1f", last$DLCO_OUT),
        `Calcium (mg/dL)` = sprintf("%.2f", last$CA_OUT),
        check.names = FALSE
      )
    })
    df_out <- do.call(rbind, rows)
    datatable(df_out, rownames=FALSE, options=list(dom='t', paging=FALSE),
              class='compact stripe') %>%
      formatStyle(columns='Granuloma',
                  backgroundColor=styleInterval(c(1.0, 1.5, 2.0),
                                                 c('#2ecc7140','#f1c40f40','#e67e2240','#e74c3c40')))
  })

  output$response_tbl <- renderDT({
    res_list <- all_scenario_results()
    if (length(res_list) == 0) return(datatable(data.frame()))

    rows <- lapply(names(res_list), function(sc_name) {
      df_sc <- res_list[[sc_name]]
      last  <- tail(df_sc, 1)
      init_gran <- df_sc$GRAN_OUT[1]
      pct_change <- (last$GRAN_OUT - init_gran) / init_gran * 100

      response <- dplyr::case_when(
        last$GRAN_OUT < 1.1 & last$FVC_OUT > input$fvc_base - 5  ~ "Remission",
        pct_change < -20                                           ~ "Partial Response",
        abs(pct_change) <= 20                                      ~ "Stable Disease",
        TRUE                                                       ~ "Progressive Disease"
      )
      color_map <- c("Remission"="#27ae60","Partial Response"="#f39c12",
                     "Stable Disease"="#3498db","Progressive Disease"="#e74c3c")

      data.frame(Scenario=sc_name,
                 `Granuloma Change`=sprintf("%+.1f%%", pct_change),
                 Response=response, check.names=FALSE)
    })
    df_out <- do.call(rbind, rows)
    datatable(df_out, rownames=FALSE, options=list(dom='t', paging=FALSE),
              class='compact stripe') %>%
      formatStyle('Response',
                  backgroundColor=styleEqual(
                    c("Remission","Partial Response","Stable Disease","Progressive Disease"),
                    c('#2ecc7140','#f1c40f40','#3498db40','#e74c3c40')
                  ))
  })

  # ---- TAB 7: Biomarker Correlation ----
  output$corr_gran_ace <- renderPlotly({
    out <- sim_result()
    if (is.null(out)) return(plot_ly())
    cor_val <- round(cor(out$GRAN_OUT, out$ACE_OUT, use="complete.obs"), 3)
    plot_ly(out, x=~GRAN_OUT, y=~ACE_OUT,
            type='scatter', mode='markers',
            marker=list(color=~week, colorscale='Viridis', size=5, showscale=TRUE,
                        colorbar=list(title="Week")),
            text=~paste0("Wk: ", round(week,1))) %>%
      layout(xaxis=list(title="Granuloma Burden"),
             yaxis=list(title="Serum ACE (U/L)"),
             title=list(text=sprintf("r = %.3f", cor_val), font=list(size=12)),
             margin=list(t=35))
  })

  output$corr_gran_fvc <- renderPlotly({
    out <- sim_result()
    if (is.null(out)) return(plot_ly())
    cor_val <- round(cor(out$GRAN_OUT, out$FVC_OUT, use="complete.obs"), 3)
    plot_ly(out, x=~GRAN_OUT, y=~FVC_OUT,
            type='scatter', mode='markers',
            marker=list(color=~week, colorscale='Plasma', size=5, showscale=TRUE,
                        colorbar=list(title="Week")),
            text=~paste0("Wk: ", round(week,1))) %>%
      layout(xaxis=list(title="Granuloma Burden"),
             yaxis=list(title="FVC %predicted"),
             title=list(text=sprintf("r = %.3f", cor_val), font=list(size=12)),
             margin=list(t=35))
  })

  output$corr_calit_ca <- renderPlotly({
    out <- sim_result()
    if (is.null(out)) return(plot_ly())
    cor_val <- round(cor(out$CALIT_OUT, out$CA_OUT, use="complete.obs"), 3)
    plot_ly(out, x=~CALIT_OUT, y=~CA_OUT,
            type='scatter', mode='markers',
            marker=list(color=~week, colorscale='RdBu', size=5, showscale=TRUE,
                        colorbar=list(title="Week")),
            text=~paste0("Wk: ", round(week,1))) %>%
      layout(xaxis=list(title="Calcitriol (pg/mL)"),
             yaxis=list(title="Serum Calcium (mg/dL)"),
             title=list(text=sprintf("r = %.3f", cor_val), font=list(size=12)),
             margin=list(t=35))
  })

  output$corr_table <- renderDT({
    out <- sim_result()
    if (is.null(out)) return(datatable(data.frame()))

    vars <- list(
      "Granuloma" = out$GRAN_OUT,
      "ACE"       = out$ACE_OUT,
      "sIL-2R"    = out$SIL2R_OUT,
      "Calcitriol"= out$CALIT_OUT,
      "Calcium"   = out$CA_OUT,
      "FVC"       = out$FVC_OUT,
      "DLCO"      = out$DLCO_OUT,
      "mMRC"      = out$mMRC_OUT
    )
    nm <- names(vars)
    mat <- matrix(NA, nrow=length(nm), ncol=length(nm), dimnames=list(nm,nm))
    for (i in seq_along(nm))
      for (j in seq_along(nm))
        mat[i,j] <- round(cor(vars[[i]], vars[[j]], use="complete.obs"), 3)

    df_cor <- as.data.frame(mat)
    df_cor <- cbind(Variable=rownames(df_cor), df_cor)
    datatable(df_cor, rownames=FALSE,
              options=list(dom='t', paging=FALSE, scrollX=TRUE),
              class='compact stripe') %>%
      formatStyle(nm, color=styleInterval(c(-0.5, 0, 0.5),
                                           c('#e74c3c','#e67e22','#7f8c8d','#27ae60')))
  })

}  # end server

# ============================================================
# Run App
# ============================================================
shinyApp(ui = ui, server = server)
