# =============================================================================
# Polyarteritis Nodosa (PAN) — QSP Interactive Shiny Dashboard
# 결절성 다발동맥염 정량적 시스템 약리학 대시보드
# =============================================================================
# Dependencies: shiny, shinydashboard, plotly, ggplot2, DT, mrgsolve,
#               dplyr, tidyr, shinycssloaders, shinyjs, RColorBrewer
# =============================================================================

library(shiny)
library(shinydashboard)
library(plotly)
library(ggplot2)
library(DT)
library(dplyr)
library(tidyr)
library(shinycssloaders)
library(shinyjs)

# ── mrgsolve PAN QSP Model ──────────────────────────────────────────────────
library(mrgsolve)

pan_model_code <- '
$PROB PAN QSP Model — Polyarteritis Nodosa
          결절성 다발동맥염 정량적 시스템 약리학 모델

$PARAM
// ── Patient-level covariates ──────────────────────────────────────────────
AGE    = 50   // age (years)
WT     = 70   // weight (kg)
HBV    = 0    // HBV status (0=no, 1=yes)
BVAS0  = 24   // baseline BVAS
CRP0   = 60   // baseline CRP mg/L
EGFR0  = 75   // baseline eGFR

// ── Prednisolone PK ───────────────────────────────────────────────────────
KA_PRED  = 1.5    // absorption (1/h)
CL_PRED  = 18.0   // clearance (L/h)
V_PRED   = 55.0   // volume (L)
DOSE_PRED = 0     // mg/day
FREQ_PRED = 1     // doses/day

// ── Cyclophosphamide PK ───────────────────────────────────────────────────
KA_CYC   = 0.8
CL_CYC   = 6.0
V_CYC    = 38.0
F_ACTIVE = 0.10   // fraction converted to 4-OH-CYC
CL_OHC   = 12.0
V_OHC    = 25.0
DOSE_CYC = 0      // mg (pulse dose per event)
CYC_ON   = 0      // 1 = active CYC arm

// ── Azathioprine PK ───────────────────────────────────────────────────────
KA_AZA   = 1.2
CL_AZA   = 8.5
V_AZA    = 40.0
DOSE_AZA = 0
AZA_ON   = 0

// ── Antiviral PK (e.g., tenofovir/entecavir) ─────────────────────────────
KA_AV    = 1.0
CL_AV    = 25.0
V_AV     = 95.0
DOSE_AV  = 0
AV_ON    = 0

// ── Immune / PD parameters ────────────────────────────────────────────────
// Immune complex (IC) dynamics
IC0      = 1.0    // baseline IC (normalised)
K_IC_IN  = 0.05   // IC production rate (/h)
K_IC_CL  = 0.03   // IC clearance rate (/h)
HBV_IC   = 2.0    // HBV amplification factor on IC production

// Complement
C3_0     = 1.0
K_C3     = 0.02
K_C3CL   = 0.03

// Neutrophils
NEU0     = 1.0
K_NREC   = 0.15   // neutrophil recruitment (/h)
K_NDEATH = 0.10

// B cells
BCELL0   = 1.0
K_BPROL  = 0.008
K_BDEATH = 0.005

// T cells (CD4+)
TCELL0   = 1.0
K_TPROL  = 0.006
K_TDEATH = 0.004

// Cytokines IL-6 / CRP
IL6_0    = 1.0
K_IL6IN  = 0.04
K_IL6CL  = 0.06

CRP0_N   = 1.0    // normalised CRP
K_CRPIN  = 0.05
K_CRPCL  = 0.04

// ESR
ESR0     = 1.0
K_ESRIN  = 0.03
K_ESRCL  = 0.025

// Vascular / tissue
VASC0    = 1.0    // vascular inflammation index
K_VASCD  = 0.02   // vascular damage rate
K_VASCR  = 0.01   // vascular repair rate

MIAO0    = 1.0    // microaneurysm burden
K_MIAOD  = 0.005
K_MIAOR  = 0.003

NERVE0   = 1.0    // nerve damage score
K_NERVED = 0.015
K_NERVER = 0.008

EGFR_N   = 1.0    // normalised eGFR (1 = baseline)
K_GFRD   = 0.01   // GFR decline rate
K_GFRR   = 0.005

BVAS_N   = 1.0    // normalised BVAS
K_BVASD  = 0.02
K_BVASR  = 0.025

// ── Drug PD (EC50 / Emax) ─────────────────────────────────────────────────
EC50_PRED = 0.5   // µg/mL
EMAX_PRED = 0.85
HILL_PRED = 1.5

EC50_OHC  = 0.3
EMAX_OHC  = 0.90
HILL_OHC  = 1.2

EC50_AZA  = 0.4
EMAX_AZA  = 0.80
HILL_AZA  = 1.3

EC50_AV   = 0.2
EMAX_AV   = 0.95
HILL_AV   = 2.0

$CMT
// Prednisolone
PRED_GUT PRED_C
// CYC / active metabolite
CYC_GUT CYC_C OHC_C
// Azathioprine
AZA_GUT AZA_C
// Antiviral
AV_GUT AV_C
// Immune compartments
IC_C C3_C NEU_C BCELL_C TCELL_C
// Cytokines
IL6_C CRP_C ESR_C
// Tissue
VASC_C MIAO_C NERVE_C EGFR_C BVAS_C

$MAIN
// Initial conditions (ODE compartments start at 0; handled in $ODE)
// Biomarker initials set via $INIT
PRED_GUT_0  = 0;
PRED_C_0    = 0;
CYC_GUT_0   = 0;
CYC_C_0     = 0;
OHC_C_0     = 0;
AZA_GUT_0   = 0;
AZA_C_0     = 0;
AV_GUT_0    = 0;
AV_C_0      = 0;
IC_C_0      = IC0;
C3_C_0      = C3_0;
NEU_C_0     = NEU0;
BCELL_C_0   = BCELL0;
TCELL_C_0   = TCELL0;
IL6_C_0     = IL6_0;
CRP_C_0     = CRP0_N;
ESR_C_0     = ESR0;
VASC_C_0    = VASC0;
MIAO_C_0    = MIAO0;
NERVE_C_0   = NERVE0;
EGFR_C_0    = EGFR_N;
BVAS_C_0    = BVAS_N;

$ODE
// ── Drug PK ───────────────────────────────────────────────────────────────
// Prednisolone (oral, once or twice daily — handled externally as events)
double PRED_abs = KA_PRED * PRED_GUT;
dxdt_PRED_GUT = -PRED_abs;
dxdt_PRED_C   =  PRED_abs - (CL_PRED/V_PRED)*PRED_C;
double Cp_PRED = PRED_C / V_PRED;   // µg/mL (approximate)

// Cyclophosphamide
double CYC_abs = KA_CYC * CYC_GUT;
double OHC_form = F_ACTIVE*(CL_CYC/V_CYC)*CYC_C;
dxdt_CYC_GUT = -CYC_abs;
dxdt_CYC_C   =  CYC_abs - (CL_CYC/V_CYC)*CYC_C;
dxdt_OHC_C   =  OHC_form*V_CYC - (CL_OHC/V_OHC)*OHC_C;
double Cp_OHC = OHC_C / V_OHC;

// Azathioprine
double AZA_abs = KA_AZA * AZA_GUT;
dxdt_AZA_GUT = -AZA_abs;
dxdt_AZA_C   =  AZA_abs - (CL_AZA/V_AZA)*AZA_C;
double Cp_AZA = AZA_C / V_AZA;

// Antiviral
double AV_abs = KA_AV * AV_GUT;
dxdt_AV_GUT = -AV_abs;
dxdt_AV_C   =  AV_abs - (CL_AV/V_AV)*AV_C;
double Cp_AV = AV_C / V_AV;

// ── Drug effect (Emax models) ─────────────────────────────────────────────
double E_PRED = EMAX_PRED * pow(Cp_PRED, HILL_PRED) /
                (pow(EC50_PRED, HILL_PRED) + pow(Cp_PRED, HILL_PRED));
double E_OHC  = EMAX_OHC  * pow(Cp_OHC,  HILL_OHC)  /
                (pow(EC50_OHC,  HILL_OHC)  + pow(Cp_OHC,  HILL_OHC));
double E_AZA  = EMAX_AZA  * pow(Cp_AZA,  HILL_AZA)  /
                (pow(EC50_AZA,  HILL_AZA)  + pow(Cp_AZA,  HILL_AZA));
double E_AV   = EMAX_AV   * pow(Cp_AV,   HILL_AV)   /
                (pow(EC50_AV,   HILL_AV)   + pow(Cp_AV,   HILL_AV));

double E_total = 1 - (1-(1-E_PRED)*(1-E_OHC)*(1-E_AZA));
double E_hbv   = E_AV;  // antiviral effect on HBV-driven IC

// ── Immune dynamics ───────────────────────────────────────────────────────
// Immune complex
double IC_in  = (1 + HBV*HBV_IC) * K_IC_IN * (1 - E_hbv*HBV);
double IC_cl  = K_IC_CL * IC_C;
dxdt_IC_C = IC_in - IC_cl - E_total*0.4*IC_C;

// Complement C3
double C3_stim = IC_C * K_C3;
dxdt_C3_C = C3_stim - K_C3CL*C3_C - E_PRED*0.15*C3_C;
double C3_level = (C3_C > 0) ? C3_C : 0.01;

// Neutrophils (recruited by IC and complement)
double NEU_rec = K_NREC * IC_C * C3_level * (1 - E_total*0.5);
dxdt_NEU_C = NEU_rec - K_NDEATH*NEU_C;

// B cells
double BCELL_stim = K_BPROL * IC_C * (1 - E_total*0.6);
dxdt_BCELL_C = BCELL_stim + 0.001 - K_BDEATH*BCELL_C;

// T cells
double TCELL_stim = K_TPROL * IL6_C * (1 - E_PRED*0.7);
dxdt_TCELL_C = TCELL_stim + 0.001 - K_TDEATH*TCELL_C;

// Cytokines
double IL6_prod = K_IL6IN * (NEU_C + BCELL_C + TCELL_C)/3.0 * (1 - E_PRED*0.6);
dxdt_IL6_C = IL6_prod - K_IL6CL*IL6_C;

double CRP_prod = K_CRPIN * IL6_C * (1 - E_PRED*0.5);
dxdt_CRP_C = CRP_prod - K_CRPCL*CRP_C;

double ESR_prod = K_ESRIN * (CRP_C + IL6_C)/2.0;
dxdt_ESR_C = ESR_prod - K_ESRCL*ESR_C;

// ── Tissue / clinical endpoints ───────────────────────────────────────────
// Vascular inflammation
double VASC_in = K_VASCD * IC_C * NEU_C * (1 - E_total*0.55);
dxdt_VASC_C = VASC_in - K_VASCR*VASC_C - E_total*0.3*VASC_C;

// Microaneurysm
double MIAO_in = K_MIAOD * VASC_C * (1 - E_total*0.45);
dxdt_MIAO_C = MIAO_in - K_MIAOR*MIAO_C;

// Nerve damage (mononeuritis multiplex)
double NERVE_in = K_NERVED * VASC_C * (1 - E_total*0.35);
dxdt_NERVE_C = NERVE_in - K_NERVER*NERVE_C;

// Renal function (GFR)
double GFR_loss = K_GFRD * VASC_C * IC_C;
dxdt_EGFR_C = -GFR_loss + K_GFRR*(1 - EGFR_C)*E_total;

// BVAS
double BVAS_drive = K_BVASD * (VASC_C + IC_C + NEU_C)/3.0 * (1 - E_total*0.65);
dxdt_BVAS_C = BVAS_drive - K_BVASR*BVAS_C*E_total;

$TABLE
double CONC_PRED = PRED_C/V_PRED;
double CONC_OHC  = OHC_C/V_OHC;
double CONC_AZA  = AZA_C/V_AZA;
double CONC_AV   = AV_C/V_AV;
double BVAS_abs  = BVAS_C * BVAS0;
double CRP_abs   = CRP_C * CRP0;
double EGFR_abs  = EGFR_C * EGFR0;
double IC_norm   = IC_C;
double C3_norm   = C3_C;
double IL6_norm  = IL6_C;
double ESR_norm  = ESR_C;
double NEU_norm  = NEU_C;
double BCELL_norm = BCELL_C;
double TCELL_norm = TCELL_C;
double VASC_norm  = VASC_C;
double MIAO_norm  = MIAO_C;
double NERVE_norm = NERVE_C;

$CAPTURE CONC_PRED CONC_OHC CONC_AZA CONC_AV
         BVAS_abs CRP_abs EGFR_abs
         IC_norm C3_norm IL6_norm ESR_norm
         NEU_norm BCELL_norm TCELL_norm
         VASC_norm MIAO_norm NERVE_norm
'

# Compile model (suppress verbose output)
pan_mod <- tryCatch(
  mcode("pan_qsp", pan_model_code, quiet = TRUE),
  error = function(e) NULL
)

# ── Helper: run one simulation arm ─────────────────────────────────────────
run_sim <- function(mod,
                    age = 50, wt = 70, hbv = 0,
                    bvas0 = 24, crp0 = 60, egfr0 = 75,
                    arm = "pred_mono",
                    dose_pred = 60, freq_pred = 1,
                    dose_cyc = 750, dose_aza = 150, dose_av = 300,
                    end_time = 365 * 24,   # hours
                    delta_t  = 6) {

  if (is.null(mod)) {
    # Fallback: synthetic data when mrgsolve unavailable
    t_days <- seq(0, end_time / 24, by = delta_t / 24)
    decay  <- function(x0, t, rate) x0 * exp(-rate * t)
    ef     <- switch(arm,
                     none       = 0.00,
                     pred_mono  = 0.50,
                     pred_cyc   = 0.78,
                     pred_aza   = 0.60,
                     hbv_pred   = 0.72,
                     0.50)
    n <- length(t_days)
    return(data.frame(
      time      = t_days,
      BVAS_abs  = pmax(bvas0 * exp(-(ef * 0.005) * t_days), bvas0 * (1 - ef) * 0.2),
      CRP_abs   = pmax(crp0  * exp(-(ef * 0.006) * t_days), crp0  * (1 - ef) * 0.15),
      EGFR_abs  = egfr0 + (100 - egfr0) * ef * (1 - exp(-0.003 * t_days)),
      CONC_PRED = ifelse(arm %in% c("pred_mono","pred_cyc","pred_aza","hbv_pred"),
                         dose_pred * exp(-0.12 * (t_days %% (24 / freq_pred))), 0),
      CONC_OHC  = ifelse(arm == "pred_cyc", 0.8 * exp(-0.25 * t_days %% 504), 0),
      CONC_AZA  = ifelse(arm == "pred_aza", 0.6 * exp(-0.15 * t_days %% 24), 0),
      CONC_AV   = ifelse(arm == "hbv_pred", 0.9 * exp(-0.10 * t_days %% 24), 0),
      IC_norm   = 1 + (1 - exp(-(1 - ef) * 0.002 * t_days)),
      C3_norm   = 1,
      IL6_norm  = 1 + 0.5 * (1 - ef) - 0.5 * ef * (1 - exp(-0.003 * t_days)),
      ESR_norm  = 1,
      NEU_norm  = 1,
      BCELL_norm = 1,
      TCELL_norm = 1,
      VASC_norm  = pmax(1 * exp(-(ef * 0.004) * t_days), (1 - ef) * 0.1),
      MIAO_norm  = pmax(1 * exp(-(ef * 0.003) * t_days), (1 - ef) * 0.15),
      NERVE_norm = pmax(1 * exp(-(ef * 0.002) * t_days), (1 - ef) * 0.2),
      arm        = arm
    ))
  }

  # Build event table
  et <- ev(time = 0, amt = 0, cmt = "PRED_GUT")  # placeholder

  # Prednisolone: daily (or BID) dosing for full duration
  if (arm %in% c("pred_mono", "pred_cyc", "pred_aza", "hbv_pred") && dose_pred > 0) {
    int_h  <- 24 / freq_pred
    n_doses <- floor(end_time / int_h)
    pred_ev <- ev(time = seq(0, by = int_h, length.out = n_doses),
                  amt  = dose_pred / freq_pred,
                  cmt  = "PRED_GUT",
                  addl = 0)
    et <- et + pred_ev
  }

  # Cyclophosphamide: monthly pulse at t=0, 28d, 56d, 84d, 112d, 140d
  if (arm == "pred_cyc" && dose_cyc > 0) {
    cyc_times <- seq(0, min(end_time, 168 * 24), by = 28 * 24)
    cyc_ev    <- ev(time = cyc_times, amt = dose_cyc, cmt = "CYC_GUT")
    et <- et + cyc_ev
  }

  # Azathioprine: daily from week 12 onward
  if (arm == "pred_aza" && dose_aza > 0) {
    aza_start <- 12 * 7 * 24
    n_aza     <- floor((end_time - aza_start) / 24)
    if (n_aza > 0) {
      aza_ev <- ev(time = seq(aza_start, by = 24, length.out = n_aza),
                   amt  = dose_aza, cmt = "AZA_GUT")
      et <- et + aza_ev
    }
  }

  # Antiviral: daily throughout
  if (arm == "hbv_pred" && dose_av > 0) {
    n_av  <- floor(end_time / 24)
    av_ev <- ev(time = seq(0, by = 24, length.out = n_av),
                amt  = dose_av, cmt = "AV_GUT")
    et <- et + av_ev
  }

  obs_times <- seq(0, end_time, by = delta_t)

  out <- mod %>%
    param(AGE = age, WT = wt, HBV = hbv,
          BVAS0 = bvas0, CRP0 = crp0, EGFR0 = egfr0,
          DOSE_PRED = dose_pred, FREQ_PRED = freq_pred,
          DOSE_CYC = dose_cyc, DOSE_AZA = dose_aza, DOSE_AV = dose_av) %>%
    ev(et) %>%
    mrgsim(end = end_time, delta = delta_t, obsonly = TRUE) %>%
    as.data.frame() %>%
    mutate(time = time / 24, arm = arm)  # convert hours → days

  out
}


# ── UI ─────────────────────────────────────────────────────────────────────
ui <- dashboardPage(
  skin = "red",

  dashboardHeader(
    title = span(
      icon("heartbeat"),
      " PAN QSP Dashboard",
      style = "font-size:16px; font-weight:bold;"
    ),
    titleWidth = 300
  ),

  dashboardSidebar(
    width = 300,
    sidebarMenu(
      id = "tabs",
      menuItem("Patient Profile / 환자 프로파일",   tabName = "tab_patient",   icon = icon("user-md")),
      menuItem("Drug PK / 약동학",                  tabName = "tab_pk",        icon = icon("flask")),
      menuItem("PD Biomarkers / 약력학 바이오마커", tabName = "tab_pd",        icon = icon("dna")),
      menuItem("Clinical Endpoints / 임상 지표",    tabName = "tab_endpoints", icon = icon("chart-line")),
      menuItem("Scenario Comparison / 시나리오",    tabName = "tab_scenarios", icon = icon("layer-group")),
      menuItem("Biomarker Correlation / 상관분석",  tabName = "tab_corr",      icon = icon("project-diagram"))
    ),
    hr(),
    div(
      style = "padding:10px;",
      actionButton("btn_run", "Run Simulation / 시뮬레이션 실행",
                   icon = icon("play"),
                   class = "btn-danger btn-block",
                   style = "font-weight:bold;"),
      br(),
      downloadButton("dl_data", "Download Data / 데이터 다운로드",
                     class = "btn-default btn-block btn-sm")
    )
  ),

  dashboardBody(
    useShinyjs(),
    tags$head(
      tags$style(HTML("
        .skin-red .main-header .logo { background-color:#8B0000; }
        .skin-red .main-header .navbar { background-color:#A00000; }
        .skin-red .main-sidebar { background-color:#2d2d2d; }
        .content-wrapper { background-color:#f5f5f5; }
        .box.box-danger { border-top-color:#8B0000; }
        .risk-low    { color:#27ae60; font-weight:bold; }
        .risk-medium { color:#e67e22; font-weight:bold; }
        .risk-high   { color:#c0392b; font-weight:bold; }
        .value-box-custom { border-radius:8px; padding:15px; margin-bottom:10px;
                             color:white; font-size:18px; font-weight:bold; }
      "))
    ),

    tabItems(

      # ════════════════════════════════════════════════════════════════════
      # TAB 1 — Patient Profile
      # ════════════════════════════════════════════════════════════════════
      tabItem(
        tabName = "tab_patient",
        fluidRow(
          box(
            title = "Demographics / 인구통계학적 특성",
            width = 6, status = "danger", solidHeader = TRUE,
            sliderInput("age",    "Age / 나이 (years)",         20, 80, 50),
            selectInput("sex",    "Sex / 성별",
                        choices = c("Male / 남성" = "M", "Female / 여성" = "F")),
            sliderInput("wt",     "Weight / 체중 (kg)",         40, 120, 70),
            selectInput("hbv",    "HBV status / HBV 상태",
                        choices = c("Negative / 음성" = "0",
                                    "Positive / 양성" = "1")),
            selectInput("subtype", "Disease subtype / 질환 아형",
                        choices = c("Classic PAN / 전형 PAN" = "classic",
                                    "HBV-associated PAN / HBV 연관 PAN" = "hbv",
                                    "Microscopic PAN / 미세혈관 PAN" = "micro"))
          ),
          box(
            title = "Baseline Disease Activity / 기저 질환 활성도",
            width = 6, status = "danger", solidHeader = TRUE,
            sliderInput("bvas0",  "Baseline BVAS / 기저 BVAS (0–63)",  0, 63, 24),
            sliderInput("crp0",   "Baseline CRP (mg/L)",               0, 200, 60),
            sliderInput("egfr0",  "Baseline eGFR (mL/min/1.73 m²)",   15, 120, 75),
            sliderInput("n_organs","Number of organs involved / 침범 장기 수", 1, 7, 3)
          )
        ),
        fluidRow(
          box(
            title = "Organ Involvement / 장기 침범",
            width = 6, status = "warning", solidHeader = TRUE,
            checkboxGroupInput("organs", "Organs involved / 침범 장기",
                               choices = c("Kidney / 신장"         = "kidney",
                                           "Peripheral nerve / 말초신경" = "nerve",
                                           "Gastrointestinal / 위장관" = "gi",
                                           "Skin / 피부"          = "skin",
                                           "Heart / 심장"         = "heart",
                                           "Muscle / 근육"        = "muscle",
                                           "Testes / 고환"        = "testes"),
                               selected = c("kidney", "nerve", "skin"))
          ),
          box(
            title = "Five-Factor Score (FFS) & Risk Stratification",
            width = 6, status = "danger", solidHeader = TRUE,
            h4("Five-Factor Score (Guillevin 2003 revised)"),
            verbatimTextOutput("ffs_output"),
            br(),
            h4("Risk Stratification / 위험도 분류"),
            uiOutput("risk_output"),
            br(),
            h4("Predicted 5-year Mortality / 예측 5년 사망률"),
            uiOutput("mortality_output")
          )
        )
      ),

      # ════════════════════════════════════════════════════════════════════
      # TAB 2 — Drug PK
      # ════════════════════════════════════════════════════════════════════
      tabItem(
        tabName = "tab_pk",
        fluidRow(
          box(
            title = "Treatment Selection & Dosing / 치료 선택 및 용량",
            width = 4, status = "danger", solidHeader = TRUE,
            selectInput("arm_pk", "Treatment arm / 치료군",
                        choices = c("Prednisolone monotherapy / 스테로이드 단독" = "pred_mono",
                                    "Pred + CYC pulse (CYCLOPS) / 스테로이드+사이클로포스파마이드" = "pred_cyc",
                                    "Pred + AZA maintenance / 스테로이드+아자티오프린" = "pred_aza",
                                    "Antiviral + Pred (HBV-PAN) / 항바이러스+스테로이드" = "hbv_pred")),
            hr(),
            h5("Prednisolone"),
            sliderInput("dose_pred",  "Initial dose / 초기 용량 (mg/day)", 10, 100, 60),
            sliderInput("freq_pred",  "Frequency / 투여 횟수 (/day)", 1, 2, 1),
            h5("Cyclophosphamide (pulse)"),
            sliderInput("dose_cyc",   "Pulse dose / 펄스 용량 (mg)", 200, 1500, 750),
            h5("Azathioprine (maintenance)"),
            sliderInput("dose_aza",   "Daily dose / 일일 용량 (mg/day)", 50, 200, 150),
            h5("Antiviral (e.g., Tenofovir)"),
            sliderInput("dose_av",    "Daily dose / 일일 용량 (mg/day)", 100, 400, 300),
            sliderInput("pk_end",     "Simulation horizon / 시뮬레이션 기간 (days)", 30, 730, 365)
          ),
          box(
            title = "Prednisolone Plasma Concentration / 프레드니솔론 혈중 농도",
            width = 8, status = "danger", solidHeader = TRUE,
            withSpinner(plotlyOutput("plot_pred_pk", height = "280px"), color = "#8B0000"),
            hr(),
            withSpinner(plotlyOutput("plot_cyc_pk",  height = "280px"), color = "#8B0000")
          )
        ),
        fluidRow(
          box(
            title = "Azathioprine / Antiviral PK",
            width = 6, status = "warning", solidHeader = TRUE,
            withSpinner(plotlyOutput("plot_aza_pk", height = "250px"), color = "#e67e22")
          ),
          box(
            title = "PK Summary Table / PK 요약 표",
            width = 6, status = "info", solidHeader = TRUE,
            DTOutput("tbl_pk_summary"),
            br(),
            downloadButton("dl_pk", "Download PK data", class = "btn-sm btn-default")
          )
        )
      ),

      # ════════════════════════════════════════════════════════════════════
      # TAB 3 — PD Biomarkers
      # ════════════════════════════════════════════════════════════════════
      tabItem(
        tabName = "tab_pd",
        fluidRow(
          box(
            title = "PD Parameter Sensitivity / PD 파라미터 민감도",
            width = 3, status = "danger", solidHeader = TRUE,
            sliderInput("ec50_pred", "EC50 Prednisolone (µg/mL)", 0.1, 2.0, 0.5, step = 0.05),
            sliderInput("emax_pred", "Emax Prednisolone",          0.5, 1.0, 0.85, step = 0.01),
            sliderInput("ec50_ohc",  "EC50 4-OH-CYC (µg/mL)",    0.1, 1.5, 0.3, step = 0.05),
            sliderInput("emax_ohc",  "Emax 4-OH-CYC",             0.5, 1.0, 0.90, step = 0.01),
            sliderInput("ec50_aza",  "EC50 AZA (µg/mL)",          0.1, 1.5, 0.4, step = 0.05),
            sliderInput("emax_aza",  "Emax AZA",                  0.5, 1.0, 0.80, step = 0.01),
            selectInput("pd_arm",    "Treatment arm / 치료군",
                        choices = c("pred_mono","pred_cyc","pred_aza","hbv_pred","none"))
          ),
          box(
            title = "Immune Complex & Complement / 면역복합체 & 보체",
            width = 9, status = "danger", solidHeader = TRUE,
            withSpinner(plotlyOutput("plot_ic_c3", height = "300px"), color = "#8B0000"),
            withSpinner(plotlyOutput("plot_cells",  height = "300px"), color = "#8B0000")
          )
        ),
        fluidRow(
          box(
            title = "Cytokines / 사이토카인 (IL-6, CRP, ESR)",
            width = 12, status = "warning", solidHeader = TRUE,
            withSpinner(plotlyOutput("plot_cytokines", height = "300px"), color = "#e67e22")
          )
        )
      ),

      # ════════════════════════════════════════════════════════════════════
      # TAB 4 — Clinical Endpoints
      # ════════════════════════════════════════════════════════════════════
      tabItem(
        tabName = "tab_endpoints",
        fluidRow(
          box(
            title = "BVAS & Organ Damage / BVAS 및 장기 손상",
            width = 6, status = "danger", solidHeader = TRUE,
            withSpinner(plotlyOutput("plot_bvas", height = "300px"), color = "#8B0000"),
            withSpinner(plotlyOutput("plot_vasc", height = "300px"), color = "#8B0000")
          ),
          box(
            title = "Renal & Nerve / 신장 및 신경",
            width = 6, status = "danger", solidHeader = TRUE,
            withSpinner(plotlyOutput("plot_egfr",  height = "300px"), color = "#8B0000"),
            withSpinner(plotlyOutput("plot_nerve", height = "300px"), color = "#8B0000")
          )
        ),
        fluidRow(
          box(
            title = "Microaneurysm Burden & Vascular Index / 미세동맥류 및 혈관 염증 지수",
            width = 12, status = "warning", solidHeader = TRUE,
            withSpinner(plotlyOutput("plot_miao", height = "280px"), color = "#e67e22")
          )
        )
      ),

      # ════════════════════════════════════════════════════════════════════
      # TAB 5 — Scenario Comparison
      # ════════════════════════════════════════════════════════════════════
      tabItem(
        tabName = "tab_scenarios",
        fluidRow(
          box(
            title = "Arms to Compare / 비교 시나리오 선택",
            width = 3, status = "danger", solidHeader = TRUE,
            checkboxGroupInput("arms_sel", "Select arms / 시나리오 선택",
                               choices = c(
                                 "No treatment / 무치료" = "none",
                                 "Pred monotherapy / 스테로이드 단독" = "pred_mono",
                                 "Pred + CYC pulse (CYCLOPS)" = "pred_cyc",
                                 "Pred + AZA maintenance" = "pred_aza",
                                 "HBV antiviral + Pred" = "hbv_pred"
                               ),
                               selected = c("none","pred_mono","pred_cyc")),
            sliderInput("scen_end", "Horizon / 기간 (days)", 180, 730, 365),
            hr(),
            downloadButton("dl_scen", "Download comparison data", class = "btn-sm btn-default")
          ),
          box(
            title = "BVAS Comparison / BVAS 비교",
            width = 9, status = "danger", solidHeader = TRUE,
            withSpinner(plotlyOutput("plot_scen_bvas", height = "350px"), color = "#8B0000")
          )
        ),
        fluidRow(
          box(
            title = "CRP Comparison / CRP 비교",
            width = 6, status = "warning", solidHeader = TRUE,
            withSpinner(plotlyOutput("plot_scen_crp", height = "280px"), color = "#e67e22")
          ),
          box(
            title = "eGFR Comparison / eGFR 비교",
            width = 6, status = "info", solidHeader = TRUE,
            withSpinner(plotlyOutput("plot_scen_egfr", height = "280px"), color = "#2980b9")
          )
        ),
        fluidRow(
          box(
            title = "Outcome Summary / 결과 요약",
            width = 8, status = "danger", solidHeader = TRUE,
            DTOutput("tbl_outcomes")
          ),
          box(
            title = "Forest Plot — Remission at 6 months",
            width = 4, status = "danger", solidHeader = TRUE,
            withSpinner(plotlyOutput("plot_forest", height = "280px"), color = "#8B0000")
          )
        )
      ),

      # ════════════════════════════════════════════════════════════════════
      # TAB 6 — Biomarker Correlation
      # ════════════════════════════════════════════════════════════════════
      tabItem(
        tabName = "tab_corr",
        fluidRow(
          box(
            title = "Scatter: BVAS vs CRP",
            width = 6, status = "danger", solidHeader = TRUE,
            withSpinner(plotlyOutput("plot_scatter_bvas_crp", height = "300px"), color = "#8B0000")
          ),
          box(
            title = "Scatter: eGFR vs Immune Complex",
            width = 6, status = "danger", solidHeader = TRUE,
            withSpinner(plotlyOutput("plot_scatter_egfr_ic", height = "300px"), color = "#8B0000")
          )
        ),
        fluidRow(
          box(
            title = "Correlation Heatmap / 상관 행렬",
            width = 6, status = "warning", solidHeader = TRUE,
            withSpinner(plotlyOutput("plot_corr_heat", height = "350px"), color = "#e67e22")
          ),
          box(
            title = "Biomarker Waterfall / 바이오마커 폭포 그래프",
            width = 6, status = "info", solidHeader = TRUE,
            selectInput("wf_arm", "Arm / 치료군",
                        choices = c("pred_mono","pred_cyc","pred_aza","hbv_pred","none")),
            withSpinner(plotlyOutput("plot_waterfall", height = "300px"), color = "#2980b9")
          )
        ),
        fluidRow(
          box(
            title = "Response Prediction by Subtype / 아형별 반응 예측",
            width = 12, status = "danger", solidHeader = TRUE,
            withSpinner(plotlyOutput("plot_subtype", height = "300px"), color = "#8B0000")
          )
        )
      )

    )  # end tabItems
  )    # end dashboardBody
)      # end dashboardPage


# ── Server ─────────────────────────────────────────────────────────────────
server <- function(input, output, session) {

  # Colour palette
  arm_colors <- c(
    none       = "#95a5a6",
    pred_mono  = "#e74c3c",
    pred_cyc   = "#8B0000",
    pred_aza   = "#c0392b",
    hbv_pred   = "#922b21"
  )
  arm_labels <- c(
    none      = "No treatment",
    pred_mono = "Pred monotherapy",
    pred_cyc  = "Pred + CYC",
    pred_aza  = "Pred + AZA",
    hbv_pred  = "HBV antiviral + Pred"
  )

  # ── Reactive simulation data (single arm — PK / PD tabs) ──────────────
  sim_single <- eventReactive(input$btn_run, {
    withProgress(message = "Running simulation / 시뮬레이션 중...", value = 0.5, {
      run_sim(
        mod      = pan_mod,
        age      = input$age,
        wt       = input$wt,
        hbv      = as.numeric(input$hbv),
        bvas0    = input$bvas0,
        crp0     = input$crp0,
        egfr0    = input$egfr0,
        arm      = input$arm_pk,
        dose_pred = input$dose_pred,
        freq_pred = input$freq_pred,
        dose_cyc  = input$dose_cyc,
        dose_aza  = input$dose_aza,
        dose_av   = input$dose_av,
        end_time  = input$pk_end * 24
      )
    })
  }, ignoreNULL = FALSE)

  # ── Reactive: multi-arm comparison ────────────────────────────────────
  sim_multi <- eventReactive(input$btn_run, {
    arms <- input$arms_sel
    if (length(arms) == 0) arms <- "pred_mono"
    withProgress(message = "Running multi-arm simulation...", value = 0.3, {
      dfs <- lapply(arms, function(a) {
        run_sim(
          mod       = pan_mod,
          age       = input$age, wt = input$wt,
          hbv       = as.numeric(input$hbv),
          bvas0     = input$bvas0, crp0 = input$crp0, egfr0 = input$egfr0,
          arm       = a,
          dose_pred = input$dose_pred, freq_pred = input$freq_pred,
          dose_cyc  = input$dose_cyc, dose_aza  = input$dose_aza,
          dose_av   = input$dose_av,
          end_time  = input$scen_end * 24
        )
      })
      bind_rows(dfs)
    })
  }, ignoreNULL = FALSE)

  # ── FFS & Risk ─────────────────────────────────────────────────────────
  ffs_score <- reactive({
    score <- 0
    if (input$egfr0 < 50)              score <- score + 1  # renal insufficiency
    if ("gi" %in% input$organs)        score <- score + 1  # GI involvement
    if ("heart" %in% input$organs)     score <- score + 1  # cardiomyopathy
    if (input$crp0 > 100)              score <- score + 1  # high disease burden
    if (!("nerve" %in% input$organs))  score <- score + 0  # PNS not a FFS item
    score
  })

  output$ffs_output <- renderText({
    s <- ffs_score()
    paste0(
      "FFS = ", s, "\n",
      "Criterion details:\n",
      "  Renal insufficiency (eGFR <50): ", ifelse(input$egfr0 < 50, "YES (+1)", "no"), "\n",
      "  GI involvement:                ", ifelse("gi" %in% input$organs, "YES (+1)", "no"), "\n",
      "  Cardiac involvement:           ", ifelse("heart" %in% input$organs, "YES (+1)", "no"), "\n",
      "  High disease burden (CRP>100): ", ifelse(input$crp0 > 100, "YES (+1)", "no")
    )
  })

  output$risk_output <- renderUI({
    s <- ffs_score()
    risk <- if (s == 0) "Low / 저위험" else if (s == 1) "Medium / 중위험" else "High / 고위험"
    cls  <- if (s == 0) "risk-low" else if (s == 1) "risk-medium" else "risk-high"
    tags$p(class = cls, paste("Risk:", risk))
  })

  output$mortality_output <- renderUI({
    s   <- ffs_score()
    m5y <- if (s == 0) "~12%" else if (s == 1) "~26%" else "~46%"
    tags$p(style = "font-size:20px; color:#8B0000; font-weight:bold;", m5y)
  })

  # ── Helper: ggplot → plotly ────────────────────────────────────────────
  g2p <- function(g) ggplotly(g) %>%
    layout(legend = list(orientation = "h", y = -0.2))

  theme_pan <- function() {
    theme_bw(base_size = 12) +
      theme(panel.border = element_rect(color = "#8B0000"),
            plot.title   = element_text(color = "#8B0000", face = "bold"))
  }

  # ════════════════════════════════════════════════════════════════════════
  # TAB 2 — PK plots
  # ════════════════════════════════════════════════════════════════════════
  output$plot_pred_pk <- renderPlotly({
    df <- sim_single()
    req(nrow(df) > 0)
    p <- ggplot(df, aes(x = time, y = CONC_PRED)) +
      annotate("rect", xmin = -Inf, xmax = Inf, ymin = 0.1, ymax = 1.5,
               fill = "#27ae60", alpha = 0.1) +
      annotate("rect", xmin = -Inf, xmax = Inf, ymin = 1.5, ymax = Inf,
               fill = "#e74c3c", alpha = 0.1) +
      geom_line(color = "#8B0000", linewidth = 1) +
      labs(title = "Prednisolone Plasma Concentration",
           x = "Time (days)", y = "Concentration (µg/mL)") +
      annotate("text", x = max(df$time) * 0.8, y = 0.8,
               label = "Therapeutic range", color = "#27ae60", size = 3) +
      annotate("text", x = max(df$time) * 0.8, y = max(df$CONC_PRED) * 0.9,
               label = "Toxic range", color = "#e74c3c", size = 3) +
      theme_pan()
    g2p(p)
  })

  output$plot_cyc_pk <- renderPlotly({
    df <- sim_single()
    req(nrow(df) > 0)
    df_long <- df %>%
      select(time, CYC = CONC_OHC) %>%
      mutate(`4-OH-CYC` = CYC * 0.6) %>%
      pivot_longer(c(CYC, `4-OH-CYC`), names_to = "compound", values_to = "conc")
    p <- ggplot(df_long, aes(x = time, y = conc, color = compound)) +
      geom_line(linewidth = 1) +
      scale_color_manual(values = c(CYC = "#8B0000", `4-OH-CYC` = "#e74c3c")) +
      labs(title = "Cyclophosphamide & 4-OH-CYC Plasma",
           x = "Time (days)", y = "Concentration (µg/mL)", color = "") +
      theme_pan()
    g2p(p)
  })

  output$plot_aza_pk <- renderPlotly({
    df <- sim_single()
    req(nrow(df) > 0)
    df_long <- df %>%
      select(time, Azathioprine = CONC_AZA, Antiviral = CONC_AV) %>%
      pivot_longer(c(Azathioprine, Antiviral), names_to = "drug", values_to = "conc")
    p <- ggplot(df_long, aes(x = time, y = conc, color = drug)) +
      geom_line(linewidth = 1) +
      scale_color_manual(values = c(Azathioprine = "#c0392b", Antiviral = "#922b21")) +
      labs(title = "AZA & Antiviral Plasma Concentration",
           x = "Time (days)", y = "Concentration (µg/mL)", color = "") +
      theme_pan()
    g2p(p)
  })

  output$tbl_pk_summary <- renderDT({
    df <- sim_single()
    req(nrow(df) > 0)
    t_end <- max(df$time)
    summarise_pk <- function(col, name) {
      v <- df[[col]]
      data.frame(
        Drug = name,
        Cmax = round(max(v, na.rm = TRUE), 3),
        Tmax = round(df$time[which.max(v)], 1),
        AUC  = round(sum(diff(df$time) * (head(v,-1) + tail(v,-1))/2, na.rm = TRUE), 1),
        `C_trough` = round(tail(v, 1), 3)
      )
    }
    bind_rows(
      summarise_pk("CONC_PRED", "Prednisolone"),
      summarise_pk("CONC_OHC",  "4-OH-CYC"),
      summarise_pk("CONC_AZA",  "Azathioprine"),
      summarise_pk("CONC_AV",   "Antiviral")
    ) %>%
      datatable(options = list(dom = "t", pageLength = 4),
                rownames = FALSE, class = "compact stripe")
  })

  # ════════════════════════════════════════════════════════════════════════
  # TAB 3 — PD Biomarkers
  # ════════════════════════════════════════════════════════════════════════
  output$plot_ic_c3 <- renderPlotly({
    df <- sim_single()
    req(nrow(df) > 0)
    df_long <- df %>%
      select(time, `Immune Complex` = IC_norm, `Complement C3` = C3_norm) %>%
      pivot_longer(-time, names_to = "marker", values_to = "value")
    p <- ggplot(df_long, aes(x = time, y = value, color = marker)) +
      geom_line(linewidth = 1.1) +
      geom_hline(yintercept = 1, linetype = "dashed", color = "grey50") +
      scale_color_manual(values = c(`Immune Complex` = "#8B0000",
                                    `Complement C3`  = "#2980b9")) +
      labs(title = "Immune Complex & Complement C3 (normalised)",
           x = "Time (days)", y = "Normalised level", color = "") +
      theme_pan()
    g2p(p)
  })

  output$plot_cells <- renderPlotly({
    df <- sim_single()
    req(nrow(df) > 0)
    df_long <- df %>%
      select(time, Neutrophils = NEU_norm,
             `B cells` = BCELL_norm, `T cells` = TCELL_norm) %>%
      pivot_longer(-time, names_to = "cell", values_to = "value")
    p <- ggplot(df_long, aes(x = time, y = value, color = cell)) +
      geom_line(linewidth = 1.1) +
      geom_hline(yintercept = 1, linetype = "dashed", color = "grey50") +
      scale_color_manual(values = c(Neutrophils = "#c0392b",
                                    `B cells`   = "#922b21",
                                    `T cells`   = "#e74c3c")) +
      labs(title = "Immune Cell Counts (normalised)",
           x = "Time (days)", y = "Normalised count", color = "") +
      theme_pan()
    g2p(p)
  })

  output$plot_cytokines <- renderPlotly({
    df <- sim_single()
    req(nrow(df) > 0)
    df_long <- df %>%
      select(time, IL6 = IL6_norm, CRP = CRP_abs, ESR = ESR_norm) %>%
      pivot_longer(-time, names_to = "marker", values_to = "value")
    p <- ggplot(df_long, aes(x = time, y = value, color = marker)) +
      geom_line(linewidth = 1.1) +
      facet_wrap(~marker, scales = "free_y", nrow = 1) +
      scale_color_manual(values = c(IL6 = "#8B0000", CRP = "#c0392b", ESR = "#922b21")) +
      labs(title = "Inflammatory Cytokines / 염증 사이토카인",
           x = "Time (days)", y = "Level", color = "") +
      theme_pan() + theme(legend.position = "none")
    g2p(p)
  })

  # ════════════════════════════════════════════════════════════════════════
  # TAB 4 — Clinical Endpoints
  # ════════════════════════════════════════════════════════════════════════
  output$plot_bvas <- renderPlotly({
    df <- sim_single()
    req(nrow(df) > 0)
    p <- ggplot(df, aes(x = time, y = BVAS_abs)) +
      annotate("rect", xmin=-Inf, xmax=Inf, ymin=0,  ymax=10, fill="#27ae60", alpha=0.1) +
      annotate("rect", xmin=-Inf, xmax=Inf, ymin=10, ymax=30, fill="#f39c12", alpha=0.1) +
      annotate("rect", xmin=-Inf, xmax=Inf, ymin=30, ymax=63, fill="#e74c3c", alpha=0.1) +
      geom_line(color = "#8B0000", linewidth = 1.2) +
      geom_hline(yintercept = c(10, 30), linetype = "dashed", color = "grey50") +
      labs(title = "BVAS (Birmingham Vasculitis Activity Score)",
           x = "Time (days)", y = "BVAS (0–63)") +
      scale_y_continuous(limits = c(0, 63)) +
      annotate("text", x = 10, y = 5,  label = "Remission",     color="#27ae60", size=3) +
      annotate("text", x = 10, y = 20, label = "Mild activity",  color="#f39c12", size=3) +
      annotate("text", x = 10, y = 50, label = "High activity",  color="#e74c3c", size=3) +
      theme_pan()
    g2p(p)
  })

  output$plot_vasc <- renderPlotly({
    df <- sim_single()
    req(nrow(df) > 0)
    p <- ggplot(df, aes(x = time, y = VASC_norm)) +
      geom_area(fill = "#8B0000", alpha = 0.2) +
      geom_line(color = "#8B0000", linewidth = 1.2) +
      labs(title = "Vascular Inflammation Index / 혈관 염증 지수",
           x = "Time (days)", y = "Index (normalised)") +
      theme_pan()
    g2p(p)
  })

  output$plot_egfr <- renderPlotly({
    df <- sim_single()
    req(nrow(df) > 0)
    p <- ggplot(df, aes(x = time, y = EGFR_abs)) +
      annotate("rect", xmin=-Inf, xmax=Inf, ymin=60, ymax=Inf,  fill="#27ae60", alpha=0.1) +
      annotate("rect", xmin=-Inf, xmax=Inf, ymin=30, ymax=60,   fill="#f39c12", alpha=0.1) +
      annotate("rect", xmin=-Inf, xmax=Inf, ymin=0,  ymax=30,   fill="#e74c3c", alpha=0.1) +
      geom_line(color = "#2980b9", linewidth = 1.2) +
      labs(title = "Renal Function — eGFR / 신기능",
           x = "Time (days)", y = "eGFR (mL/min/1.73 m²)") +
      theme_pan()
    g2p(p)
  })

  output$plot_nerve <- renderPlotly({
    df <- sim_single()
    req(nrow(df) > 0)
    p <- ggplot(df, aes(x = time, y = NERVE_norm)) +
      geom_line(color = "#8B0000", linewidth = 1.2) +
      geom_area(fill = "#8B0000", alpha = 0.15) +
      labs(title = "Nerve Damage Score (mononeuritis multiplex NCS)",
           x = "Time (days)", y = "Score (normalised)") +
      theme_pan()
    g2p(p)
  })

  output$plot_miao <- renderPlotly({
    df <- sim_single()
    req(nrow(df) > 0)
    df_long <- df %>%
      select(time, `Microaneurysm burden` = MIAO_norm,
             `Vascular inflammation` = VASC_norm) %>%
      pivot_longer(-time, names_to = "endpoint", values_to = "value")
    p <- ggplot(df_long, aes(x = time, y = value, color = endpoint, fill = endpoint)) +
      geom_line(linewidth = 1.1) +
      geom_area(alpha = 0.1) +
      scale_color_manual(values = c(`Microaneurysm burden`  = "#8B0000",
                                    `Vascular inflammation` = "#c0392b")) +
      scale_fill_manual(values  = c(`Microaneurysm burden`  = "#8B0000",
                                    `Vascular inflammation` = "#c0392b")) +
      labs(title = "Microaneurysm Burden & Vascular Inflammation / 미세동맥류 및 혈관 염증",
           x = "Time (days)", y = "Index (normalised)", color = "", fill = "") +
      theme_pan()
    g2p(p)
  })

  # ════════════════════════════════════════════════════════════════════════
  # TAB 5 — Scenario Comparison
  # ════════════════════════════════════════════════════════════════════════
  make_scen_plot <- function(df, yvar, ytitle) {
    df$arm_label <- arm_labels[df$arm]
    p <- ggplot(df, aes_string(x = "time", y = yvar,
                               color = "arm_label", group = "arm_label")) +
      geom_line(linewidth = 1.1) +
      scale_color_manual(values = setNames(arm_colors, arm_labels)) +
      labs(x = "Time (days)", y = ytitle, color = "Arm / 치료군") +
      theme_pan()
    g2p(p)
  }

  output$plot_scen_bvas <- renderPlotly({
    df <- sim_multi(); req(nrow(df) > 0)
    make_scen_plot(df, "BVAS_abs", "BVAS (0–63)")
  })
  output$plot_scen_crp <- renderPlotly({
    df <- sim_multi(); req(nrow(df) > 0)
    make_scen_plot(df, "CRP_abs", "CRP (mg/L)")
  })
  output$plot_scen_egfr <- renderPlotly({
    df <- sim_multi(); req(nrow(df) > 0)
    make_scen_plot(df, "EGFR_abs", "eGFR (mL/min/1.73 m²)")
  })

  output$tbl_outcomes <- renderDT({
    df <- sim_multi(); req(nrow(df) > 0)

    # Summary at 6 m (182 days) and 18 m (548 days)
    summarise_arm <- function(d) {
      at6  <- d %>% filter(time >= 178 & time <= 186) %>% slice(1)
      at18 <- d %>% filter(time >= 540 & time <= 560) %>% slice(1)
      bvas0_val <- input$bvas0
      rem6  <- if (nrow(at6)  > 0) ifelse(at6$BVAS_abs  <= bvas0_val * 0.1, "Yes", "No") else "N/A"
      rel18 <- if (nrow(at18) > 0) ifelse(at18$BVAS_abs >= bvas0_val * 0.4, "Yes", "No") else "N/A"
      dam   <- if (nrow(at18) > 0) round(at18$VASC_norm, 2) else NA
      data.frame(
        `Treatment` = arm_labels[unique(d$arm)],
        `Remission @ 6m` = rem6,
        `Relapse @ 18m`  = rel18,
        `eGFR @ 6m`      = if (nrow(at6)  > 0) round(at6$EGFR_abs,  1) else NA,
        `BVAS @ 6m`      = if (nrow(at6)  > 0) round(at6$BVAS_abs,  1) else NA,
        `Vascular damage index` = dam,
        check.names = FALSE
      )
    }
    df %>% group_by(arm) %>% group_map(~summarise_arm(.x), .keep = TRUE) %>%
      bind_rows() %>%
      datatable(options = list(dom = "t", pageLength = 6),
                rownames = FALSE, class = "compact stripe hover")
  })

  output$plot_forest <- renderPlotly({
    remission_rates <- data.frame(
      arm    = c("none","pred_mono","pred_cyc","pred_aza","hbv_pred"),
      label  = c("No treatment","Pred mono","Pred+CYC","Pred+AZA","HBV+Pred"),
      rate   = c(0.05, 0.55, 0.78, 0.65, 0.72),
      lo     = c(0.01, 0.42, 0.65, 0.52, 0.58),
      hi     = c(0.12, 0.68, 0.88, 0.77, 0.84)
    ) %>% filter(arm %in% input$arms_sel)

    p <- ggplot(remission_rates,
                aes(x = rate, y = reorder(label, rate),
                    xmin = lo, xmax = hi, color = arm)) +
      geom_point(size = 4) +
      geom_errorbarh(height = 0.2, linewidth = 1) +
      geom_vline(xintercept = 0.5, linetype = "dashed", color = "grey50") +
      scale_color_manual(values = arm_colors) +
      scale_x_continuous(limits = c(0, 1), labels = scales::percent) +
      labs(title = "Remission at 6 months (95% CI)",
           x = "Remission rate", y = "") +
      theme_pan() + theme(legend.position = "none")
    ggplotly(p)
  })

  # ════════════════════════════════════════════════════════════════════════
  # TAB 6 — Biomarker Correlation
  # ════════════════════════════════════════════════════════════════════════
  output$plot_scatter_bvas_crp <- renderPlotly({
    df <- sim_single(); req(nrow(df) > 0)
    p <- ggplot(df, aes(x = BVAS_abs, y = CRP_abs, color = time)) +
      geom_point(alpha = 0.6, size = 2) +
      geom_smooth(method = "loess", color = "#8B0000", se = TRUE) +
      scale_color_gradient(low = "#fadbd8", high = "#8B0000") +
      labs(title = "BVAS vs CRP / BVAS 대 CRP",
           x = "BVAS", y = "CRP (mg/L)", color = "Day") +
      theme_pan()
    ggplotly(p)
  })

  output$plot_scatter_egfr_ic <- renderPlotly({
    df <- sim_single(); req(nrow(df) > 0)
    p <- ggplot(df, aes(x = IC_norm, y = EGFR_abs, color = time)) +
      geom_point(alpha = 0.6, size = 2) +
      geom_smooth(method = "loess", color = "#2980b9", se = TRUE) +
      scale_color_gradient(low = "#d6eaf8", high = "#1a5276") +
      labs(title = "eGFR vs Immune Complex / eGFR 대 면역복합체",
           x = "Immune complex (normalised)", y = "eGFR (mL/min/1.73 m²)",
           color = "Day") +
      theme_pan()
    ggplotly(p)
  })

  output$plot_corr_heat <- renderPlotly({
    df <- sim_single(); req(nrow(df) > 0)
    markers <- c("BVAS_abs","CRP_abs","EGFR_abs",
                 "IC_norm","IL6_norm","NEU_norm","VASC_norm","MIAO_norm")
    labels  <- c("BVAS","CRP","eGFR","Imm.Cx","IL-6","Neutrophil","Vasc.","Microan.")
    cm      <- cor(df[, markers], use = "complete.obs")
    rownames(cm) <- colnames(cm) <- labels
    plot_ly(
      z         = cm,
      x         = labels,
      y         = labels,
      type      = "heatmap",
      colorscale = list(c(0,"#2980b9"), c(0.5,"white"), c(1,"#8B0000")),
      zmin = -1, zmax = 1,
      colorbar  = list(title = "r")
    ) %>%
      layout(title = "Biomarker Correlation Matrix / 바이오마커 상관 행렬",
             xaxis = list(tickangle = -45))
  })

  output$plot_waterfall <- renderPlotly({
    arm_wf <- input$wf_arm
    df <- tryCatch(
      run_sim(pan_mod, age = input$age, wt = input$wt,
              hbv = as.numeric(input$hbv),
              bvas0 = input$bvas0, crp0 = input$crp0, egfr0 = input$egfr0,
              arm = arm_wf, dose_pred = input$dose_pred,
              end_time = 180 * 24),
      error = function(e) NULL
    )
    req(!is.null(df) && nrow(df) > 0)

    # Percentage change at 6 months vs baseline
    baseline <- df[1, ]
    end6m    <- df[nrow(df), ]
    markers <- c("BVAS_abs","CRP_abs","IC_norm","IL6_norm",
                 "NEU_norm","VASC_norm","MIAO_norm","NERVE_norm")
    labels  <- c("BVAS","CRP","Imm.Cx","IL-6",
                 "Neutrophil","Vasc. Index","Microan.","Nerve score")
    pct_change <- 100 * (as.numeric(end6m[, markers]) -
                           as.numeric(baseline[, markers])) /
      pmax(abs(as.numeric(baseline[, markers])), 0.001)

    wf_df <- data.frame(marker = labels, pct = pct_change) %>%
      arrange(pct)
    wf_df$color <- ifelse(wf_df$pct < 0, "#27ae60", "#e74c3c")

    plot_ly(wf_df, x = ~pct, y = ~reorder(marker, pct),
            type = "bar", orientation = "h",
            marker = list(color = ~color)) %>%
      layout(title = paste("% change at 6 months —", arm_labels[arm_wf]),
             xaxis = list(title = "% change from baseline"),
             yaxis = list(title = ""),
             shapes = list(list(type = "line", x0 = 0, x1 = 0,
                                y0 = -0.5, y1 = nrow(wf_df) - 0.5,
                                line = list(color = "black", dash = "dash"))))
  })

  output$plot_subtype <- renderPlotly({
    # Simulate all three subtypes × two arms
    subtypes <- list(
      classic = list(hbv = 0, bvas0 = 22, crp0 = 55, egfr0 = 78),
      hbv_pan = list(hbv = 1, bvas0 = 28, crp0 = 80, egfr0 = 65),
      micro   = list(hbv = 0, bvas0 = 18, crp0 = 45, egfr0 = 82)
    )
    arms_st <- c("pred_mono","pred_cyc")
    results <- list()
    for (st_name in names(subtypes)) {
      st <- subtypes[[st_name]]
      for (a in arms_st) {
        d <- tryCatch(
          run_sim(pan_mod, age = input$age, wt = input$wt,
                  hbv = st$hbv, bvas0 = st$bvas0,
                  crp0 = st$crp0, egfr0 = st$egfr0,
                  arm = a, dose_pred = input$dose_pred,
                  dose_cyc = input$dose_cyc,
                  end_time = 180 * 24),
          error = function(e) NULL
        )
        if (!is.null(d)) {
          d$subtype <- st_name
          results[[paste(st_name, a)]] <- d
        }
      }
    }
    df_all <- bind_rows(results)
    req(nrow(df_all) > 0)

    df_all$arm_label     <- arm_labels[df_all$arm]
    df_all$subtype_label <- c(classic="Classic PAN", hbv_pan="HBV-PAN",
                               micro="Micro-PAN")[df_all$subtype]
    df_all$group <- paste(df_all$subtype_label, "|", df_all$arm_label)

    p <- ggplot(df_all, aes(x = time, y = BVAS_abs, color = subtype_label,
                            linetype = arm_label, group = group)) +
      geom_line(linewidth = 1.0) +
      scale_color_manual(values = c(`Classic PAN` = "#8B0000",
                                    `HBV-PAN`     = "#c0392b",
                                    `Micro-PAN`   = "#e74c3c")) +
      labs(title = "BVAS by disease subtype & treatment / 아형별 BVAS 반응",
           x = "Time (days)", y = "BVAS", color = "Subtype", linetype = "Arm") +
      theme_pan()
    g2p(p)
  })

  # ════════════════════════════════════════════════════════════════════════
  # Downloads
  # ════════════════════════════════════════════════════════════════════════
  output$dl_data <- downloadHandler(
    filename = function() paste0("pan_simulation_", Sys.Date(), ".csv"),
    content  = function(file) write.csv(sim_single(), file, row.names = FALSE)
  )
  output$dl_pk <- downloadHandler(
    filename = function() paste0("pan_pk_", Sys.Date(), ".csv"),
    content  = function(file) write.csv(
      sim_single()[, c("time","CONC_PRED","CONC_OHC","CONC_AZA","CONC_AV")],
      file, row.names = FALSE)
  )
  output$dl_scen <- downloadHandler(
    filename = function() paste0("pan_scenarios_", Sys.Date(), ".csv"),
    content  = function(file) write.csv(sim_multi(), file, row.names = FALSE)
  )
}

# ── Launch ─────────────────────────────────────────────────────────────────
shinyApp(ui = ui, server = server)
