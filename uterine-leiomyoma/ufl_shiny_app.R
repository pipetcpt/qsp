## ============================================================
## Uterine Leiomyoma QSP — Interactive Shiny App
## 자궁근종 정량적 시스템 약리학 인터랙티브 대시보드
##
## 6 Tabs:
##   1. 환자 프로파일 (Patient Profile)
##   2. 약물 PK (Drug Pharmacokinetics)
##   3. PD 주요 지표 (Key PD Markers — Hormones, Fibroid Volume)
##   4. 임상 엔드포인트 (Clinical Endpoints — MBL, Hgb, Symptoms)
##   5. 시나리오 비교 (Treatment Scenario Comparison)
##   6. 바이오마커 패널 (Biomarker Panel)
## ============================================================

library(shiny)
library(mrgsolve)
library(ggplot2)
library(dplyr)
library(tidyr)
library(shinydashboard)
library(DT)
library(plotly)

## ============================================================
## EMBEDDED MRGSOLVE MODEL (compact version for Shiny)
## ============================================================

ufl_model_code <- '
$PARAM
  kGnRH_pulse=0.5, pulse_freq=1.0, kGnRH_deg=2.0,
  kLH_prod=0.08, kLH_deg=0.23, kFSH_prod=0.04, kFSH_deg=0.035,
  EC50_GnRH_LH=0.5, EC50_GnRH_FSH=0.5, nH_GnRH=2.0,
  Emax_GnRH_LH=5.0, Emax_GnRH_FSH=3.0,
  IC50_E2_LH=100.0, IC50_E2_GnRH=120.0, IC50_P4_GnRH=3.0,
  Imax_E2_GnRH=0.7, Imax_P4_GnRH=0.5,
  kE2_base=5.0, kE2_LH_stim=0.8, kE2_deg=0.1,
  kP4_base=0.05, kP4_LH_lut=2.0, kP4_deg=0.2,
  V_fib_0=50.0, V_fib_max=500.0, kgrow_fib=0.003,
  EC50_E2_fib=80.0, Emax_E2_fib=0.8,
  EC50_P4_fib=2.0, Emax_P4_fib=0.6, kfib_apop=0.0005,
  kECM_syn=0.002, kECM_deg=0.001,
  MBL_base=80.0, kMBL_fib=0.5, kMBL_E2=0.3,
  Hgb_0=13.0, kHgb_loss=0.000015, kHgb_prod=0.002, Hgb_setpt=13.0,
  BMD_0=1.0, kBMD_loss=0.0002, kBMD_gain=0.0001, E2_BMD_prot=40.0,
  ka_Leu=0.010, CL_Leu=8.0, Vd_Leu=30.0,
  ka_Ela=1.5, CL_Ela=35.0, Vd_Ela=200.0, F_Ela=0.56, IC50_Ela=5.0,
  ka_Rel=0.5, CL_Rel=14.0, Vd_Rel=1200.0, F_Rel=0.12, IC50_Rel=1.5,
  ka_UPA=1.2, CL_UPA=6.2, Vd_UPA=100.0, F_UPA=0.87, IC50_UPA_PR=0.5, Emax_UPA_PR=0.75,
  Dose_E2addbk=15.0,
  use_GnRHag=0, use_Ela=0, use_Rel=0, use_UPA=0, use_addbk=0,
  PBAC_per_mL=0.8

$CMT GnRH_C, LH_C, FSH_C, E2_C, P4_C, V_fib, ECM_fib, MBL_cum, Hgb_C, BMD_C,
     Leu_depot, Leu_plasma, Ela_gut, Ela_plasma, Rel_gut, Rel_plasma, UPA_gut, UPA_plasma

$INIT
  GnRH_C=0.5, LH_C=8.0, FSH_C=5.0, E2_C=150.0, P4_C=2.0,
  V_fib=50.0, ECM_fib=20.0, MBL_cum=0.0, Hgb_C=13.0, BMD_C=1.0,
  Leu_depot=0, Leu_plasma=0, Ela_gut=0, Ela_plasma=0,
  Rel_gut=0, Rel_plasma=0, UPA_gut=0, UPA_plasma=0

$MAIN
  double C_Ela = (Ela_plasma / Vd_Ela) * 1000.0;
  double C_Rel = (Rel_plasma / Vd_Rel) * 1000.0;
  double C_Leu = (Leu_plasma / Vd_Leu) * 1000.0;
  double C_UPA = (UPA_plasma / Vd_UPA) * 1000.0;
  double E_GnRH_LH = Emax_GnRH_LH * pow(GnRH_C,nH_GnRH) / (pow(EC50_GnRH_LH,nH_GnRH) + pow(GnRH_C,nH_GnRH));
  double E_GnRH_FSH = Emax_GnRH_FSH * pow(GnRH_C,nH_GnRH) / (pow(EC50_GnRH_FSH,nH_GnRH) + pow(GnRH_C,nH_GnRH));
  double I_E2_LH = IC50_E2_LH / (IC50_E2_LH + E2_C);
  double I_E2_pulse = 1.0 - Imax_E2_GnRH * E2_C / (IC50_E2_GnRH + E2_C);
  double I_P4_pulse = 1.0 - Imax_P4_GnRH * P4_C / (IC50_P4_GnRH + P4_C);
  double GnRH_pulserate = pulse_freq * I_E2_pulse * I_P4_pulse;
  double Leu_desens = (use_GnRHag > 0.5) ? C_Leu / (C_Leu + 2.0) : 0.0;
  double Ela_block = (use_Ela > 0.5) ? C_Ela / (C_Ela + IC50_Ela) : 0.0;
  double Rel_block = (use_Rel > 0.5) ? C_Rel / (C_Rel + IC50_Rel) : 0.0;
  double GnRHant_block = 1.0 - (Ela_block > Rel_block ? Ela_block : Rel_block);
  double E2_addbk = (use_addbk > 0.5) ? Dose_E2addbk : 0.0;
  double UPA_PR_block = (use_UPA > 0.5) ? Emax_UPA_PR * C_UPA / (C_UPA + IC50_UPA_PR) : 0.0;
  double E_E2_fib = Emax_E2_fib * E2_C / (EC50_E2_fib + E2_C);
  double E_P4_fib = Emax_P4_fib * P4_C / (EC50_P4_fib + P4_C) * (1.0 - UPA_PR_block);
  double fib_growth_stim = 1.0 + E_E2_fib + E_P4_fib;
  double V_prolif = V_fib - ECM_fib;
  double fib_capacity_factor = log(V_fib_max / (V_fib + 0.001));
  double MBL_rate = (MBL_base + kMBL_fib * V_fib + kMBL_E2 * E2_C) / (24.0 * 28.0);
  double E2_BMD_eff = E2_C / (E2_BMD_prot + E2_C);
  double GnRH_drug_block = (use_Ela > 0.5 || use_Rel > 0.5) ? GnRHant_block : (1.0 - Leu_desens);
  double LH_GnRH_stim = kLH_prod * (1.0 + E_GnRH_LH) * GnRH_drug_block * I_E2_LH;
  double LH_suppression = (use_GnRHag > 0.5) ? (1.0 - Leu_desens) : 1.0;
  double E2_drug_suppress = (use_GnRHag > 0.5) ? Leu_desens : (1.0 - GnRHant_block);
  double P4_drug_suppress = (use_GnRHag > 0.5) ? Leu_desens : (1.0 - GnRHant_block);

$ODE
  dxdt_GnRH_C = kGnRH_pulse * GnRH_pulserate - kGnRH_deg * GnRH_C;
  dxdt_LH_C   = LH_GnRH_stim * LH_suppression - kLH_deg * LH_C;
  dxdt_FSH_C  = kFSH_prod * (1.0 + E_GnRH_FSH) * GnRH_drug_block - kFSH_deg * FSH_C;
  dxdt_E2_C   = (kE2_base + kE2_LH_stim * LH_C) * (1.0 - E2_drug_suppress * 0.9) - kE2_deg * E2_C + E2_addbk;
  dxdt_P4_C   = kP4_base + kP4_LH_lut * LH_C / (LH_C + 5.0) * (1.0 - P4_drug_suppress * 0.85) - kP4_deg * P4_C;
  dxdt_V_fib  = kgrow_fib * V_prolif * fib_capacity_factor * fib_growth_stim - kfib_apop * V_fib - 0.003 * V_fib * UPA_PR_block;
  dxdt_ECM_fib= kECM_syn * V_fib * (1.0 - UPA_PR_block) - kECM_deg * ECM_fib;
  dxdt_MBL_cum= MBL_rate;
  dxdt_Hgb_C  = kHgb_prod * (Hgb_setpt - Hgb_C) * (Hgb_C < Hgb_setpt ? 1.5 : 1.0) - kHgb_loss * MBL_rate * 24.0 * 28.0;
  dxdt_BMD_C  = kBMD_gain * E2_BMD_eff - kBMD_loss * (1.0 - E2_BMD_eff);
  dxdt_Leu_depot  = -ka_Leu * Leu_depot;
  dxdt_Leu_plasma =  ka_Leu * Leu_depot - (CL_Leu / Vd_Leu) * Leu_plasma;
  dxdt_Ela_gut    = -ka_Ela * Ela_gut;
  dxdt_Ela_plasma =  ka_Ela * Ela_gut * F_Ela - (CL_Ela / Vd_Ela) * Ela_plasma;
  dxdt_Rel_gut    = -ka_Rel * Rel_gut;
  dxdt_Rel_plasma =  ka_Rel * Rel_gut * F_Rel - (CL_Rel / Vd_Rel) * Rel_plasma;
  dxdt_UPA_gut    = -ka_UPA * UPA_gut;
  dxdt_UPA_plasma =  ka_UPA * UPA_gut * F_UPA - (CL_UPA / Vd_UPA) * UPA_plasma;

$TABLE
  double Conc_Ela_ngmL = (Ela_plasma / Vd_Ela) * 1000.0;
  double Conc_Rel_ngmL = (Rel_plasma / Vd_Rel) * 1000.0;
  double Conc_Leu_ngmL = (Leu_plasma / Vd_Leu) * 1000.0;
  double Conc_UPA_ngmL = (UPA_plasma / Vd_UPA) * 1000.0;
  double PBAC_score    = MBL_cum * PBAC_per_mL;
  double fib_vol_pct_chg = (V_fib - 50.0) / 50.0 * 100.0;
  double BMD_pct_chg   = (BMD_C - 1.0) * 100.0;
  double hot_flush_score = (E2_C < 30.0) ? (30.0 - E2_C) / 30.0 * 10.0 : 0.0;
  capture Conc_Ela_ngmL, Conc_Rel_ngmL, Conc_Leu_ngmL, Conc_UPA_ngmL;
  capture PBAC_score, fib_vol_pct_chg, BMD_pct_chg, hot_flush_score;
$CAPTURE E2_C, P4_C, LH_C, FSH_C, V_fib, MBL_cum, Hgb_C, BMD_C
'

## Compile model (once at startup)
mod <- mcode("ufl_shiny", ufl_model_code)

## ============================================================
## RUN SIMULATION FUNCTION
## ============================================================

run_sim <- function(scenario, duration_wk = 24, fib_vol_init = 50,
                    hgb_init = 12.0, E2_init = 150.0) {
  tmax     <- duration_wk * 7 * 24
  obs_t    <- seq(0, (duration_wk + 12) * 7 * 24, by = 12)

  # Drug events based on scenario
  ev_doses <- switch(scenario,
    "No Treatment" = ev(time = 0, amt = 0, cmt = 1),
    "Leuprolide 3.75mg" = {
      n_doses <- floor(duration_wk / 4)
      ev(time = seq(0, (n_doses - 1) * 4 * 7 * 24, by = 4 * 7 * 24),
         amt = 3750, cmt = "Leu_depot")
    },
    "Elagolix 150mg BID" = ev(
      time = seq(0, tmax - 12, by = 12), amt = 150, cmt = "Ela_gut"),
    "Elagolix 200mg BID + AB" = ev(
      time = seq(0, tmax - 12, by = 12), amt = 200, cmt = "Ela_gut"),
    "Relugolix + Add-Back" = ev(
      time = seq(0, tmax - 24, by = 24), amt = 40, cmt = "Rel_gut"),
    "UPA 5mg QD" = ev(
      time = seq(0, tmax - 24, by = 24), amt = 5, cmt = "UPA_gut")
  )

  params <- switch(scenario,
    "No Treatment"         = list(use_GnRHag=0, use_Ela=0, use_Rel=0, use_UPA=0, use_addbk=0),
    "Leuprolide 3.75mg"    = list(use_GnRHag=1, use_Ela=0, use_Rel=0, use_UPA=0, use_addbk=0),
    "Elagolix 150mg BID"   = list(use_GnRHag=0, use_Ela=1, use_Rel=0, use_UPA=0, use_addbk=0),
    "Elagolix 200mg BID + AB" = list(use_GnRHag=0, use_Ela=1, use_Rel=0, use_UPA=0, use_addbk=1),
    "Relugolix + Add-Back" = list(use_GnRHag=0, use_Ela=0, use_Rel=1, use_UPA=0, use_addbk=1),
    "UPA 5mg QD"           = list(use_GnRHag=0, use_Ela=0, use_Rel=0, use_UPA=1, use_addbk=0)
  )

  result <- mod %>%
    param(params) %>%
    init(V_fib = fib_vol_init, Hgb_C = hgb_init, E2_C = E2_init) %>%
    mrgsim(ev = ev_doses, obstime = obs_t) %>%
    as.data.frame() %>%
    mutate(time_wk = time / (7 * 24), Scenario = scenario)

  return(result)
}

## ============================================================
## UI
## ============================================================

ui <- dashboardPage(
  skin = "purple",

  dashboardHeader(
    title = HTML("🏥 Uterine Leiomyoma QSP<br><small>자궁근종 QSP 모델</small>"),
    titleWidth = 280
  ),

  dashboardSidebar(
    width = 280,
    sidebarMenu(
      menuItem("환자 프로파일", tabName = "patient",    icon = icon("user-md")),
      menuItem("약물 PK",       tabName = "pk",         icon = icon("pills")),
      menuItem("PD 주요 지표",  tabName = "pd",         icon = icon("chart-line")),
      menuItem("임상 엔드포인트",tabName = "endpoints", icon = icon("heartbeat")),
      menuItem("시나리오 비교", tabName = "scenarios",  icon = icon("balance-scale")),
      menuItem("바이오마커 패널",tabName = "biomarkers", icon = icon("vials"))
    ),

    hr(),
    h5("  치료 설정 (Treatment Settings)", style = "color: #ECF0F1; padding-left: 10px;"),

    selectInput("scenario", "치료 시나리오 선택:",
                choices = c("No Treatment",
                            "Leuprolide 3.75mg",
                            "Elagolix 150mg BID",
                            "Elagolix 200mg BID + AB",
                            "Relugolix + Add-Back",
                            "UPA 5mg QD"),
                selected = "Elagolix 200mg BID + AB"),

    sliderInput("duration_wk", "치료 기간 (주):",
                min = 12, max = 52, value = 24, step = 4),

    sliderInput("fib_vol_init", "초기 근종 용적 (cm³):",
                min = 10, max = 200, value = 50, step = 10),

    sliderInput("hgb_init", "초기 헤모글로빈 (g/dL):",
                min = 8.0, max = 14.0, value = 12.0, step = 0.5),

    sliderInput("E2_init", "초기 E2 수준 (pg/mL):",
                min = 50, max = 300, value = 150, step = 25),

    hr(),
    actionButton("run_btn", "  시뮬레이션 실행", icon = icon("play"),
                 style = "background-color: #8E44AD; color: white; width: 90%; margin: 5px 5%;")
  ),

  dashboardBody(
    tags$head(
      tags$style(HTML("
        .content-wrapper { background-color: #F5F5F5; }
        .box { border-radius: 8px; }
        .info-box { border-radius: 8px; }
        .value-box { border-radius: 8px; }
        h4 { color: #5B2C6F; }
      "))
    ),

    tabItems(
      ## -------- TAB 1: Patient Profile --------
      tabItem(tabName = "patient",
        h3("환자 프로파일 & 질환 개요 (Patient Profile & Disease Overview)"),
        fluidRow(
          valueBoxOutput("vbox_fibvol"),
          valueBoxOutput("vbox_hgb"),
          valueBoxOutput("vbox_E2")
        ),
        fluidRow(
          box(
            title = "자궁근종 발병 위험인자 (Risk Factors)", width = 6,
            solidHeader = TRUE, status = "warning",
            HTML("
              <table class='table table-striped table-sm'>
                <thead><tr><th>위험인자</th><th>상대위험도 (RR)</th></tr></thead>
                <tbody>
                  <tr><td>흑인 여성 (African-American)</td><td>2–3× ↑</td></tr>
                  <tr><td>이른 초경 (Early menarche < 10 yr)</td><td>1.5× ↑</td></tr>
                  <tr><td>미산부 (Nulliparity)</td><td>2× ↑</td></tr>
                  <tr><td>비만 BMI > 30</td><td>1.5× ↑</td></tr>
                  <tr><td>가족력 (Family history)</td><td>2.5× ↑</td></tr>
                  <tr><td>비타민D 결핍</td><td>1.4× ↑</td></tr>
                </tbody>
              </table>
            ")
          ),
          box(
            title = "치료 결정 매트릭스 (Treatment Decision Matrix)", width = 6,
            solidHeader = TRUE, status = "purple",
            HTML("
              <table class='table table-striped table-sm'>
                <thead><tr><th>상황</th><th>권고 치료</th><th>근거</th></tr></thead>
                <tbody>
                  <tr><td>HMB + 임신 원함</td><td>GnRH 길항제+AB / UPA</td><td>ELARIS UF, PEARL</td></tr>
                  <tr><td>HMB + 수술 전 전처치</td><td>GnRH 작용제 또는 길항제</td><td>Friedman 1989</td></tr>
                  <tr><td>HMB + 수술 불원</td><td>레루골릭스+AB (장기)</td><td>LIBERTY 1/2</td></tr>
                  <tr><td>근종 크기 감소만</td><td>GnRH 작용제 (단기)</td><td>Lupron PI</td></tr>
                  <tr><td>증상 없음 (관찰)</td><td>Watchful waiting</td><td>—</td></tr>
                  <tr><td>재발 근종 + 영구해결</td><td>자궁절제술</td><td>—</td></tr>
                </tbody>
              </table>
            ")
          )
        ),
        fluidRow(
          box(
            title = "현재 시뮬레이션 파라미터", width = 12,
            solidHeader = TRUE, status = "primary",
            verbatimTextOutput("params_summary")
          )
        )
      ),

      ## -------- TAB 2: Drug PK --------
      tabItem(tabName = "pk",
        h3("약물 약동학 (Drug Pharmacokinetics)"),
        fluidRow(
          box(
            title = "약물 PK 파라미터 (PK Parameters)", width = 12,
            solidHeader = TRUE, status = "primary",
            tableOutput("pk_table")
          )
        ),
        fluidRow(
          box(
            title = "혈중 약물 농도 (Drug Plasma Concentration)", width = 12,
            solidHeader = TRUE, status = "info",
            plotlyOutput("plot_pk_conc", height = "400px")
          )
        )
      ),

      ## -------- TAB 3: PD Key Markers --------
      tabItem(tabName = "pd",
        h3("PD 주요 지표 — 호르몬 및 근종 (Key PD Markers — Hormones & Fibroid)"),
        fluidRow(
          box(
            title = "혈청 에스트라디올 E2 (pg/mL)", width = 6,
            solidHeader = TRUE, status = "danger",
            plotlyOutput("plot_E2", height = "350px")
          ),
          box(
            title = "근종 용적 변화 (Fibroid Volume, cm³)", width = 6,
            solidHeader = TRUE, status = "warning",
            plotlyOutput("plot_fibvol", height = "350px")
          )
        ),
        fluidRow(
          box(
            title = "LH / FSH (IU/L)", width = 6,
            solidHeader = TRUE, status = "info",
            plotlyOutput("plot_gonadotropins", height = "350px")
          ),
          box(
            title = "프로게스테론 P4 (ng/mL)", width = 6,
            solidHeader = TRUE, status = "success",
            plotlyOutput("plot_P4", height = "350px")
          )
        )
      ),

      ## -------- TAB 4: Clinical Endpoints --------
      tabItem(tabName = "endpoints",
        h3("임상 엔드포인트 (Clinical Endpoints)"),
        fluidRow(
          valueBoxOutput("vbox_MBL_final"),
          valueBoxOutput("vbox_Hgb_final"),
          valueBoxOutput("vbox_BMD_final")
        ),
        fluidRow(
          box(
            title = "월경혈량 (MBL, mL/cycle)", width = 6,
            solidHeader = TRUE, status = "danger",
            plotlyOutput("plot_MBL", height = "350px")
          ),
          box(
            title = "헤모글로빈 Hgb (g/dL)", width = 6,
            solidHeader = TRUE, status = "warning",
            plotlyOutput("plot_Hgb", height = "350px")
          )
        ),
        fluidRow(
          box(
            title = "골밀도 BMD 변화율 (%)", width = 6,
            solidHeader = TRUE, status = "info",
            plotlyOutput("plot_BMD", height = "350px")
          ),
          box(
            title = "안면홍조 점수 (Hot Flush Score, 0–10)", width = 6,
            solidHeader = TRUE, status = "success",
            plotlyOutput("plot_hotflush", height = "350px")
          )
        )
      ),

      ## -------- TAB 5: Scenario Comparison --------
      tabItem(tabName = "scenarios",
        h3("치료 시나리오 비교 (Treatment Scenario Comparison)"),
        fluidRow(
          box(
            title = "모든 시나리오 실행 (버튼 클릭)", width = 12,
            solidHeader = TRUE, status = "primary",
            actionButton("run_all_btn", "  전체 시나리오 비교 실행",
                         icon = icon("play-circle"),
                         style = "background-color: #2980B9; color: white; font-size: 14px; padding: 10px 20px;"),
            br(), br(),
            p("6가지 치료 시나리오를 동시에 시뮬레이션하여 비교합니다.")
          )
        ),
        fluidRow(
          box(
            title = "근종 용적 비교 (All Scenarios)", width = 6,
            solidHeader = TRUE, status = "warning",
            plotlyOutput("plot_scen_fibvol", height = "350px")
          ),
          box(
            title = "E2 비교 (All Scenarios)", width = 6,
            solidHeader = TRUE, status = "danger",
            plotlyOutput("plot_scen_E2", height = "350px")
          )
        ),
        fluidRow(
          box(
            title = "월경혈량 비교 (All Scenarios)", width = 6,
            solidHeader = TRUE, status = "info",
            plotlyOutput("plot_scen_MBL", height = "350px")
          ),
          box(
            title = "골밀도 비교 (All Scenarios)", width = 6,
            solidHeader = TRUE, status = "success",
            plotlyOutput("plot_scen_BMD", height = "350px")
          )
        ),
        fluidRow(
          box(
            title = "요약 성과 테이블 (Summary Outcomes at Week 24)", width = 12,
            solidHeader = TRUE, status = "primary",
            DTOutput("summary_table")
          )
        )
      ),

      ## -------- TAB 6: Biomarker Panel --------
      tabItem(tabName = "biomarkers",
        h3("바이오마커 패널 (Biomarker Panel)"),
        fluidRow(
          box(
            title = "주요 바이오마커 (Key Biomarkers)", width = 12,
            solidHeader = TRUE, status = "purple",
            HTML("
              <table class='table table-bordered table-striped'>
                <thead style='background:#8E44AD; color:white'>
                  <tr>
                    <th>바이오마커</th><th>정상 범위</th>
                    <th>치료 목표</th><th>임상 의미</th>
                  </tr>
                </thead>
                <tbody>
                  <tr>
                    <td><b>E2 (에스트라디올)</b></td>
                    <td>20–150 pg/mL (가임기)</td>
                    <td>GnRH 길항제: 20–50 pg/mL</td>
                    <td>근종 성장 주요 드라이버</td>
                  </tr>
                  <tr>
                    <td><b>P4 (프로게스테론)</b></td>
                    <td>황체기: 6–20 ng/mL</td>
                    <td>GnRH 길항제: &lt; 1 ng/mL</td>
                    <td>ECM 축적 촉진</td>
                  </tr>
                  <tr>
                    <td><b>LH (황체화호르몬)</b></td>
                    <td>2–15 IU/L</td>
                    <td>치료 시 &lt; 2 IU/L</td>
                    <td>GnRH-R 기능 반영</td>
                  </tr>
                  <tr>
                    <td><b>FSH (난포자극호르몬)</b></td>
                    <td>3–10 IU/L</td>
                    <td>치료 시 &lt; 3 IU/L</td>
                    <td>난소 기능 반영</td>
                  </tr>
                  <tr>
                    <td><b>Hgb (헤모글로빈)</b></td>
                    <td>≥ 12.0 g/dL (여성)</td>
                    <td>치료 목표: ≥ 12.0 g/dL</td>
                    <td>AUB로 인한 빈혈</td>
                  </tr>
                  <tr>
                    <td><b>BMD (골밀도)</b></td>
                    <td>T-score ≥ -1.0 (정상)</td>
                    <td>치료 중 &lt; -3% 변화</td>
                    <td>저에스트로겐 부작용</td>
                  </tr>
                </tbody>
              </table>
            ")
          )
        ),
        fluidRow(
          box(
            title = "임상시험 주요 결과 (Key Clinical Trial Results)", width = 12,
            solidHeader = TRUE, status = "primary",
            HTML("
              <table class='table table-bordered table-striped table-sm'>
                <thead style='background:#2C3E50; color:white'>
                  <tr>
                    <th>임상시험</th><th>약물</th><th>용량</th>
                    <th>1차 엔드포인트</th><th>결과</th><th>참고문헌</th>
                  </tr>
                </thead>
                <tbody>
                  <tr>
                    <td>ELARIS UF-I</td><td>Elagolix</td><td>200mg BID+AB</td>
                    <td>Heavy MBL reduction (Week 24)</td>
                    <td>68.5% achieved (&lt;80mL & ≥50% reduction)</td>
                    <td>Simon JA et al. NEJM 2020;382:328</td>
                  </tr>
                  <tr>
                    <td>ELARIS UF-II</td><td>Elagolix</td><td>200mg BID+AB</td>
                    <td>Heavy MBL reduction (Week 24)</td>
                    <td>76.5% achieved</td>
                    <td>Schlaff WD et al. NEJM 2020;382:317</td>
                  </tr>
                  <tr>
                    <td>LIBERTY 1</td><td>Relugolix comb.</td><td>40mg QD+AB</td>
                    <td>HMB reduction (Week 24)</td>
                    <td>71.2% achieved</td>
                    <td>Lukes AS et al. NEJM 2021;384:630</td>
                  </tr>
                  <tr>
                    <td>LIBERTY 2</td><td>Relugolix comb.</td><td>40mg QD+AB</td>
                    <td>HMB reduction (Week 24)</td>
                    <td>70.6% achieved</td>
                    <td>Al-Hendy A et al. NEJM 2021;384:630</td>
                  </tr>
                  <tr>
                    <td>PRIMROSE 1</td><td>Linzagolix</td><td>200mg QD+AB</td>
                    <td>HMB reduction (Week 24)</td>
                    <td>93.9% achieved</td>
                    <td>Murji A et al. NEJM 2022;387:1767</td>
                  </tr>
                  <tr>
                    <td>PEARL I</td><td>UPA</td><td>5mg QD × 13wk</td>
                    <td>Controlled bleeding at end of treatment</td>
                    <td>91% controlled bleeding</td>
                    <td>Donnez J et al. NEJM 2012;366:409</td>
                  </tr>
                  <tr>
                    <td>Friedman 1989</td><td>Leuprolide</td><td>3.75mg depot q4w</td>
                    <td>Fibroid volume reduction</td>
                    <td>35–50% by 12 weeks</td>
                    <td>Friedman AJ. Fertil Steril 1989;51:61</td>
                  </tr>
                </tbody>
              </table>
            ")
          )
        ),
        fluidRow(
          box(
            title = "PBAC 점수 시뮬레이션 (PBAC Score)", width = 12,
            solidHeader = TRUE, status = "info",
            plotlyOutput("plot_PBAC", height = "350px")
          )
        )
      )
    )
  )
)

## ============================================================
## SERVER
## ============================================================

server <- function(input, output, session) {

  ## Reactive simulation
  sim_data <- eventReactive(input$run_btn, {
    withProgress(message = "시뮬레이션 실행 중...", value = 0.5, {
      run_sim(
        scenario    = input$scenario,
        duration_wk = input$duration_wk,
        fib_vol_init= input$fib_vol_init,
        hgb_init    = input$hgb_init,
        E2_init     = input$E2_init
      )
    })
  }, ignoreNULL = FALSE)

  ## All scenarios comparison
  all_scen_data <- eventReactive(input$run_all_btn, {
    scenarios <- c("No Treatment", "Leuprolide 3.75mg", "Elagolix 150mg BID",
                   "Elagolix 200mg BID + AB", "Relugolix + Add-Back", "UPA 5mg QD")
    withProgress(message = "전체 시나리오 실행 중...", value = 0, {
      result <- lapply(seq_along(scenarios), function(i) {
        incProgress(1 / length(scenarios), detail = scenarios[i])
        run_sim(scenarios[i],
                duration_wk  = 24,
                fib_vol_init = input$fib_vol_init,
                hgb_init     = input$hgb_init,
                E2_init      = input$E2_init)
      })
      bind_rows(result)
    })
  })

  ## ---- Value Boxes ----
  output$vbox_fibvol <- renderValueBox({
    valueBox(
      paste0(input$fib_vol_init, " cm³"),
      "초기 근종 용적", icon = icon("circle"),
      color = "orange"
    )
  })
  output$vbox_hgb <- renderValueBox({
    valueBox(
      paste0(input$hgb_init, " g/dL"),
      "초기 헤모글로빈", icon = icon("tint"),
      color = if (input$hgb_init < 12) "red" else "green"
    )
  })
  output$vbox_E2 <- renderValueBox({
    valueBox(
      paste0(input$E2_init, " pg/mL"),
      "초기 에스트라디올 (E2)", icon = icon("female"),
      color = "purple"
    )
  })

  ## ---- Params summary ----
  output$params_summary <- renderText({
    paste0(
      "선택 시나리오: ", input$scenario, "\n",
      "치료 기간: ", input$duration_wk, " 주\n",
      "초기 근종 용적: ", input$fib_vol_init, " cm³\n",
      "초기 헤모글로빈: ", input$hgb_init, " g/dL\n",
      "초기 E2: ", input$E2_init, " pg/mL\n",
      "오늘 날짜: 2026-06-25"
    )
  })

  ## ---- PK Table ----
  output$pk_table <- renderTable({
    data.frame(
      Drug = c("Leuprolide (Depot)", "Elagolix", "Relugolix", "UPA"),
      Dose = c("3.75 mg Q4W (depot)", "150/200 mg BID", "40 mg QD", "5 mg QD"),
      t_half = c("~3–4 weeks (depot release)", "4–6 h", "~60 h", "32–38 h"),
      Bioavailability = c("~95%", "56%", "12%", "87%"),
      Tmax = c("~3–4 weeks (plateau)", "~1 h", "~2 h", "~1 h"),
      CL = c("8 L/h", "35 L/h", "14 L/h", "6.2 L/h"),
      stringsAsFactors = FALSE
    )
  })

  ## ---- PK Concentration Plot ----
  output$plot_pk_conc <- renderPlotly({
    df <- sim_data()
    # Find which drug is active
    drug_col <- switch(input$scenario,
      "Leuprolide 3.75mg"     = "Conc_Leu_ngmL",
      "Elagolix 150mg BID"    = "Conc_Ela_ngmL",
      "Elagolix 200mg BID + AB" = "Conc_Ela_ngmL",
      "Relugolix + Add-Back"  = "Conc_Rel_ngmL",
      "UPA 5mg QD"            = "Conc_UPA_ngmL",
      NULL
    )
    if (is.null(drug_col) || !(drug_col %in% names(df))) {
      return(plotly_empty() %>% layout(title = "선택된 약물 없음 (No Treatment)"))
    }
    p <- ggplot(df, aes_string(x = "time_wk", y = drug_col)) +
      geom_line(color = "#8E44AD", size = 1) +
      labs(title = paste(input$scenario, "— 혈중 농도 (ng/mL)"),
           x = "Time (weeks)", y = "Drug Concentration (ng/mL)") +
      theme_bw()
    ggplotly(p)
  })

  ## ---- PD Plots ----
  make_plot <- function(df, y_col, title_str, ylab, hline_vals = NULL,
                        hline_cols = NULL, hline_labels = NULL, color = "#E74C3C") {
    p <- ggplot(df, aes_string(x = "time_wk", y = y_col)) +
      geom_line(color = color, size = 1.2) +
      labs(title = title_str, x = "Time (weeks)", y = ylab) +
      theme_bw()
    if (!is.null(hline_vals)) {
      for (i in seq_along(hline_vals)) {
        p <- p + geom_hline(yintercept = hline_vals[i], linetype = "dashed",
                             color = hline_cols[i], alpha = 0.8)
      }
    }
    ggplotly(p)
  }

  output$plot_E2 <- renderPlotly({
    make_plot(sim_data(), "E2_C",
              "혈청 에스트라디올 E2",
              "E2 (pg/mL)",
              hline_vals = c(20, 50, 150),
              hline_cols = c("red", "orange", "blue"))
  })
  output$plot_fibvol <- renderPlotly({
    make_plot(sim_data(), "V_fib",
              "근종 용적 (Fibroid Volume)",
              "Volume (cm³)", color = "#E67E22")
  })
  output$plot_P4 <- renderPlotly({
    make_plot(sim_data(), "P4_C",
              "혈청 프로게스테론 P4",
              "P4 (ng/mL)", color = "#27AE60")
  })
  output$plot_gonadotropins <- renderPlotly({
    df <- sim_data() %>%
      select(time_wk, LH_C, FSH_C) %>%
      pivot_longer(c(LH_C, FSH_C), names_to = "Hormone", values_to = "Level")
    p <- ggplot(df, aes(x = time_wk, y = Level, color = Hormone)) +
      geom_line(size = 1.2) +
      labs(title = "LH / FSH (IU/L)", x = "Time (weeks)", y = "IU/L") +
      theme_bw()
    ggplotly(p)
  })

  ## ---- Endpoint Plots ----
  output$vbox_MBL_final <- renderValueBox({
    df <- sim_data()
    last_MBL <- tail(df$MBL_cum, 1)
    valueBox(
      sprintf("%.0f mL", last_MBL),
      "최종 MBL (mL/cycle)", icon = icon("tint"),
      color = if (last_MBL > 80) "red" else "green"
    )
  })
  output$vbox_Hgb_final <- renderValueBox({
    df <- sim_data()
    last_Hgb <- tail(df$Hgb_C, 1)
    valueBox(
      sprintf("%.1f g/dL", last_Hgb),
      "최종 헤모글로빈", icon = icon("tint"),
      color = if (last_Hgb < 12) "red" else "green"
    )
  })
  output$vbox_BMD_final <- renderValueBox({
    df <- sim_data()
    last_BMD <- tail(df$BMD_pct_chg, 1)
    valueBox(
      sprintf("%.1f%%", last_BMD),
      "BMD 변화율 (%)", icon = icon("bone"),
      color = if (last_BMD < -3) "red" else "green"
    )
  })

  output$plot_MBL <- renderPlotly({
    make_plot(sim_data(), "MBL_cum",
              "월경혈량 MBL",
              "MBL (mL/cycle)",
              hline_vals = c(80), hline_cols = c("red"), color = "#C0392B")
  })
  output$plot_Hgb <- renderPlotly({
    make_plot(sim_data(), "Hgb_C",
              "헤모글로빈 Hgb",
              "Hgb (g/dL)",
              hline_vals = c(12.0, 11.0),
              hline_cols = c("orange", "red"), color = "#884EA0")
  })
  output$plot_BMD <- renderPlotly({
    make_plot(sim_data(), "BMD_pct_chg",
              "골밀도 변화율 BMD (%)",
              "BMD Change (%)",
              hline_vals = c(-3, 0),
              hline_cols = c("red", "grey"), color = "#2E86C1")
  })
  output$plot_hotflush <- renderPlotly({
    make_plot(sim_data(), "hot_flush_score",
              "안면홍조 점수 (Hot Flush Score)",
              "Score (0–10)", color = "#E74C3C")
  })
  output$plot_PBAC <- renderPlotly({
    make_plot(sim_data(), "PBAC_score",
              "PBAC 점수 (Pictorial Blood Assessment)",
              "PBAC Score",
              hline_vals = c(64), hline_cols = c("red"), color = "#8E44AD")
  })

  ## ---- Scenario Comparison Plots ----
  make_scen_plot <- function(df, y_col, title_str, ylab) {
    p <- ggplot(df, aes_string(x = "time_wk", y = y_col, color = "Scenario")) +
      geom_line(size = 1) +
      labs(title = title_str, x = "Time (weeks)", y = ylab) +
      theme_bw() +
      theme(legend.position = "bottom",
            legend.text = element_text(size = 8))
    ggplotly(p) %>% layout(legend = list(orientation = "h", y = -0.3))
  }

  output$plot_scen_fibvol <- renderPlotly({
    req(all_scen_data())
    make_scen_plot(all_scen_data(), "V_fib", "근종 용적 비교", "Volume (cm³)")
  })
  output$plot_scen_E2 <- renderPlotly({
    req(all_scen_data())
    make_scen_plot(all_scen_data(), "E2_C", "E2 비교 (pg/mL)", "E2 (pg/mL)")
  })
  output$plot_scen_MBL <- renderPlotly({
    req(all_scen_data())
    make_scen_plot(all_scen_data(), "MBL_cum", "월경혈량 비교 (mL/cycle)", "MBL (mL)")
  })
  output$plot_scen_BMD <- renderPlotly({
    req(all_scen_data())
    make_scen_plot(all_scen_data(), "BMD_pct_chg", "골밀도 변화율 비교 (%)", "BMD Change (%)")
  })

  output$summary_table <- renderDT({
    req(all_scen_data())
    df <- all_scen_data()
    wk24 <- df %>%
      filter(abs(time_wk - 24) < 0.5) %>%
      group_by(Scenario) %>%
      slice_head(n = 1) %>%
      ungroup() %>%
      mutate(
        `근종 용적 (cm³)` = round(V_fib, 1),
        `E2 (pg/mL)` = round(E2_C, 1),
        `MBL (mL/cycle)` = round(MBL_cum, 1),
        `Hgb (g/dL)` = round(Hgb_C, 2),
        `BMD 변화 (%)` = round(BMD_pct_chg, 2),
        `HMB 해소` = ifelse(MBL_cum < 80, "✓ 달성", "✗ 미달")
      ) %>%
      select(Scenario, `근종 용적 (cm³)`, `E2 (pg/mL)`,
             `MBL (mL/cycle)`, `Hgb (g/dL)`, `BMD 변화 (%)`, `HMB 해소`)

    datatable(wk24, rownames = FALSE,
              options = list(dom = "t", pageLength = 10),
              class = "table table-striped table-bordered")
  })
}

## ============================================================
## RUN APP
## ============================================================
shinyApp(ui = ui, server = server)
