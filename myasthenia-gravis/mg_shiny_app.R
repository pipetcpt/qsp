## ============================================================
## Myasthenia Gravis QSP — Interactive Shiny Dashboard
## 중증 근무력증 정량적 시스템 약리학 인터랙티브 대시보드
##
## Tabs:
##   1. 환자 프로파일 & 질환 개요 (Patient Profile & Disease Overview)
##   2. 약동학 (PK — Pyridostigmine, Prednisolone, Eculizumab, Efgartigimod)
##   3. PD 주요 지표 (AChR density, AChE inhibition, Complement, NMJ Safety Factor)
##   4. 임상 엔드포인트 (QMG, MG-ADL, MGFA)
##   5. 치료 시나리오 비교 (Scenario Comparison)
##   6. 바이오마커 패널 (Biomarker Panel)
## ============================================================

library(shiny)
library(shinydashboard)
library(mrgsolve)
library(dplyr)
library(ggplot2)
library(tidyr)
library(DT)
library(plotly)

## ============================================================
## mrgsolve model definition (inline)
## ============================================================
mg_code <- '
$PARAM @annotated
  KA_PYR   : 1.20  : Pyridostigmine absorption (1/h)
  F_PYR    : 0.18  : Pyridostigmine bioavailability
  VD_PYR   : 90.0  : Pyridostigmine Vd (L)
  KE_PYR   : 0.37  : Pyridostigmine ke (1/h)
  EMAX_PYR : 0.90  : Max AChE inhibition
  EC50_PYR : 50.0  : EC50 pyridostigmine (ng/mL)
  GAMMA_PYR: 1.50  : Hill coefficient
  KA_PRED  : 1.50  : Prednisolone absorption (1/h)
  F_PRED   : 0.82  : Prednisolone bioavailability
  VD_PRED  : 35.0  : Prednisolone Vd (L)
  CL_PRED  : 21.0  : Prednisolone CL (L/h)
  Q_PRED   : 5.0   : Pred intercompartmental CL
  VP_PRED  : 50.0  : Pred peripheral Vd
  EMAX_GR  : 0.90  : Max GR immunosuppression
  EC50_GR  : 30.0  : EC50 prednisolone (ng/mL)
  KA_AZA   : 1.00  : AZA absorption (1/h)
  F_AZA    : 0.50  : AZA bioavailability
  K_AZA2MP : 0.40  : AZA → 6-MP rate (1/h)
  K_MP2TGN : 0.08  : 6-MP → 6-TGN rate (1/h)
  KE_TGN   : 0.005 : 6-TGN elimination (1/h)
  VD_TGN   : 120.0 : 6-TGN Vd (L)
  EMAX_TGN : 0.80  : Max 6-TGN effect
  EC50_TGN : 100.0 : EC50 6-TGN
  CL_ECUL  : 0.31  : Eculizumab CL (L/h)
  VD_ECUL  : 7.7   : Eculizumab Vd (L)
  Q_ECUL   : 0.5   : Ecul intercompartmental CL
  VP_ECUL  : 3.0   : Ecul peripheral Vd
  EMAX_C5  : 0.99  : Max C5 inhibition
  EC50_C5  : 100.0 : EC50 eculizumab (ng/mL)
  CL_EFGAR : 0.15  : Efgartigimod CL (L/h)
  VD_EFGAR : 6.0   : Efgartigimod Vd (L)
  KON_FCRN : 0.002 : Efgar FcRn on-rate
  KOFF_FCRN: 0.05  : Efgar FcRn off-rate
  FCRN_TOT : 50.0  : Total FcRn (nmol/L)
  K_TFH_IN : 0.002 : Tfh input rate
  K_TFH_D  : 0.01  : Tfh death rate (1/h)
  K_GCB_IN : 0.50  : GC-B proliferation rate
  K_GCB_D  : 0.05  : GC-B death rate
  K_SLPC_F : 0.30  : SL-PC formation rate
  K_SLPC_D : 0.04  : SL-PC death rate
  K_LLPC_F : 0.10  : LL-PC formation rate
  K_LLPC_D : 0.0003: LL-PC death rate
  K_AB_SYN : 5.0   : Ab synthesis rate
  K_AB_DEG : 0.0033: IgG degradation (1/h)
  K_AB_FCRN: 0.0020: FcRn rescue rate (1/h)
  AB0      : 2.5   : Baseline AChR-Ab (nmol/L)
  ACHR0    : 1.0   : Normal AChR density
  K_ACHR_IN: 0.01  : AChR synthesis (1/h)
  K_ACHR_D : 0.01  : AChR degradation (1/h)
  K_ACHR_AB: 0.005 : Ab-mediated downregulation
  COMP0    : 1.0   : Baseline complement
  K_COMP_AB: 0.10  : Complement activation
  K_COMP_D : 0.30  : Complement decay (1/h)
  SF0      : 1.0   : Normal NMJ safety factor
  K_SF_ACHR: 0.80  : ACHR weight in SF
  K_SF_ACH : 0.20  : ACh weight in SF
  QMG0     : 30.0  : Baseline QMG score
  QMG_MAX  : 39.0  : Maximum QMG
  K_QMG_SF : 0.70  : SF weight in QMG
  K_QMG_AB : 0.30  : Ab weight in QMG

$CMT @annotated
  GUT_PYR   : Pyridostigmine gut
  CENT_PYR  : Pyridostigmine central
  GUT_PRED  : Prednisolone gut
  CENT_PRED : Prednisolone central
  PERIPH_PRED: Prednisolone peripheral
  GUT_AZA   : Azathioprine gut
  SIXMP     : 6-Mercaptopurine
  SIXTGN    : 6-TGN
  CENT_ECUL : Eculizumab central
  PERIPH_ECUL: Eculizumab peripheral
  CENT_EFGAR: Efgartigimod central
  FCRN_EFGAR: FcRn-bound Efgartigimod
  TFH       : Tfh cells
  GCB       : GC-B cells
  SLPC      : Short-lived plasma cells
  LLPC      : Long-lived plasma cells
  ACHR_AB_C : AChR antibody (nmol/L)
  ACHR_DEN  : AChR density (normalized)
  COMP_ACT  : Complement activity
  NMJ_SF    : NMJ safety factor

$MAIN
  double C_PYR   = CENT_PYR / VD_PYR;
  double C_PRED  = CENT_PRED / VD_PRED;
  double C_TGN   = SIXTGN / VD_TGN;
  double C_ECUL  = CENT_ECUL / VD_ECUL;
  double C_EFGAR = CENT_EFGAR / VD_EFGAR;
  double INH_ACHE = EMAX_PYR * pow(C_PYR, GAMMA_PYR) /
                    (pow(EC50_PYR, GAMMA_PYR) + pow(C_PYR, GAMMA_PYR));
  double ACH_AVAIL = 1.0 + INH_ACHE;
  double EFF_GR  = EMAX_GR * C_PRED / (EC50_GR + C_PRED);
  double EFF_TGN = EMAX_TGN * C_TGN / (EC50_TGN + C_TGN);
  double IMM_SUPP = 1.0 - (1.0 - EFF_GR) * (1.0 - EFF_TGN);
  double INH_C5  = EMAX_C5 * C_ECUL / (EC50_C5 + C_ECUL);
  double FCRN_FREE = FCRN_TOT - FCRN_EFGAR;
  if(FCRN_FREE < 0) FCRN_FREE = 0;
  double FCRN_OCC = FCRN_EFGAR / (FCRN_TOT + 1e-6);
  double K_AB_EFF = K_AB_DEG * (1.0 + 3.0 * FCRN_OCC);
  TFH_0 = 1.0; GCB_0 = 1.0; SLPC_0 = 1.0; LLPC_0 = 1.0;
  ACHR_AB_C_0 = AB0; ACHR_DEN_0 = ACHR0; NMJ_SF_0 = SF0; COMP_ACT_0 = COMP0;

$ODE
  double C_PYR_od   = CENT_PYR / VD_PYR;
  double INH_ACHE_od = EMAX_PYR * pow(C_PYR_od, GAMMA_PYR) /
                       (pow(EC50_PYR, GAMMA_PYR) + pow(C_PYR_od, GAMMA_PYR));
  double ACH_AVAIL_od = 1.0 + INH_ACHE_od;
  double C_PRED_od = CENT_PRED / VD_PRED;
  double C_TGN_od  = SIXTGN / VD_TGN;
  double EFF_GR_od  = EMAX_GR * C_PRED_od / (EC50_GR + C_PRED_od);
  double EFF_TGN_od = EMAX_TGN * C_TGN_od / (EC50_TGN + C_TGN_od);
  double IMM_SUPP_od = 1.0 - (1.0 - EFF_GR_od) * (1.0 - EFF_TGN_od);
  double C_ECUL_od = CENT_ECUL / VD_ECUL;
  double INH_C5_od = EMAX_C5 * C_ECUL_od / (EC50_C5 + C_ECUL_od);
  double FCRN_FREE_od = FCRN_TOT - FCRN_EFGAR;
  if(FCRN_FREE_od < 0) FCRN_FREE_od = 0;
  double FCRN_OCC_od = FCRN_EFGAR / (FCRN_TOT + 1e-6);
  double K_AB_EFF_od = K_AB_DEG * (1.0 + 3.0 * FCRN_OCC_od);
  double C_EFGAR_od = CENT_EFGAR / VD_EFGAR;
  dxdt_GUT_PYR    = -KA_PYR * GUT_PYR;
  dxdt_CENT_PYR   = KA_PYR * GUT_PYR * F_PYR - KE_PYR * CENT_PYR;
  dxdt_GUT_PRED   = -KA_PRED * GUT_PRED;
  dxdt_CENT_PRED  = KA_PRED * GUT_PRED * F_PRED
                    - (CL_PRED / VD_PRED + Q_PRED / VD_PRED) * CENT_PRED
                    + Q_PRED / VP_PRED * PERIPH_PRED;
  dxdt_PERIPH_PRED = Q_PRED / VD_PRED * CENT_PRED - Q_PRED / VP_PRED * PERIPH_PRED;
  dxdt_GUT_AZA    = -KA_AZA * GUT_AZA;
  dxdt_SIXMP      = KA_AZA * GUT_AZA * F_AZA - (K_AZA2MP + K_MP2TGN) * SIXMP;
  dxdt_SIXTGN     = K_MP2TGN * SIXMP - KE_TGN * SIXTGN;
  dxdt_CENT_ECUL  = -(CL_ECUL / VD_ECUL + Q_ECUL / VD_ECUL) * CENT_ECUL
                    + Q_ECUL / VP_ECUL * PERIPH_ECUL;
  dxdt_PERIPH_ECUL = Q_ECUL / VD_ECUL * CENT_ECUL - Q_ECUL / VP_ECUL * PERIPH_ECUL;
  dxdt_CENT_EFGAR = -(CL_EFGAR / VD_EFGAR) * CENT_EFGAR
                    - KON_FCRN * C_EFGAR_od * FCRN_FREE_od
                    + KOFF_FCRN * FCRN_EFGAR;
  dxdt_FCRN_EFGAR = KON_FCRN * C_EFGAR_od * FCRN_FREE_od
                    - KOFF_FCRN * FCRN_EFGAR
                    - (CL_EFGAR / VD_EFGAR) * FCRN_EFGAR;
  dxdt_TFH  = K_TFH_IN * (1.0 - IMM_SUPP_od) - K_TFH_D * TFH;
  dxdt_GCB  = K_GCB_IN * TFH * GCB * (1.0 - IMM_SUPP_od) - K_GCB_D * GCB;
  dxdt_SLPC = K_SLPC_F * GCB * (1.0 - IMM_SUPP_od) - K_SLPC_D * SLPC;
  dxdt_LLPC = K_LLPC_F * GCB * (1.0 - IMM_SUPP_od * 0.5) - K_LLPC_D * LLPC;
  dxdt_ACHR_AB_C = K_AB_SYN * (SLPC + LLPC) - K_AB_EFF_od * ACHR_AB_C;
  dxdt_ACHR_DEN  = K_ACHR_IN * ACHR0 - K_ACHR_D * ACHR_DEN
                  - K_ACHR_AB * ACHR_AB_C * ACHR_DEN;
  dxdt_COMP_ACT  = K_COMP_AB * ACHR_AB_C
                  - K_COMP_D * COMP_ACT
                  - K_COMP_D * INH_C5_od * COMP_ACT;
  double SF_TARG_od = K_SF_ACHR * (ACHR_DEN / ACHR0) + K_SF_ACH * ACH_AVAIL_od;
  dxdt_NMJ_SF    = 0.1 * (SF_TARG_od - NMJ_SF);

$TABLE
  double C_PYR_OUT   = CENT_PYR / VD_PYR;
  double C_PRED_OUT  = CENT_PRED / VD_PRED;
  double C_ECUL_OUT  = CENT_ECUL / VD_ECUL;
  double C_EFGAR_OUT = CENT_EFGAR / VD_EFGAR;
  double AChE_INH_PCT = 100.0 * EMAX_PYR * pow(C_PYR_OUT, GAMMA_PYR) /
                        (pow(EC50_PYR, GAMMA_PYR) + pow(C_PYR_OUT, GAMMA_PYR));
  double C5_INH_PCT  = 100.0 * EMAX_C5 * C_ECUL_OUT / (EC50_C5 + C_ECUL_OUT);
  double FCRN_OCC_PCT = 100.0 * FCRN_EFGAR / (FCRN_TOT + 1e-6);
  double IgG_Red_PCT  = 100.0 * (1.0 - ACHR_AB_C / AB0);
  double ACHR_DEN_PCT = 100.0 * ACHR_DEN / ACHR0;
  double SF_NORM = NMJ_SF / SF0;
  if(SF_NORM > 1.0) SF_NORM = 1.0;
  double QMG_OUT = QMG0 * (K_QMG_SF * (1.0 - SF_NORM + (1.0 - ACHR_DEN / ACHR0) * 0.5)
                  + K_QMG_AB * (ACHR_AB_C / AB0) * 0.3);
  if(QMG_OUT > QMG_MAX) QMG_OUT = QMG_MAX;
  if(QMG_OUT < 0.0) QMG_OUT = 0.0;
  double MGADL_OUT = QMG_OUT * (24.0 / 39.0);
  capture C_PYR_OUT, C_PRED_OUT, C_ECUL_OUT, C_EFGAR_OUT,
          AChE_INH_PCT, C5_INH_PCT, FCRN_OCC_PCT, IgG_Red_PCT,
          ACHR_DEN_PCT, COMP_ACT, NMJ_SF, ACHR_AB_C,
          TFH, GCB, SLPC, LLPC,
          QMG_OUT, MGADL_OUT;
'

## Compile model once
mod <- mcode("MG_QSP_shiny", mg_code, quiet = TRUE)

## ============================================================
## Helper: build event table from UI inputs
## ============================================================
build_events <- function(pyr_dose, pred_dose, aza_dose,
                          ecul_dose, efgar_dose,
                          sim_wk, use_ritux) {
  sim_h   <- sim_wk * 168
  ev_list <- list()

  if (pyr_dose > 0) {
    ev_list[["pyr"]] <- ev(amt = pyr_dose * 1e6, cmt = 1,
                            ii = 6, addl = floor(sim_h / 6) - 1, time = 0)
  }
  if (pred_dose > 0) {
    ev_list[["pred"]] <- ev(amt = pred_dose * 1e6, cmt = 3,
                             ii = 24, addl = floor(sim_h / 24) - 1, time = 0)
  }
  if (aza_dose > 0) {
    ev_list[["aza"]] <- ev(amt = aza_dose * 1e6, cmt = 6,
                            ii = 24, addl = floor(sim_h / 24) - 1, time = 0)
  }
  if (ecul_dose > 0) {
    ecul_times <- seq(0, sim_h - 336, by = 336)
    ev_list[["ecul"]] <- ev(amt = ecul_dose * 1e6, cmt = 9,
                             time = ecul_times)
  }
  if (efgar_dose > 0) {
    efgar_times <- seq(0, sim_h - 168, by = 168)
    efgar_nmol  <- efgar_dose * 1e6 / 50e3
    ev_list[["efgar"]] <- ev(amt = efgar_nmol, cmt = 11,
                              time = efgar_times)
  }

  if (length(ev_list) == 0) {
    return(ev(time = 0, amt = 0, cmt = 1))
  }
  Reduce(c, ev_list)
}

## ============================================================
## UI
## ============================================================
ui <- dashboardPage(
  skin = "blue",
  dashboardHeader(
    title = span("MG QSP Dashboard", style = "font-size:15px; font-weight:bold;"),
    titleWidth = 280
  ),

  dashboardSidebar(
    width = 280,
    sidebarMenu(
      menuItem("1. 환자 프로파일", tabName = "profile",
               icon = icon("user-circle")),
      menuItem("2. 약동학 (PK)", tabName = "pk",
               icon = icon("pills")),
      menuItem("3. PD 주요 지표", tabName = "pd",
               icon = icon("microscope")),
      menuItem("4. 임상 엔드포인트", tabName = "clinical",
               icon = icon("chart-line")),
      menuItem("5. 시나리오 비교", tabName = "scenario",
               icon = icon("layer-group")),
      menuItem("6. 바이오마커 패널", tabName = "biomarker",
               icon = icon("flask"))
    ),
    hr(),
    ## ---- Drug Dosing Controls ----
    h5("치료 설정 (Treatment Setup)", style = "text-align:center; color:#E0E0E0;"),
    sliderInput("pyr_dose", "Pyridostigmine (mg)",
                min = 0, max = 120, value = 60, step = 15),
    sliderInput("pred_dose", "Prednisolone (mg/d)",
                min = 0, max = 100, value = 60, step = 5),
    sliderInput("aza_dose", "Azathioprine (mg/d)",
                min = 0, max = 200, value = 0, step = 25),
    sliderInput("ecul_dose", "Eculizumab (mg q2w)",
                min = 0, max = 1200, value = 0, step = 300),
    sliderInput("efgar_dose", "Efgartigimod (mg qw)",
                min = 0, max = 1400, value = 0, step = 350),
    sliderInput("sim_wk", "시뮬레이션 기간 (weeks)",
                min = 4, max = 104, value = 52, step = 4),
    checkboxInput("use_ritux", "Rituximab (B-cell depletion)", value = FALSE),
    actionButton("run_sim", "시뮬레이션 실행",
                 icon = icon("play"),
                 style = "width:100%; background-color:#27AE60; color:white; margin-top:10px;")
  ),

  dashboardBody(
    tags$style(HTML("
      .box { border-radius:6px; }
      .info-box { border-radius:6px; }
      .small-box { border-radius:6px; }
    ")),
    tabItems(

      ## ==== Tab 1: Patient Profile & Disease Overview ====
      tabItem(
        tabName = "profile",
        fluidRow(
          box(
            title = "중증 근무력증 개요 (Myasthenia Gravis Overview)",
            status = "primary", solidHeader = TRUE, width = 12,
            fluidRow(
              column(6,
                h4("병태생리 (Pathophysiology)"),
                tags$ul(
                  tags$li("AChR (아세틸콜린 수용체) 자가항체가 NMJ 아세틸콜린 수용체를 파괴"),
                  tags$li("보체 활성화 → MAC (막공격 복합체) 형성 → NMJ 구조 손상"),
                  tags$li("AChR 항체 교차결합 → 수용체 내재화(internalization) 증가"),
                  tags$li("AChR 밀도 감소 → NMJ 안전계수(safety factor) 감소"),
                  tags$li("반복 자극 시 단연접전위(EPP) 역치 미달 → 근육 약화")
                ),
                h4("역학 (Epidemiology)"),
                tags$ul(
                  tags$li("유병률: 인구 10만명당 약 20명"),
                  tags$li("AChR-Ab (+): 약 85%; MuSK-Ab (+): 약 5-8%"),
                  tags$li("항체 음성(seronegative): 약 10%"),
                  tags$li("여성: 50세 이전 호발; 남성: 50세 이후 호발"),
                  tags$li("흉선종(thymoma) 동반: AChR-MG의 약 10-15%")
                )
              ),
              column(6,
                h4("치료 목표 (Treatment Goals)"),
                tags$ul(
                  tags$li(strong("QMG ≤ 6"), ": 최소 증상 또는 관해"),
                  tags$li(strong("MG-ADL ≤ 2"), ": 최소 증상 상태(MMS)"),
                  tags$li(strong("CSR"), ": 완전 안정 관해 (약물 없이 QMG=0, 12개월 이상)"),
                  tags$li("흉선 절제술(thymectomy): 비흉선종 AChR-MG에 효과 (MGTX 시험)"),
                  tags$li("위기(crisis) 예방: FVC, NIF 모니터링 필수")
                ),
                h4("MGFA 분류"),
                tableOutput("mgfa_table")
              )
            )
          )
        ),
        fluidRow(
          valueBoxOutput("vb_qmg",   width = 3),
          valueBoxOutput("vb_mgadl", width = 3),
          valueBoxOutput("vb_achr",  width = 3),
          valueBoxOutput("vb_sf",    width = 3)
        ),
        fluidRow(
          box(
            title = "질환 자연경과 (Natural History Simulation)",
            status = "warning", solidHeader = TRUE, width = 12,
            plotlyOutput("plot_natural_history", height = "300px")
          )
        )
      ),

      ## ==== Tab 2: PK ====
      tabItem(
        tabName = "pk",
        fluidRow(
          box(
            title = "피리도스티그민 PK (Pyridostigmine PK — 24h Profile)",
            status = "success", solidHeader = TRUE, width = 6,
            plotlyOutput("plot_pyr_pk", height = "280px")
          ),
          box(
            title = "프레드니솔론 PK (Prednisolone PK)",
            status = "info", solidHeader = TRUE, width = 6,
            plotlyOutput("plot_pred_pk", height = "280px")
          )
        ),
        fluidRow(
          box(
            title = "에쿨리주맙 PK (Eculizumab PK)",
            status = "warning", solidHeader = TRUE, width = 6,
            plotlyOutput("plot_ecul_pk", height = "280px")
          ),
          box(
            title = "에프가르티지모드 PK (Efgartigimod PK)",
            status = "danger", solidHeader = TRUE, width = 6,
            plotlyOutput("plot_efgar_pk", height = "280px")
          )
        ),
        fluidRow(
          box(
            title = "PK 파라미터 요약 (PK Parameter Summary)",
            status = "primary", solidHeader = TRUE, width = 12,
            DTOutput("pk_param_table")
          )
        )
      ),

      ## ==== Tab 3: PD Key Metrics ====
      tabItem(
        tabName = "pd",
        fluidRow(
          box(
            title = "AChE 억제율 (AChE Inhibition by Pyridostigmine)",
            status = "success", solidHeader = TRUE, width = 6,
            plotlyOutput("plot_ache_inh", height = "280px")
          ),
          box(
            title = "AChR 밀도 (NMJ AChR Density over Time)",
            status = "danger", solidHeader = TRUE, width = 6,
            plotlyOutput("plot_achr_den", height = "280px")
          )
        ),
        fluidRow(
          box(
            title = "보체 활성도 (Complement Activity: C5 Inhibition)",
            status = "warning", solidHeader = TRUE, width = 6,
            plotlyOutput("plot_complement", height = "280px")
          ),
          box(
            title = "NMJ 안전계수 (NMJ Safety Factor)",
            status = "info", solidHeader = TRUE, width = 6,
            plotlyOutput("plot_nmj_sf", height = "280px")
          )
        ),
        fluidRow(
          box(
            title = "면역세포 동태 (Immune Cell Dynamics: B Cells / Plasma Cells)",
            status = "primary", solidHeader = TRUE, width = 12,
            plotlyOutput("plot_immune_cells", height = "300px")
          )
        )
      ),

      ## ==== Tab 4: Clinical Endpoints ====
      tabItem(
        tabName = "clinical",
        fluidRow(
          box(
            title = "QMG Score 추이 (QMG Score over Time)",
            status = "danger", solidHeader = TRUE, width = 6,
            plotlyOutput("plot_qmg", height = "300px")
          ),
          box(
            title = "MG-ADL Score (Activities of Daily Living)",
            status = "warning", solidHeader = TRUE, width = 6,
            plotlyOutput("plot_mgadl", height = "300px")
          )
        ),
        fluidRow(
          box(
            title = "AChR 항체 역가 (AChR Antibody Titer over Time)",
            status = "info", solidHeader = TRUE, width = 6,
            plotlyOutput("plot_ab_titer", height = "280px")
          ),
          box(
            title = "IgG 감소율 (FcRn Blockade: IgG Reduction %)",
            status = "success", solidHeader = TRUE, width = 6,
            plotlyOutput("plot_igg_red", height = "280px")
          )
        ),
        fluidRow(
          box(
            title = "주요 시점별 임상 지표 요약",
            status = "primary", solidHeader = TRUE, width = 12,
            DTOutput("clinical_summary_table")
          )
        )
      ),

      ## ==== Tab 5: Scenario Comparison ====
      tabItem(
        tabName = "scenario",
        fluidRow(
          box(
            title = "치료 시나리오 비교 설정",
            status = "primary", solidHeader = TRUE, width = 3,
            checkboxGroupInput("scenarios_sel",
                               "비교할 시나리오 선택:",
                               choices = c("Untreated" = "s1",
                                           "Pyridostigmine" = "s2",
                                           "Pred + Pyridostigmine" = "s3",
                                           "Pred + AZA + Pyridostigmine" = "s4",
                                           "Eculizumab + Pyridostigmine" = "s5",
                                           "Efgartigimod + Pyridostigmine" = "s6",
                                           "Rituximab + Pyridostigmine" = "s7"),
                               selected = c("s1", "s2", "s3", "s5", "s6")),
            selectInput("comp_endpoint",
                        "비교 지표:",
                        choices = c("QMG Score"       = "QMG_OUT",
                                    "MG-ADL"          = "MGADL_OUT",
                                    "AChR-Ab (nmol/L)"= "ACHR_AB_C",
                                    "AChR Density (%)"= "ACHR_DEN_PCT",
                                    "NMJ Safety Factor"= "NMJ_SF",
                                    "Complement"       = "COMP_ACT",
                                    "IgG Reduction (%)"= "IgG_Red_PCT"),
                        selected = "QMG_OUT"),
            actionButton("run_comparison", "비교 실행",
                         icon = icon("balance-scale"),
                         style = "width:100%; background-color:#8E44AD; color:white;")
          ),
          box(
            title = "치료 시나리오 비교 결과",
            status = "purple", solidHeader = TRUE, width = 9,
            plotlyOutput("plot_scenario_comp", height = "420px")
          )
        ),
        fluidRow(
          box(
            title = "52주 임상 결과 비교표",
            status = "primary", solidHeader = TRUE, width = 12,
            DTOutput("scenario_result_table")
          )
        )
      ),

      ## ==== Tab 6: Biomarker Panel ====
      tabItem(
        tabName = "biomarker",
        fluidRow(
          box(
            title = "진단 바이오마커 패널 (Diagnostic Biomarker Panel)",
            status = "primary", solidHeader = TRUE, width = 5,
            tableOutput("biomarker_table")
          ),
          box(
            title = "항체 역가 변화율 (Antibody Titer Change Over Time)",
            status = "info", solidHeader = TRUE, width = 7,
            plotlyOutput("plot_biomarker_ab", height = "320px")
          )
        ),
        fluidRow(
          box(
            title = "SFEMG Jitter (단일 섬유 EMG, NMJ 전달 이상도)",
            status = "warning", solidHeader = TRUE, width = 6,
            plotlyOutput("plot_jitter", height = "280px")
          ),
          box(
            title = "형질세포 서브타입 동태 (Plasma Cell Subtype Dynamics)",
            status = "danger", solidHeader = TRUE, width = 6,
            plotlyOutput("plot_pc_dynamics", height = "280px")
          )
        ),
        fluidRow(
          box(
            title = "바이오마커 임상 해석 가이드",
            status = "success", solidHeader = TRUE, width = 12,
            DTOutput("biomarker_guide_table")
          )
        )
      )
    )  # end tabItems
  )
)

## ============================================================
## Server
## ============================================================
server <- function(input, output, session) {

  ## ---- Static tables ----
  output$mgfa_table <- renderTable({
    data.frame(
      Class = c("I", "IIa", "IIb", "IIIa", "IIIb", "IVa", "IVb", "V"),
      Description = c(
        "안구 근육만 침범, 눈꺼풀 하수 또는 복시",
        "주로 사지/몸통 침범 (경증)",
        "주로 구인두/호흡 근육 침범 (경증)",
        "주로 사지/몸통 침범 (중등도)",
        "주로 구인두/호흡 침범 (중등도)",
        "주로 사지/몸통 침범 (중증)",
        "주로 구인두/호흡 침범 (중증)",
        "삽관 필요 (위기)"
      )
    )
  }, striped = TRUE, bordered = TRUE)

  output$biomarker_guide_table <- renderDT({
    df <- data.frame(
      바이오마커 = c("Anti-AChR Ab", "Anti-MuSK Ab", "Anti-LRP4 Ab",
                   "Anti-Titin Ab", "SFEMG Jitter", "RNS 감소율",
                   "Total IgG", "sC5b-9", "Blood Tfh %", "CH50"),
      정상 범위 = c("<0.4 nmol/L (RIA)", "음성", "음성",
                   "음성", "<55 µs", "<10% at 3Hz",
                   "7-16 g/L", "<244 ng/mL", "~1-2% CD4+", "60-150 U/mL"),
      임상 의의 = c("MG 진단 및 질환 활성도", "MuSK-MG (흉부 CT, IgG4)",
                   "Seroneg MG 진단", "흉선종 관련 MG",
                   "NMJ 전달 이상 민감도 높음", "신경근육 전달 차단 객관적 지표",
                   "FcRn 치료 모니터링", "MAC 형성도 (에쿨리주맙 치료 반응)",
                   "면역 활성도", "보체 소비 (활성 MG에서 감소)")
    )
    datatable(df, options = list(pageLength = 10, scrollX = TRUE),
              rownames = FALSE)
  })

  output$pk_param_table <- renderDT({
    df <- data.frame(
      Drug = c("Pyridostigmine", "Prednisolone", "Azathioprine/6-TGN",
               "Eculizumab", "Efgartigimod alfa", "Rituximab"),
      Route = c("PO", "PO", "PO", "IV", "IV", "IV"),
      Bioavailability = c("13-25%", "82%", "50% (AZA)", "N/A", "N/A", "N/A"),
      "Vd (L)" = c("90 (1.3 L/kg)", "35 (0.5 L/kg)", "120 (RBC)", "7.7", "6.0", "8.1"),
      "t½" = c("1.9h", "2.3h", "5d (TGN)", "11d", "3.3d", "22d"),
      "Dose & Freq" = c("60mg q6h", "0.5-1.5 mg/kg/d", "2-3 mg/kg/d",
                        "900mg q2w", "10 mg/kg qw", "375mg/m² ×4"),
      "Main CYP" = c("Plasma ChE", "CYP3A4", "TPMT/HGPRT",
                     "Protein degradation", "FcRn recycling", "CD20 TMDD")
    )
    datatable(df, options = list(scrollX = TRUE), rownames = FALSE)
  })

  ## ---- Reactive simulation ----
  sim_result <- eventReactive(input$run_sim, {
    events <- build_events(
      pyr_dose  = input$pyr_dose,
      pred_dose = input$pred_dose,
      aza_dose  = input$aza_dose,
      ecul_dose = input$ecul_dose,
      efgar_dose = input$efgar_dose,
      sim_wk    = input$sim_wk,
      use_ritux = input$use_ritux
    )

    custom_params <- if (input$use_ritux) {
      list(K_GCB_IN = 0.005, K_SLPC_F = 0.05, K_LLPC_F = 0.02)
    } else {
      list()
    }

    out <- mod %>%
      param(custom_params) %>%
      mrgsim(ev = events,
             end = input$sim_wk * 168,
             delta = 24) %>%
      as.data.frame()
    out$time_wk <- out$time / 168
    out
  }, ignoreNULL = FALSE)

  ## --- Value boxes (Tab 1) ---
  get_latest <- function(col) {
    d <- sim_result()
    if (is.null(d)) return(NA)
    tail(d[[col]], 1)
  }

  output$vb_qmg <- renderValueBox({
    val <- round(get_latest("QMG_OUT"), 1)
    color <- if (!is.na(val) && val <= 6) "green" else if (!is.na(val) && val <= 15) "yellow" else "red"
    valueBox(value = val, subtitle = "QMG Score (최종, 0-39)",
             icon = icon("heartbeat"), color = color)
  })

  output$vb_mgadl <- renderValueBox({
    val <- round(get_latest("MGADL_OUT"), 1)
    color <- if (!is.na(val) && val <= 2) "green" else if (!is.na(val) && val <= 6) "yellow" else "red"
    valueBox(value = val, subtitle = "MG-ADL (최종, 0-24)",
             icon = icon("walking"), color = color)
  })

  output$vb_achr <- renderValueBox({
    val <- round(get_latest("ACHR_AB_C"), 2)
    color <- if (!is.na(val) && val <= 0.4) "green" else if (!is.na(val) && val <= 1.5) "yellow" else "red"
    valueBox(value = paste0(val, " nmol/L"), subtitle = "AChR-Ab Titer",
             icon = icon("vial"), color = color)
  })

  output$vb_sf <- renderValueBox({
    val <- round(get_latest("NMJ_SF"), 2)
    color <- if (!is.na(val) && val >= 1.0) "green" else if (!is.na(val) && val >= 0.6) "yellow" else "red"
    valueBox(value = val, subtitle = "NMJ Safety Factor",
             icon = icon("bolt"), color = color)
  })

  ## --- Plot: Natural History (Tab 1) ---
  output$plot_natural_history <- renderPlotly({
    d <- sim_result()
    if (is.null(d)) return(NULL)
    p <- ggplot(d, aes(x = time_wk, y = QMG_OUT)) +
      geom_line(color = "#E74C3C", linewidth = 1.2) +
      geom_area(fill = "#E74C3C", alpha = 0.15) +
      labs(x = "Time (weeks)", y = "QMG Score",
           title = "QMG Score under Current Treatment Settings") +
      geom_hline(yintercept = 6, linetype = "dashed", color = "darkgreen") +
      annotate("text", x = max(d$time_wk) * 0.9, y = 7, label = "MMS threshold",
               color = "darkgreen", size = 3) +
      ylim(0, 40) + theme_bw(base_size = 11)
    ggplotly(p)
  })

  ## --- PK Plots (Tab 2) ---
  output$plot_pyr_pk <- renderPlotly({
    d <- sim_result()
    if (is.null(d) || input$pyr_dose == 0) return(NULL)
    d24 <- d[d$time <= 24, ]
    p <- ggplot(d24, aes(x = time, y = C_PYR_OUT)) +
      geom_line(color = "#27AE60", linewidth = 1.3) +
      geom_vline(xintercept = c(0, 6, 12, 18), linetype = "dotted", color = "grey60") +
      labs(x = "Time (h)", y = "Conc (ng/mL)",
           title = paste0("Pyridostigmine ", input$pyr_dose, "mg q6h (24h)")) +
      theme_bw(base_size = 11)
    ggplotly(p)
  })

  output$plot_pred_pk <- renderPlotly({
    d <- sim_result()
    if (is.null(d) || input$pred_dose == 0) return(NULL)
    p <- ggplot(d, aes(x = time_wk, y = C_PRED_OUT)) +
      geom_line(color = "#2980B9", linewidth = 1.0) +
      labs(x = "Time (weeks)", y = "Conc (ng/mL)",
           title = paste0("Prednisolone ", input$pred_dose, "mg/d")) +
      theme_bw(base_size = 11)
    ggplotly(p)
  })

  output$plot_ecul_pk <- renderPlotly({
    d <- sim_result()
    if (is.null(d) || input$ecul_dose == 0) {
      return(plotly_empty() %>% layout(title = "Eculizumab: dose = 0"))
    }
    p <- ggplot(d, aes(x = time_wk, y = C_ECUL_OUT)) +
      geom_line(color = "#E67E22", linewidth = 1.0) +
      labs(x = "Time (weeks)", y = "Conc (ng/mL)",
           title = paste0("Eculizumab ", input$ecul_dose, "mg q2w")) +
      theme_bw(base_size = 11)
    ggplotly(p)
  })

  output$plot_efgar_pk <- renderPlotly({
    d <- sim_result()
    if (is.null(d) || input$efgar_dose == 0) {
      return(plotly_empty() %>% layout(title = "Efgartigimod: dose = 0"))
    }
    p <- ggplot(d, aes(x = time_wk, y = C_EFGAR_OUT)) +
      geom_line(color = "#8E44AD", linewidth = 1.0) +
      labs(x = "Time (weeks)", y = "Conc (nmol/L)",
           title = paste0("Efgartigimod ", input$efgar_dose, "mg qw")) +
      theme_bw(base_size = 11)
    ggplotly(p)
  })

  ## --- PD Plots (Tab 3) ---
  output$plot_ache_inh <- renderPlotly({
    d <- sim_result()
    if (is.null(d)) return(NULL)
    p <- ggplot(d, aes(x = time_wk, y = AChE_INH_PCT)) +
      geom_line(color = "#8E44AD", linewidth = 1.1) +
      labs(x = "Time (weeks)", y = "AChE Inhibition (%)",
           title = "AChE Inhibition (Pyridostigmine)") +
      ylim(0, 100) + theme_bw(base_size = 11)
    ggplotly(p)
  })

  output$plot_achr_den <- renderPlotly({
    d <- sim_result()
    if (is.null(d)) return(NULL)
    p <- ggplot(d, aes(x = time_wk, y = ACHR_DEN_PCT)) +
      geom_line(color = "#C0392B", linewidth = 1.1) +
      geom_hline(yintercept = 100, linetype = "dashed", color = "grey50") +
      labs(x = "Time (weeks)", y = "AChR Density (%)",
           title = "NMJ AChR Density") +
      ylim(0, 110) + theme_bw(base_size = 11)
    ggplotly(p)
  })

  output$plot_complement <- renderPlotly({
    d <- sim_result()
    if (is.null(d)) return(NULL)
    p <- ggplot(d, aes(x = time_wk, y = COMP_ACT)) +
      geom_line(color = "#E67E22", linewidth = 1.1) +
      geom_hline(yintercept = 1, linetype = "dashed", color = "grey50") +
      labs(x = "Time (weeks)", y = "Complement Activity",
           title = "Complement Activity (C5 Pathway)") +
      theme_bw(base_size = 11)
    ggplotly(p)
  })

  output$plot_nmj_sf <- renderPlotly({
    d <- sim_result()
    if (is.null(d)) return(NULL)
    p <- ggplot(d, aes(x = time_wk, y = NMJ_SF)) +
      geom_line(color = "#2980B9", linewidth = 1.1) +
      geom_hline(yintercept = 1, linetype = "dashed", color = "grey50") +
      labs(x = "Time (weeks)", y = "Safety Factor",
           title = "NMJ Safety Factor") +
      theme_bw(base_size = 11)
    ggplotly(p)
  })

  output$plot_immune_cells <- renderPlotly({
    d <- sim_result()
    if (is.null(d)) return(NULL)
    d_long <- d %>%
      select(time_wk, TFH, GCB, SLPC, LLPC) %>%
      pivot_longer(-time_wk, names_to = "Cell", values_to = "Level")
    p <- ggplot(d_long, aes(x = time_wk, y = Level, color = Cell)) +
      geom_line(linewidth = 1.0) +
      scale_color_manual(values = c(TFH = "#3498DB", GCB = "#E74C3C",
                                    SLPC = "#F39C12", LLPC = "#1ABC9C")) +
      labs(x = "Time (weeks)", y = "Relative Level",
           title = "Immune Cell Dynamics") +
      theme_bw(base_size = 11) + theme(legend.position = "right")
    ggplotly(p)
  })

  ## --- Clinical Plots (Tab 4) ---
  output$plot_qmg <- renderPlotly({
    d <- sim_result()
    if (is.null(d)) return(NULL)
    p <- ggplot(d, aes(x = time_wk, y = QMG_OUT)) +
      geom_line(color = "#E74C3C", linewidth = 1.3) +
      geom_hline(yintercept = c(6, 12, 20), linetype = "dashed",
                 color = c("darkgreen", "orange", "red")) +
      annotate("text", x = max(d$time_wk) * 0.85, y = c(5, 11, 19),
               label = c("MMS (<6)", "Mild (<12)", "Moderate (<20)"),
               size = 2.8, color = c("darkgreen", "orange", "red")) +
      labs(x = "Weeks", y = "QMG Score (0-39)", title = "QMG Score") +
      ylim(0, 39) + theme_bw(base_size = 11)
    ggplotly(p)
  })

  output$plot_mgadl <- renderPlotly({
    d <- sim_result()
    if (is.null(d)) return(NULL)
    p <- ggplot(d, aes(x = time_wk, y = MGADL_OUT)) +
      geom_line(color = "#F39C12", linewidth = 1.3) +
      geom_hline(yintercept = 2, linetype = "dashed", color = "darkgreen") +
      labs(x = "Weeks", y = "MG-ADL (0-24)", title = "MG-ADL Score") +
      ylim(0, 24) + theme_bw(base_size = 11)
    ggplotly(p)
  })

  output$plot_ab_titer <- renderPlotly({
    d <- sim_result()
    if (is.null(d)) return(NULL)
    p <- ggplot(d, aes(x = time_wk, y = ACHR_AB_C)) +
      geom_line(color = "#9B59B6", linewidth = 1.1) +
      geom_hline(yintercept = 0.4, linetype = "dashed", color = "red") +
      labs(x = "Weeks", y = "Anti-AChR IgG (nmol/L)",
           title = "AChR Antibody Titer") +
      theme_bw(base_size = 11)
    ggplotly(p)
  })

  output$plot_igg_red <- renderPlotly({
    d <- sim_result()
    if (is.null(d)) return(NULL)
    p <- ggplot(d, aes(x = time_wk, y = IgG_Red_PCT)) +
      geom_line(color = "#27AE60", linewidth = 1.1) +
      labs(x = "Weeks", y = "IgG Reduction (%)",
           title = "Total IgG Reduction (Efgartigimod / FcRn)") +
      ylim(-10, 80) + theme_bw(base_size = 11)
    ggplotly(p)
  })

  output$clinical_summary_table <- renderDT({
    d <- sim_result()
    if (is.null(d)) return(NULL)
    key_wk <- c(0, 4, 12, 26, 52)
    tbl <- d %>%
      mutate(wk = round(time_wk, 0)) %>%
      filter(wk %in% key_wk) %>%
      group_by(wk) %>%
      slice(1) %>%
      ungroup() %>%
      select(wk, QMG_OUT, MGADL_OUT, ACHR_AB_C, ACHR_DEN_PCT, NMJ_SF,
             C5_INH_PCT, IgG_Red_PCT, AChE_INH_PCT) %>%
      mutate(across(where(is.numeric), ~round(., 2))) %>%
      rename(Week = wk, QMG = QMG_OUT, "MG-ADL" = MGADL_OUT,
             "AChR-Ab(nmol/L)" = ACHR_AB_C, "AChR Den(%)" = ACHR_DEN_PCT,
             "NMJ SF" = NMJ_SF, "C5 Inh(%)" = C5_INH_PCT,
             "IgG Red(%)" = IgG_Red_PCT, "AChE Inh(%)" = AChE_INH_PCT)
    datatable(tbl, options = list(scrollX = TRUE, pageLength = 8), rownames = FALSE)
  })

  ## --- Scenario Comparison (Tab 5) ---
  scenario_sims <- eventReactive(input$run_comparison, {
    sim_h <- input$sim_wk * 168

    scenarios <- list(
      s1 = list(name = "Untreated", pyr = 0, pred = 0, aza = 0, ecul = 0, efgar = 0, ritux = FALSE),
      s2 = list(name = "Pyridostigmine", pyr = 60, pred = 0, aza = 0, ecul = 0, efgar = 0, ritux = FALSE),
      s3 = list(name = "Pred + Pyridostigmine", pyr = 60, pred = 60, aza = 0, ecul = 0, efgar = 0, ritux = FALSE),
      s4 = list(name = "Pred + AZA + Pyridostigmine", pyr = 60, pred = 60, aza = 150, ecul = 0, efgar = 0, ritux = FALSE),
      s5 = list(name = "Eculizumab + Pyridostigmine", pyr = 60, pred = 0, aza = 0, ecul = 900, efgar = 0, ritux = FALSE),
      s6 = list(name = "Efgartigimod + Pyridostigmine", pyr = 60, pred = 0, aza = 0, ecul = 0, efgar = 700, ritux = FALSE),
      s7 = list(name = "Rituximab + Pyridostigmine", pyr = 60, pred = 0, aza = 0, ecul = 0, efgar = 0, ritux = TRUE)
    )

    sel <- input$scenarios_sel
    results <- lapply(sel, function(s) {
      sc <- scenarios[[s]]
      ev_s <- build_events(sc$pyr, sc$pred, sc$aza, sc$ecul, sc$efgar,
                           input$sim_wk, sc$ritux)
      cp <- if (sc$ritux) list(K_GCB_IN = 0.005, K_SLPC_F = 0.05) else list()
      out <- mod %>%
        param(cp) %>%
        mrgsim(ev = ev_s, end = sim_h, delta = 24) %>%
        as.data.frame()
      out$scenario <- sc$name
      out$time_wk  <- out$time / 168
      out
    })
    bind_rows(results)
  }, ignoreNULL = FALSE)

  output$plot_scenario_comp <- renderPlotly({
    d <- scenario_sims()
    if (is.null(d)) return(NULL)
    ep <- input$comp_endpoint
    p <- ggplot(d, aes_string(x = "time_wk", y = ep, color = "scenario")) +
      geom_line(linewidth = 1.1) +
      labs(x = "Time (weeks)", y = ep, color = "Scenario",
           title = paste("Treatment Scenario Comparison:", ep)) +
      theme_bw(base_size = 11) +
      theme(legend.position = "bottom", legend.text = element_text(size = 8))
    ggplotly(p)
  })

  output$scenario_result_table <- renderDT({
    d <- scenario_sims()
    if (is.null(d)) return(NULL)
    tbl <- d %>%
      group_by(scenario) %>%
      filter(time_wk == max(time_wk)) %>%
      slice(1) %>%
      ungroup() %>%
      select(scenario, QMG_OUT, MGADL_OUT, ACHR_AB_C,
             ACHR_DEN_PCT, NMJ_SF, C5_INH_PCT, IgG_Red_PCT) %>%
      mutate(across(where(is.numeric), ~round(., 2))) %>%
      rename(Scenario = scenario, QMG = QMG_OUT, "MG-ADL" = MGADL_OUT,
             "AChR-Ab" = ACHR_AB_C, "AChR Den(%)" = ACHR_DEN_PCT,
             "NMJ SF" = NMJ_SF, "C5 Inh(%)" = C5_INH_PCT,
             "IgG Red(%)" = IgG_Red_PCT)
    datatable(tbl, options = list(scrollX = TRUE, pageLength = 10),
              rownames = FALSE)
  })

  ## --- Biomarker Tab ---
  output$biomarker_table <- renderTable({
    d <- sim_result()
    if (is.null(d)) return(NULL)
    latest <- tail(d, 1)
    data.frame(
      Biomarker = c("AChR-Ab titer", "AChR density", "NMJ Safety Factor",
                    "Complement activity", "IgG pool", "AChE inhibition",
                    "C5 inhibition", "Short-lived PC", "Long-lived PC"),
      Value = c(
        paste0(round(latest$ACHR_AB_C, 2), " nmol/L"),
        paste0(round(latest$ACHR_DEN_PCT, 1), "%"),
        round(latest$NMJ_SF, 2),
        round(latest$COMP_ACT, 2),
        paste0(round(latest$IgG_Red_PCT, 1), "% red."),
        paste0(round(latest$AChE_INH_PCT, 1), "%"),
        paste0(round(latest$C5_INH_PCT, 1), "%"),
        round(latest$SLPC, 2),
        round(latest$LLPC, 2)
      ),
      Status = c(
        ifelse(latest$ACHR_AB_C <= 0.4, "Normal", "Elevated"),
        ifelse(latest$ACHR_DEN_PCT >= 90, "Normal", "Reduced"),
        ifelse(latest$NMJ_SF >= 0.9, "Normal", "Reduced"),
        ifelse(latest$COMP_ACT <= 1.1, "Normal", "Elevated"),
        ifelse(latest$IgG_Red_PCT >= 0, "Reduced", "Baseline"),
        ifelse(latest$AChE_INH_PCT > 30, "Inhibited", "Minimal"),
        ifelse(latest$C5_INH_PCT > 80, "Inhibited", "Minimal"),
        round(latest$SLPC, 2),
        round(latest$LLPC, 2)
      )
    )
  }, striped = TRUE, bordered = TRUE)

  output$plot_biomarker_ab <- renderPlotly({
    d <- sim_result()
    if (is.null(d)) return(NULL)
    p <- ggplot(d, aes(x = time_wk, y = ACHR_AB_C)) +
      geom_line(color = "#9B59B6", linewidth = 1.2) +
      geom_ribbon(aes(ymin = 0, ymax = ACHR_AB_C), fill = "#9B59B6", alpha = 0.15) +
      geom_hline(yintercept = 0.4, linetype = "dashed", color = "red") +
      labs(x = "Weeks", y = "Anti-AChR IgG (nmol/L)",
           title = "AChR Antibody Titer Dynamics") +
      theme_bw(base_size = 11)
    ggplotly(p)
  })

  output$plot_jitter <- renderPlotly({
    d <- sim_result()
    if (is.null(d)) return(NULL)
    d <- d %>%
      mutate(SFEMG_jitter = 55 + (1 - NMJ_SF / max(NMJ_SF, na.rm = TRUE)) * 80)
    p <- ggplot(d, aes(x = time_wk, y = SFEMG_jitter)) +
      geom_line(color = "#E67E22", linewidth = 1.1) +
      geom_hline(yintercept = 55, linetype = "dashed", color = "darkgreen") +
      annotate("text", x = max(d$time_wk) * 0.8, y = 53,
               label = "Normal (<55 µs)", color = "darkgreen", size = 3) +
      labs(x = "Weeks", y = "SFEMG Jitter (µs, estimated)",
           title = "SFEMG Jitter (NMJ Safety Proxy)") +
      theme_bw(base_size = 11)
    ggplotly(p)
  })

  output$plot_pc_dynamics <- renderPlotly({
    d <- sim_result()
    if (is.null(d)) return(NULL)
    d_long <- d %>%
      select(time_wk, SLPC, LLPC) %>%
      pivot_longer(-time_wk, names_to = "PC_type", values_to = "Level")
    p <- ggplot(d_long, aes(x = time_wk, y = Level, color = PC_type, fill = PC_type)) +
      geom_line(linewidth = 1.1) +
      geom_area(alpha = 0.15, position = "identity") +
      scale_color_manual(values = c(SLPC = "#E74C3C", LLPC = "#3498DB"),
                         labels = c("Short-lived PC", "Long-lived PC")) +
      scale_fill_manual(values = c(SLPC = "#E74C3C", LLPC = "#3498DB"),
                        labels = c("Short-lived PC", "Long-lived PC")) +
      labs(x = "Weeks", y = "Relative Level",
           title = "Plasma Cell Subtype Dynamics",
           color = NULL, fill = NULL) +
      theme_bw(base_size = 11) + theme(legend.position = "bottom")
    ggplotly(p)
  })
}

## ============================================================
## Run
## ============================================================
shinyApp(ui = ui, server = server)
