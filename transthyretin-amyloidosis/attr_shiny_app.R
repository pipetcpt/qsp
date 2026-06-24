##############################################################################
# Transthyretin Amyloidosis (ATTR) — Interactive QSP Shiny Dashboard
# 트랜스티레틴 아밀로이드증 대화형 QSP 시뮬레이션 대시보드
#
# 탭 구성 (8개):
#   Tab 1: 환자 프로파일 설정 (Patient Profile)
#   Tab 2: TTR 약동학 & 응집 (TTR PK/Aggregation)
#   Tab 3: 심장 ATTR (Cardiac Endpoints)
#   Tab 4: 신경 ATTR (Neurological Endpoints)
#   Tab 5: 치료 시나리오 비교 (Scenario Comparison)
#   Tab 6: 바이오마커 추적 (Biomarker Tracking)
#   Tab 7: 민감도 분석 (Sensitivity Analysis)
#   Tab 8: 가상 환자 집단 (Virtual Population)
##############################################################################

library(shiny)
library(shinydashboard)
library(mrgsolve)
library(dplyr)
library(ggplot2)
library(plotly)
library(DT)
library(tidyr)
library(patchwork)
library(shinyWidgets)

# ============================================================
# mrgsolve 모델 코드 (인라인)
# ============================================================
attr_model_code <- '
$PARAM
kprod_TTR=104 kel_TTR=0.347 k12_TTR=0.1 k21_TTR=0.05
mut_factor=1.0 kdiss_base=5e-4 kconf=0.05 kolig=0.10 kfib=0.02
frac_heart=0.45 frac_nerve=0.35 frac_GI=0.20 kclear_amyl=0.003
LV_base=10 LV_kgrowth=0.8 BNP_base=80 BNP_scale=15
TROP_base=5 TROP_scale=0.8 CARD_base=1.0 CARD_decline=0.15
NIS_base=10 NIS_scale=5.0
ka_TAF=8.328 CL_TAF=7.1 Vd_TAF=18.0 F_TAF=1.0
EC50_TAF=2.0 Emax_TAF=0.80
kout_PAT=0.035 Emax_PAT=0.80 EC50_PAT=0.5
kout_VUT=0.025 Emax_VUT=0.87 EC50_VUT=0.5
kout_INO=0.045 Emax_INO=0.75 EC50_INO=0.5
SixMWT_base=420 SixMWT_k=0.5 KCCQ_base=70 KCCQ_k=2.5

$CMT
TTR_C TTR_P TTR_MF TTR_OL
AMY_H AMY_N AMY_GI
LV_THICK BNP_C TROP_C CARD_FUNC
NIS_TOT AUTO_NP
TAF_GUT TAF_C PAT_EFF VUT_EFF INO_EFF
SixMWT KCCQ_IDX

$INIT
TTR_C=300 TTR_P=150 TTR_MF=0.001 TTR_OL=0.0001
AMY_H=0 AMY_N=0 AMY_GI=0
LV_THICK=10 BNP_C=80 TROP_C=5 CARD_FUNC=1
NIS_TOT=10 AUTO_NP=0.05
TAF_GUT=0 TAF_C=0 PAT_EFF=0 VUT_EFF=0 INO_EFF=0
SixMWT=420 KCCQ_IDX=70

$MAIN
double ke_TAF = CL_TAF/Vd_TAF;
double Inh_TAF = (Emax_TAF*TAF_C)/(EC50_TAF+TAF_C+1e-10);
double kdiss_eff = kdiss_base * mut_factor * (1.0-Inh_TAF);
double Inh_PAT = (Emax_PAT*PAT_EFF)/(EC50_PAT+PAT_EFF+1e-10);
double Inh_VUT = (Emax_VUT*VUT_EFF)/(EC50_VUT+VUT_EFF+1e-10);
double Inh_INO = (Emax_INO*INO_EFF)/(EC50_INO+INO_EFF+1e-10);
double Inh_total = 1.0-(1.0-Inh_PAT)*(1.0-Inh_VUT)*(1.0-Inh_INO);
double kprod_eff = kprod_TTR*(1.0-Inh_total);
double TTR_avail = (TTR_C>0)?TTR_C:0.0;
double rate_dissoc = kdiss_eff*TTR_avail;

$ODE
dxdt_TTR_C = kprod_eff - kel_TTR*TTR_C - k12_TTR*TTR_C + k21_TTR*TTR_P - rate_dissoc;
dxdt_TTR_P = k12_TTR*TTR_C - k21_TTR*TTR_P;
dxdt_TTR_MF = kconf*rate_dissoc - kolig*TTR_MF;
dxdt_TTR_OL = kolig*TTR_MF - kfib*TTR_OL;
double rate_fib = kfib*TTR_OL;
dxdt_AMY_H  = frac_heart*rate_fib - kclear_amyl*AMY_H;
dxdt_AMY_N  = frac_nerve*rate_fib - kclear_amyl*AMY_N;
dxdt_AMY_GI = frac_GI*rate_fib   - kclear_amyl*AMY_GI;
double LV_tgt = LV_base + LV_kgrowth*AMY_H;
dxdt_LV_THICK = 0.05*(LV_tgt-LV_THICK);
double BNP_tgt = BNP_base + BNP_scale*AMY_H + 50*TTR_OL;
dxdt_BNP_C = 0.1*(BNP_tgt-BNP_C);
double TROP_tgt = TROP_base + TROP_scale*AMY_H + 0.5*TTR_OL;
dxdt_TROP_C = 0.1*(TROP_tgt-TROP_C);
double CARD_tgt = CARD_base - CARD_decline*AMY_H;
if(CARD_tgt<0.1) CARD_tgt=0.1;
dxdt_CARD_FUNC = 0.02*(CARD_tgt-CARD_FUNC);
double NIS_tgt = NIS_base + NIS_scale*AMY_N;
dxdt_NIS_TOT = 0.05*(NIS_tgt-NIS_TOT);
double AUTO_tgt = 0.05+0.3*AMY_N;
if(AUTO_tgt>1.0) AUTO_tgt=1.0;
dxdt_AUTO_NP = 0.05*(AUTO_tgt-AUTO_NP);
dxdt_TAF_GUT = -ka_TAF*TAF_GUT;
dxdt_TAF_C   = (ka_TAF*TAF_GUT*F_TAF)/Vd_TAF - ke_TAF*TAF_C;
dxdt_PAT_EFF = -kout_PAT*PAT_EFF;
dxdt_VUT_EFF = -kout_VUT*VUT_EFF;
dxdt_INO_EFF = -kout_INO*INO_EFF;
double SixMWT_tgt = SixMWT_base*CARD_FUNC - SixMWT_k*AMY_H;
if(SixMWT_tgt<50) SixMWT_tgt=50;
dxdt_SixMWT = 0.05*(SixMWT_tgt-SixMWT);
double KCCQ_tgt = KCCQ_base*CARD_FUNC*1.1 - KCCQ_k*AMY_H;
if(KCCQ_tgt<0) KCCQ_tgt=0;
if(KCCQ_tgt>100) KCCQ_tgt=100;
dxdt_KCCQ_IDX = 0.05*(KCCQ_tgt-KCCQ_IDX);

$TABLE
capture TTR_mgdL = TTR_C/10;
capture CardAmyloid = AMY_H;
capture NerveAmyloid = AMY_N;
capture LV_wall = LV_THICK;
capture NT_proBNP = BNP_C;
capture TroponinT = TROP_C;
capture CardFunc = CARD_FUNC;
capture NIS_score = NIS_TOT;
capture Walk6min = SixMWT;
capture KCCQ = KCCQ_IDX;
capture TAF_plasma = TAF_C;
capture Inh_TAF_pct = 100*(Emax_TAF*TAF_C)/(EC50_TAF+TAF_C+1e-10);
capture NAC_stage = (NT_proBNP>=3000 && TroponinT>=50)?3.0:(NT_proBNP>=3000||TroponinT>=50)?2.0:1.0;

$CAPTURE TTR_mgdL CardAmyloid NerveAmyloid LV_wall NT_proBNP TroponinT
         CardFunc NIS_score Walk6min KCCQ TAF_plasma Inh_TAF_pct NAC_stage
'

# 모델 컴파일
mod_shiny <- tryCatch(
  mcode("ATTR_Shiny", attr_model_code),
  error = function(e) { message("모델 컴파일 오류: ", e$message); NULL }
)

# ============================================================
# 시뮬레이션 실행 함수
# ============================================================
run_sim <- function(mod, params, treatment, end_years, delta_days = 14) {
  if (is.null(mod)) return(NULL)
  end_days <- end_years * 365

  # 파라미터 업데이트
  mod_run <- mod %>% param(params)

  # 투여 이벤트 설정
  events <- NULL

  if (treatment %in% c("tafamidis", "combo")) {
    taf_ev <- ev(time = 0, amt = 61, cmt = "TAF_GUT", evid = 1,
                 rate = 0, ii = 1, addl = end_days - 1)
    events <- taf_ev
  }

  if (treatment == "patisiran") {
    n_d <- floor(end_days / 21) + 1
    pat_ev <- data.frame(
      time = seq(0, by = 21, length.out = n_d),
      amt = 1.0, cmt = "PAT_EFF", evid = 1, rate = 0, ID = 1
    )
    out <- mod_run %>%
      mrgsim_df(data = pat_ev, end = end_days, delta = delta_days,
                add = list(ID = 1)) %>%
      mutate(scenario = "파티시란 Q3W IV")
    return(out)
  }

  if (treatment == "vutrisiran") {
    n_d <- floor(end_days / 91) + 1
    vut_ev <- data.frame(
      time = seq(0, by = 91, length.out = n_d),
      amt = 1.0, cmt = "VUT_EFF", evid = 1, rate = 0, ID = 1
    )
    out <- mod_run %>%
      mrgsim_df(data = vut_ev, end = end_days, delta = delta_days,
                add = list(ID = 1)) %>%
      mutate(scenario = "뷔트리시란 Q3M SC")
    return(out)
  }

  if (treatment == "inotersen") {
    n_d <- floor(end_days / 7) + 1
    ino_ev <- data.frame(
      time = seq(0, by = 7, length.out = n_d),
      amt = 1.0, cmt = "INO_EFF", evid = 1, rate = 0, ID = 1
    )
    out <- mod_run %>%
      mrgsim_df(data = ino_ev, end = end_days, delta = delta_days,
                add = list(ID = 1)) %>%
      mutate(scenario = "이노테르센 QW SC")
    return(out)
  }

  # 자연 경과 또는 타파미디스 단독
  out <- if (!is.null(events)) {
    mod_run %>%
      mrgsim(events = events, end = end_days, delta = delta_days) %>%
      as.data.frame()
  } else {
    mod_run %>%
      mrgsim(end = end_days, delta = delta_days) %>%
      as.data.frame()
  }

  out %>% mutate(scenario = switch(treatment,
    "none"       = "무치료 (자연 경과)",
    "tafamidis"  = "타파미디스 61mg QD",
    "combo"      = "타파미디스 + HF 관리",
    "unknown"
  ))
}

# ============================================================
# UI 정의
# ============================================================
ui <- dashboardPage(
  skin = "blue",
  dashboardHeader(title = "ATTR 아밀로이드증 QSP"),

  dashboardSidebar(
    sidebarMenu(
      menuItem("환자 프로파일",       tabName = "tab1", icon = icon("user")),
      menuItem("TTR PK / 응집 경로",  tabName = "tab2", icon = icon("dna")),
      menuItem("심장 ATTR",           tabName = "tab3", icon = icon("heartbeat")),
      menuItem("신경 ATTR (FAP)",     tabName = "tab4", icon = icon("brain")),
      menuItem("치료 시나리오 비교",   tabName = "tab5", icon = icon("pills")),
      menuItem("바이오마커 추적",      tabName = "tab6", icon = icon("chart-line")),
      menuItem("민감도 분석",          tabName = "tab7", icon = icon("sliders-h")),
      menuItem("가상 환자 집단",       tabName = "tab8", icon = icon("users"))
    ),
    hr(),
    # 공통 시뮬레이션 설정
    sliderInput("sim_years", "시뮬레이션 기간 (년)", 1, 10, 5, step = 0.5),
    selectInput("mut_type", "TTR 변이 유형",
      choices = c("야생형 (ATTRwt)" = "1.0",
                  "Val122I (심근병증형)" = "1.5",
                  "Val30M (다발신경병증형)" = "2.0",
                  "고위험 변이" = "3.0"),
      selected = "1.0"
    ),
    selectInput("treatment_sel", "치료제 선택",
      choices = c("무치료 (자연 경과)" = "none",
                  "타파미디스 61mg QD" = "tafamidis",
                  "파티시란 Q3W IV" = "patisiran",
                  "뷔트리시란 Q3M SC" = "vutrisiran",
                  "이노테르센 QW SC" = "inotersen",
                  "타파미디스 + HF 관리" = "combo"),
      selected = "tafamidis"
    ),
    actionButton("run_sim", "시뮬레이션 실행", icon = icon("play"),
                 style = "background-color:#27AE60; color:white; width:90%; margin:5px")
  ),

  dashboardBody(
    tags$head(tags$style(HTML("
      .box-header { background-color: #2C3E50 !important; color: white !important; }
      .content-wrapper { background-color: #ECF0F1; }
    "))),

    tabItems(
      # =========================================================
      # Tab 1: 환자 프로파일
      # =========================================================
      tabItem(tabName = "tab1",
        fluidRow(
          box(title = "환자 특성 설정", status = "primary", solidHeader = TRUE,
              width = 4,
              numericInput("age", "나이 (세)", 72, 50, 90),
              selectInput("sex", "성별", choices = c("남성" = "M", "여성" = "F")),
              numericInput("body_wt", "체중 (kg)", 70, 40, 120),
              selectInput("ethnicity", "인종",
                choices = c("아시아인" = "asian", "백인" = "white",
                            "아프리카계" = "african", "기타" = "other")),
              selectInput("diagnosis_type", "진단 유형",
                choices = c("ATTRwt 심근병증" = "attr_wt_card",
                            "ATTRv V30M 다발신경병증" = "attr_v30m",
                            "ATTRv V122I 심근병증" = "attr_v122i")),
              numericInput("diag_lv_thick", "진단시 LV 벽두께 (mm)", 14, 8, 30),
              numericInput("diag_bnp", "진단시 NT-proBNP (pg/mL)", 1500, 50, 20000),
              numericInput("diag_6mwt", "진단시 6MWT 거리 (m)", 380, 50, 600)
          ),
          box(title = "환자 프로파일 요약", status = "info", solidHeader = TRUE,
              width = 4,
              verbatimTextOutput("patient_summary"),
              hr(),
              valueBoxOutput("vbox_nac_stage", width = 12),
              valueBoxOutput("vbox_survival_est", width = 12)
          ),
          box(title = "ATTR 역학 및 임상 경로", status = "warning",
              solidHeader = TRUE, width = 4,
              h4("ATTRwt 심근병증"),
              p("• 남성에서 더 흔함 (남:여 = 4:1)"),
              p("• 진단 중앙 연령: 75세"),
              p("• 진단 후 중앙 생존기간: 2.5-5년"),
              p("• 수근관 증후군 先行: ~30%"),
              hr(),
              h4("ATTRv 다발신경병증 (FAP)"),
              p("• V30M: 포르투갈·스웨덴·일본 多"),
              p("• 발병 연령: 25-45세 (조기) / 50-70세 (만기)"),
              p("• V122I: 아프리카계 미국인 3.4% 보인자"),
              hr(),
              h4("NAC 병기 분류"),
              p("1기: NT-proBNP <3000 & TropT <50 ng/L"),
              p("2기: 하나만 초과"),
              p("3기: 둘 다 초과 → 중앙 생존 <2.5년")
          )
        ),
        fluidRow(
          box(title = "진단 알고리즘 요약", status = "success", solidHeader = TRUE,
              width = 12,
              p("1. 임상 의심: 노인 HFpEF + LV 비후 / 양측성 수근관 증후군 / 다발신경병증"),
              p("2. 초음파/CMR: LV 벽두께 ≥12mm + 심첨부 보존 GLS 패턴 + ECV↑"),
              p("3. Tc-99m DPD/PYP 신티그래피: Grade 2-3 → 생검 없이 진단 확정 (혈청 SPEP/UPEP 음성 전제)"),
              p("4. TTR 유전자형 검사: ATTRwt vs ATTRv 구분 → 치료제 선택에 중요"),
              p("5. 치료 시작: 타파미디스 (심장형) / 파티시란·뷔트리시란 (신경병증형)")
          )
        )
      ),

      # =========================================================
      # Tab 2: TTR PK / 응집
      # =========================================================
      tabItem(tabName = "tab2",
        fluidRow(
          box(title = "TTR 농도 시뮬레이션 파라미터", status = "primary",
              solidHeader = TRUE, width = 3,
              sliderInput("kprod", "TTR 생산율 (kprod, μg/mL/day)",
                          50, 200, 104, step = 5),
              sliderInput("mut_factor_2", "사량체 불안정성\n(돌연변이 인자)",
                          0.5, 4.0, 1.0, step = 0.1),
              sliderInput("kdiss", "기저 해리 속도 (×10⁻⁴/day)",
                          1, 20, 5, step = 1),
              hr(),
              h5("타파미디스 PD 파라미터"),
              sliderInput("ec50_taf2", "EC50 (μg/mL)", 0.5, 10, 2, step = 0.5),
              sliderInput("emax_taf2", "Emax (최대 억제율)", 0.3, 1.0, 0.8, step = 0.05),
              actionButton("run_tab2", "업데이트", icon = icon("refresh"),
                           style = "background-color:#3498DB; color:white")
          ),
          box(title = "혈장 TTR 농도 변화", status = "info", solidHeader = TRUE,
              width = 9,
              plotlyOutput("plt_ttr_pk", height = "350px"),
              hr(),
              plotlyOutput("plt_aggregation", height = "300px")
          )
        ),
        fluidRow(
          box(title = "TTR 응집 연쇄반응 현황", status = "warning",
              solidHeader = TRUE, width = 6,
              plotlyOutput("plt_oligomers", height = "280px")
          ),
          box(title = "약물 억제율 (시간-효과)", status = "success",
              solidHeader = TRUE, width = 6,
              plotlyOutput("plt_inhibition", height = "280px")
          )
        )
      ),

      # =========================================================
      # Tab 3: 심장 ATTR
      # =========================================================
      tabItem(tabName = "tab3",
        fluidRow(
          valueBoxOutput("vbox_lv_thick"),
          valueBoxOutput("vbox_bnp"),
          valueBoxOutput("vbox_trop")
        ),
        fluidRow(
          box(title = "LV 벽두께 및 심장 아밀로이드", status = "primary",
              solidHeader = TRUE, width = 6,
              plotlyOutput("plt_lv_cardiac", height = "350px")
          ),
          box(title = "NT-proBNP 및 트로포닌 추이", status = "danger",
              solidHeader = TRUE, width = 6,
              plotlyOutput("plt_bnp_trop", height = "350px")
          )
        ),
        fluidRow(
          box(title = "심장 기능 지수 (LVEF 대리 지표)", status = "info",
              solidHeader = TRUE, width = 6,
              plotlyOutput("plt_card_func", height = "300px")
          ),
          box(title = "ATTR-ACT 임상시험 결과 참조", status = "success",
              solidHeader = TRUE, width = 6,
              h4("ATTR-ACT (Maurer et al., NEJM 2018)"),
              p("• 대상: ATTRwt 심근병증 (N=441)"),
              p("• 치료: 타파미디스 61mg 또는 20mg QD vs 위약"),
              p("• 기간: 30개월 추적"),
              p("• 주요 결과:"),
              tags$ul(
                tags$li("모든 원인 사망: 타파미디스 29.5% vs 위약 42.9% (HR 0.70)"),
                tags$li("CV 입원: 타파미디스 0.48 vs 위약 0.70회/년 (RR 0.68)"),
                tags$li("6MWT 거리 감소 완화 (위약 대비)"),
                tags$li("KCCQ-OS 감소 완화 (위약 대비)")
              ),
              p("• 안전성: 위약과 유사 (CYP 독립 대사)"),
              hr(),
              h5("FDA 승인 (2019) 및 EMA 승인 (2019)"),
              p("Tafamidis meglumine 61mg (Vyndaqel®)"),
              p("Tafamidis 61mg (Vyndamax® — 단일 캡슐)")
          )
        )
      ),

      # =========================================================
      # Tab 4: 신경 ATTR
      # =========================================================
      tabItem(tabName = "tab4",
        fluidRow(
          box(title = "신경 병리 파라미터 설정", status = "primary",
              solidHeader = TRUE, width = 3,
              sliderInput("nis_scale", "NIS 증가 속도 (NIS/AU)", 1, 15, 5, step = 0.5),
              sliderInput("nerve_kdiss", "신경 침착 분율", 0.1, 0.6, 0.35, step = 0.05),
              selectInput("neuro_treat", "신경병증 치료제",
                choices = c("무치료" = "none",
                            "파티시란 Q3W" = "patisiran",
                            "뷔트리시란 Q3M" = "vutrisiran",
                            "이노테르센 QW" = "inotersen")),
              actionButton("run_tab4", "신경 시뮬레이션 실행",
                           style = "background-color:#8E44AD; color:white")
          ),
          box(title = "NIS / mNIS+7 점수 추이", status = "warning",
              solidHeader = TRUE, width = 9,
              plotlyOutput("plt_nis", height = "350px"),
              hr(),
              plotlyOutput("plt_auto_np", height = "250px")
          )
        ),
        fluidRow(
          box(title = "FAP 임상시험 결과 비교 (mNIS+7)", status = "info",
              solidHeader = TRUE, width = 6,
              plotOutput("plt_fap_trials", height = "320px")
          ),
          box(title = "임상시험 데이터 요약", status = "success",
              solidHeader = TRUE, width = 6,
              h4("주요 FAP 임상시험 결과"),
              DTOutput("tbl_fap_trials", height = "280px")
          )
        )
      ),

      # =========================================================
      # Tab 5: 치료 시나리오 비교
      # =========================================================
      tabItem(tabName = "tab5",
        fluidRow(
          box(title = "비교할 시나리오 선택", status = "primary",
              solidHeader = TRUE, width = 3,
              checkboxGroupInput("compare_treatments", "치료 시나리오",
                choices = c("무치료 (자연 경과)" = "none",
                            "타파미디스 61mg QD" = "tafamidis",
                            "파티시란 Q3W IV" = "patisiran",
                            "뷔트리시란 Q3M SC" = "vutrisiran",
                            "이노테르센 QW SC" = "inotersen",
                            "타파미디스 + HF 관리" = "combo"),
                selected = c("none", "tafamidis", "patisiran", "vutrisiran")
              ),
              selectInput("compare_endpoint", "비교 지표",
                choices = c("혈장 TTR (mg/dL)" = "TTR_mgdL",
                            "심장 아밀로이드 부담" = "CardAmyloid",
                            "신경 아밀로이드 부담" = "NerveAmyloid",
                            "NT-proBNP (pg/mL)" = "NT_proBNP",
                            "NIS 점수" = "NIS_score",
                            "6분보행거리 (m)" = "Walk6min",
                            "KCCQ-OS" = "KCCQ",
                            "심장 기능 지수" = "CardFunc"),
                selected = "CardAmyloid"
              ),
              actionButton("run_compare", "비교 실행",
                           style = "background-color:#E74C3C; color:white; width:90%")
          ),
          box(title = "치료 시나리오 비교 그래프", status = "danger",
              solidHeader = TRUE, width = 9,
              plotlyOutput("plt_compare", height = "450px")
          )
        ),
        fluidRow(
          box(title = "시뮬레이션 말기 요약 표", status = "info",
              solidHeader = TRUE, width = 12,
              DTOutput("tbl_compare_summary")
          )
        )
      ),

      # =========================================================
      # Tab 6: 바이오마커 추적
      # =========================================================
      tabItem(tabName = "tab6",
        fluidRow(
          box(title = "바이오마커 대시보드", status = "primary",
              solidHeader = TRUE, width = 8,
              plotlyOutput("plt_biomarker_panel", height = "600px")
          ),
          box(title = "바이오마커 임계값 알림", status = "warning",
              solidHeader = TRUE, width = 4,
              h4("NT-proBNP 상태"),
              verbatimTextOutput("biomarker_bnp_status"),
              h4("트로포닌T 상태"),
              verbatimTextOutput("biomarker_trop_status"),
              h4("LV 벽두께"),
              verbatimTextOutput("biomarker_lv_status"),
              hr(),
              h4("NAC 병기 변화"),
              plotlyOutput("plt_nac_stage", height = "200px")
          )
        )
      ),

      # =========================================================
      # Tab 7: 민감도 분석
      # =========================================================
      tabItem(tabName = "tab7",
        fluidRow(
          box(title = "민감도 분석 파라미터", status = "primary",
              solidHeader = TRUE, width = 3,
              selectInput("sens_param", "분석 파라미터",
                choices = c("타파미디스 EC50" = "EC50_TAF",
                            "타파미디스 Emax" = "Emax_TAF",
                            "사량체 해리 속도" = "kdiss_base",
                            "섬유화 속도" = "kfib",
                            "아밀로이드 제거율" = "kclear_amyl",
                            "TTR 생산율" = "kprod_TTR"),
                selected = "EC50_TAF"
              ),
              numericInput("sens_low",  "최솟값", 0.5),
              numericInput("sens_high", "최댓값", 8.0),
              numericInput("sens_n",    "분석 단계 수", 5, 3, 10, step = 1),
              selectInput("sens_endpoint", "민감도 결과 지표",
                choices = c("심장 아밀로이드" = "CardAmyloid",
                            "NT-proBNP" = "NT_proBNP",
                            "KCCQ" = "KCCQ",
                            "LV 벽두께" = "LV_wall"),
                selected = "CardAmyloid"
              ),
              actionButton("run_sens", "민감도 분석 실행",
                           style = "background-color:#F39C12; color:white; width:90%")
          ),
          box(title = "민감도 분석 결과", status = "warning",
              solidHeader = TRUE, width = 9,
              plotlyOutput("plt_sensitivity", height = "420px"),
              hr(),
              plotlyOutput("plt_tornado", height = "200px")
          )
        )
      ),

      # =========================================================
      # Tab 8: 가상 환자 집단
      # =========================================================
      tabItem(tabName = "tab8",
        fluidRow(
          box(title = "가상 환자 집단 설정", status = "primary",
              solidHeader = TRUE, width = 3,
              numericInput("vp_n", "환자 수 (N)", 50, 10, 200, step = 10),
              sliderInput("vp_cv", "파라미터 변이도 (CV%)", 5, 50, 20, step = 5),
              selectInput("vp_treatment", "치료제",
                choices = c("무치료" = "none",
                            "타파미디스 61mg QD" = "tafamidis",
                            "파티시란 Q3W IV" = "patisiran"),
                selected = "tafamidis"
              ),
              selectInput("vp_endpoint", "결과 지표",
                choices = c("심장 아밀로이드" = "CardAmyloid",
                            "NT-proBNP" = "NT_proBNP",
                            "NIS 점수" = "NIS_score",
                            "6분보행거리" = "Walk6min",
                            "KCCQ" = "KCCQ"),
                selected = "NT_proBNP"
              ),
              numericInput("vp_seed", "난수 씨앗값", 42),
              actionButton("run_vp", "가상 환자 집단 실행",
                           style = "background-color:#8E44AD; color:white; width:90%")
          ),
          box(title = "가상 환자 집단 시뮬레이션 결과", status = "info",
              solidHeader = TRUE, width = 9,
              plotlyOutput("plt_vp_spaghetti", height = "350px"),
              hr(),
              plotlyOutput("plt_vp_dist", height = "280px")
          )
        ),
        fluidRow(
          box(title = "집단 예측 구간 (VPC)", status = "success",
              solidHeader = TRUE, width = 8,
              plotlyOutput("plt_vpc", height = "320px")
          ),
          box(title = "NAC 병기 분포 (말기)", status = "warning",
              solidHeader = TRUE, width = 4,
              plotlyOutput("plt_nac_dist", height = "300px")
          )
        )
      )
    ) # end tabItems
  ) # end dashboardBody
) # end dashboardPage

# ============================================================
# Server 정의
# ============================================================
server <- function(input, output, session) {

  # 반응형 파라미터 목록
  current_params <- reactive({
    list(
      mut_factor = as.numeric(input$mut_type),
      kprod_TTR  = 104,
      kdiss_base = 5e-4
    )
  })

  # 주 시뮬레이션 (공통)
  sim_result <- eventReactive(input$run_sim, {
    req(mod_shiny)
    withProgress(message = "시뮬레이션 실행 중...", {
      run_sim(
        mod    = mod_shiny,
        params = current_params(),
        treatment = input$treatment_sel,
        end_years = input$sim_years
      )
    })
  }, ignoreNULL = FALSE)

  # ===== Tab 1: 환자 프로파일 =====
  output$patient_summary <- renderText({
    paste0(
      "나이: ", input$age, "세 | 성별: ", input$sex, "\n",
      "체중: ", input$body_wt, " kg | 인종: ", input$ethnicity, "\n",
      "진단 유형: ", input$diagnosis_type, "\n",
      "LV 벽두께: ", input$diag_lv_thick, " mm\n",
      "NT-proBNP: ", input$diag_bnp, " pg/mL\n",
      "6MWT 거리: ", input$diag_6mwt, " m\n",
      "TTR 변이: ", input$mut_type
    )
  })

  output$vbox_nac_stage <- renderValueBox({
    stage <- if (input$diag_bnp >= 3000 && input$diag_lv_thick > 15) 3 else
             if (input$diag_bnp >= 3000 || input$diag_lv_thick > 15) 2 else 1
    col <- c("green", "yellow", "red")[stage]
    valueBox(paste("NAC", stage, "기"), "추정 병기", icon = icon("stethoscope"),
             color = col)
  })

  output$vbox_survival_est <- renderValueBox({
    stage <- if (input$diag_bnp >= 3000 && input$diag_lv_thick > 15) 3 else
             if (input$diag_bnp >= 3000 || input$diag_lv_thick > 15) 2 else 1
    surv <- c("4-5년", "2.5-4년", "<2.5년")[stage]
    valueBox(surv, "추정 중앙 생존 (무치료)", icon = icon("clock"),
             color = c("green","yellow","red")[stage])
  })

  # ===== Tab 2: TTR PK/응집 =====
  sim_tab2 <- eventReactive(input$run_tab2, {
    req(mod_shiny)
    p_list <- list(
      kprod_TTR  = input$kprod,
      mut_factor = input$mut_factor_2,
      kdiss_base = input$kdiss * 1e-4,
      EC50_TAF   = input$ec50_taf2,
      Emax_TAF   = input$emax_taf2
    )
    taf_ev <- ev(time = 0, amt = 61, cmt = "TAF_GUT", evid = 1,
                 rate = 0, ii = 1, addl = 364 * 5)
    list(
      treated = mod_shiny %>% param(p_list) %>%
        mrgsim(events = taf_ev, end = 365*5, delta = 7) %>%
        as.data.frame() %>% mutate(arm = "타파미디스"),
      control = mod_shiny %>% param(p_list) %>%
        mrgsim(end = 365*5, delta = 7) %>%
        as.data.frame() %>% mutate(arm = "무치료")
    )
  }, ignoreNULL = FALSE)

  output$plt_ttr_pk <- renderPlotly({
    d <- sim_tab2()
    if (is.null(d)) return(NULL)
    df <- bind_rows(d) %>% mutate(years = time/365.25)
    p <- ggplot(df, aes(years, TTR_mgdL, color = arm)) +
      geom_line(linewidth = 0.9) +
      geom_hline(yintercept = c(20, 40), linetype = "dashed", color = "gray") +
      labs(title = "혈장 TTR 농도", x = "시간 (년)", y = "TTR (mg/dL)",
           color = "치료") +
      theme_bw()
    ggplotly(p)
  })

  output$plt_aggregation <- renderPlotly({
    d <- sim_tab2()
    if (is.null(d)) return(NULL)
    df <- bind_rows(d) %>% mutate(years = time/365.25)
    p <- ggplot(df, aes(years, CardAmyloid, color = arm)) +
      geom_line(linewidth = 0.9) +
      labs(title = "심장 아밀로이드 침착 (AU)", x = "시간 (년)", y = "아밀로이드 (AU)") +
      theme_bw()
    ggplotly(p)
  })

  output$plt_oligomers <- renderPlotly({
    d <- sim_tab2()
    if (is.null(d)) return(NULL)
    df <- bind_rows(d) %>% mutate(years = time/365.25)
    # melt manually
    df2 <- df %>% select(years, arm, TTR_C, NIS_score) %>%
      pivot_longer(c(TTR_C, NIS_score), names_to = "marker", values_to = "val")
    p <- ggplot(df2, aes(years, val, color = arm, linetype = marker)) +
      geom_line() +
      labs(title = "TTR 응집 지표", x = "시간 (년)", y = "값") +
      theme_bw()
    ggplotly(p)
  })

  output$plt_inhibition <- renderPlotly({
    d <- sim_tab2()
    if (is.null(d)) return(NULL)
    df <- bind_rows(d) %>%
      filter(arm == "타파미디스") %>%
      mutate(years = time/365.25, inh_pct = Inh_TAF_pct)
    p <- ggplot(df, aes(years, inh_pct)) +
      geom_line(color = "darkgreen", linewidth = 1) +
      scale_y_continuous(limits = c(0, 100)) +
      labs(title = "타파미디스 TTR 사량체 해리 억제율",
           x = "시간 (년)", y = "억제율 (%)") +
      theme_bw()
    ggplotly(p)
  })

  # ===== Tab 3: 심장 ATTR =====
  output$vbox_lv_thick <- renderValueBox({
    d <- sim_result()
    val <- if (!is.null(d)) round(tail(d$LV_wall, 1), 1) else 10
    valueBox(paste0(val, " mm"), "LV 벽두께 (말기)", icon = icon("heartbeat"),
             color = if (val > 15) "red" else if (val > 12) "yellow" else "green")
  })

  output$vbox_bnp <- renderValueBox({
    d <- sim_result()
    val <- if (!is.null(d)) round(tail(d$NT_proBNP, 1), 0) else 80
    valueBox(paste0(val, " pg/mL"), "NT-proBNP (말기)", icon = icon("chart-bar"),
             color = if (val > 3000) "red" else if (val > 900) "yellow" else "green")
  })

  output$vbox_trop <- renderValueBox({
    d <- sim_result()
    val <- if (!is.null(d)) round(tail(d$TroponinT, 1), 1) else 5
    valueBox(paste0(val, " ng/L"), "TroponinT (말기)", icon = icon("plus"),
             color = if (val > 50) "red" else if (val > 20) "yellow" else "green")
  })

  output$plt_lv_cardiac <- renderPlotly({
    d <- sim_result()
    if (is.null(d)) return(NULL)
    d <- d %>% mutate(years = time/365.25)
    p <- ggplot(d, aes(years)) +
      geom_line(aes(y = LV_wall, color = "LV 벽두께 (mm)"), linewidth = 1) +
      geom_line(aes(y = CardAmyloid * 5, color = "심장 아밀로이드 (AU×5)"), linewidth = 1) +
      geom_hline(yintercept = 12, linetype = "dashed", color = "blue") +
      scale_color_manual(values = c("red", "darkblue")) +
      labs(title = "심장 ATTR 부담", x = "시간 (년)", y = "값", color = NULL) +
      theme_bw() + theme(legend.position = "top")
    ggplotly(p)
  })

  output$plt_bnp_trop <- renderPlotly({
    d <- sim_result()
    if (is.null(d)) return(NULL)
    d2 <- d %>% mutate(years = time/365.25) %>%
      select(years, NT_proBNP, TroponinT) %>%
      pivot_longer(c(NT_proBNP, TroponinT), names_to = "marker", values_to = "val")
    p <- ggplot(d2, aes(years, val, color = marker)) +
      geom_line(linewidth = 1) +
      geom_hline(data = data.frame(marker = "NT_proBNP", thresh = 3000),
                 aes(yintercept = thresh), linetype = "dashed", color = "red") +
      labs(title = "심장 손상 바이오마커", x = "시간 (년)", y = "값", color = NULL) +
      theme_bw() + theme(legend.position = "top")
    ggplotly(p)
  })

  output$plt_card_func <- renderPlotly({
    d <- sim_result()
    if (is.null(d)) return(NULL)
    d <- d %>% mutate(years = time/365.25)
    p <- ggplot(d, aes(years, CardFunc)) +
      geom_line(color = "steelblue", linewidth = 1.1) +
      geom_hline(yintercept = 0.5, linetype = "dashed", color = "orange") +
      scale_y_continuous(limits = c(0, 1.1)) +
      labs(title = "심장 기능 지수 (1=정상)", x = "시간 (년)", y = "기능 지수") +
      theme_bw()
    ggplotly(p)
  })

  # ===== Tab 4: 신경 ATTR =====
  sim_neuro <- eventReactive(input$run_tab4, {
    req(mod_shiny)
    p_list <- list(
      mut_factor = as.numeric(input$mut_type),
      NIS_scale  = input$nis_scale,
      frac_nerve = input$nerve_kdiss
    )
    run_sim(mod_shiny, p_list, input$neuro_treat, input$sim_years)
  })

  output$plt_nis <- renderPlotly({
    d <- sim_neuro()
    if (is.null(d)) return(NULL)
    d <- d %>% mutate(years = time/365.25)
    p <- ggplot(d, aes(years, NIS_score)) +
      geom_line(color = "purple", linewidth = 1.1) +
      geom_hline(yintercept = 80, linetype = "dashed", color = "red") +
      annotate("text", x = 0.3, y = 85, label = "장애 임계 (80)", color = "red") +
      labs(title = "NIS 점수 (신경병증 장애 지수)", x = "시간 (년)", y = "NIS 점수") +
      theme_bw()
    ggplotly(p)
  })

  output$plt_auto_np <- renderPlotly({
    d <- sim_neuro()
    if (is.null(d)) return(NULL)
    d <- d %>% mutate(years = time/365.25)
    p <- ggplot(d, aes(years, AUTO_NP)) +
      geom_line(color = "darkorange", linewidth = 1) +
      scale_y_continuous(limits = c(0, 1)) +
      labs(title = "자율신경 병증 지수", x = "시간 (년)", y = "지수 (0=정상, 1=중증)") +
      theme_bw()
    ggplotly(p)
  })

  output$plt_fap_trials <- renderPlot({
    trial_data <- data.frame(
      trial = rep(c("APOLLO\n(파티시란)", "HELIOS-A\n(뷔트리시란)",
                    "NEURO-TTR\n(이노테르센)"), each = 2),
      arm   = rep(c("치료군", "위약"), 3),
      delta_mNIS7 = c(-0.9, 14.3, -0.7, 17.8, -6.7, 3.9)
    )
    ggplot(trial_data, aes(arm, delta_mNIS7, fill = arm)) +
      geom_col(width = 0.6) +
      geom_hline(yintercept = 0, color = "black") +
      facet_wrap(~trial) +
      scale_fill_manual(values = c("#27AE60", "#E74C3C")) +
      labs(title = "FAP 임상시험 mNIS+7 변화량 (18개월)",
           y = "mNIS+7 변화 (점수 증가 = 악화)", x = NULL, fill = NULL) +
      theme_bw() + theme(legend.position = "bottom")
  })

  output$tbl_fap_trials <- renderDT({
    data.frame(
      임상시험 = c("APOLLO", "HELIOS-A", "NEURO-TTR", "NEURO-TTRansform"),
      치료제    = c("파티시란 Q3W IV", "뷔트리시란 Q3M SC",
                   "이노테르센 QW SC", "에플론테르센 Q4W SC"),
      TTR감소   = c("80%", "87%", "75%", "82%"),
      mNIS7     = c("−0.9 vs +14.3*", "−0.7 vs +17.8*",
                    "−6.7 vs +3.9*", "−9.7 vs +15.6*"),
      기간      = c("18개월", "9개월", "15개월", "15개월"),
      연도      = c(2018, 2022, 2018, 2024)
    )
  }, options = list(dom = 't', pageLength = 10), rownames = FALSE)

  # ===== Tab 5: 비교 =====
  compare_results <- eventReactive(input$run_compare, {
    req(mod_shiny)
    treatments <- input$compare_treatments
    if (length(treatments) == 0) return(NULL)

    results_list <- lapply(treatments, function(trt) {
      tryCatch(
        run_sim(mod_shiny, current_params(), trt, input$sim_years),
        error = function(e) NULL
      )
    })
    results_list <- results_list[!sapply(results_list, is.null)]
    bind_rows(results_list) %>% mutate(years = time/365.25)
  })

  output$plt_compare <- renderPlotly({
    d <- compare_results()
    if (is.null(d)) return(NULL)
    ep <- input$compare_endpoint
    p <- ggplot(d, aes_string("years", ep, color = "scenario")) +
      geom_line(linewidth = 0.9) +
      labs(title = paste("치료 시나리오 비교 —", ep),
           x = "시간 (년)", y = ep, color = "치료") +
      theme_bw() + theme(legend.position = "bottom")
    ggplotly(p)
  })

  output$tbl_compare_summary <- renderDT({
    d <- compare_results()
    if (is.null(d)) return(NULL)
    d %>%
      group_by(scenario) %>%
      filter(time == max(time)) %>%
      slice(1) %>%
      ungroup() %>%
      select(시나리오 = scenario, TTR = TTR_mgdL, 심장아밀 = CardAmyloid,
             NT_proBNP, NIS = NIS_score, 보행거리 = Walk6min,
             KCCQ, NAC병기 = NAC_stage) %>%
      mutate(across(where(is.numeric), ~round(.x, 1)))
  }, options = list(pageLength = 10), rownames = FALSE)

  # ===== Tab 6: 바이오마커 =====
  output$plt_biomarker_panel <- renderPlotly({
    d <- sim_result()
    if (is.null(d)) return(NULL)
    d <- d %>% mutate(years = time/365.25)

    p1 <- ggplot(d, aes(years, NT_proBNP)) +
      geom_line(color="red") + geom_hline(yintercept=3000, linetype="dashed") +
      labs(title="NT-proBNP", x="", y="pg/mL") + theme_bw(base_size=9)
    p2 <- ggplot(d, aes(years, TroponinT)) +
      geom_line(color="darkred") + geom_hline(yintercept=50, linetype="dashed") +
      labs(title="TroponinT", x="", y="ng/L") + theme_bw(base_size=9)
    p3 <- ggplot(d, aes(years, LV_wall)) +
      geom_line(color="navy") + geom_hline(yintercept=12, linetype="dashed") +
      labs(title="LV 벽두께", x="", y="mm") + theme_bw(base_size=9)
    p4 <- ggplot(d, aes(years, Walk6min)) +
      geom_line(color="forestgreen") +
      labs(title="6MWT", x="시간(년)", y="m") + theme_bw(base_size=9)
    p5 <- ggplot(d, aes(years, KCCQ)) +
      geom_line(color="purple") + ylim(0,100) +
      labs(title="KCCQ-OS", x="시간(년)", y="점수") + theme_bw(base_size=9)
    p6 <- ggplot(d, aes(years, TTR_mgdL)) +
      geom_line(color="steelblue") +
      labs(title="혈장 TTR", x="시간(년)", y="mg/dL") + theme_bw(base_size=9)

    subplot(
      ggplotly(p1), ggplotly(p2),
      ggplotly(p3), ggplotly(p4),
      ggplotly(p5), ggplotly(p6),
      nrows = 3, shareX = FALSE
    )
  })

  output$biomarker_bnp_status <- renderText({
    d <- sim_result()
    val <- if (!is.null(d)) round(tail(d$NT_proBNP, 1), 0) else 80
    if (val >= 3000) paste0(val, " pg/mL — ⚠ NAC 3기 위험") else
    if (val >= 900) paste0(val, " pg/mL — ⚡ 중등도 상승") else
    paste0(val, " pg/mL — ✓ 정상 범위 내")
  })

  output$biomarker_trop_status <- renderText({
    d <- sim_result()
    val <- if (!is.null(d)) round(tail(d$TroponinT, 1), 1) else 5
    if (val >= 50) paste0(val, " ng/L — ⚠ 심근 손상 지속") else
    if (val >= 20) paste0(val, " ng/L — ⚡ 경미한 상승") else
    paste0(val, " ng/L — ✓ 정상")
  })

  output$biomarker_lv_status <- renderText({
    d <- sim_result()
    val <- if (!is.null(d)) round(tail(d$LV_wall, 1), 1) else 10
    if (val > 16) paste0(val, " mm — ⚠ 중증 비후") else
    if (val > 12) paste0(val, " mm — ⚡ 경계 비후") else
    paste0(val, " mm — ✓ 정상")
  })

  output$plt_nac_stage <- renderPlotly({
    d <- sim_result()
    if (is.null(d)) return(NULL)
    d <- d %>% mutate(years = time/365.25)
    p <- ggplot(d, aes(years, NAC_stage)) +
      geom_step(color = "darkorange", linewidth = 1.2) +
      scale_y_continuous(breaks = 1:3, limits = c(0.5, 3.5)) +
      labs(title = "NAC 병기", x = "시간 (년)", y = "병기") +
      theme_bw()
    ggplotly(p)
  })

  # ===== Tab 7: 민감도 =====
  sensitivity_data <- eventReactive(input$run_sens, {
    req(mod_shiny)
    param_name <- input$sens_param
    vals <- seq(input$sens_low, input$sens_high, length.out = input$sens_n)
    ep   <- input$sens_endpoint

    taf_ev <- ev(time = 0, amt = 61, cmt = "TAF_GUT", evid = 1,
                 rate = 0, ii = 1, addl = as.integer(input$sim_years * 365) - 1)

    res_list <- lapply(vals, function(v) {
      p_update <- setNames(list(v), param_name)
      tryCatch({
        mod_shiny %>%
          param(p_update) %>%
          mrgsim(events = taf_ev, end = input$sim_years * 365, delta = 14) %>%
          as.data.frame() %>%
          mutate(param_val = v, years = time / 365.25)
      }, error = function(e) NULL)
    })
    bind_rows(res_list[!sapply(res_list, is.null)])
  })

  output$plt_sensitivity <- renderPlotly({
    d <- sensitivity_data()
    if (is.null(d)) return(NULL)
    ep <- input$sens_endpoint
    p <- ggplot(d, aes_string("years", ep, color = "factor(param_val)")) +
      geom_line(linewidth = 0.8) +
      labs(title = paste("민감도 분석:", input$sens_param, "→", ep),
           x = "시간 (년)", y = ep, color = input$sens_param) +
      theme_bw()
    ggplotly(p)
  })

  output$plt_tornado <- renderPlotly({
    d <- sensitivity_data()
    if (is.null(d)) return(NULL)
    ep <- input$sens_endpoint
    tornado_df <- d %>%
      filter(abs(time - max(time)) < 15) %>%
      group_by(param_val) %>%
      slice(1) %>%
      ungroup()
    p <- ggplot(tornado_df, aes_string("param_val", ep)) +
      geom_point(size = 3, color = "steelblue") +
      geom_line(color = "steelblue") +
      labs(title = paste("말기", ep, "vs", input$sens_param),
           x = input$sens_param, y = ep) +
      theme_bw()
    ggplotly(p)
  })

  # ===== Tab 8: 가상 환자 집단 =====
  vp_data <- eventReactive(input$run_vp, {
    req(mod_shiny)
    n <- input$vp_n
    cv <- input$vp_cv / 100
    set.seed(input$vp_seed)

    pop <- data.frame(
      ID         = 1:n,
      kprod_TTR  = rlnorm(n, log(104),   cv),
      kel_TTR    = rlnorm(n, log(0.347), cv * 0.7),
      mut_factor = rlnorm(n, log(as.numeric(input$mut_type)), cv * 0.5),
      kdiss_base = rlnorm(n, log(5e-4),  cv),
      kfib       = rlnorm(n, log(0.02),  cv)
    )

    # 투여 데이터
    dose_df <- if (input$vp_treatment == "tafamidis") {
      do.call(rbind, lapply(1:n, function(i)
        data.frame(ID=i, time=0, amt=61, cmt="TAF_GUT", evid=1, rate=0,
                   ii=1, addl=as.integer(input$sim_years*365)-1)
      ))
    } else if (input$vp_treatment == "patisiran") {
      nd <- floor(input$sim_years*365/21)+1
      do.call(rbind, lapply(1:n, function(i)
        data.frame(ID=i, time=seq(0,by=21,length.out=nd),
                   amt=1, cmt="PAT_EFF", evid=1, rate=0)
      ))
    } else {
      data.frame(ID=1:n, time=0, amt=0, cmt="TAF_GUT", evid=0, rate=0)
    }

    out <- mod_shiny %>%
      data_set(dose_df) %>%
      idata_set(pop) %>%
      mrgsim(end = input$sim_years * 365, delta = 30) %>%
      as.data.frame() %>%
      mutate(years = time / 365.25)
    out
  })

  output$plt_vp_spaghetti <- renderPlotly({
    d <- vp_data()
    if (is.null(d)) return(NULL)
    ep <- input$vp_endpoint
    p <- ggplot(d, aes_string("years", ep, group = "ID")) +
      geom_line(alpha = 0.25, color = "steelblue") +
      stat_summary(aes(group = 1), fun = median, geom = "line",
                   color = "red", linewidth = 1.3) +
      labs(title = paste("가상 환자 집단 (N=", input$vp_n, ") —", ep),
           subtitle = "파란 선: 개별 환자 | 빨간 선: 중앙값",
           x = "시간 (년)", y = ep) +
      theme_bw()
    ggplotly(p)
  })

  output$plt_vp_dist <- renderPlotly({
    d <- vp_data()
    if (is.null(d)) return(NULL)
    ep <- input$vp_endpoint
    last_tp <- d %>% filter(abs(time - max(time)) < 31) %>%
      group_by(ID) %>% slice(1) %>% ungroup()
    p <- ggplot(last_tp, aes_string(ep)) +
      geom_histogram(bins = 20, fill = "steelblue", color = "white", alpha = 0.8) +
      labs(title = paste("말기", ep, "분포"), x = ep, y = "환자 수") +
      theme_bw()
    ggplotly(p)
  })

  output$plt_vpc <- renderPlotly({
    d <- vp_data()
    if (is.null(d)) return(NULL)
    ep <- input$vp_endpoint
    vpc_df <- d %>%
      group_by(years = round(years, 2)) %>%
      summarise(
        p05 = quantile(.data[[ep]], 0.05, na.rm = TRUE),
        p25 = quantile(.data[[ep]], 0.25, na.rm = TRUE),
        p50 = quantile(.data[[ep]], 0.50, na.rm = TRUE),
        p75 = quantile(.data[[ep]], 0.75, na.rm = TRUE),
        p95 = quantile(.data[[ep]], 0.95, na.rm = TRUE),
        .groups = "drop"
      )
    p <- ggplot(vpc_df, aes(years)) +
      geom_ribbon(aes(ymin = p05, ymax = p95), fill = "lightblue", alpha = 0.4) +
      geom_ribbon(aes(ymin = p25, ymax = p75), fill = "steelblue", alpha = 0.4) +
      geom_line(aes(y = p50), color = "darkblue", linewidth = 1.2) +
      labs(title = paste("VPC (90% PI 및 50% PI) —", ep),
           x = "시간 (년)", y = ep) +
      theme_bw()
    ggplotly(p)
  })

  output$plt_nac_dist <- renderPlotly({
    d <- vp_data()
    if (is.null(d)) return(NULL)
    last_tp <- d %>% filter(abs(time - max(time)) < 31) %>%
      group_by(ID) %>% slice(1) %>% ungroup() %>%
      mutate(NAC = factor(round(NAC_stage), levels = c(1,2,3),
                          labels = c("1기", "2기", "3기")))
    p <- ggplot(last_tp, aes(NAC, fill = NAC)) +
      geom_bar() +
      scale_fill_manual(values = c("1기" = "#27AE60", "2기" = "#F39C12", "3기" = "#E74C3C")) +
      labs(title = "말기 NAC 병기 분포", x = "병기", y = "환자 수") +
      theme_bw() + theme(legend.position = "none")
    ggplotly(p)
  })
}

# ============================================================
# 앱 실행
# ============================================================
shinyApp(ui, server)
