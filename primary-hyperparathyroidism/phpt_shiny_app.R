## ============================================================
## Primary Hyperparathyroidism (PHPT) QSP — Shiny Dashboard
## ============================================================
## Tabs:
##  1. 환자 프로파일 & 개요 (Patient Profile & Disease Overview)
##  2. 약물 PK (Drug Pharmacokinetics)
##  3. Ca/PTH/VitD 동태 (Calcium / PTH / Vitamin D Dynamics)
##  4. 골 리모델링 (Bone Remodeling & BMD)
##  5. 시나리오 비교 (Scenario Comparison)
##  6. 바이오마커 대시보드 (Biomarker Dashboard)
## ============================================================

library(shiny)
library(mrgsolve)
library(dplyr)
library(ggplot2)
library(tidyr)
library(DT)
library(plotly)
library(shinydashboard)
library(patchwork)

## ---- Embed mrgsolve model ----
phpt_code <- '
$PARAM
PHPT_severity = 1
PTx_time      = 9999
PTx_success   = 1
kel_PTH   = 10.4
PTH_ss0   = 5.5
nH_CaSR   = 3.5
Ca_set    = 1.15
Ca0       = 1.20
kel_Ca    = 0.8
PTH_Ca_bone = 0.025
PTH_Ca_TmCa = 0.003
VitD125_Ca_int = 0.0015
Ca_loss_base = 0.15
Ca_GFR_factor = 0.003
PO4_0     = 1.0
kel_PO4   = 1.2
PTH_PO4   = -0.08
VitD25_0  = 75
VitD125_0 = 90
k_CYP27B1 = 0.05
k_CYP24A1 = 0.008
k_FGF23_VD = 0.003
kel_VD25  = 0.007
kel_VD125 = 0.15
GFR_VD_factor = 0.015
OB0       = 1.0
OC0       = 1.0
kOB_form  = 0.05
kOC_form  = 0.08
kOB_death = 0.05
kOC_death = 0.08
PTH_OC    = 0.008
PTH_OB    = 0.003
kRANKL    = 0.05
kRANKL_deg= 0.05
RANKL_s0  = 1.0
BMD_LS0   = 0.960
BMD_FN0   = 0.730
kBMD_form = 0.0003
kBMD_resorb = 0.0004
GFR0      = 85
kGFR_loss = 0.00005
Ca_urine0 = 180
PTH_UCa   = 1.8
Ca_UCa    = 95
ka_cin    = 6.0
F_cin     = 0.22
Vc_cin    = 55
Vp_cin    = 180
CL_cin    = 125
Q_cin     = 167
EC50_cin  = 3.0
Emax_cin  = 0.35
Vd_deno   = 3.1
kel_deno  = 0.025
kon_deno  = 0.001
koff_deno = 0.003
kdeg_complex = 0.05
EC50_deno = 1000
Emax_deno = 0.85
kRANKL_deg_base = 0.05

$CMT
A_cin C_cin1 C_cin2
C_deno C_rankl C_complex
PTH Ca PO4 VitD25 VitD125
OB OC RANKL_s
BMD_LS BMD_FN
Ca_urine GFR

$INIT
A_cin=0 C_cin1=0 C_cin2=0
C_deno=0 C_rankl=1 C_complex=0
PTH=5.5 Ca=1.20 PO4=1.0
VitD25=75 VitD125=90
OB=1.0 OC=1.0 RANKL_s=1.0
BMD_LS=0.960 BMD_FN=0.730
Ca_urine=180 GFR=85

$MAIN
double PTH_secretion_rate;
if (SOLVERTIME < PTx_time) {
  if (PHPT_severity == 0)      PTH_secretion_rate = PTH_ss0;
  else if (PHPT_severity == 1) PTH_secretion_rate = PTH_ss0 * 5.0;
  else                          PTH_secretion_rate = PTH_ss0 * 22.0;
} else {
  PTH_secretion_rate = (PTx_success == 1) ? PTH_ss0 * 1.1 : PTH_ss0 * 4.0;
}
double Cin_conc = C_cin1 / Vc_cin;
double CaSR_shift = (Emax_cin * Cin_conc) / (EC50_cin + Cin_conc);
double Ca_set_eff = Ca_set - CaSR_shift;
double CaSR_effect = pow(Ca / Ca_set_eff, nH_CaSR) /
                     (1.0 + pow(Ca / Ca_set_eff, nH_CaSR));
double Deno_conc = C_deno / Vd_deno;
double OC_inh_deno = (Emax_deno * Deno_conc) / (EC50_deno + Deno_conc);
double VitD_Ca_absorb = VitD125_Ca_int * VitD125;
double GFR_VD_effect = GFR_VD_factor * (GFR / GFR0);
double PTH_bone_Ca = PTH_Ca_bone * PTH * OC / OB;
double PTH_TmCa = PTH_Ca_TmCa * PTH;
double Ca_in = VitD_Ca_absorb + PTH_bone_Ca + PTH_TmCa;
double Ca_out = Ca_loss_base + Ca_GFR_factor * GFR * Ca;
double RANKL_net = RANKL_s * (1.0 + PTH_OC * PTH) * (1.0 - OC_inh_deno);

$ODE
dxdt_A_cin = -ka_cin * A_cin;
dxdt_C_cin1 = ka_cin * F_cin * A_cin
              - (CL_cin/Vc_cin + Q_cin/Vc_cin) * C_cin1
              + (Q_cin/Vp_cin) * C_cin2;
dxdt_C_cin2 = (Q_cin/Vc_cin)*C_cin1 - (Q_cin/Vp_cin)*C_cin2;
dxdt_C_deno = -kel_deno*C_deno - kon_deno*C_deno*C_rankl + koff_deno*C_complex;
dxdt_C_rankl = kRANKL - kRANKL_deg_base*C_rankl
               - kon_deno*C_deno*C_rankl + koff_deno*C_complex;
dxdt_C_complex = kon_deno*C_deno*C_rankl - koff_deno*C_complex - kdeg_complex*C_complex;
double PTH_release = PTH_secretion_rate * (1.0 - CaSR_effect) * (1.0 - 0.7*CaSR_shift);
dxdt_PTH = PTH_release - kel_PTH * PTH;
dxdt_Ca  = Ca_in - Ca_out;
dxdt_PO4 = (0.3 + 0.1*VitD125/90.0) - kel_PO4*PO4 - fabs(PTH_PO4)*PTH;
double CYP27B1_activity = k_CYP27B1*PTH + GFR_VD_effect;
dxdt_VitD25  = 0.52 - kel_VD25*VitD25 - CYP27B1_activity*VitD25*0.1;
dxdt_VitD125 = CYP27B1_activity*VitD25*0.1 - kel_VD125*VitD125 - k_CYP24A1*VitD125*VitD125;
dxdt_OB = kOB_form*(1+PTH_OB*PTH*0.5) - kOB_death*OB;
dxdt_OC = kOC_form*RANKL_net - kOC_death*OC;
dxdt_RANKL_s = PTH_OC*PTH*kRANKL - kRANKL_deg_base*RANKL_s;
dxdt_BMD_LS = kBMD_form*OB - kBMD_resorb*OC;
dxdt_BMD_FN = (kBMD_form*OB - kBMD_resorb*OC)*0.85;
dxdt_Ca_urine = 0.5*(Ca_urine0 + 1.8*(PTH-5.5) + 95*(Ca-1.2) - Ca_urine);
dxdt_GFR = -kGFR_loss*(Ca_urine/Ca_urine0-1.0)*GFR;

$TABLE
double Cin_ng_mL  = C_cin1 / Vc_cin;
double Deno_ng_mL = C_deno / Vd_deno;
double Ca_total   = Ca * 2.0;
double PTH_pg_mL  = PTH * 9.425;
double OC_OB_ratio = OC / OB;
double CTX_index  = OC * 280;
double P1NP_index = OB * 55;
double T_score_LS = (BMD_LS - 0.96) / 0.12;
double T_score_FN = (BMD_FN - 0.73) / 0.10;

$CAPTURE
PTH PTH_pg_mL Ca PO4 VitD125 VitD25
OB OC RANKL_s BMD_LS BMD_FN Ca_urine GFR
Cin_ng_mL Deno_ng_mL OC_OB_ratio CTX_index P1NP_index
T_score_LS T_score_FN Ca_total
'

## Compile
mod_base <- mread_cache("phpt_shiny", tempdir(), phpt_code)

## ---- Helper: run simulation ----
run_sim <- function(severity, ptx_time = 9999, cin_dose = 0,
                    deno_dose = 0, deno_interval = 180,
                    duration = 1825, gfr0 = 85) {
  params <- list(PHPT_severity = severity, PTx_time = ptx_time,
                 GFR0 = gfr0)

  events <- NULL
  if (cin_dose > 0) {
    events <- bind_rows(events,
      tibble(time = seq(0, duration, by = 1), amt = cin_dose, cmt = 2,
             evid = 1, rate = 0))
  }
  if (deno_dose > 0) {
    deno_times <- seq(0, duration, by = deno_interval)
    events <- bind_rows(events,
      tibble(time = deno_times, amt = deno_dose * 1000, cmt = 4,
             evid = 1, rate = 0.069 * deno_dose * 1000))
  }

  m <- mod_base %>% param(params)
  if (!is.null(events)) {
    out <- m %>% data_set(events) %>%
      mrgsim(end = duration, delta = 1, obsonly = TRUE)
  } else {
    out <- m %>% mrgsim(end = duration, delta = 1, obsonly = TRUE)
  }
  as_tibble(out)
}

## ---- UI ----

ui <- dashboardPage(
  skin = "blue",

  dashboardHeader(title = "PHPT QSP Dashboard", titleWidth = 280),

  dashboardSidebar(
    width = 260,
    sidebarMenu(
      menuItem("환자 프로파일", tabName = "tab1", icon = icon("user-md")),
      menuItem("약물 PK", tabName = "tab2", icon = icon("pills")),
      menuItem("Ca/PTH/VitD 동태", tabName = "tab3", icon = icon("chart-line")),
      menuItem("골 리모델링 & BMD", tabName = "tab4", icon = icon("bone")),
      menuItem("시나리오 비교", tabName = "tab5", icon = icon("balance-scale")),
      menuItem("바이오마커 대시보드", tabName = "tab6", icon = icon("vials"))
    ),
    hr(),
    h5("질환 설정", style = "padding-left:15px; color:#ccc"),
    selectInput("severity", "PHPT 중증도",
      choices = c("정상 (Normal)" = 0,
                  "경증 PHPT (Mild)" = 1,
                  "중증 PHPT (Severe)" = 2),
      selected = 1),
    numericInput("gfr0", "기저 eGFR (mL/min)", value = 85, min = 10, max = 120, step = 5),
    numericInput("duration_yr", "시뮬레이션 기간 (년)", value = 5, min = 1, max = 10, step = 1),
    hr(),
    h5("치료 설정", style = "padding-left:15px; color:#ccc"),
    numericInput("ptx_day", "부갑상선절제술 시점 (일; 0=없음)", value = 0, min = 0, max = 3650),
    numericInput("cin_dose", "Cinacalcet 용량 (mg/day; 0=없음)", value = 0, min = 0, max = 180, step = 30),
    numericInput("deno_dose", "Denosumab 용량 (mg; 0=없음)", value = 0, min = 0, max = 120, step = 60),
    numericInput("deno_int", "Denosumab 투여 간격 (일)", value = 180, min = 90, max = 360, step = 90),
    actionButton("run_btn", "시뮬레이션 실행", class = "btn-primary", width = "100%")
  ),

  dashboardBody(
    tags$head(tags$style(HTML("
      .content-wrapper { background-color: #F8F9FA; }
      .box { border-radius: 6px; }
      .small-box { border-radius: 6px; }
    "))),

    tabItems(

      ## ---- Tab 1: Patient Profile ----
      tabItem(tabName = "tab1",
        fluidRow(
          box(title = "원발성 부갑상선 기능 항진증 (PHPT) 개요", width = 12,
              status = "primary", solidHeader = TRUE,
              tabsetPanel(
                tabPanel("질환 개요",
                  h4("병태생리"),
                  p("PHPT는 부갑상선 선종(85%), 다발성 선종(3%), 증식(12%), 암종(1%)에 의한
                    PTH 자율분비로 발생합니다. PTH과잉은 골흡수↑, 신장 Ca 재흡수↑,
                    1,25(OH)₂D 합성↑를 통해 고칼슘혈증을 유발합니다."),
                  h4("임상 분류"),
                  tableOutput("classification_tbl"),
                  h4("NIH 2022 수술 적응증"),
                  tableOutput("nih_criteria_tbl")
                ),
                tabPanel("CaSR-PTH 피드백",
                  plotlyOutput("casr_curve", height = "400px"),
                  p("정상: 혈청 Ca 상승 시 CaSR 활성화 → PTH 분비 억제.
                    PHPT: 자율분비로 인해 억제 불가.")
                ),
                tabPanel("약물 요약",
                  h4("주요 약물 PK/PD"),
                  tableOutput("drug_summary_tbl")
                )
              )
          )
        ),
        fluidRow(
          valueBoxOutput("vbox_nih_ca", width = 3),
          valueBoxOutput("vbox_nih_epi", width = 3),
          valueBoxOutput("vbox_nih_bone", width = 3),
          valueBoxOutput("vbox_nih_age", width = 3)
        )
      ),

      ## ---- Tab 2: Drug PK ----
      tabItem(tabName = "tab2",
        fluidRow(
          box(title = "Cinacalcet PK — 혈중 농도 추이", width = 6,
              status = "info", solidHeader = TRUE,
              plotlyOutput("cin_pk_plot", height = "320px")),
          box(title = "Denosumab PK — 혈중 농도 추이 (TMDD)", width = 6,
              status = "success", solidHeader = TRUE,
              plotlyOutput("deno_pk_plot", height = "320px"))
        ),
        fluidRow(
          box(title = "Cinacalcet Dose-Response (EC50 모델)", width = 6,
              status = "info", solidHeader = TRUE,
              plotlyOutput("cin_ec50_plot", height = "300px")),
          box(title = "약물 PK 파라미터 요약", width = 6,
              status = "primary", solidHeader = TRUE,
              DTOutput("pk_params_tbl"))
        )
      ),

      ## ---- Tab 3: Ca/PTH/VitD ----
      tabItem(tabName = "tab3",
        fluidRow(
          box(title = "iPTH 혈중 농도", width = 6,
              status = "danger", solidHeader = TRUE,
              plotlyOutput("pth_plot", height = "300px")),
          box(title = "혈청 칼슘 (Total)", width = 6,
              status = "warning", solidHeader = TRUE,
              plotlyOutput("ca_plot", height = "300px"))
        ),
        fluidRow(
          box(title = "혈청 인산염 (Phosphate)", width = 6,
              status = "primary", solidHeader = TRUE,
              plotlyOutput("po4_plot", height = "300px")),
          box(title = "1,25(OH)₂D (Calcitriol) & 25(OH)D", width = 6,
              status = "info", solidHeader = TRUE,
              plotlyOutput("vitd_plot", height = "300px"))
        ),
        fluidRow(
          box(title = "요중 칼슘 배설 (Urinary Ca)", width = 6,
              status = "warning", solidHeader = TRUE,
              plotlyOutput("ucal_plot", height = "300px")),
          box(title = "eGFR 추이", width = 6,
              status = "primary", solidHeader = TRUE,
              plotlyOutput("gfr_plot", height = "300px"))
        )
      ),

      ## ---- Tab 4: Bone Remodeling ----
      tabItem(tabName = "tab4",
        fluidRow(
          box(title = "요추 BMD (Lumbar Spine)", width = 6,
              status = "purple", solidHeader = TRUE,
              plotlyOutput("bmd_ls_plot", height = "300px")),
          box(title = "대퇴경부 BMD (Femoral Neck)", width = 6,
              status = "purple", solidHeader = TRUE,
              plotlyOutput("bmd_fn_plot", height = "300px"))
        ),
        fluidRow(
          box(title = "OB/OC 비율 & 골 리모델링 마커", width = 6,
              status = "danger", solidHeader = TRUE,
              plotlyOutput("oc_ob_plot", height = "300px")),
          box(title = "T-score 추이 (L1-L4)", width = 6,
              status = "warning", solidHeader = TRUE,
              plotlyOutput("tscore_plot", height = "300px"))
        ),
        fluidRow(
          box(title = "혈청 CTX & P1NP (골 회전 마커)", width = 12,
              status = "info", solidHeader = TRUE,
              plotlyOutput("btm_plot", height = "280px"))
        )
      ),

      ## ---- Tab 5: Scenario Comparison ----
      tabItem(tabName = "tab5",
        fluidRow(
          box(title = "시나리오 선택 (최대 4개)", width = 4,
              status = "primary", solidHeader = TRUE,
              checkboxGroupInput("sc_select",
                label = "비교할 시나리오 선택",
                choices = c(
                  "정상 (Normal)" = "sc0",
                  "미치료 PHPT (경증)" = "sc1",
                  "미치료 PHPT (중증)" = "sc2",
                  "Cinacalcet 60mg" = "sc3",
                  "Denosumab q6mo" = "sc4",
                  "부갑상선절제술 (day 90)" = "sc5",
                  "Cin + Deno 병용" = "sc6"
                ),
                selected = c("sc1", "sc3", "sc5")
              ),
              actionButton("run_compare", "비교 실행", class = "btn-primary")
          ),
          box(title = "Year-1 요약 테이블", width = 8,
              status = "info", solidHeader = TRUE,
              DTOutput("compare_tbl"))
        ),
        fluidRow(
          box(title = "PTH 비교", width = 6, status = "danger", solidHeader = TRUE,
              plotlyOutput("cmp_pth", height = "280px")),
          box(title = "혈청 Ca 비교", width = 6, status = "warning", solidHeader = TRUE,
              plotlyOutput("cmp_ca", height = "280px"))
        ),
        fluidRow(
          box(title = "요추 BMD 비교", width = 6, status = "success", solidHeader = TRUE,
              plotlyOutput("cmp_bmd", height = "280px")),
          box(title = "요중 Ca 비교", width = 6, status = "info", solidHeader = TRUE,
              plotlyOutput("cmp_ucal", height = "280px"))
        )
      ),

      ## ---- Tab 6: Biomarker Dashboard ----
      tabItem(tabName = "tab6",
        fluidRow(
          valueBoxOutput("vb_pth_now", width = 3),
          valueBoxOutput("vb_ca_now", width = 3),
          valueBoxOutput("vb_bmd_now", width = 3),
          valueBoxOutput("vb_gfr_now", width = 3)
        ),
        fluidRow(
          box(title = "시뮬레이션 결과 — 전체 바이오마커 테이블", width = 12,
              status = "primary", solidHeader = TRUE,
              DTOutput("biomarker_tbl"))
        ),
        fluidRow(
          box(title = "바이오마커 시계열 (Heatmap)", width = 12,
              status = "info", solidHeader = TRUE,
              plotlyOutput("heatmap_plot", height = "400px"))
        )
      )
    )
  )
)

## ---- SERVER ----

server <- function(input, output, session) {

  ## Reactive: run simulation
  sim_data <- eventReactive(input$run_btn, {
    severity <- as.integer(input$severity)
    ptx_day  <- if (input$ptx_day > 0) input$ptx_day else 9999
    dur_days <- input$duration_yr * 365
    run_sim(severity = severity, ptx_time = ptx_day,
            cin_dose = input$cin_dose, deno_dose = input$deno_dose,
            deno_interval = input$deno_int, duration = dur_days,
            gfr0 = input$gfr0)
  }, ignoreNULL = FALSE)

  ## ---- Tab 1 outputs ----

  output$classification_tbl <- renderTable({
    data.frame(
      분류 = c("무증상 PHPT", "정상칼슘 PHPT", "경증 PHPT", "중증 PHPT", "고칼슘혈증 위기"),
      `혈청 Ca` = c("정상 상한 이내", "정상", "ULN + <1 mg/dL", "> ULN + 1 mg/dL", "> 14 mg/dL"),
      `혈청 PTH` = c("↑", "↑↑", "↑", "↑↑↑", "↑↑↑"),
      빈도 = c("75-80%", "15%", "10%", "5%", "<1%"),
      stringsAsFactors = FALSE
    )
  })

  output$nih_criteria_tbl <- renderTable({
    data.frame(
      `수술 적응증 (NIH 2022)` = c(
        "혈청 Ca > 정상 상한 + 1 mg/dL (> 0.25 mmol/L)",
        "eGFR < 60 mL/min/1.73m²",
        "요로결석 또는 신석회화증",
        "요중 Ca > 400 mg/day + 고위험 결석 프로파일",
        "DXA: T-score < -2.5 (요추/대퇴/요골)",
        "척추 골절 (방사선/VFA)",
        "나이 < 50세",
        "원하지 않는 경우 (비수술 불가)"
      ),
      stringsAsFactors = FALSE
    )
  })

  output$drug_summary_tbl <- renderTable({
    data.frame(
      약물 = c("Cinacalcet (Sensipar)", "Denosumab (Prolia)", "Alendronate (Fosamax)"),
      기전 = c("CaSR allosteric potentiator → PTH↓",
               "Anti-RANKL mAb → OC 억제 → BMD↑",
               "Farnesyl-PP synthase 억제 → OC 사멸"),
      용량 = c("30-180 mg/day po", "60 mg SC q6mo", "70 mg/wk po"),
      주요효과 = c("Ca↓ ~0.5 mg/dL, PTH↓ 30-50%",
                  "CTX↓ 70-80%, BMD↑ 5%/yr",
                  "CTX↓ 50-70%, BMD 안정"),
      stringsAsFactors = FALSE
    )
  })

  ## CaSR Curve
  output$casr_curve <- renderPlotly({
    ca_seq <- seq(0.8, 1.8, by = 0.01)
    pth_normal <- sapply(ca_seq, function(x) {
      1 / (1 + (x / 1.15)^3.5) * 100
    })
    pth_cin <- sapply(ca_seq, function(x) {
      1 / (1 + (x / 0.90)^3.5) * 100
    })
    pth_phpt <- rep(100, length(ca_seq))
    df_casr <- data.frame(Ca = ca_seq, Normal = pth_normal,
                          Cinacalcet = pth_cin, PHPT = pth_phpt)
    plot_ly(df_casr, x = ~Ca) %>%
      add_lines(y = ~Normal, name = "Normal", line = list(color = "steelblue")) %>%
      add_lines(y = ~Cinacalcet, name = "Cinacalcet", line = list(color = "green", dash = "dash")) %>%
      add_lines(y = ~PHPT, name = "PHPT (autonomous)", line = list(color = "red")) %>%
      layout(title = "CaSR-PTH Sigmoidal Set-Point Curve",
             xaxis = list(title = "Ionized Ca2+ (mM)"),
             yaxis = list(title = "Relative PTH Secretion (%)"),
             legend = list(orientation = "h"))
  })

  output$vbox_nih_ca <- renderValueBox({
    valueBox("> 1 mg/dL ULN", "Ca 수술 기준", icon = icon("tint"), color = "red")
  })
  output$vbox_nih_epi <- renderValueBox({
    valueBox("< 60 mL/min", "eGFR 기준", icon = icon("kidney"), color = "yellow")
  })
  output$vbox_nih_bone <- renderValueBox({
    valueBox("T < -2.5", "BMD 기준", icon = icon("bone"), color = "orange")
  })
  output$vbox_nih_age <- renderValueBox({
    valueBox("< 50세", "연령 기준", icon = icon("user"), color = "blue")
  })

  ## ---- Tab 2: Drug PK ----

  output$cin_pk_plot <- renderPlotly({
    d <- sim_data()
    plot_ly(d, x = ~time, y = ~Cin_ng_mL, type = "scatter", mode = "lines",
            line = list(color = "#2196F3")) %>%
      add_lines(y = rep(3, nrow(d)), name = "EC50 = 3 ng/mL",
                line = list(dash = "dash", color = "red")) %>%
      layout(title = "Cinacalcet 혈중 농도",
             xaxis = list(title = "시간 (일)"),
             yaxis = list(title = "농도 (ng/mL)"))
  })

  output$deno_pk_plot <- renderPlotly({
    d <- sim_data()
    plot_ly(d, x = ~time, y = ~Deno_ng_mL, type = "scatter", mode = "lines",
            line = list(color = "#4CAF50")) %>%
      layout(title = "Denosumab 혈중 농도 (TMDD)",
             xaxis = list(title = "시간 (일)"),
             yaxis = list(title = "농도 (ng/mL)"))
  })

  output$cin_ec50_plot <- renderPlotly({
    conc <- seq(0, 50, by = 0.5)
    eff  <- 0.35 * conc / (3.0 + conc)
    df_ec <- data.frame(conc = conc, eff = eff)
    plot_ly(df_ec, x = ~conc, y = ~eff, type = "scatter", mode = "lines",
            line = list(color = "#2196F3")) %>%
      add_segments(x = 3, xend = 3, y = 0, yend = 0.175,
                   line = list(dash = "dash", color = "red")) %>%
      add_segments(x = 0, xend = 3, y = 0.175, yend = 0.175,
                   line = list(dash = "dash", color = "red")) %>%
      layout(title = "Cinacalcet Emax Model (CaSR Set-point Shift)",
             xaxis = list(title = "Cinacalcet 농도 (ng/mL)"),
             yaxis = list(title = "Set-point Shift (mM)"))
  })

  output$pk_params_tbl <- renderDT({
    df <- data.frame(
      약물 = c("Cinacalcet", "Cinacalcet", "Cinacalcet", "Cinacalcet",
               "Denosumab", "Denosumab", "Denosumab",
               "Alendronate", "Alendronate"),
      파라미터 = c("F (경구 생체이용률)", "Vc (L)", "CL (L/day)", "t½",
                 "F (SC)", "Vd (L)", "t½",
                 "F (경구)", "골 결합 t½"),
      값 = c("22% (공복)", "55", "125", "6-8h",
             "62%", "3.1", "~28일",
             "0.6%", "~10년"),
      stringsAsFactors = FALSE
    )
    datatable(df, options = list(pageLength = 10, dom = "t"),
              rownames = FALSE)
  })

  ## ---- Tab 3: Ca/PTH/VitD ----

  make_line_plot <- function(data, y_var, title, ylab, hlines = NULL, colors = "#F44336") {
    d <- data
    p <- plot_ly(d, x = ~time, y = ~get(y_var), type = "scatter", mode = "lines",
                 line = list(color = colors), name = y_var) %>%
      layout(title = title,
             xaxis = list(title = "시간 (일)"),
             yaxis = list(title = ylab))
    if (!is.null(hlines)) {
      for (hl in hlines) {
        p <- p %>% add_lines(y = rep(hl$y, nrow(d)), name = hl$label,
                             line = list(dash = "dash", color = hl$color))
      }
    }
    p
  }

  output$pth_plot <- renderPlotly({
    make_line_plot(sim_data(), "PTH_pg_mL", "iPTH 혈중 농도", "PTH (pg/mL)",
      hlines = list(list(y = 15, label = "LLN 15", color = "gray"),
                    list(y = 65, label = "ULN 65", color = "orange")))
  })

  output$ca_plot <- renderPlotly({
    make_line_plot(sim_data(), "Ca_total", "총 혈청 칼슘", "Ca (mM)",
      hlines = list(list(y = 2.2, label = "정상 하한", color = "blue"),
                    list(y = 2.6, label = "정상 상한", color = "red")),
      colors = "#FF9800")
  })

  output$po4_plot <- renderPlotly({
    make_line_plot(sim_data(), "PO4", "혈청 인산염", "PO4 (mM)",
      hlines = list(list(y = 0.8, label = "정상 하한", color = "blue"),
                    list(y = 1.5, label = "정상 상한", color = "red")),
      colors = "#2196F3")
  })

  output$vitd_plot <- renderPlotly({
    d <- sim_data()
    plot_ly(d, x = ~time) %>%
      add_lines(y = ~VitD125, name = "1,25(OH)₂D (pmol/L)", line = list(color = "#FF9800")) %>%
      add_lines(y = ~(VitD25 * 2), name = "25(OH)D ×2 (nmol/L)", line = list(color = "#4CAF50")) %>%
      layout(title = "Vitamin D 대사물", xaxis = list(title = "시간 (일)"),
             yaxis = list(title = "농도"), legend = list(orientation = "h"))
  })

  output$ucal_plot <- renderPlotly({
    make_line_plot(sim_data(), "Ca_urine", "요중 칼슘 배설량", "요중 Ca (mg/day)",
      hlines = list(list(y = 300, label = "고칼슘뇨증 기준", color = "red")),
      colors = "#9C27B0")
  })

  output$gfr_plot <- renderPlotly({
    make_line_plot(sim_data(), "GFR", "eGFR 추이", "eGFR (mL/min/1.73m²)",
      hlines = list(list(y = 60, label = "CKD G2/G3 경계", color = "orange"),
                    list(y = 30, label = "CKD G4 경계", color = "red")),
      colors = "#00838F")
  })

  ## ---- Tab 4: Bone ----

  output$bmd_ls_plot <- renderPlotly({
    make_line_plot(sim_data(), "BMD_LS", "요추 BMD (L1-L4)", "BMD (g/cm²)",
      hlines = list(list(y = 0.84, label = "골감소증 기준", color = "orange"),
                    list(y = 0.72, label = "골다공증 기준", color = "red")),
      colors = "#7B1FA2")
  })

  output$bmd_fn_plot <- renderPlotly({
    make_line_plot(sim_data(), "BMD_FN", "대퇴경부 BMD", "BMD (g/cm²)",
      hlines = list(list(y = 0.63, label = "골감소증 기준", color = "orange"),
                    list(y = 0.53, label = "골다공증 기준", color = "red")),
      colors = "#6A1B9A")
  })

  output$oc_ob_plot <- renderPlotly({
    d <- sim_data()
    plot_ly(d, x = ~time) %>%
      add_lines(y = ~OC, name = "OC (파골세포)", line = list(color = "red")) %>%
      add_lines(y = ~OB, name = "OB (조골세포)", line = list(color = "blue")) %>%
      add_lines(y = ~OC_OB_ratio, name = "OC/OB 비율", line = list(color = "purple", dash = "dash")) %>%
      layout(title = "OB/OC 세포 역학", xaxis = list(title = "시간 (일)"),
             yaxis = list(title = "상대 세포수"), legend = list(orientation = "h"))
  })

  output$tscore_plot <- renderPlotly({
    d <- sim_data()
    plot_ly(d, x = ~time) %>%
      add_lines(y = ~T_score_LS, name = "요추 T-score", line = list(color = "#7B1FA2")) %>%
      add_lines(y = ~T_score_FN, name = "대퇴경부 T-score", line = list(color = "#AB47BC")) %>%
      add_lines(y = rep(-1, nrow(d)), name = "골감소증 기준 (-1)",
                line = list(dash = "dash", color = "orange")) %>%
      add_lines(y = rep(-2.5, nrow(d)), name = "골다공증 기준 (-2.5)",
                line = list(dash = "dash", color = "red")) %>%
      layout(title = "T-score 추이", xaxis = list(title = "시간 (일)"),
             yaxis = list(title = "T-score"), legend = list(orientation = "h"))
  })

  output$btm_plot <- renderPlotly({
    d <- sim_data()
    plot_ly(d, x = ~time) %>%
      add_lines(y = ~CTX_index, name = "CTX 지표 (ng/L)", line = list(color = "red")) %>%
      add_lines(y = ~P1NP_index, name = "P1NP 지표 (µg/L)", line = list(color = "green")) %>%
      layout(title = "골 회전 마커 추이",
             xaxis = list(title = "시간 (일)"),
             yaxis = list(title = "마커 수준"),
             legend = list(orientation = "h"))
  })

  ## ---- Tab 5: Scenario Comparison ----

  sc_configs <- list(
    sc0 = list(severity = 0, ptx = 9999, cin = 0, deno = 0),
    sc1 = list(severity = 1, ptx = 9999, cin = 0, deno = 0),
    sc2 = list(severity = 2, ptx = 9999, cin = 0, deno = 0),
    sc3 = list(severity = 1, ptx = 9999, cin = 60, deno = 0),
    sc4 = list(severity = 1, ptx = 9999, cin = 0, deno = 60),
    sc5 = list(severity = 1, ptx = 90, cin = 0, deno = 0),
    sc6 = list(severity = 1, ptx = 9999, cin = 60, deno = 60)
  )
  sc_labels <- c(
    sc0 = "정상", sc1 = "미치료(경)", sc2 = "미치료(중)",
    sc3 = "Cinacalcet 60", sc4 = "Denosumab q6mo",
    sc5 = "PTx (day90)", sc6 = "Cin+Deno"
  )

  compare_data <- eventReactive(input$run_compare, {
    selected <- input$sc_select
    bind_rows(lapply(selected, function(sc) {
      cfg <- sc_configs[[sc]]
      run_sim(severity = cfg$severity, ptx_time = cfg$ptx, cin_dose = cfg$cin,
              deno_dose = cfg$deno, duration = 1825, gfr0 = input$gfr0) %>%
        mutate(scenario = sc_labels[[sc]])
    }))
  }, ignoreNULL = FALSE)

  output$compare_tbl <- renderDT({
    d <- compare_data()
    tbl <- d %>%
      filter(time == 365) %>%
      group_by(scenario) %>%
      summarise(
        `PTH (pg/mL)` = round(mean(PTH_pg_mL), 0),
        `Ca Total (mM)` = round(mean(Ca_total), 2),
        `BMD LS (g/cm²)` = round(mean(BMD_LS), 3),
        `요중 Ca (mg/d)` = round(mean(Ca_urine), 0),
        `eGFR` = round(mean(GFR), 1),
        `T-score LS` = round(mean(T_score_LS), 2),
        .groups = "drop"
      )
    datatable(tbl, rownames = FALSE,
              options = list(pageLength = 10, dom = "t"))
  })

  output$cmp_pth <- renderPlotly({
    d <- compare_data()
    plot_ly(d, x = ~time, y = ~PTH_pg_mL, color = ~scenario,
            type = "scatter", mode = "lines") %>%
      layout(title = "PTH 비교", xaxis = list(title = "일"),
             yaxis = list(title = "PTH (pg/mL)"), legend = list(orientation = "h"))
  })

  output$cmp_ca <- renderPlotly({
    d <- compare_data()
    plot_ly(d, x = ~time, y = ~Ca_total, color = ~scenario,
            type = "scatter", mode = "lines") %>%
      layout(title = "혈청 Ca 비교", xaxis = list(title = "일"),
             yaxis = list(title = "Ca (mM)"), legend = list(orientation = "h"))
  })

  output$cmp_bmd <- renderPlotly({
    d <- compare_data()
    plot_ly(d, x = ~time, y = ~BMD_LS, color = ~scenario,
            type = "scatter", mode = "lines") %>%
      layout(title = "요추 BMD 비교", xaxis = list(title = "일"),
             yaxis = list(title = "BMD (g/cm²)"), legend = list(orientation = "h"))
  })

  output$cmp_ucal <- renderPlotly({
    d <- compare_data()
    plot_ly(d, x = ~time, y = ~Ca_urine, color = ~scenario,
            type = "scatter", mode = "lines") %>%
      layout(title = "요중 Ca 비교", xaxis = list(title = "일"),
             yaxis = list(title = "mg/day"), legend = list(orientation = "h"))
  })

  ## ---- Tab 6: Biomarker Dashboard ----

  output$vb_pth_now <- renderValueBox({
    d <- tail(sim_data(), 1)
    val <- round(d$PTH_pg_mL, 0)
    color <- if (val > 65) "red" else "green"
    valueBox(paste0(val, " pg/mL"), "현재 iPTH", icon = icon("syringe"), color = color)
  })

  output$vb_ca_now <- renderValueBox({
    d <- tail(sim_data(), 1)
    val <- round(d$Ca_total, 2)
    color <- if (val > 2.6) "red" else "green"
    valueBox(paste0(val, " mM"), "혈청 총 Ca", icon = icon("tint"), color = color)
  })

  output$vb_bmd_now <- renderValueBox({
    d <- tail(sim_data(), 1)
    val <- round(d$BMD_LS, 3)
    ts  <- round(d$T_score_LS, 1)
    color <- if (ts < -2.5) "red" else if (ts < -1) "yellow" else "green"
    valueBox(paste0(val, " (T=", ts, ")"), "요추 BMD", icon = icon("bone"), color = color)
  })

  output$vb_gfr_now <- renderValueBox({
    d <- tail(sim_data(), 1)
    val <- round(d$GFR, 0)
    color <- if (val < 30) "red" else if (val < 60) "yellow" else "green"
    valueBox(paste0(val, " mL/min"), "eGFR", icon = icon("kidney"), color = color)
  })

  output$biomarker_tbl <- renderDT({
    d <- sim_data() %>%
      filter(time %% 90 == 0) %>%
      mutate(across(where(is.numeric), ~round(.x, 2))) %>%
      select(time, PTH_pg_mL, Ca_total, PO4, VitD125, VitD25,
             BMD_LS, T_score_LS, Ca_urine, GFR, OC_OB_ratio,
             CTX_index, P1NP_index)
    datatable(d, rownames = FALSE,
              options = list(pageLength = 10, scrollX = TRUE))
  })

  output$heatmap_plot <- renderPlotly({
    d <- sim_data() %>%
      filter(time %% 90 == 0) %>%
      select(time, PTH_pg_mL, Ca_total, BMD_LS, Ca_urine, GFR, VitD125) %>%
      mutate(across(-time, ~scale(.x)[, 1]))
    d_long <- pivot_longer(d, -time, names_to = "Marker", values_to = "ZScore")
    plot_ly(d_long, x = ~time, y = ~Marker, z = ~ZScore, type = "heatmap",
            colorscale = "RdBu", reversescale = TRUE) %>%
      layout(title = "바이오마커 Z-Score 히트맵",
             xaxis = list(title = "시간 (일)"),
             yaxis = list(title = ""))
  })
}

## ---- Launch ----
shinyApp(ui, server)
