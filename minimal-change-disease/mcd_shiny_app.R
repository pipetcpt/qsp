# ==============================================================================
# Minimal Change Disease (MCD) QSP Shiny Application
# Quantitative Systems Pharmacology Model
# ==============================================================================
# Dependencies: shiny, shinydashboard, shinyWidgets, plotly, DT, dplyr, tidyr
# ==============================================================================

library(shiny)
library(shinydashboard)
library(shinyWidgets)
library(plotly)
library(DT)
library(dplyr)
library(tidyr)

# ==============================================================================
# Helper Functions
# ==============================================================================

# PK Model: One-compartment oral (Prednisolone / CNI)
pk_onecomp <- function(time, dose, CL, Vd, ka, F_oral = 1, tau = 24) {
  # Steady-state superposition for oral dosing
  n_doses <- floor(max(time) / tau) + 1
  conc <- numeric(length(time))
  for (n in 0:(n_doses - 1)) {
    t_adj <- time - n * tau
    t_adj[t_adj < 0] <- NA
    ke <- CL / Vd
    c_n <- (F_oral * dose * ka) / (Vd * (ka - ke)) *
      (exp(-ke * t_adj) - exp(-ka * t_adj))
    c_n[is.na(c_n)] <- 0
    c_n[c_n < 0] <- 0
    conc <- conc + c_n
  }
  conc
}

# PK Model: Two-compartment IV (Rituximab)
pk_twocomp_iv <- function(time, dose, CL, Vc, Q, Vp) {
  k10 <- CL / Vc
  k12 <- Q / Vc
  k21 <- Q / Vp
  alpha <- 0.5 * ((k10 + k12 + k21) + sqrt((k10 + k12 + k21)^2 - 4 * k10 * k21))
  beta  <- 0.5 * ((k10 + k12 + k21) - sqrt((k10 + k12 + k21)^2 - 4 * k10 * k21))
  A <- (dose / Vc) * (alpha - k21) / (alpha - beta)
  B <- (dose / Vc) * (k21 - beta)  / (alpha - beta)
  conc <- A * exp(-alpha * time) + B * exp(-beta * time)
  conc[conc < 0] <- 0
  conc
}

# Disease progression / PD model
run_pd_model <- function(
    time_days,
    drug = "Prednisolone",
    dose = 1,          # mg/kg/day
    baseline_uprot = 8,
    baseline_alb   = 2.5,
    age            = 10,
    disease_type   = "primary",
    relapse_hist   = "first_episode"
) {
  set.seed(42)
  n <- length(time_days)

  # ---- Baseline severity factor ----
  sev <- (baseline_uprot / 20) * 0.5 + (1 - baseline_alb / 4.5) * 0.5
  relapse_factor <- switch(relapse_hist,
    first_episode   = 1.0,
    frequent_relapser = 1.3,
    steroid_dependent = 1.5,
    steroid_resistant = 2.0,
    1.0
  )

  # ---- Drug effect parameters ----
  params <- switch(drug,
    Prednisolone = list(
      t_onset = 14,  t_peak = 28, max_eff = 0.85 / relapse_factor,
      CL = 15, Vd = 50, ka = 1.2, tau = 24,
      trough_therapeutic = c(NA, NA)
    ),
    Cyclosporine = list(
      t_onset = 7,   t_peak = 21, max_eff = 0.72 / relapse_factor,
      CL = 5,  Vd = 250, ka = 0.8, tau = 12,
      trough_therapeutic = c(100, 200)   # ng/mL
    ),
    Tacrolimus = list(
      t_onset = 7,   t_peak = 21, max_eff = 0.75 / relapse_factor,
      CL = 2.5, Vd = 480, ka = 0.6, tau = 12,
      trough_therapeutic = c(5, 10)      # ng/mL
    ),
    Rituximab = list(
      t_onset = 28,  t_peak = 90, max_eff = 0.80 / relapse_factor,
      CL = 0.015, Vc = 3.1, Q = 0.3, Vp = 2.0,
      trough_therapeutic = c(NA, NA)
    )
  )

  # ---- Immune compartments ----
  # CD4 Effector T cells (relative, 1 = baseline)
  cd4_eff <- pmax(0.1,
    1 - params$max_eff * pmin(1, time_days / params$t_peak) *
    exp(-pmax(0, time_days - params$t_peak) / 60) * sev)

  # Treg cells (counter-regulatory increase)
  treg <- pmin(1.5,
    1 + 0.4 * params$max_eff * pmin(1, time_days / (params$t_peak * 1.5)))

  treg_eff_ratio <- treg / pmax(0.1, cd4_eff)

  # B cells (rituximab depletes; CNIs partially reduce; steroids mild effect)
  bcell_depletion <- if (drug == "Rituximab") {
    pmax(0.05, 1 - 0.92 * pmin(1, time_days / 60) *
           exp(-pmax(0, time_days - 90) / 120))
  } else if (drug %in% c("Cyclosporine", "Tacrolimus")) {
    pmax(0.5, 1 - 0.35 * pmin(1, time_days / 30))
  } else {
    pmax(0.7, 1 - 0.25 * pmin(1, time_days / 21))
  }

  # suPAR (soluble urokinase plasminogen activator receptor)
  supar <- pmax(0.5,
    (1 + sev * 1.5) - params$max_eff * 0.9 *
    pmin(1, time_days / params$t_peak) *
    exp(-pmax(0, time_days - params$t_peak) / 80))

  # ---- Podocyte biology ----
  nephrin_expr <- pmin(1,
    0.3 + 0.7 * params$max_eff * pmin(1, time_days / (params$t_peak * 1.2)))

  anti_nephrin_ab <- pmax(0.05,
    if (drug == "Rituximab") {
      (1 + sev * 0.8) * pmax(0.05,
        1 - 0.95 * pmin(1, time_days / 90) *
          exp(-pmax(0, time_days - 90) / 150))
    } else {
      (1 + sev * 0.8) * pmax(0.1,
        1 - 0.6 * pmin(1, time_days / params$t_peak) *
          exp(-pmax(0, time_days - params$t_peak) / 60))
    })

  podocyte_integrity <- pmin(1,
    (0.4 - sev * 0.3) + (0.6 + sev * 0.3) *
    params$max_eff * pmin(1, time_days / params$t_peak) *
    exp(-pmax(0, time_days - params$t_peak) / 90))
  podocyte_integrity <- pmax(0, pmin(1, podocyte_integrity))

  foot_process_eff <- pmax(0,
    (0.9 + sev * 0.1) - params$max_eff * 0.95 *
    pmin(1, time_days / params$t_peak) *
    exp(-pmax(0, time_days - params$t_peak) / 60))

  # ---- Clinical endpoints ----
  uprot_cr <- baseline_uprot * pmax(0.02,
    1 - params$max_eff * pmin(1, time_days / params$t_peak) *
    exp(-pmax(0, time_days - params$t_peak) / 70))

  serum_alb <- pmin(4.5,
    baseline_alb + (4.0 - baseline_alb) *
    params$max_eff * pmin(1, time_days / params$t_peak) *
    exp(-pmax(0, time_days - params$t_peak) / 80))

  cholesterol <- pmax(130,
    350 - 180 * params$max_eff * pmin(1, time_days / params$t_peak) *
    exp(-pmax(0, time_days - params$t_peak) / 75))

  edema_score <- pmax(0,
    (3 + sev) - 3.5 * params$max_eff * pmin(1, time_days / params$t_peak) *
    exp(-pmax(0, time_days - params$t_peak) / 50))
  edema_score <- pmin(4, edema_score)

  # Time to complete remission (UPCR < 0.3)
  cr_idx <- which(uprot_cr < 0.3)
  t_cr   <- if (length(cr_idx) > 0) time_days[min(cr_idx)] else NA

  # ---- Biomarkers ----
  c3 <- pmin(120,
    (70 - sev * 20) + 40 * params$max_eff * pmin(1, time_days / params$t_peak))
  c4 <- pmin(30,
    (15 - sev * 5)  + 12 * params$max_eff * pmin(1, time_days / params$t_peak))
  ige <- pmax(30,
    (200 + sev * 300) - 180 * params$max_eff * pmin(1, time_days / params$t_peak) *
    exp(-pmax(0, time_days - params$t_peak) / 90))
  eosinophil <- pmax(1,
    (5 + sev * 8) - 7 * params$max_eff * pmin(1, time_days / (params$t_peak * 0.8)))
  u_nephrin <- pmax(0,
    (sev * 3) - 2.5 * params$max_eff * pmin(1, time_days / params$t_peak) *
    exp(-pmax(0, time_days - params$t_peak) / 60))

  list(
    time             = time_days,
    cd4_eff          = cd4_eff,
    treg             = treg,
    treg_eff_ratio   = treg_eff_ratio,
    bcell            = bcell_depletion,
    supar            = supar,
    nephrin_expr     = nephrin_expr,
    anti_nephrin_ab  = anti_nephrin_ab,
    podocyte_integrity = podocyte_integrity,
    foot_process_eff = foot_process_eff,
    uprot_cr         = uprot_cr,
    serum_alb        = serum_alb,
    cholesterol      = cholesterol,
    edema_score      = edema_score,
    t_cr             = t_cr,
    c3               = c3,
    c4               = c4,
    ige              = ige,
    eosinophil       = eosinophil,
    u_nephrin        = u_nephrin,
    params           = params
  )
}

# ==============================================================================
# UI Definition
# ==============================================================================

ui <- dashboardPage(
  skin = "blue",

  dashboardHeader(
    title = span(icon("kidneys", lib = "font-awesome"),
                 "MCD QSP Model",
                 style = "font-size:16px; font-weight:bold;")
  ),

  dashboardSidebar(
    width = 280,

    # ---- Global Patient Profile ----
    tags$div(
      style = "padding: 10px 15px 5px; color: #aaa; font-size: 12px; text-transform: uppercase;",
      icon("user-circle"), " 환자 기본 설정 (Global)"
    ),
    sliderInput("baseline_uprot",
      "기저 단백뇨 (g/g Cr)",
      min = 1, max = 20, value = 8, step = 0.5),
    sliderInput("baseline_alb",
      "기저 혈청 알부민 (g/dL)",
      min = 1.5, max = 4.5, value = 2.5, step = 0.1),
    sliderInput("egfr",
      "기저 eGFR (mL/min/1.73m²)",
      min = 15, max = 120, value = 85, step = 5),
    sliderInput("age",
      "나이 (세)",
      min = 1, max = 80, value = 10, step = 1),
    radioGroupButtons("disease_type",
      "질환 유형",
      choices  = c("원발성" = "primary", "이차성" = "secondary"),
      selected = "primary",
      size     = "sm",
      justified = TRUE),
    pickerInput("relapse_hist",
      "재발 이력",
      choices = list(
        "초발"        = "first_episode",
        "빈번 재발형"  = "frequent_relapser",
        "스테로이드 의존형" = "steroid_dependent",
        "스테로이드 저항형" = "steroid_resistant"
      ),
      selected = "first_episode"),

    tags$hr(style = "border-color:#3c4b64; margin:8px 15px;"),

    # ---- Drug Settings ----
    tags$div(
      style = "padding: 5px 15px 5px; color: #aaa; font-size: 12px; text-transform: uppercase;",
      icon("pills"), " 약물 설정"
    ),
    pickerInput("drug_select",
      "치료약물 선택",
      choices = c("Prednisolone", "Cyclosporine", "Tacrolimus", "Rituximab"),
      selected = "Prednisolone"),
    conditionalPanel(
      condition = "input.drug_select != 'Rituximab'",
      sliderInput("dose_daily",
        "1일 용량 (mg/kg/day)",
        min = 0.1, max = 2.0, value = 1.0, step = 0.1)
    ),
    conditionalPanel(
      condition = "input.drug_select == 'Rituximab'",
      sliderInput("dose_rtx",
        "Rituximab 용량 (mg/m²)",
        min = 100, max = 750, value = 375, step = 25)
    ),
    sliderInput("sim_days",
      "시뮬레이션 기간 (일)",
      min = 30, max = 360, value = 180, step = 10),

    tags$hr(style = "border-color:#3c4b64; margin:8px 15px;"),
    actionBttn("run_sim",
      "시뮬레이션 실행",
      icon  = icon("play"),
      style = "material-flat",
      color = "primary",
      block = TRUE)
  ),

  # ==============================================================================
  # Dashboard Body
  # ==============================================================================
  dashboardBody(
    tags$head(
      tags$style(HTML("
        .content-wrapper, .right-side { background-color: #f4f6f9; }
        .nav-tabs-custom { margin-bottom: 0; }
        .box { border-top-color: #3c8dbc; }
        .clinical-note {
          background: #fffbea;
          border-left: 4px solid #f0ad4e;
          padding: 10px 14px;
          border-radius: 3px;
          font-size: 13px;
          margin-top: 10px;
        }
        .severity-badge {
          display: inline-block;
          padding: 4px 10px;
          border-radius: 12px;
          font-weight: bold;
          font-size: 13px;
        }
        .kpi-box {
          text-align: center;
          padding: 14px;
          background: white;
          border-radius: 8px;
          box-shadow: 0 1px 3px rgba(0,0,0,0.1);
        }
        .kpi-value { font-size: 28px; font-weight: bold; }
        .kpi-label { font-size: 12px; color: #666; margin-top: 4px; }
      "))
    ),

    tabBox(
      id = "main_tabs",
      width = 12,
      height = "auto",

      # ============================================================
      # Tab 1: Patient Profile
      # ============================================================
      tabPanel(
        title = tagList(icon("user"), " 환자 프로파일"),
        value = "tab_patient",

        fluidRow(
          box(
            title = "신증후군 중증도 분석", status = "primary",
            solidHeader = TRUE, width = 6,
            plotlyOutput("radar_severity", height = "380px"),
            tags$div(class = "clinical-note",
              icon("info-circle"),
              " 중증도 레이더 차트: 각 축은 정규화된 임상 지표(0–1)를 나타냅니다.",
              "단백뇨·저알부민혈증·부종이 치료 반응의 핵심 예측인자입니다."
            )
          ),
          box(
            title = "위험도 계층화 패널", status = "warning",
            solidHeader = TRUE, width = 6,
            uiOutput("risk_stratification_panel"),
            br(),
            plotlyOutput("patient_kpi_bars", height = "220px")
          )
        ),
        fluidRow(
          box(
            title = "병력 요약 (Clinical Summary)", status = "info",
            solidHeader = TRUE, width = 12, collapsible = TRUE,
            DT::DTOutput("patient_summary_table")
          )
        )
      ),

      # ============================================================
      # Tab 2: PK Profiles
      # ============================================================
      tabPanel(
        title = tagList(icon("chart-line"), " PK 프로파일"),
        value = "tab_pk",

        fluidRow(
          box(
            title = "혈중 농도-시간 곡선 (Concentration-Time Profile)",
            status = "primary", solidHeader = TRUE, width = 8,
            plotlyOutput("pk_conc_time", height = "420px")
          ),
          box(
            title = "PK 요약 지표", status = "info",
            solidHeader = TRUE, width = 4,
            uiOutput("pk_summary_kpi"),
            br(),
            tags$div(class = "clinical-note",
              icon("exclamation-triangle"),
              strong(" 치료 모니터링 권고"),
              br(),
              "CNI 투여 시 혈중 최저 농도(trough level)를 2주마다 모니터링하십시오.",
              br(),
              "CsA 목표 trough: 100–200 ng/mL",
              br(),
              "Tacrolimus 목표 trough: 5–10 ng/mL"
            )
          )
        ),
        fluidRow(
          box(
            title = "용량-반응 관계 (Dose-Response)", status = "success",
            solidHeader = TRUE, width = 12, collapsible = TRUE,
            plotlyOutput("pk_dose_response", height = "280px")
          )
        )
      ),

      # ============================================================
      # Tab 3: Immune Dynamics
      # ============================================================
      tabPanel(
        title = tagList(icon("shield-virus"), " 면역 역학"),
        value = "tab_immune",

        fluidRow(
          box(
            title = "CD4 효과 T세포 및 Treg 궤적",
            status = "danger", solidHeader = TRUE, width = 7,
            plotlyOutput("immune_trajectory", height = "380px")
          ),
          box(
            title = "Treg/Effector 비율 및 B세포 소진",
            status = "warning", solidHeader = TRUE, width = 5,
            plotlyOutput("treg_ratio_plot", height = "180px"),
            plotlyOutput("bcell_plot", height = "180px")
          )
        ),
        fluidRow(
          box(
            title = "면역세포 구성 (누적 면적 차트)",
            status = "primary", solidHeader = TRUE, width = 8,
            plotlyOutput("immune_stacked_area", height = "300px"),
            tags$div(class = "clinical-note",
              icon("info-circle"),
              " 음영 처리된 구간은 치료 초기 면역 전환기를 나타냅니다.",
              "Treg 비율 증가는 지속적 관해의 중요 예측인자입니다."
            )
          ),
          box(
            title = "suPAR (가용성 uPA 수용체)", status = "info",
            solidHeader = TRUE, width = 4,
            plotlyOutput("supar_plot", height = "240px"),
            tags$div(class = "clinical-note",
              "suPAR > 3 ng/mL는 족세포 손상 및 단백뇨 지속과 관련됩니다."
            )
          )
        )
      ),

      # ============================================================
      # Tab 4: Podocyte Biology
      # ============================================================
      tabPanel(
        title = tagList(icon("microscope"), " 족세포 생물학"),
        value = "tab_podocyte",

        fluidRow(
          box(
            title = "족세포 완전성 지수 및 분자 마커",
            status = "danger", solidHeader = TRUE, width = 7,
            plotlyOutput("podocyte_dynamics", height = "380px")
          ),
          box(
            title = "족돌기 소실 중증도 게이지",
            status = "warning", solidHeader = TRUE, width = 5,
            plotlyOutput("foot_process_gauge", height = "180px"),
            plotlyOutput("anti_nephrin_plot", height = "180px")
          )
        ),
        fluidRow(
          box(
            title = "족세포 분자 마커 히트맵",
            status = "primary", solidHeader = TRUE, width = 12,
            plotlyOutput("podocyte_heatmap", height = "280px"),
            tags$div(class = "clinical-note",
              icon("dna"),
              " Nephrin·Podocin 발현 감소와 Desmin·CD2AP 상향조절은",
              "족세포 탈분화 및 슬릿막 손상을 반영합니다."
            )
          )
        )
      ),

      # ============================================================
      # Tab 5: Clinical Endpoints
      # ============================================================
      tabPanel(
        title = tagList(icon("hospital"), " 임상 엔드포인트"),
        value = "tab_clinical",

        fluidRow(
          valueBoxOutput("vbox_uprot",     width = 3),
          valueBoxOutput("vbox_alb",       width = 3),
          valueBoxOutput("vbox_chol",      width = 3),
          valueBoxOutput("vbox_t_cr",      width = 3)
        ),
        fluidRow(
          box(
            title = "뇨단백/크레아티닌 비(UPCR) 및 혈청 알부민 추이",
            status = "primary", solidHeader = TRUE, width = 8,
            plotlyOutput("clinical_primary", height = "380px")
          ),
          box(
            title = "완전 관해(CR) 달성 상태", status = "success",
            solidHeader = TRUE, width = 4,
            plotlyOutput("cr_status_plot", height = "200px"),
            br(),
            tags$div(class = "clinical-note",
              icon("check-circle"),
              strong(" 완전 관해 정의"),
              tags$ul(
                tags$li("UPCR < 0.3 g/g Cr (성인)"),
                tags$li("UPCR < 0.2 g/g Cr (소아)"),
                tags$li("혈청 알부민 ≥ 3.5 g/dL")
              )
            )
          )
        ),
        fluidRow(
          box(
            title = "부종 점수 및 혈청 콜레스테롤",
            status = "warning", solidHeader = TRUE, width = 12,
            plotlyOutput("clinical_secondary", height = "260px")
          )
        )
      ),

      # ============================================================
      # Tab 6: Scenario Comparison
      # ============================================================
      tabPanel(
        title = tagList(icon("balance-scale"), " 시나리오 비교"),
        value = "tab_scenario",

        fluidRow(
          box(
            title = "6가지 치료 시나리오 비교", status = "primary",
            solidHeader = TRUE, width = 12,
            fluidRow(
              column(4, sliderInput("scen_uprot", "기저 단백뇨 (시나리오용, g/g)",
                                   min = 1, max = 20, value = 8, step = 0.5)),
              column(4, sliderInput("scen_alb",   "기저 알부민 (시나리오용, g/dL)",
                                   min = 1.5, max = 4.5, value = 2.5, step = 0.1)),
              column(4, radioGroupButtons("scen_relapse",
                "재발 유형 (시나리오)",
                choices  = c("초발" = "first_episode",
                             "빈번 재발" = "frequent_relapser",
                             "스테로이드 의존" = "steroid_dependent"),
                selected = "first_episode", size = "sm", justified = TRUE))
            )
          )
        ),
        fluidRow(
          box(
            title = "Day-90 단백뇨 감소율 (%)",
            status = "success", solidHeader = TRUE, width = 6,
            plotlyOutput("scen_bar_uprot", height = "320px")
          ),
          box(
            title = "완전 관해율 (Day-90 / Day-180)",
            status = "info", solidHeader = TRUE, width = 6,
            plotlyOutput("scen_bar_cr", height = "320px")
          )
        ),
        fluidRow(
          box(
            title = "첫 관해까지 소요 시간 (일)",
            status = "warning", solidHeader = TRUE, width = 6,
            plotlyOutput("scen_ttcr", height = "280px")
          ),
          box(
            title = "치료 효능/안전성 레이더 (종합 프로파일)",
            status = "danger", solidHeader = TRUE, width = 6,
            plotlyOutput("scen_radar", height = "280px")
          )
        ),
        fluidRow(
          box(
            title = "시나리오별 상세 비교표", status = "primary",
            solidHeader = TRUE, width = 12, collapsible = TRUE,
            DT::DTOutput("scen_table")
          )
        )
      ),

      # ============================================================
      # Tab 7: Biomarker Panel
      # ============================================================
      tabPanel(
        title = tagList(icon("vial"), " 바이오마커"),
        value = "tab_biomarker",

        fluidRow(
          box(
            title = "항-Nephrin 항체 및 뇨중 Nephrin", status = "danger",
            solidHeader = TRUE, width = 6,
            plotlyOutput("bm_nephrin", height = "300px"),
            tags$div(class = "clinical-note",
              icon("flask"),
              " 항-Nephrin IgG 항체는 MCD 활성기에 상승하며",
              "족세포 슬릿막 파괴와 직접 연관됩니다 (Beck 2023, NEJM)."
            )
          ),
          box(
            title = "보체 (C3/C4) 및 IgE / 호산구",
            status = "warning", solidHeader = TRUE, width = 6,
            plotlyOutput("bm_complement_ige", height = "300px")
          )
        ),
        fluidRow(
          box(
            title = "바이오마커 상관관계 히트맵",
            status = "primary", solidHeader = TRUE, width = 7,
            plotlyOutput("bm_correlation_heatmap", height = "320px")
          ),
          box(
            title = "바이오마커 기반 치료 결정 트리",
            status = "info", solidHeader = TRUE, width = 5,
            plotlyOutput("bm_decision_tree", height = "320px"),
            tags$div(class = "clinical-note",
              icon("sitemap"),
              " 바이오마커 프로파일에 따른 맞춤형 치료 전략 선택 가이드."
            )
          )
        ),
        fluidRow(
          box(
            title = "바이오마커 시계열 데이터", status = "success",
            solidHeader = TRUE, width = 12, collapsible = TRUE,
            DT::DTOutput("bm_table")
          )
        )
      )
    ) # end tabBox
  ) # end dashboardBody
) # end dashboardPage

# ==============================================================================
# Server Logic
# ==============================================================================

server <- function(input, output, session) {

  # ---- Reactive: Run simulation ----
  sim_data <- eventReactive(input$run_sim, {
    time_days <- seq(0, input$sim_days, by = 1)
    dose_val  <- if (input$drug_select == "Rituximab") input$dose_rtx
                 else input$dose_daily
    run_pd_model(
      time_days      = time_days,
      drug           = input$drug_select,
      dose           = dose_val,
      baseline_uprot = input$baseline_uprot,
      baseline_alb   = input$baseline_alb,
      age            = input$age,
      disease_type   = input$disease_type,
      relapse_hist   = input$relapse_hist
    )
  }, ignoreNULL = FALSE)

  # ====================================================================
  # Tab 1: Patient Profile
  # ====================================================================

  output$radar_severity <- renderPlotly({
    alb_norm   <- 1 - (input$baseline_alb - 1.5) / 3.0
    uprot_norm <- (input$baseline_uprot - 1) / 19
    egfr_norm  <- 1 - (input$egfr - 15) / 105
    age_risk   <- if (input$age < 10 || input$age > 65) 0.7 else 0.4
    relapse_n  <- switch(input$relapse_hist,
      first_episode = 0.2, frequent_relapser = 0.6,
      steroid_dependent = 0.8, steroid_resistant = 1.0)

    theta_vals <- c("단백뇨 부담", "저알부민혈증",
                    "신기능 저하", "연령 위험도",
                    "재발 이력", "단백뇨 부담")
    r_vals     <- c(uprot_norm, alb_norm, egfr_norm,
                    age_risk, relapse_n, uprot_norm)

    plot_ly(
      type = "scatterpolar",
      r    = r_vals,
      theta = theta_vals,
      fill  = "toself",
      fillcolor = "rgba(60, 141, 188, 0.35)",
      line  = list(color = "#3c8dbc", width = 2),
      name  = "환자 프로파일"
    ) %>%
      layout(
        polar = list(
          radialaxis = list(visible = TRUE, range = c(0, 1),
                            tickfont = list(size = 10))
        ),
        showlegend = FALSE,
        margin = list(t = 30, b = 30, l = 50, r = 50),
        paper_bgcolor = "white"
      )
  })

  output$risk_stratification_panel <- renderUI({
    sev_score <- ((input$baseline_uprot / 20) * 0.4 +
                  (1 - input$baseline_alb / 4.5) * 0.35 +
                  (1 - (input$egfr - 15) / 105) * 0.25)
    risk_level <- if (sev_score < 0.35) "저위험 (Low Risk)"
                  else if (sev_score < 0.60) "중등위험 (Moderate Risk)"
                  else "고위험 (High Risk)"
    risk_color <- if (sev_score < 0.35) "#27ae60"
                  else if (sev_score < 0.60) "#f39c12"
                  else "#e74c3c"

    ns_stage <- if (input$baseline_uprot > 3.5 && input$baseline_alb < 3.5)
      "완전 신증후군 (Full NS)" else "불완전 신증후군 (Incomplete NS)"

    tagList(
      tags$div(style = "margin-bottom: 12px;",
        tags$span(class = "kpi-value",
          style = paste0("color:", risk_color, "; font-size:22px;"),
          icon("exclamation-circle"), " ", risk_level)
      ),
      tags$div(class = "kpi-box", style = "margin-bottom:8px;",
        tags$div(class = "kpi-label", "신증후군 병기"),
        tags$div(class = "kpi-value", style = "font-size:18px;", ns_stage)
      ),
      fluidRow(
        column(6, tags$div(class = "kpi-box",
          tags$div(class = "kpi-label", "중증도 점수"),
          tags$div(class = "kpi-value",
            style = paste0("color:", risk_color, ";"),
            round(sev_score, 2))
        )),
        column(6, tags$div(class = "kpi-box",
          tags$div(class = "kpi-label", "KDIGO 위험군"),
          tags$div(class = "kpi-value", style = "font-size: 16px;",
            if (sev_score < 0.35) "Class I"
            else if (sev_score < 0.60) "Class II"
            else "Class III")
        ))
      )
    )
  })

  output$patient_kpi_bars <- renderPlotly({
    tgt_uprot <- input$baseline_uprot / 0.3   # ratio vs target
    tgt_alb   <- input$baseline_alb / 3.5
    tgt_egfr  <- input$egfr / 90

    plot_ly() %>%
      add_bars(
        x = c("단백뇨\n(vs 목표)", "알부민\n(vs 목표)", "eGFR\n(vs 정상)"),
        y = c(min(tgt_uprot, 5), tgt_alb, tgt_egfr),
        marker = list(
          color = c(
            ifelse(input$baseline_uprot < 0.3, "#27ae60", "#e74c3c"),
            ifelse(input$baseline_alb > 3.5, "#27ae60", "#e74c3c"),
            ifelse(input$egfr > 60, "#27ae60", "#f39c12")
          )
        )
      ) %>%
      add_segments(x = 0.5, xend = 3.5, y = 1, yend = 1,
                   line = list(dash = "dash", color = "black", width = 1)) %>%
      layout(
        xaxis = list(title = ""),
        yaxis = list(title = "목표 대비 비율"),
        showlegend = FALSE,
        margin = list(t = 10, b = 10)
      )
  })

  output$patient_summary_table <- DT::renderDT({
    df <- data.frame(
      "항목"  = c("기저 단백뇨 (g/g Cr)", "혈청 알부민 (g/dL)",
                  "eGFR (mL/min/1.73m²)", "나이 (세)", "질환 유형",
                  "재발 이력", "치료 약물", "시뮬레이션 기간 (일)"),
      "현재값" = c(input$baseline_uprot, input$baseline_alb,
                   input$egfr, input$age,
                   if (input$disease_type == "primary") "원발성 MCD" else "이차성 MCD",
                   switch(input$relapse_hist,
                     first_episode = "초발 에피소드",
                     frequent_relapser = "빈번 재발형",
                     steroid_dependent = "스테로이드 의존형",
                     steroid_resistant = "스테로이드 저항형"),
                   input$drug_select,
                   input$sim_days),
      "정상 범위" = c("< 0.3", "3.5–5.0", "≥ 60", "–", "–", "–", "–", "–"),
      "상태" = c(
        ifelse(input$baseline_uprot < 0.3, "정상", "비정상"),
        ifelse(input$baseline_alb >= 3.5, "정상", "저하"),
        ifelse(input$egfr >= 60, "정상", "저하"),
        "–", "–", "–", "–", "–"
      )
    )
    DT::datatable(df, options = list(dom = "t", pageLength = 10),
                  rownames = FALSE)
  })

  # ====================================================================
  # Tab 2: PK Profiles
  # ====================================================================

  output$pk_conc_time <- renderPlotly({
    req(sim_data())
    d   <- sim_data()
    t   <- d$time
    prm <- d$params

    dose_val <- if (input$drug_select == "Rituximab") input$dose_rtx
                else input$dose_daily

    conc <- if (input$drug_select == "Rituximab") {
      pk_twocomp_iv(t, dose_val * 1.73, prm$CL, prm$Vc, prm$Q, prm$Vp)
    } else {
      pk_onecomp(t, dose_val * 60, prm$CL, prm$Vd, prm$ka, tau = prm$tau)
    }

    unit_label <- switch(input$drug_select,
      Prednisolone = "ng/mL",
      Cyclosporine = "ng/mL (×10)",
      Tacrolimus   = "ng/mL",
      Rituximab    = "μg/mL"
    )

    p <- plot_ly() %>%
      add_lines(x = t, y = conc,
                name = input$drug_select,
                line = list(color = "#3c8dbc", width = 2.5))

    # Therapeutic window bands for CNIs
    if (!is.na(prm$trough_therapeutic[1])) {
      p <- p %>%
        add_ribbons(x = t,
                    ymin = rep(prm$trough_therapeutic[1], length(t)),
                    ymax = rep(prm$trough_therapeutic[2], length(t)),
                    name = "치료 목표 범위",
                    fillcolor = "rgba(39, 174, 96, 0.15)",
                    line = list(color = "transparent"))
    }

    p %>% layout(
      xaxis = list(title = "시간 (일)"),
      yaxis = list(title = paste0("혈중 농도 (", unit_label, ")")),
      legend = list(orientation = "h", y = -0.15),
      hovermode = "x unified",
      margin = list(t = 20)
    )
  })

  output$pk_summary_kpi <- renderUI({
    req(sim_data())
    d        <- sim_data()
    prm      <- d$params
    dose_val <- if (input$drug_select == "Rituximab") input$dose_rtx
                else input$dose_daily
    t        <- d$time

    conc <- if (input$drug_select == "Rituximab") {
      pk_twocomp_iv(t, dose_val * 1.73, prm$CL, prm$Vc, prm$Q, prm$Vp)
    } else {
      pk_onecomp(t, dose_val * 60, prm$CL, prm$Vd, prm$ka, tau = prm$tau)
    }

    cmax <- round(max(conc), 2)
    auc  <- round(sum(diff(t) * (head(conc, -1) + tail(conc, -1)) / 2), 0)

    tagList(
      fluidRow(
        column(6, tags$div(class = "kpi-box",
          tags$div(class = "kpi-label", "Cmax"),
          tags$div(class = "kpi-value", style = "font-size:20px;", cmax)
        )),
        column(6, tags$div(class = "kpi-box",
          tags$div(class = "kpi-label", "AUC (0–t)"),
          tags$div(class = "kpi-value", style = "font-size:20px;", auc)
        ))
      ),
      br(),
      fluidRow(
        column(12, tags$div(class = "kpi-box",
          tags$div(class = "kpi-label", "약물"),
          tags$div(class = "kpi-value", style = "font-size:18px;",
                   input$drug_select)
        ))
      )
    )
  })

  output$pk_dose_response <- renderPlotly({
    doses <- seq(0.1, 2.0, by = 0.1)
    prot_red <- sapply(doses, function(d) {
      max_eff <- switch(input$drug_select,
        Prednisolone = 0.85, Cyclosporine = 0.72,
        Tacrolimus = 0.75, Rituximab = 0.80)
      max_eff * d / (0.3 + d)
    })
    plot_ly(x = doses, y = prot_red * 100,
            type = "scatter", mode = "lines+markers",
            line = list(color = "#e74c3c", width = 2)) %>%
      layout(
        xaxis = list(title = "용량 (mg/kg/day)"),
        yaxis = list(title = "단백뇨 감소율 (%)"),
        hovermode = "x"
      )
  })

  # ====================================================================
  # Tab 3: Immune Dynamics
  # ====================================================================

  output$immune_trajectory <- renderPlotly({
    req(sim_data())
    d <- sim_data()
    plot_ly() %>%
      add_lines(x = d$time, y = d$cd4_eff,
                name = "CD4 효과 T세포",
                line = list(color = "#e74c3c", width = 2)) %>%
      add_lines(x = d$time, y = d$treg,
                name = "Treg 세포",
                line = list(color = "#27ae60", width = 2)) %>%
      add_ribbons(x = c(0, 14),
                  ymin = c(0, 0), ymax = c(2, 2),
                  name = "치료 초기 전환기",
                  fillcolor = "rgba(243,156,18,0.12)",
                  line = list(color = "transparent"),
                  showlegend = TRUE) %>%
      layout(
        xaxis = list(title = "시간 (일)"),
        yaxis = list(title = "상대적 세포 수 (기저치 = 1)"),
        legend = list(orientation = "h", y = -0.15),
        hovermode = "x unified"
      )
  })

  output$treg_ratio_plot <- renderPlotly({
    req(sim_data())
    d <- sim_data()
    plot_ly(x = d$time, y = d$treg_eff_ratio,
            type = "scatter", mode = "lines",
            fill = "tozeroy",
            fillcolor = "rgba(39,174,96,0.2)",
            line = list(color = "#27ae60")) %>%
      layout(
        xaxis = list(title = ""),
        yaxis = list(title = "Treg/Effector 비율"),
        showlegend = FALSE
      )
  })

  output$bcell_plot <- renderPlotly({
    req(sim_data())
    d <- sim_data()
    plot_ly(x = d$time, y = d$bcell,
            type = "scatter", mode = "lines",
            fill = "tozeroy",
            fillcolor = "rgba(52,152,219,0.2)",
            line = list(color = "#3498db")) %>%
      layout(
        xaxis = list(title = "시간 (일)"),
        yaxis = list(title = "B세포 (상대치)"),
        showlegend = FALSE
      )
  })

  output$immune_stacked_area <- renderPlotly({
    req(sim_data())
    d <- sim_data()
    nk_cells <- pmax(0.3, 0.8 - 0.3 * d$cd4_eff + 0.15 * d$treg)
    cd8      <- pmax(0.2, 0.9 - 0.4 * (1 - d$cd4_eff))

    plot_ly() %>%
      add_lines(x = d$time, y = d$cd4_eff + d$treg + d$bcell + nk_cells + cd8,
                name = "CD8 T세포", fill = "tozeroy",
                fillcolor = "rgba(231,76,60,0.6)",
                line = list(color = "rgba(231,76,60,0.8)")) %>%
      add_lines(x = d$time, y = d$cd4_eff + d$treg + d$bcell + nk_cells,
                name = "NK세포", fill = "tozeroy",
                fillcolor = "rgba(155,89,182,0.6)",
                line = list(color = "rgba(155,89,182,0.8)")) %>%
      add_lines(x = d$time, y = d$cd4_eff + d$treg + d$bcell,
                name = "B세포", fill = "tozeroy",
                fillcolor = "rgba(52,152,219,0.6)",
                line = list(color = "rgba(52,152,219,0.8)")) %>%
      add_lines(x = d$time, y = d$cd4_eff + d$treg,
                name = "Treg", fill = "tozeroy",
                fillcolor = "rgba(39,174,96,0.6)",
                line = list(color = "rgba(39,174,96,0.8)")) %>%
      add_lines(x = d$time, y = d$cd4_eff,
                name = "CD4 Effector", fill = "tozeroy",
                fillcolor = "rgba(231,76,60,0.4)",
                line = list(color = "rgba(231,76,60,0.6)")) %>%
      layout(
        xaxis = list(title = "시간 (일)"),
        yaxis = list(title = "상대적 세포 구성"),
        hovermode = "x unified",
        legend = list(orientation = "h", y = -0.18)
      )
  })

  output$supar_plot <- renderPlotly({
    req(sim_data())
    d <- sim_data()
    plot_ly(x = d$time, y = d$supar,
            type = "scatter", mode = "lines",
            fill = "tozeroy",
            fillcolor = "rgba(230,126,34,0.25)",
            line = list(color = "#e67e22", width = 2)) %>%
      add_segments(x = 0, xend = max(d$time), y = 3, yend = 3,
                   line = list(dash = "dash", color = "red", width = 1),
                   name = "위험 역치 (3 ng/mL)") %>%
      layout(
        xaxis = list(title = "시간 (일)"),
        yaxis = list(title = "suPAR (ng/mL)"),
        showlegend = FALSE
      )
  })

  # ====================================================================
  # Tab 4: Podocyte Biology
  # ====================================================================

  output$podocyte_dynamics <- renderPlotly({
    req(sim_data())
    d <- sim_data()
    plot_ly() %>%
      add_lines(x = d$time, y = d$podocyte_integrity,
                name = "족세포 완전성 지수",
                line = list(color = "#8e44ad", width = 2.5)) %>%
      add_lines(x = d$time, y = d$nephrin_expr,
                name = "Nephrin 발현",
                line = list(color = "#3498db", width = 2)) %>%
      add_lines(x = d$time, y = 1 - d$foot_process_eff,
                name = "족돌기 복원 지수",
                line = list(color = "#27ae60", width = 2)) %>%
      layout(
        xaxis = list(title = "시간 (일)"),
        yaxis = list(title = "지수 (0–1)", range = c(0, 1.1)),
        legend = list(orientation = "h", y = -0.15),
        hovermode = "x unified"
      )
  })

  output$foot_process_gauge <- renderPlotly({
    req(sim_data())
    d     <- sim_data()
    last  <- d$foot_process_eff[length(d$foot_process_eff)]
    color <- if (last < 0.2) "#27ae60" else if (last < 0.5) "#f39c12" else "#e74c3c"
    plot_ly(
      type = "indicator", mode = "gauge+number",
      value = round(last * 100, 1),
      title = list(text = "족돌기 소실 중증도 (%)", font = list(size = 13)),
      gauge = list(
        axis   = list(range = list(0, 100)),
        bar    = list(color = color),
        steps  = list(
          list(range = c(0, 20),  color = "#d5f5e3"),
          list(range = c(20, 50), color = "#fdebd0"),
          list(range = c(50, 100), color = "#fadbd8")
        ),
        threshold = list(
          line  = list(color = "red", width = 4),
          thickness = 0.75,
          value = 70
        )
      )
    ) %>% layout(margin = list(t = 40, b = 10, l = 10, r = 10))
  })

  output$anti_nephrin_plot <- renderPlotly({
    req(sim_data())
    d <- sim_data()
    plot_ly(x = d$time, y = d$anti_nephrin_ab,
            type = "scatter", mode = "lines",
            fill = "tozeroy",
            fillcolor = "rgba(142,68,173,0.2)",
            line = list(color = "#8e44ad", width = 2)) %>%
      layout(
        xaxis = list(title = "시간 (일)"),
        yaxis = list(title = "항-Nephrin 항체 (상대치)"),
        showlegend = FALSE
      )
  })

  output$podocyte_heatmap <- renderPlotly({
    req(sim_data())
    d    <- sim_data()
    tidx <- seq(1, length(d$time), length.out = 10)
    tlab <- round(d$time[tidx])

    nephrin_v  <- d$nephrin_expr[tidx]
    podocin_v  <- pmin(1, nephrin_v + 0.05)
    synapt_v   <- pmin(1, nephrin_v * 0.9)
    desmin_v   <- pmax(0, 1 - nephrin_v)
    cd2ap_v    <- pmin(1, 0.3 + 0.5 * nephrin_v)
    integrin_v <- pmin(1, 0.4 + 0.4 * d$podocyte_integrity[tidx])

    z_mat <- rbind(nephrin_v, podocin_v, synapt_v,
                   desmin_v, cd2ap_v, integrin_v)
    rownames(z_mat) <- c("Nephrin", "Podocin", "Synaptopodin",
                         "Desmin", "CD2AP", "α3β1-Integrin")

    plot_ly(
      z         = z_mat,
      x         = as.character(tlab),
      y         = rownames(z_mat),
      type      = "heatmap",
      colorscale = "RdYlGn",
      zmin       = 0, zmax = 1
    ) %>%
      layout(
        xaxis = list(title = "시간 (일)"),
        yaxis = list(title = ""),
        margin = list(l = 110)
      )
  })

  # ====================================================================
  # Tab 5: Clinical Endpoints
  # ====================================================================

  output$vbox_uprot <- renderValueBox({
    req(sim_data())
    d    <- sim_data()
    last <- round(tail(d$uprot_cr, 1), 2)
    valueBox(
      paste0(last, " g/g"),
      "최종 UPCR",
      icon  = icon("tint"),
      color = if (last < 0.3) "green" else if (last < 1.0) "yellow" else "red"
    )
  })

  output$vbox_alb <- renderValueBox({
    req(sim_data())
    d    <- sim_data()
    last <- round(tail(d$serum_alb, 1), 2)
    valueBox(
      paste0(last, " g/dL"),
      "최종 혈청 알부민",
      icon  = icon("flask"),
      color = if (last >= 3.5) "green" else if (last >= 2.5) "yellow" else "red"
    )
  })

  output$vbox_chol <- renderValueBox({
    req(sim_data())
    d    <- sim_data()
    last <- round(tail(d$cholesterol, 1), 0)
    valueBox(
      paste0(last, " mg/dL"),
      "최종 총 콜레스테롤",
      icon  = icon("heartbeat"),
      color = if (last < 200) "green" else if (last < 240) "yellow" else "red"
    )
  })

  output$vbox_t_cr <- renderValueBox({
    req(sim_data())
    d    <- sim_data()
    t_cr <- d$t_cr
    valueBox(
      if (is.na(t_cr)) "미달성" else paste0(t_cr, "일"),
      "완전 관해 달성일",
      icon  = icon("check-circle"),
      color = if (!is.na(t_cr) && t_cr <= 90) "green"
              else if (!is.na(t_cr)) "yellow"
              else "red"
    )
  })

  output$clinical_primary <- renderPlotly({
    req(sim_data())
    d <- sim_data()
    ax2_scale <- 4.5 / max(d$uprot_cr, na.rm = TRUE)

    p <- plot_ly() %>%
      add_lines(x = d$time, y = d$uprot_cr,
                name = "UPCR (g/g Cr)",
                line = list(color = "#e74c3c", width = 2.5)) %>%
      add_segments(x = 0, xend = max(d$time), y = 0.3, yend = 0.3,
                   line = list(dash = "dash", color = "red", width = 1),
                   name = "CR 목표 (0.3)") %>%
      add_lines(x = d$time, y = d$serum_alb,
                name = "혈청 알부민 (g/dL)",
                yaxis = "y2",
                line = list(color = "#3498db", width = 2.5)) %>%
      add_segments(x = 0, xend = max(d$time), y = 3.5, yend = 3.5,
                   yaxis = "y2",
                   line = list(dash = "dash", color = "blue", width = 1),
                   name = "정상 알부민 (3.5)")

    if (!is.na(d$t_cr)) {
      p <- p %>% add_segments(
        x = d$t_cr, xend = d$t_cr,
        y = 0, yend = max(d$uprot_cr),
        line = list(dash = "dot", color = "darkgreen", width = 2),
        name = paste0("CR 달성 (Day ", d$t_cr, ")")
      )
    }

    p %>% layout(
      xaxis  = list(title = "시간 (일)"),
      yaxis  = list(title = "UPCR (g/g Cr)", side = "left"),
      yaxis2 = list(title = "혈청 알부민 (g/dL)", side = "right",
                    overlaying = "y"),
      legend = list(orientation = "h", y = -0.15),
      hovermode = "x unified"
    )
  })

  output$cr_status_plot <- renderPlotly({
    req(sim_data())
    d       <- sim_data()
    cr_days <- sum(d$uprot_cr < 0.3)
    total   <- length(d$time)

    plot_ly(
      labels = c("관해 기간", "미관해 기간"),
      values = c(cr_days, total - cr_days),
      type   = "pie",
      hole   = 0.5,
      marker = list(colors = c("#27ae60", "#e74c3c"))
    ) %>%
      layout(
        showlegend = TRUE,
        legend = list(orientation = "h", y = -0.1),
        margin = list(t = 10, b = 10)
      )
  })

  output$clinical_secondary <- renderPlotly({
    req(sim_data())
    d <- sim_data()
    plot_ly() %>%
      add_lines(x = d$time, y = d$cholesterol,
                name = "총 콜레스테롤 (mg/dL)",
                line = list(color = "#e67e22", width = 2)) %>%
      add_segments(x = 0, xend = max(d$time), y = 200, yend = 200,
                   line = list(dash = "dash", color = "#e67e22"),
                   name = "정상 콜레스테롤 (200)") %>%
      add_lines(x = d$time, y = d$edema_score * 50,
                name = "부종 점수 (×50)",
                yaxis = "y2",
                line = list(color = "#9b59b6", width = 2)) %>%
      layout(
        xaxis  = list(title = "시간 (일)"),
        yaxis  = list(title = "총 콜레스테롤 (mg/dL)"),
        yaxis2 = list(title = "부종 점수 (×50)", side = "right",
                      overlaying = "y"),
        legend = list(orientation = "h", y = -0.15),
        hovermode = "x unified"
      )
  })

  # ====================================================================
  # Tab 6: Scenario Comparison
  # ====================================================================

  scenario_results <- reactive({
    t <- seq(0, 180, by = 1)
    scenarios <- list(
      list(drug = "Prednisolone", dose = 1.0,  label = "Prednisolone\n표준 용량"),
      list(drug = "Prednisolone", dose = 0.5,  label = "Prednisolone\n저용량"),
      list(drug = "Cyclosporine", dose = 0.8,  label = "Cyclosporine"),
      list(drug = "Tacrolimus",   dose = 0.1,  label = "Tacrolimus"),
      list(drug = "Rituximab",    dose = 375,  label = "Rituximab\n375 mg/m²"),
      list(drug = "Rituximab",    dose = 750,  label = "Rituximab\n고용량")
    )

    lapply(scenarios, function(s) {
      res <- run_pd_model(
        time_days      = t,
        drug           = s$drug,
        dose           = s$dose,
        baseline_uprot = input$scen_uprot,
        baseline_alb   = input$scen_alb,
        relapse_hist   = input$scen_relapse
      )
      baseline_up <- res$uprot_cr[1]
      idx90  <- which(t == 90)
      idx180 <- which(t == 180)

      list(
        label        = s$label,
        drug         = s$drug,
        prot_red_90  = (1 - res$uprot_cr[idx90]  / baseline_up) * 100,
        cr_rate_90   = as.integer(res$uprot_cr[idx90]  < 0.3) * 100,
        cr_rate_180  = as.integer(res$uprot_cr[idx180] < 0.3) * 100,
        t_cr         = if (is.na(res$t_cr)) 180 else res$t_cr,
        efficacy     = (1 - res$uprot_cr[idx180] / baseline_up) * 100,
        speed        = 100 - (if (is.na(res$t_cr)) 180 else res$t_cr) / 180 * 100,
        alb_recov    = (res$serum_alb[idx180] / 4.0) * 100,
        safety       = switch(s$drug,
          Prednisolone = 55, Cyclosporine = 65,
          Tacrolimus = 70, Rituximab = 80)
      )
    })
  })

  scen_colors <- c("#3498db", "#85c1e9", "#e74c3c", "#f39c12", "#27ae60", "#1e8449")

  output$scen_bar_uprot <- renderPlotly({
    res <- scenario_results()
    labels <- sapply(res, function(x) x$label)
    vals   <- sapply(res, function(x) x$prot_red_90)
    plot_ly(x = labels, y = vals, type = "bar",
            marker = list(color = scen_colors)) %>%
      layout(
        xaxis = list(title = "", tickangle = -20),
        yaxis = list(title = "단백뇨 감소율 (%)"),
        showlegend = FALSE
      )
  })

  output$scen_bar_cr <- renderPlotly({
    res    <- scenario_results()
    labels <- sapply(res, function(x) x$label)
    cr90   <- sapply(res, function(x) x$cr_rate_90)
    cr180  <- sapply(res, function(x) x$cr_rate_180)

    plot_ly() %>%
      add_bars(x = labels, y = cr90,  name = "Day-90 CR",
               marker = list(color = scen_colors)) %>%
      add_bars(x = labels, y = cr180, name = "Day-180 CR",
               marker = list(color = scen_colors,
                             opacity = 0.55)) %>%
      layout(
        barmode = "group",
        xaxis   = list(title = "", tickangle = -20),
        yaxis   = list(title = "완전 관해율 (%)", range = c(0, 110)),
        legend  = list(orientation = "h", y = -0.2)
      )
  })

  output$scen_ttcr <- renderPlotly({
    res    <- scenario_results()
    labels <- sapply(res, function(x) x$label)
    ttcr   <- sapply(res, function(x) x$t_cr)

    plot_ly(x = labels, y = ttcr, type = "bar",
            marker = list(color = scen_colors)) %>%
      add_segments(x = 0.5, xend = length(labels) + 0.5,
                   y = 90, yend = 90,
                   line = list(dash = "dash", color = "green", width = 1),
                   name = "90일 기준선") %>%
      layout(
        xaxis = list(title = "", tickangle = -20),
        yaxis = list(title = "관해까지 소요 시간 (일)"),
        showlegend = FALSE
      )
  })

  output$scen_radar <- renderPlotly({
    res    <- scenario_results()
    labels <- sapply(res, function(x) x$label)

    dims   <- c("치료 효능", "치료 속도", "알부민 회복", "안전성", "치료 효능")

    p <- plot_ly(type = "scatterpolar")
    for (i in seq_along(res)) {
      r_vals <- c(res[[i]]$efficacy, res[[i]]$speed,
                  res[[i]]$alb_recov, res[[i]]$safety,
                  res[[i]]$efficacy)
      p <- add_trace(p, r = r_vals, theta = dims,
                     fill = "toself", name = labels[i],
                     fillcolor = paste0(sub("^#", "rgba(", scen_colors[i]),
                                        ", 0.15)"))
    }
    p %>% layout(
      polar = list(radialaxis = list(range = c(0, 100))),
      legend = list(orientation = "h", y = -0.2)
    )
  })

  output$scen_table <- DT::renderDT({
    res <- scenario_results()
    df  <- data.frame(
      "치료 시나리오"     = sapply(res, function(x) gsub("\n", " ", x$label)),
      "90일 단백뇨 감소율(%)" = round(sapply(res, function(x) x$prot_red_90), 1),
      "90일 CR율(%)"     = sapply(res, function(x) x$cr_rate_90),
      "180일 CR율(%)"    = sapply(res, function(x) x$cr_rate_180),
      "관해 소요일"       = sapply(res, function(x) x$t_cr),
      "안전성 점수"       = sapply(res, function(x) x$safety),
      check.names = FALSE
    )
    DT::datatable(df, options = list(dom = "t"), rownames = FALSE) %>%
      DT::formatStyle("90일 CR율(%)",
                      backgroundColor = DT::styleInterval(c(0, 50, 100),
                                          c("white", "#fef9e7", "#eafaf1", "white")))
  })

  # ====================================================================
  # Tab 7: Biomarker Panel
  # ====================================================================

  output$bm_nephrin <- renderPlotly({
    req(sim_data())
    d <- sim_data()
    plot_ly() %>%
      add_lines(x = d$time, y = d$anti_nephrin_ab,
                name = "항-Nephrin 항체",
                line = list(color = "#8e44ad", width = 2.5)) %>%
      add_lines(x = d$time, y = d$u_nephrin,
                name = "뇨중 Nephrin (ng/mL)",
                yaxis = "y2",
                line = list(color = "#3498db", width = 2.5)) %>%
      layout(
        xaxis  = list(title = "시간 (일)"),
        yaxis  = list(title = "항-Nephrin 항체 (상대치)"),
        yaxis2 = list(title = "뇨중 Nephrin (ng/mL)", side = "right",
                      overlaying = "y"),
        legend = list(orientation = "h", y = -0.15),
        hovermode = "x unified"
      )
  })

  output$bm_complement_ige <- renderPlotly({
    req(sim_data())
    d <- sim_data()
    plot_ly() %>%
      add_lines(x = d$time, y = d$c3,
                name = "C3 (mg/dL)",
                line = list(color = "#27ae60", width = 2)) %>%
      add_lines(x = d$time, y = d$c4,
                name = "C4 (mg/dL)",
                line = list(color = "#2ecc71", width = 2, dash = "dash")) %>%
      add_lines(x = d$time, y = d$ige / 10,
                name = "IgE / 10 (IU/mL)",
                yaxis = "y2",
                line = list(color = "#e74c3c", width = 2)) %>%
      add_lines(x = d$time, y = d$eosinophil,
                name = "호산구 (%)",
                yaxis = "y2",
                line = list(color = "#f39c12", width = 2)) %>%
      layout(
        xaxis  = list(title = "시간 (일)"),
        yaxis  = list(title = "보체 (mg/dL)"),
        yaxis2 = list(title = "IgE/10, 호산구", side = "right",
                      overlaying = "y"),
        legend = list(orientation = "h", y = -0.2),
        hovermode = "x unified"
      )
  })

  output$bm_correlation_heatmap <- renderPlotly({
    req(sim_data())
    d     <- sim_data()
    bm_df <- data.frame(
      anti_nephrin = d$anti_nephrin_ab,
      suPAR        = d$supar,
      u_nephrin    = d$u_nephrin,
      C3           = d$c3,
      IgE          = d$ige,
      UPCR         = d$uprot_cr,
      Albumin      = d$serum_alb,
      Podocyte_Int = d$podocyte_integrity
    )
    corr_mat <- cor(bm_df, use = "pairwise.complete.obs")
    nm       <- colnames(corr_mat)

    plot_ly(
      z         = corr_mat,
      x         = nm,
      y         = nm,
      type      = "heatmap",
      colorscale = "RdBu",
      zmin       = -1, zmax = 1,
      text       = round(corr_mat, 2),
      texttemplate = "%{text}"
    ) %>%
      layout(
        margin = list(l = 110, b = 80),
        xaxis  = list(tickangle = -30)
      )
  })

  output$bm_decision_tree <- renderPlotly({
    req(sim_data())
    d <- sim_data()
    last_ab  <- tail(d$anti_nephrin_ab, 1)
    last_sup <- tail(d$supar, 1)
    last_up  <- tail(d$uprot_cr, 1)

    rec <- if (last_up < 0.3) {
      "완전 관해 달성\n현 치료 유지"
    } else if (last_ab > 1.5) {
      "고역가 항-Nephrin 항체\n→ Rituximab 추가 고려"
    } else if (last_sup > 3) {
      "고suPAR 수치\n→ CNI 전환 고려"
    } else {
      "부분 반응\n→ 용량 조정 또는 병용 요법"
    }

    plot_ly() %>%
      add_annotations(
        x = 0.5, y = 0.85,
        text = paste0("<b>항-Nephrin Ab: </b>", round(last_ab, 2),
                      "<br><b>suPAR: </b>", round(last_sup, 2),
                      "<br><b>UPCR: </b>", round(last_up, 2)),
        showarrow = FALSE, align = "center",
        font = list(size = 13)
      ) %>%
      add_annotations(
        x = 0.5, y = 0.45,
        text = paste0("<b>권고 사항:</b><br>", gsub("\n", "<br>", rec)),
        showarrow = FALSE, align = "center",
        font = list(size = 14, color = "#e74c3c")
      ) %>%
      layout(
        xaxis = list(visible = FALSE, range = c(0, 1)),
        yaxis = list(visible = FALSE, range = c(0, 1)),
        paper_bgcolor = "#f8f9fa",
        plot_bgcolor  = "#f8f9fa",
        shapes = list(
          list(type = "rect", x0 = 0.1, x1 = 0.9,
               y0 = 0.7, y1 = 1.0,
               fillcolor = "#d6eaf8", line = list(color = "#3498db")),
          list(type = "rect", x0 = 0.1, x1 = 0.9,
               y0 = 0.25, y1 = 0.60,
               fillcolor = "#fadbd8", line = list(color = "#e74c3c"))
        )
      )
  })

  output$bm_table <- DT::renderDT({
    req(sim_data())
    d     <- sim_data()
    tidx  <- round(seq(1, length(d$time), length.out = 19))
    df    <- data.frame(
      "시간 (일)"        = d$time[tidx],
      "항-Nephrin 항체"  = round(d$anti_nephrin_ab[tidx], 3),
      "suPAR (ng/mL)"   = round(d$supar[tidx], 2),
      "뇨중 Nephrin"     = round(d$u_nephrin[tidx], 3),
      "C3 (mg/dL)"      = round(d$c3[tidx], 1),
      "C4 (mg/dL)"      = round(d$c4[tidx], 1),
      "IgE (IU/mL)"     = round(d$ige[tidx], 0),
      "호산구 (%)"       = round(d$eosinophil[tidx], 1),
      check.names = FALSE
    )
    DT::datatable(df, options = list(pageLength = 10, dom = "tp"),
                  rownames = FALSE)
  })

  # ---- Auto-run simulation on startup ----
  observe({
    isolate({
      shinyjs::runjs("setTimeout(function(){
        document.getElementById('run_sim').click(); }, 500);")
    })
  })

} # end server

# ==============================================================================
# Launch
# ==============================================================================
shinyApp(ui = ui, server = server)
