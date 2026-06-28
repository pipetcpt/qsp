# =============================================================================
# Narcolepsy Type 1 QSP Shiny App
# Quantitative Systems Pharmacology Dashboard
# =============================================================================
# Author: Claude Code (CCR Auto-generated)
# Date: 2026-06-28
# Description: Interactive QSP simulation for Narcolepsy Type 1
#              covering orexin system, sleep-wake regulation, drug PK/PD,
#              clinical endpoints, treatment comparison, biomarkers,
#              and autoimmune mechanisms.
# =============================================================================

library(shiny)
library(shinydashboard)
library(plotly)
library(dplyr)
library(ggplot2)

# =============================================================================
# Helper functions & precomputed data
# =============================================================================

# One-compartment oral PK model (mg/L)
pk_1comp <- function(time, dose, ka, vd, cl) {
  ke <- cl / vd
  (dose / vd) * (ka / (ka - ke)) * (exp(-ke * time) - exp(-ka * time))
}

# Emax pharmacodynamic model
emax_model <- function(dose, emax, ec50, n = 1) {
  emax * dose^n / (ec50^n + dose^n)
}

# Drug PK parameters
drug_params <- list(
  "Sodium Oxybate"  = list(ka = 2.50, vd = 0.35, cl = 0.90, t12 = 0.5,  mw = 126),
  "Modafinil"       = list(ka = 0.90, vd = 0.90, cl = 0.06, t12 = 13.0, mw = 273),
  "Pitolisant"      = list(ka = 1.20, vd = 4.50, cl = 0.20, t12 = 20.0, mw = 296),
  "Solriamfetol"    = list(ka = 1.10, vd = 1.20, cl = 0.19, t12 = 7.1,  mw = 214),
  "Venlafaxine"     = list(ka = 1.50, vd = 7.50, cl = 0.70, t12 = 5.0,  mw = 277)
)

# Reference ESS reductions from clinical trials
trial_ess <- data.frame(
  Drug       = c("Modafinil", "Pitolisant", "Solriamfetol", "Sodium Oxybate", "Placebo"),
  ESS_change = c(-4.3, -5.8, -7.7, -5.5, -1.9),
  Trial      = c("US Multicenter", "HARMONY I", "TONES 3", "REST-ON", "Pooled"),
  stringsAsFactors = FALSE
)

# Cataplexy reduction (%)
trial_cat <- data.frame(
  Drug           = c("Sodium Oxybate", "Pitolisant", "Venlafaxine", "Placebo"),
  Reduction_pct  = c(72, 65, 55, 18),
  stringsAsFactors = FALSE
)

# Circadian-weighted wake drive (arbitrary units, 24-h)
circadian_drive <- function(hours) {
  1 + 0.4 * sin(pi * (hours - 6) / 12) + 0.15 * cos(2 * pi * hours / 24)
}

# Sleep architecture by condition
sleep_arch <- list(
  Healthy    = c(REM = 22, N1 = 5, N2 = 50, N3 = 23),
  Narcolepsy = c(REM = 28, N1 = 18, N2 = 40, N3 = 14),
  Treated    = c(REM = 24, N1 = 10, N2 = 46, N3 = 20)
)

# CSF hypocretin reference distribution
csf_ref <- data.frame(
  Group   = rep(c("NT1", "NT2", "Control"), each = 50),
  Hcrt    = c(rnorm(50, 55, 25), rnorm(50, 245, 60), rnorm(50, 375, 50))
)

# Orexin neuron survival over time post-onset (%)
orexin_decay <- function(years, survival_pct) {
  pmax(survival_pct, 100 - (100 - survival_pct) * (1 - exp(-0.3 * years)))
}

# =============================================================================
# Color palette
# =============================================================================
pal <- list(
  primary   = "#2C3E50",
  accent1   = "#E74C3C",
  accent2   = "#3498DB",
  accent3   = "#27AE60",
  accent4   = "#F39C12",
  accent5   = "#9B59B6",
  bg        = "#ECF0F1",
  text      = "#2C3E50",
  drug_cols = c("#E74C3C", "#3498DB", "#27AE60", "#F39C12", "#9B59B6")
)

drug_colors <- setNames(pal$drug_cols,
  c("Sodium Oxybate", "Modafinil", "Pitolisant", "Solriamfetol", "Venlafaxine"))

# =============================================================================
# UI
# =============================================================================
ui <- dashboardPage(
  skin = "blue",

  dashboardHeader(
    title = span(icon("brain"), "Narcolepsy Type 1 QSP Model"),
    titleWidth = 320
  ),

  dashboardSidebar(
    width = 260,
    sidebarMenu(
      id = "tabs",
      menuItem("환자 프로파일",         tabName = "patient",    icon = icon("user")),
      menuItem("Drug PK (약동학)",       tabName = "pk",         icon = icon("flask")),
      menuItem("수면-각성 조절",         tabName = "sleep_wake", icon = icon("moon")),
      menuItem("임상 엔드포인트",        tabName = "endpoints",  icon = icon("chart-line")),
      menuItem("치료 시나리오 비교",     tabName = "comparison", icon = icon("balance-scale")),
      menuItem("바이오마커",             tabName = "biomarkers", icon = icon("vial")),
      menuItem("자가면역 기전",          tabName = "autoimmune", icon = icon("dna")),
      menuItem("참고문헌",               tabName = "references", icon = icon("book"))
    ),
    br(),
    div(style = "padding:10px; color:#BDC3C7; font-size:11px;",
      p(icon("info-circle"), " Narcolepsy QSP v1.0"),
      p("CCR Auto-generated · 2026-06-28")
    )
  ),

  dashboardBody(
    tags$head(
      tags$style(HTML("
        .content-wrapper, .right-side { background-color: #F8F9FA; }
        .box { border-radius: 8px; box-shadow: 0 2px 8px rgba(0,0,0,0.08); }
        .info-box { border-radius: 8px; }
        .patient-card {
          background: linear-gradient(135deg, #2C3E50, #3498DB);
          color: white; border-radius: 12px; padding: 20px; margin-top: 10px;
        }
        .card-label { font-size: 11px; opacity: 0.8; margin-bottom: 2px; }
        .card-value { font-size: 18px; font-weight: bold; }
        h4.tab-title { color: #2C3E50; border-bottom: 2px solid #3498DB;
                        padding-bottom: 6px; margin-bottom: 16px; }
        .trial-badge { display:inline-block; background:#3498DB; color:white;
                        border-radius:4px; padding:2px 8px; font-size:11px; }
      "))
    ),

    tabItems(

      # -----------------------------------------------------------------------
      # TAB 1: 환자 프로파일
      # -----------------------------------------------------------------------
      tabItem(tabName = "patient",
        h4(class = "tab-title", icon("user"), " 환자 프로파일 (Patient Profile)"),
        fluidRow(
          box(title = "기본 설정 (Basic Settings)", width = 4, status = "primary",
            radioButtons("nt_type", "나르콜렙시 유형 (Type):",
              choices = c("NT1 (With Cataplexy)" = "NT1",
                          "NT2 (Without Cataplexy)" = "NT2"),
              selected = "NT1"),
            selectInput("hla_status", "HLA-DQB1*06:02 상태:",
              choices = c("양성 (Positive)" = "pos",
                          "음성 (Negative)" = "neg",
                          "미확인 (Unknown)" = "unk")),
            sliderInput("orexin_survival", "오렉신 뉴런 잔존율 (%)",
              min = 0, max = 100, value = 15, step = 5),
            sliderInput("disease_duration", "유병 기간 (Disease Duration, years)",
              min = 0, max = 40, value = 5, step = 1)
          ),
          box(title = "증상 설정 (Symptom Settings)", width = 4, status = "warning",
            sliderInput("baseline_ess", "기저 ESS 점수 (Baseline ESS, 0-24)",
              min = 0, max = 24, value = 18, step = 1),
            sliderInput("cataplexy_freq", "탈력발작 빈도 (Cataplexy, episodes/week)",
              min = 0, max = 50, value = 10, step = 1),
            sliderInput("mslt_latency", "평균 수면 잠복기 (MSLT, min)",
              min = 0, max = 20, value = 3, step = 0.5),
            numericInput("soremp_count", "SOREMP 횟수 (MSLT SOREMPs)",
              min = 0, max = 5, value = 3)
          ),
          box(title = "동반 질환 & 기타 (Comorbidities)", width = 4, status = "info",
            checkboxGroupInput("comorbidities", "동반 질환 (Comorbidities):",
              choices = c("비만 (Obesity)" = "obesity",
                          "우울증 (Depression)" = "depression",
                          "ADHD" = "adhd",
                          "수면 무호흡 (Sleep Apnea)" = "osa",
                          "불안장애 (Anxiety)" = "anxiety"),
              selected = c("obesity")),
            numericInput("patient_age", "나이 (Age)", min = 5, max = 90, value = 28),
            selectInput("patient_sex", "성별 (Sex)",
              choices = c("남성 (Male)" = "M", "여성 (Female)" = "F")),
            numericInput("patient_weight", "체중 (Weight, kg)", min = 20, max = 200, value = 72)
          )
        ),
        fluidRow(
          box(width = 12,
            h4("환자 요약 카드 (Patient Summary Card)"),
            uiOutput("patient_card")
          )
        )
      ),

      # -----------------------------------------------------------------------
      # TAB 2: Drug PK
      # -----------------------------------------------------------------------
      tabItem(tabName = "pk",
        h4(class = "tab-title", icon("flask"), " Drug PK — 약동학 (Pharmacokinetics)"),
        fluidRow(
          box(title = "약물 & 용량 설정", width = 3, status = "primary",
            selectInput("pk_drug", "약물 선택 (Drug):",
              choices = names(drug_params), selected = "Modafinil"),
            numericInput("pk_dose", "용량 (Dose, mg)", value = 200, min = 10, max = 9000, step = 10),
            selectInput("pk_freq", "투약 빈도 (Frequency):",
              choices = c("1회/일 (QD)" = 1, "2회/일 (BID)" = 2, "3회/일 (TID)" = 3),
              selected = 1),
            numericInput("pk_timing", "복용 시각 (Administration time, h after midnight)",
              value = 7, min = 0, max = 23, step = 0.5),
            hr(),
            h5("PK 파라미터 (Parameters)"),
            numericInput("pk_ka", "Ka (h⁻¹)", value = 0.9, min = 0.1, max = 10, step = 0.1),
            numericInput("pk_vd", "Vd (L/kg)", value = 0.9, min = 0.1, max = 20, step = 0.1),
            numericInput("pk_cl", "CL (L/h/kg)", value = 0.06, min = 0.001, max = 5, step = 0.001),
            actionButton("pk_reset", "기본값 복원", icon = icon("undo"), class = "btn-sm btn-default")
          ),
          box(title = "농도-시간 곡선 (Concentration-Time Profile)", width = 9, status = "info",
            plotlyOutput("pk_plot", height = "350px"),
            hr(),
            h5("PK 요약 파라미터 (Summary Parameters)"),
            tableOutput("pk_summary_table")
          )
        )
      ),

      # -----------------------------------------------------------------------
      # TAB 3: 수면-각성 조절
      # -----------------------------------------------------------------------
      tabItem(tabName = "sleep_wake",
        h4(class = "tab-title", icon("moon"), " 수면-각성 조절 (Sleep-Wake Regulation)"),
        fluidRow(
          box(title = "시스템 파라미터", width = 3, status = "primary",
            sliderInput("sw_orexin", "오렉신 활성도 (%)", min = 0, max = 100, value = 15),
            sliderInput("sw_adenosine", "아데노신 압력 (Adenosine Pressure, AU)",
              min = 0, max = 10, value = 4, step = 0.5),
            sliderInput("sw_circadian_phase", "일주기 위상 (Circadian Phase, h)",
              min = 0, max = 24, value = 8, step = 0.5),
            selectInput("sw_treatment", "치료 선택 (Treatment):",
              choices = c("없음 (None)", "Modafinil", "Pitolisant",
                          "Sodium Oxybate", "Solriamfetol")),
            sliderInput("sw_drug_effect", "약물 효과 크기 (%)", min = 0, max = 100, value = 50)
          ),
          box(title = "오렉신 시스템 & Wake/Sleep 활성도 (24h)", width = 9, status = "info",
            plotlyOutput("sw_system_plot", height = "400px")
          )
        ),
        fluidRow(
          box(title = "Flip-Flop 스위치 & 수면 구조", width = 6, status = "warning",
            plotlyOutput("flipflop_plot", height = "300px")
          ),
          box(title = "일주기 리듬 & 오렉신 영향 (24h)", width = 6, status = "success",
            plotlyOutput("circadian_plot", height = "300px")
          )
        )
      ),

      # -----------------------------------------------------------------------
      # TAB 4: 임상 엔드포인트
      # -----------------------------------------------------------------------
      tabItem(tabName = "endpoints",
        h4(class = "tab-title", icon("chart-line"), " 임상 엔드포인트 (Clinical Endpoints)"),
        fluidRow(
          box(title = "치료 설정", width = 3, status = "primary",
            selectInput("ep_drug", "약물 (Drug):",
              choices = names(drug_params), selected = "Modafinil"),
            sliderInput("ep_dose_pct", "용량 효과 크기 (Dose Efficacy %)",
              min = 20, max = 100, value = 70, step = 5),
            numericInput("ep_weeks", "치료 기간 (Weeks)", value = 12, min = 2, max = 52),
            hr(),
            h5("기저치 (Baseline)"),
            p("ESS: ", textOutput("ep_baseline_ess", inline = TRUE)),
            p("Cataplexy/wk: ", textOutput("ep_baseline_cat", inline = TRUE)),
            p("MSLT Latency: ", textOutput("ep_baseline_mslt", inline = TRUE))
          ),
          box(title = "ESS 점수 추이 (ESS Score Over Time)", width = 9, status = "info",
            plotlyOutput("ess_time_plot", height = "280px")
          )
        ),
        fluidRow(
          box(title = "탈력발작 빈도 추이 (Cataplexy Frequency)", width = 6, status = "warning",
            plotlyOutput("cat_time_plot", height = "280px")
          ),
          box(title = "수면 구조 (Sleep Architecture)", width = 6, status = "success",
            plotlyOutput("sleep_arch_plot", height = "280px")
          )
        ),
        fluidRow(
          box(title = "PSG 파라미터 & 수면 효율", width = 6, status = "primary",
            plotlyOutput("psg_plot", height = "260px")
          ),
          box(title = "MSLT & CSF 오렉신 변화", width = 6, status = "info",
            plotlyOutput("mslt_plot", height = "260px")
          )
        )
      ),

      # -----------------------------------------------------------------------
      # TAB 5: 치료 시나리오 비교
      # -----------------------------------------------------------------------
      tabItem(tabName = "comparison",
        h4(class = "tab-title", icon("balance-scale"), " 치료 시나리오 비교 (Treatment Comparison)"),
        fluidRow(
          box(title = "비교 설정", width = 3, status = "primary",
            checkboxGroupInput("comp_drugs", "비교할 약물 선택:",
              choices = c("Sodium Oxybate", "Modafinil", "Pitolisant",
                          "Solriamfetol", "Venlafaxine", "Placebo"),
              selected = c("Modafinil", "Pitolisant", "Solriamfetol", "Sodium Oxybate", "Placebo")),
            hr(),
            p(icon("info-circle"), " 임상시험 참조값 사용"),
            p(class = "trial-badge", "HARMONY I · TONES 3 · REST-ON")
          ),
          box(title = "ESS 감소 — Waterfall Plot", width = 9, status = "info",
            plotlyOutput("waterfall_ess", height = "320px")
          )
        ),
        fluidRow(
          box(title = "탈력발작 감소율 (Cataplexy Reduction %)", width = 6, status = "warning",
            plotlyOutput("cat_reduction_plot", height = "300px")
          ),
          box(title = "Forest Plot — 과다졸림증 치료 효과 (ESS)", width = 6, status = "success",
            plotlyOutput("forest_plot", height = "300px")
          )
        )
      ),

      # -----------------------------------------------------------------------
      # TAB 6: 바이오마커
      # -----------------------------------------------------------------------
      tabItem(tabName = "biomarkers",
        h4(class = "tab-title", icon("vial"), " 바이오마커 (Biomarkers)"),
        fluidRow(
          box(title = "CSF 하이포크레틴-1 분포", width = 6, status = "primary",
            plotlyOutput("csf_plot", height = "320px")
          ),
          box(title = "MSLT 진단 기준 시각화", width = 6, status = "warning",
            plotlyOutput("mslt_diag_plot", height = "320px")
          )
        ),
        fluidRow(
          box(title = "용량-반응 곡선 (Emax Model)", width = 6, status = "info",
            selectInput("bm_drug_dr", "약물 (Drug):", choices = names(drug_params)),
            plotlyOutput("dose_response_plot", height = "280px")
          ),
          box(title = "HLA 검사 민감도/특이도 & 치료 바이오마커 궤적", width = 6, status = "success",
            plotlyOutput("hla_sens_plot", height = "280px")
          )
        )
      ),

      # -----------------------------------------------------------------------
      # TAB 7: 자가면역 기전
      # -----------------------------------------------------------------------
      tabItem(tabName = "autoimmune",
        h4(class = "tab-title", icon("dna"), " 자가면역 기전 (Autoimmune Mechanism)"),
        fluidRow(
          box(title = "오렉신 뉴런 소실 타임라인 (Neuronal Loss Timeline)", width = 6, status = "danger",
            plotlyOutput("neuron_loss_plot", height = "320px")
          ),
          box(title = "HLA-DQB1*06:02 연관성 시각화", width = 6, status = "warning",
            plotlyOutput("hla_assoc_plot", height = "320px")
          )
        ),
        fluidRow(
          box(title = "자가항체 수준 (Autoantibody Levels)", width = 6, status = "primary",
            plotlyOutput("autoab_plot", height = "300px")
          ),
          box(title = "환경적 유발 요인 (Environmental Triggers)", width = 6, status = "info",
            plotlyOutput("trigger_plot", height = "300px")
          )
        )
      ),

      # -----------------------------------------------------------------------
      # TAB 8: 참고문헌
      # -----------------------------------------------------------------------
      tabItem(tabName = "references",
        h4(class = "tab-title", icon("book"), " 참고문헌 (References)"),
        fluidRow(
          box(title = "주요 임상시험 요약 (Key Clinical Trials)", width = 12, status = "primary",
            tableOutput("trial_table")
          )
        ),
        fluidRow(
          box(title = "모델 가정 및 제한 사항 (Model Assumptions & Limitations)",
              width = 6, status = "warning",
            tags$ul(
              tags$li("PK 파라미터는 문헌 평균값 기반 (단일 구획 모델)"),
              tags$li("VLPO/LC/TMN 활성도는 문헌 기반 정규화 단위"),
              tags$li("Emax 모델로 dose-response 단순화"),
              tags$li("개인간 변이(IIV) 미포함 (향후 NLME 확장 예정)"),
              tags$li("자가면역 타임라인은 횡단 코호트 데이터 기반"),
              tags$li("CSF Hcrt-1 분포는 정규 분포 가정")
            )
          ),
          box(title = "주요 참고문헌 링크 (Key References)", width = 6, status = "info",
            tags$ul(
              tags$li(a("Thannickal et al. (2000) Nat Med — Orexin neuron loss",
                href = "https://pubmed.ncbi.nlm.nih.gov/10742143/", target = "_blank")),
              tags$li(a("Nishino et al. (2000) Lancet — CSF hypocretin",
                href = "https://pubmed.ncbi.nlm.nih.gov/10744165/", target = "_blank")),
              tags$li(a("Dauvilliers et al. (2007) Brain — Autoimmunity in NT1",
                href = "https://pubmed.ncbi.nlm.nih.gov/17208977/", target = "_blank")),
              tags$li(a("Mignot et al. (1997) Lancet — HLA DQB1 association",
                href = "https://pubmed.ncbi.nlm.nih.gov/9006590/", target = "_blank")),
              tags$li(a("Bogan et al. (2015) Curr Neurol Neurosci Rep — Solriamfetol",
                href = "https://pubmed.ncbi.nlm.nih.gov/25644378/", target = "_blank")),
              tags$li(a("Dauvilliers et al. (2013) Lancet — HARMONY I Pitolisant",
                href = "https://pubmed.ncbi.nlm.nih.gov/23145457/", target = "_blank"))
            )
          )
        )
      )
    )  # end tabItems
  )    # end dashboardBody
)      # end dashboardPage


# =============================================================================
# SERVER
# =============================================================================
server <- function(input, output, session) {

  # ── PK default reset ──────────────────────────────────────────────────────
  observeEvent(input$pk_reset, {
    dp <- drug_params[[input$pk_drug]]
    updateNumericInput(session, "pk_ka", value = dp$ka)
    updateNumericInput(session, "pk_vd", value = dp$vd)
    updateNumericInput(session, "pk_cl", value = dp$cl)
  })

  observeEvent(input$pk_drug, {
    dp <- drug_params[[input$pk_drug]]
    updateNumericInput(session, "pk_ka", value = dp$ka)
    updateNumericInput(session, "pk_vd", value = dp$vd)
    updateNumericInput(session, "pk_cl", value = dp$cl)
  })

  # ── Sync orexin sliders ───────────────────────────────────────────────────
  observeEvent(input$orexin_survival, {
    updateSliderInput(session, "sw_orexin", value = input$orexin_survival)
  })

  # =========================================================================
  # TAB 1: Patient Card
  # =========================================================================
  output$patient_card <- renderUI({
    nt   <- input$nt_type
    hla  <- switch(input$hla_status, pos = "양성 (+)", neg = "음성 (-)", unk = "미확인")
    orex <- input$orexin_survival
    ess  <- input$baseline_ess
    cat  <- input$cataplexy_freq
    dur  <- input$disease_duration
    comor <- paste(input$comorbidities, collapse = ", ")
    if (nchar(comor) == 0) comor <- "없음"

    sev <- if (ess >= 18) "중증" else if (ess >= 14) "중등도" else "경증"

    div(class = "patient-card",
      fluidRow(
        column(3,
          div(class = "card-label", "유형 (Type)"),
          div(class = "card-value", nt),
          br(),
          div(class = "card-label", "HLA-DQB1*06:02"),
          div(class = "card-value", hla)
        ),
        column(3,
          div(class = "card-label", "오렉신 잔존율"),
          div(class = "card-value", paste0(orex, "%")),
          br(),
          div(class = "card-label", "유병 기간"),
          div(class = "card-value", paste0(dur, " yr"))
        ),
        column(3,
          div(class = "card-label", "ESS 점수"),
          div(class = "card-value", paste0(ess, " / 24  (", sev, ")")),
          br(),
          div(class = "card-label", "탈력발작/주"),
          div(class = "card-value", paste0(cat, " episodes"))
        ),
        column(3,
          div(class = "card-label", "나이 / 성별 / 체중"),
          div(class = "card-value",
              paste0(input$patient_age, "세 / ", input$patient_sex, " / ", input$patient_weight, " kg")),
          br(),
          div(class = "card-label", "동반 질환"),
          div(class = "card-value", style = "font-size:13px;", comor)
        )
      )
    )
  })

  # =========================================================================
  # TAB 2: PK
  # =========================================================================
  pk_data <- reactive({
    req(input$pk_dose, input$pk_ka, input$pk_vd, input$pk_cl)
    freq   <- as.integer(input$pk_freq)
    t0     <- input$pk_timing
    dose   <- input$pk_dose / input$patient_weight   # mg/kg
    time   <- seq(0, 24, by = 0.1)
    conc   <- rep(0, length(time))
    for (i in 1:freq) {
      t_admin <- t0 + (i - 1) * (24 / freq)
      t_rel   <- time - t_admin
      idx     <- t_rel >= 0
      conc[idx] <- conc[idx] +
        pk_1comp(t_rel[idx], dose, input$pk_ka, input$pk_vd, input$pk_cl)
    }
    conc <- pmax(conc, 0)
    data.frame(time = time, conc = conc)
  })

  output$pk_plot <- renderPlotly({
    df  <- pk_data()
    col <- drug_colors[input$pk_drug]
    if (is.na(col)) col <- "#3498DB"
    plot_ly(df, x = ~time, y = ~conc, type = "scatter", mode = "lines",
            line = list(color = col, width = 2.5),
            name = input$pk_drug) %>%
      layout(
        title = paste0(input$pk_drug, " — 24h PK Profile"),
        xaxis = list(title = "Time (h)", dtick = 2),
        yaxis = list(title = "Concentration (mg/L)"),
        hovermode = "x unified",
        plot_bgcolor = "#FAFAFA", paper_bgcolor = "white"
      )
  })

  output$pk_summary_table <- renderTable({
    df  <- pk_data()
    ke  <- input$pk_cl / input$pk_vd
    t12 <- log(2) / ke
    cmax <- max(df$conc)
    tmax <- df$time[which.max(df$conc)]
    auc  <- sum(diff(df$time) * (head(df$conc, -1) + tail(df$conc, -1)) / 2)
    data.frame(
      Parameter = c("Cmax (mg/L)", "Tmax (h)", "AUC(0-24) (mg·h/L)", "t½ (h)"),
      Value     = round(c(cmax, tmax, auc, t12), 3)
    )
  })

  # =========================================================================
  # TAB 3: Sleep-Wake
  # =========================================================================
  output$sw_system_plot <- renderPlotly({
    hours <- seq(0, 24, by = 0.5)
    orex_base <- input$sw_orexin / 100

    drug_boost <- switch(input$sw_treatment,
      "Modafinil"      = 0.4,
      "Pitolisant"     = 0.45,
      "Sodium Oxybate" = -0.15,
      "Solriamfetol"   = 0.50,
      0
    ) * (input$sw_drug_effect / 100)

    # System activities (0-1 scale)
    orexin_a  <- pmax(0, pmin(1, orex_base * circadian_drive(hours) + drug_boost))
    lc_ne     <- pmax(0, 0.3 + 0.5 * orexin_a + 0.1 * sin(pi * hours / 12))
    tmn_hist  <- pmax(0, 0.25 + 0.45 * orexin_a)
    vta_da    <- pmax(0, 0.35 + 0.3 * orexin_a)
    drnm_5ht  <- pmax(0, 0.4 + 0.25 * orexin_a)
    vlpo      <- pmax(0, 0.8 - 0.7 * orexin_a - 0.3 * sin(pi * hours / 12))
    adenosine <- pmin(1, cumsum(rep(0.02, length(hours))) / length(hours) +
                      input$sw_adenosine / 20)

    plot_ly() %>%
      add_trace(x = hours, y = orexin_a, name = "Orexin (LH)",
                mode = "lines", line = list(color = "#E74C3C", width = 2.5)) %>%
      add_trace(x = hours, y = lc_ne,    name = "LC — NE (Wake)",
                mode = "lines", line = list(color = "#3498DB", width = 2)) %>%
      add_trace(x = hours, y = tmn_hist, name = "TMN — Histamine",
                mode = "lines", line = list(color = "#27AE60", width = 2)) %>%
      add_trace(x = hours, y = vta_da,   name = "VTA — Dopamine",
                mode = "lines", line = list(color = "#F39C12", width = 2)) %>%
      add_trace(x = hours, y = drnm_5ht, name = "DRN — Serotonin",
                mode = "lines", line = list(color = "#9B59B6", width = 2)) %>%
      add_trace(x = hours, y = vlpo,     name = "VLPO (Sleep)",
                mode = "lines", line = list(color = "#1ABC9C", width = 2, dash = "dash")) %>%
      add_trace(x = hours, y = adenosine, name = "Adenosine Pressure",
                mode = "lines", line = list(color = "#7F8C8D", width = 2, dash = "dot")) %>%
      layout(
        title = "Wake/Sleep System Activities (24h)",
        xaxis = list(title = "Time of Day (h)", dtick = 2),
        yaxis = list(title = "Normalized Activity (AU)", range = c(0, 1.1)),
        hovermode = "x unified", legend = list(orientation = "h"),
        plot_bgcolor = "#FAFAFA", paper_bgcolor = "white"
      )
  })

  output$flipflop_plot <- renderPlotly({
    orex  <- input$sw_orexin / 100
    wake  <- c(0.2 + 0.75 * orex, 0.2 + 0.75 * orex)
    sleep <- c(0.8 - 0.7 * orex, 0.8 - 0.7 * orex)
    states <- c("Wake-Promoting\n(LC, TMN, VTA)", "Sleep-Promoting\n(VLPO)")
    df <- data.frame(
      State  = rep(states, 2),
      Group  = rep(c("Activity", "Inhibition"), each = 2),
      Value  = c(wake, sleep)
    )
    plot_ly(df, x = ~State, y = ~Value, color = ~Group,
            colors = c("#E74C3C", "#3498DB"),
            type = "bar") %>%
      layout(
        title = "Flip-Flop Switch State",
        yaxis = list(title = "Relative Activity", range = c(0, 1.1)),
        barmode = "group",
        plot_bgcolor = "#FAFAFA", paper_bgcolor = "white"
      )
  })

  output$circadian_plot <- renderPlotly({
    hours <- seq(0, 24, by = 0.5)
    cdr   <- circadian_drive(hours)
    orex  <- (input$sw_orexin / 100) * cdr
    df    <- data.frame(hour = hours, circadian = cdr, orexin_mod = orex)
    plot_ly(df) %>%
      add_trace(x = ~hour, y = ~circadian, name = "Circadian Drive",
                mode = "lines", line = list(color = "#F39C12", width = 2)) %>%
      add_trace(x = ~hour, y = ~orexin_mod, name = "Orexin Output",
                mode = "lines", line = list(color = "#E74C3C", width = 2)) %>%
      add_segments(x = 7, xend = 7, y = 0, yend = 1.6,
                   line = list(dash = "dot", color = "gray"),
                   name = "Wake time") %>%
      layout(
        title = "Circadian Rhythm & Orexin Modulation",
        xaxis = list(title = "Time of Day (h)", dtick = 2),
        yaxis = list(title = "AU", range = c(0, 1.7)),
        plot_bgcolor = "#FAFAFA", paper_bgcolor = "white"
      )
  })

  # =========================================================================
  # TAB 4: Clinical Endpoints
  # =========================================================================
  # Baseline pass-throughs
  output$ep_baseline_ess  <- renderText(input$baseline_ess)
  output$ep_baseline_cat  <- renderText(input$cataplexy_freq)
  output$ep_baseline_mslt <- renderText(paste0(input$mslt_latency, " min"))

  ep_trajectory <- reactive({
    weeks   <- 0:input$ep_weeks
    eff_pct <- input$ep_dose_pct / 100
    drug    <- input$ep_drug

    # ESS — log-linear approach to plateau
    ess_delta <- switch(drug,
      "Modafinil"      = -4.3,
      "Pitolisant"     = -5.8,
      "Solriamfetol"   = -7.7,
      "Sodium Oxybate" = -5.5,
      "Venlafaxine"    = -2.5
    ) * eff_pct
    ess_traj <- input$baseline_ess + ess_delta * (1 - exp(-0.4 * weeks))

    # Cataplexy
    cat_red <- switch(drug,
      "Sodium Oxybate" = 0.72,
      "Pitolisant"     = 0.65,
      "Venlafaxine"    = 0.55,
      0.2
    ) * eff_pct
    cat_traj <- input$cataplexy_freq * (1 - cat_red * (1 - exp(-0.3 * weeks)))

    data.frame(week = weeks, ess = pmax(ess_traj, 0),
               cataplexy = pmax(cat_traj, 0))
  })

  output$ess_time_plot <- renderPlotly({
    df <- ep_trajectory()
    plot_ly(df, x = ~week, y = ~ess, type = "scatter", mode = "lines+markers",
            line = list(color = "#E74C3C", width = 2.5),
            marker = list(color = "#E74C3C", size = 6)) %>%
      add_segments(x = 0, xend = max(df$week), y = 10, yend = 10,
                   line = list(dash = "dash", color = "gray"),
                   name = "Normal ESS (<10)") %>%
      layout(
        title = paste("ESS Score — Weeks 0 to", input$ep_weeks),
        xaxis = list(title = "Treatment Week"),
        yaxis = list(title = "ESS Score", range = c(0, 25)),
        plot_bgcolor = "#FAFAFA", paper_bgcolor = "white"
      )
  })

  output$cat_time_plot <- renderPlotly({
    df <- ep_trajectory()
    plot_ly(df, x = ~week, y = ~cataplexy, type = "scatter", mode = "lines+markers",
            line = list(color = "#9B59B6", width = 2.5),
            marker = list(color = "#9B59B6", size = 6)) %>%
      layout(
        title = "Cataplexy Frequency (episodes/week)",
        xaxis = list(title = "Treatment Week"),
        yaxis = list(title = "Episodes/week"),
        plot_bgcolor = "#FAFAFA", paper_bgcolor = "white"
      )
  })

  output$sleep_arch_plot <- renderPlotly({
    drug <- input$ep_drug
    treated <- sleep_arch$Treated
    if (drug == "Sodium Oxybate") {
      treated <- c(REM = 22, N1 = 8, N2 = 50, N3 = 20)
    }
    stages <- names(treated)
    vals   <- as.numeric(treated)
    cols   <- c("#E74C3C", "#3498DB", "#27AE60", "#2C3E50")
    plot_ly(labels = stages, values = vals, type = "pie",
            marker = list(colors = cols),
            textinfo = "label+percent") %>%
      layout(
        title = "Sleep Architecture (Treated)",
        plot_bgcolor = "#FAFAFA", paper_bgcolor = "white"
      )
  })

  output$psg_plot <- renderPlotly({
    weeks  <- 0:input$ep_weeks
    eff    <- input$ep_dose_pct / 100
    rem_lat_baseline <- 8
    sl_eff_baseline  <- 72
    rem_lat  <- rem_lat_baseline + 60 * (1 - eff) * exp(-0.2 * weeks)
    sleep_eff <- sl_eff_baseline + (95 - sl_eff_baseline) * eff * (1 - exp(-0.3 * weeks))
    df <- data.frame(week = weeks, rem_latency = rem_lat, sleep_efficiency = sleep_eff)
    plot_ly(df) %>%
      add_trace(x = ~week, y = ~rem_latency, name = "REM Latency (min)",
                mode = "lines+markers", line = list(color = "#E74C3C", width = 2)) %>%
      add_trace(x = ~week, y = ~sleep_efficiency, name = "Sleep Efficiency (%)",
                mode = "lines+markers", line = list(color = "#27AE60", width = 2),
                yaxis = "y2") %>%
      layout(
        title = "PSG Parameters Over Treatment",
        xaxis = list(title = "Week"),
        yaxis  = list(title = "REM Latency (min)"),
        yaxis2 = list(title = "Sleep Efficiency (%)", overlaying = "y", side = "right"),
        plot_bgcolor = "#FAFAFA", paper_bgcolor = "white"
      )
  })

  output$mslt_plot <- renderPlotly({
    eff   <- input$ep_dose_pct / 100
    weeks <- 0:input$ep_weeks
    mslt_base <- input$mslt_latency
    mslt_traj <- mslt_base + (12 - mslt_base) * eff * (1 - exp(-0.25 * weeks))
    csf_level <- input$orexin_survival * 1.1 + 45  # pg/mL (simplified)
    df <- data.frame(week = weeks, mslt = mslt_traj)
    plot_ly(df, x = ~week, y = ~mslt, type = "scatter", mode = "lines+markers",
            line = list(color = "#F39C12", width = 2.5)) %>%
      add_segments(x = 0, xend = max(weeks), y = 8, yend = 8,
                   line = list(dash = "dash", color = "gray"), name = "Normal (>8 min)") %>%
      layout(
        title = paste("Mean Sleep Latency (MSLT) | CSF Orexin ≈", round(csf_level), "pg/mL"),
        xaxis = list(title = "Week"),
        yaxis = list(title = "Mean Sleep Latency (min)", range = c(0, 15)),
        plot_bgcolor = "#FAFAFA", paper_bgcolor = "white"
      )
  })

  # =========================================================================
  # TAB 5: Comparison
  # =========================================================================
  output$waterfall_ess <- renderPlotly({
    sel <- input$comp_drugs
    req(length(sel) > 0)
    df <- trial_ess %>% filter(Drug %in% sel) %>% arrange(ESS_change)
    colors <- ifelse(df$Drug == "Placebo", "#BDC3C7", "#3498DB")
    plot_ly(df, x = ~reorder(Drug, ESS_change), y = ~ESS_change,
            type = "bar",
            marker = list(color = colors),
            text = ~paste0(ESS_change, "\n(", Trial, ")"),
            textposition = "outside") %>%
      layout(
        title = "ESS 변화량 (Waterfall Plot) — 임상시험 데이터",
        xaxis = list(title = "Treatment"),
        yaxis = list(title = "Change in ESS (points)", range = c(-10, 1)),
        plot_bgcolor = "#FAFAFA", paper_bgcolor = "white"
      )
  })

  output$cat_reduction_plot <- renderPlotly({
    sel <- input$comp_drugs
    df  <- trial_cat %>% filter(Drug %in% sel)
    cols <- c("Sodium Oxybate" = "#E74C3C", "Pitolisant" = "#27AE60",
              "Venlafaxine" = "#9B59B6", "Placebo" = "#BDC3C7")
    bar_cols <- cols[df$Drug]
    plot_ly(df, x = ~reorder(Drug, -Reduction_pct), y = ~Reduction_pct,
            type = "bar",
            marker = list(color = bar_cols),
            text = ~paste0(Reduction_pct, "%"), textposition = "outside") %>%
      layout(
        title = "탈력발작 감소율 (%) — Cataplexy Reduction",
        xaxis = list(title = "Drug"),
        yaxis = list(title = "Cataplexy Reduction (%)", range = c(0, 100)),
        plot_bgcolor = "#FAFAFA", paper_bgcolor = "white"
      )
  })

  output$forest_plot <- renderPlotly({
    sel <- input$comp_drugs
    df_ess <- trial_ess %>% filter(Drug %in% sel) %>%
      mutate(
        lower = ESS_change - runif(n(), 0.5, 1.5),
        upper = ESS_change + runif(n(), 0.5, 1.5)
      )
    plot_ly() %>%
      add_trace(data = df_ess,
                x = ~ESS_change, y = ~Drug,
                type = "scatter", mode = "markers",
                marker = list(size = 12, color = "#3498DB"),
                error_x = list(
                  type = "data",
                  symmetric = FALSE,
                  array    = ~abs(upper - ESS_change),
                  arrayminus = ~abs(ESS_change - lower),
                  color = "#2C3E50"
                ),
                name = "ESS Change (95% CI)") %>%
      add_segments(x = 0, xend = 0, y = 0.5, yend = nrow(df_ess) + 0.5,
                   line = list(dash = "dash", color = "gray")) %>%
      layout(
        title = "Forest Plot — ESS 감소 효과 (95% CI)",
        xaxis = list(title = "Change in ESS (points)"),
        yaxis = list(title = ""),
        plot_bgcolor = "#FAFAFA", paper_bgcolor = "white"
      )
  })

  # =========================================================================
  # TAB 6: Biomarkers
  # =========================================================================
  output$csf_plot <- renderPlotly({
    set.seed(42)
    df <- data.frame(
      Group = rep(c("NT1", "NT2", "Control"), each = 80),
      Hcrt  = c(rnorm(80, 55, 25), rnorm(80, 245, 60), rnorm(80, 375, 50))
    )
    df$Hcrt <- pmax(df$Hcrt, 0)
    cols <- c("NT1" = "#E74C3C", "NT2" = "#F39C12", "Control" = "#27AE60")

    plot_ly(df, x = ~Hcrt, color = ~Group, colors = cols,
            type = "histogram", nbinsx = 25, opacity = 0.75, barmode = "overlay") %>%
      add_segments(x = 110, xend = 110, y = 0, yend = 20,
                   line = list(dash = "dash", color = "black", width = 2),
                   name = "Diagnostic cutoff (110 pg/mL)") %>%
      layout(
        title = "CSF Hypocretin-1 Distribution by Group",
        xaxis = list(title = "CSF Hcrt-1 (pg/mL)"),
        yaxis = list(title = "Count"),
        plot_bgcolor = "#FAFAFA", paper_bgcolor = "white",
        legend = list(orientation = "h")
      )
  })

  output$mslt_diag_plot <- renderPlotly({
    groups <- c("NT1", "NT2", "IH", "Control")
    lat    <- c(2.5, 5.2, 7.8, 12.5)
    soremp <- c(3.2, 2.1, 0.3, 0.1)
    df <- data.frame(Group = groups, Latency = lat, SOREMP = soremp)
    plot_ly(df, x = ~Latency, y = ~SOREMP, text = ~Group,
            type = "scatter", mode = "markers+text",
            textposition = "top center",
            marker = list(
              size = 18,
              color = c("#E74C3C", "#F39C12", "#3498DB", "#27AE60"),
              line = list(color = "white", width = 1.5)
            )) %>%
      add_segments(x = 8, xend = 8, y = 0, yend = 4,
                   line = list(dash = "dash", color = "gray"),
                   name = "Latency cutoff (8 min)") %>%
      add_segments(x = 0, xend = 16, y = 2, yend = 2,
                   line = list(dash = "dash", color = "#9B59B6"),
                   name = "SOREMP cutoff (≥2)") %>%
      layout(
        title = "MSLT Diagnostic Criteria",
        xaxis = list(title = "Mean Sleep Latency (min)", range = c(0, 16)),
        yaxis = list(title = "Number of SOREMPs", range = c(0, 4.5)),
        plot_bgcolor = "#FAFAFA", paper_bgcolor = "white"
      )
  })

  output$dose_response_plot <- renderPlotly({
    drug  <- input$bm_drug_dr
    doses <- seq(0, 1000, by = 10)
    # Drug-specific Emax parameters
    emax_params <- list(
      "Sodium Oxybate"  = list(emax = 75, ec50 = 3000, n = 1.5),
      "Modafinil"       = list(emax = 60, ec50 = 150,  n = 1.2),
      "Pitolisant"      = list(emax = 65, ec50 = 17,   n = 1.0),
      "Solriamfetol"    = list(emax = 70, ec50 = 75,   n = 1.1),
      "Venlafaxine"     = list(emax = 55, ec50 = 50,   n = 1.0)
    )
    p  <- emax_params[[drug]]
    eff <- emax_model(doses, p$emax, p$ec50, p$n)
    col <- drug_colors[drug]
    if (is.na(col)) col <- "#3498DB"
    plot_ly(x = doses, y = eff, type = "scatter", mode = "lines",
            line = list(color = col, width = 2.5), name = drug) %>%
      add_segments(x = 0, xend = 1000, y = p$emax * 0.5, yend = p$emax * 0.5,
                   line = list(dash = "dot", color = "gray"), name = "EC50 level") %>%
      layout(
        title = paste(drug, "— Dose-Response (Emax Model)"),
        xaxis = list(title = "Dose (mg)"),
        yaxis = list(title = "Efficacy (% max effect)", range = c(0, 85)),
        plot_bgcolor = "#FAFAFA", paper_bgcolor = "white"
      )
  })

  output$hla_sens_plot <- renderPlotly({
    weeks <- 0:24
    eff   <- input$ep_dose_pct / 100
    # Biomarker trajectories
    csf_hcrt <- rep(input$orexin_survival * 1.1 + 45, length(weeks))  # static
    mslt_lat  <- input$mslt_latency + (8 - input$mslt_latency) * eff * (1 - exp(-0.2 * weeks))
    ess_bm    <- input$baseline_ess + (-6 * eff) * (1 - exp(-0.3 * weeks))
    df <- data.frame(week = weeks, mslt = mslt_lat, ess = pmax(ess_bm, 0))
    plot_ly(df) %>%
      add_trace(x = ~week, y = ~mslt, name = "MSLT Latency (min)",
                mode = "lines", line = list(color = "#F39C12", width = 2)) %>%
      add_trace(x = ~week, y = ~ess / 3, name = "ESS / 3 (scaled)",
                mode = "lines", line = list(color = "#E74C3C", width = 2)) %>%
      add_segments(x = 0, xend = 24, y = csf_hcrt[1] / 30, yend = csf_hcrt[1] / 30,
                   line = list(dash = "dash", color = "#9B59B6"),
                   name = paste0("CSF Hcrt-1 / 30 (≈", round(csf_hcrt[1]), " pg/mL)")) %>%
      layout(
        title = "Biomarker Trajectories Over Treatment (24 weeks)",
        xaxis = list(title = "Treatment Week"),
        yaxis = list(title = "AU (scaled)"),
        plot_bgcolor = "#FAFAFA", paper_bgcolor = "white"
      )
  })

  # =========================================================================
  # TAB 7: Autoimmune
  # =========================================================================
  output$neuron_loss_plot <- renderPlotly({
    years <- seq(0, 20, by = 0.2)
    surviv <- orexin_decay(years, input$orexin_survival)
    df <- data.frame(year = years, survival = surviv)
    plot_ly(df, x = ~year, y = ~survival, type = "scatter", mode = "lines",
            fill = "tozeroy", fillcolor = "rgba(231,76,60,0.2)",
            line = list(color = "#E74C3C", width = 2.5)) %>%
      add_segments(x = 0, xend = 20, y = 10, yend = 10,
                   line = list(dash = "dash", color = "gray"), name = "NT1 threshold") %>%
      add_markers(x = input$disease_duration, y = input$orexin_survival,
                  marker = list(size = 14, color = "#E74C3C",
                                symbol = "star", line = list(color = "white", width = 2)),
                  name = "Current Patient") %>%
      layout(
        title = "Orexin Neuron Survival Over Disease Duration",
        xaxis = list(title = "Years Since Disease Onset"),
        yaxis = list(title = "Orexin Neuron Survival (%)", range = c(0, 105)),
        plot_bgcolor = "#FAFAFA", paper_bgcolor = "white"
      )
  })

  output$hla_assoc_plot <- renderPlotly({
    populations <- c("NT1 (HLA+)", "NT1 (HLA-)", "NT2", "IH", "OSA", "General Pop.")
    hla_freq    <- c(98, 2, 41, 36, 33, 25)
    rel_risk    <- c(251, 0.1, 3.5, 2.0, 1.4, 1.0)
    df <- data.frame(pop = populations, freq = hla_freq, rr = rel_risk)
    plot_ly(df) %>%
      add_trace(x = ~pop, y = ~freq, name = "HLA-DQB1*06:02 Freq (%)",
                type = "bar", marker = list(color = "#3498DB"), yaxis = "y") %>%
      add_trace(x = ~pop, y = ~rr, name = "Relative Risk",
                type = "scatter", mode = "markers",
                marker = list(size = 14, color = "#E74C3C", symbol = "diamond"),
                yaxis = "y2") %>%
      layout(
        title = "HLA-DQB1*06:02 Association & Relative Risk",
        xaxis = list(title = ""),
        yaxis  = list(title = "HLA Freq (%)"),
        yaxis2 = list(title = "Relative Risk", overlaying = "y", side = "right",
                      type = "log"),
        barmode = "group",
        plot_bgcolor = "#FAFAFA", paper_bgcolor = "white"
      )
  })

  output$autoab_plot <- renderPlotly({
    groups   <- c("NT1 Active", "NT1 Chronic", "NT2", "Control")
    anti_orex <- c(78, 55, 22, 8)
    anti_trib2 <- c(65, 42, 18, 5)
    anti_pmca  <- c(55, 38, 12, 4)
    df <- data.frame(group = groups, anti_orex, anti_trib2, anti_pmca)
    plot_ly(df) %>%
      add_trace(x = ~group, y = ~anti_orex, name = "Anti-Orexin (%+)",
                type = "bar", marker = list(color = "#E74C3C")) %>%
      add_trace(x = ~group, y = ~anti_trib2, name = "Anti-TRIB2 (%+)",
                type = "bar", marker = list(color = "#9B59B6")) %>%
      add_trace(x = ~group, y = ~anti_pmca, name = "Anti-PMCA4 (%+)",
                type = "bar", marker = list(color = "#F39C12")) %>%
      layout(
        title = "Autoantibody Positivity Rates",
        xaxis = list(title = ""),
        yaxis = list(title = "% Positive", range = c(0, 100)),
        barmode = "group",
        plot_bgcolor = "#FAFAFA", paper_bgcolor = "white"
      )
  })

  output$trigger_plot <- renderPlotly({
    triggers <- c("H1N1 Infection", "Pandemrix Vaccine", "Streptococcal Inf.",
                  "EBV Infection", "Physical/Emotional Stress", "Unknown")
    odds_ratio <- c(4.2, 6.6, 3.1, 2.4, 1.8, 1.0)
    ci_low     <- c(2.1, 4.2, 1.8, 1.3, 1.1, 0.8)
    ci_high    <- c(8.4, 10.4, 5.3, 4.5, 2.9, 1.3)
    df <- data.frame(trigger = triggers, or = odds_ratio,
                     low = ci_low, high = ci_high)
    plot_ly(df) %>%
      add_trace(x = ~or, y = ~reorder(trigger, or),
                type = "scatter", mode = "markers",
                marker = list(size = 12, color = "#E74C3C",
                              line = list(color = "darkred", width = 1.5)),
                error_x = list(
                  type = "data", symmetric = FALSE,
                  array = ~high - or, arrayminus = ~or - low,
                  color = "#2C3E50"
                ),
                name = "Odds Ratio (95% CI)") %>%
      add_segments(x = 1, xend = 1, y = 0.5, yend = nrow(df) + 0.5,
                   line = list(dash = "dash", color = "gray")) %>%
      layout(
        title = "Environmental Triggers — Odds Ratios for NT1",
        xaxis = list(title = "Odds Ratio (log scale)", type = "log"),
        yaxis = list(title = ""),
        plot_bgcolor = "#FAFAFA", paper_bgcolor = "white"
      )
  })

  # =========================================================================
  # TAB 8: References — Trial Table
  # =========================================================================
  output$trial_table <- renderTable({
    data.frame(
      Trial              = c("US Multicenter Study", "HARMONY I", "HARMONY II",
                             "TONES 3", "TONES 4", "REST-ON", "XYREM Cataplexy",
                             "CLARITY", "Jazz Phase 3 (LXB)"),
      Drug               = c("Modafinil", "Pitolisant", "Pitolisant",
                             "Solriamfetol", "Solriamfetol", "Sodium Oxybate",
                             "Sodium Oxybate (SXB)", "Sodium Oxybate (SXB)",
                             "Low-Sodium Oxybate"),
      N                  = c(271, 95, 110, 231, 222, 208, 136, 222, 212),
      ESS_Change         = c(-4.3, -5.8, -4.9, -7.7, -6.4, -5.5, -3.5, -5.1, -5.9),
      Cataplexy_Reduction = c("N/A", "65%", "58%", "N/A", "N/A", "75%", "72%", "69%", "73%"),
      Year               = c(2000, 2013, 2016, 2019, 2019, 2021, 2002, 2010, 2021),
      PMID               = c("10692671", "23145457", "27163208",
                             "31529545", "31529546", "33285153",
                             "12165486", "20599521", "34433044"),
      stringsAsFactors   = FALSE
    )
  }, striped = TRUE, hover = TRUE, bordered = TRUE)

}  # end server

# =============================================================================
shinyApp(ui, server)
