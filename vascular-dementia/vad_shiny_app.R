## ============================================================
## Vascular Dementia (VaD) — Shiny Interactive QSP Dashboard
## 혈관성 치매 정량적 시스템 약리학 인터랙티브 대시보드
## ============================================================
## 실행 방법: shiny::runApp("vad_shiny_app.R")
## ============================================================

library(shiny)
library(shinydashboard)
library(ggplot2)
library(dplyr)
library(tidyr)
library(DT)

## ── Simulation engine (self-contained Euler integrator) ────────────
simulate_VaD <- function(
    MMSE0, WMH0, CBF0, BP0, LDL0, ACh0, Syn0,
    use_AHT, use_APT, use_STATIN, use_AChEI, use_MEM, use_CIL,
    dose_AHT, dose_APT, dose_STATIN, dose_AChEI, dose_MEM, dose_CIL,
    n_days = 730, dt = 1
) {
  # PK steady-state approximations (Cp_ss = F*Dose/CL)
  Cp_AHT    <- if (use_AHT    == 1) dose_AHT    * 0.65 / 12.0 / 24 else 0
  Cp_APT    <- if (use_APT    == 1) dose_APT    * 0.80 / 8.0  / 24 else 0
  Cp_ST     <- if (use_STATIN == 1) dose_STATIN * 0.12 / 200.0/ 24 else 0
  Cp_AChEI  <- if (use_AChEI  == 1) dose_AChEI  * 0.90 / 3.5  / 24 * 15 else 0  # *Kp_brain
  Cp_MEM    <- if (use_MEM    == 1) dose_MEM    * 0.85 / 4.5  / 24 * 8  else 0
  Cp_CIL    <- if (use_CIL    == 1) dose_CIL    * 0.90 / 35.0 / 24 else 0

  # PD effects at steady-state
  E_AHT_BP  <- use_AHT    * 25.0 * (Cp_AHT^1.5)  / (0.8^1.5  + Cp_AHT^1.5)
  E_ST_LDL  <- use_STATIN * 0.50 * Cp_ST          / (1.2      + Cp_ST)
  E_APT_SVD <- use_APT    * 0.35 * Cp_APT         / (0.5      + Cp_APT)
  E_CIL_CBF <- use_CIL    * 0.18 * Cp_CIL         / (0.6      + Cp_CIL)
  E_ST_CBF  <- use_STATIN * 0.10 * Cp_ST          / (2.0      + Cp_ST)
  AChE_inhib<- use_AChEI  * 0.75 * Cp_AChEI       / (0.05     + Cp_AChEI)
  E_MEM     <- use_MEM    * 0.60 * Cp_MEM         / (0.10     + Cp_MEM)
  E_ST_inf  <- use_STATIN * 0.30 * Cp_ST          / (2.0      + Cp_ST)
  E_NIC_CBF <- 0  # not simulated separately here

  # State variables
  BP   <- BP0; LDL  <- LDL0; CBF  <- CBF0
  WMH  <- WMH0; Inf  <- 5.0
  MG   <- 0.25; CYT  <- 3.0; ROS  <- 2.5
  ACh  <- ACh0; Syn  <- Syn0; MMSE <- MMSE0

  # Adjust BP/LDL for drug effect (immediate)
  BP  <- BP  - E_AHT_BP
  LDL <- LDL * (1 - E_ST_LDL)

  # CBF with drug effects
  CBF <- CBF * (1 + E_CIL_CBF + E_ST_CBF)

  # Storage
  out <- data.frame(
    Day = numeric(n_days + 1),
    BP = numeric(n_days + 1),   LDL = numeric(n_days + 1),
    CBF = numeric(n_days + 1),  WMH = numeric(n_days + 1),
    Inf = numeric(n_days + 1),  MG  = numeric(n_days + 1),
    CYT = numeric(n_days + 1),  ROS = numeric(n_days + 1),
    ACh = numeric(n_days + 1),  Syn = numeric(n_days + 1),
    MMSE = numeric(n_days + 1)
  )
  out[1, ] <- c(0, BP, LDL, CBF, WMH, Inf, MG, CYT, ROS, ACh, Syn, MMSE)

  for (i in seq_len(n_days)) {
    # WMH
    svd_rate <- 0.0018 * (1 + 0.0012 * max(BP - 130, 0)) *
                           (1 + 0.0006 * max(LDL - 100, 0)) *
                           (1 - E_APT_SVD) * (1 - E_AHT_BP / 25 * 0.4)
    WMH <- WMH + svd_rate * dt

    # Microinfarcts
    Inf <- Inf + 0.0005 * (WMH / WMH0) * (1 - E_APT_SVD * 0.5) * dt

    # CBF dynamics
    CBF_target <- CBF0 * (1 - 0.15 * WMH / (WMH + 15)) *
                          (1 - 0.1 * max(BP - 130, 0) / 30) *
                          (1 + E_CIL_CBF + E_ST_CBF)
    dCBF <- 0.005 * (CBF_target - CBF)
    CBF <- CBF + dCBF * dt

    # Microglia
    hypox <- max(0.5 * (1 - CBF / 55), 0)
    dMG_in  <- 0.003 * (hypox + 0.2 * max(ROS / 2.5 - 1, 0)) * (1 - MG)
    dMG_out <- 0.002 * MG * (1 + E_ST_inf)
    MG <- max(0, min(1, MG + (dMG_in - dMG_out) * dt))

    # Cytokine
    CYT <- max(0.5, CYT + (0.5 * MG - 0.003 * CYT) * dt)

    # ROS
    ros_prod  <- 0.004 * (1 + 0.5 * (CYT / 3 - 1)) * max(1 - CBF / 60, 0.01)
    ros_clear <- 0.003 * ROS * (1 + E_ST_inf * 0.3)
    ROS <- max(0.5, ROS + (ros_prod - ros_clear) * dt)

    # ACh
    ach_loss <- 0.0001 * (1 + 0.15 * max(1 - CBF / 55, 0)) * ACh
    ach_rest <- AChE_inhib * 0.30 * (1 - ACh)
    ACh <- max(0.1, min(1.0, ACh + (-ach_loss + ach_rest) * dt))

    # Synaptic density
    syn_loss <- 0.00015 * (1 + 0.10 * max(ROS / 2.5 - 1, 0) +
                                0.08 * max(CYT / 3 - 1, 0)) *
                            (1 - E_MEM * 0.3) * Syn
    syn_rest <- 0.00005 * (1 + E_ST_CBF * 0.5 + AChE_inhib * 0.2) *
                            max(0.9 - Syn, 0)
    Syn <- max(0.1, min(1.0, Syn + (-syn_loss + syn_rest) * dt))

    # MMSE
    MMSE_exp <- 8.0 * ACh + 12.0 * Syn -
                0.3 * max(WMH - WMH0, 0) - 0.4 * max(Inf - 5, 0)
    MMSE <- max(0, min(30, MMSE + 0.002 * (MMSE_exp - MMSE) * dt))

    out[i + 1, ] <- c(i, BP, LDL, CBF, WMH, Inf, MG, CYT, ROS, ACh, Syn, MMSE)
  }
  out
}

## ── UI ──────────────────────────────────────────────────────────────

ui <- dashboardPage(
  skin = "blue",
  dashboardHeader(title = "VaD QSP Dashboard — 혈관성 치매"),

  dashboardSidebar(
    sidebarMenu(
      menuItem("1. 환자 프로파일",     tabName = "profile",  icon = icon("user")),
      menuItem("2. Drug PK/PD",        tabName = "pk",       icon = icon("pills")),
      menuItem("3. 혈관·뇌관류",       tabName = "vascular", icon = icon("heart")),
      menuItem("4. 신경생물학적 기전", tabName = "neuro",    icon = icon("brain")),
      menuItem("5. 임상 엔드포인트",   tabName = "outcomes", icon = icon("chart-line")),
      menuItem("6. 시나리오 비교",     tabName = "scenario", icon = icon("sliders-h")),
      menuItem("7. 바이오마커",        tabName = "biomarker",icon = icon("flask"))
    )
  ),

  dashboardBody(
    tabItems(

      ## ── Tab 1: Patient Profile ────────────────────────────────────
      tabItem(tabName = "profile",
        fluidRow(
          box(title = "환자 기저 특성", width = 4, status = "primary", solidHeader = TRUE,
            sliderInput("MMSE0",  "기저 MMSE 점수 (0-30)",      22, min = 10, max = 28, step = 1),
            sliderInput("WMH0",   "기저 WMH 부피 (mL)",          8, min = 1,  max = 30, step = 0.5),
            sliderInput("BP0",    "기저 수축기 혈압 (mmHg)",    145, min = 110,max = 185, step = 5),
            sliderInput("CBF0",   "기저 CBF (mL/100g/min)",      50, min = 30, max = 70, step = 2),
            sliderInput("LDL0",   "기저 LDL-C (mg/dL)",        130, min = 60, max = 220, step = 5),
            sliderInput("ACh0",   "기저 ACh 톤 (rel, 0-1)",   0.60, min = 0.2, max = 1.0, step = 0.05),
            sliderInput("Syn0",   "기저 시냅스 밀도 (rel, 0-1)",0.75, min = 0.3, max = 1.0, step = 0.05)
          ),
          box(title = "위험 인자 요약", width = 4, status = "warning", solidHeader = TRUE,
            htmlOutput("riskSummary")
          ),
          box(title = "VaD 진단 기준 (VASCOG 2014)", width = 4, status = "info",
            h5("주요 진단 특징:"),
            tags$ul(
              tags$li("신경인지 장애의 증거 (MMSE / MoCA / NTB)"),
              tags$li("뇌혈관 질환의 영상 증거 (MRI)"),
              tags$li("임상적 인과 관계 (시간적·공간적 연관성)"),
              tags$li("치매 이외 원인 배제")
            ),
            h5("주요 아형:"),
            tags$ul(
              tags$li(strong("소혈관 질환형"), " — 백질변성, 소경색"),
              tags$li(strong("전략적 경색형"), " — 시상·해마·각회"),
              tags$li(strong("피질 경색형"), " — 다발성 뇌경색"),
              tags$li(strong("혼합형"), " — 알츠하이머 + 혈관성")
            )
          )
        )
      ),

      ## ── Tab 2: Drug PK/PD ─────────────────────────────────────────
      tabItem(tabName = "pk",
        fluidRow(
          box(title = "약물 선택 & 용량", width = 4, status = "success", solidHeader = TRUE,
            h5("혈관 위험인자 조절"),
            checkboxInput("use_AHT",    "강압제 (Antihypertensive)",  TRUE),
            sliderInput("dose_AHT",     "용량 (mg/day)", 10, min = 2.5, max = 20, step = 2.5),
            checkboxInput("use_APT",    "항혈소판제 (Aspirin)",  TRUE),
            sliderInput("dose_APT",     "용량 (mg/day)", 100, min = 75, max = 325, step = 25),
            checkboxInput("use_STATIN", "스타틴 (Statin)", TRUE),
            sliderInput("dose_STATIN",  "용량 (mg/day)",  40, min = 10, max = 80, step = 10),
            hr(),
            h5("증상 치료"),
            checkboxInput("use_AChEI",  "아세틸콜린에스터라제 억제제 (AChEI)", FALSE),
            sliderInput("dose_AChEI",   "용량 (mg/day)", 10, min = 5, max = 23, step = 5),
            checkboxInput("use_MEM",    "메만틴 (Memantine)", FALSE),
            sliderInput("dose_MEM",     "용량 (mg/day)", 20, min = 5, max = 20, step = 5),
            hr(),
            h5("기타"),
            checkboxInput("use_CIL",    "실로스타졸 (Cilostazol)", FALSE),
            sliderInput("dose_CIL",     "용량 (mg/day)", 200, min = 50, max = 200, step = 50)
          ),
          box(title = "약동력학: 혈중 농도 (정상상태)", width = 8, status = "success",
            plotOutput("pkPlot", height = "450px")
          )
        ),
        fluidRow(
          box(title = "약역학 효과 요약", width = 12, status = "primary",
            tableOutput("pdSummary")
          )
        )
      ),

      ## ── Tab 3: Vascular / Cerebral Perfusion ──────────────────────
      tabItem(tabName = "vascular",
        fluidRow(
          box(title = "뇌관류 & 혈관 기전 시뮬레이션", width = 12, status = "danger",
            sliderInput("sim_duration", "시뮬레이션 기간 (일)", 730, min = 90, max = 1460, step = 90),
            actionButton("run_sim", "시뮬레이션 실행", class = "btn-primary btn-lg")
          )
        ),
        fluidRow(
          box(title = "수축기 혈압 (mmHg)", width = 6, status = "danger",
            plotOutput("bpPlot", height = "280px")),
          box(title = "뇌혈류 CBF (mL/100g/min)", width = 6, status = "warning",
            plotOutput("cbfPlot", height = "280px"))
        ),
        fluidRow(
          box(title = "백질변성 WMH 진행 (mL)", width = 6, status = "warning",
            plotOutput("wmhPlot", height = "280px")),
          box(title = "미세경색 누적 (n)", width = 6, status = "danger",
            plotOutput("infPlot", height = "280px"))
        )
      ),

      ## ── Tab 4: Neurobiological Mechanisms ─────────────────────────
      tabItem(tabName = "neuro",
        fluidRow(
          box(title = "신경생물학적 기전 다이나믹스", width = 12, status = "info",
            p("혈관성 치매의 핵심 신경생물학적 메커니즘을 시각화합니다.")
          )
        ),
        fluidRow(
          box(title = "미세아교세포(M1) 활성화 & 사이토카인 지수", width = 6, status = "warning",
            plotOutput("microglia_plot", height = "280px")),
          box(title = "산화 스트레스 지수 (ROS)", width = 6, status = "success",
            plotOutput("rosPlot", height = "280px"))
        ),
        fluidRow(
          box(title = "아세틸콜린 톤 & AChE 억제 효과", width = 6, status = "info",
            plotOutput("achPlot", height = "280px")),
          box(title = "시냅스 밀도 변화", width = 6, status = "info",
            plotOutput("synPlot", height = "280px"))
        )
      ),

      ## ── Tab 5: Clinical Endpoints ──────────────────────────────────
      tabItem(tabName = "outcomes",
        fluidRow(
          box(title = "MMSE 점수 궤적", width = 8, status = "success", solidHeader = TRUE,
            plotOutput("mmsePlot", height = "400px")),
          box(title = "결과 요약 (6·12·24개월)", width = 4, status = "primary",
            tableOutput("outcomeTable"))
        ),
        fluidRow(
          box(title = "인지 도메인별 기여도", width = 6, status = "info",
            plotOutput("cogDomainPlot", height = "300px")),
          box(title = "혈관성 치매 스테이지 분류", width = 6, status = "warning",
            htmlOutput("stageHTML"))
        )
      ),

      ## ── Tab 6: Scenario Comparison ────────────────────────────────
      tabItem(tabName = "scenario",
        fluidRow(
          box(title = "6가지 치료 시나리오 비교", width = 12, status = "primary", solidHeader = TRUE,
            p("기저 환자 프로파일에서 6가지 치료 전략의 2년 결과를 비교합니다.")
          )
        ),
        fluidRow(
          box(title = "MMSE 점수 비교 (2년)", width = 6, status = "success",
            plotOutput("scenMMSE", height = "380px")),
          box(title = "WMH 진행 비교 (2년)", width = 6, status = "warning",
            plotOutput("scenWMH", height = "380px"))
        ),
        fluidRow(
          box(title = "시나리오별 2년 결과 테이블", width = 12, status = "primary",
            DTOutput("scenTable"))
        )
      ),

      ## ── Tab 7: Biomarker Panel ────────────────────────────────────
      tabItem(tabName = "biomarker",
        fluidRow(
          box(title = "바이오마커 패널 & 치료 반응 예측", width = 12, status = "info",
            p("혈관성 치매 핵심 바이오마커의 치료 전후 변화를 시각화합니다.")
          )
        ),
        fluidRow(
          box(title = "뇌 MRI 바이오마커", width = 6, status = "primary",
            plotOutput("mri_biomarker", height = "320px")),
          box(title = "혈액/CSF 바이오마커 변화", width = 6, status = "info",
            plotOutput("csf_biomarker", height = "320px"))
        ),
        fluidRow(
          box(title = "바이오마커 정상 참고치 & VaD 범위", width = 12, status = "warning",
            tableOutput("biomarkerRef"))
        )
      )

    ) # end tabItems
  ) # end dashboardBody
)

## ── Server ──────────────────────────────────────────────────────────

server <- function(input, output, session) {

  ## Reactive simulation
  sim_result <- eventReactive(input$run_sim, {
    simulate_VaD(
      MMSE0 = input$MMSE0, WMH0 = input$WMH0, CBF0 = input$CBF0,
      BP0 = input$BP0, LDL0 = input$LDL0,
      ACh0 = input$ACh0, Syn0 = input$Syn0,
      use_AHT = as.integer(input$use_AHT),
      use_APT = as.integer(input$use_APT),
      use_STATIN = as.integer(input$use_STATIN),
      use_AChEI = as.integer(input$use_AChEI),
      use_MEM = as.integer(input$use_MEM),
      use_CIL = as.integer(input$use_CIL),
      dose_AHT = input$dose_AHT, dose_APT = input$dose_APT,
      dose_STATIN = input$dose_STATIN, dose_AChEI = input$dose_AChEI,
      dose_MEM = input$dose_MEM, dose_CIL = input$dose_CIL,
      n_days = input$sim_duration
    )
  }, ignoreNULL = FALSE)

  ## ── Tab 1: Risk summary ──────────────────────────────────────────
  output$riskSummary <- renderUI({
    bp_risk  <- if (input$BP0 >= 160) "매우 높음" else if (input$BP0 >= 140) "높음" else "보통"
    mmse_cat <- if (input$MMSE0 >= 24) "경미" else if (input$MMSE0 >= 18) "경도" else "중등도"
    wmh_cat  <- if (input$WMH0 < 5) "경미 (1등급)" else if (input$WMH0 < 15) "중등도 (2-3등급)" else "심각 (4등급)"
    HTML(paste0(
      "<b>혈압:</b> ", input$BP0, " mmHg (", bp_risk, ")<br>",
      "<b>LDL-C:</b> ", input$LDL0, " mg/dL<br>",
      "<b>인지 단계:</b> ", mmse_cat, " (MMSE=", input$MMSE0, ")<br>",
      "<b>WMH 등급:</b> ", wmh_cat, "<br>",
      "<b>CBF:</b> ", input$CBF0, " mL/100g/min<br>",
      "<b>ACh 톤:</b> ", input$ACh0, " (정상 = 1.0)<br>",
      "<b>시냅스 밀도:</b> ", input$Syn0, " (정상 = 1.0)"
    ))
  })

  ## ── Tab 2: PK Plot ───────────────────────────────────────────────
  output$pkPlot <- renderPlot({
    t <- seq(0, 48, 0.5)  # 48h PK profile
    drugs <- list()
    if (input$use_AHT)    drugs[["강압제 (ARB/ACEi)"]] <-
      input$dose_AHT * 0.65 / 12 * exp(-0.25 * t) + input$dose_AHT * 0.65 / 12 * 0.3 * exp(-0.08 * t)
    if (input$use_APT)    drugs[["항혈소판제 (ASA)"]] <-
      input$dose_APT * 0.80 / 8  * exp(-0.35 * t)
    if (input$use_STATIN) drugs[["스타틴"]] <-
      input$dose_STATIN * 0.12 / 200 * exp(-0.20 * t)
    if (input$use_AChEI)  drugs[["AChEI (도네페질)"]] <-
      input$dose_AChEI * 0.90 / 3.5 * exp(-0.029 * t)
    if (input$use_MEM)    drugs[["메만틴"]] <-
      input$dose_MEM   * 0.85 / 4.5 * exp(-0.033 * t)
    if (input$use_CIL)    drugs[["실로스타졸"]] <-
      input$dose_CIL   * 0.90 / 35  * exp(-0.15 * t)
    if (length(drugs) == 0) {
      plot(0, type = "n", xlab = "Time (h)", ylab = "Cp (μg/mL)", main = "약물 없음")
      return()
    }
    df <- do.call(rbind, lapply(names(drugs), function(nm)
      data.frame(t = t, Cp = pmax(drugs[[nm]], 0), Drug = nm)))
    ggplot(df, aes(x = t, y = Cp, color = Drug)) +
      geom_line(size = 1.2) +
      labs(x = "Time (h)", y = "Cp (μg/mL)", title = "약동학: 혈중 농도-시간 프로파일 (48h)") +
      theme_bw() + theme(legend.position = "bottom", legend.text = element_text(size = 9))
  })

  output$pdSummary <- renderTable({
    Cp_AHT    <- if (input$use_AHT)    input$dose_AHT    * 0.65 / 12.0 / 24 else 0
    Cp_APT    <- if (input$use_APT)    input$dose_APT    * 0.80 / 8.0  / 24 else 0
    Cp_ST     <- if (input$use_STATIN) input$dose_STATIN * 0.12 / 200.0/ 24 else 0
    Cp_AChEI  <- if (input$use_AChEI)  input$dose_AChEI  * 0.90 / 3.5  / 24 * 15 else 0
    Cp_MEM    <- if (input$use_MEM)    input$dose_MEM    * 0.85 / 4.5  / 24 * 8  else 0
    Cp_CIL    <- if (input$use_CIL)    input$dose_CIL    * 0.90 / 35.0 / 24 else 0

    E_AHT_BP <- round(input$use_AHT * 25 * Cp_AHT^1.5 / (0.8^1.5 + Cp_AHT^1.5), 1)
    E_LDL    <- round(input$use_STATIN * 100 * 0.50 * Cp_ST / (1.2 + Cp_ST), 1)
    E_AChE   <- round(input$use_AChEI * 100 * 0.75 * Cp_AChEI / (0.05 + Cp_AChEI), 1)
    E_NMDA   <- round(input$use_MEM * 100 * 0.60 * Cp_MEM / (0.10 + Cp_MEM), 1)
    E_CBF    <- round((input$use_CIL * 18 * Cp_CIL / (0.6 + Cp_CIL)) +
                      (input$use_STATIN * 10 * Cp_ST / (2.0 + Cp_ST)), 1)

    data.frame(
      약물 = c("강압제 (ARB/ACEi)", "스타틴", "AChEI (도네페질)", "메만틴", "실로스타졸"),
      `사용` = c(input$use_AHT, input$use_STATIN, input$use_AChEI, input$use_MEM, input$use_CIL),
      `효과 지표` = c("SBP 감소 (mmHg)", "LDL 감소 (%)", "AChE 억제 (%)", "NMDA 차단 (%)", "CBF 증가 (%)"),
      `예측 효과` = c(E_AHT_BP, E_LDL, E_AChE, E_NMDA, E_CBF),
      stringsAsFactors = FALSE
    )
  })

  ## ── Vascular Tab ─────────────────────────────────────────────────
  output$bpPlot <- renderPlot({
    df <- sim_result()
    ggplot(df, aes(x = Day, y = BP)) +
      geom_line(color = "#C62828", size = 1.2) +
      geom_hline(yintercept = 130, linetype = "dashed", color = "steelblue") +
      annotate("text", x = 5, y = 131.5, label = "목표: 130 mmHg",
               hjust = 0, color = "steelblue", size = 3.5) +
      labs(x = "일(Day)", y = "SBP (mmHg)", title = "수축기 혈압 추이") +
      theme_bw()
  })

  output$cbfPlot <- renderPlot({
    df <- sim_result()
    ggplot(df, aes(x = Day, y = CBF)) +
      geom_line(color = "#1565C0", size = 1.2) +
      geom_hline(yintercept = input$CBF0, linetype = "dashed", color = "gray50") +
      labs(x = "일(Day)", y = "CBF (mL/100g/min)", title = "뇌혈류량 변화") +
      theme_bw()
  })

  output$wmhPlot <- renderPlot({
    df <- sim_result()
    ggplot(df, aes(x = Day, y = WMH)) +
      geom_area(fill = "#FF8A65", alpha = 0.4) +
      geom_line(color = "#BF360C", size = 1.2) +
      labs(x = "일(Day)", y = "WMH 부피 (mL)", title = "백질변성 (WMH) 진행") +
      theme_bw()
  })

  output$infPlot <- renderPlot({
    df <- sim_result()
    ggplot(df, aes(x = Day, y = Inf)) +
      geom_line(color = "#6A1B9A", size = 1.2) +
      labs(x = "일(Day)", y = "미세경색 수 (누적)", title = "피질 미세경색 누적") +
      theme_bw()
  })

  ## ── Neuro Tab ────────────────────────────────────────────────────
  output$microglia_plot <- renderPlot({
    df <- sim_result()
    df2 <- df %>% select(Day, MG, CYT) %>%
      pivot_longer(-Day, names_to = "Marker", values_to = "Value") %>%
      mutate(Marker = recode(Marker, "MG" = "M1 미세아교세포 (0-1)", "CYT" = "사이토카인 지수 (AU)"))
    ggplot(df2, aes(x = Day, y = Value, color = Marker)) +
      geom_line(size = 1.1) +
      facet_wrap(~Marker, scales = "free_y", ncol = 1) +
      labs(x = "일", y = "상태값", title = "신경염증 기전") +
      theme_bw() + theme(legend.position = "none")
  })

  output$rosPlot <- renderPlot({
    df <- sim_result()
    ggplot(df, aes(x = Day, y = ROS)) +
      geom_line(color = "#2E7D32", size = 1.2) +
      geom_hline(yintercept = 1.0, linetype = "dashed", color = "gray50") +
      annotate("text", x = 5, y = 1.05, label = "정상 = 1.0", hjust = 0, size = 3.5) +
      labs(x = "일(Day)", y = "ROS 지수 (AU)", title = "산화 스트레스 지수") +
      theme_bw()
  })

  output$achPlot <- renderPlot({
    df <- sim_result()
    ggplot(df, aes(x = Day, y = ACh)) +
      geom_line(color = "#1565C0", size = 1.2) +
      geom_hline(yintercept = 1.0, linetype = "dashed", color = "gray50") +
      scale_y_continuous(limits = c(0, 1.05)) +
      labs(x = "일(Day)", y = "ACh 톤 (rel, 정상=1)", title = "아세틸콜린 톤") +
      theme_bw()
  })

  output$synPlot <- renderPlot({
    df <- sim_result()
    ggplot(df, aes(x = Day, y = Syn)) +
      geom_area(fill = "#7B1FA2", alpha = 0.3) +
      geom_line(color = "#4A148C", size = 1.2) +
      scale_y_continuous(limits = c(0, 1.05)) +
      labs(x = "일(Day)", y = "시냅스 밀도 (rel)", title = "시냅스 밀도 변화") +
      theme_bw()
  })

  ## ── Outcomes Tab ─────────────────────────────────────────────────
  output$mmsePlot <- renderPlot({
    df <- sim_result()
    ggplot(df, aes(x = Day, y = MMSE)) +
      geom_ribbon(aes(ymin = MMSE - 1.5, ymax = MMSE + 1.5), fill = "#A5D6A7", alpha = 0.3) +
      geom_line(color = "#1B5E20", size = 1.5) +
      geom_hline(yintercept = input$MMSE0, linetype = "dashed", color = "gray60") +
      geom_hline(yintercept = 10, linetype = "dotted", color = "red", size = 0.8) +
      annotate("text", x = 5, y = 10.5, label = "중등도 치매 경계 (10)", hjust = 0, color = "red", size = 3) +
      scale_y_continuous(limits = c(0, 32)) +
      labs(x = "일(Day)", y = "MMSE 점수 (0-30)", title = "MMSE 점수 궤적",
           subtitle = paste0("기저값: ", input$MMSE0)) +
      theme_bw(base_size = 14)
  })

  output$outcomeTable <- renderTable({
    df <- sim_result()
    tpts <- c(180, 365, 730)
    df %>% filter(Day %in% tpts) %>%
      mutate(
        `시점` = paste0(round(Day / 30.4, 0), "개월"),
        `MMSE` = round(MMSE, 1),
        `변화` = round(MMSE - input$MMSE0, 1),
        `WMH(mL)` = round(WMH, 1),
        `CBF` = round(CBF, 1)
      ) %>%
      select(`시점`, `MMSE`, `변화`, `WMH(mL)`, `CBF`)
  })

  output$cogDomainPlot <- renderPlot({
    df <- sim_result()
    last_row <- tail(df, 1)
    domains <- data.frame(
      Domain = c("에피소드 기억", "실행 기능", "처리 속도", "언어 기능", "시공간 기능"),
      Score  = c(
        last_row$ACh * 0.9,
        last_row$Syn * 0.8 * (last_row$CBF / input$CBF0),
        last_row$CBF / input$CBF0 * 0.95,
        last_row$Syn * 0.85,
        last_row$Syn * 0.80 * (last_row$CBF / input$CBF0)
      )
    )
    ggplot(domains, aes(x = Domain, y = Score, fill = Domain)) +
      geom_col() +
      geom_hline(yintercept = 1.0, linetype = "dashed", color = "gray40") +
      scale_y_continuous(limits = c(0, 1.1)) +
      labs(x = "", y = "상대적 기능 (정상=1)", title = "2년 후 인지 도메인별 수준") +
      theme_bw() + theme(legend.position = "none",
                         axis.text.x = element_text(angle = 20, hjust = 1, size = 9))
  })

  output$stageHTML <- renderUI({
    df <- sim_result()
    last_mmse <- tail(df$MMSE, 1)
    stage <- if (last_mmse >= 24) "MCI / 경미"
             else if (last_mmse >= 18) "경도 치매"
             else if (last_mmse >= 10) "중등도 치매"
             else "중증 치매"
    color <- if (last_mmse >= 24) "success" else if (last_mmse >= 18) "info"
              else if (last_mmse >= 10) "warning" else "danger"
    HTML(paste0(
      "<div class='alert alert-", color, "'>",
      "<strong>2년 후 단계: ", stage, "</strong><br>",
      "MMSE: ", round(last_mmse, 1), " (초기: ", input$MMSE0, ")<br>",
      "변화량: ", round(last_mmse - input$MMSE0, 1), " 점<br>",
      "WMH: ", round(tail(df$WMH, 1), 1), " mL<br>",
      "CBF: ", round(tail(df$CBF, 1), 1), " mL/100g/min",
      "</div>"
    ))
  })

  ## ── Scenario Tab ─────────────────────────────────────────────────
  scenario_data <- reactive({
    base_args <- list(
      MMSE0 = input$MMSE0, WMH0 = input$WMH0, CBF0 = input$CBF0,
      BP0 = input$BP0, LDL0 = input$LDL0,
      ACh0 = input$ACh0, Syn0 = input$Syn0, n_days = 730
    )
    scen_list <- list(
      list(label = "1. 무치료", use_AHT=0,use_APT=0,use_STATIN=0,use_AChEI=0,use_MEM=0,use_CIL=0),
      list(label = "2. 강압제 단독", use_AHT=1,use_APT=0,use_STATIN=0,use_AChEI=0,use_MEM=0,use_CIL=0),
      list(label = "3. 혈관 병합\n(강압+항혈+스타틴)", use_AHT=1,use_APT=1,use_STATIN=1,use_AChEI=0,use_MEM=0,use_CIL=0),
      list(label = "4. 증상치료\n(AChEI+메만틴)", use_AHT=0,use_APT=0,use_STATIN=0,use_AChEI=1,use_MEM=1,use_CIL=0),
      list(label = "5. 포괄적 치료\n(혈관+증상+실로)", use_AHT=1,use_APT=1,use_STATIN=1,use_AChEI=1,use_MEM=1,use_CIL=1),
      list(label = "6. 최적+\n(고용량스타틴)", use_AHT=1,use_APT=1,use_STATIN=1,use_AChEI=1,use_MEM=1,use_CIL=1)
    )
    dose_args <- list(
      dose_AHT=10, dose_APT=100, dose_STATIN=40, dose_AChEI=10, dose_MEM=20, dose_CIL=200
    )

    bind_rows(lapply(scen_list, function(sc) {
      args <- c(base_args, sc[setdiff(names(sc), "label")], dose_args)
      if (sc$label == "6. 최적+\n(고용량스타틴)") args$dose_STATIN <- 80
      df <- do.call(simulate_VaD, args)
      df$Scenario <- sc$label
      df
    }))
  })

  output$scenMMSE <- renderPlot({
    df <- scenario_data()
    ggplot(df, aes(x = Day, y = MMSE, color = Scenario)) +
      geom_line(size = 1.1) +
      scale_y_continuous(limits = c(max(0, input$MMSE0 - 10), input$MMSE0 + 1)) +
      labs(x = "일(Day)", y = "MMSE 점수", title = "시나리오별 MMSE 추이") +
      theme_bw() + theme(legend.position = "bottom", legend.text = element_text(size = 8))
  })

  output$scenWMH <- renderPlot({
    df <- scenario_data()
    ggplot(df, aes(x = Day, y = WMH, color = Scenario)) +
      geom_line(size = 1.1) +
      labs(x = "일(Day)", y = "WMH 부피 (mL)", title = "시나리오별 WMH 진행") +
      theme_bw() + theme(legend.position = "bottom", legend.text = element_text(size = 8))
  })

  output$scenTable <- renderDT({
    df <- scenario_data() %>%
      filter(Day == 730) %>%
      mutate(
        `MMSE (24개월)` = round(MMSE, 1),
        `MMSE 변화`     = round(MMSE - input$MMSE0, 1),
        `WMH (24개월, mL)` = round(WMH, 1),
        `WMH 변화`      = round(WMH - input$WMH0, 1),
        `CBF (24개월)`  = round(CBF, 1),
        `SBP`           = round(BP, 1)
      ) %>%
      select(Scenario, `MMSE (24개월)`, `MMSE 변화`,
             `WMH (24개월, mL)`, `WMH 변화`, `CBF (24개월)`, `SBP`)
    datatable(df, options = list(pageLength = 10, dom = 't'), rownames = FALSE)
  })

  ## ── Biomarker Tab ────────────────────────────────────────────────
  output$mri_biomarker <- renderPlot({
    df <- sim_result()
    last <- tail(df, 1)
    bio <- data.frame(
      Marker = c("WMH 부피 (mL)", "미세경색 수\n(누적)", "CBF\n(mL/100g/min)"),
      Baseline = c(input$WMH0, 5, input$CBF0),
      End_of_study = c(last$WMH, last$Inf, last$CBF)
    ) %>%
      pivot_longer(-Marker, names_to = "Timepoint", values_to = "Value") %>%
      mutate(Timepoint = recode(Timepoint,
                                "Baseline" = "기저", "End_of_study" = "24개월 후"))
    ggplot(bio, aes(x = Marker, y = Value, fill = Timepoint)) +
      geom_col(position = "dodge") +
      scale_fill_manual(values = c("기저" = "#90CAF9", "24개월 후" = "#1565C0")) +
      labs(x = "", y = "값", title = "MRI 영상 바이오마커 (기저 vs 24개월)") +
      theme_bw() + theme(axis.text.x = element_text(angle = 10, hjust = 1))
  })

  output$csf_biomarker <- renderPlot({
    df <- sim_result()
    last <- tail(df, 1)
    bio <- data.frame(
      Marker = c("ACh 톤\n(rel)", "시냅스 밀도\n(rel)", "M1 미세아교세포\n(0-1)", "ROS 지수\n(AU)"),
      Baseline = c(input$ACh0, input$Syn0, 0.25, 2.5),
      End_of_study = c(last$ACh, last$Syn, last$MG, last$ROS)
    ) %>%
      pivot_longer(-Marker, names_to = "Timepoint", values_to = "Value") %>%
      mutate(Timepoint = recode(Timepoint,
                                "Baseline" = "기저", "End_of_study" = "24개월 후"))
    ggplot(bio, aes(x = Marker, y = Value, fill = Timepoint)) +
      geom_col(position = "dodge") +
      scale_fill_manual(values = c("기저" = "#C8E6C9", "24개월 후" = "#2E7D32")) +
      labs(x = "", y = "값", title = "신경생물학 바이오마커 (기저 vs 24개월)") +
      theme_bw() + theme(axis.text.x = element_text(angle = 10, hjust = 1))
  })

  output$biomarkerRef <- renderTable({
    data.frame(
      `바이오마커` = c(
        "WMH 부피 (MRI)",
        "뇌혈류 (CBF, ASL-MRI)",
        "CSF p-Tau 181",
        "CSF Aβ42",
        "혈청 NfL",
        "hsCRP",
        "MMSE",
        "CDR-SB"
      ),
      `정상 참고치` = c(
        "< 2 mL (Fazekas 1-2등급)",
        "50-80 mL/100g/min",
        "< 80 pg/mL",
        "> 550 pg/mL",
        "< 10 pg/mL (성인)",
        "< 1.0 mg/L",
        "27-30 (정상 인지)",
        "0 (정상 CDR)"
      ),
      `VaD 범위` = c(
        "5-40+ mL (중증)",
        "30-55 mL/100g/min",
        "80-250 pg/mL",
        "200-500 pg/mL",
        "20-80+ pg/mL",
        "2-5 mg/L",
        "10-24 (경도-중등도)",
        "2-12 (경도-중등도)"
      ),
      `임상 의의` = c(
        "WMH ↑ = 소혈관 질환, 사망위험 ↑",
        "CBF ↓ = 뇌혈류 감소, 인지 ↓",
        "p-Tau ↑ = tau 병리, 혼합형 감별",
        "Aβ42 ↓ = 아밀로이드 침착, 혼합형",
        "NfL ↑ = 신경축삭 손상 지표",
        "CRP ↑ = 전신 염증, 뇌혈관 위험",
        "MMSE ↓ = 인지 저하 정도",
        "CDR-SB ↑ = 기능 저하 중증도"
      ),
      stringsAsFactors = FALSE,
      check.names = FALSE
    )
  })

} # end server

shinyApp(ui, server)
