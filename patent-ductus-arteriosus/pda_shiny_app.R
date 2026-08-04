## =====================================================================
##  PATENT DUCTUS ARTERIOSUS OF PREMATURITY -- QSP EXPLORER (Shiny)
## =====================================================================
##
##  A dashboard over pda_mrgsolve_model.R.  The app is organised around
##  the question the field actually argues about, which is NOT "does the
##  drug close the duct" (it does) but "does closing it help" (mostly it
##  has not).  So the tabs deliberately separate three things that the
##  clinical literature keeps fusing:
##
##    * DUCTAL RESPONSE      -- tone, diameter, closure          (tab 3)
##    * HAEMODYNAMIC BURDEN  -- integral of significant shunt     (tab 4)
##    * OUTCOME              -- death / BPD / NEC / IVH           (tab 5)
##
##  Tab 8 is the one to read first: it shows the untreated duct closing on
##  its own, which is the competing process that caps every treatment
##  effect and which the negative trials were really measuring.
##
##  Run:  shiny::runApp("pda_shiny_app.R")
##  Requires: shiny, mrgsolve, dplyr, ggplot2, tidyr, DT
## =====================================================================

library(shiny)
library(mrgsolve)
library(dplyr)
library(tidyr)
library(ggplot2)

have_DT <- requireNamespace("DT", quietly = TRUE)

## ---------------------------------------------------------------------
## Load the model.  Sourcing pda_mrgsolve_model.R defines `mod`, the
## rx_* regimen builders, SCENARIOS and the helper functions, so the app
## and the batch model can never diverge.
## ---------------------------------------------------------------------
.here <- tryCatch(dirname(sys.frame(1)$ofile), error = function(e) ".")
if (is.null(.here) || !nzchar(.here)) .here <- "."
model_file <- file.path(.here, "pda_mrgsolve_model.R")
if (!file.exists(model_file)) model_file <- "pda_mrgsolve_model.R"
source(model_file, local = FALSE, chdir = TRUE)

THEME <- theme_bw(base_size = 12) +
  theme(panel.grid.minor = element_blank(),
        strip.background = element_rect(fill = "#eceff1", colour = NA),
        legend.position = "bottom")

PAL <- c("ibuprofen" = "#00796b", "indomethacin" = "#00838f",
         "acetaminophen" = "#558b2f", "expectant" = "#90a4ae")

## =====================================================================
##  UI
## =====================================================================
ui <- fluidPage(
  titlePanel("미숙아 동맥관 개존증 (PDA) QSP 시뮬레이터 · Patent Ductus Arteriosus of Prematurity"),
  tags$p(style = "color:#555; margin-top:-8px;",
    HTML("교육·연구 목적 모델입니다. 임상 의사결정에 사용할 수 없습니다. &nbsp;|&nbsp;
          Educational model — not for clinical use.")),

  sidebarLayout(
    sidebarPanel(
      width = 3,
      h4("환자 (Patient)"),
      sliderInput("GA", "재태연령 Gestational age (wk)", 23, 40, 26, step = 0.5),
      sliderInput("BW", "출생체중 Birth weight (kg)", 0.4, 3.5, 0.80, step = 0.02),
      sliderInput("PAO2", "PaO₂ (mmHg)", 30, 90, 55, step = 1),
      checkboxInput("ANTESTER", "산전 스테로이드 완료 (antenatal steroids)", TRUE),
      checkboxInput("SEPSIS",
        "융모양막염/조기 패혈증 — 관벽 peroxide 상승", FALSE),
      checkboxInput("HCORT", "조기 하이드로코르티손 병용", FALSE),

      hr(),
      h4("치료 (Treatment)"),
      selectInput("drug", "약물", c("무치료 (expectant)" = "none",
                                    "이부프로펜 IV" = "ibu",
                                    "인도메타신 IV" = "ind",
                                    "아세트아미노펜 IV" = "apap",
                                    "이부프로펜 + 아세트아미노펜" = "combo")),
      sliderInput("start_d", "치료 시작 (postnatal day)", 0.25, 21, 2, step = 0.25),
      conditionalPanel("input.drug == 'ibu' || input.drug == 'combo'",
        radioButtons("ibu_mode", "이부프로펜 용법",
          c("표준 10-5-5 mg/kg q24h" = "std",
            "고용량 20-10-10 mg/kg q24h" = "high",
            "지속주입 continuous infusion" = "inf"), selected = "std")),
      conditionalPanel("input.drug == 'ind'",
        sliderInput("ind_dose", "인도메타신 초회 용량 (mg/kg)",
                    0.05, 0.30, 0.20, step = 0.05),
        sliderInput("ind_n", "투여 횟수", 1, 6, 3, step = 1)),
      conditionalPanel("input.drug == 'apap' || input.drug == 'combo'",
        sliderInput("apap_dose", "아세트아미노펜 1회 용량 (mg/kg)",
                    7.5, 25, 15, step = 2.5),
        sliderInput("apap_days", "투여 기간 (days)", 1, 7, 3, step = 1)),
      checkboxInput("second_course", "재개통 시 2차 치료 (d+5)", FALSE),

      hr(),
      sliderInput("horizon", "관찰 기간 (days)", 20, 120, 90, step = 5),
      actionButton("go", "시뮬레이션 실행", class = "btn-primary",
                   icon = icon("play")),
      hr(),
      htmlOutput("verdict")
    ),

    mainPanel(
      width = 9,
      tabsetPanel(
        id = "tabs", type = "tabs",

        ## ---------------------------------------------------------- 1
        tabPanel("① 환자 프로파일",
          br(),
          fluidRow(
            column(6, h4("이 재태연령에서 무엇이 물리적으로 가능한가"),
                   verbatimTextOutput("profile_txt")),
            column(6, h4("관 폐쇄 천장 (occlusion ceiling)"),
                   plotOutput("plot_ceiling", height = 300))),
          hr(),
          HTML("<p><b>핵심 구조:</b> 최대 수축 가능 정도(TMAXGA)는 재태연령의
            함수이며, 이것이 남기는 잔여 내경이 폐쇄 기준(0.30 mm)보다 크면
            <b>어떤 약물 노출로도 그 관은 닫히지 않습니다.</b> 24주에서는
            천장에서도 0.4 mm 이상이 남습니다. 임상에서 관찰되는 재태연령별
            치료 성공률의 급격한 기울기는 이 구조에서 <i>도출</i>된 것이며,
            여기에 맞춰 <i>적합</i>시킨 것이 아닙니다.</p>")),

        ## ---------------------------------------------------------- 2
        tabPanel("② 약물 PK · COX 억제",
          br(),
          fluidRow(column(6, plotOutput("plot_pk", height = 300)),
                   column(6, plotOutput("plot_cox", height = 300))),
          hr(),
          fluidRow(column(6, h4("노출 요약"), tableOutput("tbl_pk")),
                   column(6, h4("두 개의 효소 부위"),
                     HTML("<p>이부프로펜·인도메타신은 <b>아라키돈산 채널</b>에서
                     기질과 경쟁하고, 아세트아미노펜은 물리적으로 분리된
                     <b>peroxidase 부위</b>에서 ferryl-protoporphyrin 라디칼을
                     환원시킵니다. 따라서 아세트아미노펜은 기질이 아니라
                     <b>peroxide와 경쟁</b>하며, 관벽 peroxide 긴장도가 올라가면
                     겉보기 IC₅₀가 상승합니다. 사이드바에서 융모양막염을
                     체크하고 두 약을 비교해 보십시오 — 이 모델의 검증 가능한
                     예측입니다.</p>")))),

        ## ---------------------------------------------------------- 3
        tabPanel("③ 관 반응 (tone · 내경 · 재형성)",
          br(),
          fluidRow(column(6, plotOutput("plot_duct", height = 300)),
                   column(6, plotOutput("plot_wall", height = 300))),
          hr(),
          fluidRow(column(12, plotOutput("plot_drive", height = 260))),
          HTML("<p><b>수축은 폐쇄가 아닙니다.</b> 영구 폐쇄에는 관 <i>벽</i>의
            저산소증이 필요합니다 — 내경이 닫히면 내강으로부터의 산소 확산이
            끊기고, 중막이 저산소가 되어 HIF-1α/VEGF/TGF-β1이 신생내막 쿠션을
            만듭니다. 미숙아의 관벽은 <i>얇아서</i> 수축 중에도 확산으로 산소가
            공급되므로 재형성 신호가 켜지지 않고, 약물이 제거되면
            <b>재개통</b>합니다.</p>")),

        ## ---------------------------------------------------------- 4
        tabPanel("④ 혈역학 · 단락 부담",
          br(),
          fluidRow(column(6, plotOutput("plot_shunt", height = 300)),
                   column(6, plotOutput("plot_press", height = 300))),
          hr(),
          fluidRow(column(6, plotOutput("plot_organ", height = 280)),
                   column(6, plotOutput("plot_burden", height = 280))),
          HTML("<p>출생 직후에는 폐혈관저항(PVR)이 높아 큰 관도 단락이 적습니다.
            PVR이 2–5일에 걸쳐 떨어지면서 단락이 <b>드러납니다</b> — 이것이
            '3일째 큰 PDA'가 생기는 이유이고, 24시간 이내 예방투여가 아직
            문제가 되지 않은 관을 치료하게 되는 이유입니다.
            저항이 내경의 <b>4승</b>에 반비례하므로 1.5 mm는 임의의 절단값이
            아니라 실제 문턱입니다.</p>")),

        ## ---------------------------------------------------------- 5
        tabPanel("⑤ 임상 엔드포인트",
          br(),
          fluidRow(column(7, plotOutput("plot_outcome", height = 320)),
                   column(5, h4("36주 교정연령 시점"),
                          tableOutput("tbl_outcome"))),
          hr(),
          HTML("<p>결과 위험을 구동하는 것은 '7일째 폐쇄 여부'가 아니라
            <b>PDA 부담</b> — 유의한 단락의 시간 적분입니다. 이것이
            Baby-OSCAR에서 관은 더 잘 닫혔는데도 사망/BPD가 개선되지 않은
            (오히려 69.4% 대 63.5%) 결과를 모델이 재현하는 방식입니다.
            재태연령 자체가 BPD 위험을 지배하고, 약물이 줄일 수 있는 부담
            성분은 그에 비해 작습니다.</p>")),

        ## ---------------------------------------------------------- 6
        tabPanel("⑥ 시나리오 비교",
          br(),
          fluidRow(column(4, actionButton("run_scen", "16개 시나리오 실행",
                                          class = "btn-warning")),
                   column(8, helpText("계산에 다소 시간이 걸립니다."))),
          br(),
          if (have_DT) DT::dataTableOutput("tbl_scen") else tableOutput("tbl_scen_plain"),
          hr(),
          plotOutput("plot_scen", height = 340)),

        ## ---------------------------------------------------------- 7
        tabPanel("⑦ 바이오마커 · 안전성",
          br(),
          fluidRow(column(6, plotOutput("plot_renal", height = 300)),
                   column(6, plotOutput("plot_safety", height = 300))),
          hr(),
          fluidRow(column(12, tableOutput("tbl_safety"))),
          HTML("<p><b>효능과 독성은 다른 장기에서 일어나는 같은 분자 사건입니다.</b>
            COX는 이 모델에서 네 번 등장합니다 — 관, 신장, 장 점막, 혈소판.
            두 NSAID를 갈라놓는 것은 COX 효능이 아니라, 인도메타신에는 있고
            이부프로펜에는 없는 <b>비-COX 혈관수축</b>입니다(뇌혈류 −20~40%
            대 거의 변화 없음).</p>")),

        ## ---------------------------------------------------------- 8
        tabPanel("⑧ 자연 폐쇄 (경쟁 위험)",
          br(),
          fluidRow(column(4, actionButton("run_spont", "재태연령 스윕 실행",
                                          class = "btn-warning")),
                   column(8, helpText("무치료로 재태연령 24–38주를 시뮬레이션합니다."))),
          br(),
          fluidRow(column(6, plotOutput("plot_spont", height = 320)),
                   column(6, tableOutput("tbl_spont"))),
          hr(),
          HTML("<p>Semberova(2017)는 비개입 정책 하에서 27주 미만 미숙아의
            94–98%가 퇴원 전 자연 폐쇄되었고, 26주 미만에서 중앙값 71일이었다고
            보고했습니다. 거의 모든 관이 결국 닫힌다면, <b>고정 시점 폐쇄
            엔드포인트에 대한 약물 효과는 그 시점의 미폐쇄 분율로 상한이
            정해집니다.</b> 이 천장이 — 약효가 아니라 — 음성 시험들이 실제로
            측정한 것입니다.</p>")),

        ## ---------------------------------------------------------- 9
        tabPanel("⑨ 임상시험 재현",
          br(),
          h4("모델이 맞춘 것과 예측한 것"),
          tableOutput("tbl_trials"),
          hr(),
          HTML("<p><b>적합(fitted):</b> 8개 재태연령의 자연 폐쇄 시점 →
            NET50 · NETW · TAUSYN0 · KTAUSYN · KINVGA. 약물별 7일 폐쇄율 →
            Kᵢ·IC₅₀. Baby-OSCAR 1차 결과 양 군 → 위험 계수 3개.<br>
            <b>예측(predicted, 이후 어떤 파라미터도 조정하지 않음):</b>
            재개통률과 그 재태연령 의존성 · peroxide 상승 시 아세트아미노펜
            효력 소실 · 후기 구제요법의 고용량 필요성 · BeNeDuctus 복합
            결과 · TIPP의 IVH 감소와 결과 무개선 · 두 NSAID의 신장/뇌/장
            혈류 분리 · 표적 치료 전략 전체.</p>"))
      )
    )
  )
)

## =====================================================================
##  SERVER
## =====================================================================
server <- function(input, output, session) {

  build_rx <- reactive({
    s <- input$start_d * 24
    e <- NULL
    if (input$drug %in% c("ibu", "combo")) {
      e <- if (input$ibu_mode == "inf") rx_ibuprofen_inf(s)
           else rx_ibuprofen(s, high = (input$ibu_mode == "high"))
    }
    if (input$drug == "ind")
      e <- rx_indomethacin(s, dose = input$ind_dose,
                           maint = min(input$ind_dose, 0.1),
                           n = input$ind_n)
    if (input$drug %in% c("apap", "combo")) {
      a <- rx_acetaminophen(s, dose = input$apap_dose, days = input$apap_days)
      e <- if (is.null(e)) a else ev_bind(e, a)
    }
    if (!is.null(e) && isTRUE(input$second_course)) {
      s2 <- s + 5 * 24
      e2 <- if (input$drug == "ind")
              rx_indomethacin(s2, dose = input$ind_dose,
                              maint = min(input$ind_dose, 0.1), n = input$ind_n)
            else if (input$drug == "apap")
              rx_acetaminophen(s2, dose = input$apap_dose, days = input$apap_days)
            else rx_ibuprofen(s2, high = TRUE)
      e <- ev_bind(e, e2)
    }
    e
  })

  sim <- eventReactive(input$go, ignoreNULL = FALSE, {
    withProgress(message = "시뮬레이션 중...", value = 0.4, {
      as.data.frame(
        sim_pda(doses = build_rx(), days = input$horizon, delta = 1,
                GA = input$GA, BW = input$BW, PAO2 = input$PAO2,
                SEPSIS = as.numeric(input$SEPSIS),
                HCORT = as.numeric(input$HCORT),
                ANTESTER = as.numeric(input$ANTESTER)))
    })
  })

  ## reference: identical patient, no drug -- every "effect" is a difference
  sim_ref <- eventReactive(input$go, ignoreNULL = FALSE, {
    as.data.frame(
      sim_pda(doses = NULL, days = input$horizon, delta = 1,
              GA = input$GA, BW = input$BW, PAO2 = input$PAO2,
              SEPSIS = as.numeric(input$SEPSIS),
              HCORT = as.numeric(input$HCORT),
              ANTESTER = as.numeric(input$ANTESTER)))
  })

  cd <- function(d) { i <- which(d$DDUCT < 0.30); if (length(i)) d$DAY[i[1]] else NA }

  ## -------------------------------------------------------------- verdict
  output$verdict <- renderUI({
    d <- sim(); r <- sim_ref()
    c1 <- cd(d); c0 <- cd(r)
    reop <- { i <- which(d$DDUCT < 0.30)
              if (length(i)) any(d$DDUCT[i[1]:nrow(d)] > 0.48) else FALSE }
    b1 <- d$PDABUR[nrow(d)]; b0 <- r$PDABUR[nrow(r)]
    p1 <- 100 * d$PCOMP[nrow(d)]; p0 <- 100 * r$PCOMP[nrow(r)]
    HTML(sprintf(
      "<div style='font-size:12px;line-height:1.5'>
       <b>폐쇄일</b>: %s (무치료 %s)<br>
       <b>가속</b>: %s일<br>
       <b>재개통</b>: %s<br>
       <b>PDA 부담</b>: %.0f (무치료 %.0f, %+.0f%%)<br>
       <b>사망/BPD</b>: %.1f%% (무치료 %.1f%%)</div>",
      ifelse(is.na(c1), "미폐쇄", sprintf("d%.1f", c1)),
      ifelse(is.na(c0), "미폐쇄", sprintf("d%.1f", c0)),
      ifelse(is.na(c1) || is.na(c0), "—", sprintf("%.1f", c0 - c1)),
      ifelse(reop, "<span style='color:#c62828'>예</span>", "아니오"),
      b1, b0, ifelse(b0 > 0, 100 * (b1 - b0) / b0, 0), p1, p0))
  })

  ## -------------------------------------------------------------- tab 1
  output$profile_txt <- renderPrint({
    ga <- input$GA
    tmax <- 1 / (1 + exp(-(ga - 21.5) * 0.55))
    dmax <- max(0.6, -1.35 + 0.135 * ga)
    dmin <- dmax * (1 - 0.965 * tmax)
    p50  <- 30 * exp(0.130 * (28 - ga))
    wall <- max(0.05, -1.20 + 0.075 * ga)
    cat(sprintf("재태연령                       %.1f wk\n", ga))
    cat(sprintf("비수축 관 내경  d_max          %.2f mm\n", dmax))
    cat(sprintf("최대 수축 가능도 TMAXGA        %.3f\n", tmax))
    cat(sprintf("천장에서의 잔여 내경           %.3f mm\n", dmin))
    cat(sprintf("폐쇄 기준                      %.2f mm\n", 0.30))
    cat(sprintf("→ 약물로 폐쇄 가능?            %s\n",
                ifelse(dmin < 0.30, "가능", "불가능 (구조적 한계)")))
    cat(sprintf("\nO2 감수성 P50                  %.1f mmHg  (PaO2 %.0f)\n",
                p50, input$PAO2))
    cat(sprintf("관벽 두께 지수                 %.2f\n", wall))
    cat(sprintf("→ 수축 시 벽 저산소 가능?      %s\n",
                ifelse(wall > 0.75, "예 → 영구 재형성", "약함 → 재개통 위험")))
  })

  output$plot_ceiling <- renderPlot({
    ga <- seq(23, 40, 0.1)
    tmax <- 1 / (1 + exp(-(ga - 21.5) * 0.55))
    dmax <- pmax(0.6, -1.35 + 0.135 * ga)
    df <- data.frame(ga, resid = dmax * (1 - 0.965 * tmax))
    ggplot(df, aes(ga, resid)) +
      geom_hline(yintercept = 0.30, linetype = 2, colour = "#c62828") +
      geom_line(linewidth = 1.1, colour = "#3949ab") +
      geom_vline(xintercept = input$GA, colour = "#00796b", linetype = 3) +
      annotate("text", x = 36, y = 0.33, label = "폐쇄 기준 0.30 mm",
               colour = "#c62828", size = 3.4) +
      labs(x = "재태연령 (wk)", y = "최대 수축 시 잔여 내경 (mm)",
           title = "구조적 폐쇄 한계") + THEME
  })

  ## -------------------------------------------------------------- tab 2
  output$plot_pk <- renderPlot({
    d <- sim()
    df <- d %>% select(DAY, CIBU, CIND, CAPAP) %>%
      pivot_longer(-DAY) %>%
      mutate(name = recode(name, CIBU = "ibuprofen", CIND = "indomethacin",
                           CAPAP = "acetaminophen")) %>%
      filter(value > 1e-6)
    if (!nrow(df)) return(ggplot() + annotate("text", 1, 1, label = "무치료") +
                          THEME + theme_void())
    ggplot(df, aes(DAY, value, colour = name)) +
      geom_line(linewidth = 0.9) + scale_colour_manual(values = PAL) +
      coord_cartesian(xlim = c(0, min(30, input$horizon))) +
      labs(x = "생후 일수", y = "총 혈중 농도 (mg/L)", colour = NULL,
           title = "약물 농도") + THEME
  })

  output$plot_cox <- renderPlot({
    d <- sim()
    df <- d %>% select(DAY, ICOXD, ICHAN, IPEROX, ICOXK, ICOXG, ICOXPLT) %>%
      pivot_longer(-DAY) %>%
      mutate(name = factor(name,
        levels = c("ICOXD", "ICHAN", "IPEROX", "ICOXK", "ICOXG", "ICOXPLT"),
        labels = c("관 (net)", "채널 성분", "peroxidase 성분",
                   "신장", "장 점막", "혈소판")))
    ggplot(df, aes(DAY, 100 * value, colour = name)) +
      geom_line(linewidth = 0.85) +
      coord_cartesian(xlim = c(0, min(30, input$horizon))) +
      labs(x = "생후 일수", y = "COX 억제 (%)", colour = NULL,
           title = "하나의 효소, 네 개의 장기") + THEME
  })

  output$tbl_pk <- renderTable({
    d <- sim()
    f <- function(x, nm, u) {
      if (max(x) < 1e-6) return(NULL)
      data.frame(항목 = nm, Cmax = sprintf("%.2f", max(x)),
                 단위 = u, stringsAsFactors = FALSE)
    }
    rows <- rbind(f(d$CIBU, "이부프로펜 (총)", "mg/L"),
                  f(d$UIBU, "이부프로펜 (유리)", "uM"),
                  f(d$CIND, "인도메타신 (총)", "mg/L"),
                  f(d$UIND, "인도메타신 (유리)", "uM"),
                  f(d$CAPAP, "아세트아미노펜 (총)", "mg/L"),
                  f(d$UAPAP, "아세트아미노펜 (유리)", "uM"))
    if (is.null(rows)) rows <- data.frame(항목 = "무치료", Cmax = "-", 단위 = "-")
    rbind(rows, data.frame(항목 = "최대 관 COX 억제",
                           Cmax = sprintf("%.1f", 100 * max(d$ICOXD)),
                           단위 = "%"))
  }, digits = 3)

  ## -------------------------------------------------------------- tab 3
  output$plot_duct <- renderPlot({
    d <- sim(); r <- sim_ref()
    df <- rbind(data.frame(DAY = d$DAY, v = d$DDUCT, arm = "치료"),
                data.frame(DAY = r$DAY, v = r$DDUCT, arm = "무치료"))
    ggplot(df, aes(DAY, v, colour = arm)) +
      geom_hline(yintercept = 0.30, linetype = 2, colour = "#c62828") +
      geom_line(linewidth = 1) +
      scale_colour_manual(values = c("치료" = "#00796b", "무치료" = "#90a4ae")) +
      labs(x = "생후 일수", y = "관 내경 (mm)", colour = NULL,
           title = "관 내경 — 치료 대 자연 경과") + THEME
  })

  output$plot_wall <- renderPlot({
    d <- sim()
    df <- d %>% select(DAY, TONE, REMOD, WALLO2) %>%
      mutate(WALLO2 = WALLO2 / 60) %>% pivot_longer(-DAY) %>%
      mutate(name = recode(name, TONE = "관 tone (0-1)",
                           REMOD = "신생내막 쿠션 (0-1)",
                           WALLO2 = "관벽 PO2 / 60 mmHg"))
    ggplot(df, aes(DAY, value, colour = name)) +
      geom_line(linewidth = 0.95) +
      coord_cartesian(xlim = c(0, min(40, input$horizon))) +
      labs(x = "생후 일수", y = NULL, colour = NULL,
           title = "수축 · 벽 산소 · 영구 재형성") + THEME
  })

  output$plot_drive <- renderPlot({
    d <- sim()
    df <- d %>% select(DAY, RELAXP, RELAXN, O2GAIN, NETD, TONETGT) %>%
      pivot_longer(-DAY) %>%
      mutate(name = recode(name, RELAXP = "PGE2/EP4 이완",
                           RELAXN = "NO 이완", O2GAIN = "O2 수축 gain",
                           NETD = "순 수축 구동 (net drive)",
                           TONETGT = "목표 tone"))
    ggplot(df, aes(DAY, value, colour = name)) +
      geom_line(linewidth = 0.85) +
      coord_cartesian(xlim = c(0, min(40, input$horizon))) +
      labs(x = "생후 일수", y = NULL, colour = NULL,
           title = "tone 은 gain 의 곱이 아니라 순 구동의 시그모이드") + THEME
  })

  ## -------------------------------------------------------------- tab 4
  output$plot_shunt <- renderPlot({
    d <- sim(); r <- sim_ref()
    df <- rbind(data.frame(DAY = d$DAY, v = d$QSH, arm = "치료"),
                data.frame(DAY = r$DAY, v = r$QSH, arm = "무치료"))
    ggplot(df, aes(DAY, v, colour = arm)) +
      geom_hline(yintercept = 100, linetype = 2, colour = "#c62828") +
      geom_line(linewidth = 1) +
      scale_colour_manual(values = c("치료" = "#c62828", "무치료" = "#90a4ae")) +
      coord_cartesian(xlim = c(0, min(40, input$horizon))) +
      annotate("text", x = 30, y = 112, size = 3.2, colour = "#c62828",
               label = "유의성 문턱") +
      labs(x = "생후 일수", y = "좌-우 단락 (mL/min/kg)", colour = NULL,
           title = "단락 — PVR 이 떨어지며 드러난다") + THEME
  })

  output$plot_press <- renderPlot({
    d <- sim()
    df <- d %>% select(DAY, PAOM, PPA, PDIA, QPQS) %>%
      mutate(QPQS = QPQS * 10) %>% pivot_longer(-DAY) %>%
      mutate(name = recode(name, PAOM = "평균 대동맥압 (mmHg)",
                           PPA = "평균 폐동맥압 (mmHg)",
                           PDIA = "확장기압 (mmHg)",
                           QPQS = "Qp:Qs x 10"))
    ggplot(df, aes(DAY, value, colour = name)) + geom_line(linewidth = 0.85) +
      coord_cartesian(xlim = c(0, min(40, input$horizon))) +
      labs(x = "생후 일수", y = NULL, colour = NULL, title = "압력과 유량비") +
      THEME
  })

  output$plot_organ <- renderPlot({
    d <- sim()
    df <- d %>% select(DAY, QCERREL, QMESREL, QRENREL) %>%
      pivot_longer(-DAY) %>%
      mutate(name = recode(name, QCERREL = "뇌", QMESREL = "장",
                           QRENREL = "신장"))
    ggplot(df, aes(DAY, 100 * value, colour = name)) +
      geom_line(linewidth = 0.9) +
      coord_cartesian(xlim = c(0, min(20, input$horizon))) +
      labs(x = "생후 일수", y = "기저치 대비 혈류 (%)", colour = NULL,
           title = "확장기 유출 + 약물 혈관수축") + THEME
  })

  output$plot_burden <- renderPlot({
    d <- sim(); r <- sim_ref()
    df <- rbind(data.frame(DAY = d$DAY, v = d$PDABUR, arm = "치료"),
                data.frame(DAY = r$DAY, v = r$PDABUR, arm = "무치료"))
    ggplot(df, aes(DAY, v, colour = arm)) + geom_line(linewidth = 1) +
      scale_colour_manual(values = c("치료" = "#37474f", "무치료" = "#90a4ae")) +
      labs(x = "생후 일수", y = "PDA 부담 (누적)", colour = NULL,
           title = "결과가 보는 변수는 폐쇄가 아니라 부담") + THEME
  })

  ## -------------------------------------------------------------- tab 5
  output$plot_outcome <- renderPlot({
    d <- sim(); r <- sim_ref()
    mk <- function(x, arm) x %>% select(DAY, PBPD, PNEC, PIVH, PDTH, PCOMP) %>%
      pivot_longer(-DAY) %>% mutate(arm = arm)
    df <- rbind(mk(d, "치료"), mk(r, "무치료")) %>%
      mutate(name = recode(name, PBPD = "중등도-중증 BPD", PNEC = "NEC >= 2기",
                           PIVH = "중증 IVH", PDTH = "사망",
                           PCOMP = "사망 또는 BPD"))
    ggplot(df, aes(DAY, 100 * value, colour = name, linetype = arm)) +
      geom_line(linewidth = 0.85) +
      labs(x = "생후 일수", y = "누적 확률 (%)", colour = NULL, linetype = NULL,
           title = "결과 위험") + THEME
  })

  output$tbl_outcome <- renderTable({
    d <- sim(); r <- sim_ref()
    g <- function(x, nm) x[[nm]][nrow(x)] * 100
    data.frame(
      엔드포인트 = c("사망 또는 중등도-중증 BPD", "중등도-중증 BPD", "사망",
                     "NEC >= 2기", "중증 IVH", "장 자연천공", "PDA 부담"),
      치료 = c(g(d, "PCOMP"), g(d, "PBPD"), g(d, "PDTH"), g(d, "PNEC"),
               g(d, "PIVH"), g(d, "PSIP"), d$PDABUR[nrow(d)]),
      무치료 = c(g(r, "PCOMP"), g(r, "PBPD"), g(r, "PDTH"), g(r, "PNEC"),
                 g(r, "PIVH"), g(r, "PSIP"), r$PDABUR[nrow(r)]),
      stringsAsFactors = FALSE)
  }, digits = 1)

  ## -------------------------------------------------------------- tab 6
  scen <- eventReactive(input$run_scen, {
    withProgress(message = "16개 시나리오 실행 중...", value = 0.2,
                 run_scenarios(days = 90))
  })
  if (have_DT) {
    output$tbl_scen <- DT::renderDataTable({
      DT::datatable(scen(), options = list(pageLength = 16, dom = "t",
                                          scrollX = TRUE), rownames = FALSE)
    })
  } else {
    output$tbl_scen_plain <- renderTable({ scen() }, digits = 2)
  }

  output$plot_scen <- renderPlot({
    s <- scen()
    s$scenario <- factor(s$id, levels = rev(s$id))
    ggplot(s, aes(P_comp, scenario)) +
      geom_col(aes(fill = burden), width = 0.72) +
      scale_fill_gradient(low = "#b2dfdb", high = "#c62828") +
      labs(x = "사망 또는 중등도-중증 BPD (%)", y = NULL, fill = "PDA 부담",
           title = "부담이 낮은 시나리오가 결과도 좋다 — 폐쇄 자체가 아니라") +
      THEME
  })

  ## -------------------------------------------------------------- tab 7
  output$plot_renal <- renderPlot({
    d <- sim(); r <- sim_ref()
    mk <- function(x, arm) data.frame(DAY = x$DAY, SCr = x$SCR, UO = x$UO,
                                      GFR = x$GFR, arm = arm)
    df <- rbind(mk(d, "치료"), mk(r, "무치료")) %>% pivot_longer(-c(DAY, arm))
    ggplot(df, aes(DAY, value, colour = name, linetype = arm)) +
      geom_line(linewidth = 0.9) +
      coord_cartesian(xlim = c(0, min(20, input$horizon))) +
      labs(x = "생후 일수", y = NULL, colour = NULL, linetype = NULL,
           title = "신장: SCr (mg/dL) · 요량 (mL/kg/h) · GFR") + THEME
  })

  output$plot_safety <- renderPlot({
    d <- sim()
    df <- d %>% select(DAY, ALT, TXA2, BTIME, BFREE) %>%
      mutate(ALT = ALT / 20) %>% pivot_longer(-DAY) %>%
      mutate(name = recode(name, ALT = "ALT / 20 U/L",
                           TXA2 = "혈소판 TXA2 (분율)",
                           BTIME = "출혈시간 지수",
                           BFREE = "유리 빌리루빈 지수"))
    ggplot(df, aes(DAY, value, colour = name)) + geom_line(linewidth = 0.9) +
      coord_cartesian(xlim = c(0, min(20, input$horizon))) +
      labs(x = "생후 일수", y = NULL, colour = NULL, title = "안전성 지표") +
      THEME
  })

  output$tbl_safety <- renderTable({
    d <- sim(); r <- sim_ref()
    e <- d[d$DAY <= 12, ]; e0 <- r[r$DAY <= 12, ]
    data.frame(
      지표 = c("최고 SCr (mg/dL)", "요량 최저 (mL/kg/h)", "GFR 최저 (mL/min/kg)",
               "뇌혈류 최저 (%)", "장혈류 최저 (%)", "최고 ALT (U/L)",
               "출혈시간 지수 최고", "유리 빌리루빈 지수 최고"),
      치료 = c(max(e$SCR), min(e$UO), min(e$GFR), 100 * min(e$QCERREL),
               100 * min(e$QMESREL), max(e$ALT), max(e$BTIME), max(e$BFREE)),
      무치료 = c(max(e0$SCR), min(e0$UO), min(e0$GFR), 100 * min(e0$QCERREL),
                 100 * min(e0$QMESREL), max(e0$ALT), max(e0$BTIME),
                 max(e0$BFREE)), stringsAsFactors = FALSE)
  }, digits = 2)

  ## -------------------------------------------------------------- tab 8
  spont <- eventReactive(input$run_spont, {
    withProgress(message = "재태연령 스윕...", value = 0.2,
                 spontaneous_closure_table())
  })
  output$tbl_spont <- renderTable({ spont() }, digits = 2)
  output$plot_spont <- renderPlot({
    s <- spont()
    tgt <- data.frame(GA = c(24, 25, 26, 27, 28, 30, 32, 38),
                      obs = c(85, 71, 48, 32, 21, 10, 5, 1.5))
    ggplot(s, aes(GA, closure_day)) +
      geom_line(linewidth = 1.1, colour = "#3949ab") +
      geom_point(size = 2, colour = "#3949ab") +
      geom_point(data = tgt, aes(GA, obs), colour = "#c62828",
                 shape = 4, size = 3, stroke = 1.2) +
      scale_y_log10() +
      labs(x = "재태연령 (wk)", y = "자연 폐쇄일 (log)",
           title = "모델(선) 대 관찰 목표(×)") + THEME
  })

  ## -------------------------------------------------------------- tab 9
  output$tbl_trials <- renderTable({
    data.frame(
      시험 = c("Semberova 2017", "Baby-OSCAR 2024", "BeNeDuctus 2023",
               "TIPP 2001", "PDA-TOLERATE 2019", "Cochrane / NMA"),
      관찰 = c("<27주 94-98% 자연폐쇄, <26주 중앙값 71일",
               "조기 이부프로펜 대 기대: 사망/BPD 69.4% 대 63.5%",
               "기대요법 비열등: 복합 46.5% 대 63.5%",
               "예방적 인도메타신: 중증 IVH 9% 대 13%, 18개월 결과 무변화",
               "6-14일 조기 치료의 이득 없음",
               "폐쇄율 3약제 유사; 이부프로펜이 NEC·AKI 적음"),
      역할 = c("적합 (involution 파라미터)", "적합 (위험 계수 3개)",
               "예측", "예측", "예측", "적합 (Ki/IC50) + 예측 (독성 분리)"),
      stringsAsFactors = FALSE)
  })
}

shinyApp(ui, server)
