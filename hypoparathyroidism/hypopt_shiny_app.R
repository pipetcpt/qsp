## ============================================================================
##  Chronic Hypoparathyroidism (HypoPT) — QSP Interactive Dashboard
##  ============================================================================
##  Front-end for hypopt_mrgsolve_model.R.
##
##  Run:
##      shiny::runApp("hypopt_shiny_app.R")
##  or from this directory:
##      R -e 'shiny::runApp(".", launch.browser = TRUE)'
##
##  Nine tabs:
##      1. 환자 프로파일   — aetiology, chief-cell mass, CaSR set-point, labs
##      2. PK              — PTH molecules, calcitriol, thiazide, encaleret
##      3. 칼슘·인 항상성  — the homeostatic core (Ca, iCa, Pi, Mg, 1,25D, FGF23)
##      4. 신장 안전성     — urine Ca, FE_Ca, TmP/GFR, nephrocalcinosis, eGFR
##      5. 뼈              — turnover markers, osteoblast/osteoclast, BMD Z
##      6. 임상 엔드포인트 — symptoms, QTc, QoL, triple composite response
##      7. 시나리오 비교   — the 15 prebuilt scenarios side by side
##      8. 바이오마커 요약 — one table of every steady-state readout
##      9. 문헌·도움말     — model provenance and how to read it
##
##  DISCLAIMER: educational / research tool. Not for clinical use.
## ============================================================================

suppressPackageStartupMessages({
  library(shiny)
  library(mrgsolve)
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(DT)
})

source("hypopt_mrgsolve_model.R", local = TRUE, chdir = TRUE)

MOD <- HYPOPT_build()

## ---------------------------------------------------------------------------
##  Plot furniture
## ---------------------------------------------------------------------------
theme_hp <- function() {
  theme_minimal(base_size = 12) +
    theme(panel.grid.minor = element_blank(),
          strip.text  = element_text(face = "bold"),
          legend.position = "bottom",
          plot.title  = element_text(face = "bold"))
}

PAL <- c("#2471A3", "#C0392B", "#1E8449", "#B9770E", "#7D3C98",
         "#117A65", "#CA6F1E", "#34495E", "#CB4335", "#5D6D7E")

## reference bands drawn on the calcium panels
band <- function(ymin, ymax, fill = "#1E8449", alpha = 0.08) {
  annotate("rect", xmin = -Inf, xmax = Inf, ymin = ymin, ymax = ymax,
           fill = fill, alpha = alpha)
}

lineplot <- function(d, y, ylab, title, xdays = TRUE, bands = NULL) {
  d$xx <- if (xdays) d$time/24 else d$time
  p <- ggplot(d, aes(x = xx, y = .data[[y]], colour = scenario))
  if (!is.null(bands)) p <- p + band(bands[1], bands[2])
  p + geom_line(linewidth = 0.7) +
    scale_colour_manual(values = rep(PAL, 4)) +
    labs(x = if (xdays) "시간 (일)" else "시간 (시간)", y = ylab,
         title = title, colour = NULL) +
    theme_hp()
}

## ---------------------------------------------------------------------------
##  UI
## ---------------------------------------------------------------------------
ui <- fluidPage(
  titlePanel("만성 부갑상선기능저하증 (Chronic Hypoparathyroidism) — QSP 시뮬레이터"),
  tags$p(style = "color:#7F8C8D;margin-top:-8px;",
         "교육·연구용 정량적 시스템 약리학 모델. 임상 의사결정에 사용하지 마십시오."),

  sidebarLayout(
    sidebarPanel(
      width = 3,

      h4("환자 · 질환"),
      selectInput("aetiology", "병인 (Aetiology)",
                  choices = c("수술 후 (Post-surgical)"        = "surg",
                              "자가면역 (Autoimmune)"          = "ai",
                              "ADH1 (CASR 기능획득)"           = "adh1",
                              "저마그네슘혈증 (기능성)"        = "mglow",
                              "정상 대조군"                    = "healthy"),
                  selected = "surg"),
      sliderInput("ptmass", "잔존 주세포 기능 (0-1)", 0, 1, 0.05, step = 0.01),
      sliderInput("setpt", "CaSR 설정점 (mmol/L)", 0.70, 1.40, 1.20, step = 0.01),
      sliderInput("alb", "혈청 알부민 (g/dL)", 2.0, 5.0, 4.0, step = 0.1),
      sliderInput("ph", "동맥혈 pH", 7.25, 7.60, 7.40, step = 0.01),
      sliderInput("gfr0", "기저 eGFR (mL/min/1.73m²)", 20, 130, 120, step = 5),
      sliderInput("dietca", "식이 칼슘 (mg/일)", 300, 1600, 1000, step = 50),
      sliderInput("dietpi", "식이 인 (mg/일)", 400, 2000, 1200, step = 50),
      sliderInput("dietmg", "식이 마그네슘 (mg/일)", 20, 500, 300, step = 10),

      hr(),
      h4("기존 치료 (Conventional)"),
      sliderInput("casup", "경구 칼슘 1회 용량 (mg 원소칼슘)", 0, 1500, 500, step = 50),
      selectInput("caii", "칼슘 투여 간격", c("1일 1회" = 24, "1일 2회" = 12,
                                              "1일 3회" = 8, "1일 4회" = 6),
                  selected = 8),
      sliderInput("ctr", "칼시트리올 1회 용량 (µg)", 0, 2, 0.25, step = 0.125),
      selectInput("ctrii", "칼시트리올 간격", c("1일 1회" = 24, "1일 2회" = 12),
                  selected = 12),
      sliderInput("d3", "콜레칼시페롤 (IU/일)", 0, 4000, 1000, step = 200),
      checkboxInput("thz", "티아지드 (HCTZ 25 mg BID)", FALSE),
      checkboxInput("binder", "인결합제 / 저인 식이", FALSE),
      sliderInput("mgsup", "마그네슘 보충 (mg/일)", 0, 1200, 0, step = 50),

      hr(),
      h4("PTH 대체요법"),
      selectInput("pthdrug", "약물",
                  choices = c("없음"                                = "none",
                              "rhPTH(1-84) SC QD"                   = "rhpth84",
                              "테리파라타이드 PTH(1-34) SC"          = "terip",
                              "팔로페그테리파라타이드 (TransCon PTH)" = "transcon",
                              "에네보파라타이드 (RG 선택적)"          = "eneb"),
                  selected = "none"),
      sliderInput("pthdose", "1회 용량 (µg)", 0, 120, 18, step = 2),
      selectInput("pthii", "투여 간격", c("1일 1회" = 24, "1일 2회" = 12,
                                          "1일 3회" = 8), selected = 24),

      hr(),
      h4("CaSR 표적치료"),
      sliderInput("enc", "엔칼레렛 1회 용량 (mg, BID)", 0, 120, 0, step = 10),

      hr(),
      sliderInput("days", "시뮬레이션 기간 (일)", 14, 3650, 180, step = 14),
      actionButton("go", "시뮬레이션 실행", class = "btn-primary btn-block")
    ),

    mainPanel(
      width = 9,
      tabsetPanel(
        id = "tabs",

        tabPanel("1. 환자 프로파일",
                 br(),
                 fluidRow(
                   column(6, h4("기저 (치료 전) 정상상태"), tableOutput("tbl_base")),
                   column(6, h4("치료 후 정상상태"),        tableOutput("tbl_treat"))),
                 hr(),
                 h4("칼슘 설정점 곡선 — 이 환자의 Ca–PTH 관계"),
                 plotOutput("p_setpoint", height = "330px"),
                 helpText("점선은 정상 설정점(1.20 mmol/L), 실선은 현재 설정과",
                          "잔존 주세포 기능을 반영한 곡선입니다.",
                          "ADH1에서는 곡선이 왼쪽으로 이동하여 낮은 칼슘에서도",
                          "PTH가 '부적절하게 정상'으로 보입니다.")),

        tabPanel("2. PK",
                 br(),
                 fluidRow(
                   column(6, plotOutput("p_pth", height = "300px")),
                   column(6, plotOutput("p_pthzoom", height = "300px"))),
                 fluidRow(
                   column(6, plotOutput("p_d125pk", height = "300px")),
                   column(6, plotOutput("p_drugpk", height = "300px"))),
                 helpText("오른쪽 위 그래프는 마지막 7일의 확대입니다.",
                          "rhPTH(1-84)의 이상성 흡수는 뾰족한 곡선을,",
                          "TransCon PTH의 서방 방출은 거의 평탄한 곡선을 만듭니다.")),

        tabPanel("3. 칼슘·인 항상성",
                 br(),
                 fluidRow(
                   column(6, plotOutput("p_ca",   height = "290px")),
                   column(6, plotOutput("p_ica",  height = "290px"))),
                 fluidRow(
                   column(6, plotOutput("p_pi",   height = "290px")),
                   column(6, plotOutput("p_mg",   height = "290px"))),
                 fluidRow(
                   column(6, plotOutput("p_d125", height = "290px")),
                   column(6, plotOutput("p_fgf",  height = "290px")))),

        tabPanel("4. 신장 안전성",
                 br(),
                 fluidRow(
                   column(6, plotOutput("p_uca",  height = "290px")),
                   column(6, plotOutput("p_feca", height = "290px"))),
                 fluidRow(
                   column(6, plotOutput("p_tmp",  height = "290px")),
                   column(6, plotOutput("p_caxp", height = "290px"))),
                 fluidRow(
                   column(6, plotOutput("p_nc",   height = "290px")),
                   column(6, plotOutput("p_egfr", height = "290px"))),
                 helpText("이 탭이 이 모델의 핵심입니다. 기존 치료는 혈청 칼슘을",
                          "목표에 넣지만 요중 칼슘은 낮추지 못합니다.",
                          "PTH 대체요법만이 같은 혈청 칼슘에서 요중 칼슘을 낮춥니다.")),

        tabPanel("5. 뼈",
                 br(),
                 fluidRow(
                   column(6, plotOutput("p_turn", height = "290px")),
                   column(6, plotOutput("p_cells", height = "290px"))),
                 fluidRow(
                   column(6, plotOutput("p_markers", height = "290px")),
                   column(6, plotOutput("p_bmd",     height = "290px")))),

        tabPanel("6. 임상 엔드포인트",
                 br(),
                 fluidRow(
                   column(6, plotOutput("p_sym", height = "290px")),
                   column(6, plotOutput("p_qtc", height = "290px"))),
                 fluidRow(
                   column(6, plotOutput("p_qol", height = "290px")),
                   column(6, plotOutput("p_bg",  height = "290px"))),
                 hr(),
                 h4("삼중 복합 반응 (PaTHway/PARALLAX 정의)"),
                 verbatimTextOutput("txt_resp3"),
                 helpText("알부민 보정 칼슘 8.0–10.6 mg/dL AND 경구 칼슘 ≤600 mg/일",
                          "AND 활성 비타민 D 중단 — 세 조건을 모두 만족해야 합니다.")),

        tabPanel("7. 시나리오 비교",
                 br(),
                 checkboxGroupInput(
                   "scn", "비교할 시나리오",
                   choices = c(
                     "01 정상 대조군" = 1,  "02 무치료" = 2,
                     "03 기존 치료" = 3,    "04 기존+티아지드" = 4,
                     "05 과치료" = 5,       "06 rhPTH(1-84)" = 6,
                     "07 테리파라타이드" = 7, "08 TransCon 18 µg" = 8,
                     "09 TransCon 30 µg" = 9, "10 ADH1 무치료" = 10,
                     "11 ADH1+엔칼레렛" = 11, "12 저Mg → 보충" = 12,
                     "13 급성 위기 (IV Ca)" = 13,
                     "14 10년 기존치료" = 14, "15 10년 TransCon" = 15),
                   selected = c(2, 3, 8), inline = TRUE),
                 actionButton("goscn", "시나리오 실행", class = "btn-success"),
                 br(), br(),
                 fluidRow(
                   column(6, plotOutput("s_ca",  height = "300px")),
                   column(6, plotOutput("s_uca", height = "300px"))),
                 fluidRow(
                   column(6, plotOutput("s_pi",   height = "300px")),
                   column(6, plotOutput("s_egfr", height = "300px"))),
                 hr(),
                 h4("시나리오 요약"),
                 DT::dataTableOutput("s_tbl"),
                 hr(),
                 h4("투여간격 내 혈청 칼슘 변동폭 (Ca swing)"),
                 DT::dataTableOutput("s_swing")),

        tabPanel("8. 바이오마커 요약",
                 br(),
                 h4("현재 설정의 정상상태 바이오마커"),
                 DT::dataTableOutput("b_tbl"),
                 hr(),
                 h4("기존 치료 용량 격자 — 목표 칼슘 달성의 신장 비용"),
                 helpText("칼슘/칼시트리올 조합별로 도달 혈청 칼슘과 그 대가로",
                          "발생하는 24시간 요중 칼슘을 계산합니다 (120일 정상상태)."),
                 actionButton("gogrid", "격자 계산 (수십 초 소요)"),
                 br(), br(),
                 DT::dataTableOutput("b_grid")),

        tabPanel("9. 문헌 · 도움말",
                 br(),
                 uiOutput("refs"))
      )
    )
  )
)

## ---------------------------------------------------------------------------
##  Server
## ---------------------------------------------------------------------------
server <- function(input, output, session) {

  ## keep the aetiology preset and the sliders in sync
  observeEvent(input$aetiology, {
    p <- switch(input$aetiology,
                surg    = list(pt = 0.05, sp = 1.20, mg = 300),
                ai      = list(pt = 0.10, sp = 1.20, mg = 300),
                adh1    = list(pt = 1.00, sp = 0.85, mg = 300),
                mglow   = list(pt = 1.00, sp = 1.20, mg = 40),
                healthy = list(pt = 1.00, sp = 1.20, mg = 300))
    updateSliderInput(session, "ptmass", value = p$pt)
    updateSliderInput(session, "setpt",  value = p$sp)
    updateSliderInput(session, "dietmg", value = p$mg)
  })

  ## ---- parameter block for the current UI state --------------------------
  base_pars <- reactive({
    list(PTMASS = input$ptmass, SETPT = input$setpt,
         ALB = input$alb, PHART = input$ph, GFRC0 = input$gfr0,
         DIETCA = input$dietca, DIETPI = input$dietpi, DIETMG = input$dietmg,
         MGSUP = 0, BINDER = 0)
  })

  treat_pars <- reactive({
    p <- base_pars()
    p$MGSUP  <- input$mgsup
    p$BINDER <- as.numeric(input$binder)
    p$DCASUP <- input$casup*(24/as.numeric(input$caii))
    p$DCTR   <- input$ctr*(24/as.numeric(input$ctrii))
    if (input$pthdrug != "none")
      p <- utils::modifyList(p, HYPOPT_pth_params[[input$pthdrug]])
    p
  })

  events <- reactive({
    d  <- input$days
    ev <- NULL
    add <- function(e) if (is.null(ev)) e else c(ev, e)
    if (input$casup > 0)
      ev <- add(hypopt_calcium(input$casup, ii = as.numeric(input$caii), days = d))
    if (input$ctr > 0)
      ev <- add(hypopt_calcitriol(input$ctr, ii = as.numeric(input$ctrii), days = d))
    if (input$d3 > 0)
      ev <- add(hypopt_vitd3(input$d3, ii = 24, days = d))
    if (isTRUE(input$thz))
      ev <- add(hypopt_thiazide(25, ii = 12, days = d))
    if (input$enc > 0)
      ev <- add(hypopt_encaleret(input$enc, ii = 12, days = d))
    if (input$pthdrug == "transcon" && input$pthdose > 0)
      ev <- add(hypopt_transcon(input$pthdose, ii = as.numeric(input$pthii), days = d))
    if (input$pthdrug %in% c("rhpth84", "terip", "eneb") && input$pthdose > 0)
      ev <- add(hypopt_pth_sc(input$pthdose, ii = as.numeric(input$pthii), days = d))
    ev
  })

  ## ---- the two runs: untreated baseline, then therapy --------------------
  baseline_state <- eventReactive(input$go, {
    withProgress(message = "질환 정상상태 계산 중...", value = 0.3, {
      hypopt_baseline(MOD, base_pars(), days = 400)
    })
  }, ignoreNULL = FALSE)

  sim <- eventReactive(input$go, {
    withProgress(message = "시뮬레이션 실행 중...", value = 0.6, {
      init <- hypopt_baseline(MOD, base_pars(), days = 400)
      m    <- update(MOD, param = treat_pars(), init = as.list(init))
      d    <- input$days
      delta <- if (d <= 60) 0.5 else if (d <= 400) 2 else 24*7
      ev   <- events()
      out  <- if (is.null(ev)) {
        mrgsim_df(m, end = d*24, delta = delta, hmax = 2)
      } else {
        mrgsim_df(m, events = ev, end = d*24, delta = delta, hmax = 2)
      }
      out$scenario <- "현재 설정"
      out
    })
  }, ignoreNULL = FALSE)

  ss <- function(d, col, frac = 0.1) {
    n <- nrow(d); k <- max(1, floor(n*(1 - frac)))
    mean(d[[col]][k:n], na.rm = TRUE)
  }

  ## ---- Tab 1 -------------------------------------------------------------
  output$tbl_base <- renderTable({
    iv <- baseline_state()
    m  <- update(MOD, param = base_pars(), init = as.list(iv))
    d  <- mrgsim_df(m, end = 48, delta = 2)
    data.frame(
      항목 = c("알부민 보정 칼슘 (mg/dL)", "이온화 칼슘 (mmol/L)",
               "혈청 인 (mg/dL)", "PTH (pg/mL)", "1,25(OH)₂D (pg/mL)",
               "24h 요중 칼슘 (mg)", "골회전 지수", "BMD Z-score"),
      값 = round(c(ss(d, "CACORRD"), ss(d, "CAIONMM"), ss(d, "PISER"),
                   ss(d, "PTHPG"), ss(d, "D125PG"), ss(d, "UCA24"),
                   ss(d, "TURNOV"), ss(d, "BMDZS")), 2))
  }, striped = TRUE)

  output$tbl_treat <- renderTable({
    d <- sim()
    data.frame(
      항목 = c("알부민 보정 칼슘 (mg/dL)", "이온화 칼슘 (mmol/L)",
               "혈청 인 (mg/dL)", "PTH 신호 (pg/mL eq)", "1,25(OH)₂D (pg/mL)",
               "24h 요중 칼슘 (mg)", "골회전 지수", "BMD Z-score",
               "eGFR (mL/min/1.73m²)", "증상 점수 (0-100)"),
      값 = round(c(ss(d, "CACORRD"), ss(d, "CAIONMM"), ss(d, "PISER"),
                   ss(d, "PTHPG"), ss(d, "D125PG"), ss(d, "UCA24"),
                   ss(d, "TURNOV"), ss(d, "BMDZS"), ss(d, "EGFR"),
                   ss(d, "SYMSCOR")), 2))
  }, striped = TRUE)

  output$p_setpoint <- renderPlot({
    ica  <- seq(0.60, 1.60, by = 0.005)
    curve <- function(sp, mass) {
      sec <- 0.10 + 0.90/(1 + (ica/sp)^6)
      mass*934*sec*0.823*0.835/10.4
    }
    df <- rbind(
      data.frame(ica = ica, pth = curve(1.20, 1.00), grp = "정상 (설정점 1.20)"),
      data.frame(ica = ica, pth = curve(input$setpt, input$ptmass),
                 grp = sprintf("환자 (설정점 %.2f, 주세포 %.2f)",
                               input$setpt, input$ptmass)))
    ggplot(df, aes(ica, pth, colour = grp, linetype = grp)) +
      geom_line(linewidth = 0.9) +
      annotate("rect", xmin = 1.15, xmax = 1.30, ymin = -Inf, ymax = Inf,
               fill = "#1E8449", alpha = 0.08) +
      scale_colour_manual(values = c("#5D6D7E", "#C0392B")) +
      scale_linetype_manual(values = c("dashed", "solid")) +
      labs(x = "이온화 칼슘 (mmol/L)", y = "PTH (pg/mL)",
           title = "Ca–PTH 설정점 곡선", colour = NULL, linetype = NULL) +
      theme_hp()
  })

  ## ---- Tab 2 : PK --------------------------------------------------------
  output$p_pth <- renderPlot(
    lineplot(sim(), "PTHPG", "PTH 신호 (pg/mL eq)", "총 생물학적 PTH 노출"))
  output$p_pthzoom <- renderPlot({
    d <- sim(); d <- d[d$time >= max(0, max(d$time) - 7*24), ]
    lineplot(d, "PTHPG", "PTH 신호 (pg/mL eq)", "마지막 7일 확대")
  })
  output$p_d125pk <- renderPlot(
    lineplot(sim(), "D125PG", "1,25(OH)₂D (pg/mL)", "칼시트리올 노출"))
  output$p_drugpk <- renderPlot({
    d <- sim() %>% select(time, PTHEXOG, D125PG, D25NG) %>%
      pivot_longer(-time)
    ggplot(d, aes(time/24, value, colour = name)) +
      geom_line(linewidth = 0.7) +
      facet_wrap(~name, scales = "free_y", ncol = 1) +
      scale_colour_manual(values = PAL) +
      labs(x = "시간 (일)", y = NULL, title = "약물·대사물 농도", colour = NULL) +
      theme_hp() + theme(legend.position = "none")
  })

  ## ---- Tab 3 : homeostasis ----------------------------------------------
  output$p_ca <- renderPlot(
    lineplot(sim(), "CACORRD", "보정 칼슘 (mg/dL)",
             "알부민 보정 혈청 칼슘", bands = c(8.0, 9.0)))
  output$p_ica <- renderPlot(
    lineplot(sim(), "CAIONMM", "이온화 칼슘 (mmol/L)",
             "이온화 칼슘", bands = c(1.15, 1.30)))
  output$p_pi <- renderPlot(
    lineplot(sim(), "PISER", "인 (mg/dL)", "혈청 인", bands = c(2.5, 4.5)))
  output$p_mg <- renderPlot(
    lineplot(sim(), "MGSER", "마그네슘 (mg/dL)", "혈청 마그네슘",
             bands = c(1.7, 2.4)))
  output$p_d125 <- renderPlot(
    lineplot(sim(), "D125PG", "1,25(OH)₂D (pg/mL)", "활성 비타민 D",
             bands = c(20, 60)))
  output$p_fgf <- renderPlot(
    lineplot(sim(), "FGF23RU", "FGF23 (상대단위)", "FGF23"))

  ## ---- Tab 4 : renal -----------------------------------------------------
  output$p_uca <- renderPlot(
    lineplot(sim(), "UCA24", "24h 요중 칼슘 (mg/일)",
             "요중 칼슘 — 이 질환의 진짜 대가", bands = c(0, 300)))
  output$p_feca <- renderPlot(
    lineplot(sim(), "FEXCA", "FE_Ca (%)", "칼슘 분획 배설률"))
  output$p_tmp <- renderPlot(
    lineplot(sim(), "TMPPERG", "TmP/GFR (mg/dL)", "신장 인 역치",
             bands = c(2.5, 4.2)))
  output$p_caxp <- renderPlot(
    lineplot(sim(), "CAPROD", "Ca × P (mg²/dL²)", "칼슘-인 곱",
             bands = c(0, 45)))
  output$p_nc <- renderPlot(
    lineplot(sim(), "NCBURD", "신석회화증 부담 (단위)", "신석회화증 축적"))
  output$p_egfr <- renderPlot(
    lineplot(sim(), "EGFR", "eGFR (mL/min/1.73m²)", "신기능 궤적"))

  ## ---- Tab 5 : bone ------------------------------------------------------
  output$p_turn <- renderPlot(
    lineplot(sim(), "TURNOV", "골회전 지수 (1 = 정상)", "골재형성 활성"))
  output$p_cells <- renderPlot({
    d <- sim() %>% select(time, OB, OC) %>% pivot_longer(-time)
    ggplot(d, aes(time/24, value, colour = name)) +
      geom_line(linewidth = 0.8) +
      scale_colour_manual(values = c(OB = "#7D3C98", OC = "#C0392B"),
                          labels = c(OB = "조골세포", OC = "파골세포")) +
      labs(x = "시간 (일)", y = "활성 (정상 = 1)",
           title = "조골세포 · 파골세포", colour = NULL) + theme_hp()
  })
  output$p_markers <- renderPlot({
    d <- sim() %>% select(time, P1NP, CTX, BSAP) %>% pivot_longer(-time)
    ggplot(d, aes(time/24, value, colour = name)) +
      geom_line(linewidth = 0.8) +
      facet_wrap(~name, scales = "free_y", ncol = 1) +
      scale_colour_manual(values = PAL) +
      labs(x = "시간 (일)", y = NULL, title = "골표지자", colour = NULL) +
      theme_hp() + theme(legend.position = "none")
  })
  output$p_bmd <- renderPlot(
    lineplot(sim(), "BMDZS", "요추 BMD Z-score", "골밀도 궤적"))

  ## ---- Tab 6 : endpoints -------------------------------------------------
  output$p_sym <- renderPlot(
    lineplot(sim(), "SYMSCOR", "증상 점수 (0-100)", "신경근 증상"))
  output$p_qtc <- renderPlot(
    lineplot(sim(), "QTCMS", "QTc (ms)", "QTc 간격", bands = c(350, 450)))
  output$p_qol <- renderPlot(
    lineplot(sim(), "QOLSC", "QoL 점수 (0-100)", "삶의 질"))
  output$p_bg <- renderPlot(
    lineplot(sim(), "BGCALCI", "기저핵 석회화 (단위)", "기저핵 석회화 축적"))

  output$txt_resp3 <- renderPrint({
    d  <- sim()
    ca <- ss(d, "CACORRD")
    tp <- treat_pars()
    cat(sprintf("알부민 보정 칼슘   : %.2f mg/dL  (%s)\n", ca,
                if (ca >= 8.0 && ca <= 10.6) "충족" else "미충족"))
    cat(sprintf("경구 칼슘          : %.0f mg/일   (%s)\n", tp$DCASUP,
                if (tp$DCASUP <= 600) "충족" else "미충족"))
    cat(sprintf("활성 비타민 D      : %.3f µg/일  (%s)\n", tp$DCTR,
                if (tp$DCTR <= 0) "충족" else "미충족"))
    cat(sprintf("\n삼중 복합 반응     : %s\n",
                if (mean(d$RESP3) > 0.5) "달성" else "미달성"))
    cat(sprintf("관찰 기간 중 달성 비율: %.0f%%\n", 100*mean(d$RESP3)))
  })

  ## ---- Tab 7 : scenarios -------------------------------------------------
  scn <- eventReactive(input$goscn, {
    withProgress(message = "시나리오 실행 중 (수십 초 소요)...", value = 0.5, {
      HYPOPT_simulate_scenarios(MOD, which = as.integer(input$scn))
    })
  })

  output$s_ca   <- renderPlot(lineplot(scn(), "CACORRD", "보정 칼슘 (mg/dL)",
                                       "혈청 칼슘", bands = c(8.0, 9.0)))
  output$s_uca  <- renderPlot(lineplot(scn(), "UCA24", "24h 요중 칼슘 (mg/일)",
                                       "요중 칼슘", bands = c(0, 300)))
  output$s_pi   <- renderPlot(lineplot(scn(), "PISER", "인 (mg/dL)", "혈청 인",
                                       bands = c(2.5, 4.5)))
  output$s_egfr <- renderPlot(lineplot(scn(), "EGFR", "eGFR", "신기능"))

  output$s_tbl   <- DT::renderDataTable(
    DT::datatable(HYPOPT_summary(scn()), options = list(pageLength = 15,
                                                        scrollX = TRUE)))
  output$s_swing <- DT::renderDataTable({
    d <- scn()
    from <- max(0, floor(max(d$time)/24) - 30)
    DT::datatable(HYPOPT_swing(d, from_day = from),
                  options = list(pageLength = 15, scrollX = TRUE))
  })

  ## ---- Tab 8 : biomarkers ------------------------------------------------
  output$b_tbl <- DT::renderDataTable({
    d <- sim()
    cols <- c("CACORRD","CAIONMM","PISER","MGSER","PTHPG","PTHENDO","PTHEXOG",
              "D125PG","D25NG","FGF23RU","TMPPERG","FRCAPCT","FEXCA","UCA24",
              "UPI24","UCAMGKG","CAABS","FABSPCT","CAPROD","EGFR","NCBURD",
              "BMDZS","P1NP","CTX","BSAP","TURNOV","SYMSCOR","QTCMS","QOLSC")
    lab <- c("보정 칼슘 (mg/dL)","이온화 칼슘 (mmol/L)","인 (mg/dL)",
             "마그네슘 (mg/dL)","PTH 신호 (pg/mL eq)","내인성 PTH (pg/mL)",
             "외인성 PTH (pg/mL)","1,25D (pg/mL)","25D (ng/mL)","FGF23 (RU)",
             "TmP/GFR (mg/dL)","세뇨관 Ca 재흡수 (%)","FE_Ca (%)",
             "24h 요중 Ca (mg)","24h 요중 P (mg)","요중 Ca (mg/kg/일)",
             "흡수 칼슘 (mg/일)","칼슘 흡수분율 (%)","Ca×P (mg²/dL²)",
             "eGFR","신석회화증 (단위)","BMD Z-score","P1NP (µg/L)",
             "CTX (ng/mL)","BSAP (µg/L)","골회전 지수","증상 점수","QTc (ms)",
             "QoL 점수")
    DT::datatable(
      data.frame(지표 = lab,
                 정상상태 = round(vapply(cols, function(c) ss(d, c), 0), 3)),
      options = list(pageLength = 30), rownames = FALSE)
  })

  grid <- eventReactive(input$gogrid, {
    withProgress(message = "용량 격자 계산 중...", value = 0.5, {
      HYPOPT_titrate_conventional(MOD)
    })
  })
  output$b_grid <- DT::renderDataTable(
    DT::datatable(grid(), options = list(pageLength = 16, scrollX = TRUE)))

  ## ---- Tab 9 : references ------------------------------------------------
  output$refs <- renderUI({
    f <- "hypopt_references.md"
    if (file.exists(f)) includeMarkdown(f) else p("참고문헌 파일을 찾을 수 없습니다.")
  })
}

shinyApp(ui, server)
