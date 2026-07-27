## ============================================================================
## Urea Cycle Disorders (UCD) QSP — Shiny Dashboard
## ----------------------------------------------------------------------------
## 8 tabs:
##   1 Patient profile   — build the virtual UCD patient (locus, residual
##                         activity, diet) and read the analytically derived
##                         untreated baseline BEFORE any drug is given
##   2 Drug PK           — phenylbutyrate -> phenylacetate -> PAGN and the
##                         benzoate -> hippurate axis; NaPBA vs GPB profiles
##   3 Ammonia & N flux  — plasma ammonia, glutamine and the complete
##                         nitrogen-disposal ledger (urea / PAGN / hippurate /
##                         renal), in umol N/h and as % of the daily load
##   4 Brain             — brain ammonia, the astrocytic glutamine trap,
##                         myo-inositol depletion, brain water, ICP, HE stage
##   5 Clinical endpoints— coma hours, cumulative injury, modelled IQ and the
##                         natural-protein tolerance that patients actually feel
##   6 Scenario compare  — all 14 prebuilt arms side by side + summary table
##   7 Biomarkers        — U-PAGN, PAA:PAGN ratio, orotate, citrulline,
##                         arginine, glycine, BCAA: the monitoring panel
##   8 Safety            — phenylacetate neurotoxicity, sodium load, BCAA and
##                         glycine depletion, dialysis rebound
##
## Dependencies: shiny, shinydashboard, mrgsolve, dplyr, tidyr, ggplot2, DT
## Run with:  shiny::runApp("ucd_shiny_app.R")
##            (keep ucd_mrgsolve_model.R in the same directory)
##
## EDUCATIONAL / RESEARCH USE ONLY — not validated for clinical decisions.
## ----------------------------------------------------------------------------
suppressPackageStartupMessages({
  library(shiny)
  library(shinydashboard)
  library(mrgsolve)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(DT)
})

## ---------------------------------------------------------------- model load
UCD_SPEC <- "ucd_mrgsolve_model.R"
if (!file.exists(UCD_SPEC)) {
  stop("ucd_mrgsolve_model.R not found — keep it in the same directory as the app.")
}
source(UCD_SPEC, local = FALSE)
MOD <- UCD_build()

THEME <- theme_bw(base_size = 12) +
  theme(panel.grid.minor = element_blank(),
        legend.position  = "bottom",
        legend.title     = element_blank(),
        plot.title       = element_text(face = "bold", size = 12))

PAL <- c("#1b6ca8", "#c0392b", "#27865c", "#b8860b", "#7d3c98",
         "#16a085", "#d35400", "#2c3e50")

## band helper — draws a reference band behind a line plot
band <- function(lo, hi, fill = "#2ecc71", alpha = 0.10) {
  annotate("rect", xmin = -Inf, xmax = Inf, ymin = lo, ymax = hi,
           fill = fill, alpha = alpha)
}

## ============================================================================
## UI
## ============================================================================
ui <- dashboardPage(
  skin = "blue",
  dashboardHeader(title = "UCD QSP 대시보드", titleWidth = 300),

  dashboardSidebar(
    width = 300,
    sidebarMenu(
      id = "tabs",
      menuItem("1. 환자 프로파일",      tabName = "patient",  icon = icon("user")),
      menuItem("2. 약물 PK",            tabName = "pk",       icon = icon("pills")),
      menuItem("3. 암모니아·질소 균형", tabName = "nitrogen", icon = icon("scale-balanced")),
      menuItem("4. 뇌 (별아교세포)",    tabName = "brain",    icon = icon("brain")),
      menuItem("5. 임상 엔드포인트",    tabName = "endpoint", icon = icon("stethoscope")),
      menuItem("6. 시나리오 비교",      tabName = "compare",  icon = icon("layer-group")),
      menuItem("7. 바이오마커",         tabName = "biomark",  icon = icon("vial")),
      menuItem("8. 안전성",             tabName = "safety",   icon = icon("triangle-exclamation"))
    ),
    hr(),
    h4("환자 (Patient)", style = "padding-left:15px;"),
    selectInput("locus", "결손 효소 (Enzyme block)",
                choices = c("OTC (X-linked)"          = "OTC",
                            "CPS1"                    = "CPS1",
                            "NAGS"                    = "NAGS",
                            "ASS1 / ASL (distal arm)" = "DISTAL"),
                selected = "OTC"),
    sliderInput("act",  "잔존 효소 활성 Residual activity (%)",
                min = 0, max = 100, value = 32, step = 1),
    sliderInput("wt",   "체중 Weight (kg)", min = 3, max = 100, value = 70, step = 0.5),
    sliderInput("vbr",  "뇌 용적 Brain volume (L)", min = 0.35, max = 1.5,
                value = 1.35, step = 0.05),
    sliderInput("prot", "천연 단백 섭취 Natural protein (g/kg/day)",
                min = 0, max = 2.5, value = 0.60, step = 0.05),

    hr(),
    h4("치료 (Therapy)", style = "padding-left:15px;"),
    sliderInput("gpb",  "글리세롤 페닐부티레이트 GPB (mL/day)",
                min = 0, max = 25, value = 0, step = 0.5),
    sliderInput("napba","나트륨 페닐부티레이트 NaPBA (g/day)",
                min = 0, max = 30, value = 0, step = 0.5),
    sliderInput("bz",   "벤조산나트륨 Na-benzoate (mg/kg/day)",
                min = 0, max = 500, value = 0, step = 25),
    sliderInput("cit",  "L-시트룰린 L-citrulline (mg/kg/day)",
                min = 0, max = 250, value = 0, step = 10),
    sliderInput("ncg",  "카글룸산 Carglumic acid (mg/kg/day)",
                min = 0, max = 250, value = 0, step = 10),

    hr(),
    h4("스트레스·구조 (Stress / rescue)", style = "padding-left:15px;"),
    sliderInput("cat",  "이화 자극 강도 Catabolic amplitude",
                min = 0, max = 4, value = 0, step = 0.1),
    sliderInput("catt", "이화 시작 시각 (h)", min = 0, max = 168, value = 24, step = 2),
    sliderInput("catd", "이화 지속 (h)", min = 6, max = 168, value = 48, step = 6),
    sliderInput("anab", "동화 구조 Anabolic rescue (0-1)",
                min = 0, max = 1, value = 0, step = 0.1),
    sliderInput("hd",   "투석 청소율 Dialysis CL (L/h)",
                min = 0, max = 25, value = 0, step = 0.5),
    sliderInput("hdt",  "투석 시작 (h)", min = 0, max = 168, value = 36, step = 2),
    sliderInput("hdd",  "투석 지속 (h)", min = 1, max = 48, value = 6, step = 1),

    hr(),
    sliderInput("days", "시뮬레이션 기간 Horizon (days)",
                min = 2, max = 180, value = 14, step = 1),
    actionButton("go", "시뮬레이션 실행 (Run)", icon = icon("play"),
                 class = "btn-primary", style = "margin:10px 15px;width:80%;")
  ),

  dashboardBody(
    tags$head(tags$style(HTML(".content-wrapper{background-color:#f7f9fb;}"))),
    tabItems(

      ## ---------------------------------------------------------- 1 patient
      tabItem(
        "patient",
        fluidRow(
          box(width = 12, title = "무엇을 보는 화면인가 (What this shows)",
              status = "primary", solidHeader = TRUE,
              HTML("모델은 환자 기술자(결손 효소·잔존 활성·체중·단백 섭취)로부터
                    <b>치료 전 정상상태를 해석적으로 계산</b>하여 모든 구획의
                    초기값으로 사용합니다. 따라서 약을 주지 않으면 곡선은
                    완전히 평평하며, 시뮬레이션에서 보이는 모든 변화는
                    <b>기전에서 나온 것이지 초기값 오차가 아닙니다.</b><br>
                    The untreated steady state is solved from the patient
                    descriptors and used as the initial condition, so an
                    untreated arm is exactly stationary."))),
        fluidRow(
          valueBoxOutput("vb_nh3",  width = 3),
          valueBoxOutput("vb_gln",  width = 3),
          valueBoxOutput("vb_tol",  width = 3),
          valueBoxOutput("vb_cap",  width = 3)),
        fluidRow(
          box(width = 6, title = "잔존 활성 vs 기저 암모니아 (Activity–ammonia cliff)",
              status = "info", solidHeader = TRUE, plotOutput("p_cliff", height = 330)),
          box(width = 6, title = "단백 섭취 vs 기저 암모니아 (Protein–ammonia cliff)",
              status = "info", solidHeader = TRUE, plotOutput("p_protcliff", height = 330))),
        fluidRow(
          box(width = 12, title = "유도된 기저 상태 (Derived baseline)",
              status = "info", solidHeader = TRUE, DTOutput("t_base")))),

      ## --------------------------------------------------------------- 2 PK
      tabItem(
        "pk",
        fluidRow(
          box(width = 6, title = "페닐부티레이트 / 페닐아세테이트 (PBA · PAA)",
              status = "primary", solidHeader = TRUE, plotOutput("p_pba", height = 320)),
          box(width = 6, title = "페닐아세틸글루타민 (PAGN)",
              status = "primary", solidHeader = TRUE, plotOutput("p_pagn", height = 320))),
        fluidRow(
          box(width = 6, title = "벤조산 · 히푸르산 (Benzoate · hippurate)",
              status = "info", solidHeader = TRUE, plotOutput("p_bz", height = 300)),
          box(width = 6, title = "요중 배설 (Urinary PAGN & hippurate)",
              status = "info", solidHeader = TRUE, plotOutput("p_uexc", height = 300))),
        fluidRow(
          box(width = 12, title = "NaPBA vs GPB — 같은 PBA 몰수, 다른 프로파일",
              status = "warning", solidHeader = TRUE,
              plotOutput("p_pbavsgpb", height = 320),
              HTML("<small>동일한 PBA 몰수(NaPBA 20 g/day ≈ GPB 17.5 mL/day ≈ 107-109 mmol/day)를
                    투여해도 GPB는 췌장 리파아제 의존적 서방출 때문에 <b>Cmax가 낮고 24시간
                    프로파일이 평탄</b>합니다 — Diaz 2013 Hepatology 57:2171의 교차설계 결과.</small>")))),

      ## --------------------------------------------------------- 3 nitrogen
      tabItem(
        "nitrogen",
        fluidRow(
          box(width = 8, title = "혈장 암모니아 (Plasma ammonia)",
              status = "primary", solidHeader = TRUE, plotOutput("p_nh3", height = 340)),
          box(width = 4, title = "질소 처리 경로 비중 (Nitrogen disposal share)",
              status = "primary", solidHeader = TRUE, plotOutput("p_nshare", height = 340))),
        fluidRow(
          box(width = 6, title = "혈장 글루타민 — 조기 경보 지표",
              status = "info", solidHeader = TRUE, plotOutput("p_gln", height = 300),
              HTML("<small>글루타민은 암모니아보다 <b>먼저</b> 오릅니다. 1000 µmol/L 초과는
                    임상적으로 위기 경보로 취급됩니다.</small>")),
          box(width = 6, title = "심부 조직 글루타민 저장고 (Deep reservoir)",
              status = "info", solidHeader = TRUE, plotOutput("p_glndeep", height = 300),
              HTML("<small>투석 후 반동(rebound)의 원천입니다.</small>"))),
        fluidRow(
          box(width = 12, title = "질소 수지 (Nitrogen ledger, µmol N/h)",
              status = "info", solidHeader = TRUE, plotOutput("p_nledger", height = 320)))),

      ## ------------------------------------------------------------ 4 brain
      tabItem(
        "brain",
        fluidRow(
          box(width = 12, title = "삼투압 가설 (The osmotic cascade)",
              status = "primary", solidHeader = TRUE,
              HTML("혈중 NH<sub>3</sub> → BBB 확산 → <b>별아교세포 glutamine synthetase</b> →
                    글루타민 축적 → 삼투질로 작용 → 세포 팽창. 방어기전은 <b>myo-inositol 배출</b>인데
                    <b>만성 고암모니아혈증에서는 이미 고갈</b>되어 있어, 급성 상승 시 완충이 없습니다.
                    이 상호작용이 UCD에서 가장 중요한 비직관적 거동입니다."))),
        fluidRow(
          box(width = 6, title = "뇌 암모니아 (Brain ammonia)",
              status = "info", solidHeader = TRUE, plotOutput("p_nh3b", height = 300)),
          box(width = 6, title = "뇌 글루타민 (Brain glutamine, MRS 2.05-2.45 ppm)",
              status = "info", solidHeader = TRUE, plotOutput("p_glnb", height = 300))),
        fluidRow(
          box(width = 6, title = "Myo-inositol (삼투질 예비능)",
              status = "warning", solidHeader = TRUE, plotOutput("p_mins", height = 300)),
          box(width = 6, title = "뇌 부종 & 두개내압 (Brain water · ICP)",
              status = "danger", solidHeader = TRUE, plotOutput("p_bw", height = 300))),
        fluidRow(
          box(width = 12, title = "뇌병증 단계 (Encephalopathy stage 0-4)",
              status = "danger", solidHeader = TRUE, plotOutput("p_he", height = 300)))),

      ## -------------------------------------------------------- 5 endpoints
      tabItem(
        "endpoint",
        fluidRow(
          valueBoxOutput("vb_coma", width = 3),
          valueBoxOutput("vb_iq",   width = 3),
          valueBoxOutput("vb_peak", width = 3),
          valueBoxOutput("vb_auc",  width = 3)),
        fluidRow(
          box(width = 6, title = "혼수 범위 누적 시간 (Cumulative coma-range hours)",
              status = "danger", solidHeader = TRUE, plotOutput("p_coma", height = 300)),
          box(width = 6, title = "누적 신경 손상 지수 & 예측 IQ",
              status = "danger", solidHeader = TRUE, plotOutput("p_iq", height = 300))),
        fluidRow(
          box(width = 6, title = "암모니아 AUC>100 (µmol/L·h)",
              status = "warning", solidHeader = TRUE, plotOutput("p_auc", height = 300)),
          box(width = 6, title = "천연 단백 내성 (Natural-protein tolerance)",
              status = "success", solidHeader = TRUE, plotOutput("p_tol", height = 300),
              HTML("<small>환자가 실제로 체감하는 치료 효과 —
                    ‘무엇을 얼마나 먹을 수 있는가’.</small>")))),

      ## --------------------------------------------------------- 6 compare
      tabItem(
        "compare",
        fluidRow(
          box(width = 12, title = "사전 정의 시나리오 (14 prebuilt arms)",
              status = "primary", solidHeader = TRUE,
              checkboxGroupInput("scen", NULL, inline = TRUE,
                                 choices  = setNames(1:14, paste0("S", 1:14)),
                                 selected = c(1, 3, 4, 5, 6, 12)),
              actionButton("goscen", "시나리오 실행 (Run scenarios)",
                           icon = icon("play"), class = "btn-primary"))),
        fluidRow(
          box(width = 6, title = "암모니아", status = "info", solidHeader = TRUE,
              plotOutput("c_nh3", height = 330)),
          box(width = 6, title = "글루타민", status = "info", solidHeader = TRUE,
              plotOutput("c_gln", height = 330))),
        fluidRow(
          box(width = 6, title = "단백 내성", status = "info", solidHeader = TRUE,
              plotOutput("c_tol", height = 330)),
          box(width = 6, title = "뇌 글루타민", status = "info", solidHeader = TRUE,
              plotOutput("c_glnb", height = 330))),
        fluidRow(
          box(width = 12, title = "요약표 (Scenario summary)",
              status = "primary", solidHeader = TRUE, DTOutput("c_tab")))),

      ## -------------------------------------------------------- 7 biomarker
      tabItem(
        "biomark",
        fluidRow(
          box(width = 6, title = "요중 PAGN — 용량 적정성·복약순응도",
              status = "primary", solidHeader = TRUE, plotOutput("b_upagn", height = 300)),
          box(width = 6, title = "PAA : PAGN 비 — 포합 포화의 지표",
              status = "warning", solidHeader = TRUE, plotOutput("b_ratio", height = 300),
              HTML("<small>µg/mL 기준 2.5 초과 = 포합 포화, PAA가 오를 것이라는 경고.</small>"))),
        fluidRow(
          box(width = 4, title = "오로트산 (Orotic acid)", status = "info",
              solidHeader = TRUE, plotOutput("b_orot", height = 280),
              HTML("<small>낮은 시트룰린 + <b>높은 오로트산</b> = OTC 결핍;
                    낮은 시트룰린 + 정상 오로트산 = CPS1/NAGS 결핍.</small>")),
          box(width = 4, title = "시트룰린 (Citrulline)", status = "info",
              solidHeader = TRUE, plotOutput("b_cit", height = 280)),
          box(width = 4, title = "아르기닌 (Arginine)", status = "info",
              solidHeader = TRUE, plotOutput("b_arg", height = 280))),
        fluidRow(
          box(width = 6, title = "글리신 (Glycine) — 벤조산이 소모",
              status = "info", solidHeader = TRUE, plotOutput("b_gly", height = 280)),
          box(width = 6, title = "분지사슬아미노산 (BCAA) — 페닐부티레이트가 소모",
              status = "info", solidHeader = TRUE, plotOutput("b_bcaa", height = 280)))),

      ## ----------------------------------------------------------- 8 safety
      tabItem(
        "safety",
        fluidRow(
          box(width = 12, title = "안전성 신호 (Safety signals)",
              status = "danger", solidHeader = TRUE,
              HTML("① <b>페닐아세테이트 신경독성</b> — 500 µg/mL 초과에서 졸림·오심·두통.
                    신부전·간부전·고용량에서 위험이 커집니다.<br>
                    ② <b>나트륨 부하</b> — NaPBA와 벤조산나트륨은 상당한 Na를 함께 투여합니다.<br>
                    ③ <b>BCAA·글리신 고갈</b> — 각각 페닐부티레이트와 벤조산의 특징적 부작용.<br>
                    ④ <b>투석 후 반동</b> — 심부 글루타민 저장고에서 질소가 재분포합니다.<br>
                    ⑤ <b>발프로산</b>은 NAGS/CPS1을 억제하므로 UCD에서 금기에 가깝습니다."))),
        fluidRow(
          box(width = 6, title = "혈장 PAA vs 500 µg/mL 독성 문턱",
              status = "danger", solidHeader = TRUE, plotOutput("s_paa", height = 320)),
          box(width = 6, title = "누적 PAA 초과 노출 (µg/mL·h)",
              status = "danger", solidHeader = TRUE, plotOutput("s_paatox", height = 320))),
        fluidRow(
          box(width = 6, title = "나트륨 부하 추정 (mmol/day)",
              status = "warning", solidHeader = TRUE, plotOutput("s_na", height = 300)),
          box(width = 6, title = "투석 후 반동 (Post-dialysis rebound)",
              status = "warning", solidHeader = TRUE, plotOutput("s_rebound", height = 300))),
        fluidRow(
          box(width = 12, title = "면책 (Disclaimer)", status = "danger",
              solidHeader = TRUE,
              HTML("<b>본 대시보드는 교육·연구 목적의 반정량적 QSP 모델입니다.</b>
                    임상 의사결정, 처방, 규제 제출에 사용해서는 안 됩니다."))))
    )
  )
)

## ============================================================================
## Server
## ============================================================================
server <- function(input, output, session) {

  ## ---------------------------------------------------- parameter assembly
  pars <- reactive({
    a <- input$act/100
    p <- list(WT = input$wt, VBR = input$vbr, PROT = input$prot,
              OTCACT = 1, CPS1A = 1, NAGSA = 1, DISTALA = 1,
              CATAMP = input$cat, CATT0 = input$catt, CATDUR = input$catd,
              ANAB = input$anab,
              CLHD = input$hd, HDT0 = ifelse(input$hd > 0, input$hdt, 1e9),
              HDDUR = input$hdd,
              BSA = round(0.024265 * (input$wt^0.5378) * (170^0.3964), 3))
    switch(input$locus,
           OTC    = { p$OTCACT  <- a },
           CPS1   = { p$CPS1A   <- a },
           NAGS   = { p$NAGSA   <- a },
           DISTAL = { p$DISTALA <- a })
    ## citrulline is diagnostic: low in proximal blocks, high in distal ones
    p$CITBASE <- if (input$locus %in% c("OTC", "CPS1", "NAGS")) 12 else 180
    p
  })

  evts <- reactive({
    d  <- input$days
    e  <- NULL
    add <- function(x) if (is.null(e)) e <<- x else e <<- c(e, x)
    if (input$gpb   > 0) add(ucd_gpb(input$gpb,  n = 3, ii = 8, addl = ucd_addl(d, 8)))
    if (input$napba > 0) add(ucd_napba(input$napba, n = 3, ii = 8, addl = ucd_addl(d, 8)))
    if (input$bz    > 0) add(ucd_benzoate(input$bz, input$wt, n = 3, ii = 8,
                                          addl = ucd_addl(d, 8)))
    if (input$cit   > 0) add(ucd_citrulline(input$cit, input$wt, n = 3, ii = 8,
                                            addl = ucd_addl(d, 8)))
    if (input$ncg   > 0) add(ucd_carglumic(input$ncg, input$wt, n = 2, ii = 12,
                                           addl = ucd_addl(d, 12)))
    e
  })

  sim <- eventReactive(input$go, {
    m <- do.call(mrgsolve::param, c(list(MOD), pars()))
    dl <- max(0.1, input$days*24/1500)
    as.data.frame(mrgsolve::mrgsim(m, events = evts(), end = input$days*24,
                                   delta = dl, hmax = 0.25))
  }, ignoreNULL = FALSE)

  base <- reactive({
    m <- do.call(mrgsolve::param, c(list(MOD), pars()))
    as.data.frame(mrgsolve::mrgsim(m, end = 1, delta = 1))[1, ]
  })

  tt <- function(d) d$time/24

  ## -------------------------------------------------------------- 1 patient
  output$vb_nh3 <- renderValueBox({
    v <- round(base()$NH3, 1)
    valueBox(paste0(v, " µmol/L"), "기저 혈장 암모니아 (baseline NH3)",
             icon = icon("droplet"),
             color = if (v < 50) "green" else if (v < 100) "yellow" else "red")
  })
  output$vb_gln <- renderValueBox({
    v <- round(base()$GLNP)
    valueBox(paste0(v, " µmol/L"), "기저 혈장 글루타민 (baseline Gln)",
             icon = icon("flask"),
             color = if (v < 800) "green" else if (v < 1000) "yellow" else "red")
  })
  output$vb_tol <- renderValueBox({
    v <- round(base()$PROTTOL, 2)
    valueBox(paste0(v, " g/kg/d"), "천연 단백 내성 (protein tolerance)",
             icon = icon("drumstick-bite"),
             color = if (v > 1.0) "green" else if (v > 0.6) "yellow" else "red")
  })
  output$vb_cap <- renderValueBox({
    v <- round(100*base()$ACTEFF, 1)
    valueBox(paste0(v, " %"), "유효 요소회로 활성 (effective activity)",
             icon = icon("gauge-high"),
             color = if (v > 50) "green" else if (v > 20) "yellow" else "red")
  })

  output$p_cliff <- renderPlot({
    p <- pars()
    grid <- seq(0.01, 1, by = 0.01)
    y <- vapply(grid, function(a) {
      q <- p
      switch(input$locus, OTC = {q$OTCACT <- a}, CPS1 = {q$CPS1A <- a},
             NAGS = {q$NAGSA <- a}, DISTAL = {q$DISTALA <- a})
      m <- do.call(mrgsolve::param, c(list(MOD), q))
      as.data.frame(mrgsolve::mrgsim(m, end = 1, delta = 1))$NH3[1]
    }, numeric(1))
    ggplot(data.frame(a = 100*grid, y = y), aes(a, y)) +
      band(0, 35) + geom_line(colour = PAL[1], linewidth = 1.1) +
      geom_vline(xintercept = input$act, linetype = 2, colour = PAL[2]) +
      geom_hline(yintercept = 100, linetype = 3) +
      scale_y_log10() +
      labs(x = "잔존 효소 활성 (%)", y = "기저 암모니아 (µmol/L, log)",
           title = "질환이 경사가 아니라 절벽인 이유") + THEME
  })

  output$p_protcliff <- renderPlot({
    p <- pars()
    grid <- seq(0.1, 2.5, by = 0.05)
    y <- vapply(grid, function(g) {
      q <- p; q$PROT <- g
      m <- do.call(mrgsolve::param, c(list(MOD), q))
      as.data.frame(mrgsolve::mrgsim(m, end = 1, delta = 1))$NH3[1]
    }, numeric(1))
    ggplot(data.frame(g = grid, y = y), aes(g, y)) +
      band(0, 35) + geom_line(colour = PAL[3], linewidth = 1.1) +
      geom_vline(xintercept = input$prot, linetype = 2, colour = PAL[2]) +
      geom_hline(yintercept = 100, linetype = 3) +
      scale_y_log10() +
      labs(x = "천연 단백 섭취 (g/kg/day)", y = "기저 암모니아 (µmol/L, log)",
           title = "단백 처방이 왜 그렇게 까다로운가") + THEME
  })

  output$t_base <- renderDT({
    b <- base()
    out <- data.frame(
      item = c("혈장 암모니아 (µmol/L)", "혈장 글루타민 (µmol/L)",
               "뇌 암모니아 (µmol/L)", "뇌 글루타민 (mmol/L)",
               "뇌 myo-inositol (mmol/L)", "뇌 수분 초과 (%)",
               "두개내압 (mmHg)", "혈장 시트룰린 (µmol/L)",
               "혈장 아르기닌 (µmol/L)", "요소 질소 처리 (µmol N/h)",
               "총 질소 처리 (µmol N/h)", "천연 단백 내성 (g/kg/d)"),
      value = round(c(b$NH3, b$GLNP, b$NH3BR, b$GLNBRAIN, b$MINSMM, b$BW, b$ICP,
                      b$CITP, b$ARGP, b$NUREA, b$NTOTAL, b$PROTTOL), 2),
      stringsAsFactors = FALSE)
    colnames(out) <- c("\uD56D\uBAA9 (Item)", "\uAC12 (Value)")
    out
  }, options = list(dom = "t", pageLength = 12), rownames = FALSE)

  ## ------------------------------------------------------------------ 2 PK
  output$p_pba <- renderPlot({
    d <- sim()
    dd <- d |> select(time, PBAUG, PAAUG) |>
      pivot_longer(-time, names_to = "k", values_to = "v")
    ggplot(dd, aes(time/24, v, colour = k)) + geom_line(linewidth = 0.9) +
      scale_colour_manual(values = PAL[c(1, 2)],
                          labels = c("PAA (활성체)", "PBA (전구체)")) +
      labs(x = "Days", y = "µg/mL", title = "PBA → PAA") + THEME
  })
  output$p_pagn <- renderPlot({
    ggplot(sim(), aes(time/24, PGNUG)) +
      geom_line(colour = PAL[3], linewidth = 0.9) +
      labs(x = "Days", y = "PAGN (µg/mL)",
           title = "페닐아세틸글루타민 — 질소 2개를 싣고 나간다") + THEME
  })
  output$p_bz <- renderPlot({
    dd <- sim() |> select(time, BZUG, HIPUG) |>
      pivot_longer(-time, names_to = "k", values_to = "v")
    ggplot(dd, aes(time/24, v, colour = k)) + geom_line(linewidth = 0.9) +
      scale_colour_manual(values = PAL[c(4, 5)],
                          labels = c("벤조산", "히푸르산")) +
      labs(x = "Days", y = "µg/mL", title = "벤조산 → 히푸르산 (질소 1개)") + THEME
  })
  output$p_uexc <- renderPlot({
    dd <- sim() |> select(time, UPAGN, UHIPP) |>
      pivot_longer(-time, names_to = "k", values_to = "v")
    ggplot(dd, aes(time/24, v, colour = k)) + geom_line(linewidth = 0.9) +
      scale_colour_manual(values = PAL[c(3, 5)],
                          labels = c("요중 PAGN", "요중 히푸르산")) +
      labs(x = "Days", y = "mmol/day", title = "요중 배설 속도") + THEME
  })
  output$p_pbavsgpb <- renderPlot({
    m <- do.call(mrgsolve::param, c(list(MOD), pars()))
    a <- as.data.frame(mrgsolve::mrgsim(
      m, events = ucd_napba(20, n = 3, ii = 8, addl = 8), end = 72,
      delta = 0.1, hmax = 0.25)); a$arm <- "NaPBA 20 g/day TID"
    b <- as.data.frame(mrgsolve::mrgsim(
      m, events = ucd_gpb(17.5, n = 3, ii = 8, addl = 8), end = 72,
      delta = 0.1, hmax = 0.25)); b$arm <- "GPB 17.5 mL/day TID"
    ggplot(rbind(a, b), aes(time, PAAUG, colour = arm)) +
      geom_line(linewidth = 0.9) +
      scale_colour_manual(values = PAL[c(2, 1)]) +
      labs(x = "Hours", y = "혈장 PAA (µg/mL)",
           title = "동일 PBA 몰수, 다른 Cmax") + THEME
  })

  ## ------------------------------------------------------------ 3 nitrogen
  output$p_nh3 <- renderPlot({
    ggplot(sim(), aes(time/24, NH3)) +
      band(0, 35) +
      geom_hline(yintercept = c(100, 200, 360), linetype = c(3, 2, 1),
                 colour = c("#b8860b", "#d35400", "#c0392b")) +
      geom_line(colour = PAL[1], linewidth = 1.0) +
      labs(x = "Days", y = "혈장 암모니아 (µmol/L)",
           title = "100 = 치료, 200 = 혼수 위험, 360 = 투석 고려") + THEME
  })
  output$p_nshare <- renderPlot({
    d <- sim(); k <- d[nrow(d), ]
    df <- data.frame(
      route = c("요소 (urea)", "PAGN", "히푸르산", "신장", "시트룰린 경로"),
      v = c(k$NUREA, k$NPAGN, k$NHIPP, k$NRENL, k$NCITR))
    df$route <- factor(df$route, levels = df$route)
    ggplot(df, aes(route, v, fill = route)) +
      geom_col(width = 0.7) + scale_fill_manual(values = PAL[1:5]) +
      labs(x = NULL, y = "µmol N/h", title = "최종 시점 질소 처리 경로") +
      THEME + theme(legend.position = "none",
                    axis.text.x = element_text(angle = 20, hjust = 1))
  })
  output$p_gln <- renderPlot({
    ggplot(sim(), aes(time/24, GLNP)) +
      band(450, 750) + geom_hline(yintercept = 1000, linetype = 2,
                                  colour = PAL[2]) +
      geom_line(colour = PAL[3], linewidth = 1.0) +
      labs(x = "Days", y = "혈장 글루타민 (µmol/L)", title = NULL) + THEME
  })
  output$p_glndeep <- renderPlot({
    dd <- sim() |> select(time, GLNP, GLNDEEP) |>
      pivot_longer(-time, names_to = "k", values_to = "v")
    ggplot(dd, aes(time/24, v, colour = k)) + geom_line(linewidth = 0.9) +
      scale_colour_manual(values = PAL[c(6, 3)],
                          labels = c("심부 저장고", "혈장")) +
      labs(x = "Days", y = "µmol/L", title = NULL) + THEME
  })
  output$p_nledger <- renderPlot({
    dd <- sim() |> select(time, NUREA, NPAGN, NHIPP, NRENL, NCITR) |>
      pivot_longer(-time, names_to = "k", values_to = "v")
    dd$k <- factor(dd$k, levels = c("NUREA", "NPAGN", "NHIPP", "NCITR", "NRENL"))
    ggplot(dd, aes(time/24, v, fill = k)) +
      geom_area(alpha = 0.85) +
      scale_fill_manual(values = PAL[1:5],
                        labels = c("요소", "PAGN", "히푸르산", "시트룰린", "신장")) +
      labs(x = "Days", y = "µmol N/h", title = "질소 처리 경로의 시간 변화") + THEME
  })

  ## --------------------------------------------------------------- 4 brain
  output$p_nh3b <- renderPlot({
    ggplot(sim(), aes(time/24, NH3BR)) +
      geom_line(colour = PAL[1], linewidth = 1.0) +
      labs(x = "Days", y = "뇌 암모니아 (µmol/L)", title = NULL) + THEME
  })
  output$p_glnb <- renderPlot({
    ggplot(sim(), aes(time/24, GLNBRAIN)) +
      band(4.0, 6.0) + geom_line(colour = PAL[2], linewidth = 1.0) +
      labs(x = "Days", y = "뇌 글루타민 (mmol/L)", title = NULL) + THEME
  })
  output$p_mins <- renderPlot({
    ggplot(sim(), aes(time/24, MINSMM)) +
      band(4.0, 6.0) + geom_line(colour = PAL[4], linewidth = 1.0) +
      labs(x = "Days", y = "뇌 myo-inositol (mmol/L)",
           title = "소진된 완충능 = 다음 위기가 치명적인 이유") + THEME
  })
  output$p_bw <- renderPlot({
    d <- sim()
    ggplot(d, aes(time/24)) +
      geom_line(aes(y = BW), colour = PAL[2], linewidth = 1.0) +
      geom_line(aes(y = ICP/25), colour = PAL[1], linewidth = 0.8, linetype = 2) +
      scale_y_continuous(name = "뇌 수분 초과 (%)",
                         sec.axis = sec_axis(~.*25, name = "ICP (mmHg)")) +
      labs(x = "Days", title = "실선 = 뇌 수분, 점선 = ICP") + THEME
  })
  output$p_he <- renderPlot({
    ggplot(sim(), aes(time/24, HESTAGE)) +
      geom_line(colour = PAL[2], linewidth = 1.1) +
      scale_y_continuous(limits = c(0, 4)) +
      labs(x = "Days", y = "뇌병증 단계 (0-4)", title = NULL) + THEME
  })

  ## ----------------------------------------------------------- 5 endpoints
  output$vb_coma <- renderValueBox({
    v <- round(sim()$COMAH[nrow(sim())], 1)
    valueBox(paste0(v, " h"), "혼수 범위 (>200 µmol/L) 누적", icon = icon("bed"),
             color = if (v < 1) "green" else if (v < 24) "yellow" else "red")
  })
  output$vb_iq <- renderValueBox({
    v <- round(sim()$IQEST[nrow(sim())], 1)
    valueBox(v, "모델 예측 전체 IQ", icon = icon("brain"),
             color = if (v > 90) "green" else if (v > 75) "yellow" else "red")
  })
  output$vb_peak <- renderValueBox({
    v <- round(max(sim()$NH3), 0)
    valueBox(paste0(v, " µmol/L"), "최고 암모니아", icon = icon("arrow-up"),
             color = if (v < 100) "green" else if (v < 360) "yellow" else "red")
  })
  output$vb_auc <- renderValueBox({
    v <- round(sim()$NAUC[nrow(sim())], 0)
    valueBox(format(v, big.mark = ","), "AUC>100 (µmol/L·h)",
             icon = icon("chart-area"),
             color = if (v < 500) "green" else if (v < 5000) "yellow" else "red")
  })
  output$p_coma <- renderPlot({
    ggplot(sim(), aes(time/24, COMAH)) +
      geom_line(colour = PAL[2], linewidth = 1.0) +
      labs(x = "Days", y = "누적 시간 (h)", title = NULL) + THEME
  })
  output$p_iq <- renderPlot({
    d <- sim()
    ggplot(d, aes(time/24)) +
      geom_line(aes(y = IQEST), colour = PAL[1], linewidth = 1.0) +
      geom_line(aes(y = 100 - NEURO), colour = PAL[2], linewidth = 0.7,
                linetype = 2) +
      labs(x = "Days", y = "IQ (실선) / 100−손상지수 (점선)", title = NULL) + THEME
  })
  output$p_auc <- renderPlot({
    ggplot(sim(), aes(time/24, NAUC)) +
      geom_line(colour = PAL[4], linewidth = 1.0) +
      labs(x = "Days", y = "µmol/L·h", title = NULL) + THEME
  })
  output$p_tol <- renderPlot({
    ggplot(sim(), aes(time/24, PROTTOL)) +
      geom_line(colour = PAL[3], linewidth = 1.0) +
      geom_hline(yintercept = input$prot, linetype = 2, colour = PAL[2]) +
      labs(x = "Days", y = "g/kg/day",
           title = "점선 = 현재 처방된 섭취량") + THEME
  })

  ## ------------------------------------------------------------- 6 compare
  scen <- eventReactive(input$goscen, {
    req(length(input$scen) > 0)
    UCD_simulate_scenarios(MOD, which = as.integer(input$scen))
  }, ignoreNULL = FALSE)

  cplot <- function(v, lab, logy = FALSE) {
    renderPlot({
      d <- scen(); req(nrow(d) > 0)
      g <- ggplot(d, aes(time/24, .data[[v]], colour = scenario)) +
        geom_line(linewidth = 0.8) +
        labs(x = "Days", y = lab, title = NULL) + THEME +
        theme(legend.text = element_text(size = 7))
      if (logy) g <- g + scale_y_log10()
      g
    })
  }
  output$c_nh3  <- cplot("NH3",      "혈장 암모니아 (µmol/L)", TRUE)
  output$c_gln  <- cplot("GLNP",     "혈장 글루타민 (µmol/L)")
  output$c_tol  <- cplot("PROTTOL",  "단백 내성 (g/kg/day)")
  output$c_glnb <- cplot("GLNBRAIN", "뇌 글루타민 (mmol/L)")
  output$c_tab  <- renderDT({
    d <- scen(); req(nrow(d) > 0)
    datatable(UCD_summarise(d), options = list(scrollX = TRUE, pageLength = 14),
              rownames = FALSE)
  })

  ## ----------------------------------------------------------- 7 biomarker
  simple <- function(v, lab, hline = NULL, bandlo = NULL, bandhi = NULL,
                     col = PAL[1]) {
    renderPlot({
      g <- ggplot(sim(), aes(time/24, .data[[v]]))
      if (!is.null(bandlo)) g <- g + band(bandlo, bandhi)
      g <- g + geom_line(colour = col, linewidth = 0.9) +
        labs(x = "Days", y = lab, title = NULL) + THEME
      if (!is.null(hline)) g <- g + geom_hline(yintercept = hline, linetype = 2,
                                               colour = PAL[2])
      g
    })
  }
  output$b_upagn <- simple("UPAGN", "요중 PAGN (mmol/day)", col = PAL[3])
  output$b_ratio <- simple("PARATIO", "PAA : PAGN (µg/mL 기준)", hline = 2.5,
                           col = PAL[4])
  output$b_orot  <- simple("OROTC", "오로트산 풀 (µmol)", col = PAL[5])
  output$b_cit   <- simple("CITP", "시트룰린 (µmol/L)", bandlo = 20, bandhi = 45)
  output$b_arg   <- simple("ARGP", "아르기닌 (µmol/L)", bandlo = 50, bandhi = 130,
                           col = PAL[6])
  output$b_gly   <- simple("GLYC", "글리신 (µmol/L)", bandlo = 180, bandhi = 320,
                           col = PAL[7])
  output$b_bcaa  <- simple("BCAC", "BCAA (µmol/L)", bandlo = 300, bandhi = 500,
                           col = PAL[8])

  ## -------------------------------------------------------------- 8 safety
  output$s_paa <- renderPlot({
    ggplot(sim(), aes(time/24, PAAUG)) +
      band(0, 100) +
      geom_hline(yintercept = 500, linetype = 1, colour = PAL[2], linewidth = 0.8) +
      geom_line(colour = PAL[1], linewidth = 1.0) +
      labs(x = "Days", y = "혈장 PAA (µg/mL)",
           title = "적색선 = 500 µg/mL 신경독성 문턱") + THEME
  })
  output$s_paatox <- renderPlot({
    ggplot(sim(), aes(time/24, PAATOX)) +
      geom_line(colour = PAL[2], linewidth = 1.0) +
      labs(x = "Days", y = "µg/mL·h (>500 초과분)", title = NULL) + THEME
  })
  output$s_na <- renderPlot({
    ## NaPBA 186.2 g/mol and Na-benzoate 144.1 g/mol each carry 1 mol Na
    na <- input$napba/186.2*1000 + input$bz*input$wt/1000/144.11*1000
    df <- data.frame(source = c("NaPBA", "Na-benzoate"),
                     mmol = c(input$napba/186.2*1000,
                              input$bz*input$wt/1000/144.11*1000))
    ggplot(df, aes(source, mmol, fill = source)) +
      geom_col(width = 0.6) +
      geom_hline(yintercept = 100, linetype = 2, colour = PAL[2]) +
      scale_fill_manual(values = PAL[c(1, 4)]) +
      labs(x = NULL, y = "mmol Na/day",
           title = paste0("총 ", round(na), " mmol Na/day (점선 = 100)")) +
      THEME + theme(legend.position = "none")
  })
  output$s_rebound <- renderPlot({
    d <- sim()
    g <- ggplot(d, aes(time/24, NH3)) +
      geom_line(colour = PAL[1], linewidth = 1.0) +
      labs(x = "Days", y = "혈장 암모니아 (µmol/L)",
           title = "투석 종료 후 심부 글루타민에서 질소가 재분포") + THEME
    if (input$hd > 0) {
      g <- g + annotate("rect", xmin = input$hdt/24,
                        xmax = (input$hdt + input$hdd)/24,
                        ymin = -Inf, ymax = Inf, fill = PAL[3], alpha = 0.15)
    }
    g
  })
}

shinyApp(ui, server)
