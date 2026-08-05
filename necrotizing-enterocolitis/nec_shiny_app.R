## =============================================================================
##  nec_shiny_app.R
##  Necrotising Enterocolitis (NEC) — interactive QSP dashboard
##  신생아 괴사성 장염 QSP 대시보드
##
##  Requires nec_mrgsolve_model.R in the same directory.
##      shiny::runApp("nec_shiny_app.R")
##
##  DESIGN INTENT
##  -------------
##  The app is organised around the ONE STRUCTURAL CLAIM of the model, not
##  around the list of state variables.  Every tab answers one question:
##
##    ① 환자      who is this infant, and how far is their own threshold?
##    ② 핵심 고리  is the positive loop closing?  (E -> Pb -> Jtr -> INJ -> E)
##    ③ 장내 생태  which of the two ecological basins did the gut land in?
##    ④ 임계 부하  where are B_lo and B_hi for THIS infant, and where is B now?
##    ⑤ 병변·병기  has the lesion passed NECa_crit = INJth/wNEC?
##    ⑥ 임상 종점  Bell stage, surgery, growth, cholestasis, NDI index
##    ⑦ 시나리오   arm-vs-arm comparison under identical provocation
##    ⑧ 바이오마커 CRP / platelets / lactate / I-FABP-like readouts and their
##                 lead time relative to Bell II — the point being that they
##                 mostly do NOT lead
##    ⑨ 약동학     plasma AND luminal exposure, because only the luminal
##                 concentration reshapes the microbiome
##
##  DISCLAIMER: teaching / research tool.  Not for clinical use.
## =============================================================================

library(shiny)
library(mrgsolve)
library(dplyr)
library(tidyr)
library(ggplot2)

source("nec_mrgsolve_model.R", local = TRUE)

PAL <- c("#2f6fb0", "#c0392b", "#2e8b57", "#b8860b", "#6a51a3", "#0e7c86")

thin <- function(d, every = 4) d[seq(1, nrow(d), by = every), , drop = FALSE]

long_plot <- function(d, vars, labs_map, title, ylab = NULL, free = TRUE) {
  d %>%
    select(time, all_of(vars)) %>%
    pivot_longer(-time) %>%
    mutate(name = factor(labs_map[name], levels = unname(labs_map[vars]))) %>%
    ggplot(aes(time, value)) +
    geom_line(colour = PAL[1], linewidth = 0.65) +
    facet_wrap(~name, scales = if (free) "free_y" else "fixed") +
    labs(x = "출생 후 일수 (postnatal day)", y = ylab, title = title) +
    theme_bw(base_size = 11) +
    theme(strip.background = element_rect(fill = "#eef2f7"),
          panel.grid.minor = element_blank())
}

## =============================================================================
##  UI
## =============================================================================
ui <- fluidPage(
  tags$head(tags$style(HTML("
    .well { background:#f7f9fc; border-color:#dbe3ee; }
    h4 { color:#1f456f; margin-top:0; }
    .keynum { font-size:22px; font-weight:700; color:#1f456f; }
    .keylab { font-size:11px; color:#667; }
    .claim { background:#fff6f5; border-left:4px solid #c0392b;
             padding:8px 12px; margin-bottom:10px; font-size:12px; }
  "))),

  titlePanel("신생아 괴사성 장염 (NEC) — QSP 시뮬레이터 · Necrotising Enterocolitis QSP Simulator"),

  div(class = "claim",
      HTML("<b>이 모델의 구조적 주장:</b> NEC는 장세포 완전성 <i>E</i> 하나에
      걸린 닫힌 양성 되먹임 고리이며, 그 고리의 이득은 장내 병원균 부하
      <i>B</i>가 정한다. &nbsp;
      <code>E → Pb → Jtr → TLR4s → (세포소실↑, 증식↓) → E</code> &nbsp;
      <b>B_lo</b>에서 NEC가 <i>가능</i>해지고 <b>B_hi</b>에서 <i>필연</i>이 되며,
      그 사이에서는 부하가 아니라 <b>침습 삽화</b>가 운명을 정한다.")),

  sidebarLayout(
    sidebarPanel(
      width = 3,
      h4("① 환자 (Patient)"),
      sliderInput("GA", "재태연령 GA (주)", 23, 34, 27, step = 0.5),
      sliderInput("BW", "출생체중 (kg)", 0.40, 2.20, 0.92, step = 0.02),
      radioButtons("binfant", "HMO 대사 균주 보유 (B. infantis)",
                   c("보유 (carrier, ~35%)" = "1",
                     "미보유 (non-carrier, ~65%)" = "0.05"),
                   selected = "1"),
      sliderInput("dysbio", "이상균총 압력 (환경 유입 배수)",
                  0.5, 3.0, 1.5, step = 0.1),

      h4("② 영양 (Nutrition)"),
      selectInput("milk", "수유 종류",
                  c("친모 모유 (MOM)" = "MOM",
                    "기증모유 저온살균 (Donor)" = "DONOR",
                    "미숙아 분유 (Formula)" = "FORMULA",
                    "혼합 (Mixed)" = "MIXED"), selected = "MOM"),
      sliderInput("feedrate", "증량 속도 (mL/kg/d per day)", 0, 40, 20, step = 5),
      sliderInput("feedmax", "목표 수유량 (mL/kg/d)", 60, 200, 160, step = 10),
      checkboxInput("npo", "완전 금식 유지 (NPO + TPN)", FALSE),

      h4("③ 약물 (Drugs)"),
      sliderInput("abx", "경험적 암피실린+젠타마이신 (일)", 0, 14, 3, step = 1),
      sliderInput("mtz", "메트로니다졸 (일)", 0, 14, 0, step = 1),
      sliderInput("ind", "인도메타신 PDA 치료 (일)", 0, 3, 0, step = 1),
      sliderInput("ibu", "이부프로펜 PDA 치료 (일)", 0, 3, 0, step = 1),
      sliderInput("prob", "프로바이오틱스 일일 용량 (1e9)", 0, 2, 0, step = 0.1),
      sliderInput("thr", "문턱 상승 약제 (thr_boost)", 1.0, 2.5, 1.0, step = 0.1),

      h4("④ 침습 삽화 (Precipitating insult)"),
      sliderInput("insday", "깊은 저혈압 삽화 발생일", 3, 35, 14, step = 1),
      sliderInput("insdep", "삽화 깊이 (장혈류 하강분)", 0, 0.6, 0.40, step = 0.02),
      sliderInput("insdur", "삽화 지속 (일)", 0.1, 1.5, 0.6, step = 0.1),
      checkboxInput("tx", "수혈 2회 (10일 · 20일)", TRUE),

      h4("⑤ 치료 반응 (Response)"),
      sliderInput("dxlag", "진단–치료 지연 dx_lag (일)", 0, 3, 1.2, step = 0.1),
      checkboxInput("treat", "Bell ≥ II 에서 금식 + 광범위 항생제", TRUE),
      hr(),
      sliderInput("tend", "시뮬레이션 기간 (일)", 21, 70, 56, step = 7),
      helpText(HTML("<small>이 앱은 교육·연구 목적입니다. 임상 판단에
                     사용하지 마십시오.</small>"))
    ),

    mainPanel(
      width = 9,
      fluidRow(
        column(2, div(class = "keylab", "최고 Bell 병기"),
                  div(class = "keynum", textOutput("kBell", inline = TRUE))),
        column(2, div(class = "keylab", "확진일 (Bell II)"),
                  div(class = "keynum", textOutput("kOnset", inline = TRUE))),
        column(2, div(class = "keylab", "B_lo (가능)"),
                  div(class = "keynum", textOutput("kBlo", inline = TRUE))),
        column(2, div(class = "keylab", "B_hi (필연)"),
                  div(class = "keynum", textOutput("kBhi", inline = TRUE))),
        column(2, div(class = "keylab", "최대 B"),
                  div(class = "keynum", textOutput("kBmax", inline = TRUE))),
        column(2, div(class = "keylab", "회복불가점 통과"),
                  div(class = "keynum", textOutput("kCrit", inline = TRUE)))
      ),
      hr(),

      tabsetPanel(
        id = "tabs", type = "tabs",

        tabPanel("① 환자 프로파일",
          br(),
          fluidRow(column(6, plotOutput("pMat", height = 300)),
                   column(6, plotOutput("pFeed", height = 300))),
          h4("성숙도가 정하는 것들 (what maturation sets)"),
          tableOutput("tMat"),
          helpText(HTML("교정연령이 올라가면 TLR4 과발현이 줄고, PAF-AH·IL-10·
          디펜신·밀착연접이 함께 올라옵니다. 이 표의 값들이 <b>같은 균 부하에
          대해 이 아이가 얼마나 크게 반응할지</b>를 정합니다."))
        ),

        tabPanel("② 핵심 되먹임 고리",
          br(),
          plotOutput("pLoop", height = 620),
          helpText(HTML("위에서 아래로 고리를 한 바퀴 읽으십시오:
          <code>E</code> → 투과성 <code>Pb</code> → 전위 플럭스 <code>Jtr</code>
          → NF-κB → 총 손상 <code>INJ</code> → 다시 <code>E</code>.
          손상이 <code>INJth</code>를 넘는 구간에서만 병변이 자랍니다."))
        ),

        tabPanel("③ 장내 생태",
          br(),
          plotOutput("pEco", height = 460),
          h4("어느 basin에 떨어졌는가"),
          tableOutput("tEco"),
          helpText(HTML("HMO는 상재균만 쓸 수 있는 <b>사유 기질</b>입니다.
          그래서 HMO를 대사할 균주가 없으면(<i>non-carrier</i>) 모유를 먹여도
          Enterobacteriaceae 우세 상태로 떨어질 수 있습니다."))
        ),

        tabPanel("④ 임계 부하 · 분수령",
          br(),
          fluidRow(column(7, plotOutput("pField", height = 420)),
                   column(5, plotOutput("pWindow", height = 420))),
          helpText(HTML("왼쪽: 축소된 1차원 장(場) <code>g(E)</code>. 근(root)이
          곧 운명이며, 가운데 불안정근이 <b>분수령 E*</b>입니다.
          오른쪽: 이 아이의 <code>B_lo</code>–<code>B_hi</code> 창과 실제 부하의
          궤적. 창 안에 있으면 부하가 아니라 삽화가 결정합니다."))
        ),

        tabPanel("⑤ 병변 · 병기",
          br(),
          plotOutput("pLesion", height = 480),
          helpText(HTML("괴사 조직 자체가 염증 자극(DAMPs)이므로 병변은 자기가
          만든 손상을 유지합니다. 그 결과 <b>회복 불가점
          NECa_crit = INJth / wNEC</b>이 존재하고, 진단–치료 지연
          <code>dx_lag</code> 동안 이 선을 넘느냐가 내과적 NEC와 수술적 NEC를
          가릅니다. 장벽내 공기(pneumatosis)는 점막이 이미 뚫린 뒤에만
          생깁니다 (<code>Pb &gt; PbPNE</code>)."))
        ),

        tabPanel("⑥ 임상 종점",
          br(),
          fluidRow(column(6, plotOutput("pBell", height = 300)),
                   column(6, plotOutput("pOrgan", height = 300))),
          h4("종점 요약"),
          tableOutput("tEnd")
        ),

        tabPanel("⑦ 시나리오 비교",
          br(),
          checkboxGroupInput("arms", "비교할 아암 (동일한 침습 삽화 조건)",
            choices = c("친모 모유" = "MOM",
                        "기증모유" = "DONOR",
                        "분유" = "FORMULA",
                        "분유 + 프로바이오틱스" = "FORM_PROB",
                        "분유 + 문턱 상승 약제" = "FORM_THR",
                        "분유 + 둘 다" = "FORM_BOTH",
                        "완전 금식 + TPN" = "NPO"),
            selected = c("MOM", "FORMULA", "FORM_PROB", "FORM_BOTH"),
            inline = TRUE),
          actionButton("go", "시나리오 실행", class = "btn-primary"),
          br(), br(),
          plotOutput("pArms", height = 420),
          tableOutput("tArms"),
          helpText(HTML("<b>문턱을 옮기는 약</b>과 <b>상태를 옮기는 약</b>은
          더해지지 않고 곱해집니다. 모유는 둘 다이므로, 단일표적 약 하나가 아니라
          <b>둘을 합친 것</b>이 공정한 비교 대상입니다."))
        ),

        tabPanel("⑧ 바이오마커 · 선행 시간",
          br(),
          plotOutput("pBio", height = 420),
          h4("Bell II 대비 선행 시간 (음수 = 뒤늦게 움직임)"),
          tableOutput("tLead"),
          helpText(HTML("여기서 볼 것은 표지자가 <b>잘 오른다</b>는 사실이 아니라,
          대부분 <b>스위치가 넘어간 뒤에</b> 오른다는 사실입니다. 앞서 움직이는
          것은 상태변수(<code>E</code>, <code>Jtr</code>)이며, 그것이 예방이
          조기진단보다 앞서는 구조적 이유입니다."))
        ),

        tabPanel("⑨ 약동학 (혈장 vs 장내)",
          br(),
          plotOutput("pPK", height = 460),
          h4("장내 노출 (미생물이 실제로 보는 농도)"),
          tableOutput("tPK"),
          helpText(HTML("미생물 군집을 바꾸는 것은 혈장 농도가 아니라
          <b>장내 농도</b>입니다. 암피실린(f<sub>lum</sub>≈0.45)과
          메트로니다졸(≈0.80)은 장내에 도달하지만 정맥 아미노글리코사이드
          (≈0.05)는 거의 도달하지 않습니다. 경험적 암피실린+젠타마이신이
          혐기균만 초토화하고 Enterobacteriaceae는 비교적 온존시키는 비대칭이
          여기서 나옵니다."))
        )
      )
    )
  )
)

## =============================================================================
##  SERVER
## =============================================================================
server <- function(input, output, session) {

  ## ---- insult schedule from the sidebar ----------------------------------
  insults <- reactive({
    z <- list(INS1A = input$insday, INS1B = input$insday + input$insdur,
              INS1DEP = input$insdep, INS1NFK = 0.010)
    if (isTRUE(input$tx)) {
      z <- c(z, list(INS2A = 10.0, INS2B = 10.5, INS2DEP = 0.16, INS2NFK = 0.030,
                     INS3A = 20.0, INS3B = 20.5, INS3DEP = 0.16, INS3NFK = 0.030))
    } else {
      z <- c(z, list(INS2A = 1e6, INS2B = 1e6, INS2DEP = 0, INS2NFK = 0,
                     INS3A = 1e6, INS3B = 1e6, INS3DEP = 0, INS3NFK = 0))
    }
    z
  })

  ## ---- the index simulation ---------------------------------------------
  sim <- reactive({
    run_nec(GA = input$GA, BW = input$BW, milk = input$milk,
            feedrate = input$feedrate, feedmax = input$feedmax,
            abx_days = input$abx, mtz_days = input$mtz,
            ind_days = input$ind, ibu_days = input$ibu,
            prob_dose = input$prob, thrboost = input$thr,
            dysbio = input$dysbio,
            binfant = as.numeric(input$binfant),
            NPO_ON = as.numeric(isTRUE(input$npo)),
            dx_lag = input$dxlag,
            npo_on_nec = as.numeric(isTRUE(input$treat)),
            insults = insults(), end = input$tend)
  })

  onset  <- reactive({ d <- sim(); if (any(d$BellStage >= 2))
                       min(d$time[d$BellStage >= 2]) else NA_real_ })
  critNa <- reactive({ first(sim()$NECa_crit) })

  loads <- reactive({
    pp  <- par_list()
    pp$Ktlr0 <- pp$Ktlr0 * input$thr
    ctx <- nec_context(pp, milk = input$milk, GA = input$GA,
                       feed = input$feedmax, probiotic = input$prob > 0)
    bifurcation_loads(pp, ctx)
  })

  ## ---- key numbers -------------------------------------------------------
  output$kBell  <- renderText(as.character(max(sim()$BellStage)))
  output$kOnset <- renderText(ifelse(is.na(onset()), "—",
                                     sprintf("d%.1f", onset())))
  output$kBlo   <- renderText(sprintf("%.1f", loads()[["B_lo"]]))
  output$kBhi   <- renderText(sprintf("%.1f", loads()[["B_hi"]]))
  output$kBmax  <- renderText(sprintf("%.1f", max(sim()$B)))
  output$kCrit  <- renderText(ifelse(any(sim()$NECA > critNa()), "예 (YES)", "아니오"))

  ## ---- ① patient ---------------------------------------------------------
  output$pMat <- renderPlot({
    d <- thin(sim())
    pp <- par_list()
    d %>% mutate(TLR4expr = pp$TLR4max*(1 - pp$phiTLR*Maturation),
                 PAF_AH   = pp$dPAF0 + pp$dPAFm*Maturation) %>%
      select(time, Maturation, TLR4expr, PAF_AH) %>%
      pivot_longer(-time) %>%
      ggplot(aes(time, value, colour = name)) +
      geom_line(linewidth = 0.7) +
      scale_colour_manual(values = PAL[1:3], name = NULL) +
      labs(x = "출생 후 일수", y = NULL, title = "성숙도와 그것이 정하는 계수") +
      theme_bw(base_size = 11)
  })

  output$pFeed <- renderPlot({
    d <- thin(sim())
    d %>% select(time, FeedVol, SUB, Ischaemia) %>%
      pivot_longer(-time) %>%
      ggplot(aes(time, value)) + geom_line(colour = PAL[4], linewidth = 0.7) +
      facet_wrap(~name, scales = "free_y", ncol = 1) +
      labs(x = "출생 후 일수", y = NULL,
           title = "수유는 양쪽 부호: 영양 · 기질 · 산소요구") +
      theme_bw(base_size = 11)
  })

  output$tMat <- renderTable({
    d <- sim(); pp <- par_list()
    m0 <- first(d$Maturation); m1 <- last(d$Maturation)
    data.frame(
      항목 = c("성숙도 MAT", "TLR4 과발현 계수", "PAF-AH 용량",
               "신호 문턱 Ktlr", "영양 자극 Ftroph(최대)"),
      시작 = c(sprintf("%.3f", m0),
               sprintf("%.3f", pp$TLR4max*(1 - pp$phiTLR*m0)),
               sprintf("%.2f", pp$dPAF0 + pp$dPAFm*m0),
               sprintf("%.1f", first(d$Threshold)), "—"),
      종료 = c(sprintf("%.3f", m1),
               sprintf("%.3f", pp$TLR4max*(1 - pp$phiTLR*m1)),
               sprintf("%.2f", pp$dPAF0 + pp$dPAFm*m1),
               sprintf("%.1f", max(d$Threshold)), "—"))
  })

  ## ---- ② the loop --------------------------------------------------------
  output$pLoop <- renderPlot({
    d <- thin(sim())
    lm <- c(E = "① 장세포 완전성 E", Permeability = "② 투과성 Pb",
            Translocat = "③ 전위 플럭스 Jtr", TLRS = "④ NF-κB 활성",
            Injury = "⑤ 총 손상 INJ", NECA = "⑥ 전층 괴사 NECa")
    pp <- par_list()
    p <- d %>% select(time, E, Permeability, Translocat, TLRS, Injury, NECA) %>%
      pivot_longer(-time) %>%
      mutate(name = factor(lm[name], levels = unname(lm))) %>%
      ggplot(aes(time, value)) +
      geom_line(colour = PAL[2], linewidth = 0.7) +
      facet_wrap(~name, scales = "free_y", ncol = 1) +
      labs(x = "출생 후 일수", y = NULL,
           title = "닫힌 양성 되먹임 고리를 한 바퀴") +
      theme_bw(base_size = 11)
    if (!is.na(onset())) p <- p + geom_vline(xintercept = onset(),
                                            linetype = 2, colour = "#888")
    p
  })

  ## ---- ③ ecology --------------------------------------------------------
  output$pEco <- renderPlot({
    d <- thin(sim())
    lm <- c(B = "병원균 B (Enterobacteriaceae)", C = "상재 혐기균 C",
            P = "프로바이오틱스 P", SCFA = "SCFA (mM)",
            HMO = "HMO 풀 (g/kg)", GAS = "장내 가스")
    long_plot(d, c("B", "C", "P", "SCFA", "HMO", "GAS"), lm,
              "장내 생태 — 두 basin 중 어디에 떨어졌는가")
  })

  output$tEco <- renderTable({
    d <- sim(); L <- loads()
    data.frame(
      지표 = c("최종 B", "최종 C", "최종 SCFA", "B_lo", "B_hi",
               "최대 B의 창 내 위치"),
      값 = c(sprintf("%.2f", last(d$B)), sprintf("%.2f", last(d$C)),
             sprintf("%.1f mM", last(d$SCFA)),
             sprintf("%.1f", L[["B_lo"]]), sprintf("%.1f", L[["B_hi"]]),
             sprintf("%.0f %%", 100*(max(d$B) - L[["B_lo"]]) /
                                (L[["B_hi"]] - L[["B_lo"]]))))
  })

  ## ---- ④ field and window ------------------------------------------------
  output$pField <- renderPlot({
    pp <- par_list(); pp$Ktlr0 <- pp$Ktlr0 * input$thr
    ctx <- nec_context(pp, milk = input$milk, GA = input$GA,
                       feed = input$feedmax, probiotic = input$prob > 0)
    L <- loads()
    Bs <- round(c(0.5*L[["B_lo"]], L[["B_lo"]],
                  0.5*(L[["B_lo"]] + L[["B_hi"]]), L[["B_hi"]]), 1)
    grid <- expand.grid(E = seq(0.02, 0.99, by = 0.005), B = Bs)
    grid$g <- mapply(reduced_field, grid$E, grid$B,
                     MoreArgs = list(pp = pp, ctx = ctx))
    ggplot(grid, aes(E, g, colour = factor(B))) +
      geom_hline(yintercept = 0, linetype = 2, colour = "#777") +
      geom_line(linewidth = 0.75) +
      scale_colour_manual(values = PAL, name = "B (1e9/g)") +
      labs(x = "장세포 완전성 E", y = "g(E) = dE/dt",
           title = "축소된 1차원 장: 근(root)이 곧 운명") +
      theme_bw(base_size = 11)
  })

  output$pWindow <- renderPlot({
    d <- thin(sim()); L <- loads()
    ggplot(d, aes(time, B)) +
      annotate("rect", xmin = -Inf, xmax = Inf,
               ymin = L[["B_lo"]], ymax = L[["B_hi"]],
               fill = "#f6d9d5", alpha = 0.7) +
      annotate("text", x = max(d$time)*0.5,
               y = (L[["B_lo"]] + L[["B_hi"]])/2,
               label = "이중안정 창: 삽화가 운명을 정한다",
               size = 3.4, colour = "#8a2418") +
      geom_hline(yintercept = L[["B_lo"]], linetype = 2, colour = "#c0392b") +
      geom_hline(yintercept = L[["B_hi"]], linetype = 2, colour = "#8a2418") +
      geom_line(colour = PAL[1], linewidth = 0.8) +
      labs(x = "출생 후 일수", y = "병원균 부하 B (1e9 CFU/g)",
           title = "이 아이의 임계 부하와 실제 부하") +
      theme_bw(base_size = 11)
  })

  ## ---- ⑤ lesion ---------------------------------------------------------
  output$pLesion <- renderPlot({
    d <- thin(sim()); cr <- critNa()
    lm <- c(Injury = "총 손상 INJ", NECA = "전층 괴사 NECa",
            PNEU = "장벽내 공기 PNEU", PERF = "장혈류 PERF")
    p <- d %>% select(time, Injury, NECA, PNEU, PERF) %>%
      pivot_longer(-time) %>%
      mutate(name = factor(lm[name], levels = unname(lm))) %>%
      ggplot(aes(time, value)) +
      geom_line(colour = PAL[2], linewidth = 0.7) +
      facet_wrap(~name, scales = "free_y") +
      labs(x = "출생 후 일수", y = NULL,
           title = sprintf("병변과 회복 불가점 (NECa_crit = %.3f)", cr)) +
      theme_bw(base_size = 11)
    p
  })

  ## ---- ⑥ endpoints ------------------------------------------------------
  output$pBell <- renderPlot({
    d <- thin(sim())
    ggplot(d, aes(time, BellStage)) +
      geom_step(colour = PAL[2], linewidth = 0.8) +
      scale_y_continuous(breaks = 0:3, limits = c(0, 3)) +
      labs(x = "출생 후 일수", y = "Bell 병기", title = "임상 병기 경과") +
      theme_bw(base_size = 11)
  })

  output$pOrgan <- renderPlot({
    d <- thin(sim())
    lm <- c(BILI = "담즙정체 (mg/dL)", NIN = "신경염증 지수",
            KIN = "세뇨관 손상 지수", WtGain_g = "체중 증가 (g)")
    long_plot(d, c("BILI", "NIN", "KIN", "WtGain_g"), lm, "장 외 결과")
  })

  output$tEnd <- renderTable({
    d <- sim()
    data.frame(
      종점 = c("최고 Bell 병기", "확진일 (Bell II)", "Bell III 도달",
               "회복 불가점 통과", "최저 E", "최종 담즙 빌리루빈",
               "체중 증가", "신경염증 지수", "세뇨관 손상 지수"),
      값 = c(as.character(max(d$BellStage)),
             ifelse(is.na(onset()), "—", sprintf("d%.1f", onset())),
             ifelse(any(d$BellStage >= 3),
                    sprintf("d%.1f", min(d$time[d$BellStage >= 3])), "—"),
             ifelse(any(d$NECA > critNa()), "예", "아니오"),
             sprintf("%.3f", min(d$E)),
             sprintf("%.2f mg/dL", last(d$BILI)),
             sprintf("%.0f g", last(d$WtGain_g)),
             sprintf("%.3f", last(d$NIN)),
             sprintf("%.3f", last(d$KIN))))
  })

  ## ---- ⑦ scenario comparison --------------------------------------------
  armdefs <- list(
    MOM       = list(label = "친모 모유", args = list(milk = "MOM")),
    DONOR     = list(label = "기증모유", args = list(milk = "DONOR")),
    FORMULA   = list(label = "분유", args = list(milk = "FORMULA")),
    FORM_PROB = list(label = "분유 + 프로바이오틱스",
                     args = list(milk = "FORMULA", prob_dose = 1.10)),
    FORM_THR  = list(label = "분유 + 문턱 상승 약제",
                     args = list(milk = "FORMULA", thrboost = 1.6)),
    FORM_BOTH = list(label = "분유 + 둘 다",
                     args = list(milk = "FORMULA", prob_dose = 1.10,
                                 thrboost = 1.6)),
    NPO       = list(label = "완전 금식 + TPN",
                     args = list(milk = "MOM", NPO_ON = 1))
  )

  arms_sim <- eventReactive(input$go, {
    req(length(input$arms) > 0)
    withProgress(message = "시나리오 실행 중...", {
      dplyr::bind_rows(lapply(input$arms, function(k) {
        incProgress(1/length(input$arms), detail = armdefs[[k]]$label)
        a <- c(list(GA = input$GA, BW = input$BW, feedrate = input$feedrate,
                    feedmax = input$feedmax, abx_days = input$abx,
                    mtz_days = input$mtz, ind_days = input$ind,
                    dysbio = input$dysbio,
                    binfant = as.numeric(input$binfant),
                    dx_lag = input$dxlag,
                    insults = insults(), end = input$tend),
               armdefs[[k]]$args)
        do.call(run_nec, a) %>% mutate(arm = armdefs[[k]]$label)
      }))
    })
  }, ignoreNULL = FALSE)

  output$pArms <- renderPlot({
    d <- arms_sim(); req(nrow(d) > 0)
    d <- d[seq(1, nrow(d), by = 4), ]
    lm <- c(E = "장세포 완전성 E", B = "병원균 부하 B",
            Injury = "총 손상 INJ", NECA = "전층 괴사 NECa")
    d %>% select(time, arm, E, B, Injury, NECA) %>%
      pivot_longer(-c(time, arm)) %>%
      mutate(name = factor(lm[name], levels = unname(lm))) %>%
      ggplot(aes(time, value, colour = arm)) +
      geom_line(linewidth = 0.7) +
      facet_wrap(~name, scales = "free_y") +
      scale_colour_manual(values = PAL, name = NULL) +
      labs(x = "출생 후 일수", y = NULL,
           title = "동일한 침습 삽화 아래 아암 비교") +
      theme_bw(base_size = 11) + theme(legend.position = "bottom")
  })

  output$tArms <- renderTable({
    d <- arms_sim(); req(nrow(d) > 0)
    d %>% group_by(arm) %>%
      summarise(`최고 Bell` = max(BellStage),
                `확진일` = ifelse(any(BellStage >= 2),
                                  sprintf("d%.1f", min(time[BellStage >= 2])), "—"),
                `회복불가점 통과` = ifelse(any(NECA > first(NECa_crit)), "예", "아니오"),
                `최저 E` = sprintf("%.3f", min(E)),
                `최대 B` = sprintf("%.1f", max(B)),
                `최종 SCFA` = sprintf("%.1f", last(SCFA)),
                `체중 증가 g` = sprintf("%.0f", last(WtGain_g)),
                .groups = "drop")
  })

  ## ---- ⑧ biomarkers -----------------------------------------------------
  output$pBio <- renderPlot({
    d <- thin(sim())
    lm <- c(CRP = "CRP (mg/L)", PLT = "혈소판 (1e9/L)", LAC = "젖산 (mmol/L)",
            LPSP = "내독소 지수", PNEU = "장벽내 공기", Translocat = "Jtr (기전량)")
    p <- long_plot(d, c("CRP", "PLT", "LAC", "LPSP", "PNEU", "Translocat"), lm,
                   "임상 표지자와 기전량")
    if (!is.na(onset())) p <- p + geom_vline(xintercept = onset(),
                                            linetype = 2, colour = "#c0392b")
    p
  })

  output$tLead <- renderTable({
    d <- sim(); o <- onset()
    cross <- function(v, thr, below = FALSE) {
      hit <- if (below) which(v < thr) else which(v > thr)
      if (!length(hit)) NA_real_ else d$time[hit[1]]
    }
    mk <- list("CRP > 12"      = cross(d$CRP, 12),
               "혈소판 < 100"  = cross(d$PLT, 100, TRUE),
               "젖산 > 4"      = cross(d$LAC, 4),
               "장벽내 공기"   = cross(d$PNEU, 0.30),
               "Jtr > 3 (기전량)" = cross(d$Translocat, 3))
    data.frame(
      표지자 = names(mk),
      교차일 = vapply(mk, function(x) ifelse(is.na(x), "—",
                                             sprintf("d%.2f", x)), ""),
      `선행시간_시간` = vapply(mk, function(x)
        ifelse(is.na(x) || is.na(o), "—", sprintf("%+.1f", (o - x)*24)), ""),
      row.names = NULL)
  })

  ## ---- ⑨ PK -------------------------------------------------------------
  output$pPK <- renderPlot({
    d <- thin(sim(), 2); pp <- par_list()
    d %>% transmute(time,
                    `암피실린 혈장` = AMPC,
                    `암피실린 장내` = AMPC*pp$fl_amp,
                    `젠타마이신 혈장` = GENC,
                    `젠타마이신 장내` = GENC*pp$fl_gen,
                    `메트로니다졸 혈장` = MTZC,
                    `메트로니다졸 장내` = MTZC*pp$fl_mtz) %>%
      pivot_longer(-time) %>%
      mutate(drug = sub(" .*", "", name),
             site = ifelse(grepl("장내", name), "장내 (lumen)", "혈장 (plasma)")) %>%
      ggplot(aes(time, value, colour = site)) +
      geom_line(linewidth = 0.65) +
      facet_wrap(~drug, scales = "free_y") +
      scale_colour_manual(values = c(`장내 (lumen)` = PAL[2],
                                     `혈장 (plasma)` = PAL[1]), name = NULL) +
      labs(x = "출생 후 일수", y = "농도 (mg/L)",
           title = "혈장과 장내는 다른 노출이다") +
      theme_bw(base_size = 11) + theme(legend.position = "bottom")
  })

  output$tPK <- renderTable({
    d <- sim(); pp <- par_list()
    data.frame(
      약물 = c("암피실린", "젠타마이신", "메트로니다졸",
               "인도메타신", "이부프로펜"),
      `f_lum` = c(sprintf("%.2f", pp$fl_amp), sprintf("%.2f", pp$fl_gen),
                  sprintf("%.2f", pp$fl_mtz), "—", "—"),
      `혈장 Cmax` = c(sprintf("%.1f", max(d$AMPC)), sprintf("%.2f", max(d$GENC)),
                      sprintf("%.1f", max(d$MTZC)), sprintf("%.3f", max(d$INDC)),
                      sprintf("%.2f", max(d$IBUC))),
      `장내 Cmax` = c(sprintf("%.1f", max(d$AMPC)*pp$fl_amp),
                      sprintf("%.2f", max(d$GENC)*pp$fl_gen),
                      sprintf("%.1f", max(d$MTZC)*pp$fl_mtz), "—", "—"),
      `주 표적` = c("상재 혐기균 C (해로운 쪽)", "병원균 B (이로운 쪽)",
                    "상재 혐기균 C (해로운 쪽)", "프로스타노이드 → 장혈류",
                    "프로스타노이드 (영향 더 작음)"))
  })
}

shinyApp(ui, server)
