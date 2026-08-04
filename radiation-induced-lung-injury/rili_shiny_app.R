## =====================================================================
##  rili_shiny_app.R
##  Radiation-Induced Lung Injury (RILI) — interactive QSP dashboard
##  방사선 유발 폐손상 QSP 대시보드
##
##  The app is built around the one thing this model does that a
##  mean-dose model cannot: the treatment plan is a DOSE-VOLUME
##  HISTOGRAM, so the user edits a vector of six volumes, not a dose.
##  Tab 2 exists to let you hold the mean lung dose fixed while you
##  reshape the histogram, and watch pneumonitis and fibrosis move in
##  opposite directions.
##
##  Ten tabs
##    1  환자·플랜 (Patient & Plan)          DVH editor, MLD/V5/V20/V40
##    2  등가 평균선량 실험 (Iso-MLD)         same MLD, different shape
##    3  약동학 (Pharmacokinetics)            8 drugs
##    4  빠른 고리 (Fast Loop)                DAMP → cytokine → oedema
##    5  느린 고리 (Slow Loop)                TGF-β → myofibroblast → collagen
##    6  이중안정 스위치 (Bistable Switch)    phase plane + separatrix
##    7  폐기능 (Pulmonary Function)          DLCO, FVC, per-bin map
##    8  임상 엔드포인트 (Clinical Endpoints)  NTCP, CTCAE grade, latency
##    9  치료비 (Therapeutic Ratio)           TCP vs NTCP vs UCP
##   10  시나리오 비교 (Scenario Comparison)   the 40 canned scenarios
##
##  RUN
##    shiny::runApp("rili_shiny_app.R")
##  The model file rili_mrgsolve_model.R must sit in the same directory.
## =====================================================================

library(shiny)
library(mrgsolve)
library(dplyr)
library(tidyr)
library(ggplot2)

MODEL_DIR <- "."
mod <- mread_cache("rili_mrgsolve_model", MODEL_DIR)

## ---- plan library ---------------------------------------------------
PLANS <- list(
  "IMRT 60 Gy / 30 fx"       = list(DRX = 60, NFX = 30,
                                    v = c(.40, .20, .12, .11, .09, .08)),
  "3D-CRT 60 Gy / 30 fx"     = list(DRX = 60, NFX = 30,
                                    v = c(.52, .13, .08, .08, .09, .10)),
  "Proton 60 Gy / 30 fx"     = list(DRX = 60, NFX = 30,
                                    v = c(.58, .16, .09, .07, .055, .045)),
  "IMRT 74 Gy / 37 fx"       = list(DRX = 74, NFX = 37,
                                    v = c(.34, .20, .13, .12, .11, .10)),
  "Hypofx 60 Gy / 8 fx"      = list(DRX = 60, NFX = 8,
                                    v = c(.700, .130, .070, .045, .033,
                                          .022)),
  "SBRT 54 Gy / 3 fx"        = list(DRX = 54, NFX = 3,
                                    v = c(.860, .075, .032, .018, .010,
                                          .005)),
  "Low-dose bath (V5 high)"  = list(DRX = 60, NFX = 30,
                                    v = c(.06, .55, .29, .05, .03, .02)),
  "Focal hot spot (V40 high)" = list(DRX = 60, NFX = 30,
                                     v = c(.600, .110, .050, .050, .070,
                                           .120))
)
EDGEF <- c(0, .08, .24, .44, .66, .88, 1.05)
REPF  <- (EDGEF[-7] + EDGEF[-1]) / 2

mld_of <- function(v, DRX) sum(v / sum(v) * DRX * REPF)
vx_of <- function(v, DRX, x) {
  vv <- v / sum(v)
  tot <- 0
  for (j in seq_len(6)) {
    lo <- DRX * EDGEF[j]; hi <- DRX * EDGEF[j + 1]
    if (hi <= x) next
    tot <- tot + if (lo >= x) vv[j] else vv[j] * (hi - x) / (hi - lo)
  }
  tot
}
geud_of <- function(v, DRX, a) {
  vv <- v / sum(v); sum(vv * (DRX * REPF) ^ a) ^ (1 / a)
}

THEME <- theme_minimal(base_size = 12) +
  theme(panel.grid.minor = element_blank(),
        legend.position = "bottom",
        plot.title = element_text(face = "bold", size = 12))

## =====================================================================
##  UI
## =====================================================================
ui <- fluidPage(
  titlePanel(HTML(paste0(
    "<b>방사선 유발 폐손상 QSP 모델</b> ",
    "<span style='font-size:70%'>Radiation-Induced Lung Injury &middot; ",
    "79 ODEs &middot; 6 DVH bins &times; 10 states</span>"))),
  sidebarLayout(
    sidebarPanel(
      width = 3,
      selectInput("plan", "플랜 프리셋 (Plan preset)",
                  choices = names(PLANS), selected = "IMRT 60 Gy / 30 fx"),
      fluidRow(
        column(6, numericInput("DRX", "처방선량 D_RX (Gy)", 60, 20, 90,
                               step = 2)),
        column(6, numericInput("NFX", "분할 횟수 NFX", 30, 1, 45, step = 1))
      ),
      tags$hr(),
      tags$b("선량-체적 히스토그램 (DVH)"),
      helpText(HTML(paste0(
        "구간별 폐 체적 분율. 합이 1이 아니어도 자동 정규화됩니다.<br>",
        "<i>이 여섯 개 숫자가 이 모델의 '처방'입니다.</i>"))),
      fluidRow(
        column(6, numericInput("v1", "bin1 (~4% D_RX)", .40, 0, 1,
                               step = .01)),
        column(6, numericInput("v2", "bin2 (~16%)", .20, 0, 1, step = .01))
      ),
      fluidRow(
        column(6, numericInput("v3", "bin3 (~34%)", .12, 0, 1, step = .01)),
        column(6, numericInput("v4", "bin4 (~55%)", .11, 0, 1, step = .01))
      ),
      fluidRow(
        column(6, numericInput("v5", "bin5 (~77%)", .09, 0, 1, step = .01)),
        column(6, numericInput("v6", "bin6 (~97%)", .08, 0, 1, step = .01))
      ),
      tags$hr(),
      tags$b("환자 예비능 (Patient reserve)"),
      sliderInput("DLCO0", "기저 DLCO (%pred)", 30, 110, 85, step = 5),
      sliderInput("FVC0", "기저 FVC (%pred)", 40, 120, 95, step = 5),
      checkboxInput("emph", "폐기종 관류 가중 (emphysema perfusion)",
                    FALSE),
      tags$hr(),
      tags$b("약물 (Drugs)"),
      tags$div(style = "font-size:88%",
        tags$span(style = "color:#B8860B", "● 빠른 고리 (fast loop)")),
      fluidRow(
        column(6, numericInput("rdex", "프레드니솔론 mg/d", 0, 0, 100,
                               step = 10)),
        column(6, numericInput("tdex", "시작일 (day)", 56, 0, 400,
                               step = 7))
      ),
      checkboxInput("durv", "더발루맙 1500 mg q4w (PACIFIC)", FALSE),
      tags$div(style = "font-size:88%",
        tags$span(style = "color:#9B4FA8", "● 느린 고리 (slow loop)")),
      fluidRow(
        column(6, numericInput("rpir", "피르페니돈 mg/d", 0, 0, 2403,
                               step = 801)),
        column(6, numericInput("tpir", "시작일 (day)", 0, 0, 400,
                               step = 7))
      ),
      fluidRow(
        column(6, numericInput("rnin", "닌테다닙 mg/d", 0, 0, 300,
                               step = 150)),
        column(6, numericInput("race", "리시노프릴 mg/d", 0, 0, 40,
                               step = 10))
      ),
      tags$div(style = "font-size:88%",
        tags$span(style = "color:#E08B2E", "● 초기 손상 (initial damage)")),
      checkboxInput("ami", "아미포스틴 (amifostine)", FALSE),
      sliderInput("amidel", "조사 전 간격 (min before beam)", 0, 90, 20,
                  step = 5),
      helpText(HTML(paste0("WR-1065 반감기 8분 — 보호계수는 ",
                           "exp(&minus;ln2&middot;&Delta;t/8) 로 감쇠합니다."))),
      tags$hr(),
      sliderInput("tend", "시뮬레이션 기간 (days)", 200, 1095, 730,
                  step = 30)
    ),
    mainPanel(
      width = 9,
      tabsetPanel(
        id = "tabs", type = "tabs",
        tabPanel("1 환자·플랜", br(),
                 fluidRow(column(7, plotOutput("dvhPlot", height = 300)),
                          column(5, tableOutput("planTab"))),
                 hr(), plotOutput("bedPlot", height = 260)),
        tabPanel("2 등가 평균선량 실험", br(),
                 htmlOutput("isoNote"),
                 plotOutput("isoPlot", height = 380),
                 hr(), tableOutput("isoTab")),
        tabPanel("3 약동학", br(),
                 plotOutput("pkPlot", height = 460),
                 hr(), htmlOutput("protNote")),
        tabPanel("4 빠른 고리", br(),
                 plotOutput("fastPlot", height = 500),
                 hr(), htmlOutput("gainNote")),
        tabPanel("5 느린 고리", br(),
                 plotOutput("slowPlot", height = 500)),
        tabPanel("6 이중안정 스위치", br(),
                 htmlOutput("switchNote"),
                 plotOutput("phasePlot", height = 420),
                 hr(), tableOutput("fpTab")),
        tabPanel("7 폐기능", br(),
                 plotOutput("pftPlot", height = 320),
                 hr(), plotOutput("binMap", height = 300)),
        tabPanel("8 임상 엔드포인트", br(),
                 fluidRow(column(6, plotOutput("ntcpPlot", height = 320)),
                          column(6, plotOutput("gradePlot", height = 320))),
                 hr(), tableOutput("endTab")),
        tabPanel("9 치료비", br(),
                 plotOutput("trPlot", height = 400),
                 hr(), tableOutput("trTab")),
        tabPanel("10 시나리오 비교", br(),
                 checkboxGroupInput(
                   "scen", "비교할 시나리오 (choose scenarios)",
                   choices = names(PLANS),
                   selected = c("SBRT 54 Gy / 3 fx", "IMRT 60 Gy / 30 fx",
                                "IMRT 74 Gy / 37 fx"),
                   inline = TRUE),
                 plotOutput("scenPlot", height = 420),
                 hr(), tableOutput("scenTab"))
      )
    )
  )
)

## =====================================================================
##  SERVER
## =====================================================================
server <- function(input, output, session) {

  ## keep the DVH boxes in sync with the preset selector
  observeEvent(input$plan, {
    p <- PLANS[[input$plan]]
    updateNumericInput(session, "DRX", value = p$DRX)
    updateNumericInput(session, "NFX", value = p$NFX)
    for (j in 1:6)
      updateNumericInput(session, paste0("v", j), value = p$v[j])
  })

  vvec <- reactive({
    v <- c(input$v1, input$v2, input$v3, input$v4, input$v5, input$v6)
    v[is.na(v) | v < 0] <- 0
    if (sum(v) <= 0) v <- c(1, 0, 0, 0, 0, 0)
    v / sum(v)
  })

  parlist <- reactive({
    v <- vvec()
    pf <- if (isTRUE(input$emph)) c(1, 1, .8, .6, .5, .4) else rep(1, 6)
    lst <- list(DRX = input$DRX, NFX = input$NFX,
                V1 = v[1], V2 = v[2], V3 = v[3],
                V4 = v[4], V5 = v[5], V6 = v[6],
                PF1 = pf[1], PF2 = pf[2], PF3 = pf[3],
                PF4 = pf[4], PF5 = pf[5], PF6 = pf[6],
                DLCO0 = input$DLCO0, FVC0 = input$FVC0,
                RDEX = input$rdex, TDEX0 = if (input$rdex > 0)
                  input$tdex else 1e9,
                RPIR = input$rpir, TPIR0 = if (input$rpir > 0)
                  input$tpir else 1e9,
                TPIR1 = if (input$rpir > 0) input$tpir + 365 else 1e9,
                RNIN = input$rnin, TNIN0 = if (input$rnin > 0) 0 else 1e9,
                TNIN1 = if (input$rnin > 0) 365 else 1e9,
                RACE = input$race, TACE0 = if (input$race > 0) 0 else 1e9,
                TACE1 = if (input$race > 0) 730 else 1e9,
                AMION = as.numeric(isTRUE(input$ami)),
                AMIDEL = input$amidel)
    lst
  })

  dosing <- reactive({
    if (isTRUE(input$durv))
      ev(amt = 1500, cmt = "DURC", time = 56, ii = 28, addl = 12)
    else ev(amt = 0, cmt = "DURC", time = 0)
  })

  sim <- reactive({
    mod %>% param(parlist()) %>% ev(dosing()) %>%
      mrgsim(end = input$tend, delta = 1) %>% as_tibble()
  })

  simplan <- function(nm, extra = list()) {
    p <- PLANS[[nm]]
    v <- p$v / sum(p$v)
    pr <- c(list(DRX = p$DRX, NFX = p$NFX, V1 = v[1], V2 = v[2],
                 V3 = v[3], V4 = v[4], V5 = v[5], V6 = v[6],
                 DLCO0 = input$DLCO0, FVC0 = input$FVC0), extra)
    mod %>% param(pr) %>% mrgsim(end = input$tend, delta = 2) %>%
      as_tibble() %>% mutate(plan = nm)
  }

  tcourse <- reactive(input$NFX / 5 * 7)

  ## ---- 1  plan ------------------------------------------------------
  output$dvhPlot <- renderPlot({
    v <- vvec(); D <- input$DRX
    d <- tibble(bin = factor(1:6), dose = D * REPF, vol = 100 * v)
    ggplot(d, aes(dose, vol)) +
      geom_col(fill = "#4A90D9", width = D * .08) +
      geom_text(aes(label = sprintf("%.0f%%", vol)), vjust = -.4,
                size = 3.4) +
      geom_vline(xintercept = 42, linetype = "dashed",
                 colour = "#9B4FA8", linewidth = .9) +
      annotate("text", x = 42, y = max(100 * v) * .95,
               label = "섬유화 문턱 42 Gy\n(fibrotic threshold)",
               hjust = -.05, size = 3.2, colour = "#9B4FA8") +
      labs(title = "선량-체적 히스토그램 (DVH)",
           subtitle = "보라색 선 = 이 모델이 산출한 섬유화 전환 문턱",
           x = "구간 대표선량 (Gy)", y = "폐 체적 (%)") + THEME
  })

  output$planTab <- renderTable({
    v <- vvec(); D <- input$DRX
    tibble(
      `지표 (metric)` = c("평균 폐선량 MLD (Gy)", "V5 (%)", "V20 (%)",
                          "V40 (%)", "gEUD a=1 (Gy)", "gEUD a=0.7 (Gy)",
                          "분할선량 d (Gy)", "치료 기간 (days)",
                          "종양 BED10 (Gy)", "최대 구간 폐 BED3 (Gy)"),
      `값 (value)` = c(
        sprintf("%.2f", mld_of(v, D)),
        sprintf("%.1f", 100 * vx_of(v, D, 5)),
        sprintf("%.1f", 100 * vx_of(v, D, 20)),
        sprintf("%.1f", 100 * vx_of(v, D, 40)),
        sprintf("%.2f", geud_of(v, D, 1)),
        sprintf("%.2f", geud_of(v, D, 0.7)),
        sprintf("%.2f", D / input$NFX),
        sprintf("%.0f", tcourse()),
        sprintf("%.1f", D * (1 + (D / input$NFX) / 10)),
        sprintf("%.1f", D * REPF[6] * (1 + (D * REPF[6] / input$NFX) / 3))))
  }, striped = TRUE)

  output$bedPlot <- renderPlot({
    v <- vvec(); D <- input$DRX; n <- input$NFX
    db <- D * REPF; dfx <- db / n
    d <- tibble(bin = factor(1:6), dose = db,
                BED3 = db * (1 + dfx / 3),
                BED10 = db * (1 + dfx / 10)) %>%
      pivot_longer(c(BED3, BED10))
    ggplot(d, aes(dose, value, colour = name)) +
      geom_line(linewidth = 1) + geom_point(size = 2.4) +
      scale_colour_manual(values = c(BED3 = "#C0392B", BED10 = "#4CA85F"),
                          labels = c("폐 BED₃ (α/β=3)",
                                     "종양 BED₁₀ (α/β=10)")) +
      labs(title = "같은 물리선량, 다른 생물학적 등가선량",
           subtitle = paste0("분할선량 ", sprintf("%.2f", D / n),
                             " Gy — 분할선량이 커질수록 두 곡선이 벌어진다"),
           x = "물리선량 (Gy)", y = "BED (Gy)", colour = NULL) + THEME
  })

  ## ---- 2  iso-MLD experiment ---------------------------------------
  output$isoNote <- renderUI({
    HTML(paste0(
      "<div style='background:#F2F8FD;padding:10px;border-left:4px solid ",
      "#4A90D9;margin-bottom:8px'><b>이 탭이 이 모델의 핵심 주장입니다.</b> ",
      "아래 두 플랜은 <b>평균 폐선량이 거의 같지만</b>(15.5 vs 15.4 Gy) ",
      "V40 이 3.8배 다릅니다. 폐렴과 섬유화는 두 플랜의 순위를 ",
      "<b>서로 반대로</b> 매깁니다 — 폐렴은 병렬 장기에 대한 체적 가중 ",
      "<i>합</i>이고, 섬유화는 국소 <i>문턱</i>이기 때문입니다.</div>"))
  })

  isoSim <- reactive({
    bind_rows(simplan("Low-dose bath (V5 high)"),
              simplan("Focal hot spot (V40 high)"))
  })

  output$isoPlot <- renderPlot({
    d <- isoSim() %>%
      select(time, plan, NTCP, VFIB, PNIE, FIB) %>%
      mutate(NTCP = 100 * NTCP, VFIB = 100 * VFIB) %>%
      pivot_longer(c(NTCP, VFIB, PNIE, FIB))
    lab <- c(NTCP = "NTCP grade≥2 (%)", VFIB = "섬유화 전환 체적 (%)",
             PNIE = "폐렴 지수 PNIe", FIB = "섬유화 점수 FIB")
    ggplot(d, aes(time, value, colour = plan)) +
      geom_line(linewidth = 1) +
      facet_wrap(~ name, scales = "free_y",
                 labeller = labeller(name = lab)) +
      scale_colour_manual(values = c("#4A90D9", "#C0392B")) +
      labs(title = "평균선량은 같고 히스토그램 모양만 다른 두 플랜",
           x = "시간 (days)", y = NULL, colour = NULL) + THEME
  })

  output$isoTab <- renderTable({
    isoSim() %>% group_by(plan) %>%
      summarise(MLD = first(MLD), V5 = 100 * first(V5G),
                V40 = 100 * first(V40G),
                `NTCP peak %` = 100 * max(NTCP),
                `Vfib final %` = 100 * last(VFIB),
                `DLCO final` = last(DLCO), .groups = "drop")
  }, digits = 2, striped = TRUE)

  ## ---- 3  PK --------------------------------------------------------
  output$pkPlot <- renderPlot({
    d <- sim() %>% select(time, CDEXo, CPIRo, CNINo, CDURo, CACEo) %>%
      pivot_longer(-time) %>% filter(value > 1e-9)
    if (!nrow(d)) {
      return(ggplot() + annotate("text", 0, 0, size = 5,
        label = "약물이 선택되지 않았습니다\n(no drug selected)") +
        theme_void())
    }
    lab <- c(CDEXo = "프레드니솔론 (mg/L)", CPIRo = "피르페니돈 (mg/L)",
             CNINo = "닌테다닙 (mg/L)", CDURo = "더발루맙 (mg/L)",
             CACEo = "리시노프릴 (mg/L)")
    ggplot(d, aes(time, value, colour = name)) +
      geom_line(linewidth = .9) +
      geom_vline(xintercept = tcourse(), linetype = "dotted") +
      facet_wrap(~ name, scales = "free_y",
                 labeller = labeller(name = lab)) +
      labs(title = "약동학 (점선 = 방사선치료 종료일)",
           x = "시간 (days)", y = "농도", colour = NULL) +
      THEME + theme(legend.position = "none")
  })

  output$protNote <- renderUI({
    if (!isTRUE(input$ami))
      return(HTML("<i>아미포스틴이 꺼져 있습니다.</i>"))
    pf <- 0.65 * exp(-log(2) * input$amidel / 8)
    HTML(sprintf(paste0(
      "<div style='background:#FFF6EC;padding:10px;border-left:4px solid ",
      "#E08B2E'><b>시간적 동시성 (temporal coincidence).</b> 조사 ",
      "%d분 전 투여 → 보호계수 <b>%.3f</b> (최대 0.650). ",
      "WR-1065 반감기가 8분이므로 30–60분 지연되면 보호는 거의 사라집니다. ",
      "이것이 아미포스틴 시험 결과가 서로 엇갈리는 기전적 이유에 대한 ",
      "이 모델의 설명입니다.</div>"), input$amidel, pf))
  })

  ## ---- 4  fast loop -------------------------------------------------
  output$fastPlot <- renderPlot({
    s <- sim()
    d <- s %>% select(time, DOM_6, CYT_6, PRM_6, SUR_6, CYTS, CYT_1) %>%
      pivot_longer(-time)
    lab <- c(DOM_6 = "치사손상 저장고 DOOM (bin6)",
             CYT_6 = "국소 사이토카인 (bin6)",
             PRM_6 = "폐포 부종 PERM (bin6)",
             SUR_6 = "계면활성제 SURF (bin6)",
             CYTS  = "전신 사이토카인 (systemic)",
             CYT_1 = "사이토카인 (bin1, 저선량)")
    pk <- s$time[which.max(s$PNIE)]
    ggplot(d, aes(time, value)) +
      geom_line(linewidth = .9, colour = "#D6A81E") +
      geom_vline(xintercept = tcourse(), linetype = "dotted") +
      geom_vline(xintercept = pk, linetype = "dashed",
                 colour = "#C7538C") +
      facet_wrap(~ name, scales = "free_y",
                 labeller = labeller(name = lab)) +
      labs(title = "빠른 고리 — DAMP → NF-κB → 사이토카인 → 부종",
           subtitle = paste0("점선 = RT 종료 (", round(tcourse()),
                             "일), 분홍 파선 = 폐렴 지수 정점 (",
                             round(pk), "일)"),
           x = "시간 (days)", y = NULL) + THEME
  })

  output$gainNote <- renderUI({
    g <- if (isTRUE(input$durv)) 1.35 else 1.00
    budget <- g * (0.025 + 2.0 * 0.015 / 1.0) + 0.10 * 0.05 / 0.30
    HTML(sprintf(paste0(
      "<div style='background:#FFFBEA;padding:10px;border-left:4px solid ",
      "#D6A81E'><b>이득 예산 (gain budget).</b> ",
      "GIMM(KAMP + KDAMP·KDCYT/KMCYT) + KSYSIN·KSOUT/KCLS = ",
      "<b>%.4f</b> &lt; KCE = 0.150. 빠른 고리가 임계 이하이므로 폐렴은 ",
      "스스로 소실됩니다. 더발루맙은 항을 더하지 않고 이득 GIMM 을 ",
      "%.2f 배로 올려 유효 시간상수를 12.8일에서 %.1f일로 늘립니다 — ",
      "그래서 초과 위험이 평균 폐선량에 따라 <i>커집니다</i>.</div>"),
      budget, g, 1 / (0.15 - budget)))
  })

  ## ---- 5  slow loop -------------------------------------------------
  output$slowPlot <- renderPlot({
    d <- sim() %>%
      select(time, TGF_6, MFB_6, COL_6, XCL_6, TGF_4, COL_4) %>%
      pivot_longer(-time)
    lab <- c(TGF_6 = "활성 TGF-β1 (bin6)", MFB_6 = "근섬유아세포 (bin6)",
             COL_6 = "가용성 콜라겐 COL (bin6)",
             XCL_6 = "가교결합 콜라겐 XCOL (bin6)",
             TGF_4 = "활성 TGF-β1 (bin4)",
             COL_4 = "가용성 콜라겐 COL (bin4)")
    ggplot(d, aes(time, value)) +
      geom_line(linewidth = .9, colour = "#9B4FA8") +
      geom_hline(data = tibble(name = c("COL_6", "COL_4"), y = 1.586),
                 aes(yintercept = y), linetype = "dashed",
                 colour = "#C0392B") +
      facet_wrap(~ name, scales = "free_y",
                 labeller = labeller(name = lab)) +
      labs(title = "느린 고리 — TGF-β1 → 근섬유아세포 → 콜라겐 → 강성",
           subtitle = paste0("빨간 파선 = 분리선 COL 1.586. ",
                             "이 선을 넘으면 되돌아오지 않는다."),
           x = "시간 (days)", y = NULL) + THEME
  })

  ## ---- 6  bistable switch -------------------------------------------
  output$switchNote <- renderUI({
    HTML(paste0(
      "<div style='background:#FAF0FB;padding:10px;border-left:4px solid ",
      "#9B4FA8'><b>느린 고리는 이중안정입니다.</b> 세 고정점은 ",
      "파라미터에서 대수적으로 계산되며 시뮬레이션이 이를 재현합니다. ",
      "섬유화는 <i>연속적인 용량-반응</i>이 아니라 <b>문턱</b>입니다 — ",
      "40 Gy 에서 콜라겐 초과가 0.09, 42 Gy 에서 2.85 입니다. ",
      "임상에서 방사선 섬유화가 등선량선으로 <b>날카롭게 경계지어지는</b> ",
      "이유가 여기 있습니다.</div>"))
  })

  output$phasePlot <- renderPlot({
    s <- sim()
    d <- bind_rows(
      tibble(bin = "bin6 (최고선량)", TGF = s$TGF_6,
             COLT = s$COL_6 + s$XCL_6, time = s$time),
      tibble(bin = "bin5", TGF = s$TGF_5, COLT = s$COL_5 + s$XCL_5,
             time = s$time),
      tibble(bin = "bin4", TGF = s$TGF_4, COLT = s$COL_4 + s$XCL_4,
             time = s$time),
      tibble(bin = "bin3", TGF = s$TGF_3, COLT = s$COL_3 + s$XCL_3,
             time = s$time))
    fp <- tibble(TGF = c(0.0401, 0.2115, 0.7767),
                 COLT = c(1.000, 1.586, 2.940),
                 lab = c("정상 (stable)", "분리선 (UNSTABLE)",
                         "섬유화 (stable)"))
    ggplot(d, aes(TGF, COLT, colour = bin)) +
      geom_path(linewidth = .9) +
      geom_point(data = fp, aes(TGF, COLT), inherit.aes = FALSE,
                 size = 4, shape = c(19, 1, 19), colour = "#7D3C98") +
      geom_text(data = fp, aes(TGF, COLT, label = lab), inherit.aes = FALSE,
                hjust = -.12, size = 3.3, colour = "#7D3C98") +
      geom_hline(yintercept = 1.586, linetype = "dashed",
                 colour = "#C0392B") +
      labs(title = "위상 평면 — 활성 TGF-β1 대 총 콜라겐",
           subtitle = "각 궤적은 하나의 선량 구간. 분리선을 넘은 구간만 섬유화한다.",
           x = "활성 TGF-β1", y = "총 콜라겐 (가용성 + 가교결합)",
           colour = NULL) + THEME
  })

  output$fpTab <- renderTable({
    tibble(
      `고정점 (fixed point)` = c("정상 (healthy)", "분리선 (separatrix)",
                                 "섬유화 (fibrotic)"),
      `활성 TGF-β1` = c(0.0401, 0.2115, 0.7767),
      `근섬유아세포` = c(0.00057, 0.3563, 0.9692),
      `총 콜라겐` = c(1.000, 1.586, 2.940),
      `안정성` = c("안정 (stable)", "불안정 (UNSTABLE)", "안정 (stable)"))
  }, digits = 4, striped = TRUE)

  ## ---- 7  pulmonary function ---------------------------------------
  output$pftPlot <- renderPlot({
    d <- sim() %>% select(time, DLCO, FVC) %>% pivot_longer(-time)
    ggplot(d, aes(time, value, colour = name)) +
      geom_line(linewidth = 1) +
      geom_vline(xintercept = tcourse(), linetype = "dotted") +
      scale_colour_manual(values = c(DLCO = "#4CA85F", FVC = "#3FA3AD")) +
      labs(title = "폐기능 경과 (DLCO · FVC)",
           subtitle = paste0("영구 성분은 모세혈관 소실 천장 ",
                             "EC_max = 1/(1+KECIRR·BED) 에서 온다 — ",
                             "섬유화하지 않은 구간에도 존재한다."),
           x = "시간 (days)", y = "% predicted", colour = NULL) + THEME
  })

  output$binMap <- renderPlot({
    s <- tail(sim(), 1)
    v <- vvec()
    d <- tibble(
      bin = factor(1:6),
      dose = input$DRX * REPF,
      vol = 100 * v,
      AT2 = c(s$AT2_1, s$AT2_2, s$AT2_3, s$AT2_4, s$AT2_5, s$AT2_6),
      EC = c(s$EC_1, s$EC_2, s$EC_3, s$EC_4, s$EC_5, s$EC_6),
      COLT = c(s$COL_1 + s$XCL_1, s$COL_2 + s$XCL_2, s$COL_3 + s$XCL_3,
               s$COL_4 + s$XCL_4, s$COL_5 + s$XCL_5,
               s$COL_6 + s$XCL_6)) %>%
      pivot_longer(c(AT2, EC, COLT))
    lab <- c(AT2 = "생존 AT2", EC = "내피 온전성 EC",
             COLT = "총 콜라겐")
    ggplot(d, aes(dose, value, size = vol)) +
      geom_point(colour = "#6A75C9", alpha = .8) +
      facet_wrap(~ name, scales = "free_y",
                 labeller = labeller(name = lab)) +
      scale_size_continuous(range = c(2, 9), name = "폐 체적 (%)") +
      labs(title = paste0("구간별 최종 상태 (day ", input$tend, ")"),
           subtitle = "점 크기 = 그 구간이 차지하는 폐 체적",
           x = "구간 대표선량 (Gy)", y = NULL) + THEME
  })

  ## ---- 8  clinical endpoints ---------------------------------------
  output$ntcpPlot <- renderPlot({
    s <- sim()
    ggplot(s, aes(time, 100 * NTCP)) +
      geom_line(linewidth = 1, colour = "#C7538C") +
      geom_vline(xintercept = tcourse(), linetype = "dotted") +
      geom_hline(yintercept = 20, linetype = "dashed",
                 colour = "#888888") +
      annotate("text", x = input$tend * .7, y = 21,
               label = "QUANTEC 20% (MLD 20 Gy)", size = 3,
               colour = "#888888") +
      labs(title = "NTCP — CTCAE grade ≥ 2 방사선 폐렴 위험",
           x = "시간 (days)", y = "위험 (%)") + THEME
  })

  output$gradePlot <- renderPlot({
    s <- sim()
    ggplot(s, aes(time, GRADE)) +
      geom_step(linewidth = 1, colour = "#E06AA0") +
      geom_vline(xintercept = tcourse(), linetype = "dotted") +
      scale_y_continuous(breaks = 0:4, limits = c(0, 4)) +
      labs(title = "CTCAE 등급 (연속 상태에서 산출)",
           subtitle = "등급은 입력이 아니라 출력이다",
           x = "시간 (days)", y = "grade") + THEME
  })

  output$endTab <- renderTable({
    s <- sim(); tc <- tcourse()
    pk <- s[which.max(s$PNIE), ]
    d365 <- s[which.min(abs(s$time - 365)), ]
    tibble(
      `엔드포인트 (endpoint)` = c(
        "폐렴 지수 정점 PNIe", "정점 시점 (RT 종료 후, 일)",
        "정점 시점 (주)", "NTCP grade≥2 (%)", "최고 CTCAE 등급",
        "DLCO 최저 (%pred)", "DLCO 365일 변화 (%)",
        "섬유화 전환 체적 (%)", "섬유화 점수 FIB 365일"),
      `값 (value)` = c(
        sprintf("%.4f", max(s$PNIE)),
        sprintf("%.0f", pk$time - tc), sprintf("%.1f", (pk$time - tc) / 7),
        sprintf("%.1f", 100 * max(s$NTCP)),
        sprintf("%.0f", max(s$GRADE)),
        sprintf("%.2f", min(s$DLCO)),
        sprintf("%.2f", 100 * (d365$DLCO - input$DLCO0) / input$DLCO0),
        sprintf("%.1f", 100 * last(s$VFIB)),
        sprintf("%.4f", d365$FIB)))
  }, striped = TRUE)

  ## ---- 9  therapeutic ratio ----------------------------------------
  trData <- reactive({
    bind_rows(lapply(names(PLANS), function(nm) {
      s <- simplan(nm)
      tibble(plan = nm, MLD = first(s$MLD),
             dfx = PLANS[[nm]]$DRX / PLANS[[nm]]$NFX,
             TCP = max(s$TCP), NTCP = 100 * max(s$NTCP),
             UCP = max(s$TCP) / 100 * (1 - max(s$NTCP)) * 100)
    }))
  })

  output$trPlot <- renderPlot({
    d <- trData()
    ggplot(d, aes(NTCP, TCP, label = plan)) +
      geom_point(aes(size = dfx, colour = MLD)) +
      geom_text(hjust = -.08, size = 3.2) +
      scale_colour_gradient(low = "#4CA85F", high = "#C0392B",
                            name = "MLD (Gy)") +
      scale_size_continuous(range = c(3, 10), name = "분할선량 (Gy)") +
      expand_limits(x = c(0, 45)) +
      labs(title = "치료비 — 종양 제어확률 대 폐렴 위험",
           subtitle = paste0("좌상단이 좋다. SBRT 가 두 축 모두에서 이기는 ",
                             "것은 방사선생물학이 아니라 기하학 때문이다 ",
                             "(MLD 17.8 → 4.3 Gy)."),
           x = "NTCP grade≥2 (%)", y = "TCP (%)") + THEME
  })

  output$trTab <- renderTable({
    trData() %>% arrange(desc(UCP)) %>%
      rename(`플랜` = plan, `분할선량 (Gy)` = dfx,
             `MLD (Gy)` = MLD, `TCP (%)` = TCP,
             `NTCP (%)` = NTCP, `UCP (%)` = UCP)
  }, digits = 2, striped = TRUE)

  ## ---- 10  scenario comparison -------------------------------------
  scenData <- reactive({
    req(length(input$scen) > 0)
    bind_rows(lapply(input$scen, simplan))
  })

  output$scenPlot <- renderPlot({
    d <- scenData() %>%
      select(time, plan, DLCO, PNIE, NTCP, VFIB) %>%
      mutate(NTCP = 100 * NTCP, VFIB = 100 * VFIB) %>%
      pivot_longer(c(DLCO, PNIE, NTCP, VFIB))
    lab <- c(DLCO = "DLCO (%pred)", PNIE = "폐렴 지수 PNIe",
             NTCP = "NTCP grade≥2 (%)", VFIB = "섬유화 전환 체적 (%)")
    ggplot(d, aes(time, value, colour = plan)) +
      geom_line(linewidth = .9) +
      facet_wrap(~ name, scales = "free_y",
                 labeller = labeller(name = lab)) +
      labs(title = "시나리오 비교", x = "시간 (days)", y = NULL,
           colour = NULL) + THEME
  })

  output$scenTab <- renderTable({
    scenData() %>% group_by(plan) %>%
      summarise(`MLD (Gy)` = first(MLD), `V20 (%)` = 100 * first(V20G),
                `V40 (%)` = 100 * first(V40G),
                `PNIe peak` = max(PNIE),
                `NTCP (%)` = 100 * max(NTCP),
                `DLCO nadir` = min(DLCO),
                `Vfib (%)` = 100 * last(VFIB), .groups = "drop") %>%
      arrange(`MLD (Gy)`)
  }, digits = 2, striped = TRUE)
}

shinyApp(ui, server)
