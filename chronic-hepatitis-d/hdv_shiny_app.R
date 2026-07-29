## =============================================================================
##  CHRONIC HEPATITIS D (HDV) -- QSP DASHBOARD  (Shiny)
##  hdv_shiny_app.R      companion to hdv_mrgsolve_model.R
## =============================================================================
##
##  RUN
##      install.packages(c("shiny", "mrgsolve", "ggplot2", "dplyr", "tidyr", "DT"))
##      shiny::runApp("hdv_shiny_app.R")
##
##  WHAT THIS APP IS FOR
##  --------------------
##  Not "here are some curves".  Each tab answers one question that the model
##  exists to answer, and the tab titles say which:
##
##   1  Patient       who is this patient, and where does the floor sit for them
##   2  Drug PK       exposure and, for bulevirtide, DERIVED NTCP occupancy
##   3  Virology      serum HDV RNA + the intracellular pool, side by side --
##                    this is where the lonafarnib paradox becomes visible
##   4  Endpoints     HDV RNA response, ALT normalisation, combined response
##   5  Decoupling    why ALT and HDV RNA disagree, plotted against each other
##   6  The FLOOR     decomposition of the inflow maintaining Id, and the
##                    ceiling of the entire entry-inhibitor class
##   7  Bile acids    total bile acids as an NTCP-occupancy read-out, and the
##                    benefit-saturates-while-cost-does-not dose curve
##   8  Compare       any number of regimens overlaid on any output
##   9  Outcomes      fibrosis, platelets, HCC over 5 years
##
##  The numbers the app should reproduce are in hdv_model_report.txt.
## =============================================================================

library(shiny)
library(mrgsolve)
library(ggplot2)
library(dplyr)
library(tidyr)

MODEL_FILE <- "hdv_mrgsolve_model.R"
mod0 <- mread_cache("hdv", MODEL_FILE)

MW_BLV <- 5398.9

## ---------------------------------------------------------------- dosing ----
ev_blv <- function(mg, dur) if (mg <= 0 || dur <= 0) NULL else
  ev(time = 0, amt = mg * 1e-3 / MW_BLV * 1e9, cmt = "ASC", ii = 1, addl = dur - 1)
ev_ifn <- function(ug, dur) if (ug <= 0 || dur < 7) NULL else
  ev(time = 0, amt = ug, cmt = "ISC", ii = 7, addl = floor(dur / 7) - 1)
ev_lnf <- function(mg, rtv, dur) {
  if (mg <= 0 || dur <= 0) return(NULL)
  e <- ev(time = 0, amt = mg * 1e-3 / 638.6 * 1e6, cmt = "LGUT",
          ii = 0.5, addl = dur * 2 - 1)
  if (rtv > 0) e <- e + ev(time = 0, amt = 0.7 * rtv * 1e-3 / 720.9 * 1e6,
                           cmt = "RCEN", ii = 0.5, addl = dur * 2 - 1)
  e
}
ev_sir <- function(mg, dur) if (mg <= 0 || dur < 28) NULL else
  ev(time = 0, amt = mg * 1e-3 / 16000 * 1e9, cmt = "QCEN",
     ii = 28, addl = floor(dur / 28) - 1)

bind_ev <- function(...) {
  es <- Filter(Negate(is.null), list(...))
  if (!length(es)) return(NULL)
  Reduce(`+`, es)
}

## ------------------------------------------------------------ simulation ----
simulate_regimen <- function(pars, horizon_wk, blv, ifn, lam, lnf, rtv, sir,
                             tx_wk, occfix = NA) {
  dur <- tx_wk * 7
  m <- param(mod0, pars)
  if (!is.na(occfix)) m <- param(m, list(OCCFIX = occfix))
  m <- param(m, list(LAMBDA = if (lam) 1 else 0))
  e <- bind_ev(ev_blv(blv, dur), ev_ifn(ifn, dur),
               ev_lnf(lnf, rtv, dur), ev_sir(sir, dur))
  m <- update(m, end = horizon_wk * 7, delta = 1)
  d <- as.data.frame(if (is.null(e)) mrgsim(m) else mrgsim(m, events = e))
  d$week <- d$time / 7
  d
}

## per-patient parameter overrides from the sidebar
patient_pars <- function(input) {
  list(
    DIMM  = input$dimm,
    DHEP  = input$dhep,
    KCC   = 0.0137499 * input$ccmult,
    PHINT = input$phint,
    KCURE = input$kcure,
    ALTBASE = input$altbase,
    NUC   = if (input$nuc) 1 else 0
  )
}

## initial conditions from the sidebar (baseline severity)
patient_init <- function(m, input) {
  init(m, list(VD = 10^input$lgvd0, SSER = 10^input$lgsag0,
               ID = input$id0, ALT = input$alt0, FIB = input$fib0))
}

simulate_full <- function(input, blv, ifn, lam, lnf, rtv, sir, tx_wk,
                          horizon_wk, occfix = NA) {
  m <- param(mod0, patient_pars(input))
  if (!is.na(occfix)) m <- param(m, list(OCCFIX = occfix))
  m <- param(m, list(LAMBDA = if (lam) 1 else 0))
  m <- patient_init(m, input)
  dur <- tx_wk * 7
  e <- bind_ev(ev_blv(blv, dur), ev_ifn(ifn, dur),
               ev_lnf(lnf, rtv, dur), ev_sir(sir, dur))
  m <- update(m, end = horizon_wk * 7, delta = 1)
  d <- as.data.frame(if (is.null(e)) mrgsim(m) else mrgsim(m, events = e))
  d$week <- d$time / 7
  d$dlog <- d$LGVD - input$lgvd0
  d
}

THEME <- theme_minimal(base_size = 12) +
  theme(panel.grid.minor = element_blank(),
        plot.title = element_text(face = "bold", size = 12),
        plot.subtitle = element_text(size = 10, colour = "grey30"),
        legend.position = "bottom")

## ==================================================================== UI ====
ui <- fluidPage(
  titlePanel("만성 D형 간염 (Chronic Hepatitis D) — QSP 대시보드"),
  tags$p(style = "color:#555;margin-top:-8px",
         HTML("HDV를 <b>두 기질 조립 라인</b>으로 놓은 모델: 유전체는 숙주 RNA Pol II가, ",
              "외피는 HBV cccDNA가 공급한다. 모든 약은 <b>어떤 플럭스를 끊는지</b>로 분류된다 — ",
              "진입(불레비티드) · 외피화(로나파닙·siRNA) · 세포내 복제와 감염세포 소실(인터페론).")),

  sidebarLayout(
    sidebarPanel(
      width = 3,
      h4("① 환자 프로파일"),
      sliderInput("lgvd0", "기저 HDV RNA (log10 IU/mL)", 3.0, 8.0, 5.5, 0.1),
      sliderInput("lgsag0", "기저 HBsAg (log10 IU/mL)", 2.0, 5.0, 3.85, 0.05),
      sliderInput("id0", "HDAg 양성 간세포 분율", 0.01, 0.30, 0.10, 0.01),
      sliderInput("alt0", "기저 ALT (U/L)", 40, 400, 110, 5),
      sliderInput("altbase", "비-HDV ALT 성분 (U/L)", 15, 80, 42, 1),
      sliderInput("fib0", "기저 섬유화 (Ishak 0–6)", 0, 6, 2, 0.5),
      checkboxInput("nuc", "핵산유사체(NUC) 병용", TRUE),

      hr(), h4("② 숙주 세포생물학 — 바닥(floor)"),
      helpText(HTML("<small>이 네 개가 진입억제제의 한계를 결정한다. ",
                    "표적 친화도가 아니다.</small>")),
      sliderInput("dimm", "감염세포 살해율 DIMM (1/day)", 0.005, 0.08, 0.02636, 0.001),
      sliderInput("dhep", "간세포 전환율 DHEP (1/day)", 0.001, 0.012, 0.004, 0.0005),
      sliderInput("ccmult", "세포간 전파 배수 (× KCC)", 0.0, 3.0, 1.0, 0.1),
      sliderInput("phint", "세포간 전파의 NTCP 의존 분율", 0.0, 1.0, 0.5, 0.05),
      sliderInput("kcure", "비세포용해성 세포내 제거 (1/day)", 0.0, 0.01, 0.002, 0.0005),

      hr(), h4("③ 치료"),
      sliderInput("blv", "불레비티드 (mg SC qd)", 0, 20, 2, 0.5),
      sliderInput("ifn", "Peg-IFN (µg SC qw)", 0, 240, 0, 10),
      checkboxInput("lam", "Peg-IFN 람다로 전환 (혈액학적 독성 없음, 간독성 있음)", FALSE),
      sliderInput("lnf", "로나파닙 (mg PO BID)", 0, 100, 0, 5),
      sliderInput("rtv", "리토나비르 부스터 (mg PO BID)", 0, 200, 100, 50),
      sliderInput("sir", "HBsAg siRNA (mg SC q4w)", 0, 400, 0, 50),
      sliderInput("txwk", "치료 기간 (주)", 4, 144, 48, 4),
      sliderInput("horiz", "관찰 기간 (주)", 12, 260, 96, 4),
      checkboxInput("perfect", "가상의 '완전' 진입억제제 (점유율 100%) 함께 표시", TRUE)
    ),

    mainPanel(
      width = 9,
      tabsetPanel(
        type = "tabs",

        tabPanel("1 · 환자·바닥", br(),
          fluidRow(column(6, h5("이 환자의 Id 유지 플럭스 분해"),
                          plotOutput("floorbar", height = "300px")),
                   column(6, h5("간세포 풀"), plotOutput("pools", height = "300px"))),
          hr(), verbatimTextOutput("patient_txt")),

        tabPanel("2 · 약물 PK와 표적 점유율", br(),
          fluidRow(column(6, plotOutput("pk_blv", height = "280px")),
                   column(6, plotOutput("pk_occ", height = "280px"))),
          fluidRow(column(6, plotOutput("pk_ifn", height = "260px")),
                   column(6, plotOutput("pk_lnf", height = "260px"))),
          hr(), verbatimTextOutput("occ_txt")),

        tabPanel("3 · 바이러스학 (혈중 vs 세포내)", br(),
          plotOutput("vir_serum", height = "300px"),
          plotOutput("vir_intra", height = "280px"),
          helpText(HTML("<b>로나파닙 역설.</b> 조립을 막으면 유전체가 세포 안에 갇힌다. ",
                        "혈중 HDV RNA는 <i>조립 플럭스</i>를 재는 것이고 저수지 크기를 재는 것이 아니다 — ",
                        "그래서 중단 시 반동이 빠르다."))),

        tabPanel("4 · 임상 엔드포인트", br(),
          plotOutput("ep_plot", height = "320px"),
          hr(), h5("주요 시점 요약"), tableOutput("ep_tab"),
          helpText(HTML("반응 정의는 MYR301과 동일: HDV RNA 검출한계 미만 또는 2 log10 이상 감소, ",
                        "그리고 ALT 정상화. 복합반응은 두 조건 모두."))),

        tabPanel("5 · ALT–RNA 해리", br(),
          fluidRow(column(7, plotOutput("dec_time", height = "320px")),
                   column(5, plotOutput("dec_scatter", height = "320px"))),
          hr(), verbatimTextOutput("dec_txt")),

        tabPanel("6 · 바닥과 천장", br(),
          plotOutput("ceiling", height = "320px"),
          hr(), h5("바닥 파라미터 민감도 (48주 log10 감소량)"),
          tableOutput("sens_tab"),
          helpText(HTML("<b>이 모델의 반증 가능한 핵심 주장:</b> 48주 반응은 ",
                        "표적 친화도(KDNTCP)보다 숙주 세포생물학(DHEP, KCC, PHINT)에 ",
                        "더 민감하다."))),

        tabPanel("7 · 담즙산 = 점유율 계측기", br(),
          fluidRow(column(6, plotOutput("ba_time", height = "300px")),
                   column(6, plotOutput("ba_dose", height = "300px"))),
          hr(), tableOutput("ba_tab"),
          helpText(HTML("NTCP는 HDV 수용체이면서 담즙산 수송체다. 그래서 총담즙산 상승은 ",
                        "<b>표적 점유율의 직접 판독값</b>이 된다. 효능은 포화하지만 ",
                        "비용은 포화하지 않는다."))),

        tabPanel("8 · 시나리오 비교", br(),
          checkboxGroupInput("cmp", "비교할 요법",
            choices = c("무치료 (NUC only)" = "none",
                        "BLV 2 mg" = "blv2",
                        "BLV 10 mg" = "blv10",
                        "Peg-IFN alfa 180 µg" = "ifn",
                        "Peg-IFN lambda 180 µg" = "lam",
                        "BLV 2 mg + Peg-IFN" = "blvifn",
                        "LNF 50 + RTV 100 BID" = "lnf",
                        "LNF/RTV + Peg-IFN" = "lnfifn",
                        "BLV 2 mg + HBsAg siRNA" = "blvsir",
                        "완전 진입 차단 (가상)" = "perfect"),
            selected = c("none", "blv2", "blv10", "ifn", "blvifn", "lnf"),
            inline = TRUE),
          selectInput("cmp_y", "표시할 변수",
            choices = c("Δlog10 HDV RNA" = "dlog", "ALT (U/L)" = "ALT",
                        "HDAg 양성 간세포 (%)" = "HDAGPC",
                        "세포내 HDV RNA (배수)" = "RGFOLD",
                        "혈중 HBsAg (log10)" = "LGSAG",
                        "총담즙산 (배수)" = "TBAFLD",
                        "탈진 수준" = "EXH",
                        "섬유화 (Ishak)" = "FIB",
                        "혈소판 (10^9/L)" = "PLT",
                        "호중구 (10^9/L)" = "NEU"),
            selected = "dlog"),
          plotOutput("cmp_plot", height = "420px"),
          hr(), tableOutput("cmp_tab")),

        tabPanel("9 · 장기 결과 (5년)", br(),
          fluidRow(column(6, plotOutput("out_fib", height = "290px")),
                   column(6, plotOutput("out_hcc", height = "290px"))),
          fluidRow(column(6, plotOutput("out_plt", height = "260px")),
                   column(6, plotOutput("out_alt", height = "260px"))),
          helpText(HTML("섬유화는 HDV RNA가 아니라 <b>ALT</b>를 따라간다. ",
                        "그래서 요법의 장기 결과 순위는 바이러스학이 아니라 ",
                        "염증 정상화로 결정된다."))),

        tabPanel("ℹ 모델 정보", br(),
          verbatimTextOutput("about"))
      )
    )
  )
)

## ================================================================ SERVER ====
server <- function(input, output, session) {

  cur <- reactive({
    simulate_full(input, input$blv, input$ifn, input$lam, input$lnf,
                  input$rtv, input$sir, input$txwk, input$horiz)
  })
  perf <- reactive({
    if (!isTRUE(input$perfect)) return(NULL)
    simulate_full(input, 0, 0, FALSE, 0, 0, 0, input$txwk, input$horiz,
                  occfix = 1.0)
  })
  untx <- reactive({
    simulate_full(input, 0, 0, FALSE, 0, 0, 0, input$txwk, input$horiz)
  })

  ## ---- 1. patient & floor ----
  output$floorbar <- renderPlot({
    d <- cur()
    p <- unlist(param(mod0))
    ID <- input$id0; IB <- 0.35
    occ <- d$OCCUP[which.min(abs(d$week - min(input$txwk, input$horiz)))]
    fentry <- p[["BETAD"]] * 10^input$lgvd0 * IB
    fcc    <- 0.0137499 * input$ccmult * ID * IB
    dthId  <- input$dhep + p[["DINN"]] + input$dimm
    Hd     <- (input$dhep + p[["DHBV"]]) * IB + dthId * ID
    fdiv   <- p[["LOCREN"]] * Hd * ID / (IB + ID)
    df <- data.frame(
      flux = factor(c("(E) 진입/재감염", "(C) 세포간 전파", "(D) 분열매개 전파"),
                    levels = c("(E) 진입/재감염", "(C) 세포간 전파", "(D) 분열매개 전파")),
      untreated = c(fentry, fcc, fdiv),
      treated = c(fentry * (1 - occ), fcc * (1 - input$phint * occ), fdiv)) %>%
      pivot_longer(-flux, names_to = "state", values_to = "v")
    df$state <- factor(df$state, levels = c("untreated", "treated"),
                       labels = c("무치료", sprintf("치료 중 (점유율 %.2f)", occ)))
    ggplot(df, aes(flux, v, fill = state)) +
      geom_col(position = "dodge") +
      scale_fill_manual(values = c("#b04a4a", "#2f7a2f"), name = NULL) +
      labs(x = NULL, y = "유입 플럭스 (분율/일)",
           title = "감염세포 풀을 유지하는 세 플럭스",
           subtitle = "(D)는 NTCP 비의존 — 진입억제제가 절대 닿지 못하는 바닥") +
      THEME
  })

  output$pools <- renderPlot({
    d <- cur() %>% select(week, T, IB, ID) %>%
      pivot_longer(-week, names_to = "pool", values_to = "v")
    d$pool <- factor(d$pool, levels = c("T", "IB", "ID"),
                     labels = c("T · HBsAg 음성", "Ib · HBsAg+ HDV−(외피 공급)",
                                "Id · HBsAg+ HDV+"))
    ggplot(d, aes(week, v, colour = pool)) + geom_line(linewidth = .8) +
      scale_colour_manual(values = c("#556699", "#c07a30", "#b04a4a"), name = NULL) +
      labs(x = "주", y = "간세포 분율", title = "간세포 구획") + THEME
  })

  output$patient_txt <- renderPrint({
    d <- cur(); u <- untx()
    i48 <- which.min(abs(d$week - 48))
    cat(sprintf(
"기저: HDV RNA %.2f log10 IU/mL · HBsAg %.2f log10 · HDAg+ 간세포 %.0f%% · ALT %.0f U/L · Ishak %.1f
48주: HDV RNA 변화 %+.2f log10 · ALT %.0f U/L (%s) · HDAg+ %.2f%% · 세포내 RNA %.2f배
     NTCP 점유율 %.3f · 자유 NTCP %.1f%% · 총담즙산 %.1f배
무치료 대조 48주: HDV RNA 변화 %+.2f log10 · ALT %.0f U/L
",
      input$lgvd0, input$lgsag0, 100 * input$id0, input$alt0, input$fib0,
      d$dlog[i48], d$ALT[i48],
      ifelse(d$ALT[i48] <= 40, "정상화", "비정상"),
      d$HDAGPC[i48], d$RGFOLD[i48], d$OCCUP[i48], 100 * d$FREENT[i48],
      d$TBAFLD[i48], u$dlog[i48], u$ALT[i48]))
  })

  ## ---- 2. PK ----
  output$pk_blv <- renderPlot({
    d <- cur() %>% filter(week <= min(8, input$horiz))
    ggplot(d, aes(week * 7, CBLV)) + geom_line(colour = "#2f7a2f", linewidth = .8) +
      labs(x = "일", y = "불레비티드 (nM)", title = "불레비티드 혈중 농도",
           subtitle = "포화성 표적매개 소실 → 용량 초비례적 노출") + THEME
  })
  output$pk_occ <- renderPlot({
    d <- cur()
    ggplot(d, aes(week, OCCUP)) + geom_line(colour = "#2f7a2f", linewidth = .9) +
      geom_hline(yintercept = 1, linetype = 3) +
      scale_y_continuous(limits = c(0, 1)) +
      labs(x = "주", y = "NTCP 점유율", title = "NTCP 점유율 (도출값, 가정값 아님)",
           subtitle = "담즙산 상승 배수 2개로부터 역산된 Kd를 통해 계산") + THEME
  })
  output$pk_ifn <- renderPlot({
    d <- cur()
    ggplot(d, aes(week)) +
      geom_line(aes(y = CIFN, colour = "Peg-IFN (ng/mL)"), linewidth = .8) +
      geom_line(aes(y = 20 * ISG, colour = "ISG × 20"), linewidth = .8) +
      geom_line(aes(y = 20 * SOCS, colour = "SOCS1 × 20"), linewidth = .6, linetype = 2) +
      scale_colour_manual(values = c("#2f6a9a", "#3a8a8a", "#7d4a9a"), name = NULL) +
      labs(x = "주", y = NULL, title = "인터페론 노출과 ISG 유도",
           subtitle = "SOCS1 음성 피드백 = 내성(tachyphylaxis)") + THEME
  })
  output$pk_lnf <- renderPlot({
    d <- cur() %>% filter(week <= min(4, input$horiz))
    ggplot(d, aes(week * 7)) +
      geom_line(aes(y = CLNF, colour = "로나파닙 (µM)"), linewidth = .8) +
      geom_line(aes(y = CRTV, colour = "리토나비르 (µM)"), linewidth = .8) +
      geom_line(aes(y = IFTC, colour = "FTase 억제율"), linewidth = .9, linetype = 2) +
      scale_colour_manual(values = c("#7d4a9a", "#a04a70", "#333333"), name = NULL) +
      labs(x = "일", y = NULL, title = "로나파닙 / 리토나비르",
           subtitle = "CYP3A4 억제 부스팅으로 노출 4–5배") + THEME
  })
  output$occ_txt <- renderPrint({
    p <- unlist(param(mod0))
    css <- function(mg) {
      inp <- p[["FBLV"]] * (mg * 1e-3 / MW_BLV * 1e9)
      f <- function(c) p[["CLBLV"]] * c + p[["VMBLV"]] * c / (p[["KMBLV"]] + c) - inp
      uniroot(f, c(0, 1e5))$root
    }
    o <- function(mg) { c <- css(mg); c / (c + p[["KDNTCP"]]) }
    cat(sprintf(
"적합된 Kd_NTCP = %.3f nM  (담즙산 상승 배수 2개만으로 역산)
  2 mg : Css %.2f nM → 점유율 %.4f → 자유 NTCP %.1f%%  → 담즙산 %.1f배
 10 mg : Css %.2f nM → 점유율 %.4f → 자유 NTCP %.1f%%  → 담즙산 %.1f배
잔여 진입 플럭스 비 (2 mg / 10 mg) = %.1f배

핵심: 승인 용량 2 mg는 NTCP를 포화시키지 않는다. 그런데도 10 mg가 복합반응을
거의 올리지 못한다 — 진입이 Id 유지 유입의 일부일 뿐이기 때문이다.
",
      p[["KDNTCP"]], css(2), o(2), 100 * (1 - o(2)),
      (1 + p[["ROATP"]]) / ((1 - o(2)) + p[["ROATP"]]),
      css(10), o(10), 100 * (1 - o(10)),
      (1 + p[["ROATP"]]) / ((1 - o(10)) + p[["ROATP"]]),
      (1 - o(2)) / (1 - o(10))))
  })

  ## ---- 3. virology ----
  output$vir_serum <- renderPlot({
    d <- cur(); u <- untx(); pf <- perf()
    df <- rbind(data.frame(week = d$week, v = d$dlog, arm = "선택한 요법"),
                data.frame(week = u$week, v = u$dlog, arm = "무치료"))
    if (!is.null(pf)) df <- rbind(df, data.frame(week = pf$week, v = pf$dlog,
                                                 arm = "완전 진입 차단 (가상)"))
    ggplot(df, aes(week, v, colour = arm)) + geom_line(linewidth = .9) +
      geom_hline(yintercept = -2, linetype = 2, colour = "grey40") +
      annotate("text", x = 2, y = -2.15, label = "2 log10 반응 기준",
               hjust = 0, size = 3, colour = "grey35") +
      geom_vline(xintercept = input$txwk, linetype = 3) +
      scale_colour_manual(values = c("#b04a4a", "#2f7a2f", "#888888"), name = NULL) +
      labs(x = "주", y = expression(Delta*log[10]*" HDV RNA"),
           title = "혈중 HDV RNA",
           subtitle = "점선 수직선 = 치료 종료. 진입억제제는 1상(첫 주 급감)이 없다") + THEME
  })
  output$vir_intra <- renderPlot({
    d <- cur()
    ggplot(d, aes(week)) +
      geom_line(aes(y = RGFOLD, colour = "세포내 HDV RNA (배수)"), linewidth = .9) +
      geom_line(aes(y = 10^(d$dlog), colour = "혈중 HDV RNA (배수)"), linewidth = .9) +
      geom_hline(yintercept = 1, linetype = 3) +
      scale_y_log10() +
      scale_colour_manual(values = c("#c07a30", "#b04a4a"), name = NULL) +
      labs(x = "주", y = "기저 대비 배수 (log 축)",
           title = "세포내 vs 혈중 — 같은 방향으로 가지 않는다",
           subtitle = "로나파닙을 켜면 두 곡선이 갈라진다 (조립 차단으로 유전체가 갇힌다)") + THEME
  })

  ## ---- 4. endpoints ----
  output$ep_plot <- renderPlot({
    d <- cur()
    df <- data.frame(week = rep(d$week, 3),
                     v = c(d$VR, d$ALTN, d$COMB),
                     ep = rep(c("HDV RNA 반응", "ALT 정상화", "복합반응"),
                              each = nrow(d)))
    ggplot(df, aes(week, v, colour = ep)) +
      geom_step(linewidth = .9) +
      geom_vline(xintercept = input$txwk, linetype = 3) +
      scale_y_continuous(breaks = c(0, 1), labels = c("미달", "달성")) +
      scale_colour_manual(values = c("#2f6a9a", "#b04a4a", "#111111"), name = NULL) +
      labs(x = "주", y = NULL, title = "엔드포인트 달성 여부 (이 환자)",
           subtitle = "집단 반응률은 hdv_reference_model.py의 가상 집단에서 계산된다") + THEME
  })
  output$ep_tab <- renderTable({
    d <- cur()
    wks <- c(4, 12, 24, 48, 96, input$horiz)
    wks <- unique(wks[wks <= input$horiz])
    do.call(rbind, lapply(wks, function(wk) {
      i <- which.min(abs(d$week - wk))
      data.frame(주 = wk,
                 `Δlog10 HDV` = sprintf("%+.2f", d$dlog[i]),
                 `ALT (U/L)` = sprintf("%.0f", d$ALT[i]),
                 `HDAg+ (%)` = sprintf("%.2f", d$HDAGPC[i]),
                 `세포내 RNA (배)` = sprintf("%.2f", d$RGFOLD[i]),
                 `HBsAg (log10)` = sprintf("%.2f", d$LGSAG[i]),
                 `담즙산 (배)` = sprintf("%.1f", d$TBAFLD[i]),
                 `RNA 반응` = ifelse(d$VR[i] > .5, "O", "-"),
                 `ALT 정상` = ifelse(d$ALTN[i] > .5, "O", "-"),
                 `복합` = ifelse(d$COMB[i] > .5, "O", "-"),
                 check.names = FALSE)
    }))
  })

  ## ---- 5. decoupling ----
  output$dec_time <- renderPlot({
    d <- cur()
    ggplot(d, aes(week)) +
      geom_line(aes(y = 100 * ALT / input$alt0, colour = "ALT (기저 %)"),
                linewidth = .9) +
      geom_line(aes(y = 100 * 10^dlog, colour = "HDV RNA (기저 %)"),
                linewidth = .9) +
      scale_colour_manual(values = c("#b04a4a", "#2f6a9a"), name = NULL) +
      labs(x = "주", y = "기저 대비 %",
           title = "ALT가 먼저 떨어진다",
           subtitle = "ALT는 새로 감염되는 세포의 '유입'을 읽고, HDV RNA는 감염 풀의 '크기'를 읽는다") +
      THEME
  })
  output$dec_scatter <- renderPlot({
    d <- cur()
    ggplot(d, aes(dlog, ALT, colour = week)) + geom_path(linewidth = .9) +
      geom_hline(yintercept = 40, linetype = 2) +
      geom_vline(xintercept = -2, linetype = 2) +
      annotate("text", x = -0.2, y = 38, label = "ALT 정상화만 달성",
               hjust = 1, size = 3, colour = "grey30") +
      scale_colour_viridis_c(name = "주") +
      labs(x = expression(Delta*log[10]*" HDV RNA"), y = "ALT (U/L)",
           title = "두 엔드포인트의 궤적",
           subtitle = "왼쪽 아래 사분면 = 복합반응") + THEME
  })
  output$dec_txt <- renderPrint({
    d <- cur(); i4 <- which.min(abs(d$week - 4)); i48 <- which.min(abs(d$week - 48))
    cat(sprintf(
"4주:  ALT는 이미 %.0f%% 감소, HDV RNA는 %.0f%%만 감소 (%+.2f log10)
48주: ALT %.0f%% 감소, HDV RNA %.0f%% 감소 (%+.2f log10)

이 해리는 모델의 '가설'이다: 손상 항에 kappa × (진입 플럭스)가 들어 있다.
반증 조건 — 바이러스학적 무반응자에서도 ALT가 (1 − 점유율)에 비례해 떨어져야 한다.
환자별로 ALT 감소가 HDV RNA 감소를 따라간다면 이 항은 틀렸다.
",
      100 * (1 - d$ALT[i4] / input$alt0), 100 * (1 - 10^d$dlog[i4]), d$dlog[i4],
      100 * (1 - d$ALT[i48] / input$alt0), 100 * (1 - 10^d$dlog[i48]), d$dlog[i48]))
  })

  ## ---- 6. ceiling & sensitivity ----
  output$ceiling <- renderPlot({
    doses <- c(0, 0.5, 1, 2, 5, 10, 20)
    df <- do.call(rbind, lapply(doses, function(mg) {
      d <- simulate_full(input, mg, 0, FALSE, 0, 0, 0, 48, 48)
      data.frame(dose = mg, occ = tail(d$OCCUP, 1), dlog = tail(d$dlog, 1))
    }))
    pf <- simulate_full(input, 0, 0, FALSE, 0, 0, 0, 48, 48, occfix = 1)
    ggplot(df, aes(occ, dlog)) +
      geom_line(colour = "#2f7a2f", linewidth = .9) +
      geom_point(size = 2, colour = "#2f7a2f") +
      geom_hline(yintercept = tail(pf$dlog, 1), linetype = 2, colour = "#b04a4a") +
      annotate("text", x = 0.05, y = tail(pf$dlog, 1) + 0.12,
               label = "완전 차단 = 이 계열의 천장", hjust = 0,
               size = 3.2, colour = "#b04a4a") +
      geom_text(aes(label = paste0(dose, " mg")), vjust = -0.8, size = 3) +
      labs(x = "NTCP 점유율", y = expression("48주 "*Delta*log[10]*" HDV RNA"),
           title = "점유율을 올려도 어디까지인가",
           subtitle = "천장은 약이 아니라 숙주 세포생물학이 정한다") + THEME
  })
  output$sens_tab <- renderTable({
    base <- tail(simulate_full(input, 2, 0, FALSE, 0, 0, 0, 48, 48)$dlog, 1)
    pars <- c(DIMM = "dimm", DHEP = "dhep", KCC = "ccmult",
              PHINT = "phint", KCURE = "kcure")
    rows <- lapply(names(pars), function(pn) {
      key <- pars[[pn]]
      v0 <- input[[key]]
      got <- sapply(c(0.7, 1.3), function(k) {
        inp <- reactiveValuesToList(input)
        inp[[key]] <- v0 * k
        m <- param(mod0, patient_pars(inp))
        m <- patient_init(m, inp)
        m <- update(m, end = 336, delta = 14)
        d <- as.data.frame(mrgsim(m, events = ev_blv(2, 336)))
        tail(d$LGVD, 1) - input$lgvd0
      })
      data.frame(파라미터 = pn, `-30%` = sprintf("%+.2f", got[1]),
                 기준 = sprintf("%+.2f", base),
                 `+30%` = sprintf("%+.2f", got[2]),
                 `평균 |영향|` = sprintf("%.3f", mean(abs(got - base))),
                 check.names = FALSE)
    })
    kd <- unlist(param(mod0))[["KDNTCP"]]
    gotkd <- sapply(c(0.7, 1.3), function(k) {
      m <- param(mod0, patient_pars(input)); m <- param(m, list(KDNTCP = kd * k))
      m <- patient_init(m, input); m <- update(m, end = 336, delta = 14)
      tail(as.data.frame(mrgsim(m, events = ev_blv(2, 336)))$LGVD, 1) - input$lgvd0
    })
    rows[[length(rows) + 1]] <- data.frame(
      파라미터 = "KDNTCP (표적 친화도)", `-30%` = sprintf("%+.2f", gotkd[1]),
      기준 = sprintf("%+.2f", base), `+30%` = sprintf("%+.2f", gotkd[2]),
      `평균 |영향|` = sprintf("%.3f", mean(abs(gotkd - base))), check.names = FALSE)
    do.call(rbind, rows)
  })

  ## ---- 7. bile acids ----
  output$ba_time <- renderPlot({
    d <- cur()
    ggplot(d, aes(week, TBAFLD)) + geom_line(colour = "#a08a20", linewidth = .9) +
      geom_hline(yintercept = 1, linetype = 3) +
      geom_vline(xintercept = input$txwk, linetype = 3) +
      labs(x = "주", y = "총담즙산 (기저 배수)",
           title = "총담즙산 — 무증상이지만 정보가 많다",
           subtitle = "가려움·담즙정체 없이 상승. 표적 점유율의 직접 판독값") + THEME
  })
  output$ba_dose <- renderPlot({
    doses <- c(0.5, 1, 2, 5, 10, 20)
    df <- do.call(rbind, lapply(doses, function(mg) {
      d <- simulate_full(input, mg, 0, FALSE, 0, 0, 0, 48, 48)
      data.frame(dose = mg, ba = tail(d$TBAFLD, 1),
                 eff = -tail(d$dlog, 1))
    }))
    ggplot(df, aes(dose)) +
      geom_line(aes(y = eff, colour = "효능 (−Δlog10 HDV, 48주)"), linewidth = .9) +
      geom_point(aes(y = eff, colour = "효능 (−Δlog10 HDV, 48주)"), size = 2) +
      geom_line(aes(y = ba / 5, colour = "비용 (담즙산 배수 ÷ 5)"), linewidth = .9) +
      geom_point(aes(y = ba / 5, colour = "비용 (담즙산 배수 ÷ 5)"), size = 2) +
      geom_vline(xintercept = 2, linetype = 2, colour = "grey40") +
      annotate("text", x = 2.3, y = 0.3, label = "승인 용량", hjust = 0, size = 3) +
      scale_x_log10(breaks = doses) +
      scale_colour_manual(values = c("#2f7a2f", "#a08a20"), name = NULL) +
      labs(x = "불레비티드 용량 (mg/일, log 축)", y = NULL,
           title = "효능은 포화하고 비용은 포화하지 않는다",
           subtitle = "2 mg가 무릎(knee)에 있다 — 표적을 포화시키지 않은 채로") + THEME
  })
  output$ba_tab <- renderTable({
    doses <- c(0.5, 1, 2, 5, 10, 20)
    do.call(rbind, lapply(doses, function(mg) {
      d <- simulate_full(input, mg, 0, FALSE, 0, 0, 0, 48, 48)
      i <- nrow(d)
      data.frame(`용량 (mg)` = mg,
                 `Css (nM)` = sprintf("%.2f", mean(tail(d$CBLV, 48))),
                 `NTCP 점유율` = sprintf("%.4f", d$OCCUP[i]),
                 `자유 NTCP (%)` = sprintf("%.1f", 100 * d$FREENT[i]),
                 `담즙산 (배)` = sprintf("%.1f", d$TBAFLD[i]),
                 `48주 Δlog10` = sprintf("%+.2f", d$dlog[i]),
                 `48주 ALT` = sprintf("%.0f", d$ALT[i]),
                 check.names = FALSE)
    }))
  })

  ## ---- 8. comparison ----
  cmp_data <- reactive({
    defs <- list(
      none    = list(0, 0, FALSE, 0, 0, 0, NA),
      blv2    = list(2, 0, FALSE, 0, 0, 0, NA),
      blv10   = list(10, 0, FALSE, 0, 0, 0, NA),
      ifn     = list(0, 180, FALSE, 0, 0, 0, NA),
      lam     = list(0, 180, TRUE, 0, 0, 0, NA),
      blvifn  = list(2, 180, FALSE, 0, 0, 0, NA),
      lnf     = list(0, 0, FALSE, 50, 100, 0, NA),
      lnfifn  = list(0, 180, FALSE, 50, 100, 0, NA),
      blvsir  = list(2, 0, FALSE, 0, 0, 200, NA),
      perfect = list(0, 0, FALSE, 0, 0, 0, 1.0)
    )
    labs <- c(none = "무치료", blv2 = "BLV 2 mg", blv10 = "BLV 10 mg",
              ifn = "Peg-IFN alfa", lam = "Peg-IFN lambda",
              blvifn = "BLV 2 + Peg-IFN", lnf = "LNF/RTV",
              lnfifn = "LNF/RTV + Peg-IFN", blvsir = "BLV 2 + siRNA",
              perfect = "완전 진입 차단")
    sel <- input$cmp
    if (!length(sel)) return(NULL)
    do.call(rbind, lapply(sel, function(k) {
      a <- defs[[k]]
      d <- simulate_full(input, a[[1]], a[[2]], a[[3]], a[[4]], a[[5]], a[[6]],
                         input$txwk, input$horiz, occfix = a[[7]])
      d$arm <- labs[[k]]
      d
    }))
  })

  output$cmp_plot <- renderPlot({
    d <- cmp_data(); if (is.null(d)) return(NULL)
    d$y <- d[[input$cmp_y]]
    ggplot(d, aes(week, y, colour = arm)) + geom_line(linewidth = .85) +
      geom_vline(xintercept = input$txwk, linetype = 3) +
      labs(x = "주", y = names(which(c(
        "Δlog10 HDV RNA" = "dlog", "ALT (U/L)" = "ALT",
        "HDAg 양성 간세포 (%)" = "HDAGPC", "세포내 HDV RNA (배수)" = "RGFOLD",
        "혈중 HBsAg (log10)" = "LGSAG", "총담즙산 (배수)" = "TBAFLD",
        "탈진 수준" = "EXH", "섬유화 (Ishak)" = "FIB",
        "혈소판 (10^9/L)" = "PLT", "호중구 (10^9/L)" = "NEU") == input$cmp_y)),
        title = "요법 비교", subtitle = "점선 = 치료 종료") +
      THEME + theme(legend.position = "right")
  })
  output$cmp_tab <- renderTable({
    d <- cmp_data(); if (is.null(d)) return(NULL)
    d %>% group_by(arm) %>%
      summarise(
        `치료종료 Δlog10` = sprintf("%+.2f", dlog[which.min(abs(week - input$txwk))]),
        `치료종료 ALT` = sprintf("%.0f", ALT[which.min(abs(week - input$txwk))]),
        `관찰종료 Δlog10` = sprintf("%+.2f", dlog[n()]),
        `관찰종료 ALT` = sprintf("%.0f", ALT[n()]),
        `관찰종료 Ishak` = sprintf("%.2f", FIB[n()]),
        `담즙산 (배)` = sprintf("%.1f", TBAFLD[which.min(abs(week - input$txwk))]),
        `호중구` = sprintf("%.2f", NEU[which.min(abs(week - input$txwk))]),
        .groups = "drop") %>% as.data.frame()
  })

  ## ---- 9. long-term outcomes ----
  long5 <- reactive({
    defs <- list(
      list("무치료", 0, 0, FALSE, 0, 0, 0),
      list("BLV 2 mg 지속", 2, 0, FALSE, 0, 0, 0),
      list("BLV 10 mg 지속", 10, 0, FALSE, 0, 0, 0),
      list("BLV 2 + Peg-IFN 48주", 2, 180, FALSE, 0, 0, 0)
    )
    do.call(rbind, lapply(defs, function(a) {
      d <- simulate_full(input, a[[2]], a[[3]], a[[4]], a[[5]], a[[6]], a[[7]],
                         48, 260)
      d$arm <- a[[1]]; d
    }))
  })
  output$out_fib <- renderPlot({
    ggplot(long5(), aes(week / 52, FIB, colour = arm)) + geom_line(linewidth = .85) +
      geom_hline(yintercept = 5, linetype = 2, colour = "grey40") +
      annotate("text", x = 0.1, y = 5.15, label = "간경변", hjust = 0, size = 3) +
      labs(x = "년", y = "Ishak 섬유화 단계", title = "섬유화 진행/회귀") + THEME
  })
  output$out_hcc <- renderPlot({
    ggplot(long5(), aes(week / 52, HCCINC, colour = arm)) +
      geom_line(linewidth = .85) +
      labs(x = "년", y = "누적 HCC 발생률 (%)", title = "간세포암 위험") + THEME
  })
  output$out_plt <- renderPlot({
    ggplot(long5(), aes(week / 52, PLT, colour = arm)) + geom_line(linewidth = .85) +
      labs(x = "년", y = "혈소판 (10^9/L)",
           title = "혈소판 = 문맥압 대리지표") + THEME
  })
  output$out_alt <- renderPlot({
    ggplot(long5(), aes(week / 52, ALT, colour = arm)) + geom_line(linewidth = .85) +
      geom_hline(yintercept = 40, linetype = 2, colour = "grey40") +
      labs(x = "년", y = "ALT (U/L)",
           title = "섬유화를 움직이는 것은 ALT다") + THEME
  })

  output$about <- renderPrint({
    cat("
만성 D형 간염 (HDV) QSP 모델 — 33개 ODE 구획
==============================================================================
중심 아이디어
  HDV는 '죽여야 할 바이러스'가 아니라 두 기질 조립 라인이다.
    기질 1  HDV 유전체  — 숙주 RNA Pol II가 만든다 (HBV 복제와 무관)
    기질 2  HBsAg 외피  — HBV cccDNA가 만든다 (같은 세포 또는 주변 HBsAg+ 풀)
  그래서 약은 '어떤 플럭스를 끊는가'로 분류된다.
    (E) 진입/재감염     NTCP 의존       불레비티드
    (A) 외피화/방출     FTase, HBsAg    로나파닙, siRNA/NAP
    (C) 세포내 + 세포소실 Pol II, 면역   Peg-IFN alfa/lambda
  그리고 감염세포 풀은 NTCP와 무관한 두 경로로도 유지된다.
    분열매개 전파  — 감염세포가 분열하면 딸세포 둘 다 감염된 상태다
    세포간 전파    — 일부만 NTCP 의존
  이 둘이 진입억제제가 넘을 수 없는 '바닥'이다.

적합된 파라미터는 4개뿐이다
  r_oatp, Kd_ntcp   총담즙산 상승 배수 2개 (2 mg ×3.2, 10 mg ×13.0)
                    → NTCP 점유율은 '도출'된다: 2 mg 0.710, 10 mg 0.954
  dth_Id_immune     MYR301 2 mg 48주 바이러스 반응률 71% (적합값 0.02636/day)
  ALT_base          MYR301 2 mg 48주 ALT 정상화율 51%
나머지는 문헌값이거나 무치료 정상상태 조건으로부터 역산된 값이다.

예측(적합에 쓰이지 않은 것): 10 mg arm 전체, 무치료 대조군, 96주, Peg-IFN의
치료 중/종료 후 반응, 병용 시너지, 로나파닙/리토나비르, Peg-IFN 람다,
진입억제 계열의 천장, 담즙산 치료지수 곡선, 섬유화·HCC 5년 추정.

솔직한 불일치: 담즙산 앵커에 묶인 상태에서 이 모델은 바이러스학적
엔드포인트에서 2 mg와 10 mg를 실제 시험보다 크게 벌린다. 세 가지 가능한
해석과 그중 검증 가능한 하나가 hdv_model_report.txt A13절에 있다.

파일
  hdv_qsp_model.dot/.svg/.png  기계론적 지도 (18 클러스터, 157 노드)
  hdv_mrgsolve_model.R         이 앱이 쓰는 ODE 모델 (33 구획, 14 시나리오)
  hdv_reference_model.py       의존성 없는 파이썬 참조 구현 + 가상 집단
  hdv_model_report.txt         위 모든 수치의 계산 결과 (A0–A14)
  hdv_references.md            파라미터별 출처 (PubMed 확인 완료 항목 표시)

교육·연구 목적의 모델이며 임상 의사결정에 사용할 수 없다.
")
  })
}

shinyApp(ui, server)
