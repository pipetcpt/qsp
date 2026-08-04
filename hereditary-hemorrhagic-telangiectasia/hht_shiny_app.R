## ============================================================================
##  HHT (유전성 출혈성 모세혈관확장증) QSP — Shiny 대시보드
##  Hereditary Hemorrhagic Telangiectasia · interactive QSP explorer
##
##  10 tabs:
##   1. 환자 프로파일     patient profile & natural history
##   2. 약동학 (PK)       five drugs, including the nasal route
##   3. 병변 생물학       pSMAD · VEGF · mural coverage · lesion number and calibre
##   4. ★ 전단응력 되먹임  the sign flip, the feeder ceiling, escape thresholds
##   5. 출혈 · 코출혈      episodes, duration, rate, monthly duration
##   6. ★ 철 · 빈혈        the iron ceiling — where clinical benefit is threshold-like
##   7. 심장 · 션트        cardiac index, shunt flow, the anaemia feedback loop
##   8. 임상 엔드포인트    ESS, HHT-QOL, transfusion burden, iron independence
##   9. 시나리오 비교      18 treatment scenarios side by side
##  10. 안전성            drug-specific adverse effects
##
##  실행 (run):
##    setwd("hereditary-hemorrhagic-telangiectasia"); shiny::runApp("hht_shiny_app.R")
##  필요 패키지: shiny, mrgsolve, ggplot2, dplyr, tidyr, DT
## ============================================================================

library(shiny)
library(mrgsolve)
library(ggplot2)
library(dplyr)
library(tidyr)

## the model, phenotypes, scenarios and helpers all live in the model file
source("hht_mrgsolve_model.R")

THEME <- theme_bw(base_size = 12) +
  theme(panel.grid.minor = element_blank(),
        strip.background = element_rect(fill = "grey93", colour = NA),
        legend.position = "bottom")

PAL <- c("#2471a3", "#c0392b", "#1a7f37", "#8e44ad", "#e67e22",
         "#16a085", "#7f8c8d", "#d35400")

## ---------------------------------------------------------------------------
ui <- fluidPage(
  titlePanel("HHT (유전성 출혈성 모세혈관확장증) QSP 대시보드"),
  tags$p(style = "color:#555;margin-top:-8px;",
         "ALK1/ENG/SMAD4는 내피 전단응력 설정점 조절자다 — 이를 잃으면 되먹임의 부호가 뒤집힌다. ",
         "임상적 이득은 연속적이지 않고 철 흡수 상한을 넘느냐로 갈린다."),
  sidebarLayout(
    sidebarPanel(
      width = 3,
      h4("환자 (Patient)"),
      selectInput("pheno", "표현형 (phenotype)",
                  choices = c("경증 mild" = "mild",
                              "중등증 moderate (PATH-HHT)" = "moderate",
                              "중증·수혈의존 severe (Parambil)" = "severe",
                              "간 AVM 고박출 hepatic (Dupuis-Girod)" = "hepatic",
                              "HHT1 폐 AVM" = "hht1_pav",
                              "JP-HHT (SMAD4)" = "jphht"),
                  selected = "moderate"),
      sliderInput("age", "치료 시작 연령 (years)", 20, 75, 52, step = 1),
      sliderInput("wt", "체중 (kg)", 45, 120, 70, step = 1),
      hr(),
      h4("표현형 조절 (phenotype knobs)"),
      sliderInput("sev",  "체성 타격률 SEV (병변 개수)", 0.1, 2.5, 1.05, step = 0.05),
      sliderInput("angf", "혈관신생 구동 ANGF (병변 크기)", 0.4, 1.6, 1.02, step = 0.02),
      sliderInput("gif",  "위장 병변 부하 GIF", 0, 50, 0.55, step = 0.25),
      sliderInput("hepf", "간 침범 HEPF", 0, 1, 0.20, step = 0.05),
      sliderInput("pulf", "폐 침범 PULF", 0, 1, 0.00, step = 0.05),
      hr(),
      h4("치료 (Treatment)"),
      checkboxGroupInput(
        "drugs", NULL,
        choices = c("베바시주맙 IV (유지)"    = "bev",
                    "베바시주맙 비강분무"      = "bevnas",
                    "파조파닙"                = "paz",
                    "포말리도마이드"           = "pom",
                    "트라넥삼산 (경구)"        = "txa",
                    "ALK1 경로 회복 (가설)"    = "alk1"),
        selected = character(0)),
      conditionalPanel("input.drugs.indexOf('paz') > -1",
        sliderInput("pazdose", "파조파닙 (mg/일)", 25, 800, 150, step = 25)),
      conditionalPanel("input.drugs.indexOf('pom') > -1",
        sliderInput("pomdose", "포말리도마이드 (mg/일)", 1, 4, 4, step = 0.5)),
      conditionalPanel("input.drugs.indexOf('bevnas') > -1",
        sliderInput("bevnasdose", "비강 분무 1회 용량 (mg)", 25, 1000, 75, step = 25)),
      selectInput("iron", "철 보충 (iron support)",
                  choices = c("없음 none" = "none",
                              "경구 65 mg/일" = "oral65",
                              "경구 200 mg/일" = "oral200",
                              "경구 + IV 1 g/8주" = "iv8",
                              "경구 + IV 1 g/4주" = "iv4",
                              "경구 + IV 1 g/2주" = "iv2"),
                  selected = "oral200"),
      sliderInput("dur", "시뮬레이션 기간 (일)", 90, 1095, 365, step = 15),
      actionButton("go", "실행 (Run)", class = "btn-primary", width = "100%")
    ),
    mainPanel(
      width = 9,
      tabsetPanel(
        id = "tabs",
        tabPanel("1. 환자 프로파일", br(),
                 plotOutput("p_nat", height = "420px"),
                 h5("자연 경과 — 병변은 수십 년에 걸쳐 축적되고 50-60대에 정체됩니다."),
                 tableOutput("t_base")),
        tabPanel("2. 약동학 (PK)", br(),
                 plotOutput("p_pk", height = "440px"),
                 helpText("베바시주맙은 mg/kg 정맥, 파조파닙·포말리도마이드·트라넥삼산은 경구. ",
                          "비강 경로는 점액섬모청소(t½ ~15분) 때문에 최고 농도는 높아도 ",
                          "역치 이상 시간이 무시할 수준입니다.")),
        tabPanel("3. 병변 생물학", br(),
                 plotOutput("p_bio", height = "480px"),
                 helpText("pSMAD가 사라지면 VEGF·AKT가 올라가고 벽세포 피복이 무너집니다. ",
                          "취약성 FRAG = (1 − 피복률)^1.55 이 출혈 빈도를 결정합니다.")),
        tabPanel("4. ★ 전단응력 되먹임", br(),
                 fluidRow(column(6, plotOutput("p_wss", height = "360px")),
                          column(6, plotOutput("p_remod", height = "360px"))),
                 h5("부호 반전 (the sign flip)"),
                 helpText("REMOD = KSH · g(WSS) · (1 − 2·pSMAD).  pSMAD ≈ 1 이면 음수(안쪽 재형성, ",
                          "자기제한). pSMAD ≈ 0 이면 양수(바깥쪽 재형성, 자기증폭). ",
                          "전단응력은 내강 크기에 단조증가하지 않습니다 — 공급혈관이 저항을 지배하면 ",
                          "다시 떨어집니다. 그래서 코 병변은 안정 내강을 갖고 간·폐 AVM은 갖지 않습니다."),
                 tableOutput("t_escape")),
        tabPanel("5. 출혈 · 코출혈", br(),
                 plotOutput("p_bleed", height = "460px"),
                 helpText("트라넥삼산은 응괴에만 작용하므로 지속시간은 줄이지만 횟수는 바꾸지 못합니다 ",
                          "(ATERO 시험에서 관찰된 해리 — 이 모델에서는 적합이 아니라 구조적 예측입니다).")),
        tabPanel("6. ★ 철 · 빈혈", br(),
                 plotOutput("p_iron", height = "380px"),
                 h5("철 흡수 상한 — 임상적 문턱값"),
                 tableOutput("t_ceiling"),
                 helpText("정상 상태는 (흡수 − 의무손실) = 0.0347 × Hb × 출혈률 을 요구합니다. ",
                          "헵시딘이 하루 흡수량을 ~10 mg으로 묶기 때문에 경구 철분 200 mg/일은 ",
                          "65 mg/일보다 나은 것이 없습니다. 이 선을 넘느냐가 수혈 의존을 가릅니다.")),
        tabPanel("7. 심장 · 션트", br(),
                 plotOutput("p_card", height = "440px"),
                 helpText("빈혈은 결과이자 입력입니다: Hb↓ → 산소운반 방어를 위한 심박출↑ → ",
                          "관류압·전단응력↑ → 병변 확대 → 출혈↑. 이 고리 때문에 빈혈 교정은 ",
                          "단순한 보존치료가 아니라 병변 성장의 기계적 구동을 걷어내는 조치입니다.")),
        tabPanel("8. 임상 엔드포인트", br(),
                 plotOutput("p_end", height = "440px"),
                 tableOutput("t_end"),
                 helpText("ESS의 MCID는 0.71점입니다 (PMID 26393959). ",
                          "그러나 연속 점수는 철 자립이라는 문턱 사건을 보여주지 못합니다.")),
        tabPanel("9. 시나리오 비교", br(),
                 actionButton("go_scen", "18개 시나리오 실행", class = "btn-success"),
                 br(), br(),
                 plotOutput("p_scen", height = "520px"),
                 tableOutput("t_scen")),
        tabPanel("10. 안전성", br(),
                 plotOutput("p_ae", height = "440px"),
                 helpText("항VEGF: 고혈압·단백뇨·상처치유 지연. 파조파닙: 간효소 상승. ",
                          "포말리도마이드: 호중구감소·변비·발진·정맥혈전. ",
                          "모델은 파조파닙 고용량에서 PDGFRβ 차단에 의한 벽세포 손실로 ",
                          "순효과가 역전되는 것을 예측합니다."))
      )
    )
  )
)

## ---------------------------------------------------------------------------
server <- function(input, output, session) {

  observeEvent(input$pheno, {
    p <- PHENO[[input$pheno]]
    updateSliderInput(session, "sev",  value = p$SEV)
    updateSliderInput(session, "angf", value = p$ANGF)
    updateSliderInput(session, "gif",  value = p$GIF)
    updateSliderInput(session, "hepf", value = p$HEPF)
    updateSliderInput(session, "pulf", value = if (is.null(p$PULF)) 0 else p$PULF)
  })

  pars <- reactive({
    list(SEV = input$sev, ANGF = input$angf, GIF = input$gif,
         HEPF = input$hepf, PULF = input$pulf, WT = input$wt)
  })

  iron_par <- reactive({
    switch(input$iron,
           none    = list(FE_ORAL = 0,   FE_IV_RATE = 0),
           oral65  = list(FE_ORAL = 65,  FE_IV_RATE = 0),
           oral200 = list(FE_ORAL = 200, FE_IV_RATE = 0),
           iv8     = list(FE_ORAL = 200, FE_IV_RATE = 1000 / 56),
           iv4     = list(FE_ORAL = 200, FE_IV_RATE = 1000 / 28),
           iv2     = list(FE_ORAL = 200, FE_IV_RATE = 1000 / 14))
  })

  ## ---- natural history to the chosen age ---------------------------------
  natural <- eventReactive(input$go, {
    mod %>% param(pars()) %>%
      mrgsim(end = 365 * input$age, delta = 365, atol = 1e-8, rtol = 1e-6) %>%
      as.data.frame()
  }, ignoreNULL = FALSE)

  base_state <- reactive({
    df <- natural()
    as.list(df[nrow(df), names(init(mod))])
  })

  ## ---- treatment events ---------------------------------------------------
  events <- reactive({
    d <- input$drugs
    e <- NULL
    add <- function(x) if (is.null(e)) x else c(e, x)
    n <- input$dur
    if ("bev"    %in% d) e <- add(ev(amt = 5 * input$wt, cmt = 1, ii = 28, addl = floor(n / 28)))
    if ("bevnas" %in% d) e <- add(ev(amt = input$bevnasdose, cmt = 3, ii = 14, addl = floor(n / 14)))
    if ("paz"    %in% d) e <- add(ev(amt = input$pazdose, cmt = 5, ii = 1, addl = n - 1))
    if ("pom"    %in% d) e <- add(ev(amt = input$pomdose, cmt = 7, ii = 1, addl = n - 1))
    if ("txa"    %in% d) e <- add(ev(amt = 1000, cmt = 9, ii = 1/3, addl = 3 * n - 1))
    if (is.null(e)) ev(amt = 0, cmt = 1, time = 0) else e
  })

  sim <- eventReactive(input$go, {
    pp <- c(pars(), iron_par())
    if ("alk1" %in% input$drugs) pp$ALK1AG <- 0.60
    mod %>% param(pp) %>% init(base_state()) %>%
      mrgsim(events = events(), end = input$dur, delta = 1,
             atol = 1e-10, rtol = 1e-8) %>% as.data.frame()
  }, ignoreNULL = FALSE)

  ## ---- 1. natural history --------------------------------------------------
  output$p_nat <- renderPlot({
    df <- natural() %>% mutate(age = time / 365) %>%
      select(age, `병변 개수 N_n` = N_n, `평균 내강 S_n` = S_n,
             `헤모글로빈 Hb` = HB, `출혈률 mL/일` = BLRo,
             `ESS` = ESSo, `심계수 CI` = CIo) %>%
      pivot_longer(-age)
    ggplot(df, aes(age, value)) +
      geom_line(colour = PAL[1], linewidth = 0.9) +
      facet_wrap(~name, scales = "free_y") +
      labs(x = "연령 (years)", y = NULL) + THEME
  })

  output$t_base <- renderTable({
    df <- natural(); r <- df[nrow(df), ]
    data.frame(
      `지표` = c("ESS", "월간 코출혈 (분)", "월 코출혈 횟수", "사건당 지속 (분)",
                 "총 출혈률 (mL/일)", "헤모글로빈 (g/dL)", "심계수", "SpO2"),
      `값` = c(round(r$ESSo, 2), round(r$MEDo, 1), round(r$EPMo, 1),
               round(r$dur_evo, 1), round(r$BLRo, 1), round(r$HB, 2),
               round(r$CIo, 2), round(r$SPO2o, 3)),
      check.names = FALSE)
  })

  ## ---- 2. PK ---------------------------------------------------------------
  output$p_pk <- renderPlot({
    df <- sim() %>%
      select(time, `베바시주맙 (µg/mL)` = C_BEVo,
             `파조파닙 유리 (ng/mL)` = CU_PAZo,
             `포말리도마이드 (ng/mL)` = C_POMo,
             `트라넥삼산 (µg/mL)` = C_TXAo,
             `비강 점막하 (mg)` = BEVN_T) %>%
      pivot_longer(-time) %>% filter(value > 0)
    validate(need(nrow(df) > 0, "선택된 약물이 없습니다 — 왼쪽에서 치료를 선택하고 실행하세요."))
    ggplot(df, aes(time, value, colour = name)) +
      geom_line(linewidth = 0.9) +
      facet_wrap(~name, scales = "free_y") +
      scale_colour_manual(values = PAL, guide = "none") +
      labs(x = "시간 (일)", y = "농도") + THEME
  })

  ## ---- 3. lesion biology ---------------------------------------------------
  output$p_bio <- renderPlot({
    df <- sim() %>%
      select(time, `pSMAD1/5 (비강)` = PS_n, `ID1 전사체` = ID1_n,
             `조직 VEGF-A (pg/mL)` = VEGF_n, `PI3K/AKT` = AKT_n,
             `PDGF-B` = PDGFB_n, `벽세포 피복 MU_n` = MU_n,
             `병변 개수 N_n` = N_n, `평균 내강 S_n` = S_n) %>%
      pivot_longer(-time)
    ggplot(df, aes(time, value)) +
      geom_line(colour = PAL[3], linewidth = 0.9) +
      facet_wrap(~name, scales = "free_y", ncol = 4) +
      labs(x = "시간 (일)", y = NULL) + THEME
  })

  ## ---- 4. shear feedback ---------------------------------------------------
  output$p_wss <- renderPlot({
    p <- as.list(param(mod))
    S <- seq(0.5, 50, length.out = 400)
    beds <- data.frame(
      bed = rep(c("비강 nasal", "위장 GI", "간 hepatic", "폐 pulmonary"), each = length(S)),
      S = rep(S, 4),
      w0 = rep(c(p$WSS0_n, p$WSS0_g, p$WSS0_h, p$WSS0_p), each = length(S)),
      ss = rep(c(p$SSAT_n, p$SSAT_g, p$SSAT_h, p$SSAT_p), each = length(S)))
    beds$WSS <- beds$w0 * beds$S / (1 + (beds$S / beds$ss)^4)
    ggplot(beds, aes(S, WSS, colour = bed)) +
      geom_line(linewidth = 1) +
      scale_x_log10() +
      scale_colour_manual(values = PAL) +
      labs(title = "전단응력은 내강에 단조증가하지 않는다",
           subtitle = "WSS ∝ S/(1+(S/S_feeder)⁴) — 공급혈관이 저항을 지배하면 다시 떨어진다",
           x = "내강 지수 S (log)", y = "벽면 전단응력") + THEME
  })

  output$p_remod <- renderPlot({
    p <- as.list(param(mod))
    ps <- seq(0, 1, length.out = 200)
    wss <- c(0.2, 0.5, 1.0, 3.0)
    df <- expand.grid(pSMAD = ps, WSS = wss)
    df$REMOD <- p$KSH * (df$WSS / (p$KWSS + df$WSS)) * (1 - 2 * df$pSMAD)
    df$WSSlab <- factor(paste0("WSS = ", df$WSS))
    ggplot(df, aes(pSMAD, REMOD, colour = WSSlab)) +
      geom_hline(yintercept = 0, linetype = 2, colour = "grey40") +
      geom_vline(xintercept = 0.5, linetype = 3, colour = "grey40") +
      geom_line(linewidth = 1) +
      annotate("text", x = 0.15, y = max(df$REMOD) * 0.85,
               label = "바깥쪽 (outward)\n자기증폭", size = 3.4, colour = "#c0392b") +
      annotate("text", x = 0.85, y = min(df$REMOD) * 0.85,
               label = "안쪽 (inward)\n자기제한", size = 3.4, colour = "#1a7f37") +
      scale_colour_manual(values = PAL, name = NULL) +
      labs(title = "되먹임의 부호 반전",
           subtitle = "REMOD = KSH · g(WSS) · (1 − 2·pSMAD)",
           x = "상대 pSMAD1/5 활성", y = "재형성 항") + THEME
  })

  output$t_escape <- renderTable({
    p <- as.list(param(mod))
    ps <- p$GD_LES * p$BMP9 / (p$KD9 + p$BMP9)
    beds <- list(n = c(p$WSS0_n, p$SSAT_n, 0.52), g = c(p$WSS0_g, p$SSAT_g, 0.52),
                 h = c(p$WSS0_h, p$SSAT_h, 0), p = c(p$WSS0_p, p$SSAT_p, 0),
                 c = c(p$WSS0_c, p$SSAT_c, 0))
    nm <- c(n = "비강", g = "위장", h = "간", p = "폐", c = "뇌")
    do.call(rbind, lapply(c(1.0, 1.45, 2.0, 2.6), function(co) {
      row <- sapply(names(beds), function(b) {
        v <- beds[[b]]
        Sp <- v[2] / 3^0.25                       # calibre where shear peaks
        w <- v[1] * co * Sp / (1 + (Sp / v[2])^4)
        remod <- p$KSH * (w / (p$KWSS + w)) * (1 - 2 * ps)
        sink <- p$KMAT * v[3] + p$KDEC
        sprintf("%s (%.2f)", ifelse(remod > sink, "탈출", "억제"), remod / sink)
      })
      setNames(data.frame(sprintf("%.2f", co), t(row), stringsAsFactors = FALSE),
               c("CO_rel", nm[names(beds)]))
    }))
  }, caption = "완전한 항혈관신생 약물로 구동을 모두 제거했을 때, 전단응력만으로 병변이 유지되는가? (비 > 1 = 약물 불응)",
     caption.placement = "top")

  ## ---- 5. bleeding ---------------------------------------------------------
  output$p_bleed <- renderPlot({
    df <- sim() %>%
      select(time, `월간 코출혈 (분)` = MEDo, `월 코출혈 횟수` = EPMo,
             `사건당 지속 (분)` = dur_evo, `사건당 속도 (mL/분)` = rate_evo,
             `비강 출혈률 (mL/일)` = BLR_No, `위장 출혈률 (mL/일)` = BLR_Go) %>%
      pivot_longer(-time)
    ggplot(df, aes(time, value)) +
      geom_line(colour = PAL[2], linewidth = 0.9) +
      facet_wrap(~name, scales = "free_y") +
      labs(x = "시간 (일)", y = NULL) + THEME
  })

  ## ---- 6. iron -------------------------------------------------------------
  output$p_iron <- renderPlot({
    df <- sim() %>%
      select(time, `헤모글로빈 (g/dL)` = HB, `철 저장 (mg)` = FES,
             `헵시딘` = HEP, `출혈 철 손실 (mg/일)` = FE_LOSSo,
             `흡수 철 (mg/일)` = FE_ABSo,
             `지속가능 최대 출혈 (mL/일)` = FE_CEILo,
             `실제 출혈률 (mL/일)` = BLRo) %>%
      pivot_longer(-time)
    ggplot(df, aes(time, value)) +
      geom_line(colour = PAL[5], linewidth = 0.9) +
      facet_wrap(~name, scales = "free_y") +
      labs(x = "시간 (일)", y = NULL) + THEME
  })

  output$t_ceiling <- renderTable({ iron_ceiling() },
    caption = "각 철 보충 전략이 감당할 수 있는 최대 출혈률 (Hb 14.6 g/dL 기준)",
    caption.placement = "top")

  ## ---- 7. cardiac ----------------------------------------------------------
  output$p_card <- renderPlot({
    df <- sim() %>%
      select(time, `심계수 (L/분/m²)` = CIo, `심박출량 (L/분)` = CO,
             `간 션트 지수 S_h` = S_h, `폐 션트 지수 S_p` = S_p,
             `SpO2` = SPO2o, `좌심실 재형성` = LVR,
             `우심방압 지수` = RAP, `헤모글로빈` = HB) %>%
      pivot_longer(-time)
    ggplot(df, aes(time, value)) +
      geom_line(colour = PAL[4], linewidth = 0.9) +
      facet_wrap(~name, scales = "free_y", ncol = 4) +
      labs(x = "시간 (일)", y = NULL) + THEME
  })

  ## ---- 8. endpoints --------------------------------------------------------
  output$p_end <- renderPlot({
    df <- sim() %>%
      select(time, `ESS (0-10)` = ESSo, `HHT 삶의 질 (0-16)` = QOL,
             `누적 출혈 (mL)` = CUM_BL, `누적 수혈 (단위 등가)` = CUM_TX,
             `철 자립 (1=예)` = IRONINDo, `헤모글로빈` = HB) %>%
      pivot_longer(-time)
    ggplot(df, aes(time, value)) +
      geom_line(colour = PAL[1], linewidth = 0.9) +
      facet_wrap(~name, scales = "free_y") +
      labs(x = "시간 (일)", y = NULL) + THEME
  })

  output$t_end <- renderTable({
    d <- sim(); a <- d[1, ]; b <- d[nrow(d), ]
    data.frame(
      `지표` = c("ESS", "HHT-QOL", "월간 코출혈 (분)", "총 출혈률 (mL/일)",
                 "헤모글로빈 (g/dL)", "심계수", "철 자립"),
      `기저` = c(round(a$ESSo, 2), round(a$QOL, 2), round(a$MEDo, 1),
                 round(a$BLRo, 1), round(a$HB, 2), round(a$CIo, 2),
                 ifelse(a$IRONINDo > 0, "예", "아니오")),
      `종료` = c(round(b$ESSo, 2), round(b$QOL, 2), round(b$MEDo, 1),
                 round(b$BLRo, 1), round(b$HB, 2), round(b$CIo, 2),
                 ifelse(b$IRONINDo > 0, "예", "아니오")),
      `변화` = c(round(b$ESSo - a$ESSo, 2), round(b$QOL - a$QOL, 2),
                 round(b$MEDo - a$MEDo, 1), round(b$BLRo - a$BLRo, 1),
                 round(b$HB - a$HB, 2), round(b$CIo - a$CIo, 2), ""),
      check.names = FALSE)
  })

  ## ---- 9. scenarios --------------------------------------------------------
  scen <- eventReactive(input$go_scen, {
    withProgress(message = "18개 시나리오 실행 중…", value = 0.5, {
      run_scenarios(input$pheno, input$age, input$dur)
    })
  })

  output$p_scen <- renderPlot({
    df <- scen() %>% mutate(scenario = reorder(scenario, dESS))
    ggplot(df, aes(dESS, scenario, fill = dESS)) +
      geom_col() +
      geom_vline(xintercept = -0.71, linetype = 2, colour = "#c0392b") +
      annotate("text", x = -0.71, y = 1, label = "MCID 0.71", hjust = -0.08,
               size = 3.2, colour = "#c0392b") +
      scale_fill_gradient(low = "#1a7f37", high = "#c0392b", guide = "none") +
      labs(x = "ESS 변화 (음수 = 개선)", y = NULL,
           title = "치료 시나리오 비교", subtitle = "붉은 점선 = 최소 임상적 의미 차이 0.71점") +
      THEME
  })

  output$t_scen <- renderTable({ scen() })

  ## ---- 10. safety ----------------------------------------------------------
  output$p_ae <- renderPlot({
    df <- sim() %>%
      select(time, `수축기 혈압 (mmHg)` = SBP, `단백뇨 (UPCR)` = UPCR,
             `ALT (U/L)` = ALT, `호중구 (10⁹/L)` = ANC,
             `누적 VTE 위험` = VTEH, `벽세포 피복 MU_n` = MU_n) %>%
      pivot_longer(-time)
    ggplot(df, aes(time, value)) +
      geom_line(colour = PAL[8], linewidth = 0.9) +
      facet_wrap(~name, scales = "free_y") +
      labs(x = "시간 (일)", y = NULL) + THEME
  })
}

shinyApp(ui, server)
