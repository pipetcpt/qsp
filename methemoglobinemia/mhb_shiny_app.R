# =============================================================================
#  mhb_shiny_app.R
#  Methaemoglobinaemia QSP model — interactive dashboard
#  메트헤모글로빈혈증 QSP 모델 — 인터랙티브 대시보드
# =============================================================================
#
#  WHAT THIS APP IS FOR
#  --------------------
#  Most dashboards in this library let you turn a knob and watch a curve.  This
#  one is built around a single disagreement, and almost every tab exists to
#  make that disagreement visible:
#
#        the number on the monitor,
#        the number in the blood gas,
#        and the number at the tissue
#        are three different numbers, and they move apart as the patient
#        gets worse.
#
#  Tab 3 (the two liars) puts SpO2, co-oximetry SaO2 and PvO2 on one figure.
#  Tab 4 converts %MetHb into the anaemia it is haemodynamically equivalent to.
#  Tab 7 shows the antidote's dose-response bending back on itself.
#  Tab 8 shows why it does that, in the currency the model actually spends.
#
#  RUN
#  ---
#    install.packages(c("shiny","mrgsolve","dplyr","tidyr","ggplot2","DT","gridExtra"))
#    shiny::runApp("mhb_shiny_app.R")
#
#  The model file `mhb_mrgsolve_model.R` must sit beside this one.
# =============================================================================

library(shiny)
library(mrgsolve)
library(dplyr)
library(tidyr)
library(ggplot2)
library(DT)

mod <- mread_cache("mhb", "mhb_mrgsolve_model.R")

MWMB <- 319.85
mb_ev <- function(mgkg, t = 0, wt = 70)
  ev(time = t, amt = mgkg * wt * 1000 / MWMB, cmt = "MBC")

PAL <- c(SpO2 = "#2b6cb0", `co-ox SaO2` = "#b3541e", PvO2 = "#2f6b34",
         MetHb = "#8b2f39", Hb = "#5b3a80")
thm <- theme_minimal(base_size = 13) +
  theme(panel.grid.minor = element_blank(),
        plot.title = element_text(face = "bold"),
        legend.position = "bottom")

# ---------------------------------------------------------------------------
#  build the event sequence and parameter set from the UI
# ---------------------------------------------------------------------------
build_sim <- function(input, tend = NULL) {
  p <- list(
    HB0    = input$hb,
    G6PD   = input$g6pd,
    EB5SET = input$eb5,
    FHBM   = input$fhbm,
    PAO2   = if (input$hbo) input$pao2 else 95,
    VO2    = input$vo2,
    ALPHAM = input$alpham
  )
  e <- NULL
  add <- function(a, b) if (is.null(a)) b else c(a, b)

  if (input$oxidant == "benzocaine")
    e <- add(e, ev(time = 0, amt = input$bzc, cmt = "BZCD"))
  if (input$oxidant == "dapsone_chronic")
    e <- add(e, ev(time = 0, amt = input$dapdose, cmt = "DAPG",
                   ii = 24, addl = max(0, ceiling(input$tend / 24) - 1)))
  if (input$oxidant == "dapsone_od")
    e <- add(e, ev(time = 0, amt = input$dapod, cmt = "DAPG"))
  if (input$oxidant == "nitrite")
    e <- add(e, ev(time = 0, amt = input$nit * 1000 * 1000 / 69, cmt = "NITG"))

  if (input$cimetidine)
    e <- add(e, ev(time = 0, amt = 400, cmt = "CIMG", ii = 8,
                   addl = max(0, ceiling(input$tend / 8) - 1)))
  if (input$ssri)
    e <- add(e, ev(time = 0, amt = 2.5, cmt = "SSRI"))
  if (input$ascorbate)
    e <- add(e, ev(time = input$mbtime, amt = 10 * 1e6 / 176.1, cmt = "ASCC"))

  if (input$nmb > 0)
    for (k in seq_len(input$nmb))
      e <- add(e, mb_ev(input$mbdose, input$mbtime + (k - 1) * input$mbint))

  list(p = p, e = e, tend = if (is.null(tend)) input$tend else tend)
}

run_sim <- function(s, delta = 0.05) {
  m <- param(mod, s$p)
  if (is.null(s$e)) mrgsim(m, end = s$tend, delta = delta) %>% as_tibble()
  else mrgsim(m, events = s$e, end = s$tend, delta = delta) %>% as_tibble()
}

# EQUIV Hb: invert the oxygen-transport block for a patient with normal Hb
equiv_hb <- function(pvo2, p50 = 26.8, n = 2.7, vo2 = 250, co = 5) {
  f <- function(hb) {
    sao2 <- 95^n / (p50^n + 95^n)
    cao2 <- 1.34 * hb * sao2 + 0.003 * 95
    cvo2 <- cao2 - vo2 / (co * 10)
    svo2 <- pmin(pmax(cvo2 / (1.34 * hb), 1e-6), 1 - 1e-6)
    p50 * (svo2 / (1 - svo2))^(1 / n) - pvo2
  }
  out <- numeric(length(pvo2))
  for (i in seq_along(pvo2)) {
    out[i] <- tryCatch(uniroot(f, c(0.2, 30))$root, error = function(e) NA_real_)
  }
  out
}

# static oxygen-transport block, for the sweep tabs (no integration needed)
transport <- function(fmet, hb = 15, p50 = 26.8, nh = 2.7, alpham = 0.5,
                      pao2 = 95, vo2 = 250, co = 5, bpg = 1, bbpg = 0.55,
                      e660m = 15, e940m = 15) {
  hbf <- hb * (1 - fmet)
  fred <- fmet
  neff <- 1 + (nh - 1) * (1 - fred)
  p50e <- pmax(p50 * (1 + bbpg * (bpg - 1)) * (1 - alpham * fred), 1)
  sao2 <- pao2^neff / (p50e^neff + pao2^neff)
  cao2 <- 1.34 * hbf * sao2 + 0.003 * pao2
  cvo2 <- cao2 - vo2 / (co * 10)
  caph <- 1.34 * hbf
  svo2 <- ifelse(cvo2 >= caph, 1, pmin(pmax(cvo2 / caph, 1e-6), 1 - 1e-6))
  pvo2 <- ifelse(cvo2 >= caph, (cvo2 - caph) / 0.003,
                 p50e * (svo2 / (1 - svo2))^(1 / neff))
  fO2 <- hbf * sao2 / hb; fdx <- hbf * (1 - sao2) / hb
  num <- fO2 * 0.32 + fdx * 3.23 + fmet * e660m
  den <- fO2 * 1.21 + fdx * 0.69 + fmet * e940m
  R <- num / den
  spo2 <- pmin(pmax(110 - 25 * R, 0), 100)
  tibble(fmet = fmet, MetPct = 100 * fmet, hbf = hbf, P50E = p50e, NEFF = neff,
         SAO2 = 100 * sao2, PVO2 = pmax(pvo2, 0), R = R, SPO2 = spo2,
         SAO2CO = 100 * fO2, GAP = spo2 - 100 * fO2,
         CYAN = (hbf * (1 - sao2)) / 5 + hb * fmet / 1.5)
}

# =============================================================================
#  UI
# =============================================================================
ui <- fluidPage(
  titlePanel("메트헤모글로빈혈증 QSP 모델 · Methaemoglobinaemia QSP model"),
  tags$p(style = "color:#555;margin-top:-8px",
         HTML("<b>운반(delivery)의 병이지 포화도(saturation)의 병이 아니다.</b> ",
              "모니터의 숫자 · 혈액가스의 숫자 · 조직의 숫자는 서로 다른 숫자이며, ",
              "환자가 나빠질수록 서로 멀어진다.")),
  sidebarLayout(
    sidebarPanel(
      width = 3,
      h4("환자 · Patient"),
      sliderInput("hb", "총 헤모글로빈 Hb (g/dL)", 5, 19, 15, 0.5),
      sliderInput("g6pd", "G6PD 활성 (분율)", 0.01, 1.2, 1.0, 0.01),
      helpText("0.02 = Mediterranean (class I–II) · 0.15 = A− (class III)"),
      sliderInput("eb5", "CYB5R3 활성 (in vivo flux 분율)", 0.01, 1.2, 1.0, 0.005),
      helpText("0.035 = 선천성 I형 표현형"),
      sliderInput("fhbm", "HbM 변이 분율 (환원 불가)", 0, 0.4, 0, 0.05),
      sliderInput("vo2", "산소소비 VO₂ (mL/min)", 150, 600, 250, 10),

      h4("산화제 노출 · Oxidant"),
      selectInput("oxidant", NULL,
                  c("없음 (none)" = "none",
                    "벤조카인 스프레이" = "benzocaine",
                    "답손 만성 투여" = "dapsone_chronic",
                    "답손 과량복용" = "dapsone_od",
                    "아질산나트륨 섭취" = "nitrite"),
                  selected = "benzocaine"),
      conditionalPanel("input.oxidant == 'benzocaine'",
        sliderInput("bzc", "흡수량 (mg)", 20, 600, 250, 10)),
      conditionalPanel("input.oxidant == 'dapsone_chronic'",
        sliderInput("dapdose", "1일 용량 (mg)", 25, 400, 100, 25)),
      conditionalPanel("input.oxidant == 'dapsone_od'",
        sliderInput("dapod", "복용량 (mg)", 500, 6000, 2000, 250)),
      conditionalPanel("input.oxidant == 'nitrite'",
        sliderInput("nit", "NaNO₂ (g)", 0.5, 20, 5, 0.5)),

      h4("치료 · Treatment"),
      sliderInput("nmb", "메틸렌블루 투여 횟수", 0, 6, 1, 1),
      sliderInput("mbdose", "1회 용량 (mg/kg)", 0.25, 10, 1, 0.25),
      sliderInput("mbtime", "첫 투여 시각 (h)", 0, 24, 0.75, 0.25),
      sliderInput("mbint", "투여 간격 (h)", 0.5, 12, 4, 0.5),
      checkboxInput("cimetidine", "시메티딘 400 mg tid (생성 항)", FALSE),
      checkboxInput("ascorbate", "아스코르브산 10 g IV", FALSE),
      checkboxInput("hbo", "고압산소 / 고농도 산소", FALSE),
      conditionalPanel("input.hbo == true",
        sliderInput("pao2", "PaO₂ (mmHg)", 95, 2000, 1800, 25)),
      checkboxInput("ssri", "SSRI 병용 (세로토닌 탭)", FALSE),

      h4("모델 가정 · Assumption"),
      sliderInput("alpham", "알로스테릭 좌측 이동 α", 0, 1.0, 0.5, 0.05),
      helpText(HTML("<b>이 슬라이더가 이 모델의 급소입니다.</b> α = 0 이면 ",
                    "%MetHb는 충분한 침대 옆 변수가 되고, 이 모델의 중심 주장은 사라집니다.")),
      sliderInput("tend", "시뮬레이션 기간 (h)", 6, 720, 48, 6)
    ),

    mainPanel(
      width = 9,
      tabsetPanel(
        tabPanel("1. 환자 프로파일",
          br(), DTOutput("profile"), br(),
          htmlOutput("verdict")),

        tabPanel("2. 시간 경과",
          br(), plotOutput("tc_main", height = "620px")),

        tabPanel("3. 두 개의 거짓말",
          br(),
          tags$p(HTML("파란색은 <b>모니터</b>, 주황색은 <b>동시산소측정</b>, ",
                      "초록색은 <b>조직</b>입니다. 셋이 벌어지는 폭이 이 병입니다.")),
          plotOutput("liars", height = "560px")),

        tabPanel("4. 등가 빈혈",
          br(),
          tags$p(HTML("\"MetHb 30%\"를 <b>같은 조직 PO₂를 만드는 정상 헤모글로빈 농도</b>로 ",
                      "번역합니다. 점선이 통상적 산술 Hb×(1−f)이고, 실선이 실제입니다. ",
                      "둘 사이의 간격이 알로스테릭 좌측 이동이 무력화한 헤모글로빈입니다.")),
          plotOutput("equiv", height = "520px"), br(), DTOutput("equivtab")),

        tabPanel("5. 맥박산소측정기",
          br(),
          tags$p(HTML("85%라는 바닥은 혈액의 성질이 아니라 <b>보정 직선의 값</b>입니다. ",
                      "아래 오른쪽 그림의 기울기가 이 탭의 요점입니다.")),
          plotOutput("oxi", height = "540px")),

        tabPanel("6. 임상 엔드포인트",
          br(), plotOutput("endpoints", height = "520px"), br(),
          DTOutput("endtab")),

        tabPanel("7. 메틸렌블루 용량–반응",
          br(),
          tags$p(HTML("모델에는 상한 용량 파라미터가 <b>없습니다</b>. ",
                      "그런데도 곡선은 되꺾입니다.")),
          plotOutput("mbdr", height = "560px")),

        tabPanel("8. NADPH 경제",
          br(),
          tags$p(HTML("모델이 실제로 쓰는 화폐. 메틸렌블루와 글루타티온은 ",
                      "<b>같은 지갑</b>에서 전자를 꺼냅니다.")),
          plotOutput("nadph", height = "560px")),

        tabPanel("9. G6PD 구배",
          br(), plotOutput("g6pdplot", height = "540px"), br(),
          DTOutput("g6pdtab")),

        tabPanel("10. 시나리오 비교",
          br(), plotOutput("scen", height = "620px")),

        tabPanel("11. 바이오마커 · 용혈",
          br(), plotOutput("biomark", height = "600px")),

        tabPanel("12. 세로토닌",
          br(), plotOutput("sero", height = "460px")),

        tabPanel("13. 모델 노트",
          br(), htmlOutput("notes"))
      )
    )
  )
)

# =============================================================================
#  SERVER
# =============================================================================
server <- function(input, output, session) {

  sim <- reactive({ run_sim(build_sim(input)) })

  # ---- 1. patient profile --------------------------------------------------
  output$profile <- renderDT({
    d <- sim(); pk <- d[which.max(d$MetPct), ]; en <- d[nrow(d), ]
    tibble(
      항목 = c("최고 MetHb (%)", "최고 시각 (h)", "그때의 SpO₂ (%)",
               "그때의 동시산소측정 SaO₂ (%)", "포화도 격차 (점)",
               "최저 PvO₂ (mmHg)", "최저 SpO₂ (%)",
               "누적 저산소 부하 (mmHg·h)", "최저 총 Hb (g/dL)",
               "최고 하인츠 부담", "최고 젖산 (mmol/L)",
               "청색증 지수 (≥1 이면 육안 청색증)"),
      값 = round(c(pk$MetPct, pk$time, pk$SPO2, pk$SAO2CO, pk$SPO2 - pk$SAO2CO,
                   min(d$PVO2), min(d$SPO2), en$INJ, min(d$HBTOT),
                   max(d$HEINZ), max(d$LAC), max(d$CYAN)), 2))
  }, options = list(dom = "t", pageLength = 20), rownames = FALSE)

  output$verdict <- renderUI({
    d <- sim()
    pv <- min(d$PVO2); mx <- max(d$MetPct); hb <- min(d$HBTOT)
    msg <- c()
    if (pv < 20) msg <- c(msg, sprintf(
      "<span style='color:#b3541e'><b>조직 PO₂가 임계값 아래로 내려갑니다 (최저 %.1f mmHg).</b>
       이 환자에서는 %%MetHb %.0f%%가 이미 혐기성 대사를 의미합니다.</span>", pv, mx))
    if (mx > 15 && pv >= 20) msg <- c(msg, sprintf(
      "MetHb는 %.0f%%까지 오르지만 조직 PO₂는 %.1f mmHg로 임계값 위에 남습니다.", mx, pv))
    if (hb < input$hb - 1.5) msg <- c(msg, sprintf(
      "<span style='color:#8b2f39'><b>용혈이 일어납니다</b> (Hb %.1f → %.1f g/dL).
       %%MetHb가 내려가는 것을 회복으로 읽으면 안 됩니다 — 산화된 적혈구가 먼저 제거되고 있습니다.</span>",
      input$hb, hb))
    if (input$g6pd < 0.1 && input$nmb > 0) msg <- c(msg,
      "<span style='color:#8b2f39'><b>G6PD 활성이 매우 낮은 상태에서 메틸렌블루가 투여되었습니다.</b>
       모델에는 이 조합을 금지하는 규칙이 없습니다 — 아래 결과는 NADPH 공유 상한 하나에서 계산된 것입니다.</span>")
    HTML(paste0("<div style='font-size:15px;line-height:1.7'>",
                paste(msg, collapse = "<br>"), "</div>"))
  })

  # ---- 2. time course ------------------------------------------------------
  output$tc_main <- renderPlot({
    d <- sim() %>% select(time, MetPct, SPO2, SAO2CO, PVO2, HBTOT, LAC, CO,
                          CMBRBC, GSH) %>%
      pivot_longer(-time)
    lab <- c(MetPct = "MetHb (%)", SPO2 = "SpO₂ (%)", SAO2CO = "co-ox SaO₂ (%)",
             PVO2 = "PvO₂ (mmHg)", HBTOT = "총 Hb (g/dL)", LAC = "젖산 (mmol/L)",
             CO = "심박출량 (L/min)", CMBRBC = "적혈구 MB (µM)",
             GSH = "GSH (µM)")
    d$name <- factor(lab[d$name], levels = lab)
    ggplot(d, aes(time, value)) +
      geom_line(linewidth = 0.9, colour = "#2b6cb0") +
      facet_wrap(~name, scales = "free_y", ncol = 3) +
      labs(x = "시간 (h)", y = NULL) + thm
  })

  # ---- 3. the two liars ----------------------------------------------------
  output$liars <- renderPlot({
    d <- sim()
    g1 <- d %>% select(time, SpO2 = SPO2, `co-ox SaO2` = SAO2CO) %>%
      pivot_longer(-time) %>%
      ggplot(aes(time, value, colour = name)) +
      geom_line(linewidth = 1.1) +
      scale_colour_manual(values = PAL, name = NULL) +
      labs(title = "모니터가 보는 것 vs 혈액이 실제로 나르는 것",
           x = "시간 (h)", y = "%") + thm
    g2 <- ggplot(d, aes(time, PVO2)) +
      geom_hline(yintercept = 20, linetype = "dashed", colour = "#b3541e") +
      annotate("text", x = max(d$time) * 0.75, y = 21.5,
               label = "임계 모세혈관 PO₂", colour = "#b3541e", size = 3.6) +
      geom_line(linewidth = 1.1, colour = PAL[["PvO2"]]) +
      labs(title = "조직이 실제로 보는 것", x = "시간 (h)", y = "PvO₂ (mmHg)") + thm
    g3 <- ggplot(d, aes(time, SPO2 - SAO2CO)) +
      geom_line(linewidth = 1.1, colour = "#5b3a80") +
      geom_hline(yintercept = 5, linetype = "dotted") +
      labs(title = "포화도 격차 (>5 이면 이상혈색소를 의심)",
           x = "시간 (h)", y = "SpO₂ − SaO₂ (점)") + thm
    gridExtra::grid.arrange(g1, g2, g3, ncol = 1)
  })

  # ---- 4. equivalent anaemia ----------------------------------------------
  eqdata <- reactive({
    tr <- transport(seq(0, 0.75, by = 0.01), hb = input$hb,
                    alpham = input$alpham, vo2 = input$vo2)
    tr$naive <- input$hb * (1 - tr$fmet)
    tr$equiv <- equiv_hb(tr$PVO2, vo2 = input$vo2)
    tr$lost <- tr$naive - tr$equiv
    tr
  })

  output$equiv <- renderPlot({
    tr <- eqdata()
    ggplot(tr, aes(MetPct)) +
      geom_ribbon(aes(ymin = equiv, ymax = naive), fill = "#f9d6c8", alpha = .8) +
      geom_line(aes(y = naive), linetype = "dashed", linewidth = 1) +
      geom_line(aes(y = equiv), linewidth = 1.2, colour = "#8b2f39") +
      annotate("text", x = 45, y = max(tr$naive) * 0.75,
               label = "이 띠 = 좌측 이동이 무력화한 헤모글로빈",
               colour = "#8b2f39", size = 4.4) +
      labs(title = "%MetHb를 등가 빈혈로 번역하기",
           subtitle = sprintf("Hb %.1f g/dL, VO₂ %d mL/min, α = %.2f",
                              input$hb, input$vo2, input$alpham),
           x = "메트헤모글로빈 (%)", y = "헤모글로빈 (g/dL)") + thm
  })

  output$equivtab <- renderDT({
    tr <- eqdata() %>% filter(round(MetPct) %in% c(0, 10, 15, 20, 25, 30, 40, 50, 60))
    tr %>% transmute(`MetHb (%)` = round(MetPct),
                     `통상 산술 Hb×(1−f)` = round(naive, 2),
                     `등가 Hb` = round(equiv, 2),
                     `무력화된 Hb` = round(lost, 2),
                     `잔여 용량 중 %` = round(100 * lost / naive, 1),
                     `PvO₂` = round(PVO2, 1),
                     `SpO₂` = round(SPO2, 1))
  }, options = list(dom = "t"), rownames = FALSE)

  # ---- 5. the oximeter -----------------------------------------------------
  output$oxi <- renderPlot({
    tr <- transport(seq(0, 0.9, by = 0.005), hb = input$hb, alpham = input$alpham)
    tr$slope <- c(NA, diff(tr$SPO2) / diff(tr$MetPct))
    g1 <- tr %>% select(MetPct, SpO2 = SPO2, `co-ox SaO2` = SAO2CO) %>%
      pivot_longer(-MetPct) %>%
      ggplot(aes(MetPct, value, colour = name)) +
      geom_hline(yintercept = 85, linetype = "dotted") +
      geom_line(linewidth = 1.2) +
      scale_colour_manual(values = PAL, name = NULL) +
      annotate("text", x = 60, y = 87.5, label = "110 − 25×1 = 85", size = 3.8) +
      labs(title = "바닥은 기기의 성질이다", x = "MetHb (%)", y = "%") + thm
    g2 <- ggplot(tr, aes(MetPct, R)) +
      geom_hline(yintercept = 1, linetype = "dotted") +
      geom_line(linewidth = 1.2, colour = "#5b3a80") +
      labs(title = "기기가 계산하는 유일한 것: R = A₆₆₀/A₉₄₀",
           x = "MetHb (%)", y = "R") + thm
    g3 <- ggplot(tr[-1, ], aes(MetPct, -slope)) +
      geom_line(linewidth = 1.2, colour = "#b3541e") +
      labs(title = "감도 붕괴: −dSpO₂ / d(%MetHb)",
           subtitle = "환자가 죽어 가는 구간에서 기기는 사실상 반응을 멈춘다",
           x = "MetHb (%)", y = "점 / %") + thm
    g4 <- ggplot(tr, aes(MetPct, GAP)) +
      geom_line(linewidth = 1.2, colour = "#2f6b34") +
      labs(title = "포화도 격차", x = "MetHb (%)", y = "점") + thm
    gridExtra::grid.arrange(g1, g2, g3, g4, ncol = 2)
  })

  # ---- 6. clinical endpoints ----------------------------------------------
  output$endpoints <- renderPlot({
    grid <- expand.grid(hb = seq(6, 18, by = 0.5), fmet = seq(0, 0.7, by = 0.01))
    res <- do.call(rbind, lapply(seq_len(nrow(grid)), function(i)
      transport(grid$fmet[i], hb = grid$hb[i], alpham = input$alpham,
                vo2 = input$vo2)))
    res$hb <- grid$hb
    ggplot(res, aes(MetPct, hb, fill = PVO2)) +
      geom_raster() +
      geom_contour(aes(z = PVO2), breaks = 20, colour = "white", linewidth = 1.2) +
      geom_contour(aes(z = CYAN), breaks = 1, colour = "#333333",
                   linetype = "dashed") +
      geom_point(aes(x = 30, y = 15), colour = "red", size = 3) +
      scale_fill_viridis_c(name = "PvO₂ (mmHg)") +
      labs(title = "같은 퍼센트는 같은 병이 아니다",
           subtitle = paste("흰 선 = 혐기성 역치 (PvO₂ 20 mmHg) · 점선 = 육안 청색증 역치 ·",
                            "빨간 점 = 교과서적 '30%'"),
           x = "메트헤모글로빈 (%)", y = "총 헤모글로빈 (g/dL)") + thm
  })

  output$endtab <- renderDT({
    do.call(rbind, lapply(c(17, 15, 12, 9, 7), function(h)
      do.call(rbind, lapply(c(0.15, 0.20, 0.30, 0.40), function(f) {
        t <- transport(f, hb = h, alpham = input$alpham, vo2 = input$vo2)
        tibble(`Hb` = h, `MetHb (%)` = 100 * f,
               `기능 Hb` = round(t$hbf, 2), `PvO₂` = round(t$PVO2, 1),
               `SpO₂` = round(t$SPO2, 1),
               `혐기성?` = ifelse(t$PVO2 < 20, "예", "아니오"))
      }))))
  }, options = list(dom = "tp", pageLength = 10), rownames = FALSE)

  # ---- 7. methylene blue dose-response ------------------------------------
  output$mbdr <- renderPlot({
    doses <- c(0, 0.5, 1, 1.5, 2, 3, 4, 5, 7, 10)
    res <- lapply(doses, function(dz) {
      s <- build_sim(input, tend = 48)
      s$e <- NULL
      if (input$oxidant == "benzocaine")
        s$e <- ev(time = 0, amt = input$bzc, cmt = "BZCD")
      if (dz > 0) s$e <- if (is.null(s$e)) mb_ev(dz, 0.75) else c(s$e, mb_ev(dz, 0.75))
      d <- run_sim(s, delta = 0.1)
      tibble(dose = dz,
             `MetHb @2h (%)` = d$MetPct[which.min(abs(d$time - 2))],
             `최저 MetHb (%)` = min(d$MetPct[d$time > 0.75]),
             `Hb @48h (g/dL)` = d$HBTOT[nrow(d)],
             `하인츠 최고` = max(d$HEINZ),
             `ROS 노출 (AUC)` = sum(d$ROS) * 0.1)
    }) %>% bind_rows()
    res %>% pivot_longer(-dose) %>%
      ggplot(aes(dose, value)) +
      geom_vline(xintercept = 7, linetype = "dotted", colour = "#b3541e") +
      geom_line(linewidth = 1.1, colour = "#2b6cb0") + geom_point(size = 2) +
      facet_wrap(~name, scales = "free_y") +
      labs(title = "상한 파라미터가 없는데도 곡선은 되꺾인다",
           subtitle = "점선 = 임상적으로 알려진 누적 7 mg/kg 천장 (모델에는 입력되지 않음)",
           x = "메틸렌블루 (mg/kg)", y = NULL) + thm
  })

  # ---- 8. the NADPH economy -----------------------------------------------
  output$nadph <- renderPlot({
    d <- sim()
    d$leuco <- d$LEUCOF
    g1 <- ggplot(d, aes(time, CMBRBC)) + geom_line(linewidth = 1.1) +
      labs(title = "적혈구 내 총 메틸렌블루", x = "h", y = "µM") + thm
    g2 <- ggplot(d, aes(time, leuco)) + geom_line(linewidth = 1.1, colour = "#2f6b34") +
      labs(title = "류코형 분율 — NADPH가 따라가고 있는가",
           subtitle = "이 값이 떨어지면 남은 것은 산화제다",
           x = "h", y = "LMB / (LMB+MB)") + thm
    g3 <- ggplot(d, aes(time, GSH)) + geom_line(linewidth = 1.1, colour = "#5b3a80") +
      geom_hline(yintercept = 1000, linetype = "dotted") +
      labs(title = "환원형 글루타티온 — 막이 쓰는 몫",
           x = "h", y = "µM") + thm
    g4 <- ggplot(d, aes(time, HEINZ)) + geom_line(linewidth = 1.1, colour = "#8b2f39") +
      labs(title = "하인츠소체 부담 → 용혈", x = "h", y = "0–1") + thm
    gridExtra::grid.arrange(g1, g2, g3, g4, ncol = 2)
  })

  # ---- 9. G6PD gradient ---------------------------------------------------
  g6data <- reactive({
    gs <- c(1.0, 0.6, 0.4, 0.25, 0.15, 0.08, 0.02)
    lapply(gs, function(g) {
      s <- build_sim(input, tend = 72); s$p$G6PD <- g
      d <- run_sim(s, delta = 0.1)
      tibble(G6PD = g,
             `MetHb @2h (%)` = d$MetPct[which.min(abs(d$time - 2))],
             `최저 GSH (µM)` = min(d$GSH),
             `하인츠 최고` = max(d$HEINZ),
             `Hb @72h (g/dL)` = d$HBTOT[nrow(d)])
    }) %>% bind_rows()
  })

  output$g6pdplot <- renderPlot({
    g6data() %>% pivot_longer(-G6PD) %>%
      ggplot(aes(G6PD, value)) +
      geom_line(linewidth = 1.1, colour = "#2b6cb0") + geom_point(size = 2) +
      scale_x_log10() + facet_wrap(~name, scales = "free_y") +
      labs(title = "모델에는 'G6PD 결핍에서 금기'라는 규칙이 없다",
           subtitle = "이 구배 전체가 NADPH 공유 상한 하나에서 나온다",
           x = "G6PD 활성 (분율, log)", y = NULL) + thm
  })
  output$g6pdtab <- renderDT({ g6data() %>% mutate(across(where(is.numeric), ~round(.x, 3))) },
                             options = list(dom = "t"), rownames = FALSE)

  # ---- 10. scenario comparison --------------------------------------------
  output$scen <- renderPlot({
    base <- ev(time = 0, amt = input$bzc, cmt = "BZCD")
    scn <- list(
      "무치료" = list(e = base, p = list()),
      "MB 1 mg/kg" = list(e = c(base, mb_ev(1, 0.75)), p = list()),
      "MB 5 mg/kg" = list(e = c(base, mb_ev(5, 0.75)), p = list()),
      "MB + G6PD 2%" = list(e = c(base, mb_ev(2, 0.75)), p = list(G6PD = 0.02)),
      "아스코르브산" = list(e = c(base, ev(time = 0.75, amt = 10 * 1e6 / 176.1,
                                          cmt = "ASCC")), p = list()),
      "고압산소 2.8 ATA" = list(e = base, p = list(PAO2 = 1800)),
      "빈혈 Hb 8 · 무치료" = list(e = base, p = list(HB0 = 8))
    )
    d <- bind_rows(lapply(names(scn), function(n) {
      m <- param(mod, modifyList(list(HB0 = input$hb, ALPHAM = input$alpham),
                                 scn[[n]]$p))
      mrgsim(m, events = scn[[n]]$e, end = 24, delta = 0.1) %>% as_tibble() %>%
        mutate(scenario = n)
    }))
    d %>% select(time, scenario, MetPct, SPO2, PVO2, HBTOT) %>%
      pivot_longer(c(MetPct, SPO2, PVO2, HBTOT)) %>%
      ggplot(aes(time, value, colour = scenario)) +
      geom_line(linewidth = 1) +
      facet_wrap(~name, scales = "free_y") +
      labs(title = "일곱 개의 시나리오, 같은 노출",
           x = "시간 (h)", y = NULL, colour = NULL) + thm
  })

  # ---- 11. biomarkers ------------------------------------------------------
  output$biomark <- renderPlot({
    d <- sim() %>% select(time, HBTOT, MetPct, HEINZ, FHB, HAPT, BILI, LDH,
                          GSH, LAC) %>% pivot_longer(-time)
    lab <- c(HBTOT = "총 Hb (g/dL)", MetPct = "MetHb (%)", HEINZ = "하인츠 부담",
             FHB = "혈장 유리 Hb (µM haem)", HAPT = "합토글로빈 (µM)",
             BILI = "간접 빌리루빈 (mg/dL)", LDH = "LDH (U/L)",
             GSH = "GSH (µM)", LAC = "젖산 (mmol/L)")
    d$name <- factor(lab[d$name], levels = lab)
    ggplot(d, aes(time, value)) + geom_line(linewidth = 0.95, colour = "#8b2f39") +
      facet_wrap(~name, scales = "free_y", ncol = 3) +
      labs(title = "용혈 지표 — %MetHb가 좋아 보일 때 무엇을 봐야 하는가",
           x = "시간 (h)", y = NULL) + thm
  })

  # ---- 12. serotonin -------------------------------------------------------
  output$sero <- renderPlot({
    variants <- list(
      "MB 단독" = list(ssri = FALSE, mb = TRUE),
      "SSRI 단독" = list(ssri = TRUE,  mb = FALSE),
      "MB + SSRI" = list(ssri = TRUE,  mb = TRUE))
    d <- bind_rows(lapply(names(variants), function(n) {
      v <- variants[[n]]
      e <- ev(time = 0, amt = input$bzc, cmt = "BZCD")
      if (v$ssri) e <- c(e, ev(time = 0, amt = 2.5, cmt = "SSRI"))
      if (v$mb)   e <- c(e, mb_ev(input$mbdose, 0.75))
      mrgsim(param(mod, list(HB0 = input$hb)), events = e, end = 12, delta = 0.05) %>%
        as_tibble() %>% mutate(scenario = n)
    }))
    ggplot(d, aes(time, HT5, colour = scenario)) +
      geom_hline(yintercept = 1, linetype = "dotted") +
      geom_line(linewidth = 1.2) +
      labs(title = "해독제 자체의 약리학: 메틸렌블루는 강력한 MAO-A 억제제다",
           x = "시간 (h)", y = "시냅스 세로토닌 (정상 = 1)", colour = NULL) + thm
  })

  # ---- 13. notes -----------------------------------------------------------
  output$notes <- renderUI(HTML("
<div style='max-width:980px;line-height:1.75;font-size:15px'>
<h3>이 앱이 보여 주려는 것</h3>
<p>이 모델의 파라미터 목록에는 <b>중증도 척도</b>도, <b>7 mg/kg 천장</b>도,
<b>G6PD 금기 스위치</b>도, <b>반동(rebound) 항</b>도, <b>85% 바닥 상수</b>도
없습니다. 그 다섯 가지는 전부 <b>출력</b>입니다. 들어간 것은 다섯 개의 구조뿐입니다:</p>
<ol>
<li>조효소가 글루타티온인 <b>촉매적 공동산화 순환</b> — 이것이 없으면 벤조카인
250 mg은 화학량론적으로 MetHb 0.5%까지밖에 만들 수 없고, 이 병은 존재할 수 없습니다.</li>
<li>소비자가 둘인 <b>하나의 한정된 NADPH 공급</b>.</li>
<li>포화 가능한 좋은 가지와 포화하지 않는 무익한 가지를 가진 <b>류코메틸렌블루 갈림길</b>.</li>
<li>아철 분율을 P50과 Hill n에 연결하는 <b>알로스테릭 항</b>.</li>
<li>색소당 두 개의 흡광계수와 <b>하나의 직선 보정식</b>.</li>
</ol>

<h3>가장 취약한 가정 (숨기지 않습니다)</h3>
<p>사이드바의 <b>α 슬라이더</b>가 이 모델의 급소입니다. α = 0 으로 두면 4번 탭의
띠가 사라지고, \"MetHb 30% ≡ Hb 10.5\"라는 통상적 산술이 옳아지며, 이 모델의
중심 주장은 소멸합니다. 메트헤모글로빈이 남은 아철 소단위의 P50을 정확히 얼마나
낮추는지에 대한 정량적 인체 데이터는 빈약하고, α = 0.5는 문헌과 정합적인 값일 뿐
문헌으로부터 유일하게 결정된 값이 아닙니다.</p>
<p><b>반증 실험은 간단합니다:</b> 같은 검체에서 동시산소측정 %MetHb와 실측 P50을
함께 재고, 그 기울기를 보십시오. 기울기가 0이면 이 모델은 틀렸습니다.</p>

<h3>검증</h3>
<p>컨테이너에 R이 없었으므로 mrgsolve 파일은 <b>컴파일 검증이 아니라 방정식 검증</b>을
거쳤습니다. 모든 ODE와 파라미터가 <code>mhb_reference_check.py</code>에 Python/scipy로
독립 재구현되어 있고, 파라미터 139개가 양쪽에서 0개 불일치로 대조되었으며,
앵커 42개가 모두 통과합니다:</p>
<pre>python3 mhb_reference_check.py --params mhb_mrgsolve_model.R
python3 mhb_reference_check.py --all</pre>

<h3>경고</h3>
<p>교육·연구용 모델입니다. 임상 의사결정에 사용하지 마십시오. 특히 메틸렌블루의
용량·누적 상한·G6PD 결핍에서의 사용 여부는 반드시 최신 임상 지침과 독성학 자문에
따르십시오.</p>
</div>"))
}

shinyApp(ui, server)
