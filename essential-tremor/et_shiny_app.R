# =============================================================================
#  Essential Tremor QSP — Shiny dashboard
#  ---------------------------------------------------------------------------
#  10 tabs.  The organising idea of the UI is the same as the model's: the user
#  does not set "tremor severity", they set LOOP GAIN, and then watch amplitude,
#  frequency and every clinical scale fall out of it.
#
#    1  환자 프로파일   Patient profile — G0, a_O, w_P, effector shares
#    2  진동자 지도     Oscillator — the bifurcation diagram, live
#    3  약물 PK         Plasma/brain concentrations, receptor occupancy
#    4  이득 분해       Gain decomposition: which phi is doing the work
#    5  떨림 진폭/주파수 Amplitude and frequency, with the waveform
#    6  임상 척도       TETRAS / FTM / spiral / QUEST, and the log artefact
#    7  시나리오 비교   Side-by-side regimens
#    8  수술·신경조절   Lesion volume and DBS frequency sweeps
#    9  감별진단        Weight-loading test: ET vs enhanced physiological tremor
#   10  장기계·이상반응 HR, FEV1, sedation, ataxia, grip, HCO3, BMD
#
#  Run:  shiny::runApp("et_shiny_app.R")
#  Requires: shiny, mrgsolve, dplyr, tidyr, ggplot2, DT
#  The model code is sourced from et_mrgsolve_model.R (its scenario block is
#  skipped by the ET_SHINY guard at the bottom of that file's helpers).
# =============================================================================

library(shiny); library(mrgsolve); library(dplyr); library(tidyr)
library(ggplot2); library(DT)

# ---- load the model without executing the scenario block --------------------
src   <- readLines("et_mrgsolve_model.R")
cut   <- grep("^#  SCENARIOS", src)
if (length(cut)) src <- src[1:(cut[1] - 3)]
eval(parse(text = paste(src, collapse = "\n")))

W  <- 24 * 7
TH <- theme_minimal(base_size = 12) +
  theme(panel.grid.minor = element_blank(),
        plot.title = element_text(face = "bold", size = 13),
        legend.position = "bottom")
PAL <- c("#2c5aa0", "#c1121f", "#2d6a4f", "#e07b39", "#7b2cbf",
         "#0b7285", "#b08900", "#6a4c93")

# =============================================================================
#  UI
# =============================================================================
ui <- fluidPage(
  titlePanel("본태성 떨림 QSP 모델 — Essential Tremor (oscillator formulation)"),
  tags$p(style = "color:#555;margin-top:-8px",
         HTML("떨림은 <b>양(level)</b>이 아니라 <b>한계순환(limit cycle)</b>입니다. ",
              "진폭은 루프 이득 G, 주파수는 루프 지연 &tau;, 치료의 천장은 루프 위상(topology)이 정합니다.")),
  sidebarLayout(
    sidebarPanel(
      width = 3,
      h4("환자 (Patient)"),
      sliderInput("G0", "루프 이득 G0 (1.15 경증 · 1.6 중등도 · 6 중증)",
                  1.02, 12, 1.60, 0.02),
      sliderInput("AO", "올리브 분율 a_O  (T형 차단제의 천장)", 0.05, 1.0, 0.35, 0.05),
      sliderInput("WP", "말초 루프 비중 w_P (β₂ 차단의 상한)", 0.0, 0.8, 0.40, 0.05),
      sliderInput("AGE", "연령 (주파수만 바꿈)", 30, 90, 60, 1),
      sliderInput("HDG", "경부 효과기 이득 분율 HDG", 0.30, 1.20, 0.55, 0.05),
      sliderInput("STRESS", "아드레날린 구동 배수 (스트레스·카페인)", 0, 9, 0, 0.5),
      checkboxInput("ASTHMA", "천식 동반 (β₂ 금기 계산)", FALSE),
      hr(),
      h4("약물 (Therapy)"),
      selectInput("beta", "β 차단제", c("없음", "프로프라놀롤", "아테놀롤", "나돌롤")),
      conditionalPanel("input.beta == '프로프라놀롤'",
                       sliderInput("dprp", "프로프라놀롤 LA mg/일", 0, 320, 160, 20)),
      sliderInput("dprm", "프리미돈 mg/일 (q8h)", 0, 750, 0, 25),
      sliderInput("dtop", "토피라메이트 mg/일 (q12h)", 0, 400, 0, 25),
      sliderInput("dttb", "T형 칼슘차단제 mg/일", 0, 400, 0, 25),
      sliderInput("drinks", "에탄올 표준잔 (단회)", 0, 4, 0, 1),
      hr(),
      h4("수술·기기 (Device)"),
      sliderInput("vles", "Vim 병변 부피 mm³ (0 = 없음)", 0, 400, 0, 10),
      checkboxInput("dbs", "Vim DBS 켜기", FALSE),
      conditionalPanel("input.dbs == true",
                       sliderInput("fstim", "자극 주파수 Hz", 10, 250, 130, 5),
                       sliderInput("vta", "활성 조직 부피 mm³", 30, 400, 250, 10)),
      sliderInput("btx", "보툴리눔 U (전완)", 0, 150, 0, 10),
      sliderInput("fspill", "확산 분율 f_spill (0.15 유도 · 0.45 비유도)",
                  0.05, 0.60, 0.15, 0.05),
      hr(),
      sliderInput("mload", "손목 부하 질량 kg (감별검사)", 0, 1.0, 0, 0.1),
      sliderInput("weeks", "시뮬레이션 기간 (주)", 1, 104, 24, 1),
      actionButton("go", "시뮬레이션 실행", class = "btn-primary btn-block")
    ),
    mainPanel(
      width = 9,
      tabsetPanel(
        id = "tabs",
        tabPanel("① 환자 프로파일", br(),
                 fluidRow(column(6, tableOutput("prof")),
                          column(6, plotOutput("profplot", height = 300))),
                 hr(), htmlOutput("profnote")),
        tabPanel("② 진동자 지도", br(),
                 plotOutput("bif", height = 380),
                 htmlOutput("bifnote"), hr(),
                 plotOutput("sqrtplot", height = 280)),
        tabPanel("③ 약물 PK", br(),
                 plotOutput("pk", height = 340), hr(),
                 plotOutput("occ", height = 280)),
        tabPanel("④ 이득 분해", br(),
                 plotOutput("phibar", height = 340),
                 htmlOutput("phinote"), hr(),
                 plotOutput("phitime", height = 280)),
        tabPanel("⑤ 진폭 · 주파수", br(),
                 plotOutput("amp", height = 320), hr(),
                 fluidRow(column(6, plotOutput("freq", height = 280)),
                          column(6, plotOutput("wave", height = 280)))),
        tabPanel("⑥ 임상 척도", br(),
                 plotOutput("scales", height = 320), hr(),
                 htmlOutput("lognote"), plotOutput("logmap", height = 300)),
        tabPanel("⑦ 시나리오 비교", br(),
                 DTOutput("cmp"), hr(), plotOutput("cmpplot", height = 340)),
        tabPanel("⑧ 수술 · 신경조절", br(),
                 fluidRow(column(6, plotOutput("lesion", height = 320)),
                          column(6, plotOutput("dbsf", height = 320))),
                 hr(), htmlOutput("surgnote"), plotOutput("habit", height = 280)),
        tabPanel("⑨ 감별진단", br(),
                 plotOutput("loadtest", height = 340),
                 htmlOutput("dxnote"), hr(), DTOutput("dxtab")),
        tabPanel("⑩ 장기계 · 이상반응", br(),
                 plotOutput("organs", height = 380), hr(),
                 fluidRow(column(6, plotOutput("btxwin", height = 300)),
                          column(6, tableOutput("adrtab"))))
      )
    )
  )
)

# =============================================================================
#  SERVER
# =============================================================================
server <- function(input, output, session) {

  pars <- reactive({
    p <- list(G0 = input$G0, AO = input$AO, AR = 1 - input$AO,
              WP = input$WP, WC = 1 - input$WP, AGE = input$AGE,
              HDG = input$HDG, STRESS = input$STRESS,
              ASTHMA = as.numeric(input$ASTHMA), MLOAD = input$mload,
              VLES = input$vles, LESION = as.numeric(input$vles > 0),
              DBSON = as.numeric(input$dbs), FSPILL = input$fspill)
    if (input$dbs) { p$FSTIM <- input$fstim; p$VTA <- input$vta }
    if (input$beta == "프로프라놀롤") p$DPRP <- input$dprp
    p
  })

  evs <- reactive({
    e <- NULL; wk <- input$weeks
    add <- function(x) if (is.null(e)) x else c(e, x)
    if (input$beta == "프로프라놀롤" && input$dprp > 0)
      e <- add(ev(amt = input$dprp, cmt = "A_PRPG", ii = 24, addl = wk*7 - 1))
    if (input$beta == "아테놀롤")
      e <- add(ev(amt = 100, cmt = "A_ATNG", ii = 24, addl = wk*7 - 1))
    if (input$beta == "나돌롤")
      e <- add(ev(amt = 120, cmt = "A_NADG", ii = 24, addl = wk*7 - 1))
    if (input$dprm > 0)
      e <- add(ev(amt = input$dprm/3, cmt = "A_PRMG", ii = 8, addl = wk*21 - 1))
    if (input$dtop > 0)
      e <- add(ev(amt = input$dtop/2, cmt = "A_TOPG", ii = 12, addl = wk*14 - 1))
    if (input$dttb > 0)
      e <- add(ev(amt = input$dttb, cmt = "A_TTBG", ii = 24, addl = wk*7 - 1))
    if (input$drinks > 0)
      e <- add(ev(amt = 14*input$drinks, cmt = "A_ETHG", time = 24))
    if (input$btx > 0) {
      e <- add(ev(amt = input$btx*(1 - input$fspill), cmt = "A_BTXT"))
      e <- add(ev(amt = input$btx*input$fspill,       cmt = "A_BTXG"))
    }
    e
  })

  sim <- eventReactive(input$go, {
    m <- do.call(param, c(list(mod), pars()))
    dlt <- if (input$weeks <= 4) 0.25 else 2
    mrgsim(m, events = evs(), end = input$weeks*W, delta = dlt) %>% as.data.frame()
  }, ignoreNULL = FALSE)

  base_sim <- reactive({
    m <- do.call(param, c(list(mod), pars()))
    mrgsim(m, end = 8*W, delta = 4) %>% as.data.frame() %>% tail(1)
  })
  fin <- reactive({ d <- sim(); tail(d, min(25, nrow(d))) %>%
      summarise(across(where(is.numeric), mean)) })

  # ---- ① profile -----------------------------------------------------------
  output$prof <- renderTable({
    b <- base_sim()
    data.frame(항목 = c("루프 이득 G", "분기 파라미터 μ = G−1", "한계순환 진폭 r*",
                        "상지 진폭 A (cm)", "생리적 떨림 바닥 (cm)",
                        "중추 진동 주파수 (Hz)", "사지 기계 공명 f₀ (Hz)",
                        "관측 우세 주파수 (Hz)", "TETRAS 수행 (0-64)",
                        "경부 μ_HD (음수 = 두부 떨림 없음)"),
               값 = sprintf("%.3f", c(b$GTOT, b$MU, sqrt(pmax(b$MU, 0)),
                                      b$A_UL, b$A_PHYS, b$FNEUR, b$F0M,
                                      b$F_OBS, b$TETRAS_PS, b$MU_HD)))
  }, striped = TRUE)

  output$profplot <- renderPlot({
    g <- seq(1.0, max(3, input$G0*1.3), length.out = 300)
    data.frame(G = g, A = ifelse(g > 1, sqrt(g - 1), 0)) %>%
      ggplot(aes(G, A)) + geom_line(linewidth = 1.2, colour = PAL[1]) +
      geom_vline(xintercept = 1, linetype = 2, colour = PAL[2]) +
      geom_point(data = data.frame(G = input$G0, A = sqrt(max(input$G0 - 1, 0))),
                 size = 4, colour = PAL[2]) +
      labs(title = "r* = √(G−1): 이 환자가 곡선의 어디에 있는가",
           x = "루프 이득 G", y = "한계순환 진폭 r*") + TH
  })
  output$profnote <- renderUI(HTML(
    "<div style='background:#eef3ff;padding:10px;border-left:4px solid #2c5aa0'>
     <b>G0</b>가 이 모델의 유일한 중증도 축입니다. 곡선이 G=1에서 <i>수직 접선</i>을 갖기
     때문에, G=1 근처의 환자는 같은 약으로 떨림이 <b>사라지고</b> 오른쪽 끝의 환자는
     <b>거의 변하지 않습니다</b>. 시험에서 반복 관찰되는 반응자/비반응자 분리는 두 종류의
     생물학이 아니라 하나의 방정식이 문턱을 넘느냐의 문제입니다.</div>"))

  # ---- ② bifurcation ------------------------------------------------------
  output$bif <- renderPlot({
    b <- base_sim(); f <- fin()
    g <- seq(0.85, max(3.5, input$G0*1.25), length.out = 400)
    df <- data.frame(G = g, r = ifelse(g > 1, sqrt(g - 1), 0))
    ggplot(df, aes(G, r)) +
      geom_line(linewidth = 1.2, colour = PAL[1]) +
      geom_hline(yintercept = 0, colour = "grey70") +
      geom_vline(xintercept = 1, linetype = 2, colour = PAL[2]) +
      annotate("text", x = 0.93, y = max(df$r)*0.9, label = "μ<0\n한계순환 없음",
               colour = PAL[2], size = 3.5) +
      geom_segment(aes(x = b$GTOT, xend = f$GTOT,
                       y = sqrt(max(b$MU, 0)), yend = sqrt(max(f$MU, 0))),
                   arrow = arrow(length = unit(0.22, "cm")), linewidth = 1,
                   colour = PAL[3]) +
      geom_point(data = data.frame(G = c(b$GTOT, f$GTOT),
                                   r = sqrt(pmax(c(b$MU, f$MU), 0)),
                                   w = c("치료 전", "치료 후")),
                 aes(G, r, colour = w), size = 4) +
      scale_colour_manual(values = c("치료 전" = PAL[2], "치료 후" = PAL[3]), name = NULL) +
      labs(title = "초임계 Hopf 분기: 치료는 환자를 곡선 위에서 왼쪽으로 밀 뿐이다",
           x = "루프 이득 G", y = "진폭 포락선 r*") + TH
  })
  output$bifnote <- renderUI({
    b <- base_sim(); f <- fin()
    HTML(sprintf("<div style='background:#f7f9fc;padding:8px'>G %.3f → %.3f
      (Δ %.1f%%) · μ %.3f → %.3f · 진폭 %.3f → %.3f cm (<b>%+.1f%%</b>) ·
      TETRAS %.1f → %.1f (<b>%+.2f점</b>)</div>",
      b$GTOT, f$GTOT, 100*(f$GTOT-b$GTOT)/b$GTOT, b$MU, f$MU,
      b$A_UL, f$A_UL, 100*(f$A_UL-b$A_UL)/b$A_UL,
      b$TETRAS_PS, f$TETRAS_PS, f$TETRAS_PS-b$TETRAS_PS))
  })
  output$sqrtplot <- renderPlot({
    # what the SAME gain reduction does across severities
    fr <- (fin()$GTOT/base_sim()$GTOT)
    g0 <- c(1.05, 1.15, 1.3, 1.6, 2.4, 4, 6, 9, 12)
    data.frame(G0 = g0,
               dA = 100*(sqrt(pmax(g0*fr - 1, 0))/sqrt(pmax(g0 - 1, 1e-9)) - 1)) %>%
      ggplot(aes(factor(G0), dA)) +
      geom_col(fill = PAL[1]) + geom_hline(yintercept = -50, linetype = 2, colour = PAL[2]) +
      labs(title = sprintf("현재 처방이 만드는 이득비 %.3f 를 여러 중증도에 적용", fr),
           x = "환자의 G0", y = "진폭 변화 (%)") + TH
  })

  # ---- ③ PK ---------------------------------------------------------------
  output$pk <- renderPlot({
    sim() %>% select(time, C_PRP, C_PRM, C_PB, C_TOP, C_ETH) %>%
      pivot_longer(-time) %>% filter(value > 1e-9) %>%
      ggplot(aes(time/24, value, colour = name)) + geom_line(linewidth = 0.9) +
      scale_colour_manual(values = PAL, name = NULL) +
      facet_wrap(~name, scales = "free_y") +
      labs(title = "혈장 농도 (mg/L; 에탄올은 g/L)", x = "일 (day)", y = NULL) + TH
  })
  output$occ <- renderPlot({
    sim() %>% select(time, OCCB1, OCCB2, BLK_TT, P_EFF, REB) %>%
      pivot_longer(-time) %>%
      ggplot(aes(time/24, value, colour = name)) + geom_line(linewidth = 0.9) +
      scale_colour_manual(values = PAL, name = NULL) +
      labs(title = "수용체 점유·차단·GABA-A 증강·반동 인자",
           subtitle = "OCCB2 하나가 떨림 효능과 기관지수축을 동시에 만든다",
           x = "일 (day)", y = "분율") + TH
  })

  # ---- ④ gain decomposition ----------------------------------------------
  output$phibar <- renderPlot({
    b <- base_sim(); f <- fin()
    data.frame(term = rep(c("φ_olive", "φ_cblthal", "φ_thal", "φ_ctx",
                            "φ_spindle", "φ_nmj"), 2),
               when = rep(c("치료 전", "치료 후"), each = 6),
               v = c(b$PHI_OL, b$PHI_CBL, b$PHI_TH, 1, b$PHI_SPIN, b$PHI_NMJ,
                     f$PHI_OL, f$PHI_CBL, f$PHI_TH, 1, f$PHI_SPIN, f$PHI_NMJ)) %>%
      ggplot(aes(term, v, fill = when)) +
      geom_col(position = "dodge") + geom_hline(yintercept = 1, linetype = 2) +
      scale_fill_manual(values = c(PAL[2], PAL[3]), name = NULL) +
      labs(title = "어느 항이 실제로 일을 하고 있는가", x = NULL, y = "인자 (1 = 무변화)") + TH
  })
  output$phinote <- renderUI(HTML(
    "<div style='background:#fff2f2;padding:10px;border-left:4px solid #c1121f'>
     <b>φ_olive</b>는 <i>병렬</i> 가지이므로 0까지 내려도 Φ_C는 a_R 아래로 못 갑니다
     (T형 차단제의 천장). <b>φ_thal</b>은 <i>직렬</i> 인자이므로 중추 루프 전체를 곱합니다
     (수술이 약을 이기는 이유는 효력이 아니라 위상입니다).</div>"))
  output$phitime <- renderPlot({
    sim() %>% select(time, PHI_OL, PHI_CBL, PHI_TH, PHI_SPIN, PHI_NMJ, GTOT) %>%
      pivot_longer(-time) %>%
      ggplot(aes(time/24, value, colour = name)) + geom_line(linewidth = 0.9) +
      scale_colour_manual(values = PAL, name = NULL) +
      labs(title = "이득 항의 시간 경과", x = "일 (day)", y = NULL) + TH
  })

  # ---- ⑤ amplitude / frequency / waveform --------------------------------
  output$amp <- renderPlot({
    sim() %>% select(time, A_UL, A_LC, A_PHYS, A_HD) %>% pivot_longer(-time) %>%
      ggplot(aes(time/24, value, colour = name)) + geom_line(linewidth = 0.9) +
      scale_colour_manual(values = PAL, name = NULL,
                          labels = c("총 관측 진폭 (cm)", "한계순환 성분 (cm)",
                                     "생리적 바닥 (cm)", "두부 (deg)")) +
      labs(title = "진폭: 한계순환 + 생리적 바닥 (직교 합)",
           x = "일 (day)", y = NULL) + TH
  })
  output$freq <- renderPlot({
    sim() %>% select(time, FNEUR, F0M, F_OBS) %>% pivot_longer(-time) %>%
      ggplot(aes(time/24, value, colour = name)) + geom_line(linewidth = 1) +
      scale_colour_manual(values = PAL, name = NULL) +
      labs(title = "주파수: 약물은 여기 나타나지 않는다",
           subtitle = "f는 루프 지연이 정하고, 지연은 연령·Purkinje 손실이 정한다",
           x = "일 (day)", y = "Hz") + TH
  })
  output$wave <- renderPlot({
    f <- fin(); fr <- f$F_OBS; amp <- f$A_UL
    t <- seq(0, 3, by = 0.002)
    data.frame(t = t, x = amp*sin(2*pi*fr*t)) %>%
      ggplot(aes(t, x)) + geom_line(colour = PAL[1]) +
      labs(title = sprintf("파형 (%.2f Hz, %.2f cm)", fr, amp),
           x = "초 (s)", y = "cm") + TH
  })

  # ---- ⑥ scales ----------------------------------------------------------
  output$scales <- renderPlot({
    b <- base_sim(); f <- fin()
    data.frame(scale = rep(c("TETRAS-PS/64", "TETRAS-ADL/48", "FTM/144",
                             "나선/4", "Bain-Findley/10", "QUEST/100"), 2),
               when = rep(c("치료 전", "치료 후"), each = 6),
               v = c(b$TETRAS_PS/64, b$TETRAS_ADL/48, b$FTM/144, b$SPIRAL/4,
                     b$BAINF/10, b$QUEST/100,
                     f$TETRAS_PS/64, f$TETRAS_ADL/48, f$FTM/144, f$SPIRAL/4,
                     f$BAINF/10, f$QUEST/100)) %>%
      ggplot(aes(scale, 100*v, fill = when)) + geom_col(position = "dodge") +
      scale_fill_manual(values = c(PAL[2], PAL[3]), name = NULL) +
      labs(title = "임상 척도 (각 척도 최대값에 대한 %)", x = NULL, y = "%") + TH
  })
  output$lognote <- renderUI({
    b <- base_sim(); f <- fin()
    HTML(sprintf("<div style='background:#f2ecff;padding:10px;border-left:4px solid #5a189a'>
      가속도계 <b>%+.1f%%</b> 인데 TETRAS는 <b>%+.2f점</b>(기저 %.1f점의 %+.1f%%)입니다.
      이것은 불일치가 아니라 <b>로그</b>입니다. 등급 = 2 + 2·log₁₀(A) 이므로 1점을 얻으려면
      진폭이 <b>3.16배</b> 줄어야 합니다.</div>",
      100*(f$A_UL-b$A_UL)/b$A_UL, f$TETRAS_PS-b$TETRAS_PS, b$TETRAS_PS,
      100*(f$TETRAS_PS-b$TETRAS_PS)/b$TETRAS_PS))
  })
  output$logmap <- renderPlot({
    a <- 10^seq(-1.5, 1.4, length.out = 300)
    data.frame(A = a, T = pmin(4, pmax(0, 2 + 2*log10(a)))) %>%
      ggplot(aes(A, T)) + geom_line(linewidth = 1.2, colour = PAL[5]) +
      scale_x_log10() +
      geom_point(data = data.frame(A = c(base_sim()$A_UL, fin()$A_UL),
                                   T = c(base_sim()$T_R, fin()$T_R),
                                   w = c("전", "후")),
                 aes(colour = w), size = 4) +
      scale_colour_manual(values = c(PAL[2], PAL[3]), name = NULL) +
      labs(title = "Elble 로그 변환: 등급 = 2 + 2·log₁₀(진폭 cm)",
           x = "진폭 (cm, 로그축)", y = "항목 등급 (0-4)") + TH
  })

  # ---- ⑦ scenario comparison ---------------------------------------------
  scen <- reactive({
    base_p <- pars(); base_p$DPRP <- 0
    mk <- function(nm, pp, e) {
      m <- do.call(param, c(list(mod), pp))
      d <- mrgsim(m, events = e, end = 8*W, delta = 4) %>% as.data.frame() %>% tail(7) %>%
        summarise(across(where(is.numeric), mean))
      data.frame(요법 = nm, G = d$GTOT, `진폭_cm` = d$A_UL, TETRAS = d$TETRAS_PS,
                 QUEST = d$QUEST, HR = d$HR, 졸림 = d$SED, 실조 = d$ATAX,
                 악력 = d$GRIP, check.names = FALSE)
    }
    R <- function(a, c_, ii) ev(amt = a, cmt = c_, ii = ii, addl = ceiling(8*W/ii) - 1)
    bind_rows(
      mk("무치료", base_p, NULL),
      mk("프로프라놀롤 160", c(base_p, list(DPRP = 160)), R(160, "A_PRPG", 24)),
      mk("아테놀롤 100", base_p, R(100, "A_ATNG", 24)),
      mk("나돌롤 120", base_p, R(120, "A_NADG", 24)),
      mk("프리미돈 750", base_p, R(250, "A_PRMG", 8)),
      mk("병용 (프로프+프리미돈)", c(base_p, list(DPRP = 160)),
         c(R(160, "A_PRPG", 24), R(250, "A_PRMG", 8))),
      mk("토피라메이트 400", base_p, R(200, "A_TOPG", 12)),
      mk("T형 차단제 100", base_p, R(100, "A_TTBG", 24)),
      mk("MRgFUS 120 mm³", c(base_p, list(LESION = 1, VLES = 120)), NULL),
      mk("Vim DBS 130 Hz", c(base_p, list(DBSON = 1, FSTIM = 130)), NULL)
    ) %>% mutate(`진폭변화_%` = 100*(진폭_cm/진폭_cm[1] - 1),
                 `TETRAS변화` = TETRAS - TETRAS[1])
  })
  output$cmp <- renderDT(scen() %>% mutate(across(where(is.numeric), ~round(., 2))),
                         options = list(pageLength = 12, dom = "t"))
  output$cmpplot <- renderPlot({
    scen() %>% filter(요법 != "무치료") %>%
      ggplot(aes(reorder(요법, `진폭변화_%`), `진폭변화_%`)) +
      geom_col(fill = PAL[1]) + coord_flip() +
      geom_text(aes(label = sprintf("%+.1f%% / %+.1f점", `진폭변화_%`, `TETRAS변화`)),
                hjust = 1.05, colour = "white", size = 3.2) +
      labs(title = "진폭 변화 (%)와 같은 결과의 TETRAS 점수 변화",
           x = NULL, y = "%") + TH
  })

  # ---- ⑧ surgery ----------------------------------------------------------
  output$lesion <- renderPlot({
    v <- seq(10, 400, by = 15)
    do.call(rbind, lapply(v, function(x) {
      d <- mrgsim(do.call(param, c(list(mod), pars(), list(LESION = 1, VLES = x))),
                  end = 4*W, delta = 24) %>% as.data.frame() %>% tail(1)
      data.frame(V = x, dA = 100*(d$A_UL/base_sim()$A_UL - 1), ATAX = d$ATAX)
    })) %>% pivot_longer(-V) %>%
      ggplot(aes(V, abs(value), colour = name)) + geom_line(linewidth = 1.1) +
      scale_colour_manual(values = c(PAL[3], PAL[2]), name = NULL,
                          labels = c("실조 지수", "진폭 감소 |%|")) +
      labs(title = "병변 부피: 이득은 포화하고 실조는 포화하지 않는다",
           subtitle = "V50 효능 45 mm³ vs V50 실조 260 mm³ — 치료창은 유도된 것",
           x = "병변 부피 (mm³)", y = NULL) + TH
  })
  output$dbsf <- renderPlot({
    f <- c(seq(10, 100, 10), 130, 160, 185, 220, 250)
    do.call(rbind, lapply(f, function(x) {
      d <- mrgsim(do.call(param, c(list(mod), pars(),
                                   list(DBSON = 1, FSTIM = x, VTA = input$vta))),
                  end = 4*W, delta = 24) %>% as.data.frame() %>% tail(1)
      data.frame(f = x, dA = 100*(d$A_UL/base_sim()$A_UL - 1))
    })) %>% ggplot(aes(f, dA)) + geom_line(linewidth = 1.1, colour = PAL[1]) +
      geom_hline(yintercept = 0, linetype = 2, colour = PAL[2]) +
      geom_vline(xintercept = 100, linetype = 3) +
      labs(title = "DBS 주파수는 다이얼이 아니라 스위치",
           subtitle = "Hill 4, f₅₀ 80 Hz — 저주파에서는 동조로 오히려 악화",
           x = "자극 주파수 (Hz)", y = "진폭 변화 (%)") + TH
  })
  output$surgnote <- renderUI(HTML(
    "<div style='background:#eefaf1;padding:10px;border-left:4px solid #2d6a4f'>
     같은 병변 부피가 <b>치료 효과</b>와 <b>보행 실조</b>를 동시에 만듭니다. 효능은 45 mm³에서
     반포화되지만 실조는 260 mm³까지 계속 오르므로, 부피를 키워 얻는 추가 이득은 거의 없고
     추가 실조는 계속 붙습니다. 습관화는 <b>기간</b>이 아니라 <b>재배선 용량</b>의 문제이며,
     문턱을 넘는 환자에서만 (그리고 비교적 갑작스럽게) 떨림이 돌아옵니다.</div>"))
  output$habit <- renderPlot({
    do.call(rbind, lapply(c(0.45, 0.75, 0.95, 1.00), function(rm) {
      do.call(rbind, lapply(c(0.25, 1, 2, 3, 4, 5), function(yr) {
        d <- mrgsim(do.call(param, c(list(mod), pars(),
                    list(LESION = 1, VLES = 120, REROUTE_MAX = rm))),
                    end = yr*8760, delta = 168) %>% as.data.frame() %>% tail(1)
        data.frame(rm = factor(rm), yr = yr,
                   dA = 100*(d$A_UL/base_sim()$A_UL - 1))
      }))
    })) %>% ggplot(aes(yr, dA, colour = rm)) +
      geom_line(linewidth = 1) + geom_point() +
      scale_colour_manual(values = PAL, name = "재배선 용량") +
      labs(title = "MRgFUS 후 습관화: 문턱을 넘는 환자만 재발한다",
           x = "년 (year)", y = "진폭 변화 (%)") + TH
  })

  # ---- ⑨ differential diagnosis -----------------------------------------
  dxd <- reactive({
    mk <- function(nm, pp, ml) {
      d <- mrgsim(do.call(param, c(list(mod), pars(), pp, list(MLOAD = ml))),
                  end = 4*W, delta = 24) %>% as.data.frame() %>% tail(1)
      data.frame(조건 = nm, 부하_kg = ml, `f0_Hz` = d$F0M, `f중추_Hz` = d$FNEUR,
                 `f관측_Hz` = d$F_OBS, `A_cm` = d$A_UL,
                 한계순환 = d$A_LC, 생리적 = d$A_PHYS, check.names = FALSE)
    }
    bind_rows(
      mk("ET", list(), 0), mk("ET", list(), 0.5), mk("ET", list(), 1.0),
      mk("강화 생리적 떨림 (EPT)", list(G0 = 0.42, THYRO = 9), 0),
      mk("강화 생리적 떨림 (EPT)", list(G0 = 0.42, THYRO = 9), 0.5),
      mk("강화 생리적 떨림 (EPT)", list(G0 = 0.42, THYRO = 9), 1.0))
  })
  output$loadtest <- renderPlot({
    dxd() %>% pivot_longer(c(`f0_Hz`, `f관측_Hz`)) %>%
      ggplot(aes(부하_kg, value, colour = 조건, linetype = name)) +
      geom_line(linewidth = 1.1) + geom_point(size = 2.5) +
      scale_colour_manual(values = c(PAL[1], PAL[2]), name = NULL) +
      scale_linetype_manual(values = c(2, 1), name = NULL) +
      labs(title = "질량 부하 시험 — 모델이 진단 검사 자체를 재현한다",
           subtitle = "EPT의 피크는 기계적 공명이라 1/√J로 내려가고, ET의 피크는 중추 지연이라 거의 안 움직인다",
           x = "손목 부하 질량 (kg)", y = "주파수 (Hz)") + TH
  })
  output$dxnote <- renderUI(HTML(
    "<div style='background:#f0f3f6;padding:10px;border-left:4px solid #495867'>
     같은 방정식 하나로 두 질환이 나옵니다. <b>ET</b>는 μ&gt;0인 한계순환이고 주파수는
     1/τ_loop, <b>EPT</b>는 μ&lt;0이라 한계순환이 없고 잡음이 기계 공명을 통과한 것이므로
     주파수가 f₀입니다. 그래서 <b>추만 달아도</b> 둘이 갈립니다.</div>"))
  output$dxtab <- renderDT(dxd() %>% mutate(across(where(is.numeric), ~round(., 3))),
                           options = list(dom = "t"))

  # ---- ⑩ organs -----------------------------------------------------------
  output$organs <- renderPlot({
    sim() %>% select(time, HR, SBP, FEV1, SED, ATAX, GRIP, HCO3, COG) %>%
      pivot_longer(-time) %>%
      ggplot(aes(time/24, value, colour = name)) + geom_line(linewidth = 0.9) +
      facet_wrap(~name, scales = "free_y", ncol = 4) +
      scale_colour_manual(values = rep(PAL, 2), guide = "none") +
      labs(title = "장기계 지표", x = "일 (day)", y = NULL) + TH
  })
  output$btxwin <- renderPlot({
    do.call(rbind, lapply(c(0.05, 0.15, 0.25, 0.35, 0.45, 0.60), function(fs) {
      d <- mrgsim(do.call(param, c(list(mod), pars(), list(FSPILL = fs))),
                  events = c(ev(amt = max(input$btx, 100)*(1 - fs), cmt = "A_BTXT"),
                             ev(amt = max(input$btx, 100)*fs,       cmt = "A_BTXG")),
                  end = 6*W, delta = 24) %>% as.data.frame()
      k <- which.min(d$A_UL)
      data.frame(f_spill = fs, `떨림감소` = -100*(d$A_UL[k]/base_sim()$A_UL - 1),
                 `악력손실` = 100 - d$GRIP[k], check.names = FALSE)
    })) %>% pivot_longer(-f_spill) %>%
      ggplot(aes(f_spill, value, colour = name)) +
      geom_line(linewidth = 1.1) + geom_point() +
      scale_colour_manual(values = c(PAL[2], PAL[3]), name = NULL) +
      labs(title = "보툴리눔: 치료창은 용량이 아니라 정밀도가 만든다",
           subtitle = "용량을 줄이면 둘이 같이 줄고, f_spill을 줄이면 둘이 갈린다",
           x = "확산 분율 f_spill", y = "%") + TH
  })
  output$adrtab <- renderTable({
    f <- fin()
    data.frame(지표 = c("심박수 (bpm)", "수축기 혈압 (mmHg)", "FEV₁ (L)",
                        "졸림 (0-100)", "실조 (0-100)", "취기 (0-100)",
                        "악력 (%)", "HCO₃⁻ (mmol/L)", "체중 (kg)",
                        "어휘인출 장애 (0-100)", "골밀도 지수", "ALT (U/L)"),
               값 = sprintf("%.2f", c(f$HR, f$SBP, f$FEV1, f$SED, f$ATAX, f$INTOX,
                                      f$GRIP, f$HCO3, f$BW, f$COG, f$BMD, f$ALT)))
  }, striped = TRUE)
}

shinyApp(ui, server)
