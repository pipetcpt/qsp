## =============================================================================
##  HIGH-ALTITUDE ILLNESS QSP — Shiny dashboard
##  고산병 (AMS · HACE · HAPE) 인터랙티브 대시보드
## =============================================================================
##
##  실행:
##    source("hai_mrgsolve_model.R")   # 모델 + 시나리오 정의
##    shiny::runApp("hai_shiny_app.R")
##
##  12 tabs:
##    ① 환자 프로파일   ② 등반 프로파일   ③ 가스교환      ④ 환기·산-염기
##    ⑤ 수면·CO₂ 예비량 ⑥ 뇌 (AMS/HACE)  ⑦ 폐 (HAPE)     ⑧ 임상 엔드포인트
##    ⑨ 시나리오 비교   ⑩ 하산 등가       ⑪ 바이오마커    ⑫ 이분기·민감도
##
##  이 대시보드의 편집 방침: 모든 탭은 "무엇이 무엇을 움직이는가"를 하나씩
##  보여준다.  특히 ⑩번 탭(하산 등가)은 모든 처치를 단 하나의 통화 —
##  "몇 미터의 하산에 해당하는가" — 로 환산해서 나란히 놓는다.  덱사메타손이
##  0 m 로 나오는 것이 이 앱에서 가장 중요한 화면이다.
## =============================================================================

library(shiny)
library(mrgsolve)
library(dplyr)
library(tidyr)
library(ggplot2)

if (!exists("mod")) source("hai_mrgsolve_model.R")

THEME <- theme_bw(base_size = 12) +
  theme(strip.background = element_rect(fill = "#eceff1"),
        panel.grid.minor = element_blank(),
        legend.position = "bottom")

PAL <- c("#1565c0", "#c62828", "#2e7d32", "#ef6c00", "#6a1b9a", "#00838f",
         "#5d4037", "#455a64")

## =============================================================================
##  UI
## =============================================================================
ui <- fluidPage(
  titlePanel("고산병 QSP 모델 — AMS · HACE · HAPE (High-Altitude Illness)"),
  tags$p(style = "color:#555;margin-top:-8px",
         HTML("하나의 흡입 산소분압, 세 개의 방어기전, 세 개의 시간상수. ",
              "<b>P<sub>IO2</sub> = F<sub>IO2</sub> × (P<sub>B</sub>(h) − 47)</b>")),

  sidebarLayout(
    sidebarPanel(
      width = 3,

      h4("① 개인 표현형"),
      selectInput("pheno", "표현형 프리셋",
                  c("전형적 트레커" = "typ",
                    "HAPE 감수성" = "hapes",
                    "Tight-fit (낮은 두개내 순응도)" = "tight",
                    "엘리트 등반가 (높은 HVR)" = "elite",
                    "직접 설정" = "custom")),
      conditionalPanel(
        "input.pheno == 'custom'",
        sliderInput("hvr",  "HVR 배수 (저산소 환기반응)", 0.4, 2.5, 1.0, 0.05),
        sliderInput("lam",  "λ_max — 최대 HPV 수축 인자", 1, 14, 5, 0.5),
        sliderInput("ahet", "a — HPV 반응 혈관상 분율 (불균일성)", 0.05, 0.92, 0.50, 0.01),
        sliderInput("pvi",  "PVI — 두개척수 압력-용적 지수 (mL)", 10, 40, 25, 1),
        sliderInput("hbsl", "해수면 Hb (g/dL)", 11, 19, 15, 0.1)
      ),
      helpText(HTML("<small>모세혈관 과관류의 <b>상한은 1/[1−a(1−κ)]</b> 이다. ",
                    "HPV 의 세기는 도달 속도만 정한다.</small>")),

      hr(),
      h4("② 등반 프로파일"),
      selectInput("profile", "프로파일",
                  c("급속 상승 4559 m (Capanna Margherita)" = "rapid",
                    "단계 상승 4559 m (~400 m/일)"          = "graded",
                    "높이 오르고 낮게 자기"                 = "chsl",
                    "4000 m 취침"                           = "sleep4000",
                    "3800 m 3주 순응"                       = "accl3800",
                    "에베레스트 등정일 (8848 m)"            = "everest",
                    "직접 설정 (일정 고도)"                 = "flat")),
      conditionalPanel("input.profile == 'flat'",
                       sliderInput("flatalt", "고도 (m)", 0, 8848, 4500, 50)),
      sliderInput("tend", "시뮬레이션 기간 (h)", 24, 504, 120, 12),
      sliderInput("exer", "주간 운동강도 (× 안정시 V̇O₂)", 1, 4, 1, 0.1),
      helpText(HTML("<small>운동은 P<sub>cap</sub> 항 전체에 <b>곱해진다</b>. ",
                    "임상의가 실제로 조절할 수 있는 분기 파라미터.</small>")),

      hr(),
      h4("③ 약물 / 처치"),
      checkboxInput("acz", "아세타졸아미드", FALSE),
      conditionalPanel("input.acz",
                       selectInput("aczdose", "용량", c("125 mg bid" = 125,
                                                        "250 mg bid" = 250))),
      checkboxInput("dex", "덱사메타손 4 mg q12h", FALSE),
      checkboxInput("nif", "니페디핀 SR 30 mg bid", FALSE),
      checkboxInput("tad", "타다라필 10 mg bid", FALSE),
      checkboxInput("sal", "살메테롤 125 µg bid", FALSE),
      checkboxInput("ibu", "이부프로펜 600 mg tid", FALSE),
      sliderInput("fio2", "흡입 산소분율 F_IO₂", 0.2094, 0.60, 0.2094, 0.005),
      sliderInput("bag",  "가모우백 가압 (mmHg)", 0, 220, 0, 5),

      hr(),
      actionButton("go", "시뮬레이션 실행", class = "btn-primary btn-block")
    ),

    mainPanel(
      width = 9,
      tabsetPanel(
        id = "tabs", type = "tabs",

        ## ------------------------------------------------------------------
        tabPanel("① 환자 프로파일",
          h4("표현형 요약"),
          tableOutput("phenotab"),
          h4("이 표현형이 무엇을 바꾸는가"),
          plotOutput("phenoplot", height = "380px"),
          helpText(HTML(
            "왼쪽: 개인의 HVR 이 고도별 정상상태 S<sub>aO2</sub> 를 어디에 놓는지. ",
            "오른쪽: 불균일성 a 가 모세혈관 압력의 <b>천장</b>을 어디에 두는지 — ",
            "λ→∞ 에서 증폭 상한은 1/[1−a(1−κ)] 로 수렴한다 (κ = 0 일 때만 1/(1−a))."))
        ),

        ## ------------------------------------------------------------------
        tabPanel("② 등반 프로파일",
          plotOutput("profplot", height = "260px"),
          h4("고도가 하는 일은 이것뿐이다"),
          plotOutput("pio2plot", height = "300px"),
          tableOutput("statictab")
        ),

        ## ------------------------------------------------------------------
        tabPanel("③ 가스교환",
          plotOutput("gasplot", height = "460px"),
          h4("약물 농도"),
          plotOutput("pkplot", height = "240px")
        ),

        ## ------------------------------------------------------------------
        tabPanel("④ 환기 · 산-염기",
          plotOutput("ventplot", height = "460px"),
          helpText(HTML(
            "중추 무호흡 역치 <b>B<sub>c</sub></b> 는 CSF 중탄산이 정한다. ",
            "순응이란 이 역치가 내려가는 과정이고, 아세타졸아미드는 그것을 ",
            "신장 34시간 지연을 건너뛰고 곧바로 하는 약이다."))
        ),

        ## ------------------------------------------------------------------
        tabPanel("⑤ 수면 · CO₂ 예비량",
          fluidRow(
            column(6, plotOutput("resplot", height = "320px")),
            column(6, plotOutput("ahiplot", height = "320px"))),
          h4("CO₂ 예비량 = P_aCO₂ − B_c"),
          tableOutput("restab"),
          helpText(HTML(
            "아세타졸아미드는 호흡을 자극해서 주기성 호흡을 없애는 것이 아니다. ",
            "<b>작동점보다 바닥을 더 빨리 내려서</b> 예비량을 넓힌다. ",
            "이 표의 두 행을 비교하면 그것이 산수로 보인다."))
        ),

        ## ------------------------------------------------------------------
        tabPanel("⑥ 뇌 (AMS / HACE)",
          plotOutput("brainplot", height = "460px"),
          fluidRow(
            column(6, plotOutput("icpcurve", height = "300px")),
            column(6, tableOutput("braintab"))),
          helpText(HTML(
            "Monro–Kellie: ICP = ICP₀·10^(ΔV/PVI). 같은 ΔV 라도 PVI 가 ",
            "16 mL 인 사람과 30 mL 인 사람의 ICP 는 완전히 다르다 — ",
            "AMS 개인차의 대부분이 여기에 있다."))
        ),

        ## ------------------------------------------------------------------
        tabPanel("⑦ 폐 (HAPE)",
          plotOutput("lungplot", height = "440px"),
          h4("두 구획 모델: 심박출량과 불균일성의 등고선"),
          plotOutput("pcapcontour", height = "360px"),
          helpText(HTML(
            "P<sub>cap,open</sub> = P<sub>LA</sub> + Q̇·(r<sub>v</sub>/µ)/[a/λ+(1−a)/µ]. ",
            "Q̇ 는 <b>곱셈</b>으로 들어간다 — 그래서 HAPE 는 고도의 병이 아니라 ",
            "고도 × 운동의 병이다."))
        ),

        ## ------------------------------------------------------------------
        tabPanel("⑧ 임상 엔드포인트",
          plotOutput("endplot", height = "420px"),
          tableOutput("endtab"),
          helpText("AMS 정의(2018 Lake Louise): 두통 ≥1 이고 총점 ≥3.")
        ),

        ## ------------------------------------------------------------------
        tabPanel("⑨ 시나리오 비교",
          checkboxGroupInput("cmparms", "비교할 arm",
                             choices = c("예방 없음"        = "none",
                                         "아세타졸아미드 250" = "acz",
                                         "덱사메타손 4 bid"  = "dex",
                                         "단계 상승"          = "graded",
                                         "이부프로펜"        = "ibu"),
                             selected = c("none", "acz", "dex"), inline = TRUE),
          selectInput("cmpvar", "변수",
                      c("LLS", "SaO2", "PaCO2", "HCO3", "ICP", "AHI",
                        "CO2res", "EVLW", "mPAP")),
          plotOutput("cmpplot", height = "400px"),
          h4("★ 비대칭"),
          tableOutput("asymtab"),
          helpText(HTML(
            "두 약이 증상 점수를 비슷하게 낮춘다. 그런데 <b>산소를 움직인 것은 ",
            "하나뿐이다.</b> 「증상이 지속되면 하산하라」는 모든 지침은 증상을 ",
            "산소 게이지로 쓰고 있고, 덱사메타손은 그 게이지를 부수는 약이다."))
        ),

        ## ------------------------------------------------------------------
        tabPanel("⑩ 하산 등가",
          h4("모든 처치를 하나의 통화로"),
          sliderInput("deqalt", "기준 고도 (m)", 2500, 8000, 4559, 50),
          plotOutput("deqplot", height = "380px"),
          tableOutput("deqtab"),
          helpText(HTML(
            "덱사메타손이 <b>0 m</b> 로 나오는 것은 비판이 아니라 설명이다. ",
            "가스교환 사슬의 어느 항도 건드리지 않기 때문이며, 그것이 바로 ",
            "이 약이 위험한 이유다."))
        ),

        ## ------------------------------------------------------------------
        tabPanel("⑪ 바이오마커",
          plotOutput("biomarker", height = "460px"),
          h4("1주차의 Hct 상승은 적혈구가 아니라 물이다"),
          plotOutput("hbsplit", height = "280px")
        ),

        ## ------------------------------------------------------------------
        tabPanel("⑫ 이분기 · 민감도",
          h4("HAPE 임계 고도"),
          fluidRow(
            column(6, sliderInput("bifq", "심박출량 Q̇ (L/min)", 5, 22, 6, 0.5)),
            column(6, sliderInput("bifa", "불균일성 a", 0.1, 0.92, 0.5, 0.01))),
          plotOutput("bifplot", height = "340px"),
          h4("최적 헤마토크리트 — 고도와 무관하다"),
          plotOutput("hctplot", height = "300px"),
          helpText(HTML(
            "Ḋ<sub>O2</sub> ∝ Hct·exp(−k·γ·Hct) ⇒ Hct* = 1/(k·γ). ",
            "S<sub>aO2</sub> 가 소거되므로 <b>최적값은 고도에 따라 움직이지 ",
            "않는다.</b> 안데스 주민의 Hct 55 % 를 정당화하려면 순환이 점도 ",
            "부담의 대부분을 흡수해야 한다(γ ≤ 0.79)."))
        )
      )
    )
  )
)

## =============================================================================
##  SERVER
## =============================================================================
server <- function(input, output, session) {

  pheno <- reactive({
    switch(input$pheno,
      typ   = TYPICAL,
      hapes = HAPE_SUSC,
      tight = TIGHT_FIT,
      elite = ELITE,
      custom = phenotype("직접 설정", HVRSUB = input$hvr, LAMMAX = input$lam,
                         AHET = input$ahet, PVI = input$pvi, HBSL = input$hbsl))
  })

  alt_fun <- reactive({
    switch(input$profile,
      rapid     = RAPID,
      graded    = GRADED,
      chsl      = function(t) {
                    b <- approx(c(0, 6, 24, 30, 48, 54, 72),
                                c(1130, 2400, 2400, 3000, 3000, 3600, 3600),
                                xout = t, rule = 2)$y
                    d <- t %% 24; if (d >= 10 && d <= 16) b + 700 else b },
      sleep4000 = ramp(list(c(0, 1130), c(6, 4000))),
      accl3800  = ramp(list(c(0, 200), c(10, 3800))),
      everest   = ramp(list(c(0, 5300), c(240, 5300), c(300, 6400),
                            c(400, 7100), c(450, 7900), c(470, 8848),
                            c(476, 7900))),
      flat      = function(t) input$flatalt)
  })

  doses <- reactive({
    n <- ceiling(input$tend/12) + 1
    d <- NULL
    add <- function(d, rows) if (is.null(d)) rows else dplyr::bind_rows(d, rows)
    if (input$acz) d <- add(d, dose_rows(0, 12, n, "ACZA", as.numeric(input$aczdose)))
    if (input$dex) d <- add(d, dose_rows(0, 12, n, "DEXA", 4))
    if (input$nif) d <- add(d, dose_rows(0, 12, n, "NIFA", 30))
    if (input$tad) d <- add(d, dose_rows(0, 12, n, "TADA", 10))
    if (input$sal) d <- add(d, dose_rows(0, 12, n, "SALE", 1.4))
    if (input$ibu) d <- add(d, dose_rows(0,  8, ceiling(input$tend/8)+1, "IBUA", 600))
    d
  })

  sim <- eventReactive(input$go, {
    withProgress(message = "적분 중...", {
      run_scenario(pheno(), input$tend, alt_fun(),
                   exer_fun = function(t) {
                     d <- t %% 24
                     if (d >= 9 && d <= 15) input$exer else 1 },
                   fio2_fun = function(t) input$fio2,
                   bag_fun  = function(t) input$bag,
                   doses = doses())
    })
  }, ignoreNULL = FALSE)

  long <- function(df, vars) df %>% select(time, all_of(vars)) %>%
    pivot_longer(-time)

  ## ---- ① phenotype --------------------------------------------------------
  output$phenotab <- renderTable({
    p <- pheno()$par
    data.frame(항목 = c("표현형", "HVR 배수", "λ_max", "불균일성 a",
                        "과관류 천장 1/[1−a(1−κ)]", "PVI (mL)", "해수면 Hb"),
               값 = c(pheno()$label, sprintf("%.2f", p$HVRSUB),
                      sprintf("%.1f", p$LAMMAX), sprintf("%.2f", p$AHET),
                      sprintf("%.2f", 1/(1 - p$AHET*(1 - 0.08))),
                      sprintf("%.0f", p$PVI), sprintf("%.1f", p$HBSL)))
  }, striped = TRUE)

  output$phenoplot <- renderPlot({
    alts <- seq(0, 8500, 250)
    p <- pheno()$par
    d1 <- do.call(rbind, lapply(c(0.6, 1.0, 1.5), function(h)
      data.frame(alt = alts, HVR = factor(h),
                 SaO2 = 100*severinghaus(pmax(1, pao2_alt(alts, 40 - 0.0028*alts*h) - 5)))))
    g1 <- ggplot(d1, aes(alt, SaO2, colour = HVR)) + geom_line(linewidth = 1) +
      geom_vline(xintercept = 4559, linetype = 2, colour = "grey50") +
      scale_colour_manual(values = PAL) +
      labs(x = "고도 (m)", y = "S_aO2 (%)", title = "HVR 이 정상상태를 어디에 놓는가") + THEME
    lam <- 10^seq(0, 3, length.out = 120)
    d2 <- do.call(rbind, lapply(c(0.25, 0.50, 0.75, 0.85), function(a)
      data.frame(lambda = lam, a = factor(a),
                 Pcap = vapply(lam, function(l) pcap_two_bed(l, a), numeric(1)),
                 ceiling = 8 + 6*0.5/(1 - a))))
    g2 <- ggplot(d2, aes(lambda, Pcap, colour = a)) + geom_line(linewidth = 1) +
      geom_hline(yintercept = 19.5, linetype = 2, colour = "#c62828") +
      scale_x_log10() + scale_colour_manual(values = PAL) +
      labs(x = "λ (HPV 수축 인자, 로그)", y = "P_cap,open (mmHg)",
           title = "천장은 a 가 정한다 (붉은 선 = 응력파괴 역치)") + THEME
    gridExtra::grid.arrange(g1, g2, ncol = 2)
  })

  ## ---- ② profile ----------------------------------------------------------
  output$profplot <- renderPlot({
    tt <- seq(0, input$tend, 0.25)
    data.frame(time = tt, alt = vapply(tt, alt_fun(), numeric(1))) %>%
      ggplot(aes(time, alt)) + geom_line(linewidth = 1, colour = PAL[1]) +
      labs(x = "시간 (h)", y = "고도 (m)") + THEME
  })

  output$pio2plot <- renderPlot({
    a <- seq(0, 8848, 50)
    data.frame(alt = a,
               `P_B`    = pb_west(a),
               `P_IO2`  = pio2_alt(a),
               `P_AO2 (PaCO2 40)` = pao2_alt(a, 40),
               check.names = FALSE) %>%
      pivot_longer(-alt) %>%
      ggplot(aes(alt, value, colour = name)) + geom_line(linewidth = 1) +
      scale_colour_manual(values = PAL) +
      labs(x = "고도 (m)", y = "mmHg", colour = NULL) + THEME
  })

  output$statictab <- renderTable({
    a <- c(0, 1600, 2500, 3500, 4559, 5300, 6400, 7100, 8000, 8848)
    data.frame(`고도 (m)` = a, `P_B` = round(pb_west(a), 1),
               `P_IO2` = round(pio2_alt(a), 1),
               `해수면 대비 (%)` = round(100*pio2_alt(a)/pio2_alt(0), 1),
               check.names = FALSE)
  }, striped = TRUE)

  ## ---- ③ gas exchange -----------------------------------------------------
  output$gasplot <- renderPlot({
    long(sim(), c("ALT", "PIO2", "PAO2", "PaO2", "PaCO2", "SaO2",
                  "AaDO2", "CaO2", "DO2")) %>%
      ggplot(aes(time, value)) + geom_line(linewidth = 0.8, colour = PAL[1]) +
      facet_wrap(~name, scales = "free_y", ncol = 3) +
      labs(x = "시간 (h)", y = NULL) + THEME
  })

  output$pkplot <- renderPlot({
    long(sim(), c("CACZ", "CDEX")) %>%
      ggplot(aes(time, value, colour = name)) + geom_line(linewidth = 0.9) +
      scale_colour_manual(values = PAL) +
      labs(x = "시간 (h)", y = "mg/L", colour = NULL) + THEME
  })

  ## ---- ④ ventilation / acid-base ------------------------------------------
  output$ventplot <- renderPlot({
    long(sim(), c("VE", "VA", "PaCO2", "HCO3", "pHa", "Bc", "Bp", "P50a", "SaO2")) %>%
      ggplot(aes(time, value)) + geom_line(linewidth = 0.8, colour = PAL[3]) +
      facet_wrap(~name, scales = "free_y", ncol = 3) +
      labs(x = "시간 (h)", y = NULL) + THEME
  })

  ## ---- ⑤ sleep -------------------------------------------------------------
  output$resplot <- renderPlot({
    d <- sim()
    ggplot(d, aes(time)) +
      geom_line(aes(y = PaCO2, colour = "P_aCO2"), linewidth = 0.9) +
      geom_line(aes(y = Bc,    colour = "B_c (무호흡 역치)"), linewidth = 0.9) +
      geom_ribbon(aes(ymin = pmin(Bc, PaCO2), ymax = pmax(Bc, PaCO2)),
                  fill = "#90caf9", alpha = 0.35) +
      scale_colour_manual(values = PAL) +
      labs(x = "시간 (h)", y = "mmHg", colour = NULL,
           title = "음영 = CO₂ 예비량") + THEME
  })

  output$ahiplot <- renderPlot({
    sim() %>% mutate(night = floor(time/24),
                     sl = (time %% 24) >= 22 | (time %% 24) < 6) %>%
      filter(sl) %>% group_by(night) %>%
      summarise(AHI = mean(AHI), .groups = "drop") %>%
      ggplot(aes(factor(night + 1), AHI)) +
      geom_col(fill = PAL[5]) +
      labs(x = "밤 (night)", y = "예측 AHI (/h)") + THEME
  })

  output$restab <- renderTable({
    alts <- c(0, 2500, 3500, 4000, 4559, 5300)
    do.call(rbind, lapply(alts, function(a) {
      data.frame(`고도 (m)` = a, check.names = FALSE,
                 `대략적 CO₂ 예비량 (mmHg)` =
                   round(4.11 - (4.11 - 1.53)*(a/5300), 2))
    }))
  }, striped = TRUE)

  ## ---- ⑥ brain -------------------------------------------------------------
  output$brainplot <- renderPlot({
    long(sim(), c("ICP", "LLS", "SaO2", "PaCO2", "HACER", "ALT")) %>%
      ggplot(aes(time, value)) + geom_line(linewidth = 0.8, colour = PAL[7]) +
      facet_wrap(~name, scales = "free_y", ncol = 3) +
      labs(x = "시간 (h)", y = NULL) + THEME
  })

  output$icpcurve <- renderPlot({
    dv <- seq(0, 22, 0.2)
    do.call(rbind, lapply(c(16, 20, 25, 30), function(pv)
      data.frame(dV = dv, PVI = factor(pv), ICP = pmin(80, 10*10^(dv/pv)))) ) %>%
      ggplot(aes(dV, ICP, colour = PVI)) + geom_line(linewidth = 1) +
      geom_hline(yintercept = 22, linetype = 2, colour = "#c62828") +
      scale_colour_manual(values = PAL) +
      labs(x = "두개내 용적 증가 ΔV (mL)", y = "ICP (mmHg)",
           title = "같은 부종, 다른 사람") + THEME
  })

  output$braintab <- renderTable({
    d <- sim()
    data.frame(지표 = c("최고 ICP (mmHg)", "최고 LLS", "AMS 시간 (h)",
                        "최대 HACE 위험", "최저 S_aO2 (%)"),
               값 = c(sprintf("%.1f", max(d$ICP)), sprintf("%.2f", max(d$LLS)),
                      sprintf("%.0f", sum(d$AMS)*(d$time[2]-d$time[1])),
                      sprintf("%.3f", max(d$HACER)), sprintf("%.1f", min(d$SaO2))))
  }, striped = TRUE)

  ## ---- ⑦ lung --------------------------------------------------------------
  output$lungplot <- renderPlot({
    long(sim(), c("mPAP", "PCAPOP", "LAMBDA", "EVLW", "SaO2", "QC")) %>%
      ggplot(aes(time, value)) + geom_line(linewidth = 0.8, colour = PAL[6]) +
      facet_wrap(~name, scales = "free_y", ncol = 3) +
      labs(x = "시간 (h)", y = NULL) + THEME
  })

  output$pcapcontour <- renderPlot({
    grid <- expand.grid(Q = seq(4, 22, 0.5), a = seq(0.1, 0.9, 0.02))
    grid$Pcap <- mapply(function(Q, a) pcap_two_bed(6, a, Q = Q), grid$Q, grid$a)
    ggplot(grid, aes(Q, a, fill = Pcap)) + geom_raster() +
      geom_contour(aes(z = Pcap), breaks = 19.5, colour = "white", linewidth = 1.2) +
      scale_fill_viridis_c(option = "magma") +
      labs(x = "심박출량 Q̇ (L/min)", y = "불균일성 a",
           fill = "P_cap\n(mmHg)",
           title = "흰 선 = 응력파괴 역치 19.5 mmHg (λ = 6 고정)") + THEME
  })

  ## ---- ⑧ endpoints ---------------------------------------------------------
  output$endplot <- renderPlot({
    long(sim(), c("LLS", "AMS", "HACER", "EVLW", "SaO2", "AHI")) %>%
      ggplot(aes(time, value)) + geom_line(linewidth = 0.8, colour = PAL[4]) +
      facet_wrap(~name, scales = "free_y", ncol = 3) +
      labs(x = "시간 (h)", y = NULL) + THEME
  })

  output$endtab <- renderTable({ summarise_run(sim()) }, digits = 2, striped = TRUE)

  ## ---- ⑨ comparison --------------------------------------------------------
  arms <- reactive({
    n <- ceiling(120/12) + 1
    L <- list()
    if ("none"   %in% input$cmparms) L[["예방 없음"]] <- run_scenario(pheno(), 120, RAPID)
    if ("acz"    %in% input$cmparms) L[["아세타졸아미드 250"]] <-
      run_scenario(pheno(), 120, RAPID, doses = dose_rows(0, 12, n, "ACZA", 250))
    if ("dex"    %in% input$cmparms) L[["덱사메타손 4 bid"]] <-
      run_scenario(pheno(), 120, RAPID, doses = dose_rows(0, 12, n, "DEXA", 4))
    if ("graded" %in% input$cmparms) L[["단계 상승"]] <-
      run_scenario(pheno(), 144, GRADED)
    if ("ibu"    %in% input$cmparms) L[["이부프로펜"]] <-
      run_scenario(pheno(), 120, RAPID, doses = dose_rows(0, 8, 16, "IBUA", 600))
    L
  })

  output$cmpplot <- renderPlot({
    L <- arms(); validate(need(length(L) > 0, "arm 을 하나 이상 선택하세요"))
    bind_rows(lapply(names(L), function(n)
      data.frame(time = L[[n]]$time, value = L[[n]][[input$cmpvar]], arm = n))) %>%
      ggplot(aes(time, value, colour = arm)) + geom_line(linewidth = 0.9) +
      scale_colour_manual(values = PAL) +
      labs(x = "시간 (h)", y = input$cmpvar, colour = NULL) + THEME
  })

  output$asymtab <- renderTable({
    L <- arms(); validate(need(length(L) > 0, ""))
    do.call(rbind, lapply(names(L), function(n) {
      d <- L[[n]]
      data.frame(arm = n,
                 `최저 S_aO2 (%)` = round(min(d$SaO2), 1),
                 `최종 HCO3⁻`     = round(tail(d$HCO3, 1), 1),
                 `최고 LLS`       = round(max(d$LLS), 2),
                 `AMS 시간`       = round(sum(d$AMS)*(d$time[2]-d$time[1])),
                 check.names = FALSE)
    }))
  }, striped = TRUE)

  ## ---- ⑩ descent equivalents ----------------------------------------------
  deq <- reactive({
    a0 <- input$deqalt
    ## 정상상태 근사: 급성 (해수면 산-염기) 조건에서 SaO2(고도)
    sao2_of <- function(a, fio2 = 0.2094, bag = 0) {
      paco2 <- 40 - 0.0026*a
      pao2  <- pio2_alt(a, fio2, bag) - paco2*(fio2 + (1 - fio2)/0.85)
      severinghaus(max(pao2 - 6, 1))
    }
    base <- sao2_of(a0)
    opts <- list(
      `아무것도 안 함`              = base,
      `아세타졸아미드 250 bid`      = sao2_of(a0) + 0.035,
      `완전 순응 (약 1주)`          = sao2_of(a0) + 0.060,
      `산소 F_IO2 0.28`             = sao2_of(a0, fio2 = 0.28),
      `산소 F_IO2 0.35`             = sao2_of(a0, fio2 = 0.35),
      `가모우백 2 psi (+105 mmHg)`  = sao2_of(a0, bag = 105),
      `가모우백 4 psi (+207 mmHg)`  = sao2_of(a0, bag = 207),
      `덱사메타손`                  = base)
    data.frame(
      처치 = names(opts),
      `S_aO2 (%)` = round(100*unlist(opts), 1),
      `하산 등가 (m)` = vapply(unlist(opts), function(s)
        tryCatch(descent_equivalent(a0, s, sao2_of), error = function(e) NA_real_),
        numeric(1)),
      check.names = FALSE)
  })

  output$deqtab <- renderTable({ deq() }, digits = 0, striped = TRUE)

  output$deqplot <- renderPlot({
    d <- deq()
    d$처치 <- factor(d$처치, levels = d$처치[order(d$`하산 등가 (m)`)])
    ggplot(d, aes(처치, `하산 등가 (m)`)) +
      geom_col(fill = PAL[3]) + coord_flip() +
      labs(x = NULL, y = "몇 미터의 하산에 해당하는가") + THEME
  })

  ## ---- ⑪ biomarkers --------------------------------------------------------
  output$biomarker <- renderPlot({
    long(sim(), c("HB", "HCT", "DO2", "CaO2", "mPAP", "AHI")) %>%
      ggplot(aes(time, value)) + geom_line(linewidth = 0.8, colour = PAL[2]) +
      facet_wrap(~name, scales = "free_y", ncol = 3) +
      labs(x = "시간 (h)", y = NULL) + THEME
  })

  output$hbsplit <- renderPlot({
    d <- sim()
    data.frame(time = d$time,
               `실제 측정되는 Hb` = d$HB,
               `혈장량이 변하지 않았다면` = d$HB[1]) %>%
      pivot_longer(-time) %>%
      ggplot(aes(time, value, colour = name)) + geom_line(linewidth = 1) +
      scale_colour_manual(values = PAL) +
      labs(x = "시간 (h)", y = "Hb (g/dL)", colour = NULL) + THEME
  })

  ## ---- ⑫ bifurcation -------------------------------------------------------
  output$bifplot <- renderPlot({
    a <- seq(1500, 7500, 100)
    paco2 <- 40 - 0.0026*a
    pao2  <- pio2_alt(a) - paco2*(0.2094 + (1 - 0.2094)/0.85)
    lam   <- 1 + input$bifa*0 + pheno()$par$LAMMAX/(1 + (pmax(pao2, 1)/45)^4)
    pc    <- mapply(function(l) pcap_two_bed(l, input$bifa, Q = input$bifq), lam)
    data.frame(alt = a, Pcap = pc) %>%
      ggplot(aes(alt, Pcap)) + geom_line(linewidth = 1.1, colour = PAL[6]) +
      geom_hline(yintercept = 19.5, linetype = 2, colour = "#c62828") +
      labs(x = "고도 (m)", y = "P_cap,open (mmHg)",
           title = "붉은 선을 넘는 고도 = 이 표현형의 HAPE 임계 고도") + THEME
  })

  output$hctplot <- renderPlot({
    hct <- seq(0.20, 0.80, 0.005)
    do.call(rbind, lapply(c(1.0, 0.9, 0.8, 0.7), function(g)
      data.frame(Hct = hct, gamma = factor(g),
                 DO2 = hct*exp(-2.31*g*hct)/max(hct*exp(-2.31*g*hct))))) %>%
      ggplot(aes(100*Hct, DO2, colour = gamma)) + geom_line(linewidth = 1) +
      geom_vline(xintercept = 43.3, linetype = 2, colour = "grey40") +
      scale_colour_manual(values = PAL) +
      labs(x = "헤마토크리트 (%)", y = "상대 Ḋ_O2", colour = "γ") + THEME
  })
}

shinyApp(ui, server)
