## =====================================================================
##  G6PD DEFICIENCY — interactive QSP dashboard (Shiny)
##  포도당-6-인산 탈수소효소 결핍증 · 인터랙티브 시뮬레이터
##
##  Run with:
##      shiny::runApp("g6pd_shiny_app.R")
##  (the app sources g6pd_mrgsolve_model.R from the same directory)
##
##  The app is organised around ONE control that matters more than all the
##  others put together: the variant selector, which sets exactly two
##  numbers, E0 and TAU. Tab 1 draws the curve those two numbers make.
##  Every other tab is a consequence of it.
##
##  Tabs
##    1  환자·변이 (Patient & variant)   the age-activity curve, and a*
##    2  약동학 (Pharmacokinetics)       drug concentrations, oxidant flux
##    3  적혈구 연령 구조 (Age structure) the histogram with a* drawn on it
##    4  혈액학 (Haematology)             Hb, reticulocytes, EPO
##    5  시나리오 비교 (Scenarios)         side-by-side, incl. the 1-vs-2 case
##    6  바이오마커 (Biomarkers)           MetHb, bilirubin, haptoglobin, AKI
##    7  진단의 함정 (The assay trap)      what the laboratory would report
##    8  신생아 황달 (Neonatal jaundice)   the G6PD x UGT1A1 product
## =====================================================================

library(shiny)
library(mrgsolve)
library(ggplot2)
library(dplyr)
library(tidyr)

## NOTE ON MASKING. mrgsolve exports req(), param(), plot() and filter(), and
## it is attached AFTER shiny, so a bare req() in the server would resolve to
## mrgsolve::req and fail at runtime (a parse check does not catch this).
## dplyr is attached after mrgsolve, so filter() correctly resolves to dplyr.
## Anything ambiguous below is namespace-qualified.

source("g6pd_mrgsolve_model.R", local = TRUE)

THEME <- theme_bw(base_size = 13) +
  theme(panel.grid.minor = element_blank(),
        strip.background = element_rect(fill = "grey93"),
        legend.position = "bottom")

VARIANT_CHOICES <- setNames(names(VARIANTS),
                            vapply(VARIANTS, function(v) v$label, ""))

DRUG_CHOICES <- c("없음 (none)"                       = "none",
                  "프리마퀸 매일 (primaquine daily)"   = "pq",
                  "프리마퀸 주 1회 (primaquine weekly)"= "pqw",
                  "타페노퀸 단회 (tafenoquine)"        = "tq",
                  "답손 (dapsone)"                     = "dap",
                  "라스부리카제 (rasburicase, TLS)"    = "rbx",
                  "잠두 (fava bean meal)"              = "fava")

## ---------------------------------------------------------------------
##  UI
## ---------------------------------------------------------------------
ui <- fluidPage(
  titlePanel("G6PD 결핍증 QSP 시뮬레이터 — 효소는 숫자가 아니라 세포 연령의 함수다"),
  tags$p(style = "color:#666; margin-top:-8px;",
         "Glucose-6-phosphate dehydrogenase deficiency · age-structured ",
         "erythrocyte redox model · 교육/연구용, 임상 사용 불가"),

  sidebarLayout(
    sidebarPanel(
      width = 3,

      h4("1. 유전형 (the whole disease)"),
      selectInput("variant", "G6PD 변이", VARIANT_CHOICES, selected = "Aminus"),
      helpText(HTML("변이가 바꾸는 것은 <b>딱 두 숫자</b>: 새 적혈구의 활성 E0와,
                     세포 안에서 효소가 사라지는 반감기 TAU. 나머지는 전부
                     여기서 파생된다.")),
      checkboxInput("custom", "E0 / TAU 직접 지정", FALSE),
      conditionalPanel("input.custom",
        sliderInput("E0",  "E0 (새 적혈구 활성, 정상=1)", 0.005, 1.0, 0.55, 0.005),
        sliderInput("TAU", "TAU (세포 내 효소 반감기, 일)", 3, 80, 13, 1)),

      hr(),
      h4("2. 환자"),
      numericInput("wt", "체중 (kg)", 70, 3, 150, 1),
      selectInput("ugt", "UGT1A1 프로모터 (TA)n", c("6/6", "6/7", "7/7"), "6/6"),
      sliderInput("cyp2d6", "CYP2D6 활성 점수", 0, 2, 1, 0.25),
      selectInput("nat2", "NAT2 아세틸화", c("fast", "slow"), "fast"),
      sliderInput("marrow", "골수 예비능 (1=정상, 0.15=재생불량 위기)",
                  0.05, 1.5, 1.0, 0.05),
      checkboxInput("splenx", "비장절제 상태", FALSE),
      checkboxInput("infect", "동반 감염 (산화 스트레스 추가)", FALSE),

      hr(),
      h4("3. 노출"),
      selectInput("drug", "약물 / 유발인자", DRUG_CHOICES, selected = "pq"),
      conditionalPanel("input.drug == 'pq'",
        sliderInput("pq_dose", "프리마퀸 (mg/일)", 7.5, 60, 30, 7.5),
        sliderInput("pq_days", "투여 기간 (일)", 1, 60, 60, 1)),
      conditionalPanel("input.drug == 'pqw'",
        sliderInput("pqw_dose", "프리마퀸 (mg, 주 1회)", 15, 60, 45, 7.5),
        sliderInput("pqw_wk", "주 수", 1, 12, 8, 1)),
      conditionalPanel("input.drug == 'tq'",
        sliderInput("tq_dose", "타페노퀸 단회 (mg)", 50, 600, 300, 50)),
      conditionalPanel("input.drug == 'dap'",
        sliderInput("dap_dose", "답손 (mg/일)", 25, 300, 100, 25),
        checkboxInput("mb", "메틸렌블루 2 mg/kg 투여 (30일째)", FALSE)),
      conditionalPanel("input.drug == 'rbx'",
        sliderInput("urate0", "초기 요산 (mg/dL)", 4, 30, 20, 1),
        sliderInput("rbx_days", "라스부리카제 투여일 수", 1, 5, 3, 1)),
      conditionalPanel("input.drug == 'fava'",
        sliderInput("fava_amt", "흡수 가능 아글리콘 (mg)", 100, 3000, 1000, 100)),

      hr(),
      sliderInput("days", "관찰 기간 (일)", 20, 180, 60, 5),
      actionButton("go", "시뮬레이션 실행", class = "btn-primary btn-block")
    ),

    mainPanel(
      width = 9,
      tabsetPanel(
        id = "tabs", type = "tabs",

        ## ---------------------------------------------------------------
        tabPanel("1 · 환자·변이",
          br(),
          fluidRow(
            column(7, h4("E(a) = E0 · exp(−ln2 · a / TAU)"),
                      plotOutput("p_curve", height = 330)),
            column(5, h4("임계 연령 a*"),
                      plotOutput("p_astar", height = 330))),
          hr(),
          h4("이 두 곡선이 이 앱의 전부"),
          HTML("<p>왼쪽: 성숙 적혈구는 핵도 리보솜도 없어 효소를 새로 만들 수
                없다. 효소는 한 번 실리고 그 뒤로는 <b>줄기만</b> 한다. 그래서
                활성은 세포 연령의 지수함수가 되고, 변이란 그 <b>기울기</b>다.</p>
                <p>오른쪽: 산화제 부하 OX가 주어지면 세포는
                <code>VMAXOX·E(a) &gt; OX</code>일 때만 살아남는다. 그 경계가
                임계 연령 <b>a*</b>이고, 닫힌 형태로 풀린다:</p>
                <p style='margin-left:2em'><code>a* = (TAU / ln2) ·
                ln(E0 · VMAXOX / OX<sub>eff</sub>)</code>&nbsp;&nbsp;&nbsp;
                <small>(OX<sub>eff</sub> = OX − KPIT·HZ50, 제거를 유발하지
                않고 견딜 수 있는 작은 부하를 뺀 값)</small></p>
                <p>이 식에 <b>적합된 용혈 속도상수가 없다</b>는 점이 요점이다.
                들어 있는 것은 유전형(E0, TAU), 생화학적 용량(VMAXOX),
                그리고 노출(OX)뿐이다.</p>
                <p>용혈은 pool에 곱해지는 속도상수가 아니라 <b>연령 히스토그램을
                지나가는 칼날</b>이다. A−는 기울기가 가팔라서 a*가 분포
                <i>안쪽</i>에 떨어지고 젊고 저항성 있는 생존자를 남긴다 —
                그래서 약을 계속 먹어도 회복한다. Mediterranean은 곡선이
                납작하고 낮아서 남길 생존자가 없다.</p>"),
          hr(),
          tableOutput("t_variants")),

        ## ---------------------------------------------------------------
        tabPanel("2 · 약동학",
          br(),
          fluidRow(column(6, plotOutput("p_pk", height = 300)),
                   column(6, plotOutput("p_ox", height = 300))),
          hr(),
          HTML("<p><b>반감기가 곧 독성학이다.</b> 프리마퀸은 반감기 약 7시간이라
                용혈이 시작되면 중단이 실제로 작동한다. 타페노퀸은 약 15일이라
                <i>중단할 것이 없다</i> — 그리고 망상적혈구 구조 요청은 도착에
                4–7일이 걸린다. 같은 총 노출량을 빠르게 끄는 약으로 주느냐
                느리게 끄는 약으로 주느냐에 따라 최저 혈색소가 달라지는 이유가
                이것이고, FDA가 타페노퀸에만 정량 검사를 요구하는 이유이기도
                하다.</p>"),
          plotOutput("p_astar_t", height = 280)),

        ## ---------------------------------------------------------------
        tabPanel("3 · 적혈구 연령 구조",
          br(),
          sliderInput("snap", "관찰 시점 (일)", 0, 60, 8, 1, width = "70%"),
          plotOutput("p_hist", height = 340),
          hr(),
          plotOutput("p_heat", height = 300),
          HTML("<p>위: 그 시점의 적혈구 연령 히스토그램. 붉은 세로선이 a*이며,
                그 오른쪽의 세포들이 지금 파괴되고 있는 세포다.
                아래: 시간 × 연령 평면에서 각 연령대 세포 수의 변화. A−에서는
                오래된 연령대만 비워지고 젊은 쪽이 다시 채워지는 것이 보인다.</p>")),

        ## ---------------------------------------------------------------
        tabPanel("4 · 혈액학",
          br(),
          fluidRow(column(6, plotOutput("p_hb",   height = 300)),
                   column(6, plotOutput("p_ret",  height = 300))),
          fluidRow(column(6, plotOutput("p_epo",  height = 280)),
                   column(6, plotOutput("p_lys",  height = 280))),
          hr(), tableOutput("t_summary")),

        ## ---------------------------------------------------------------
        tabPanel("5 · 시나리오 비교",
          br(),
          checkboxGroupInput("cmp", "비교할 변이 (같은 노출을 각각에 적용)",
            choices = VARIANT_CHOICES,
            selected = c("normal", "Aminus", "Mediterranean"), inline = TRUE),
          actionButton("go_cmp", "비교 실행", class = "btn-primary"),
          br(), br(),
          plotOutput("p_cmp", height = 380),
          hr(),
          tableOutput("t_cmp"),
          HTML("<p><b>이 탭이 모델의 핵심 결과다.</b> 같은 약, 같은 용량,
                둘 다 '중증 결핍'이라는 같은 라벨. A−는 바닥을 친 뒤
                <i>약을 먹는 중에</i> 올라오고, Mediterranean은 올라오지 않는다.
                차이를 만드는 파라미터는 붕괴 상수 하나뿐이다.</p>")),

        ## ---------------------------------------------------------------
        tabPanel("6 · 바이오마커",
          br(),
          fluidRow(column(6, plotOutput("p_met",  height = 280)),
                   column(6, plotOutput("p_bili", height = 280))),
          fluidRow(column(6, plotOutput("p_hapto",height = 280)),
                   column(6, plotOutput("p_aki",  height = 280))),
          hr(),
          HTML("<p><b>메틸렌블루가 왜 금기인가는 외울 사실이 아니라 계산
                결과다.</b> NADPH 풀 하나가 두 소비자를 먹인다: 글루타티온
                환원효소(산화 방어)와 NADPH-메트헤모글로빈 환원효소(메틸렌블루를
                류코메틸렌블루로 바꾸는 경로). 결핍에서는 두 번째가 돌지 않으니
                메틸렌블루는 <i>해독제로 실패하고 동시에 새 산화제가 된다</i>.
                답손 탭에서 메틸렌블루를 켜고 정상형과 A−를 비교해 보라.</p>")),

        ## ---------------------------------------------------------------
        tabPanel("7 · 진단의 함정",
          br(),
          plotOutput("p_assay", height = 360),
          hr(),
          HTML("<p>정량 G6PD 검사는 순환 중인 적혈구 <b>전체의 연령 가중
                평균</b>을 보고한다. 그런데 용혈은 방금 늙고 효소가 없는 세포를
                지워버렸고, 그 자리를 <b>젊은 세포 활성의 1.5배</b>를 지닌
                망상적혈구가 채웠다. 그래서 검사 수치는 용혈 중과 직후에
                <b>올라간다</b>.</p>
                <p>위 그래프에서 진짜 A− 환자가 최저점 부근에서 정상 범위로
                판독되는 것을 볼 수 있다. 그래서 급성 삽화 중의 정상 결과는
                결핍을 배제하지 못하며, 재검은 <b>3개월 뒤</b> — 적혈구 수명
                한 바퀴 — 해야 한다.</p>"),
          tableOutput("t_assay")),

        ## ---------------------------------------------------------------
        tabPanel("8 · 신생아 황달",
          br(),
          p("3 kg 만삭아, 생후 12시간부터 12일간. G6PD와 UGT1A1은 각각",
            "생성(분자)과 포합(분모)을 건드리므로 효과는 더해지는 것이 아니라",
            "곱해진다. 네 조합을 동시에 돌린다."),
          actionButton("go_neo", "신생아 4군 시뮬레이션", class = "btn-primary"),
          br(), br(),
          plotOutput("p_neo", height = 380),
          hr(),
          tableOutput("t_neo"),
          HTML("<p>둘 중 하나만 있으면 광선치료 문턱 근처에서 멈추고, 둘이
                함께 있을 때만 교환수혈 영역으로 넘어간다. G6PD 결핍이 전 세계
                핵황달의 주요 원인인 이유는 결핍 자체가 특별히 강력해서가
                아니라 <b>흔한 두 번째 결함과 곱해지기</b> 때문이다.</p>"))
      )
    )
  )
)

## ---------------------------------------------------------------------
##  SERVER
## ---------------------------------------------------------------------
server <- function(input, output, session) {

  patient <- reactive({
    p <- g6pd_patient(variant = input$variant, wt = input$wt, ugt = input$ugt,
                      cyp2d6 = input$cyp2d6, nat2 = input$nat2,
                      spleen = if (input$splenx) 0 else 1,
                      marrow = input$marrow)
    if (isTRUE(input$custom)) { p$E0 <- input$E0; p$TAU <- input$TAU }
    if (isTRUE(input$infect)) p$INFON <- 1
    p
  })

  regimen <- reactive({
    switch(input$drug,
      none  = NULL,
      pq    = rx_primaquine(input$pq_dose, days = input$pq_days),
      pqw   = rx_pq_weekly(input$pqw_dose, weeks = input$pqw_wk),
      tq    = rx_tafenoquine(input$tq_dose),
      dap   = if (isTRUE(input$mb))
                c(rx_dapsone(input$dap_dose, days = input$days),
                  rx_methyleneblue(2 * input$wt, start = 30))
              else rx_dapsone(input$dap_dose, days = input$days),
      rbx   = rx_rasburicase(0.2 * input$wt, days = input$rbx_days),
      fava  = rx_fava(input$fava_amt))
  })

  extra <- reactive({
    if (input$drug == "rbx") {
      ## urate mg/dL -> mmol in a 0.5 L/kg volume
      vur <- 0.5 * input$wt
      list(VUR = vur, URATE_0 = input$urate0 * 10 * vur / 168.1,
           KGENUR = 3 * input$wt / 70 * (input$urate0 / 6))
    } else list()
  })

  sim <- eventReactive(input$go, {
    withProgress(message = "시뮬레이션 중...", {
      run_g6pd(patient = patient(), rx = regimen(), days = input$days,
               extra = extra(), delta = 0.05,
               name = VARIANTS[[input$variant]]$label)
    })
  }, ignoreNULL = FALSE)

  observe({
    updateSliderInput(session, "snap", max = input$days,
                      value = min(8, input$days))
  })

  ## ---- tab 1 -------------------------------------------------------
  output$p_curve <- renderPlot({
    a <- seq(0, 120, 1)
    d <- do.call(rbind, lapply(names(VARIANTS), function(v) {
      p <- VARIANTS[[v]]
      data.frame(age = a, E = p$E0 * exp(-log(2) * a / p$TAU),
                 variant = p$label, sel = (v == input$variant))
    }))
    ggplot(d, aes(age, E, colour = variant, linewidth = sel)) +
      geom_line() +
      scale_linewidth_manual(values = c(`FALSE` = 0.5, `TRUE` = 1.8),
                             guide = "none") +
      scale_y_log10(labels = scales::percent) +
      labs(x = "적혈구 연령 (일)", y = "G6PD 활성 (정상 젊은 세포 대비, 로그)",
           colour = NULL,
           title = "변이 = 이 곡선의 기울기") +
      THEME
  })

  output$p_astar <- renderPlot({
    p  <- patient()
    ox <- 10^seq(-3, 2, length.out = 400)
    ## the same closed form the model reports as ASTAR, offset included
    oxe <- ox - 0.5 * 0.10
    astar <- ifelse(oxe <= 0, 120,
                    pmin(120, pmax(0, (p$TAU / log(2)) *
                                     log(p$E0 * 60 / pmax(oxe, 1e-12)))))
    ggplot(data.frame(ox, astar), aes(ox, astar)) +
      geom_line(linewidth = 1.2, colour = "#a02c2c") +
      geom_hline(yintercept = 120, linetype = 3) +
      annotate("text", x = 0.0015, y = 116, hjust = 0, size = 3.4,
               label = "위험한 세포 없음") +
      scale_x_log10() +
      labs(x = "산화제 부하 OX (mmol H2O2 / L RBC / 일, 로그)",
           y = "임계 연령 a* (일)",
           title = "a*보다 늙은 세포는 용혈한다") +
      THEME
  })

  output$t_variants <- renderTable(g6pd_age_table(), digits = 3)

  ## ---- tab 2 -------------------------------------------------------
  output$p_pk <- renderPlot({
    d <- sim() %>% filter(time >= 0) %>%
      select(time, CPQ_o, CTQ_o, CDP_o) %>%
      pivot_longer(-time) %>% filter(value > 1e-8)
    if (!nrow(d)) return(ggplot() + THEME + labs(title = "약물 없음"))
    ggplot(d, aes(time, value, colour = name)) + geom_line(linewidth = 0.9) +
      scale_colour_manual(values = c(CPQ_o = "#2f6fb5", CTQ_o = "#a02c2c",
                                     CDP_o = "#2c7a2c"),
        labels = c(CPQ_o = "프리마퀸", CTQ_o = "타페노퀸", CDP_o = "답손")) +
      labs(x = "시간 (일)", y = "혈장 농도 (mg/L)", colour = NULL,
           title = "약동학") + THEME
  })

  output$p_ox <- renderPlot({
    ggplot(filter(sim(), time >= 0), aes(time, OXFLUX)) +
      geom_line(colour = "#b5651d", linewidth = 1) +
      labs(x = "시간 (일)", y = "OX (mmol H2O2 / L RBC / 일)",
           title = "총 산화제 부하 — 모든 유발인자가 만나는 한 지점") + THEME
  })

  output$p_astar_t <- renderPlot({
    ggplot(filter(sim(), time >= 0), aes(time, ASTAR)) +
      geom_line(colour = "#a02c2c", linewidth = 1.1) +
      geom_hline(yintercept = 120, linetype = 3) +
      labs(x = "시간 (일)", y = "임계 연령 a* (일)",
           title = "a*가 연령 분포를 얼마나 깊이 쓸고 지나가는가") +
      ylim(0, 125) + THEME
  })

  ## ---- tab 3 -------------------------------------------------------
  age_long <- reactive({
    s <- sim()
    rn <- paste0("R", 1:12)
    s %>% select(time, all_of(rn)) %>%
      pivot_longer(-time, names_to = "bin", values_to = "cells") %>%
      mutate(bin = as.integer(sub("R", "", bin)),
             age = 10 * (bin - 0.5))
  })

  output$p_hist <- renderPlot({
    s  <- sim()
    tt <- s$time[which.min(abs(s$time - input$snap))]
    d  <- filter(age_long(), abs(time - tt) < 1e-6)
    as <- s$ASTAR[which.min(abs(s$time - tt))]
    ggplot(d, aes(age, cells)) +
      geom_col(fill = "#a9cdee", colour = "#25588c", width = 9) +
      geom_vline(xintercept = as, colour = "#a02c2c", linewidth = 1.3) +
      annotate("text", x = as, y = Inf, vjust = 1.6, hjust = -0.08,
               colour = "#a02c2c", size = 4.2,
               label = sprintf("a* = %.0f일", as)) +
      labs(x = "적혈구 연령 (일)", y = "세포 수 (10^12/L)",
           title = sprintf("연령 히스토그램, %.0f일째", tt)) + THEME
  })

  output$p_heat <- renderPlot({
    ggplot(filter(age_long(), time >= 0), aes(time, age, fill = cells)) +
      geom_raster() +
      geom_line(data = filter(sim(), time >= 0),
                aes(time, ASTAR), inherit.aes = FALSE,
                colour = "white", linewidth = 1.1) +
      scale_fill_viridis_c(option = "mako", direction = -1) +
      labs(x = "시간 (일)", y = "적혈구 연령 (일)", fill = "10^12/L",
           title = "흰 선이 a* — 그 위쪽이 비워지는 것이 용혈이다") + THEME
  })

  ## ---- tab 4 -------------------------------------------------------
  output$p_hb <- renderPlot({
    s <- sim(); b <- s$HB[which.min(abs(s$time))]
    ggplot(s, aes(time, HB)) +
      geom_vline(xintercept = 0, linetype = 3) +
      geom_hline(yintercept = b * 0.75, linetype = 2, colour = "#a02c2c") +
      annotate("text", x = max(s$time), y = b * 0.75, vjust = -0.5, hjust = 1,
               size = 3.3, colour = "#a02c2c", label = "기저 대비 -25%") +
      geom_line(linewidth = 1.1, colour = "#25588c") +
      labs(x = "시간 (일)", y = "혈색소 (g/dL)", title = "혈색소") + THEME
  })

  output$p_ret <- renderPlot({
    ggplot(sim(), aes(time, RETPCT)) +
      geom_line(linewidth = 1.1, colour = "#2c7a2c") +
      labs(x = "시간 (일)", y = "망상적혈구 (%)",
           title = "구조 요청 — 도착까지 4–7일") + THEME
  })

  output$p_epo <- renderPlot({
    ggplot(sim(), aes(time, EPOOUT)) +
      geom_line(linewidth = 1, colour = "#6b4a99") +
      scale_y_log10() +
      labs(x = "시간 (일)", y = "EPO (mIU/mL, 로그)", title = "에리스로포이에틴") +
      THEME
  })

  output$p_lys <- renderPlot({
    ggplot(filter(sim(), time >= 0), aes(time, ATRISK)) +
      geom_area(fill = "#f2b3ae", colour = "#963029") +
      labs(x = "시간 (일)", y = "a*보다 늙은 적혈구 비율 (%)",
           title = "지금 위험에 놓인 적혈구 질량") + THEME
  })

  output$t_summary <- renderTable(summarise_g6pd(sim()))

  ## ---- tab 5 -------------------------------------------------------
  cmp <- eventReactive(input$go_cmp, {
    shiny::req(length(input$cmp) > 0)
    withProgress(message = "변이별 비교 중...", {
      do.call(rbind, lapply(input$cmp, function(v) {
        p <- g6pd_patient(v, wt = input$wt, ugt = input$ugt,
                          cyp2d6 = input$cyp2d6, nat2 = input$nat2,
                          spleen = if (input$splenx) 0 else 1,
                          marrow = input$marrow)
        if (isTRUE(input$infect)) p$INFON <- 1
        run_g6pd(patient = p, rx = regimen(), days = input$days,
                 extra = extra(), delta = 0.05,
                 name = VARIANTS[[v]]$label)
      }))
    })
  })

  output$p_cmp <- renderPlot({
    d <- cmp()
    d %>% group_by(scenario) %>%
      mutate(pct = 100 * HB / HB[which.min(abs(time))]) %>%
      ggplot(aes(time, pct, colour = scenario)) +
      geom_hline(yintercept = 100, linetype = 3) +
      geom_hline(yintercept = 75, linetype = 2, colour = "#a02c2c") +
      geom_line(linewidth = 1.1) +
      labs(x = "시간 (일)", y = "혈색소 (기저 대비 %)", colour = NULL,
           title = "같은 약, 같은 용량 — 다른 붕괴 상수") + THEME
  })

  output$t_cmp <- renderTable({
    d <- cmp()
    do.call(rbind, lapply(split(d, d$scenario), summarise_g6pd))
  })

  ## ---- tab 6 -------------------------------------------------------
  output$p_met <- renderPlot({
    ggplot(sim(), aes(time, METPCT)) +
      geom_hline(yintercept = c(15, 30), linetype = 2, colour = "#a02c2c") +
      geom_line(linewidth = 1, colour = "#2f4a91") +
      labs(x = "시간 (일)", y = "메트헤모글로빈 (%)",
           title = "MetHb — 15% 청색증, 30% 증상") + THEME
  })

  output$p_bili <- renderPlot({
    ggplot(sim(), aes(time, TSB)) +
      geom_line(linewidth = 1, colour = "#8f7212") +
      labs(x = "시간 (일)", y = "총 빌리루빈 (mg/dL)", title = "빌리루빈") + THEME
  })

  output$p_hapto <- renderPlot({
    s <- sim()
    ggplot(s, aes(time)) +
      geom_line(aes(y = HPOUT, colour = "합토글로빈 (g/L)"), linewidth = 1) +
      geom_line(aes(y = FHBOUT * 5, colour = "유리 Hb x5 (g/L)"), linewidth = 1) +
      scale_colour_manual(values = c("#2b6f68", "#963029")) +
      labs(x = "시간 (일)", y = "농도", colour = NULL,
           title = "혈관내 용혈 — 합토글로빈이 먼저 0이 된다") + THEME
  })

  output$p_aki <- renderPlot({
    ggplot(sim(), aes(time, CREA)) +
      geom_line(linewidth = 1, colour = "#963029") +
      labs(x = "시간 (일)", y = "혈청 크레아티닌 (mg/dL)",
           title = "여과된 유리 헤모글로빈에 의한 세뇨관 손상") + THEME
  })

  ## ---- tab 7 -------------------------------------------------------
  output$p_assay <- renderPlot({
    s <- sim()
    ggplot(s, aes(time)) +
      annotate("rect", xmin = -Inf, xmax = Inf, ymin = 70, ymax = Inf,
               fill = "#e8f6e8") +
      annotate("rect", xmin = -Inf, xmax = Inf, ymin = 30, ymax = 70,
               fill = "#fdf6ee") +
      annotate("rect", xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = 30,
               fill = "#fbe0e0") +
      geom_hline(yintercept = c(30, 70), linetype = 2) +
      geom_line(aes(y = G6PDPCT), linewidth = 1.3, colour = "#8c4c11") +
      geom_line(aes(y = 100 * HB / max(HB)), linewidth = 0.8,
                linetype = 2, colour = "#25588c") +
      annotate("text", x = max(s$time), y = 74, hjust = 1, size = 3.4,
               label = "70% — 타페노퀸 허용선") +
      annotate("text", x = max(s$time), y = 33, hjust = 1, size = 3.4,
               label = "30% — 매일 프리마퀸 허용선") +
      labs(x = "시간 (일)",
           y = "검사가 보고할 G6PD 활성 (정상 %) · 점선은 혈색소(상대값)",
           title = "실선이 검사 수치 — 용혈 중에 올라간다") + THEME
  })

  output$t_assay <- renderTable({
    s <- filter(sim(), time >= 0)
    nadir_t <- s$time[which.min(s$HB)]
    j14 <- which.min(abs(s$time - (nadir_t + 14)))
    ## column names are given as STRINGS, not identifiers: R will not parse a
    ## Hangul identifier under a C locale, which is what a bare container has
    data.frame(
      "시점"        = c("노출 전", "혈색소 최저점", "최저점 +14일", "관찰 종료"),
      "일"          = round(c(0, nadir_t, nadir_t + 14, max(s$time)), 1),
      "검사값 (%)"  = round(c(s$G6PDPCT[1], s$G6PDPCT[which.min(s$HB)],
                              s$G6PDPCT[j14], tail(s$G6PDPCT, 1)), 1),
      "혈색소"      = round(c(s$HB[1], min(s$HB), s$HB[j14],
                              tail(s$HB, 1)), 2),
      check.names = FALSE)
  })

  ## ---- tab 8 -------------------------------------------------------
  neo <- eventReactive(input$go_neo, {
    withProgress(message = "신생아 4군 시뮬레이션...", {
      grid <- list(c("normal", "6/6", "둘 다 없음"),
                   c("Mediterranean", "6/6", "G6PD 결핍만"),
                   c("normal", "7/7", "UGT1A1 7/7만"),
                   c("Mediterranean", "7/7", "둘 다"))
      do.call(rbind, lapply(grid, function(g)
        run_g6pd(patient = g6pd_patient(g[1], wt = 3, ugt = g[2],
                                        neonate = TRUE, pna = 0.5),
                 rx = NULL, days = 12, lead_in = 0, delta = 0.02,
                 name = g[3])))
    })
  })

  output$p_neo <- renderPlot({
    ggplot(neo(), aes(time, TSB, colour = scenario)) +
      geom_hline(yintercept = 15, linetype = 2, colour = "#c9a227") +
      geom_hline(yintercept = 20, linetype = 2, colour = "#a02c2c") +
      annotate("text", x = 0.2, y = 15.6, hjust = 0, size = 3.4,
               colour = "#8f7212", label = "광선치료 문턱 (대략)") +
      annotate("text", x = 0.2, y = 20.6, hjust = 0, size = 3.4,
               colour = "#a02c2c", label = "교환수혈 문턱 (대략)") +
      geom_line(linewidth = 1.2) +
      labs(x = "생후 일수", y = "총 혈청 빌리루빈 (mg/dL)", colour = NULL,
           title = "G6PD × UGT1A1 — 더해지는 것이 아니라 곱해진다") + THEME
  })

  output$t_neo <- renderTable({
    d <- neo()
    do.call(rbind, lapply(split(d, d$scenario), function(x)
      data.frame("군"                  = x$scenario[1],
                 "최고 TSB (mg/dL)"    = round(max(x$TSB), 2),
                 "최고일"              = round(x$time[which.max(x$TSB)], 1),
                 "최고 자유빌리루빈 (nM)" = round(max(x$BFREE), 1),
                 "최저 혈색소 (g/dL)"  = round(min(x$HB), 2),
                 check.names = FALSE)))
  })
}

shinyApp(ui, server)
