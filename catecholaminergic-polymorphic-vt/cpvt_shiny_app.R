## =============================================================================
##  cpvt_shiny_app.R
##  카테콜아민성 다형성 심실빈맥 (CPVT) QSP 인터랙티브 대시보드
##
##  10 tabs. The app is organised around the model's central identity rather
##  than around drug classes, because the whole clinical argument of CPVT is
##  visible only when the two crossings are shown separately:
##
##      PVCR = W( CAJSRF / CTHR )  x  P( VDAD / VREQ )
##
##  Tab 1  환자 프로파일    — genotype, triggers, electrolytes, devices
##  Tab 2  두 교차 (핵심)   — load vs threshold AND DAD vs requirement, side by side
##  Tab 3  약물 PK          — plasma concentrations over the dosing interval
##  Tab 4  β1 차단 최저치   — the trough argument: nadolol vs metoprolol
##  Tab 5  Ca 순환          — SERCA/NCX/RyR2 fluxes and SR content
##  Tab 6  부정맥 지표      — PVC rate, VT burden, VE score, onset heart rate
##  Tab 7  임상 엔드포인트  — cumulative hazard, syncope, ICD shocks
##  Tab 8  시나리오 비교    — run several arms and rank them
##  Tab 9  복약 순응도      — what one missed dose costs, per drug
##  Tab 10 안전성/바이오마커— QRS width, fatigue, CaMKII, hypokalaemia
##
##  Run with:  shiny::runApp("cpvt_shiny_app.R")
##  Requires cpvt_mrgsolve_model.R in the same directory.
## =============================================================================

library(shiny)
library(mrgsolve)
library(dplyr)
library(ggplot2)
library(tidyr)
library(gridExtra)

source("cpvt_mrgsolve_model.R")

PAL <- c("#C0392B", "#1976D2", "#00838F", "#F9A825", "#7E57C2",
         "#2E7D32", "#D81B60", "#5D4037")
th <- theme_bw(base_size = 12) +
  theme(strip.background = element_rect(fill = "#F0F0F0"),
        legend.position = "bottom", plot.title = element_text(face = "bold"))

## ---------------------------------------------------------------------------
##  UI
## ---------------------------------------------------------------------------
ui <- fluidPage(
  titlePanel("CPVT QSP 시뮬레이터 — 두 교차의 곱으로서의 촉발 박동"),
  tags$p(style = "color:#555;margin-top:-8px",
    "촉발 박동 = W(SR 부하 / SOICR 문턱) × P(DAD 진폭 / 필요 탈분극량). ",
    "β차단제는 첫 번째 교차만, 플레카이니드는 주로 두 번째 교차를 움직인다."),
  sidebarLayout(
    sidebarPanel(
      width = 3,
      h4("환자 · 유전형"),
      selectInput("geno", "유전형",
        choices = c("야생형 (WT)" = "WT",
                    "CPVT1 전형 (RyR2)" = "CPVT1",
                    "CPVT1 중증" = "CPVT1s",
                    "CPVT1 경증 침투" = "CPVT1m",
                    "CPVT2 (CASQ2)" = "CPVT2"),
        selected = "CPVT1"),
      sliderInput("ko", "혈청 K+ (mmol/L)", 2.5, 5.5, 4.2, 0.1),
      checkboxInput("lcsd", "좌심장 교감신경 절제술 (LCSD)", FALSE),
      checkboxInput("icd", "ICD 삽입", FALSE),
      checkboxInput("shkfb", "쇼크 → 카테콜아민 되먹임 포함", TRUE),

      hr(), h4("부하 프로토콜"),
      radioButtons("prot", NULL,
        c("운동부하 검사" = "ex", "에피네프린 부하" = "epi",
          "안정 시 놀람 자극" = "st"), selected = "ex"),
      conditionalPanel("input.prot == 'ex'",
        sliderInput("peak", "최대 운동 강도 (1.0 = 최대 운동)", 0.2, 1.0, 1.0, 0.05)),

      hr(), h4("β차단제"),
      selectInput("bb", "약물",
        c("없음" = "none", "나돌롤" = "nad", "메토프롤롤" = "met",
          "비소프롤롤" = "bis"), selected = "nad"),
      conditionalPanel("input.bb != 'none'",
        numericInput("bbdose", "1회 용량 (mg)", 80, 5, 400, 5),
        selectInput("bbtau", "투여 간격", c("1일 1회" = 24, "1일 2회" = 12),
                    selected = 24),
        sliderInput("bboff", "마지막 복용 후 경과 시간 (h)", 0.5, 72, 24, 0.5),
        helpText("이 슬라이더가 이 모델의 핵심입니다. 최고 농도가 아니라 ",
                 "최저 농도에서 무슨 일이 일어나는지가 결과를 결정합니다.")),
      conditionalPanel("input.bb == 'met'",
        selectInput("cyp", "CYP2D6 표현형",
          c("정상 대사자 (EM)" = 1, "초고속 대사자 (UM)" = 2.5,
            "저대사자 (PM)" = 0.3), selected = 1)),

      hr(), h4("항부정맥제"),
      checkboxInput("useflec", "플레카이니드", FALSE),
      conditionalPanel("input.useflec == true",
        numericInput("fledose", "1회 용량 (mg, 1일 2회)", 100, 25, 200, 25),
        sliderInput("fleoff", "마지막 복용 후 경과 (h)", 0.5, 24, 12, 0.5),
        radioButtons("flearm", "작용 팔 (기전 해부용)",
          c("두 팔 모두" = "both", "RyR2 팔만" = "ryr", "Na 팔만" = "na"),
          selected = "both")),
      checkboxInput("useverap", "베라파밀 240 mg qd", FALSE),
      hr(),
      actionButton("go", "시뮬레이션 실행", class = "btn-primary btn-block")
    ),

    mainPanel(
      width = 9,
      tabsetPanel(
        id = "tabs",
        tabPanel("① 환자 프로파일", br(),
          fluidRow(column(6, h4("설정 요약"), tableOutput("profile")),
                   column(6, h4("정상상태 혈중 농도"), tableOutput("concs"))),
          hr(), h4("이 환자에서 두 교차의 안정 시 위치"),
          plotOutput("restbars", height = "260px"),
          helpText("안정 시 X = 부하/문턱 이 1보다 작으면 심전도·심초음파는 ",
                   "완전히 정상입니다. CPVT가 '구조적으로 정상인 심장의 병'인 ",
                   "이유가 이 막대 하나에 있습니다.")),

        tabPanel("② 두 교차 (핵심)", br(),
          plotOutput("crossings", height = "620px"),
          helpText("위: 교차 1 — 유리 접합부 SR Ca(부하)와 SOICR 문턱. 두 선이 ",
                   "교차하는 순간부터만 파동이 발생합니다. 아래: 교차 2 — DAD ",
                   "진폭과 촉발에 필요한 탈분극량.")),

        tabPanel("③ 약물 PK", br(),
          plotOutput("pkplot", height = "420px"),
          hr(), h4("투여 간격 전체에 걸친 경쟁 이동계수 1 + C/EC50"),
          plotOutput("shiftcurve", height = "300px")),

        tabPanel("④ β1 차단 최저치", br(),
          h4("동일한 수용체, 동일한 경쟁식 — 차이는 반감기에서만 나온다"),
          tableOutput("shifttab"),
          plotOutput("troughplot", height = "400px"),
          helpText("메토프롤롤은 배용량의 최고 농도에서 겨우 나돌롤의 최저 ",
                   "농도에 해당하는 차단을 달성합니다. CPVT 사건은 예고 없이 ",
                   "일어나므로 재판정에 서는 시간은 최저 농도 시점입니다.")),

        tabPanel("⑤ Ca 순환", br(), plotOutput("caplot", height = "560px")),

        tabPanel("⑥ 부정맥 지표", br(),
          fluidRow(column(4, h4("운동부하 심실부정맥 점수"), tableOutput("vetab")),
                   column(8, plotOutput("arrplot", height = "380px"))),
          hr(), plotOutput("onsetplot", height = "280px")),

        tabPanel("⑦ 임상 엔드포인트", br(),
          plotOutput("endplot", height = "480px"),
          hr(), tableOutput("endtab")),

        tabPanel("⑧ 시나리오 비교", br(),
          checkboxGroupInput("arms", "비교할 치료군",
            choices = c("무치료" = "none",
                        "나돌롤 최저 (24 h)" = "nad_tr",
                        "나돌롤 최고 (2 h)" = "nad_pk",
                        "메토프롤롤 100 bid 최저" = "met_tr",
                        "메토프롤롤 200 bid 최고" = "met2_pk",
                        "비소프롤롤 최저" = "bis_tr",
                        "나돌롤 + 플레카이니드" = "nad_fle",
                        "플레카이니드 단독" = "fle_only",
                        "나돌롤 + LCSD" = "nad_lcsd",
                        "LCSD 단독" = "lcsd"),
            selected = c("none", "nad_tr", "met_tr", "nad_fle"), inline = TRUE),
          actionButton("gocmp", "비교 실행", class = "btn-primary"),
          br(), br(), plotOutput("cmpplot", height = "400px"),
          hr(), tableOutput("cmptab")),

        tabPanel("⑨ 복약 순응도", br(),
          h4("한 번 빠뜨린 복용의 대가는 약마다 다르다"),
          actionButton("gomiss", "누락 시나리오 실행", class = "btn-primary"),
          br(), br(), plotOutput("missplot", height = "400px"),
          hr(), tableOutput("misstab"),
          helpText("LCSD를 병용한 군을 함께 보십시오. LCSD의 임상적 가치는 ",
                   "최대 효과가 아니라 복약 순응도와 무관하게 유지되는 ",
                   "부분에 있습니다.")),

        tabPanel("⑩ 안전성 · 바이오마커", br(),
          plotOutput("safeplot", height = "520px"),
          hr(), h4("혈청 칼륨과 재분극 예비력"),
          plotOutput("kplot", height = "280px"))
      )
    )
  )
)

## ---------------------------------------------------------------------------
##  SERVER
## ---------------------------------------------------------------------------
server <- function(input, output, session) {

  drug_list <- reactive({
    d <- list()
    if (input$bb != "none")
      d <- c(d, list(list(drug = input$bb, dose = input$bbdose,
                          tau = as.numeric(input$bbtau), offset = input$bboff)))
    if (isTRUE(input$useflec))
      d <- c(d, list(list(drug = "fle", dose = input$fledose, tau = 12,
                          offset = input$fleoff)))
    if (isTRUE(input$useverap))
      d <- c(d, list(list(drug = "ver", dose = 240, tau = 24, offset = 4)))
    d
  })

  pmod_list <- reactive({
    p <- list(KO = input$ko,
              FLCSD = as.numeric(isTRUE(input$lcsd)),
              ICD = as.numeric(isTRUE(input$icd)))
    if (!isTRUE(input$shkfb)) p$ESHK <- 0
    if (input$bb == "met") p$CYP2D6 <- as.numeric(input$cyp)
    if (isTRUE(input$useflec)) {
      if (input$flearm == "ryr") p$KONNA <- 0
      if (input$flearm == "na")  p$EFLR  <- 0
    }
    p
  })

  protocol <- reactive({
    switch(input$prot,
           ex  = exercise_data(peak = input$peak),
           epi = epi_data(),
           st  = startle_data())
  })

  sim <- eventReactive(input$go, {
    run_cpvt(genotype = input$geno, drugs = drug_list(),
             pmod = pmod_list(), protocol = protocol(),
             label = input$geno)
  }, ignoreNULL = FALSE)

  ## ---- Tab 1 -------------------------------------------------------------
  output$profile <- renderTable({
    data.frame(항목 = c("유전형", "혈청 K+", "LCSD", "ICD", "부하 프로토콜",
                        "β차단제", "플레카이니드"),
               값 = c(input$geno, paste(input$ko, "mmol/L"),
                      ifelse(input$lcsd, "시행", "미시행"),
                      ifelse(input$icd, "삽입", "없음"),
                      switch(input$prot, ex = "운동부하", epi = "에피네프린",
                             st = "놀람 자극"),
                      ifelse(input$bb == "none", "없음",
                             sprintf("%s %g mg q%sh, 복용 후 %g h",
                                     input$bb, input$bbdose, input$bbtau,
                                     input$bboff)),
                      ifelse(input$useflec,
                             sprintf("%g mg bid, 복용 후 %g h",
                                     input$fledose, input$fleoff), "없음")))
  })

  output$concs <- renderTable({
    cs <- attr(sim(), "conc")
    if (is.null(cs) || !length(cs)) return(data.frame(메시지 = "투여 약물 없음"))
    data.frame(약물 = names(cs), `농도_ng_mL` = round(as.numeric(cs), 1))
  })

  output$restbars <- renderPlot({
    s <- sim()
    d <- data.frame(
      항목 = c("유리 SR Ca (부하)", "SOICR 문턱", "DAD 진폭", "필요 탈분극량"),
      값 = c(s$CAJSRFo[1], s$CTHRo[1], max(s$VDAD[1], 0.1), 22),
      군 = c("교차 1", "교차 1", "교차 2", "교차 2"))
    ggplot(d, aes(항목, 값, fill = 군)) +
      geom_col(width = 0.6) + facet_wrap(~군, scales = "free") +
      scale_fill_manual(values = PAL[c(1, 2)]) +
      labs(x = NULL, y = NULL, title = "안정 시 두 교차의 위치") + th
  })

  ## ---- Tab 2: the two crossings -----------------------------------------
  output$crossings <- renderPlot({
    s <- sim()
    p1 <- s %>% select(time, `유리 SR Ca` = CAJSRFo, `SOICR 문턱` = CTHRo) %>%
      pivot_longer(-time) %>%
      ggplot(aes(time, value, colour = name)) +
      geom_line(linewidth = 1) + scale_colour_manual(values = PAL[c(1, 2)]) +
      labs(x = NULL, y = "µM", colour = NULL,
           title = "교차 1 — 부하가 문턱을 넘는 순간부터만 파동이 생긴다") + th
    p2 <- ggplot(s, aes(time, Xo)) +
      geom_hline(yintercept = 1, linetype = 2, colour = "#C0392B") +
      geom_line(linewidth = 1, colour = PAL[5]) +
      labs(x = NULL, y = "X = 부하/문턱", title = NULL) + th
    p3 <- s %>% mutate(`필요 탈분극량` = 22 / NAAV^1.5) %>%
      select(time, `DAD 진폭` = VDAD, `필요 탈분극량`) %>%
      pivot_longer(-time) %>%
      ggplot(aes(time, value, colour = name)) +
      geom_line(linewidth = 1) + scale_colour_manual(values = PAL[c(1, 3)]) +
      labs(x = NULL, y = "mV", colour = NULL,
           title = "교차 2 — DAD가 필요 탈분극량을 넘어야 실제 박동이 된다") + th
    p4 <- ggplot(s, aes(time, PVCR)) +
      geom_line(linewidth = 1, colour = PAL[1]) +
      labs(x = "시간 (min)", y = "촉발 박동 (1/min)",
           title = "두 교차의 곱") + th
    gridExtra::grid.arrange(p1, p2, p3, p4, ncol = 1, heights = c(3, 2, 3, 2))
  })

  ## ---- Tab 3: PK ---------------------------------------------------------
  output$pkplot <- renderPlot({
    s <- sim()
    s %>% select(time, 나돌롤 = CNADo, 메토프롤롤 = CMETo,
                 비소프롤롤 = CBISo, 플레카이니드 = CFLEo, 베라파밀 = CVERo) %>%
      pivot_longer(-time) %>% filter(value > 1e-6) %>%
      ggplot(aes(time, value, colour = name)) +
      geom_line(linewidth = 1) + scale_colour_manual(values = PAL) +
      labs(x = "시간 (min)", y = "혈중 농도 (ng/mL)", colour = NULL,
           title = "시뮬레이션 구간 중 혈중 농도") + th
  })

  output$shiftcurve <- renderPlot({
    grid <- expand.grid(off = seq(0.5, 48, by = 0.5),
                        reg = c("나돌롤 80 qd", "메토프롤롤 100 bid",
                                "메토프롤롤 200 bid", "비소프롤롤 10 qd"),
                        stringsAsFactors = FALSE)
    spec <- list("나돌롤 80 qd"       = c("nad", 80, 24, 20),
                 "메토프롤롤 100 bid" = c("met", 100, 12, 45),
                 "메토프롤롤 200 bid" = c("met", 200, 12, 45),
                 "비소프롤롤 10 qd"   = c("bis", 10, 24, 12))
    grid$shift <- mapply(function(o, r) {
      sp <- spec[[r]]
      cc <- attr(pk_steady(sp[1], as.numeric(sp[2]), as.numeric(sp[3]), o), "conc")
      1 + cc / as.numeric(sp[4])
    }, grid$off, grid$reg)
    ggplot(grid, aes(off, shift, colour = reg)) +
      geom_line(linewidth = 1) +
      scale_colour_manual(values = PAL) +
      labs(x = "마지막 복용 후 경과 시간 (h)", y = "경쟁 이동계수 1 + C/EC50",
           colour = NULL,
           title = "복용 간격 전체에 걸친 β1 차단 — 최저점이 결과를 결정한다") + th
  })

  ## ---- Tab 4: the trough argument ---------------------------------------
  output$shifttab <- renderTable({ shift_table() }, digits = 2)

  output$troughplot <- renderPlot({
    arms <- list(
      list(l = "무치료", d = list(), p = list()),
      list(l = "나돌롤 최저 (24 h)",
           d = list(list(drug = "nad", dose = 80, tau = 24, offset = 24))),
      list(l = "메토프롤롤 100 최저 (12 h)",
           d = list(list(drug = "met", dose = 100, tau = 12, offset = 12))),
      list(l = "메토프롤롤 100 최고 (2 h)",
           d = list(list(drug = "met", dose = 100, tau = 12, offset = 2))),
      list(l = "메토프롤롤 200 최고 (2 h)",
           d = list(list(drug = "met", dose = 200, tau = 12, offset = 2))))
    bind_rows(lapply(arms, function(a) {
      s <- run_cpvt(input$geno, a$d, list(KO = input$ko),
                    exercise_data(), label = a$l)
      data.frame(arm = a$l, time = s$time, PVCR = s$PVCR)
    })) %>%
      ggplot(aes(time, PVCR, colour = arm)) +
      geom_line(linewidth = 0.9) + scale_colour_manual(values = PAL) +
      labs(x = "시간 (min)", y = "촉발 박동 (1/min)", colour = NULL,
           title = "최저 농도에서의 β1 차단이 운동 중 부정맥을 결정한다") + th
  })

  ## ---- Tab 5: Ca cycling -------------------------------------------------
  output$caplot <- renderPlot({
    sim() %>%
      select(time, `세포질 Ca (µM)` = CAI, `NSR Ca (µM)` = CANSR,
             `JSR 총 Ca (µM)` = CAJSRT, `유리 JSR Ca (µM)` = CAJSRFo,
             `P-PLB 분율` = PPLB, `P-RyR2 분율` = PRYR,
             `RyR2 불응 분율` = RREF, `미토콘드리아 Ca` = CAMT) %>%
      pivot_longer(-time) %>%
      ggplot(aes(time, value)) + geom_line(colour = PAL[1], linewidth = 0.8) +
      facet_wrap(~name, scales = "free_y", ncol = 2) +
      labs(x = "시간 (min)", y = NULL, title = "Ca 순환 상태변수") + th
  })

  ## ---- Tab 6: arrhythmia -------------------------------------------------
  output$vetab <- renderTable({
    e <- endpoints(sim())
    data.frame(지표 = c("VE 점수 (0-4)", "최고 심박수", "PVC 출현 심박수",
                        "VT 출현 심박수", "최고 PVC 발생률", "총 PVC 수",
                        "VT 지속시간 (min)", "X 최대값"),
               값 = c(e$VE_score, round(e$HR_peak), round(e$HR_at_PVC),
                      round(e$HR_at_VT), round(e$PVC_peak, 1),
                      round(e$PVC_count), round(e$VT_minutes, 2),
                      round(e$X_peak, 3)))
  })

  output$arrplot <- renderPlot({
    sim() %>% select(time, `심박수` = HR, `촉발 박동 (1/min)` = PVCR,
                     `VT 부담` = VTB, `X = 부하/문턱` = Xo) %>%
      pivot_longer(-time) %>%
      ggplot(aes(time, value)) + geom_line(colour = PAL[1], linewidth = 0.9) +
      facet_wrap(~name, scales = "free_y", ncol = 2) +
      labs(x = "시간 (min)", y = NULL) + th
  })

  output$onsetplot <- renderPlot({
    s <- sim()
    ggplot(s, aes(HR, PVCR)) +
      geom_path(colour = PAL[1], linewidth = 0.9) +
      labs(x = "심박수 (bpm)", y = "촉발 박동 (1/min)",
           title = "심박수-부정맥 이력 (운동 중 상승과 회복이 겹치지 않는다)") + th
  })

  ## ---- Tab 7: endpoints --------------------------------------------------
  output$endplot <- renderPlot({
    sim() %>% select(time, `누적 부정맥 위험` = HZRD, `누적 실신 위험` = SYNCH,
                     `누적 PVC 수` = NPVC, `ICD 쇼크 수` = NSHOCK) %>%
      pivot_longer(-time) %>%
      ggplot(aes(time, value)) + geom_line(colour = PAL[2], linewidth = 0.9) +
      facet_wrap(~name, scales = "free_y", ncol = 2) +
      labs(x = "시간 (min)", y = NULL) + th
  })
  output$endtab <- renderTable({ endpoints(sim()) }, digits = 3)

  ## ---- Tab 8: scenario comparison ---------------------------------------
  ARMS <- list(
    none     = list(l = "무치료", d = list(), p = list()),
    nad_tr   = list(l = "나돌롤 최저 (24 h)",
                    d = list(list(drug = "nad", dose = 80, tau = 24, offset = 24))),
    nad_pk   = list(l = "나돌롤 최고 (2 h)",
                    d = list(list(drug = "nad", dose = 80, tau = 24, offset = 2))),
    met_tr   = list(l = "메토프롤롤 100 bid 최저",
                    d = list(list(drug = "met", dose = 100, tau = 12, offset = 12))),
    met2_pk  = list(l = "메토프롤롤 200 bid 최고",
                    d = list(list(drug = "met", dose = 200, tau = 12, offset = 2))),
    bis_tr   = list(l = "비소프롤롤 최저",
                    d = list(list(drug = "bis", dose = 10, tau = 24, offset = 24))),
    nad_fle  = list(l = "나돌롤 + 플레카이니드",
                    d = list(list(drug = "nad", dose = 80, tau = 24, offset = 24),
                             list(drug = "fle", dose = 100, tau = 12, offset = 12))),
    fle_only = list(l = "플레카이니드 단독",
                    d = list(list(drug = "fle", dose = 100, tau = 12, offset = 12))),
    nad_lcsd = list(l = "나돌롤 + LCSD",
                    d = list(list(drug = "nad", dose = 80, tau = 24, offset = 24)),
                    p = list(FLCSD = 1)),
    lcsd     = list(l = "LCSD 단독", d = list(), p = list(FLCSD = 1)))

  cmp <- eventReactive(input$gocmp, {
    sel <- if (length(input$arms)) input$arms else c("none", "nad_tr")
    bind_rows(lapply(sel, function(k) {
      a <- ARMS[[k]]
      pm <- c(list(KO = input$ko), if (is.null(a$p)) list() else a$p)
      s <- run_cpvt(input$geno, a$d, pm, exercise_data(peak = input$peak), a$l)
      cbind(endpoints(s), arm = a$l)
    }))
  })

  output$cmptab <- renderTable({
    cmp() %>% select(arm, VE_score, HR_peak, HR_at_PVC, PVC_peak, PVC_count,
                     VT_minutes, hazard, QRS_max) %>% arrange(PVC_count)
  }, digits = 3)

  output$cmpplot <- renderPlot({
    cmp() %>% mutate(arm = reorder(arm, PVC_count)) %>%
      ggplot(aes(arm, PVC_count, fill = arm)) +
      geom_col(width = 0.65, show.legend = FALSE) +
      scale_fill_manual(values = rep(PAL, 3)) + coord_flip() +
      labs(x = NULL, y = "운동부하 검사 중 총 촉발 박동 수",
           title = "치료군 순위 (낮을수록 좋음)") + th
  })

  ## ---- Tab 9: adherence --------------------------------------------------
  miss <- eventReactive(input$gomiss, {
    rows <- list(
      list(l = "나돌롤, 정시 (24 h)",  d = "nad", dose = 80, tau = 24, off = 24, lc = 0),
      list(l = "나돌롤, 1회 누락 (48 h)", d = "nad", dose = 80, tau = 24, off = 48, lc = 0),
      list(l = "나돌롤, 2회 누락 (72 h)", d = "nad", dose = 80, tau = 24, off = 72, lc = 0),
      list(l = "메토프롤롤, 정시 (12 h)", d = "met", dose = 100, tau = 12, off = 12, lc = 0),
      list(l = "메토프롤롤, 1회 누락 (24 h)", d = "met", dose = 100, tau = 12, off = 24, lc = 0),
      list(l = "비소프롤롤, 1회 누락 (48 h)", d = "bis", dose = 10, tau = 24, off = 48, lc = 0),
      list(l = "LCSD + 나돌롤 1회 누락", d = "nad", dose = 80, tau = 24, off = 48, lc = 1),
      list(l = "LCSD + 나돌롤 2회 누락", d = "nad", dose = 80, tau = 24, off = 72, lc = 1))
    bind_rows(lapply(rows, function(r) {
      s <- run_cpvt(input$geno,
                    list(list(drug = r$d, dose = r$dose, tau = r$tau, offset = r$off)),
                    list(KO = input$ko, FLCSD = r$lc), exercise_data(), r$l)
      cbind(endpoints(s), arm = r$l,
            conc = round(as.numeric(attr(s, "conc")[1]), 1))
    }))
  })

  output$misstab <- renderTable({
    miss() %>% select(arm, conc, shift_min, VE_score, PVC_count, VT_minutes)
  }, digits = 3)

  output$missplot <- renderPlot({
    miss() %>% mutate(arm = factor(arm, levels = arm)) %>%
      ggplot(aes(arm, PVC_count, fill = arm)) +
      geom_col(width = 0.65, show.legend = FALSE) +
      scale_fill_manual(values = rep(PAL, 2)) + coord_flip() +
      labs(x = NULL, y = "총 촉발 박동 수",
           title = "반감기는 순응도와 같은 축이다") + th
  })

  ## ---- Tab 10: safety ----------------------------------------------------
  output$safeplot <- renderPlot({
    sim() %>% select(time, `QRS 폭 (ms)` = QRSD, `Na 통로 가용도` = NAAV,
                     `피로 지수` = FATG, `CaMKII 활성` = CAMK,
                     `누적 Ca 과부하` = CAINT, `섬유화 지수` = FIBR) %>%
      pivot_longer(-time) %>%
      ggplot(aes(time, value)) + geom_line(colour = PAL[4], linewidth = 0.9) +
      facet_wrap(~name, scales = "free_y", ncol = 2) +
      labs(x = "시간 (min)", y = NULL,
           title = "안전성 지표 및 만성 되먹임") + th
  })

  output$kplot <- renderPlot({
    ks <- seq(2.8, 5.4, by = 0.2)
    bind_rows(lapply(ks, function(k) {
      s <- run_cpvt(input$geno, drug_list(),
                    c(pmod_list()[setdiff(names(pmod_list()), "KO")], list(KO = k)),
                    exercise_data(peak = input$peak))
      cbind(endpoints(s), K = k)
    })) %>%
      ggplot(aes(K, PVC_count)) +
      geom_line(colour = PAL[1], linewidth = 1) + geom_point(colour = PAL[1]) +
      labs(x = "혈청 K+ (mmol/L)", y = "총 촉발 박동 수",
           title = "저칼륨혈증은 I_K1 sink를 줄여 같은 파동을 더 큰 DAD로 만든다") + th
  })
}

shinyApp(ui, server)
