## =====================================================================
##  sal_shiny_app.R
##  Salicylate (Aspirin) Poisoning — interactive QSP dashboard
##  11 tabs · driven by sal_mrgsolve_model.R
##
##  살리실산(아스피린) 중독 — QSP 대시보드
## ---------------------------------------------------------------------
##  The app is built around ONE decomposition, which every tab returns to:
##
##      C_brain = C_total x fu x f_n(pH_plasma)/f_n(pH_brain)
##
##  Tab 3 ("두 배수") shows the three factors side by side over time, so the
##  user can watch the reported number and the brain concentration come
##  apart.  Tab 8 lets the user move the ventilator and see it happen.
##
##  RUN
##    R -e "shiny::runApp('sal_shiny_app.R', port = 8080)"
##  (sal_mrgsolve_model.R must be in the same directory; it is sourced on
##   startup, which compiles the model once.)
## =====================================================================

library(shiny)
library(mrgsolve)
library(dplyr)
library(ggplot2)
library(tidyr)

source("sal_mrgsolve_model.R")

THEME <- theme_bw(base_size = 12) +
  theme(panel.grid.minor = element_blank(),
        strip.background = element_rect(fill = "#ECEFF1", colour = NA),
        legend.position = "bottom")

PAL <- c(plasma = "#1565C0", free = "#7986CB", brain = "#00838F",
         tissue = "#8D6E63", target = "#C62828", alt = "#2E7D32")

## ---------------------------------------------------------------------
##  UI
## ---------------------------------------------------------------------
ui <- fluidPage(
  titlePanel("살리실산(아스피린) 중독 QSP 대시보드 — Salicylate Poisoning QSP Dashboard"),
  tags$p(style = "color:#555;margin-top:-8px",
         HTML("측정되는 값(총 혈장 살리실산)과 사람을 죽이는 값(뇌 세포 내 살리실산)은
               <b>유리분율</b>과 <b>pH 분배</b>라는 두 번의 곱셈으로 갈라진다.
               이 앱은 그 두 배수를 계속 함께 보여준다.")),

  sidebarLayout(
    sidebarPanel(
      width = 3,
      h4("환자 (Patient)"),
      numericInput("bw", "체중 Body weight (kg)", 70, 30, 150, 5),
      sliderInput("alb", "알부민 Albumin (g/L)", 15, 50, 40, 1),
      sliderInput("gfr", "GFR (mL/min)", 10, 140, 120, 5),
      checkboxInput("aged", "고령 (폐손상 감수성 ↑)", FALSE),

      h4("노출 (Exposure)"),
      selectInput("mode", "형태 Pattern",
                  c("급성 단회 Acute single" = "acute",
                    "장용정/결석 Enteric-coated / concretion" = "ec",
                    "만성 과다 Chronic overshoot" = "chronic")),
      conditionalPanel("input.mode != 'chronic'",
        sliderInput("dose", "아스피린 용량 Dose (g)", 1, 60, 30, 1)),
      conditionalPanel("input.mode == 'chronic'",
        sliderInput("cdose", "1회 용량 (mg), q6h", 325, 1950, 1300, 325),
        sliderInput("cdays", "투여 기간 (일)", 3, 21, 14, 1)),

      h4("치료 (Treatment)"),
      sliderInput("trx", "치료 시작 시각 (h)", 0, 24, 4, 1),
      checkboxInput("fluid", "수액 crystalloid", TRUE),
      checkboxInput("bic", "탄산수소나트륨 NaHCO3", TRUE),
      sliderInput("kcl", "칼륨 보충 K (mmol/h)", 0, 30, 10, 1),
      checkboxInput("dex", "포도당 dextrose (혈당 11 mmol/L)", FALSE),
      checkboxInput("ac", "활성탄 MDAC", FALSE),
      checkboxInput("hd", "혈액투석 haemodialysis", FALSE),
      conditionalPanel("input.hd == true",
        sliderInput("thd", "투석 시작 (h)", 2, 36, 6, 1)),

      h4("호흡 (Ventilation)"),
      radioButtons("vent", "",
                   c("자발호흡 Spontaneous" = "spont",
                     "삽관 · 통상 PaCO2 Intubated, conventional" = "conv",
                     "삽관 · 분당환기량 일치 Intubated, matched" = "match",
                     "진정만 Sedation only (45%)" = "sed")),
      conditionalPanel("input.vent != 'spont'",
        sliderInput("tvent", "삽관/진정 시각 (h)", 1, 24, 8, 1)),

      hr(),
      sliderInput("tmax", "시뮬레이션 기간 (h)", 12, 340, 72, 4),
      actionButton("go", "다시 계산 Run", class = "btn-primary")
    ),

    mainPanel(
      width = 9,
      tabsetPanel(
        id = "tabs",
        tabPanel("1. 환자 프로파일", br(),
                 fluidRow(column(6, h4("요약 Summary"), tableOutput("summary")),
                          column(6, h4("경보 Alerts"), uiOutput("alerts"))),
                 hr(), h4("주요 지표 한눈에"), plotOutput("overview", height = "420px")),

        tabPanel("2. 약동학 PK", br(),
                 p("총 혈장농도(측정되는 값), 유리농도, 조직농도, 그리고 뇌농도."),
                 plotOutput("pk", height = "380px"),
                 h4("제거 경로 Elimination pathways (mg/h)"),
                 plotOutput("clr", height = "300px")),

        tabPanel("3. 두 배수 (The two multipliers)", br(),
                 p(HTML("<b>C_brain = C_total × fu × Kp</b>. 세 항을 따로 그린다.
                         셋 중 두 개는 검사실에서 보이지 않는다.")),
                 plotOutput("mult", height = "420px"),
                 h4("총 혈장농도 대비 뇌농도 비 (brain : reported level)"),
                 plotOutput("ratio", height = "260px")),

        tabPanel("4. 산-염기 Acid-base", br(),
                 plotOutput("ab", height = "440px"),
                 h4("환기 구동 Ventilatory drive"),
                 plotOutput("vent", height = "260px")),

        tabPanel("5. 신장 · 전해질 Renal", br(),
                 plotOutput("renal", height = "440px"),
                 p(HTML("<b>요 pH가 오르지 않으면 알칼리화는 아무 일도 하지 않는다.</b>
                         그리고 요 pH는 혈장 중탄산이 신세뇨관 역치를 넘어야 오르며,
                         저칼륨혈증은 그 역치를 올린다."))),

        tabPanel("6. 임상 엔드포인트", br(),
                 plotOutput("endp", height = "440px"),
                 h4("누적 중추신경 노출 (duration, not just peak)"),
                 plotOutput("cnsi", height = "240px")),

        tabPanel("7. 시나리오 비교", br(),
                 p("모델에 내장된 14개 시나리오. 두 개를 골라 겹쳐 본다."),
                 fluidRow(
                   column(5, selectInput("sc1", "시나리오 A",
                                         setNames(names(SCEN),
                                                  sapply(SCEN, `[[`, "label")),
                                         selected = "S07")),
                   column(5, selectInput("sc2", "시나리오 B",
                                         setNames(names(SCEN),
                                                  sapply(SCEN, `[[`, "label")),
                                         selected = "S08")),
                   column(2, br(), actionButton("gosc", "비교 Compare"))),
                 plotOutput("scen", height = "480px"),
                 tableOutput("scentab")),

        tabPanel("8. 인공호흡기 실험", br(),
                 h4("스위치 하나만 바꾼다 — One switch, nothing else"),
                 p(HTML("과호흡하고 있는 환자를 삽관하고 <b>통상적인 분당환기량</b>으로
                         맞추면 PaCO2가 오르고 혈장 pH가 떨어진다. 세포 내 pH는 완충되어
                         거의 움직이지 않으므로 <b>분배계수가 커지고 약물이 뇌로 들어간다</b>.
                         같은 순간 검사실 수치는 <b>내려간다</b>.")),
                 plotOutput("ventexp", height = "460px"),
                 tableOutput("venttab")),

        tabPanel("9. Done 노모그램", br(),
                 p(HTML("노모그램은 <b>입력이 아니라 출력</b>으로 취급한다.
                         시뮬레이션에서 읽은 (시간, 혈장농도)를 노모그램에 넣어 나온 등급과,
                         모델이 계산한 뇌농도 등급을 나란히 비교한다.")),
                 plotOutput("nomo", height = "440px"),
                 tableOutput("nomotab")),

        tabPanel("10. 투석 결정 EXTRIP", br(),
                 uiOutput("extrip"),
                 plotOutput("hdplot", height = "400px")),

        tabPanel("11. 검사 패널 Biomarkers", br(),
                 plotOutput("labs", height = "520px"),
                 tableOutput("labtab")),

        tabPanel("정보 About", br(),
                 htmlOutput("about"))
      )
    )
  )
)

## ---------------------------------------------------------------------
##  Scenario construction from the sidebar
## ---------------------------------------------------------------------
build_stages <- function(inp) {
  base <- list(GFR0 = inp$gfr * 60 / 1000, ALB = inp$alb, BW = inp$bw,
               AGEF = if (isTRUE(inp$aged)) 1 else 0)
  if (inp$mode == "acute") base <- c(base, MASSIVE)
  if (inp$mode == "ec")    base <- c(base, ENTERIC)

  stages <- list(list(t = 0, par = base))
  rx <- list()
  if (isTRUE(inp$fluid)) rx <- c(rx, list(RFLU = 0.25))
  if (isTRUE(inp$bic))   rx <- c(rx, list(RBIC = 75.0))
  rx <- c(rx, list(RKCL = inp$kcl))
  if (isTRUE(inp$dex))   rx <- c(rx, list(GLCP = 11.0))
  if (length(rx)) {
    stages <- c(stages, list(list(t = inp$trx, par = rx)))
    if (isTRUE(inp$bic))
      stages <- c(stages, list(list(t = inp$trx + 4, par = list(RBIC = 37.5))))
  }
  if (isTRUE(inp$hd)) {
    stages <- c(stages, list(list(t = inp$thd, par = list(CLHD = 6.0, HDBIC = 45))),
                list(list(t = inp$thd + 4, par = list(CLHD = 0, HDBIC = 0))))
  }
  if (inp$vent != "spont") {
    vp <- switch(inp$vent,
                 conv  = list(VENT = 1, VASET = 1.0),
                 match = list(VENT = 1, VASET = 2.3),
                 sed   = list(SED = 0.45))
    stages <- c(stages, list(list(t = inp$tvent, par = vp)))
  }
  stages[order(sapply(stages, `[[`, "t"))]
}

build_events <- function(inp) {
  ev <- if (inp$mode == "chronic") {
    dose_seq(inp$cdose, ii = 6, n = ceiling(inp$cdays * 4))
  } else if (inp$mode == "ec") {
    rbind(asa_dose(inp$dose * 0.55), asa_dose(inp$dose * 0.45, cmt = "ACONC"))
  } else {
    asa_dose(inp$dose)
  }
  if (isTRUE(inp$bic)) ev <- rbind(ev, bic_bolus(100, time = inp$trx))
  if (isTRUE(inp$ac))
    ev <- rbind(ev, charcoal(50, max(0.5, inp$trx - 2)),
                charcoal(50, inp$trx + 2), charcoal(50, inp$trx + 6))
  ev[order(ev$time), , drop = FALSE]
}

## ---------------------------------------------------------------------
##  SERVER
## ---------------------------------------------------------------------
server <- function(input, output, session) {

  sim <- eventReactive(input$go, ignoreNULL = FALSE, {
    stages <- build_stages(input)
    events <- build_events(input)
    run_stages(mod, stages, end = input$tmax, delta = 0.1, events = events)
  })

  long <- function(d, cols, labels = cols) {
    d %>% select(time, all_of(cols)) %>%
      pivot_longer(-time) %>%
      mutate(name = factor(name, levels = cols, labels = labels))
  }

  ## ---- 1. profile ----
  output$summary <- renderTable({
    d <- sim()
    data.frame(
      지표 = c("최고 총 혈장농도 (mg/dL)", "최고 유리농도 (mg/L)",
               "최고 뇌농도 (mg/L)", "최저 동맥혈 pH", "최저 HCO3 (mmol/L)",
               "최저 PaCO2 (mmHg)", "최고 요 pH", "최고 신장청소율 (mL/min)",
               "최저 K (mmol/L)", "최고 체온 (°C)", "최저 뇌 포도당 (mmol/L)",
               "최고 CNS 점수", "누적 CNS 손상"),
      값 = sprintf("%.2f", c(max(d$CTOTD), max(d$CFREE), max(d$CCNS),
                             min(d$PH), min(d$HCO3C), min(d$PACO2C),
                             max(d$UPH), max(d$CLREN), min(d$KPL),
                             max(d$TCORE), min(d$GLUB), max(d$CNS),
                             max(d$CNSI))))
  }, striped = TRUE)

  output$alerts <- renderUI({
    d <- sim(); a <- list()
    if (max(d$CTOTD) > 90) a <- c(a, "혈장농도 > 90 mg/dL — EXTRIP 투석 적응증")
    if (min(d$PH) < 7.20) a <- c(a, "산혈증 pH < 7.20 — 뇌 분배계수 급증, 투석 적응증")
    if (max(d$CNS) > 50) a <- c(a, "중추신경 억제 — 의식 변화, 투석 적응증")
    if (min(d$KPL) < 3.5) a <- c(a, "저칼륨혈증 — 요 알칼리화가 실패한다")
    if (max(d$UPH) < 7.0 && input$bic) a <- c(a, "요 pH가 오르지 않음 — K와 용량 확인")
    if (min(d$GLUB) < 0.8) a <- c(a, "뇌 저포도당 — 혈당이 정상이어도 포도당 투여")
    if (max(d$TCORE) > 38.5) a <- c(a, "고체온 — 탈공역의 표시, 예후 불량")
    if (input$vent == "conv") a <- c(a, "통상 환기 설정으로 삽관됨 — 이 모델에서 가장 위험한 조작")
    if (!length(a)) a <- "특이 경보 없음"
    HTML(paste0("<ul>", paste0("<li>", a, "</li>", collapse = ""), "</ul>"))
  })

  output$overview <- renderPlot({
    d <- sim()
    p <- long(d, c("CTOT", "CCNS", "PH", "UPH", "KPL", "CNS"),
              c("총 혈장 (mg/L)", "뇌 (mg/L)", "동맥혈 pH", "요 pH",
                "혈청 K (mmol/L)", "CNS 점수"))
    ggplot(p, aes(time, value)) + geom_line(linewidth = 0.9, colour = PAL["plasma"]) +
      facet_wrap(~name, scales = "free_y") + labs(x = "시간 (h)", y = NULL) + THEME
  })

  ## ---- 2. PK ----
  output$pk <- renderPlot({
    d <- sim()
    p <- long(d, c("CTOT", "CFREE", "CTIS", "CCNS"),
              c("총 혈장 (측정값)", "유리 혈장", "조직", "뇌"))
    ggplot(p, aes(time, value, colour = name)) + geom_line(linewidth = 1) +
      scale_colour_manual(values = unname(PAL[c("plasma","free","tissue","brain")])) +
      labs(x = "시간 (h)", y = "농도 (mg/L)", colour = NULL) + THEME
  })

  output$clr <- renderPlot({
    d <- sim()
    d2 <- d %>% mutate(
      salicylurate = c(0, diff(ASU) / diff(time)),
      glucuronide  = c(0, diff(APG + AAG) / diff(time)),
      renal        = c(0, diff(AUR) / diff(time)),
      dialysis     = c(0, diff(AHD) / diff(time)))
    p <- long(d2, c("salicylurate", "glucuronide", "renal", "dialysis"))
    ggplot(p, aes(time, value, colour = name)) + geom_line(linewidth = 1) +
      labs(x = "시간 (h)", y = "제거 속도 (mg/h)", colour = NULL) + THEME
  })

  ## ---- 3. the two multipliers ----
  output$mult <- renderPlot({
    d <- sim()
    p <- long(d, c("CTOT", "FU", "KPBR", "CCNS"),
              c("① C_total  (검사실 수치, mg/L)", "② fu  유리분율",
                "③ Kp  뇌:유리혈장 분배", "= C_brain  (mg/L)"))
    ggplot(p, aes(time, value)) +
      geom_line(linewidth = 1, colour = PAL["brain"]) +
      facet_wrap(~name, scales = "free_y", ncol = 2) +
      labs(x = "시간 (h)", y = NULL) + THEME
  })

  output$ratio <- renderPlot({
    d <- sim()
    ggplot(d, aes(time, RATIO)) + geom_line(linewidth = 1, colour = PAL["target"]) +
      labs(x = "시간 (h)", y = "뇌농도 / 총 혈장농도") + THEME
  })

  ## ---- 4. acid-base ----
  output$ab <- renderPlot({
    d <- sim()
    p <- long(d, c("PH", "PHB", "PACO2C", "HCO3C", "AGAP", "UNC"),
              c("동맥혈 pH", "뇌 세포내 pH", "PaCO2 (mmHg)",
                "HCO3 (mmol/L)", "음이온차 (mmol/L)", "탈공역 정도"))
    ggplot(p, aes(time, value)) + geom_line(linewidth = 0.9, colour = PAL["plasma"]) +
      facet_wrap(~name, scales = "free_y") + labs(x = "시간 (h)", y = NULL) + THEME
  })

  output$vent <- renderPlot({
    d <- sim()
    ggplot(d, aes(time, VASP)) + geom_line(linewidth = 1, colour = PAL["alt"]) +
      geom_hline(yintercept = 1, linetype = 2) +
      labs(x = "시간 (h)",
           y = "환자 자신의 폐포환기 (정상 대비 배수)") + THEME
  })

  ## ---- 5. renal ----
  output$renal <- renderPlot({
    d <- sim()
    p <- long(d, c("UPH", "CLREN", "UO", "KPL", "GFRC", "HCO3C"),
              c("요 pH", "신장 살리실산 청소율 (mL/min)", "요량 (mL/h)",
                "혈청 K (mmol/L)", "GFR (mL/min)", "혈장 HCO3 (mmol/L)"))
    ggplot(p, aes(time, value)) + geom_line(linewidth = 0.9, colour = PAL["alt"]) +
      facet_wrap(~name, scales = "free_y") + labs(x = "시간 (h)", y = NULL) + THEME
  })

  ## ---- 6. endpoints ----
  output$endp <- renderPlot({
    d <- sim()
    p <- long(d, c("CNS", "TINN", "SEV", "TCORE", "GLUB", "LUNG"),
              c("CNS 억제 점수", "이명 점수", "중증도 지수",
                "체온 (°C)", "뇌 포도당 (mmol/L)", "폐손상 (0-1)"))
    ggplot(p, aes(time, value)) + geom_line(linewidth = 0.9, colour = PAL["target"]) +
      facet_wrap(~name, scales = "free_y") + labs(x = "시간 (h)", y = NULL) + THEME
  })

  output$cnsi <- renderPlot({
    d <- sim()
    ggplot(d, aes(time, CNSI)) + geom_area(fill = "#EF9A9A", alpha = 0.6) +
      geom_line(linewidth = 1, colour = PAL["target"]) +
      labs(x = "시간 (h)", y = "누적 CNS 노출 ∫max(0, C_brain − 45)dt") + THEME
  })

  ## ---- 7. scenario comparison ----
  scen2 <- eventReactive(input$gosc, ignoreNULL = FALSE, {
    bind_rows(run_scenario(input$sc1, delta = 0.1),
              run_scenario(input$sc2, delta = 0.1))
  })

  output$scen <- renderPlot({
    d <- scen2()
    p <- d %>% select(time, scenario, CTOT, CCNS, PH, UPH, CLREN, CNSI) %>%
      pivot_longer(-c(time, scenario))
    ggplot(p, aes(time, value, colour = scenario)) + geom_line(linewidth = 1) +
      facet_wrap(~name, scales = "free_y") +
      labs(x = "시간 (h)", y = NULL, colour = NULL) + THEME
  })

  output$scentab <- renderTable({
    d <- scen2()
    bind_rows(lapply(split(d, d$scenario), summarise_scenario))
  })

  ## ---- 8. the ventilator experiment ----
  ventexp <- reactive({
    bind_rows(run_scenario("S07", delta = 0.05),
              run_scenario("S08", delta = 0.05)) %>%
      mutate(arm = ifelse(scenario == "S07", "통상 PaCO2로 삽관",
                          "분당환기량 일치"))
  })

  output$ventexp <- renderPlot({
    d <- ventexp() %>%
      select(time, arm, CTOT, CCNS, PH, PACO2C, KPBR, RATIO) %>%
      pivot_longer(-c(time, arm)) %>%
      mutate(name = recode(name,
        CTOT = "총 혈장 (검사실 수치, mg/L)", CCNS = "뇌 (mg/L)",
        PH = "동맥혈 pH", PACO2C = "PaCO2 (mmHg)",
        KPBR = "Kp 뇌:유리혈장", RATIO = "뇌 / 총 혈장"))
    ggplot(d, aes(time, value, colour = arm)) + geom_line(linewidth = 1) +
      geom_vline(xintercept = 8, linetype = 2, colour = "grey40") +
      facet_wrap(~name, scales = "free_y") +
      scale_colour_manual(values = c(PAL[["target"]], PAL[["alt"]])) +
      labs(x = "시간 (h)", y = NULL, colour = NULL) + THEME
  })

  output$venttab <- renderTable({
    d <- ventexp()
    at <- function(a, v, t) approx(d$time[d$arm == a], d[[v]][d$arm == a], t)$y
    ts <- c(8, 9, 12, 24, 36)
    data.frame(
      `시각 h` = ts,
      `혈장 A 통상` = round(sapply(ts, function(t) at("통상 PaCO2로 삽관", "CTOT", t)), 1),
      `혈장 B 일치` = round(sapply(ts, function(t) at("분당환기량 일치", "CTOT", t)), 1),
      `뇌 A` = round(sapply(ts, function(t) at("통상 PaCO2로 삽관", "CCNS", t)), 1),
      `뇌 B` = round(sapply(ts, function(t) at("분당환기량 일치", "CCNS", t)), 1),
      `뇌 A/B` = round(sapply(ts, function(t)
        at("통상 PaCO2로 삽관", "CCNS", t) / at("분당환기량 일치", "CCNS", t)), 3),
      check.names = FALSE)
  })

  ## ---- 9. Done nomogram ----
  done_zone <- function(c_mgdl, t_h) {
    t_h <- pmin(pmax(t_h, 6), 60)
    sev <- 130 * 0.5 ^ ((t_h - 6) / 20)
    mod <- 65 * 0.5 ^ ((t_h - 6) / 20)
    mild <- 32.5 * 0.5 ^ ((t_h - 6) / 20)
    ifelse(c_mgdl >= sev, "severe",
           ifelse(c_mgdl >= mod, "moderate",
                  ifelse(c_mgdl >= mild, "mild", "asymptomatic")))
  }
  brain_zone <- function(cb)
    ifelse(cb > 90, "severe", ifelse(cb > 55, "moderate",
                                     ifelse(cb > 25, "mild", "asymptomatic")))

  output$nomo <- renderPlot({
    d <- sim() %>% filter(time >= 6)
    bands <- data.frame(t = seq(6, 60, 0.5)) %>%
      mutate(severe = 130 * 0.5 ^ ((t - 6) / 20),
             moderate = 65 * 0.5 ^ ((t - 6) / 20),
             mild = 32.5 * 0.5 ^ ((t - 6) / 20)) %>%
      pivot_longer(-t)
    ggplot() +
      geom_line(data = bands, aes(t, value, linetype = name), colour = "grey40") +
      geom_path(data = d, aes(time, CTOTD, colour = brain_zone(CCNS)),
                linewidth = 1.4) +
      scale_y_log10() +
      labs(x = "섭취 후 시간 (h)", y = "혈장 살리실산 (mg/dL, log)",
           colour = "모델이 계산한 뇌농도 등급", linetype = "Done 경계") +
      THEME
  })

  output$nomotab <- renderTable({
    d <- sim()
    ts <- c(6, 12, 24, 36)
    ts <- ts[ts <= max(d$time)]
    ap <- function(v, t) approx(d$time, d[[v]], t)$y
    data.frame(
      `시각 h` = ts,
      `혈장 mg/dL` = round(sapply(ts, ap, v = "CTOTD"), 1),
      `Done 등급` = done_zone(sapply(ts, ap, v = "CTOTD"), ts),
      `뇌 mg/L` = round(sapply(ts, ap, v = "CCNS"), 1),
      `모델 등급` = brain_zone(sapply(ts, ap, v = "CCNS")),
      일치 = ifelse(done_zone(sapply(ts, ap, v = "CTOTD"), ts) ==
                      brain_zone(sapply(ts, ap, v = "CCNS")), "✓", "✗ 불일치"),
      check.names = FALSE)
  })

  ## ---- 10. EXTRIP ----
  output$extrip <- renderUI({
    d <- sim()
    crit <- c(
      sprintf("혈장 살리실산 > 90 mg/dL (신기능 정상): 최고 %.1f mg/dL", max(d$CTOTD)),
      sprintf("의식 변화 (CNS 점수 > 50): 최고 %.0f", max(d$CNS)),
      sprintf("동맥혈 pH ≤ 7.20: 최저 %.2f", min(d$PH)),
      sprintf("급성 폐손상: 최고 %.2f", max(d$LUNG)),
      sprintf("신기능 저하 (GFR < 45 mL/min): 최저 %.0f", min(d$GFRC)))
    met <- c(max(d$CTOTD) > 90, max(d$CNS) > 50, min(d$PH) <= 7.20,
             max(d$LUNG) > 0.3, min(d$GFRC) < 45)
    HTML(paste0(
      "<h4>EXTRIP 기준 대조 (Juurlink 2015)</h4><ul>",
      paste0("<li>", ifelse(met, "<b style='color:#C62828'>충족</b> — ", "미충족 — "),
             crit, "</li>", collapse = ""),
      "</ul><p>", if (any(met))
        "<b style='color:#C62828'>체외 제거를 고려할 상태입니다.</b>" else
        "현재 설정에서는 EXTRIP 기준을 충족하지 않습니다.",
      " 이 모델에서 투석 청소율(6 L/h)은 최적 알칼리뇨 청소율의 약 5.7배입니다.</p>"))
  })

  output$hdplot <- renderPlot({
    d <- bind_rows(run_scenario("S09", delta = 0.1),
                   run_scenario("S10", delta = 0.1)) %>%
      mutate(arm = ifelse(scenario == "S09", "투석 6시간", "투석 18시간"))
    p <- d %>% select(time, arm, CTOT, CCNS, CNSI, PH) %>%
      pivot_longer(-c(time, arm))
    ggplot(p, aes(time, value, colour = arm)) + geom_line(linewidth = 1) +
      facet_wrap(~name, scales = "free_y") +
      labs(x = "시간 (h)", y = NULL, colour = NULL) + THEME
  })

  ## ---- 11. laboratory panel ----
  output$labs <- renderPlot({
    d <- sim()
    p <- long(d, c("CTOTD", "PH", "HCO3C", "PACO2C", "AGAP", "KPL",
                   "GFRC", "UPH", "GLUB"),
              c("살리실산 (mg/dL)", "pH", "HCO3 (mmol/L)", "PaCO2 (mmHg)",
                "음이온차", "K (mmol/L)", "GFR (mL/min)", "요 pH",
                "뇌 포도당 (mmol/L)"))
    ggplot(p, aes(time, value)) + geom_line(linewidth = 0.9) +
      facet_wrap(~name, scales = "free_y") + labs(x = "시간 (h)", y = NULL) + THEME
  })

  output$labtab <- renderTable({
    d <- sim()
    ts <- unique(pmin(c(0, 2, 4, 6, 12, 24, 48, 72), max(d$time)))
    ap <- function(v) round(approx(d$time, d[[v]], ts)$y, 2)
    data.frame(`시각 h` = ts, `살리실산 mg/dL` = ap("CTOTD"), pH = ap("PH"),
               HCO3 = ap("HCO3C"), PaCO2 = ap("PACO2C"), `음이온차` = ap("AGAP"),
               K = ap("KPL"), `요 pH` = ap("UPH"), `뇌 mg/L` = ap("CCNS"),
               check.names = FALSE)
  })

  ## ---- about ----
  output$about <- renderUI(HTML('
    <h4>이 모델이 주장하는 것</h4>
    <p><code>C_brain = C_total × fu(C_total, pH) × f_n(pH_plasma)/f_n(pH_brain)</code></p>
    <p>검사실은 첫 항만 보고한다. 두 번째 항(유리분율)은 알부민이 포화되면서
       150 mg/L에서 8.5%, 800 mg/L에서 43%로 5배 커지고, 세 번째 항(pH 분배)은
       <b>혈장 pH는 움직이는데 뇌 세포내 pH는 완충되어 잘 움직이지 않기 때문에</b>
       커진다. PaCO2를 25에서 60 mmHg로 올리면 혈장 pH는 0.38 움직이지만
       뇌 세포내 pH는 0.14만 움직이고, 분배계수는 1.7배가 된다.</p>
    <h4>구성</h4>
    <ul>
      <li>28개 ODE · 14개 시나리오 · 11개 탭</li>
      <li>기계론적 지도: <code>sal_qsp_model.dot / .svg / .png</code> (125 노드, 16 클러스터)</li>
      <li>모델: <code>sal_mrgsolve_model.R</code></li>
      <li>독립 검증: <code>sal_verify_python.py</code> → <code>sal_verification_output.txt</code>
          (Python/scipy 재구현, 54개 문헌 기준 중 54개 통과)</li>
      <li>문헌: <code>sal_references.md</code> (PMID는 모두 PubMed에서 조회 확인)</li>
    </ul>
    <h4>한계</h4>
    <p>이 모델의 핵심 변수인 뇌 농도는 살아 있는 사람에게서 측정된 적이 없다.
       동물 분배 실험(Hill 1971)과 약산의 pH 물리학으로 계산했을 뿐이며,
       임상 결과에 맞춰 보정하지 않았다. 반증 가능한 예측은
       <b>PaCO2 변화 뒤 혈장과 뇌가 서로 반대 방향으로 움직인다</b>는 것이다.</p>
    <p style="color:#B71C1C"><b>교육·연구용 모델입니다. 임상 판단에 사용하지 마십시오.</b></p>'))
}

shinyApp(ui, server)
