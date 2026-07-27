##  Chronic Insomnia Disorder (CID) — QSP Simulator (Shiny)
##  ============================================================================
##  Front-end for ins_mrgsolve_model.R (63-ODE chronic insomnia QSP model).
##
##  Run with:
##      shiny::runApp("ins_shiny_app.R")
##  from inside the chronic-insomnia-disorder/ directory, or
##      shiny::runApp("chronic-insomnia-disorder/ins_shiny_app.R")
##  from the repository root.
##
##  Requires: shiny, mrgsolve, dplyr, tidyr, ggplot2, DT
##
##  Design intent
##  -------------
##  Insomnia is the one sleep disorder where the patient's own behaviour is
##  part of the pathophysiology, so the app is laid out to keep the SCHEDULE
##  (time in bed, rise time, light) visible next to the pharmacology on every
##  tab. Two things the app is built to make obvious:
##
##    1. A single night and twelve weeks are different questions. The
##       "Two-process & switch" tab shows why a patient falls asleep tonight;
##       the "Perpetuating loop" and "Endpoints" tabs show why they will still
##       be an insomniac in three months. Hypnotics change the first plot and
##       not the second; CBT-I does the reverse.
##
##    2. Every arm is compared against its OWN untreated control, run from the
##       same lead-in state, because the phenotypes have different baselines.
##
##  Tabs
##    1  환자 프로파일        Patient & phenotype — arousal decomposition
##    2  두 과정 & 스위치     Process S x C, wake drive, one-night sleep switch
##    3  일주기 & 멜라토닌    Circadian phase, light, melatonin, temperature
##    4  약물 PK / 표적점유   Concentrations and receptor occupancy
##    5  야간 수면 지표       SOL / WASO / TST / SE trajectories
##    6  수면 구조            N3, REM, awakenings, hypnogram surrogate
##    7  유지 고리 (3-P)      COND / SEFF / DBAS / TIB and the CBT-I response
##    8  임상 엔드포인트      ISI, subjective sleep, remission, safety
##    9  시나리오 비교        Several arms side by side with a summary table
##  ============================================================================

## NOTE: mrgsolve exports its own req(), which masks shiny::req() because it
## is attached later. Shiny's version is called with its namespace below.
library(shiny)
library(mrgsolve)
library(dplyr)
library(tidyr)
library(ggplot2)
library(DT)

source("ins_mrgsolve_model.R")

## ---------------------------------------------------------------- helpers --
PHENOS <- c(
  "좋은 수면자 (good sleeper)"                = "good_sleeper",
  "만성 불면장애 (chronic insomnia)"           = "chronic_insomnia",
  "중증 불면장애 (severe)"                     = "severe_insomnia",
  "입면곤란 우세형 (sleep-onset)"              = "sleep_onset",
  "수면유지곤란 우세형 (sleep-maintenance)"    = "sleep_maintenance",
  "고령 (elderly)"                             = "elderly",
  "우울 동반 (comorbid depression)"            = "comorbid_depression",
  "수면위상지연 (delayed phase)"               = "delayed_phase",
  "교대근무자 (shift worker)"                  = "shift_worker"
)

DRUGS <- list(
  "없음 (none)"                 = list(cmt = NA,     dose = 0),
  "졸피뎀 Zolpidem"             = list(cmt = "ZOLD", dose = 10),
  "에스조피클론 Eszopiclone"    = list(cmt = "ESZD", dose = 3),
  "수보렉산트 Suvorexant"       = list(cmt = "SUVD", dose = 20),
  "렘보렉산트 Lemborexant"      = list(cmt = "LEMD", dose = 10),
  "다리도렉산트 Daridorexant"   = list(cmt = "DARD", dose = 50),
  "라멜테온 Ramelteon"          = list(cmt = "RAMD", dose = 8),
  "멜라토닌 Melatonin"          = list(cmt = "MELD", dose = 2),
  "저용량 독세핀 Doxepin"       = list(cmt = "DOXD", dose = 6),
  "트라조돈 Trazodone"          = list(cmt = "TRZD", dose = 50)
)

theme_ins <- function() {
  theme_minimal(base_size = 12) +
    theme(panel.grid.minor = element_blank(),
          strip.text = element_text(face = "bold"),
          legend.position = "bottom",
          plot.title = element_text(face = "bold"))
}

night_shade <- function(days, waket, tib) {
  bed <- (waket + 24 - tib) %% 24
  data.frame(
    xmin = sapply(0:(days - 1), function(d) d*24 + bed - ifelse(bed > waket, 0, 24)),
    xmax = sapply(0:(days - 1), function(d) d*24 + waket)
  )
}

## -------------------------------------------------------------------- UI ---
ui <- fluidPage(
  titlePanel("만성 불면장애 QSP 시뮬레이터 / Chronic Insomnia Disorder QSP Simulator"),
  tags$p(style = "color:#666;margin-top:-8px",
         "Two-process (S x C) core · VLPO/arousal flip-flop · 24-h hyperarousal · ",
         "Spielman 3-P perpetuating loop · BzRA / DORA / melatonergic / H1 pharmacology"),

  sidebarLayout(
    sidebarPanel(
      width = 3,
      h4("① 환자 / Patient"),
      selectInput("pheno", "표현형 Phenotype", choices = PHENOS,
                  selected = "chronic_insomnia"),
      sliderInput("A0", "구성적 과각성 A0", 0.20, 1.20, 0.74, step = 0.02),
      sliderInput("TIB0", "습관적 침상시간 TIB (h)", 5.0, 10.0, 8.8, step = 0.1),
      sliderInput("WAKET", "기상 시각 (clock h)", 4, 11, 7, step = 0.25),
      sliderInput("TSTNEED", "개인 수면 요구량 (min)", 360, 540, 450, step = 10),

      h4("② 동반질환 / Comorbidity"),
      sliderInput("PAIN", "만성 통증 부담 (0-1)", 0, 1, 0, step = 0.05),
      sliderInput("VMS",  "혈관운동증상 부담 (0-1)", 0, 1, 0, step = 0.05),
      sliderInput("STR0", "촉발 스트레스 잔여 (0-1)", 0, 1, 0, step = 0.05),

      h4("③ 약물 / Drug"),
      selectInput("drug", "야간 수면제", choices = names(DRUGS),
                  selected = "없음 (none)"),
      numericInput("dose", "1회 용량 (mg)", 10, min = 0, max = 200, step = 0.5),
      sliderInput("txweeks", "투약 기간 (주)", 0, 12, 12, step = 1),
      sliderInput("bedhour", "복용 시각 (clock h)", 19, 24, 22, step = 0.25),
      sliderInput("CYP3A", "CYP3A4 활성 배수", 0.2, 2.0, 1.0, step = 0.05),

      h4("④ 비약물 / Non-pharmacological"),
      checkboxInput("cbti", "CBT-I 시행", FALSE),
      sliderInput("cbtiday", "CBT-I 시작일", 1, 42, 1, step = 1),
      checkboxInput("blt", "아침 고조도 광치료", FALSE),
      sliderInput("caf", "카페인 (mg, 16시)", 0, 400, 0, step = 25),
      sliderInput("etoh", "취침 전 알코올 (g)", 0, 60, 0, step = 5),
      sliderInput("luxeve", "저녁 조도 (lux)", 20, 600, 120, step = 10),

      hr(),
      sliderInput("weeks", "관찰 기간 (주)", 2, 12, 12, step = 1),
      actionButton("run", "시뮬레이션 실행 / Run", class = "btn-primary btn-block"),
      tags$small(style = "color:#888",
        "각 실행은 150일 무치료 lead-in으로 환자를 만성 평형 상태에 도달시킨 뒤 ",
        "무작위배정 시점부터의 결과만 표시합니다. 대조군은 동일한 lead-in에서 ",
        "출발한 무치료 아암입니다.")
    ),

    mainPanel(
      width = 9,
      tabsetPanel(
        id = "tabs",

        tabPanel("1 환자 프로파일",
          fluidRow(
            column(6, plotOutput("p_arousal", height = 320)),
            column(6, plotOutput("p_hpa", height = 320))),
          fluidRow(column(12, h4("각성 구성요소 분해 / Arousal decomposition"),
                          DTOutput("t_profile")))),

        tabPanel("2 두 과정 & 스위치",
          fluidRow(column(12, plotOutput("p_night", height = 380))),
          fluidRow(
            column(6, plotOutput("p_switch", height = 300)),
            column(6, plotOutput("p_drivecomp", height = 300))),
          tags$p(style = "color:#666",
            "위: 대표적인 하룻밤. 수면 확률 SLP는 sigmoid((S·wS + 아데노신 − C(t) − 각성 ",
            "− 상행각성구동 + 약물효과 − θ)/k)이며, 각성계 톤이 SLP에 의해 다시 억제되므로 ",
            "이 스위치는 쌍안정(bistable)입니다. 입면과 각성이 매끄럽지 않고 '넘어가는' ",
            "이유입니다.")),

        tabPanel("3 일주기 & 멜라토닌",
          fluidRow(
            column(6, plotOutput("p_phi", height = 300)),
            column(6, plotOutput("p_mel", height = 300))),
          fluidRow(
            column(6, plotOutput("p_temp", height = 280)),
            column(6, plotOutput("p_light", height = 280)))),

        tabPanel("4 약물 PK / 표적점유",
          fluidRow(column(12, plotOutput("p_pk", height = 330))),
          fluidRow(column(12, plotOutput("p_occ", height = 330)))),

        tabPanel("5 야간 수면 지표",
          fluidRow(column(12, plotOutput("p_metrics", height = 420))),
          fluidRow(column(12, DTOutput("t_weekly")))),

        tabPanel("6 수면 구조",
          fluidRow(
            column(6, plotOutput("p_stage", height = 320)),
            column(6, plotOutput("p_stagepct", height = 320))),
          fluidRow(column(12, plotOutput("p_naw", height = 260)))),

        tabPanel("7 유지 고리 (3-P)",
          fluidRow(column(12, plotOutput("p_loop", height = 400))),
          fluidRow(column(12, plotOutput("p_tib", height = 260))),
          tags$p(style = "color:#666",
            "이 탭이 이 모델의 요점입니다. 수면제는 이 네 상태변수를 거의 건드리지 ",
            "않으므로 중단하면 효과가 사라지고, CBT-I는 이 네 변수를 직접 표적하므로 ",
            "치료 종료 후에도 효과가 남습니다.")),

        tabPanel("8 임상 엔드포인트",
          fluidRow(
            column(6, plotOutput("p_isi", height = 320)),
            column(6, plotOutput("p_subj", height = 320))),
          fluidRow(
            column(6, plotOutput("p_safety", height = 300)),
            column(6, plotOutput("p_conseq", height = 300))),
          fluidRow(column(12, DTOutput("t_endpoint")))),

        tabPanel("9 시나리오 비교",
          fluidRow(column(12,
            checkboxGroupInput("scen", "비교할 시나리오",
              choices = names(ins_scenarios()),
              selected = names(ins_scenarios())[c(1, 2, 7, 13)], inline = FALSE),
            actionButton("runscen", "선택 시나리오 실행", class = "btn-success"))),
          fluidRow(column(12, plotOutput("p_scen", height = 420))),
          fluidRow(column(12, DTOutput("t_scen"))))
      )
    )
  )
)

## ----------------------------------------------------------------- SERVER --
server <- function(input, output, session) {

  observeEvent(input$pheno, {
    p <- ins_phenotype(input$pheno)
    if (!is.null(p$A0))   updateSliderInput(session, "A0",   value = p$A0)
    if (!is.null(p$TIB0)) updateSliderInput(session, "TIB0", value = p$TIB0)
    if (!is.null(p$WAKET))updateSliderInput(session, "WAKET",value = p$WAKET)
  })

  observeEvent(input$drug, {
    d <- DRUGS[[input$drug]]
    updateNumericInput(session, "dose", value = d$dose)
  })

  build_events <- reactive({
    ev <- NULL
    d <- DRUGS[[input$drug]]
    nd <- input$txweeks*7
    if (!is.na(d$cmt) && input$dose > 0 && nd > 0)
      ev <- ev_nightly(d$cmt, input$dose, days = nd, bedhour = input$bedhour)
    if (input$caf > 0)
      ev <- ev %+% ev_daily("CAFD", input$caf, days = input$weeks*7, hour = 16)
    if (input$etoh > 0)
      ev <- ev %+% ev_nightly("ETHD", input$etoh, days = input$weeks*7,
                              bedhour = max(19, input$bedhour - 0.5))
    ev
  })

  user_param <- reactive({
    list(A0 = input$A0, TIB0 = input$TIB0, WAKET = input$WAKET,
         TSTNEED = input$TSTNEED, PAIN = input$PAIN, VMS = input$VMS,
         LUXEVE = input$luxeve)
  })

  treat_param <- reactive({
    p <- list(CYP3A = input$CYP3A)
    if (input$blt) p <- c(p, list(BLTON = 1, BLTLUX = 10000,
                                  BLTT0 = input$WAKET, BLTDUR = 0.75))
    p
  })

  sim <- eventReactive(input$run, {
    withProgress(message = "시뮬레이션 중 (lead-in 150일 + 관찰기간)...", value = 0.3, {
      act <- ins_simulate(phenotype = input$pheno, weeks = input$weeks,
                          param = user_param(), treat = treat_param(),
                          events = build_events(),
                          cbti_start_day = if (input$cbti) input$cbtiday else NA_real_)
      incProgress(0.5, detail = "대조군")
      ctl <- ins_simulate(phenotype = input$pheno, weeks = input$weeks,
                          param = user_param(), treat = list(),
                          events = NULL, cbti_start_day = NA_real_)
      list(act = act, ctl = ctl,
           nact = ins_nightly(act), nctl = ins_nightly(ctl))
    })
  }, ignoreNULL = FALSE)

  ## a representative night, taken from the middle of the observation window
  night <- reactive({
    s <- sim()$act
    d0 <- floor(input$weeks*7/2)
    s[s$time >= (d0*24 + 16) & s$time <= (d0*24 + 40), ]
  })

  long_state <- function(df, vars, labs) {
    out <- df[, c("time", vars)]
    names(out) <- c("time", labs)
    tidyr::pivot_longer(out, -time, names_to = "state", values_to = "value")
  }

  ## ---- tab 1 --------------------------------------------------------------
  output$p_arousal <- renderPlot({
    s <- sim()$act
    d <- long_state(s, c("AROU", "COND", "SEFF", "DBAS"),
                    c("과각성 AROU", "조건화각성 COND", "수면노력 SEFF", "신념 DBAS"))
    ggplot(d, aes(time/24, value, colour = state)) + geom_line(linewidth = 0.6) +
      labs(x = "일 (day)", y = "상태값", colour = NULL,
           title = "과각성과 그 구동 요소") + theme_ins()
  })

  output$p_hpa <- renderPlot({
    s <- night()
    d <- long_state(s, c("CORT", "CRH", "ACTH"),
                    c("코르티솔 (ug/dL)", "CRH", "ACTH"))
    ggplot(d, aes(time %% 24, value, colour = state)) + geom_line(linewidth = 0.7) +
      labs(x = "시각 (clock h)", y = NULL, colour = NULL,
           title = "HPA 축 — 야간 코르티솔 nadir 상승") + theme_ins()
  })

  output$t_profile <- renderDT({
    s <- sim()$act; n <- nrow(s)
    a <- sim()$nact; b <- sim()$nctl
    wk <- function(x) x[x$night > max(x$night) - 7, ]
    out <- data.frame(
      item = c("과각성 AROU", "조건화각성 COND", "수면노력 SEFF", "역기능적 신념 DBAS",
               "침상시간 TIB (h)", "일주기 위상 PHI (h)", "SOL (min)", "WASO (min)",
               "TST (min)", "SE (%)", "ISI"),
      active = round(c(s$AROU[n], s$COND[n], s$SEFF[n], s$DBAS[n], s$TIBS[n], s$PHI[n],
                       mean(wk(a)$SOL), mean(wk(a)$WASO), mean(wk(a)$TST),
                       mean(wk(a)$SE), mean(wk(a)$ISI)), 2),
      control = round(c(sim()$ctl$AROU[n], sim()$ctl$COND[n], sim()$ctl$SEFF[n],
                        sim()$ctl$DBAS[n], sim()$ctl$TIBS[n], sim()$ctl$PHI[n],
                        mean(wk(b)$SOL), mean(wk(b)$WASO), mean(wk(b)$TST),
                        mean(wk(b)$SE), mean(wk(b)$ISI)), 2)
    )
    names(out) <- c("\uD56D\uBAA9", "\uCE58\uB8CC\uAD70", "\uB300\uC870\uAD70")
    out
  }, options = list(dom = "t", pageLength = 12))

  ## ---- tab 2 --------------------------------------------------------------
  output$p_night <- renderPlot({
    s <- night()
    d <- long_state(s, c("S", "Cwake", "AROU", "Wdrv", "drive", "SLP"),
                    c("Process S", "C(t) 각성신호", "과각성 A", "상행각성구동",
                      "순 수면구동 drive", "수면확률 SLP"))
    ggplot(d, aes(time %% 24, value, colour = state)) +
      geom_hline(yintercept = 0, colour = "grey70", linetype = 2) +
      geom_line(linewidth = 0.7) +
      scale_x_continuous(breaks = seq(0, 24, 3)) +
      labs(x = "시각 (clock h)", y = NULL, colour = NULL,
           title = "대표적인 하룻밤 — 두 과정, 각성, 스위치") + theme_ins()
  })

  output$p_switch <- renderPlot({
    s <- night()
    ggplot(s, aes(drive, SLP)) + geom_path(linewidth = 0.7, colour = "#3949ab") +
      labs(x = "순 수면구동", y = "수면확률",
           title = "스위치 궤적 (히스테리시스)") + theme_ins()
  })

  output$p_drivecomp <- renderPlot({
    s <- night()
    d <- long_state(s, c("adoEff", "Ebz", "Emt", "oxBlk", "h1blk"),
                    c("유효 아데노신", "BzRA 효과", "멜라토닌계 효과",
                      "오렉신 차단", "H1 차단"))
    ggplot(d, aes(time %% 24, value, colour = state)) + geom_line(linewidth = 0.7) +
      labs(x = "시각 (clock h)", y = NULL, colour = NULL,
           title = "약리학적 구동 성분") + theme_ins()
  })

  ## ---- tab 3 --------------------------------------------------------------
  output$p_phi <- renderPlot({
    s <- sim()$act; c0 <- sim()$ctl
    d <- rbind(data.frame(time = s$time, PHI = s$PHI, arm = "치료군"),
               data.frame(time = c0$time, PHI = c0$PHI, arm = "대조군"))
    ggplot(d, aes(time/24, PHI, colour = arm)) + geom_line(linewidth = 0.6) +
      labs(x = "일", y = "위상 이동 Φ (h, + = 지연)", colour = NULL,
           title = "일주기 위상 — 광 PRC로 동조된 결과") + theme_ins()
  })

  output$p_mel <- renderPlot({
    s <- night()
    ggplot(s, aes(time %% 24)) +
      geom_line(aes(y = MELP, colour = "내인성 (pg/mL)"), linewidth = 0.8) +
      geom_line(aes(y = CMLX, colour = "외인성 (pg/mL)"), linewidth = 0.8) +
      geom_line(aes(y = 100*occMT, colour = "MT 점유율 (%)"), linetype = 2) +
      labs(x = "시각 (clock h)", y = NULL, colour = NULL,
           title = "멜라토닌 — 내인성 · 외인성 · 수용체 점유") + theme_ins()
  })

  output$p_temp <- renderPlot({
    s <- night()
    ggplot(s, aes(time %% 24, TEMPC)) + geom_line(colour = "#c62828", linewidth = 0.8) +
      labs(x = "시각 (clock h)", y = "심부체온 (°C)",
           title = "심부체온 — 최저점 부근에서 수면 유지가 가장 안정") + theme_ins()
  })

  output$p_light <- renderPlot({
    s <- night()
    ggplot(s, aes(time %% 24, lux + 1)) + geom_line(colour = "#f9a825", linewidth = 0.8) +
      scale_y_log10() +
      labs(x = "시각 (clock h)", y = "조도 (lux, log)",
           title = "광 노출 프로파일") + theme_ins()
  })

  ## ---- tab 4 --------------------------------------------------------------
  output$p_pk <- renderPlot({
    s <- night()
    cs <- c("CZOL", "CESZ", "CSUV", "CLEM", "CDAR", "CRAM", "CRM2",
            "CDOX", "CTRZ", "CCAF", "CETH")
    keep <- cs[sapply(cs, function(v) max(s[[v]], na.rm = TRUE) > 1e-9)]
    if (!length(keep)) return(ggplot() + labs(title = "투여된 약물 없음") + theme_ins())
    d <- long_state(s, keep, keep)
    ggplot(d, aes(time %% 24, value, colour = state)) + geom_line(linewidth = 0.8) +
      labs(x = "시각 (clock h)", y = "농도 (mg/L; 알코올 g/L)", colour = NULL,
           title = "혈장 농도 — 대표적인 하룻밤") + theme_ins()
  })

  output$p_occ <- renderPlot({
    s <- night()
    d <- long_state(s, c("occBZ", "oxBlk", "occMT", "h1blk", "e5ht"),
                    c("GABA-A BZ 부위", "오렉신 수용체 차단", "MT1/MT2 점유",
                      "H1 차단", "5-HT2A 차단"))
    ggplot(d, aes(time %% 24, 100*value, colour = state)) + geom_line(linewidth = 0.8) +
      labs(x = "시각 (clock h)", y = "점유/차단 (%)", colour = NULL,
           title = "표적 점유율") + theme_ins()
  })

  ## ---- tab 5 --------------------------------------------------------------
  metrics_long <- reactive({
    a <- sim()$nact; b <- sim()$nctl
    f <- function(d, arm) data.frame(
      night = d$night, arm = arm,
      SOL = d$SOL, WASO = d$WASO, TST = d$TST, SE = d$SE)
    tidyr::pivot_longer(rbind(f(a, "치료군"), f(b, "대조군")),
                        c(SOL, WASO, TST, SE), names_to = "metric")
  })

  output$p_metrics <- renderPlot({
    ggplot(metrics_long(), aes(night, value, colour = arm)) +
      geom_line(linewidth = 0.6) +
      facet_wrap(~metric, scales = "free_y") +
      labs(x = "밤 (night)", y = NULL, colour = NULL,
           title = "야간 수면 지표 — 치료군 대 무치료 대조군") + theme_ins()
  })

  output$t_weekly <- renderDT({
    a <- sim()$nact; b <- sim()$nctl
    wk <- function(d) d %>% mutate(week = ceiling(night/7)) %>%
      group_by(week) %>% summarise(across(c(SOL, WASO, TST, SE, ISI), mean),
                                   .groups = "drop")
    A <- wk(a); B <- wk(b)
    out <- data.frame(week = A$week,
      SOL = round(A$SOL), SOLd = round(A$SOL - B$SOL),
      WASO = round(A$WASO), WASOd = round(A$WASO - B$WASO),
      TST = round(A$TST), TSTd = round(A$TST - B$TST),
      SEpct = round(A$SE, 1), ISI = round(A$ISI, 1),
      ISId = round(A$ISI - B$ISI, 1))
    names(out) <- c("\uC8FC\uCC28", "SOL", "SOL \u0394", "WASO", "WASO \u0394",
                    "TST", "TST \u0394", "SE %", "ISI", "ISI \u0394")
    datatable(out, options = list(dom = "tp", pageLength = 12), rownames = FALSE)
  })

  ## ---- tab 6 --------------------------------------------------------------
  output$p_stage <- renderPlot({
    a <- sim()$nact
    d <- tidyr::pivot_longer(a[, c("night", "N3", "REM")], -night)
    ggplot(d, aes(night, value, colour = name)) + geom_line(linewidth = 0.7) +
      labs(x = "밤", y = "분/밤", colour = NULL,
           title = "수면 단계 — 절대 시간") + theme_ins()
  })

  output$p_stagepct <- renderPlot({
    a <- sim()$nact
    d <- data.frame(night = a$night,
                    "N3 %" = 100*a$N3/pmax(a$TST, 1),
                    "REM %" = 100*a$REM/pmax(a$TST, 1), check.names = FALSE)
    d <- tidyr::pivot_longer(d, -night)
    ggplot(d, aes(night, value, colour = name)) + geom_line(linewidth = 0.7) +
      labs(x = "밤", y = "% of TST", colour = NULL,
           title = "수면 단계 — TST 대비 비율") + theme_ins()
  })

  output$p_naw <- renderPlot({
    a <- sim()$nact; b <- sim()$nctl
    d <- rbind(data.frame(night = a$night, NAW = a$NAW, arm = "치료군"),
               data.frame(night = b$night, NAW = b$NAW, arm = "대조군"))
    ggplot(d, aes(night, NAW, colour = arm)) + geom_line(linewidth = 0.6) +
      labs(x = "밤", y = "각성 횟수", colour = NULL,
           title = "야간 각성 횟수") + theme_ins()
  })

  ## ---- tab 7 --------------------------------------------------------------
  output$p_loop <- renderPlot({
    s <- sim()$act; c0 <- sim()$ctl
    mk <- function(df, arm) long_state(df, c("COND", "SEFF", "DBAS", "MISP"),
      c("조건화각성 COND", "수면노력 SEFF", "역기능적 신념 DBAS", "수면 오지각 MISP")) %>%
      mutate(arm = arm)
    d <- rbind(mk(s, "치료군"), mk(c0, "대조군"))
    ggplot(d, aes(time/24, value, colour = arm)) + geom_line(linewidth = 0.6) +
      facet_wrap(~state, scales = "free_y") +
      labs(x = "일", y = NULL, colour = NULL,
           title = "3-P 유지 고리 상태변수") + theme_ins()
  })

  output$p_tib <- renderPlot({
    s <- sim()$act; c0 <- sim()$ctl
    d <- rbind(data.frame(time = s$time, TIB = s$TIBS, SE = 100*s$SEBAR, arm = "치료군"),
               data.frame(time = c0$time, TIB = c0$TIBS, SE = 100*c0$SEBAR, arm = "대조군"))
    ggplot(d, aes(time/24)) +
      geom_line(aes(y = TIB, colour = arm), linewidth = 0.7) +
      geom_line(aes(y = SE/10, colour = arm), linetype = 2) +
      scale_y_continuous("침상시간 TIB (h)",
                         sec.axis = sec_axis(~.*10, name = "수면효율 (%)")) +
      labs(x = "일", colour = NULL,
           title = "침상시간 적정 — 실선 TIB, 점선 수면효율") + theme_ins()
  })

  ## ---- tab 8 --------------------------------------------------------------
  output$p_isi <- renderPlot({
    s <- sim()$act; c0 <- sim()$ctl
    d <- rbind(data.frame(time = s$time, ISI = s$ISI, arm = "치료군"),
               data.frame(time = c0$time, ISI = c0$ISI, arm = "대조군"))
    ggplot(d, aes(time/24, ISI, colour = arm)) + geom_line(linewidth = 0.7) +
      geom_hline(yintercept = 8, linetype = 2, colour = "grey50") +
      annotate("text", x = 1, y = 8.8, label = "관해 기준 ISI<8",
               hjust = 0, size = 3.2, colour = "grey40") +
      labs(x = "일", y = "ISI (0-28)", colour = NULL,
           title = "불면 중증도 지수") + theme_ins()
  })

  output$p_subj <- renderPlot({
    a <- sim()$nact
    d <- data.frame(night = a$night, "objective TST" = a$TST,
                    "subjective TST" = a$sTST, check.names = FALSE)
    d <- tidyr::pivot_longer(d, -night)
    ggplot(d, aes(night, value, colour = name)) + geom_line(linewidth = 0.7) +
      labs(x = "밤", y = "분", colour = NULL,
           title = "객관 대 주관 수면시간 — 오지각의 크기") + theme_ins()
  })

  output$p_safety <- renderPlot({
    s <- sim()$act
    d <- long_state(s, c("TOLB", "WDR", "RESBAR"),
                    c("BzRA 내성", "금단/반동", "주간 잔류진정"))
    ggplot(d, aes(time/24, value, colour = state)) + geom_line(linewidth = 0.7) +
      labs(x = "일", y = NULL, colour = NULL,
           title = "내성 · 반동 · 잔류 진정") + theme_ins()
  })

  output$p_conseq <- renderPlot({
    s <- sim()$act
    d <- long_state(s, c("DEP", "IL6", "GLYM"),
                    c("우울 부담", "IL-6 (pg/mL)", "글림프 청소 부채"))
    ggplot(d, aes(time/24, value, colour = state)) + geom_line(linewidth = 0.7) +
      labs(x = "일", y = NULL, colour = NULL,
           title = "장기 결과 대리지표") + theme_ins()
  })

  output$t_endpoint <- renderDT({
    a <- sim()$nact; b <- sim()$nctl
    l <- function(d) d[d$night > max(d$night) - 7, ]
    A <- l(a); B <- l(b)
    out <- data.frame(
      endpoint = c("SOL (min)", "WASO (min)", "TST (min)", "SE (%)",
                     "주관 TST (min)", "ISI", "N3 (min)", "REM (min)",
                     "숙취 지수 (0-100)", "관해 (ISI<8 & SE>85%)"),
      active = c(round(mean(A$SOL)), round(mean(A$WASO)), round(mean(A$TST)),
                 round(mean(A$SE), 1), round(mean(A$sTST)), round(mean(A$ISI), 1),
                 round(mean(A$N3)), round(mean(A$REM)), round(mean(A$HANG)),
                 ifelse(mean(A$ISI) < 8 && mean(A$SE) > 85, "예", "아니오")),
      control = c(round(mean(B$SOL)), round(mean(B$WASO)), round(mean(B$TST)),
                  round(mean(B$SE), 1), round(mean(B$sTST)), round(mean(B$ISI), 1),
                  round(mean(B$N3)), round(mean(B$REM)), round(mean(B$HANG)),
                  ifelse(mean(B$ISI) < 8 && mean(B$SE) > 85, "\uC608", "\uC544\uB2C8\uC624"))
    )
    names(out) <- c("\uC5D4\uB4DC\uD3EC\uC778\uD2B8", "\uCE58\uB8CC\uAD70",
                    "\uB300\uC870\uAD70")
    out
  }, options = list(dom = "t", pageLength = 12))

  ## ---- tab 9 --------------------------------------------------------------
  scen <- eventReactive(input$runscen, {
    shiny::req(length(input$scen) > 0)   # mrgsolve also exports req()
    withProgress(message = "시나리오 실행 중...", value = 0, {
      out <- lapply(seq_along(input$scen), function(i) {
        incProgress(1/length(input$scen), detail = input$scen[i])
        r <- ins_run_scenario(input$scen[i], weeks = input$weeks)
        r$nightly$scenario <- input$scen[i]
        r$nightly
      })
      do.call(rbind, out)
    })
  })

  output$p_scen <- renderPlot({
    d <- scen()
    dl <- tidyr::pivot_longer(d[, c("night", "scenario", "SOL", "WASO", "TST", "ISI")],
                              c(SOL, WASO, TST, ISI), names_to = "metric")
    ggplot(dl, aes(night, value, colour = scenario)) + geom_line(linewidth = 0.6) +
      facet_wrap(~metric, scales = "free_y") +
      labs(x = "밤", y = NULL, colour = NULL, title = "시나리오 비교") +
      theme_ins() + theme(legend.text = element_text(size = 8))
  })

  output$t_scen <- renderDT({
    d <- scen()
    s <- d %>% group_by(scenario) %>%
      filter(night > max(night) - 7) %>%
      summarise(SOL = round(mean(SOL)), WASO = round(mean(WASO)),
                TST = round(mean(TST)), `SE %` = round(mean(SE), 1),
                ISI = round(mean(ISI), 1), sTST = round(mean(sTST)),
                N3 = round(mean(N3)), REM = round(mean(REM)),
                HANG = round(mean(HANG)), .groups = "drop")
    datatable(s, options = list(dom = "t"), rownames = FALSE)
  })
}

shinyApp(ui, server)
