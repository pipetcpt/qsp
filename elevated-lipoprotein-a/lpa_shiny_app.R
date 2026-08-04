## =====================================================================
##  lpa_shiny_app.R
##  Elevated Lipoprotein(a) — QSP interactive dashboard
##  10 tabs · drives the 50-ODE model in lpa_mrgsolve_model.R
##
##  고지단백(a)혈증 — QSP 대시보드
## ---------------------------------------------------------------------
##  The app is organised around the model's thesis rather than around the
##  compartment list.  Each tab answers ONE question:
##
##   1  환자 프로파일   Who is this patient, and what does their KIV-2
##                      repeat number do to everything downstream?
##   2  생산 · 조립     Where in the production chain does each drug act?
##   3  약동학          PK of eight agents on one axis.
##   4  Lp(a) 반응      The headline PD curve, four assays at once.
##   5  ★ 측정의 함정   The tab the model exists for: the SAME patient
##                      reported four ways, and who gets missed.
##   6  지질 · apoB     LDL-C contamination and the Dahlen correction.
##   7  염증 축         OxPL -> monocyte -> IL-6 -> CRP, and the loop gain.
##   8  혈관 · 위험     Plaque, vulnerability, and the MR-vs-trial split.
##   9  대동맥판막      Why starting age, not potency, decides this one.
##  10  시나리오 비교   All 20 scenarios side by side.
##
##  Requires: shiny, mrgsolve, dplyr, tidyr, ggplot2, DT
##  Run with:  shiny::runApp("lpa_shiny_app.R")
## =====================================================================

library(shiny)
library(mrgsolve)
library(dplyr)
library(tidyr)
library(ggplot2)
library(DT)

## Load the model from the companion file (keeps ONE source of truth for
## the equations; the app never redefines them).
LPA_MODEL_ONLY <- TRUE
source("lpa_mrgsolve_model.R", local = FALSE, echo = FALSE)
## `mod` is now available.

YR <- 365.25

DRUGS <- c("없음 (None)"                 = "none",
           "고강도 스타틴 (Statin)"       = "statin",
           "에볼로쿠맙 (Evolocumab)"      = "evolocumab",
           "펠라카르센 ASO (Pelacarsen)"   = "pelacarsen",
           "올파시란 siRNA (Olpasiran)"    = "olpasiran",
           "레포디시란 siRNA (Lepodisiran)"= "lepodisiran",
           "무발라플린 경구 (Muvalaplin)"  = "muvalaplin",
           "니아신 (Niacin ER)"           = "niacin",
           "오비세트라핍 (Obicetrapib)"    = "obicetrapib",
           "질티베키맙 (Ziltivekimab)"     = "ziltivekimab",
           "지단백 성분채집술 (Apheresis)"  = "apheresis")

build_ev <- function(drugs, start_day = 0) {
  e <- NULL
  add <- function(x) if (is.null(e)) x else c(e, x)
  for (d in drugs) {
    e <- switch(d,
      statin       = add(ev(amt =   20, cmt = 12, ii = 1,   addl = 99999, time = start_day)),
      evolocumab   = add(ev(amt =  420, cmt =  9, ii = 28,  addl = 999,   time = start_day)),
      pelacarsen   = add(ev(amt =   80, cmt =  1, ii = 28,  addl = 999,   time = start_day)),
      olpasiran    = add(ev(amt =   75, cmt =  4, ii = 84,  addl = 999,   time = start_day)),
      lepodisiran  = add(ev(amt =  400, cmt =  4, ii = 168, addl = 999,   time = start_day)),
      muvalaplin   = add(ev(amt =  240, cmt =  7, ii = 1,   addl = 99999, time = start_day)),
      niacin       = add(ev(amt = 2000, cmt = 14, ii = 1,   addl = 99999, time = start_day)),
      obicetrapib  = add(ev(amt =   10, cmt = 17, ii = 1,   addl = 99999, time = start_day)),
      ziltivekimab = add(ev(amt =   30, cmt = 15, ii = 28,  addl = 999,   time = start_day)),
      e)
  }
  if (is.null(e)) ev(amt = 0, cmt = 1, time = 0) else e
}

THEME <- theme_minimal(base_size = 13) +
  theme(panel.grid.minor = element_blank(),
        plot.title = element_text(face = "bold"),
        legend.position = "bottom")

## =====================================================================
##  UI
## =====================================================================
ui <- fluidPage(
  titlePanel("고지단백(a)혈증 QSP 대시보드 — Elevated Lipoprotein(a)"),
  tags$p(style = "color:#555;margin-top:-10px",
         "하나의 입자 · 하나의 율속 단계 · 세 개의 작용 팔 · 두 개의 눈금 ",
         tags$em("(one particle, one rate-limiting step, three arms, two scales)")),

  sidebarLayout(
    sidebarPanel(
      width = 3,
      h4("환자 (Patient)"),
      sliderInput("LPA0", "기저 Lp(a) 입자 (nmol/L)", 5, 500, 250, step = 5),
      sliderInput("NKIV2", "KIV-2 반복수 (발현 우세 대립유전자)", 5, 40, 12, step = 1),
      helpText("반복수가 작을수록 분비 효율이 높아 Lp(a)가 높고, ",
               "동시에 질량 검사에서 과소 보고됩니다."),
      sliderInput("LDLP0", "기저 LDL 입자 (nmol/L)", 500, 2200, 1200, step = 50),
      sliderInput("ERISK", "혈관 위험 부담 (흡연·고혈압·당뇨)", 0, 1.5, 0, step = 0.1),
      sliderInput("IL6EXO", "동반 염증 IL-6 (pg/mL)", 0, 40, 0, step = 1),
      selectInput("FLDLRFN", "LDL 수용체 기능",
                  c("정상 (Normal)" = 1, "이형접합 FH" = 0.5, "동형접합 FH" = 0.05)),

      hr(), h4("치료 (Therapy)"),
      checkboxGroupInput("drugs", NULL, choices = DRUGS, selected = "pelacarsen"),
      sliderInput("start_yr", "치료 시작 (모델 연차)", 0, 50, 0, step = 1),
      sliderInput("horizon", "시뮬레이션 기간 (년)", 1, 62, 2, step = 1),

      hr(), h4("검사 가정 (Assay)"),
      sliderInput("NCAL", "질량 검사 캘리브레이터 아이소폼", 8, 40, 22, step = 1),
      sliderInput("RFREE", "전통적 apo(a) 검사의 유리 apo(a) 몰당 반응", 1, 6, 1, step = 0.1),
      sliderInput("FCHOL", "Lp(a) 질량 중 콜레스테롤 분율", 0.15, 0.40, 0.30, step = 0.01),
      actionButton("go", "실행 (Run)", class = "btn-primary btn-block")
    ),

    mainPanel(
      width = 9,
      tabsetPanel(
        id = "tabs",
        tabPanel("1 · 환자 프로파일",
                 br(), fluidRow(column(6, tableOutput("iso_tbl")),
                                column(6, plotOutput("iso_plot", height = 300))),
                 hr(), verbatimTextOutput("iso_note")),
        tabPanel("2 · 생산 · 조립",  br(), plotOutput("prod_plot", height = 560)),
        tabPanel("3 · 약동학 (PK)",  br(), plotOutput("pk_plot",   height = 560)),
        tabPanel("4 · Lp(a) 반응",   br(), plotOutput("pd_plot",   height = 420),
                 hr(), tableOutput("pd_tbl")),
        tabPanel("5 · ★ 측정의 함정", br(),
                 h4("같은 환자를 네 가지 방법으로 보고하면"),
                 plotOutput("assay_plot", height = 340),
                 hr(), h4("같은 TRUE 질량 50 mg/dL, 아이소폼만 다른 환자들"),
                 DTOutput("miss_tbl")),
        tabPanel("6 · 지질 · apoB",  br(), plotOutput("lipid_plot", height = 480),
                 hr(), verbatimTextOutput("lipid_note")),
        tabPanel("7 · 염증 축",      br(), plotOutput("infl_plot", height = 480),
                 hr(), verbatimTextOutput("loop_note")),
        tabPanel("8 · 혈관 · 위험",  br(), plotOutput("risk_plot", height = 480),
                 hr(), verbatimTextOutput("risk_note")),
        tabPanel("9 · 대동맥판막",   br(), plotOutput("valve_plot", height = 480),
                 hr(), verbatimTextOutput("valve_note")),
        tabPanel("10 · 시나리오 비교", br(), DTOutput("scen_tbl"))
      )
    )
  )
)

## =====================================================================
##  SERVER
## =====================================================================
server <- function(input, output, session) {

  pars <- reactive({
    list(LPA0 = input$LPA0, NKIV2 = input$NKIV2, LDLP0 = input$LDLP0,
         ERISK = input$ERISK, IL6EXO = input$IL6EXO,
         FLDLRFN = as.numeric(input$FLDLRFN), NCAL = input$NCAL,
         RFREE = input$RFREE, FCHOL = input$FCHOL,
         APHON = as.numeric("apheresis" %in% input$drugs))
  })

  sim <- eventReactive(input$go, {
    m  <- param(mod, pars())
    e  <- build_ev(setdiff(input$drugs, c("none", "apheresis")),
                   start_day = input$start_yr*YR)
    end <- max(input$horizon, input$start_yr + 1)*YR
    dl  <- if (end > 10*YR) 30 else 1
    mrgsim_df(m, events = e, end = end, delta = dl)
  }, ignoreNULL = FALSE)

  ## ---- tab 1 : patient / isoform ------------------------------------
  output$iso_tbl <- renderTable({
    d <- sim()[1, ]
    data.frame(
      항목 = c("KIV-2 반복수", "apo(a) 분자량 (kDa)", "입자 분자량 (MDa)",
               "분비 효율 SECEFF", "환산계수 (nmol/L per mg/dL)",
               "질량검사 편향 (x TRUE)", "TRUE 질량 (mg/dL)",
               "질량검사 보고값 (mg/dL)", "입자 농도 (nmol/L)"),
      값 = c(sprintf("%d", input$NKIV2),
             sprintf("%.0f", d$MWAPOAO/1000),
             sprintf("%.2f", (3.30e6 + d$MWAPOAO)/1e6),
             sprintf("%.3f", d$SECEFFO),
             sprintf("%.3f", d$CONVF),
             sprintf("%.3f", d$EPITB),
             sprintf("%.1f", d$LPA_MASS_T),
             sprintf("%.1f", d$ASSAY_MASS),
             sprintf("%.1f", d$LPA_NMOL))
    )
  })

  output$iso_plot <- renderPlot({
    nk <- 5:40
    conv <- 1e7/(3.30e6 + 14000*(nk + 10) + 35000)
    df <- data.frame(nk,
                     `환산계수 (nmol/L per mg/dL)` = conv,
                     `질량검사 편향` = (nk + 10)/(input$NCAL + 10),
                     `분비 효율` = 22^3/(22^3 + nk^3), check.names = FALSE)
    pivot_longer(df, -nk) %>%
      ggplot(aes(nk, value)) +
      geom_line(linewidth = 1.1, colour = "#0D47A1") +
      geom_vline(xintercept = input$NKIV2, linetype = 2, colour = "#B71C1C") +
      facet_wrap(~name, scales = "free_y") +
      labs(x = "KIV-2 반복수", y = NULL,
           title = "아이소폼 크기가 바꾸는 세 가지 — 그리고 그 크기 차이") + THEME
  })

  output$iso_note <- renderText({
    d <- sim()[1, ]
    paste0(
      "화학적 환산계수는 전체 아이소폼 범위에서 2.52-2.79 로 약 9% 만 변합니다.\n",
      "같은 범위에서 항체 편향은 0.56-1.41 로 150% 변합니다.\n",
      "즉 단위 환산이 아니라 항체가 문제입니다.\n\n",
      sprintf("이 환자: TRUE %.1f mg/dL 이지만 질량 검사는 %.1f mg/dL 로 보고하고,\n",
              d$LPA_MASS_T, d$ASSAY_MASS),
      sprintf("입자 농도는 %.0f nmol/L 입니다.  50 mg/dL 역치 -> %s / 125 nmol/L 역치 -> %s",
              d$LPA_NMOL,
              ifelse(d$ASSAY_MASS >= 50, "양성", "음성 (놓침)"),
              ifelse(d$LPA_NMOL  >= 125, "양성", "음성")))
  })

  ## ---- tab 2 : production / assembly chain --------------------------
  output$prod_plot <- renderPlot({
    sim() %>%
      select(time, MRNA, APOA_ER, APOA_FR, LPA_P, ASSEM, KNOCK, FMUV, TRANS) %>%
      pivot_longer(-time) %>%
      mutate(name = factor(name,
        levels = c("TRANS","MRNA","KNOCK","APOA_ER","APOA_FR","FMUV","ASSEM","LPA_P"),
        labels = c("전사 TRANS","LPA mRNA","mRNA 분해 배수 KNOCK","ER apo(a) 풀",
                   "유리 apo(a)","조립 잔여분율 FMUV","조립 flux ASSEM","Lp(a) 입자"))) %>%
      ggplot(aes(time/YR, value)) +
      geom_line(linewidth = 1, colour = "#00695C") +
      facet_wrap(~name, scales = "free_y", ncol = 3) +
      labs(x = "시간 (년)", y = NULL,
           title = "생산 사슬을 따라간 약물 작용 지점") + THEME
  })

  ## ---- tab 3 : PK ----------------------------------------------------
  output$pk_plot <- renderPlot({
    sim() %>%
      select(time, CPELL, CMUV, CEVO, CSTA, CZIL, COBI, SIR_RISC) %>%
      pivot_longer(-time) %>% filter(value > 0) %>%
      mutate(name = recode(name,
        CPELL = "펠라카르센 간 (mg/L)", CMUV = "무발라플린 혈장 (mg/L)",
        CEVO = "에볼로쿠맙 혈장 (mg/L)", CSTA = "스타틴 혈장 (mg/L)",
        CZIL = "질티베키맙 혈장 (mg/L)", COBI = "오비세트라핍 혈장 (mg/L)",
        SIR_RISC = "siRNA RISC 적재량 (mg)")) %>%
      ggplot(aes(time/YR, value)) +
      geom_line(linewidth = 0.9, colour = "#1B5E20") +
      facet_wrap(~name, scales = "free_y", ncol = 2) +
      labs(x = "시간 (년)", y = NULL,
           title = "약동학 — RISC 구획의 반감기가 PD 지속을 만든다") + THEME
  })

  ## ---- tab 4 : Lp(a) response ---------------------------------------
  output$pd_plot <- renderPlot({
    sim() %>%
      ggplot(aes(time/YR, LPA_NMOL)) +
      geom_hline(yintercept = c(75, 125, 180), linetype = 3, colour = "grey50") +
      annotate("text", x = 0, y = c(75, 125, 180), hjust = 0, vjust = -0.4,
               size = 3, colour = "grey40",
               label = c("75 정상 상한", "125 위험", "180 매우 높음")) +
      geom_line(linewidth = 1.2, colour = "#B71C1C") +
      labs(x = "시간 (년)", y = "Lp(a) 입자 (nmol/L)",
           title = "Lp(a) 입자 농도") + THEME
  })

  output$pd_tbl <- renderTable({
    d <- sim(); b <- d[1, ]; e <- d[nrow(d), ]
    data.frame(
      지표 = c("Lp(a) (nmol/L)", "Lp(a) TRUE 질량 (mg/dL)",
               "질량 검사 보고값 (mg/dL)", "전통 apo(a) 검사 (nmol/L)",
               "유리 apo(a) (nmol/L)", "절대 감소량 (nmol/L)"),
      기저 = c(b$LPA_NMOL, b$LPA_MASS_T, b$ASSAY_MASS, b$ASSAY_APOA,
               b$ASSAY_APOA - b$LPA_NMOL, NA),
      종료 = c(e$LPA_NMOL, e$LPA_MASS_T, e$ASSAY_MASS, e$ASSAY_APOA,
               e$ASSAY_APOA - e$LPA_NMOL, NA),
      변화율 = c(100*(e$LPA_NMOL/b$LPA_NMOL - 1),
                 100*(e$LPA_MASS_T/b$LPA_MASS_T - 1),
                 100*(e$ASSAY_MASS/b$ASSAY_MASS - 1),
                 100*(e$ASSAY_APOA/b$ASSAY_APOA - 1),
                 100*((e$ASSAY_APOA - e$LPA_NMOL)/(b$ASSAY_APOA - b$LPA_NMOL) - 1),
                 e$LPA_NMOL - b$LPA_NMOL)
    )
  }, digits = 1)

  ## ---- tab 5 : the measurement trap ---------------------------------
  output$assay_plot <- renderPlot({
    sim() %>%
      transmute(time,
                `입자 몰농도 검사 (nmol/L)` = ASSAY_NMOL,
                `intact-Lp(a) 검사 (nmol/L)` = ASSAY_INTACT,
                `전통 apo(a) 검사 (nmol/L)` = ASSAY_APOA,
                `질량 검사 x 환산 (nmol/L 환산)` = ASSAY_MASS*CONVF) %>%
      pivot_longer(-time) %>%
      ggplot(aes(time/YR, value, colour = name)) +
      geom_line(linewidth = 1.1) +
      scale_colour_manual(values = c("#FF6F00", "#00838F", "#6A1B9A", "#B71C1C")) +
      labs(x = "시간 (년)", y = "보고값 (nmol/L 기준)", colour = NULL,
           title = "네 개의 검사, 하나의 환자 — 조립 억제제에서 가장 크게 갈린다") +
      THEME
  })

  output$miss_tbl <- renderDT({
    nk <- c(6, 8, 12, 16, 22, 28, 35, 40)
    conv <- 1e7/(3.30e6 + 14000*(nk + 10) + 35000)
    true_nm <- 50*conv
    bias <- (nk + 10)/(input$NCAL + 10)
    data.frame(
      `KIV-2` = nk,
      `apo(a) MW (kDa)` = round((14000*(nk + 10) + 35000)/1000),
      `TRUE 질량 (mg/dL)` = 50,
      `TRUE 입자 (nmol/L)` = round(true_nm, 1),
      `환산계수` = round(conv, 3),
      `항체 편향` = round(bias, 3),
      `보고 질량 (mg/dL)` = round(50*bias, 1),
      `질량 50 역치` = ifelse(50*bias >= 50, "양성", "음성 (놓침)"),
      `몰 125 역치` = ifelse(true_nm >= 125, "양성", "음성"),
      check.names = FALSE
    )
  }, options = list(dom = "t", pageLength = 10), rownames = FALSE)

  ## ---- tab 6 : lipids / apoB ----------------------------------------
  output$lipid_plot <- renderPlot({
    sim() %>%
      select(time, LDLC_TRUE, LDLC_MEAS, LDLC_CORR, LPAC, APOB_TOT, LPA_APOB_PCT) %>%
      pivot_longer(-time) %>%
      mutate(name = recode(name,
        LDLC_TRUE = "TRUE LDL-C (mg/dL)", LDLC_MEAS = "보고된 LDL-C (mg/dL)",
        LDLC_CORR = "Dahlen 보정 LDL-C (mg/dL)", LPAC = "Lp(a)-C (mg/dL)",
        APOB_TOT = "총 apoB (mg/dL)", LPA_APOB_PCT = "apoB 입자 중 Lp(a) 비율 (%)")) %>%
      ggplot(aes(time/YR, value)) +
      geom_line(linewidth = 1, colour = "#AD1457") +
      facet_wrap(~name, scales = "free_y", ncol = 3) +
      labs(x = "시간 (년)", y = NULL,
           title = "'LDL-C' 안에 들어 있는 것") + THEME
  })

  output$lipid_note <- renderText({
    b <- sim()[1, ]
    paste0(
      sprintf("보고되는 LDL-C %.1f mg/dL 중 %.1f mg/dL (%.0f%%) 은 LDL 이 아니라 Lp(a) 콜레스테롤입니다.\n",
              b$LDLC_MEAS, b$LPAC, 100*b$LPAC/b$LDLC_MEAS),
      sprintf("Dahlen 보정을 (편향된) 질량 검사값으로 수행하면 %.1f mg/dL 이 되어\n",
              b$LDLC_CORR),
      sprintf("TRUE LDL-C %.1f mg/dL 대비 %+.1f mg/dL 만큼 어긋납니다.\n", b$LDLC_TRUE,
              b$LDLC_CORR - b$LDLC_TRUE),
      sprintf("Lp(a) 는 apoB 입자의 %.1f%% 에 불과하지만 내막 저류 flux 의 %.1f%% 를 차지합니다.",
              b$LPA_APOB_PCT, b$RETSHARE))
  })

  ## ---- tab 7 : inflammation -----------------------------------------
  output$infl_plot <- renderPlot({
    sim() %>%
      select(time, OXPL, MONO, IL6, IL6EFF, CRP, FIL6) %>%
      pivot_longer(-time) %>%
      mutate(name = recode(name,
        OXPL = "OxPL-apoB (상대값)", MONO = "훈련된 단핵구 지수",
        IL6 = "IL-6 (pg/mL)", IL6EFF = "신호전달 가능 IL-6 (pg/mL)",
        CRP = "hsCRP (mg/L)", FIL6 = "LPA 전사에 대한 IL-6 배수")) %>%
      ggplot(aes(time/YR, value)) +
      geom_line(linewidth = 1, colour = "#E65100") +
      facet_wrap(~name, scales = "free_y", ncol = 3) +
      labs(x = "시간 (년)", y = NULL,
           title = "OxPL -> 단핵구 -> IL-6 -> LPA 프로모터") + THEME
  })

  output$loop_note <- renderText({
    paste0(
      "피드포워드 루프는 실재하지만 약합니다. 개방 루프 이득 g 는 정상 염증 상태에서\n",
      "약 0.02 이며, 증폭은 1/(1-g) = 1.02 에 불과합니다.\n\n",
      "약하지 않은 것은 '기저 IL-6 톤' 입니다: 기저 IL-6 만으로 LPA 전사의 약 25% 가\n",
      "설명되며, 그래서 IL-6 차단이 비염증 환자에서 Lp(a) 를 약 25% 낮추고\n",
      "류마티스 관절염에서는 약 35-40% 낮춥니다 — 같은 방정식, 바뀐 것은 IL-6 뿐입니다.")
  })

  ## ---- tab 8 : plaque / risk ----------------------------------------
  output$risk_plot <- renderPlot({
    sim() %>%
      select(time, INT_LPA, FOAM, NECRO, PLAQUE, VULN, HR_MACE, HR_SLOW, HR_FAST) %>%
      pivot_longer(-time) %>%
      mutate(name = recode(name,
        INT_LPA = "내막 저류 Lp(a)", FOAM = "거품세포 부담",
        NECRO = "괴사핵 지수", PLAQUE = "죽종 용적 PAV (%)",
        VULN = "취약성 지수", HR_MACE = "MACE 위험비 HR",
        HR_SLOW = "HR — 느린 성분 (죽종)", HR_FAST = "HR — 빠른 성분 (취약성)")) %>%
      ggplot(aes(time/YR, value)) +
      geom_line(linewidth = 1, colour = "#3E2723") +
      facet_wrap(~name, scales = "free_y", ncol = 4) +
      labs(x = "시간 (년)", y = NULL,
           title = "느린 성분과 빠른 성분 — MR 과 임상시험이 갈리는 이유") + THEME
  })

  output$risk_note <- renderText({
    paste0(
      "위험은 시간상수가 다른 두 성분이 나릅니다.\n",
      "  느림 : 죽종 부담, 반감기 약 14년 → 초과 위험의 약 70%. 멘델 무작위화가 보는 것.\n",
      "  빠름 : OxPL/IL-6 매개 취약성, 반감기 약 2개월 → 약 30%. 5년 시험이 움직일 수 있는 것.\n\n",
      "그 결과 평생 노출 효과 / 5년 시험 효과 비율이 약 2.6-2.7배로 '자동으로' 나옵니다.\n",
      "이 비율은 보정계수로 넣은 값이 아니라 두 시간상수에서 유도된 값입니다.")
  })

  ## ---- tab 9 : aortic valve -----------------------------------------
  output$valve_plot <- renderPlot({
    sim() %>%
      select(time, ATXV, LYSOPA, VIC_OST, VCALC, AVA, MGRAD) %>%
      pivot_longer(-time) %>%
      mutate(name = recode(name,
        ATXV = "판막 오토탁신", LYSOPA = "판막 LysoPA",
        VIC_OST = "골형성 VIC 분율", VCALC = "판막 칼슘 AVC (AU)",
        AVA = "대동맥판 면적 (cm2)", MGRAD = "평균 압력차 (mmHg)")) %>%
      ggplot(aes(time/YR, value)) +
      geom_line(linewidth = 1, colour = "#283593") +
      facet_wrap(~name, scales = "free_y", ncol = 3) +
      labs(x = "시간 (년)", y = NULL,
           title = "판막 — 시작 시점이 약효보다 중요한 유일한 엔드포인트") + THEME
  })

  output$valve_note <- renderText({
    paste0(
      "판막 석회화는 문턱(Hill n=4)을 넘으면 Lp(a) 와 무관하게 자기영속합니다.\n",
      "62년 시뮬레이션에서 (Lp(a) 250 nmol/L, 펠라카르센):\n",
      "  30세 시작 → AVC 약 206 AU, AVA 2.45 cm2   (칼슘의 87% 예방)\n",
      "  60세 시작 → AVC 약 1327 AU, AVA 0.93 cm2  (20년 치료로 칼슘 15% 감소)\n",
      "  무치료    → AVC 약 1558 AU, AVA 0.82 cm2\n\n",
      "SEAS · ASTRONOMER · SALTIRE 가 모두 음성이었던 이유이며, 동시에\n",
      "판막 적응증은 아무리 강력한 약이라도 이차예방 설계로는 성공할 수 없다는 예측입니다.")
  })

  ## ---- tab 10 : scenario comparison ---------------------------------
  output$scen_tbl <- renderDT({
    withProgress(message = "20개 시나리오 실행 중...", value = 0, {
      rows <- lapply(seq_along(sc), function(i) {
        incProgress(1/length(sc))
        s <- sc[[i]]
        m <- if (length(s$p)) param(mod, s$p) else mod
        o <- mrgsim_df(m, events = s$ev, end = 2*YR, delta = 7)
        b <- o[1, ]; e <- o[nrow(o), ]
        data.frame(
          시나리오 = names(sc)[i],
          `기저 nmol/L` = round(b$LPA_NMOL, 0),
          `2년 nmol/L` = round(e$LPA_NMOL, 1),
          `변화 %` = round(100*(e$LPA_NMOL/b$LPA_NMOL - 1), 1),
          `절대 Δ nmol/L` = round(e$LPA_NMOL - b$LPA_NMOL, 1),
          `LDL-C(보고) %` = round(100*(e$LDLC_MEAS/b$LDLC_MEAS - 1), 1),
          `apoB %` = round(100*(e$APOB_TOT/b$APOB_TOT - 1), 1),
          `hsCRP %` = round(100*(e$CRP/b$CRP - 1), 1),
          check.names = FALSE)
      })
      do.call(rbind, rows)
    })
  }, options = list(dom = "t", pageLength = 25, scrollX = TRUE), rownames = FALSE)
}

shinyApp(ui, server)
