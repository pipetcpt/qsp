## =============================================================================
## Anthracycline-Induced Cardiotoxicity (AIC) — Shiny dashboard
## 안트라사이클린 심장독성 QSP 모델 인터랙티브 대시보드
##
## Eight tabs, each built around one thing the model claims:
##   1. 환자 프로파일 (Patient)      risk phenotype -> which parameters move
##   2. PK · 심근 노출 (Exposure)    the two cardiac pools and their metrics
##   3. 손상 축 (Injury arms)        Top2b/DSB/p53 vs iron/ROS/mitochondria
##   4. 두 결손과 가면 (Deficits)     MYO vs FUNC, and LVEF with/without the mask
##   5. 임상 엔드포인트 (Endpoints)   LVEF, GLS, CTRCD, NYHA
##   6. 시나리오 비교 (Scenarios)     same cumulative dose, different schedules
##   7. 바이오마커 (Biomarkers)       troponin/NT-proBNP lead times
##   8. 가역성 창 (Window)            when HF therapy is started
##
## Run:  shiny::runApp("aic_shiny_app.R")
## Requires: shiny, mrgsolve, dplyr, tidyr, ggplot2, DT
## =============================================================================

library(shiny)
library(mrgsolve)
library(dplyr)
library(tidyr)
library(ggplot2)
library(DT)

MODEL_FILE <- "aic_mrgsolve_model.R"
mod <- mread_cache("aic", MODEL_FILE)

BSA_DEFAULT <- 1.8
WT_DEFAULT  <- 70

## ---------------------------------------------------------------------------
## dosing helpers (mirror the R driver in aic_mrgsolve_model.R)
## ---------------------------------------------------------------------------
dox_ev <- function(dose_mg_m2, ncyc, interval, bsa, tinf = 0, lipo = FALSE,
                   start = 0) {
  ev(amt = dose_mg_m2 * bsa, cmt = if (lipo) 4 else 1, time = start,
     ii = interval, addl = ncyc - 1,
     rate = if (tinf > 0) dose_mg_m2 * bsa / tinf else 0)
}
dex_ev <- function(dose_mg_m2, ncyc, interval, bsa, ratio = 10) {
  ev(amt = ratio * dose_mg_m2 * bsa, cmt = 10, time = 0, ii = interval,
     addl = ncyc - 1)
}
tras_ev <- function(start, wt, weeks = 52) {
  c(ev(amt = 8 * wt, cmt = 12, time = start),
    ev(amt = 6 * wt, cmt = 12, time = start + 21, ii = 21,
       addl = floor(weeks / 3) - 2))
}

SCHEDULES <- c("IV push (q3w)"          = "push",
               "72-h continuous infusion" = "inf72",
               "Weekly fractionated"    = "weekly",
               "Pegylated liposomal (PLD)" = "pld")

theme_aic <- function() {
  theme_minimal(base_size = 13) +
    theme(panel.grid.minor = element_blank(),
          legend.position = "bottom",
          plot.title = element_text(face = "bold"))
}

## ---------------------------------------------------------------------------
## UI
## ---------------------------------------------------------------------------
ui <- fluidPage(
  titlePanel(paste("안트라사이클린 심장독성 QSP 모델 —",
                   "Anthracycline-Induced Cardiotoxicity")),
  helpText(paste("누적용량은 노출 지표가 아니다: 같은 mg/m² 라도 첨두농도를",
                 "바꾸면 손상이 바뀐다. LVEF는 보상성 비대에 가려진",
                 "후행지표이며, 회복 가능성은 섬유화 시계가 결정한다.")),

  sidebarLayout(
    sidebarPanel(
      width = 3,
      h4("환자 (Patient)"),
      sliderInput("age", "나이 (years)", 20, 85, 55, 1),
      sliderInput("bsa", "체표면적 BSA (m²)", 1.3, 2.3, BSA_DEFAULT, 0.05),
      sliderInput("ef0", "기저 LVEF (%)", 45, 72, 62, 1),
      checkboxInput("htn", "고혈압 (hypertension)", FALSE),
      checkboxInput("rt", "종격부 방사선 병력 (prior mediastinal RT)", FALSE),
      selectInput("cbr1", "CBR1/AKR1C3 대사형 (doxorubicinol FM)",
                  c("Slow (0.12)" = 0.12, "Intermediate (0.18)" = 0.18,
                    "Typical (0.25)" = 0.25, "Fast (0.34)" = 0.34,
                    "Very fast (0.45)" = 0.45), selected = 0.25),

      hr(), h4("항암 요법 (Chemotherapy)"),
      selectInput("sched", "스케줄 · 제형", SCHEDULES, selected = "push"),
      sliderInput("dose", "회당 용량 (mg/m²)", 10, 90, 60, 5),
      sliderInput("ncyc", "주기 수 (cycles)", 1, 24, 8, 1),
      sliderInput("ii", "주기 간격 (days)", 7, 28, 21, 7),
      helpText(textOutput("cumtxt")),

      hr(), h4("심장보호 (Cardioprotection)"),
      checkboxInput("dex", "덱스라족산 10:1 (dexrazoxane)", FALSE),
      checkboxInput("sta", "아토르바스타틴 40 mg", FALSE),
      sliderInput("sta_day", "스타틴 시작일 (day)", 0, 400, 0, 10),
      checkboxInput("ace", "ACE 억제제 / ARB", FALSE),
      checkboxInput("bb", "베타차단제 (카르베딜롤)", FALSE),
      checkboxInput("arni", "사쿠비트릴/발사르탄 (ARNI)", FALSE),
      checkboxInput("sglt", "SGLT2 억제제", FALSE),
      sliderInput("hf_day", "ACEi/BB/ARNI 시작일 (day)", 0, 540, 0, 30),

      hr(), h4("HER2 표적 (Trastuzumab)"),
      selectInput("tras", "트라스투주맙",
                  c("없음 (none)" = "none",
                    "순차 (sequential, day 84)" = "seq",
                    "동시 (concurrent, day 0)" = "conc"), selected = "none"),

      hr(), h4("시뮬레이션"),
      sliderInput("tend", "관찰기간 (days)", 365, 1460, 730, 30),
      sliderInput("nid", "가상 환자 수 (virtual subjects)", 1, 500, 100, 1),
      actionButton("go", "실행 (Run)", class = "btn-primary")
    ),

    mainPanel(
      width = 9,
      tabsetPanel(
        id = "tabs",
        tabPanel(
          "1. 환자 프로파일",
          h4("위험 표현형이 어떤 파라미터를 움직이는가"),
          tableOutput("risk_tbl"),
          helpText(paste("고령·고혈압·방사선 병력은 별개의 '위험점수'가 아니라",
                         "근세포 사멸 감수성(KDROS)·섬유화 성향(KFIBIN)·재생능",
                         "(KREG)·기저 예비력(EF0)을 각각 다르게 움직인다.")),
          h4("보상 예비력 (masking reserve)"),
          plotOutput("reserve_plot", height = "300px"),
          helpText(paste("HYPMAX = 0.32 이므로 유효 수축단위 결손이 약 20%를",
                         "넘어서면 보상이 포화되고 LVEF가 급격히 떨어진다.",
                         "기저 LVEF가 낮으면 이 여유가 작다."))
        ),

        tabPanel(
          "2. PK · 심근 노출",
          h4("혈장 · 대사체 농도"),
          plotOutput("pk_plot", height = "260px"),
          h4("심근 내 두 개의 풀 — 서로 다른 노출 지표"),
          plotOutput("cardiac_plot", height = "300px"),
          tableOutput("exposure_tbl"),
          helpText(paste("빠른 핵 풀(CHF)은 첨두를 추종하고 느린 잔류 풀(CH,",
                         "CHM)은 AUC를 적분한다. 스케줄을 바꾸면 앞의 것만",
                         "바뀐다 — 이것이 지속주입·리포조말이 작동하는 이유다."))
        ),

        tabPanel(
          "3. 손상 축",
          fluidRow(
            column(6, h4("Top2b / DSB / p53 (첨두 구동)"),
                   plotOutput("arm_nuc", height = "300px")),
            column(6, h4("철 · ROS · 미토콘드리아 (AUC 구동)"),
                   plotOutput("arm_redox", height = "300px"))),
          h4("Top2b 단백 수준 (덱스라족산 효과)"),
          plotOutput("t2b_plot", height = "220px"),
          helpText(paste("두 축은 사멸률에서 곱으로 결합한다:",
                         "kdeath ∝ Hill(ROS;n=3) × (1 + 1.2·p53).",
                         "따라서 한쪽만 제거해도 이득이 지분보다 크다."))
        ),

        tabPanel(
          "4. 두 결손과 가면",
          h4("가역적 결손 FUNC vs 비가역적 근세포 소실"),
          plotOutput("deficit_plot", height = "300px"),
          h4("가면의 크기: 보상성 비대를 제거하면 LVEF는 어디에 있었을까"),
          plotOutput("mask_plot", height = "300px"),
          tableOutput("mask_tbl"),
          helpText(paste("LVEF(실제)와 LVEF_NOMASK의 차이가 보상성 비대가",
                         "가려주고 있는 점수이고, LVEF_IRREV와의 차이가 아직",
                         "되찾을 수 있는 가역 성분이다."))
        ),

        tabPanel(
          "5. 임상 엔드포인트",
          fluidRow(
            column(6, h4("LVEF (%)"), plotOutput("ef_plot", height = "300px")),
            column(6, h4("|GLS| (%)"), plotOutput("gls_plot", height = "300px"))),
          h4("사건 발생률 (virtual population)"),
          DTOutput("endpoint_tbl"),
          helpText(paste("CTRCD = LVEF 10점 이상 하락 & <50% (ASE/ESC).",
                         "GLS 상대 15% 감소는 같은 손상을 더 이른 시점에",
                         "포착한다."))
        ),

        tabPanel(
          "6. 시나리오 비교",
          h4("동일 누적용량, 다른 스케줄 · 제형 · 보호제"),
          actionButton("go_cmp", "시나리오 비교 실행", class = "btn-success"),
          plotOutput("cmp_plot", height = "340px"),
          DTOutput("cmp_tbl"),
          helpText(paste("모든 행의 누적용량은 동일하다. 차이는 전부",
                         "첨두농도와 Top2b 축에서 온다 — 항종양 AUC는",
                         "유지된다는 점이 치료지수 개선의 근거다."))
        ),

        tabPanel(
          "7. 바이오마커",
          h4("hs-cTnI — 손상률(rate)의 지표"),
          plotOutput("tni_plot", height = "280px"),
          h4("NT-proBNP — 벽응력의 지표"),
          plotOutput("bnp_plot", height = "250px"),
          h4("선행시간 (lead time)"),
          tableOutput("lead_tbl"),
          helpText(paste("트로포닌 반감기는 약 12시간이므로 누적 손상이 아니라",
                         "'지금 죽고 있는 속도'를 읽는다. 사이클마다 음성인",
                         "환자에서 CTRCD가 드물다는 점(높은 음성예측도)이",
                         "트로포닌 유도 전략의 근거다."))
        ),

        tabPanel(
          "8. 가역성 창",
          h4("심부전 치료를 언제 시작하는가"),
          actionButton("go_win", "가역성 창 계산", class = "btn-warning"),
          plotOutput("win_plot", height = "320px"),
          DTOutput("win_tbl"),
          helpText(paste("회복 가능성은 시간 자체가 아니라 그 시점의 섬유화",
                         "수준을 따라간다: FUNC 회복속도가 1/(1+3·FIB)로",
                         "느려지기 때문이다. 창을 닫는 것은 섬유화 시계다."))
        )
      )
    )
  )
)

## ---------------------------------------------------------------------------
## server
## ---------------------------------------------------------------------------
server <- function(input, output, session) {

  output$cumtxt <- renderText({
    sprintf("누적용량 = %d mg/m²  (총 %.0f mg, BSA %.2f m²)",
            input$dose * input$ncyc, input$dose * input$ncyc * input$bsa,
            input$bsa)
  })

  ## risk phenotype -> parameter multipliers
  risk_par <- reactive({
    p <- list(EF0 = input$ef0, FM = as.numeric(input$cbr1))
    kd <- 1; kf <- 1; kfib <- 1; kreg <- 1
    if (input$age >= 70) { kd <- kd * 1.6; kreg <- kreg * 0.4 }
    else if (input$age >= 60) { kd <- kd * 1.2; kreg <- kreg * 0.7 }
    if (input$age < 18) { kreg <- kreg * 1.4 }
    if (input$htn) { kfib <- kfib * 1.3; kf <- kf * 1.2 }
    if (input$rt)  { kfib <- kfib * 1.6; kd <- kd * 1.2 }
    p$KDROS  <- 4.8e-4 * kd
    p$KFIN   <- 1.0e-3 * kf
    p$KFIBIN <- 8.0e-3 * kfib
    p$KREG   <- 2.0e-5 * kreg
    p
  })

  drug_par <- reactive({
    p <- list()
    p$TGT_STA  <- if (input$sta)  0.75 else 0
    p$TON_STA  <- input$sta_day
    p$TGT_ACE  <- if (input$ace)  0.80 else 0
    p$TGT_BB   <- if (input$bb)   0.80 else 0
    p$TGT_ARNI <- if (input$arni) 0.90 else 0
    p$TGT_SGLT <- if (input$sglt) 0.70 else 0
    p$TON_ACE  <- input$hf_day
    p$TON_BB   <- input$hf_day
    p$TON_ARNI <- if (input$arni) input$hf_day else 1e6
    p$TON_SGLT <- if (input$sglt) 0 else 1e6
    p
  })

  regimen <- reactive({
    s <- input$sched
    e <- switch(
      s,
      push   = dox_ev(input$dose, input$ncyc, input$ii, input$bsa),
      inf72  = dox_ev(input$dose, input$ncyc, input$ii, input$bsa, tinf = 3),
      weekly = dox_ev(input$dose, input$ncyc, 7, input$bsa),
      pld    = dox_ev(input$dose, input$ncyc, 28, input$bsa,
                      lipo = TRUE))
    if (input$dex && s != "pld")
      e <- c(e, dex_ev(input$dose, input$ncyc, input$ii, input$bsa))
    if (input$tras == "seq")  e <- c(e, tras_ev(84, WT_DEFAULT))
    if (input$tras == "conc") e <- c(e, tras_ev(0,  WT_DEFAULT))
    e
  })

  sim <- eventReactive(input$go, {
    m <- param(mod, c(risk_par(), drug_par()))
    n <- input$nid
    out <- m %>% ev(regimen()) %>%
      mrgsim(end = input$tend, delta = 1, nid = n) %>% as_tibble()
    out
  }, ignoreNULL = FALSE)

  idx <- reactive({
    m <- param(mod, c(risk_par(), drug_par()))
    m %>% ev(regimen()) %>% mrgsim(end = input$tend, delta = 1) %>% as_tibble()
  })

  ## ---- tab 1 ------------------------------------------------------------
  output$risk_tbl <- renderTable({
    p <- risk_par()
    tibble(파라미터 = c("KDROS (사멸 감수성)", "KFIN (기능손상 감수성)",
                        "KFIBIN (섬유화 성향)", "KREG (근세포 재생)",
                        "EF0 (기저 LVEF)", "FM (독소루비시놀 생성분율)"),
           기준값 = c("4.8e-4", "1.0e-3", "8.0e-3", "2.0e-5", "62", "0.25"),
           현재값 = c(sprintf("%.2e", p$KDROS), sprintf("%.2e", p$KFIN),
                      sprintf("%.2e", p$KFIBIN), sprintf("%.2e", p$KREG),
                      sprintf("%.0f", p$EF0), sprintf("%.2f", p$FM)))
  })

  output$reserve_plot <- renderPlot({
    d <- tibble(deficit = seq(0, 0.5, 0.005)) %>%
      mutate(HYP = pmin(0.32, 1.6 * deficit),
             CONT = (1 - deficit) * (1 + HYP),
             LVEF = input$ef0 * CONT^0.9,
             LVEF_nomask = input$ef0 * (1 - deficit)^0.9)
    ggplot(d, aes(100 * deficit)) +
      geom_line(aes(y = LVEF, colour = "보상 후 (관측되는 LVEF)"), linewidth = 1.2) +
      geom_line(aes(y = LVEF_nomask, colour = "보상 없다면"),
                linewidth = 1.0, linetype = 2) +
      geom_vline(xintercept = 100 * 0.32 / 1.6, linetype = 3) +
      annotate("text", x = 20.5, y = input$ef0 * 0.55,
               label = "예비력 포화 (HYPMAX)", hjust = 0, size = 4) +
      labs(x = "유효 수축단위 결손 (%)", y = "LVEF (%)", colour = NULL,
           title = "가면(mask)의 대수적 형태") +
      theme_aic()
  })

  ## ---- tab 2 ------------------------------------------------------------
  output$pk_plot <- renderPlot({
    d <- idx() %>% select(time, CPLASMA, CDOXOL) %>%
      pivot_longer(-time)
    ggplot(d, aes(time, value, colour = name)) + geom_line(linewidth = 0.8) +
      scale_y_log10() +
      labs(x = "day", y = "concentration (mg/L, log)", colour = NULL,
           title = "독소루비신 · 독소루비시놀") + theme_aic()
  })

  output$cardiac_plot <- renderPlot({
    d <- idx() %>% select(time, CHF, CH, CHM) %>% pivot_longer(-time)
    lab <- c(CHF = "CHF 빠른 핵 풀 (첨두 추종)",
             CH = "CH 느린 잔류 풀 (AUC 적분)",
             CHM = "CHM 심근 독소루비시놀 (축적)")
    ggplot(d, aes(time, value, colour = lab[name])) +
      geom_line(linewidth = 0.9) +
      facet_wrap(~lab[name], scales = "free_y", ncol = 1) +
      labs(x = "day", y = "cardiac pool", colour = NULL) +
      theme_aic() + theme(legend.position = "none")
  })

  output$exposure_tbl <- renderTable({
    d <- idx()
    tibble(지표 = c("peak CHF (핵 첨두)", "peak CH (잔류)",
                    "peak CHM (대사체)", "AUC_CH (누적 잔류 노출)",
                    "peak p53 (유전독성 기억)"),
           값 = sprintf("%.3f", c(max(d$CHF), max(d$CH), max(d$CHM),
                                  max(d$AUCH), max(d$P53))))
  })

  ## ---- tab 3 ------------------------------------------------------------
  output$arm_nuc <- renderPlot({
    d <- idx() %>% select(time, DSB, P53) %>% pivot_longer(-time)
    ggplot(d, aes(time, value, colour = name)) + geom_line(linewidth = 0.9) +
      labs(x = "day", y = "level", colour = NULL) + theme_aic()
  })
  output$arm_redox <- renderPlot({
    d <- idx() %>% select(time, ROS, MITOD, LIP) %>% pivot_longer(-time)
    ggplot(d, aes(time, value, colour = name)) + geom_line(linewidth = 0.9) +
      labs(x = "day", y = "level", colour = NULL) + theme_aic()
  })
  output$t2b_plot <- renderPlot({
    ggplot(idx(), aes(time, TOP2B)) + geom_line(linewidth = 0.9) +
      ylim(0, 1.05) +
      labs(x = "day", y = "Top2b (1 = 정상)",
           title = "덱스라족산은 Top2b 단백을 분해한다 (재합성 t½ 2 d)") +
      theme_aic()
  })

  ## ---- tab 4 ------------------------------------------------------------
  output$deficit_plot <- renderPlot({
    d <- idx() %>%
      transmute(time, `비가역 근세포 소실 (1-MYO)` = 1 - MYOCYTES,
                `가역 기능 결손 FUNC` = FUNCDEF,
                `섬유화 FIB` = FIBROSIS, `보상성 비대 HYP` = MASK) %>%
      pivot_longer(-time)
    ggplot(d, aes(time, value, colour = name)) + geom_line(linewidth = 1) +
      labs(x = "day", y = "fraction", colour = NULL) + theme_aic()
  })

  output$mask_plot <- renderPlot({
    d <- idx() %>% select(time, LVEF, LVEF_NOMASK, LVEF_IRREV) %>%
      pivot_longer(-time)
    lab <- c(LVEF = "LVEF (관측)", LVEF_NOMASK = "보상성 비대 제거 시",
             LVEF_IRREV = "가역 성분 회복 시 (상한)")
    ggplot(d, aes(time, value, colour = lab[name])) +
      geom_line(linewidth = 1) + geom_hline(yintercept = 50, linetype = 3) +
      labs(x = "day", y = "LVEF (%)", colour = NULL) + theme_aic()
  })

  output$mask_tbl <- renderTable({
    idx() %>% filter(time %in% c(0, 60, 120, 170, 240, 365, 540, 730)) %>%
      transmute(day = time, LVEF = round(LVEF, 2),
                `가면 점수` = round(LVEF - LVEF_NOMASK, 2),
                `가역 성분 점수` = round(LVEF_IRREV - LVEF, 2),
                `근세포 (%)` = round(100 * MYOCYTES, 1),
                `FIB` = round(FIBROSIS, 3))
  })

  ## ---- tab 5 ------------------------------------------------------------
  output$ef_plot <- renderPlot({
    d <- sim()
    q <- d %>% group_by(time) %>%
      summarise(med = median(LVEF), lo = quantile(LVEF, .05),
                hi = quantile(LVEF, .95), .groups = "drop")
    ggplot(q, aes(time)) +
      geom_ribbon(aes(ymin = lo, ymax = hi), alpha = .2) +
      geom_line(aes(y = med), linewidth = 1.1) +
      geom_hline(yintercept = 50, linetype = 3) +
      geom_hline(yintercept = 40, linetype = 2, colour = "red") +
      labs(x = "day", y = "LVEF (%)",
           title = "중앙값과 90% 구간") + theme_aic()
  })

  output$gls_plot <- renderPlot({
    d <- sim()
    q <- d %>% group_by(time) %>%
      summarise(med = median(GLS), lo = quantile(GLS, .05),
                hi = quantile(GLS, .95), .groups = "drop")
    ggplot(q, aes(time)) +
      geom_ribbon(aes(ymin = lo, ymax = hi), alpha = .2) +
      geom_line(aes(y = med), linewidth = 1.1) +
      labs(x = "day", y = "|GLS| (%)",
           title = "GLS는 비대 보상을 받지 않는다") + theme_aic()
  })

  output$endpoint_tbl <- renderDT({
    d <- sim()
    s <- d %>% group_by(ID) %>%
      summarise(EF0 = first(LVEF), nadir = min(LVEF),
                CTRCD = any(LVEF < 50 & (first(LVEF) - LVEF) >= 10),
                HF = any(LVEF < 40),
                GLS15 = any((first(GLS) - GLS) / first(GLS) >= 0.15),
                .groups = "drop")
    datatable(tibble(
      지표 = c("CTRCD (LVEF ≥10점↓ & <50%)", "증상성 심부전 대리 (LVEF <40%)",
               "GLS 상대 15% 이상 감소", "평균 LVEF 최저치",
               "평균 LVEF 하락 (최저)"),
      값 = c(sprintf("%.1f %%", 100 * mean(s$CTRCD)),
             sprintf("%.1f %%", 100 * mean(s$HF)),
             sprintf("%.1f %%", 100 * mean(s$GLS15)),
             sprintf("%.2f %%", mean(s$nadir)),
             sprintf("%.2f points", mean(s$EF0 - s$nadir)))),
      options = list(dom = "t"), rownames = FALSE)
  })

  ## ---- tab 6 ------------------------------------------------------------
  cmp <- eventReactive(input$go_cmp, {
    cum <- input$dose * input$ncyc
    base <- param(mod, c(risk_par(), list(TGT_STA = 0, TGT_ACE = 0,
                                          TGT_BB = 0, TGT_ARNI = 0,
                                          TGT_SGLT = 0)))
    defs <- list(
      list(k = "IV push", e = dox_ev(60, round(cum / 60), 21, input$bsa), p = list()),
      list(k = "72-h infusion",
           e = dox_ev(60, round(cum / 60), 21, input$bsa, tinf = 3), p = list()),
      list(k = "weekly 20 mg/m²",
           e = dox_ev(20, round(cum / 20), 7, input$bsa), p = list()),
      list(k = "liposomal (PLD)",
           e = dox_ev(40, round(cum / 40), 28, input$bsa, lipo = TRUE), p = list()),
      list(k = "push + dexrazoxane",
           e = c(dox_ev(60, round(cum / 60), 21, input$bsa),
                 dex_ev(60, round(cum / 60), 21, input$bsa)), p = list()),
      list(k = "push + statin",
           e = dox_ev(60, round(cum / 60), 21, input$bsa),
           p = list(TGT_STA = 0.75, TON_STA = 0)),
      list(k = "push + ACEi/BB",
           e = dox_ev(60, round(cum / 60), 21, input$bsa),
           p = list(TGT_ACE = 0.8, TON_ACE = 0, TGT_BB = 0.8, TON_BB = 0)))
    bind_rows(lapply(defs, function(d) {
      m <- if (length(d$p)) param(base, d$p) else base
      o <- m %>% ev(d$e) %>%
        mrgsim(end = input$tend, delta = 1, nid = input$nid) %>% as_tibble()
      i <- m %>% ev(d$e) %>% mrgsim(end = input$tend, delta = 1) %>%
        as_tibble()
      o %>% group_by(ID) %>%
        summarise(CTRCD = any(LVEF < 50 & (first(LVEF) - LVEF) >= 10),
                  HF = any(LVEF < 40), nadir = min(LVEF), .groups = "drop") %>%
        summarise(scenario = d$k, `CTRCD %` = 100 * mean(CTRCD),
                  `LVEF<40 %` = 100 * mean(HF),
                  `mean nadir` = mean(nadir)) %>%
        mutate(`peak nuclear CHF` = max(i$CHF), `peak p53` = max(i$P53),
               `AUC retained` = max(i$AUCH))
    }))
  })

  output$cmp_tbl <- renderDT({
    datatable(cmp() %>% mutate(across(where(is.numeric), ~round(.x, 2))),
              options = list(dom = "t", pageLength = 10), rownames = FALSE)
  })

  output$cmp_plot <- renderPlot({
    d <- cmp()
    ggplot(d, aes(reorder(scenario, `CTRCD %`), `CTRCD %`)) +
      geom_col(fill = "#c0392b", alpha = .85) +
      geom_text(aes(label = sprintf("%.1f%%", `CTRCD %`)), hjust = -0.15) +
      coord_flip() + ylim(0, max(d$`CTRCD %`) * 1.2) +
      labs(x = NULL, y = "CTRCD (%)",
           title = sprintf("누적용량 %d mg/m² 고정 — 스케줄·제형·보호제만 다름",
                           input$dose * input$ncyc)) +
      theme_aic()
  })

  ## ---- tab 7 ------------------------------------------------------------
  output$tni_plot <- renderPlot({
    d <- sim()
    q <- d %>% group_by(time) %>%
      summarise(med = median(CTNI), hi = quantile(CTNI, .95), .groups = "drop")
    ggplot(q, aes(time)) +
      geom_ribbon(aes(ymin = 0, ymax = hi), alpha = .15) +
      geom_line(aes(y = med), linewidth = 1) +
      geom_hline(yintercept = 14, linetype = 3) +
      labs(x = "day", y = "hs-cTnI (ng/L)",
           title = "점선 = 14 ng/L (99th percentile)") + theme_aic()
  })

  output$bnp_plot <- renderPlot({
    d <- sim()
    q <- d %>% group_by(time) %>%
      summarise(med = median(NTBNP), hi = quantile(NTBNP, .95),
                .groups = "drop")
    ggplot(q, aes(time)) +
      geom_ribbon(aes(ymin = 0, ymax = hi), alpha = .15) +
      geom_line(aes(y = med), linewidth = 1) +
      labs(x = "day", y = "NT-proBNP (pg/mL)") + theme_aic()
  })

  output$lead_tbl <- renderTable({
    d <- sim()
    f <- function(cond_expr) {
      d %>% group_by(ID) %>% filter({{ cond_expr }}) %>%
        summarise(day = min(time), .groups = "drop")
    }
    tni <- f(CTNI > 14)
    gls <- d %>% group_by(ID) %>%
      filter((first(GLS) - GLS) / first(GLS) >= 0.15) %>%
      summarise(day = min(time), .groups = "drop")
    ef <- d %>% group_by(ID) %>%
      filter(LVEF < 50, (first(LVEF) - LVEF) >= 10) %>%
      summarise(day = min(time), .groups = "drop")
    nid <- n_distinct(d$ID)
    tibble(지표 = c("hs-cTnI > 14 ng/L", "GLS 상대 15% 감소", "CTRCD (LVEF)"),
           `도달 환자 비율` = sprintf("%.1f %%",
                                      100 * c(nrow(tni), nrow(gls),
                                              nrow(ef)) / nid),
           `중앙 도달일` = c(if (nrow(tni)) median(tni$day) else NA,
                             if (nrow(gls)) median(gls$day) else NA,
                             if (nrow(ef)) median(ef$day) else NA))
  })

  ## ---- tab 8 ------------------------------------------------------------
  win <- eventReactive(input$go_win, {
    starts <- c(60, 90, 120, 180, 270, 365, 540)
    bind_rows(lapply(starts, function(t0) {
      m <- param(mod, c(risk_par(), list(TGT_ACE = 0.8, TON_ACE = t0,
                                         TGT_BB = 0.8, TON_BB = t0,
                                         TGT_STA = 0, TGT_ARNI = 0,
                                         TGT_SGLT = 0)))
      o <- m %>% ev(regimen()) %>%
        mrgsim(end = t0 + 365, delta = 1, nid = input$nid) %>% as_tibble()
      o %>% group_by(ID) %>%
        summarise(had = any(LVEF < 50 & (first(LVEF) - LVEF) >= 10 &
                              time <= 365),
                  EF_base = first(LVEF),
                  EF_start = LVEF[which.min(abs(time - t0))],
                  FIB_start = FIBROSIS[which.min(abs(time - t0))],
                  EF_end = last(LVEF), .groups = "drop") %>%
        filter(had) %>%
        summarise(`시작일` = t0, n = n(),
                  `시작 시 LVEF` = mean(EF_start),
                  `12개월 후 LVEF` = mean(EF_end),
                  `회복 비율 %` = 100 * mean(EF_end >= 50 |
                                               (EF_base - EF_end) < 5),
                  `시작 시 FIB` = mean(FIB_start))
    }))
  })

  output$win_tbl <- renderDT({
    datatable(win() %>% mutate(across(where(is.numeric), ~round(.x, 2))),
              options = list(dom = "t"), rownames = FALSE)
  })

  output$win_plot <- renderPlot({
    d <- win()
    ggplot(d, aes(`시작일`)) +
      geom_line(aes(y = `회복 비율 %`), linewidth = 1.2, colour = "#c0392b") +
      geom_point(aes(y = `회복 비율 %`), size = 3, colour = "#c0392b") +
      geom_line(aes(y = 100 * `시작 시 FIB`), linewidth = 1,
                linetype = 2, colour = "#616161") +
      labs(x = "심부전 치료 시작일 (day)",
           y = "회복 비율 (%) / 100 × FIB",
           title = "회복은 시간이 아니라 섬유화를 따라간다 (점선 = FIB)") +
      theme_aic()
  })
}

shinyApp(ui, server)
