## =====================================================================
##  ed_shiny_app.R
##  Erectile Dysfunction QSP model — interactive dashboard
##  발기부전 QSP 모델 — 인터랙티브 대시보드
##
##  14 tabs:
##   1  환자 프로파일          patient profile and derived baseline
##   2  약동학                 PK of the four oral PDE5 inhibitors
##   3  발기 역학              intracavernosal pressure / rigidity time course
##   4  신호전달 캐스케이드    NO -> cGMP -> PKG -> Ca -> MLC -> relaxation
##   5  혈역학                 inflow, outflow, veno-occlusion, P-V loop
##   6  용량-반응              dose-response and the saturation of the gain
##   7  작용시간 창            time window after a single dose (the 36 h claim)
##   8  임상 엔드포인트        IIEF-EF, SEP2/SEP3, EHS, MCID
##   9  시나리오 비교          the 16 therapeutic scenarios side by side
##  10  구조 리모델링          the slow arm: SMI, collagen, smooth muscle
##  11  음경 재활              REACTT-style rehabilitation with washout
##  12  안전성                 MAP, nitrate interaction, PDE6/PDE11 occupancy
##  13  가상 집단              virtual population and responder fraction
##  14  바이오마커             ADMA, ROS, endocrine axis, Doppler proxies
##
##  Run:  shiny::runApp("ed_shiny_app.R")
##  Requires: shiny, mrgsolve, dplyr, tidyr, ggplot2, DT
## =====================================================================

library(shiny)
library(mrgsolve)
library(dplyr)
library(tidyr)
library(ggplot2)

source("ed_mrgsolve_model.R")     # provides ed_mod, ED_DRUGS, ED_ARCH, helpers

THEME <- theme_bw(base_size = 12) +
  theme(panel.grid.minor = element_blank(),
        strip.background = element_rect(fill = "#ECEFF1"),
        legend.position = "bottom")

PAL <- c("#1E88E5", "#D81B60", "#43A047", "#FB8C00", "#8E24AA", "#00897B",
         "#6D4C41", "#546E7A")

## ============================== UI ==================================
ui <- fluidPage(
  titlePanel("발기부전 (Erectile Dysfunction) QSP 모델 — 45 ODE / 4 PDE5 억제제"),
  tags$p(style = "color:#546E7A;margin-top:-8px",
         paste("발기는 포화형 증폭기에 놓인 문턱 판독이다:",
               "질산화물 펄스 × cGMP 이득 × 정맥폐쇄 천장 → 강직도 문턱.",
               "PDE5 억제제는 이득만 곱한다.")),
  sidebarLayout(
    sidebarPanel(
      width = 3,
      selectInput("arch", "환자 원형 (archetype)",
                  choices = names(ED_ARCH), selected = "vasculo"),
      selectInput("drug", "경구 PDE5 억제제",
                  choices = names(ED_DRUGS), selected = "sildenafil"),
      sliderInput("dose", "용량 (mg)", 0, 200, 100, step = 5),
      sliderInput("lead", "복용 → 시도 간격 (h)", 0.25, 48, 1, step = 0.25),
      hr(),
      h5("동반질환 / 위험인자"),
      sliderInput("agey", "연령 (년)", 20, 85, 60),
      sliderInput("hba1c", "HbA1c (%)", 4.5, 13, 5.6, step = 0.1),
      sliderInput("bmi", "BMI (kg/m²)", 18, 45, 29, step = 0.5),
      sliderInput("ldl", "LDL-C (mg/dL)", 50, 250, 145, step = 5),
      sliderInput("smoke", "흡연 (0 = 비흡연, 1 = 현재)", 0, 1, 0.5, step = 0.1),
      sliderInput("raas", "RAAS 활성 (-)", 1, 1.8, 1.25, step = 0.05),
      hr(),
      h5("신경 · 내분비 · 심리"),
      sliderInput("nrv", "동체신경 온전성 NRV (0–1)", 0.05, 1, 0.43, step = 0.01),
      sliderInput("nrvmax", "신경보존 회복 천장 NRVMAX", 0.10, 1, 1.0, step = 0.05),
      sliderInput("kstl", "Leydig 기능 KSTL", 40, 260, 229.6, step = 5),
      sliderInput("stress", "심리적 스트레스 STRESS", 0, 1.5, 0, step = 0.05),
      hr(),
      h5("보조 · 병용"),
      sliderInput("pge", "알프로스타딜 ICI (µg)", 0, 40, 0, step = 2),
      sliderInput("pap", "파파베린 (mg)", 0, 40, 0, step = 2),
      sliderInput("phen", "펜톨라민 (mg)", 0, 2, 0, step = 0.1),
      sliderInput("gtn", "니트로글리세린 (ng/mL 상당)", 0, 3, 0, step = 0.05),
      sliderInput("cact", "sGC 활성제 (cinaciguat 상당)", 0, 3, 0, step = 0.1),
      sliderInput("cstim", "sGC 자극제 (riociguat 상당)", 0, 3, 0, step = 0.1),
      sliderInput("rocki", "ROCK 억제제", 0, 3, 0, step = 0.1),
      checkboxInput("exer", "운동 프로그램", FALSE),
      checkboxInput("statin", "스타틴", FALSE),
      actionButton("go", "시뮬레이션 실행", class = "btn-primary")
    ),
    mainPanel(
      width = 9,
      tabsetPanel(
        id = "tabs", type = "tabs",
        tabPanel("1 환자 프로파일",
                 h4("유도된 기저 상태 (모든 값은 방정식에서 계산됨)"),
                 tableOutput("t_base"),
                 plotOutput("p_base", height = "330px"),
                 helpText(paste("SMI = 평활근/(평활근+콜라겐), 정맥폐쇄 천장을",
                                "결정한다. STEN은 축적된 산화 부하에서 계산된",
                                "협착 분율이다."))),
        tabPanel("2 약동학",
                 plotOutput("p_pk", height = "330px"),
                 plotOutput("p_pk2", height = "300px"),
                 tableOutput("t_pk"),
                 helpText("Cu/IC50가 PDE5 잔존활성 1/(1+Cu/IC50)을 결정한다.")),
        tabPanel("3 발기 역학",
                 plotOutput("p_icp", height = "380px"),
                 fluidRow(column(6, tableOutput("t_att")),
                          column(6, plotOutput("p_ehs", height = "260px")))),
        tabPanel("4 신호전달",
                 plotOutput("p_sig", height = "460px"),
                 helpText(paste("NO는 신경이 발화하는 동안에만 생성된다.",
                                "PDE5 억제제는 NO를 만들지 않고 cGMP 이득만",
                                "곱한다."))),
        tabPanel("5 혈역학",
                 fluidRow(column(7, plotOutput("p_flow", height = "340px")),
                          column(5, plotOutput("p_pv", height = "340px"))),
                 plotOutput("p_vocc", height = "280px")),
        tabPanel("6 용량-반응",
                 plotOutput("p_dr", height = "360px"),
                 tableOutput("t_dr"),
                 helpText(paste("증폭기가 포화하므로 100 mg 이상에서 거의",
                                "변화가 없다 — 라벨 최대용량의 구조적 이유."))),
        tabPanel("7 작용시간 창",
                 plotOutput("p_win", height = "380px"),
                 tableOutput("t_win")),
        tabPanel("8 임상 엔드포인트",
                 fluidRow(column(6, plotOutput("p_end", height = "330px")),
                          column(6, plotOutput("p_mcid", height = "330px"))),
                 tableOutput("t_end")),
        tabPanel("9 시나리오 비교",
                 plotOutput("p_scn", height = "420px"),
                 tableOutput("t_scn")),
        tabPanel("10 구조 리모델링",
                 sliderInput("months", "관찰 기간 (개월)", 1, 36, 12),
                 plotOutput("p_struct", height = "420px"),
                 tableOutput("t_struct")),
        tabPanel("11 음경 재활",
                 helpText(paste("REACTT 설계: 9개월 tadalafil 5 mg 매일 /",
                                "20 mg 온디맨드 / 위약 → 6주 약물중단 세척기 후",
                                "무보조 발기 평가.")),
                 plotOutput("p_rehab", height = "400px"),
                 tableOutput("t_rehab")),
        tabPanel("12 안전성",
                 plotOutput("p_map", height = "330px"),
                 fluidRow(column(6, plotOutput("p_iso", height = "300px")),
                          column(6, tableOutput("t_safe")))),
        tabPanel("13 가상 집단",
                 sliderInput("npop", "가상 환자 수", 11, 81, 31, step = 10),
                 sliderInput("sigma", "NO 생성능 로그정규 SD", 0.2, 1.0, 0.55,
                             step = 0.05),
                 plotOutput("p_pop", height = "380px"),
                 tableOutput("t_pop")),
        tabPanel("14 바이오마커",
                 plotOutput("p_bio", height = "420px"),
                 tableOutput("t_bio"))
      )
    )
  )
)

## ============================ server ================================
server <- function(input, output, session) {

  ## ---- assemble the patient ----------------------------------------
  patient <- eventReactive(input$go, {
    A <- ED_ARCH[[input$arch]]
    par <- modifyList(A$par, list(
      AGEY = input$agey, HBA1C = input$hba1c, BMI = input$bmi,
      LDL = input$ldl, SMOKE = input$smoke, RAAS = input$raas,
      KSTL = input$kstl, STRESS = input$stress, NRVMAX = input$nrvmax,
      EXER = as.numeric(input$exer), STATIN = as.numeric(input$statin),
      CACT = input$cact, CSTIM = input$cstim, CROCKI = input$rocki))
    arch2 <- A; arch2$par <- par; arch2$nrv <- input$nrv
    ED_ARCH[["_ui"]] <<- arch2
    ed_baseline(ed_mod, "_ui", input$drug)
  }, ignoreNULL = FALSE)

  attempt <- reactive({
    ed_attempt(patient(), dose = input$dose, lead = input$lead,
               pge = input$pge, pap = input$pap, phen = input$phen,
               gtn = input$gtn, window = input$lead + 4)
  })

  ## ---- 1 patient profile -------------------------------------------
  output$t_base <- renderTable({
    y <- patient()$ini
    d <- as.data.frame(mrgsim(patient()$mod, end = 24, delta = 0.02,
                              hmax = 0.02))
    data.frame(
      `지표` = c("ROS (배)", "eNOS 결합 분율", "ADMA (µmol/L)", "산화 sGC 분율",
               "협착 STEN", "평활근 SM", "콜라겐 COL", "SMI",
               "정맥폐쇄 천장", "동체 pO₂ (mmHg)", "nNOS", "PDE5A 발현",
               "총 테스토스테론 (ng/dL)", "LH (IU/L)", "수행불안 PA",
               "야간발기 최고 ICP (mmHg)"),
      `값` = sprintf("%.3g", c(y$ROS, y$ECPL, y$ADMA, y$SGCOX, y$STEN, y$SM,
                             y$COL, y$SM / (y$SM + y$COL) * 2, max(d$VOCMAXo),
                             y$PO2, y$NNOS, y$PDE5E, y$TT, y$LH,
                             tail(d$PA, 1), max(d$ICPo))))
  })
  output$p_base <- renderPlot({
    y <- patient()$ini
    ref <- c(ROS = 1, ECPL = 0.714, ADMA = 0.45, SM = 1, COL = 1,
             NNOS = 1, PDE5E = 1, PO2 = 95)
    v <- c(ROS = y$ROS, ECPL = y$ECPL, ADMA = y$ADMA, SM = y$SM, COL = y$COL,
           NNOS = y$NNOS, PDE5E = y$PDE5E, PO2 = y$PO2)
    ggplot(data.frame(k = names(v), rel = as.numeric(v / ref[names(v)])),
           aes(reorder(k, rel), rel, fill = rel > 1)) +
      geom_col(width = .65) + geom_hline(yintercept = 1, linetype = 2) +
      coord_flip() + scale_fill_manual(values = c("#1E88E5", "#D81B60")) +
      labs(x = NULL, y = "건강 기준값 대비 비율", fill = "> 정상") + THEME
  })

  ## ---- 2 PK ---------------------------------------------------------
  pkall <- reactive({
    do.call(rbind, lapply(names(ED_DRUGS), function(dg) {
      ds <- c(sildenafil = 100, tadalafil = 20, vardenafil = 20,
              avanafil = 200)[[dg]]
      m <- ed_drug(ed_mod, dg) %>% param(NPTON = 0, ATTT = -1)
      d <- as.data.frame(mrgsim(m, events = ev(amt = ds, cmt = "AGUT"),
                                end = 60, delta = 0.05))
      data.frame(drug = dg, dose = ds, time = d$time, Cp = d$CP1 * 1000,
                 Cu = d$CUo, ratio = d$RATIO, res = d$P5RESo,
                 p6 = d$P6INH, p11 = d$P11INH)
    }))
  })
  output$p_pk <- renderPlot({
    ggplot(pkall(), aes(time, Cu, colour = drug)) + geom_line(linewidth = 1) +
      scale_y_log10() + scale_colour_manual(values = PAL) +
      labs(x = "시간 (h)", y = "비결합 혈장 농도 Cu (nM)",
           title = "라벨 최대용량에서의 비결합 노출") + THEME
  })
  output$p_pk2 <- renderPlot({
    ggplot(pkall(), aes(time, res, colour = drug)) + geom_line(linewidth = 1) +
      scale_colour_manual(values = PAL) + ylim(0, 1) +
      labs(x = "시간 (h)", y = "PDE5 잔존활성 1/(1+Cu/IC50)",
           title = "네 약제 모두 포화 구간에서 출발한다") + THEME
  })
  output$t_pk <- renderTable({
    pkall() %>% group_by(drug, dose) %>%
      summarise(`Cmax (ng/mL)` = max(Cp), `Cu,max (nM)` = max(Cu),
                `Cu/IC50` = max(ratio), `최소 잔존활성` = min(res),
                `PDE6 점유 (%)` = max(p6), `PDE11 점유 (%)` = max(p11),
                .groups = "drop")
  }, digits = 3)

  ## ---- 3 erection dynamics -----------------------------------------
  output$p_icp <- renderPlot({
    d <- attempt()$sim
    d %>% select(time, ICPo, RIGo, VSIN) %>%
      pivot_longer(-time) %>%
      mutate(name = recode(name, ICPo = "동체내압 ICP (mmHg)",
                           RIGo = "축방향 강직도 (%)",
                           VSIN = "동체 혈액량 (mL)")) %>%
      ggplot(aes(time, value)) + geom_line(colour = "#1E88E5", linewidth = 1) +
      geom_hline(data = data.frame(name = "축방향 강직도 (%)", y = 60),
                 aes(yintercept = y), linetype = 2, colour = "#D81B60") +
      facet_wrap(~name, scales = "free_y", ncol = 1) +
      labs(x = "시간 (h)", y = NULL) + THEME
  })
  output$t_att <- renderTable({
    r <- attempt()
    data.frame(`지표` = c("최고 ICP (mmHg)", "최고 강직도 (%)",
                        "문턱 초과 시간 (분)", "최고 cGMP", "최고 NO (nM)",
                        "최고 이완 R", "Cu (nM)", "Cu/IC50",
                        "SEP2 확률", "SEP3 확률", "IIEF-EF",
                        "최저 MAP (mmHg)"),
               `값` = sprintf("%.3g", c(r$ICPpk, r$RIGpk, r$TAEmin, r$CGMPpk,
                                      r$NOpk, r$Rpk, r$CUpk, r$RATIO,
                                      r$PS2, r$PS3, r$IIEF, r$MAPmin)))
  })
  output$p_ehs <- renderPlot({
    d <- attempt()$sim
    ggplot(d, aes(time, EHS)) + geom_step(colour = "#43A047", linewidth = 1) +
      ylim(1, 4) + labs(x = "시간 (h)", y = "발기 경도 점수 EHS") + THEME
  })

  ## ---- 4 signalling -------------------------------------------------
  output$p_sig <- renderPlot({
    attempt()$sim %>%
      select(time, NO, CGMP, CAMP, CAI, MLCP, Ro) %>%
      pivot_longer(-time) %>%
      mutate(name = factor(name, c("NO", "CGMP", "CAMP", "CAI", "MLCP", "Ro"),
                           c("NO (nM)", "cGMP (배)", "cAMP (배)",
                             "세포내 Ca²⁺ (배)", "인산화 MLC20 분율",
                             "이완 R"))) %>%
      ggplot(aes(time, value)) + geom_line(colour = "#8E24AA", linewidth = 1) +
      facet_wrap(~name, scales = "free_y", ncol = 2) +
      labs(x = "시간 (h)", y = NULL) + THEME
  })

  ## ---- 5 haemodynamics ---------------------------------------------
  output$p_flow <- renderPlot({
    attempt()$sim %>% select(time, QINo, QOUTo) %>% pivot_longer(-time) %>%
      mutate(name = recode(name, QINo = "동맥 유입 Q_in",
                           QOUTo = "정맥 유출 Q_out")) %>%
      ggplot(aes(time, value, colour = name)) + geom_line(linewidth = 1) +
      scale_colour_manual(values = PAL) +
      labs(x = "시간 (h)", y = "유량 (mL/h)", colour = NULL) + THEME
  })
  output$p_pv <- renderPlot({
    ggplot(attempt()$sim, aes(VSIN, ICPo, colour = time)) +
      geom_path(linewidth = 1) + scale_colour_viridis_c(option = "C") +
      labs(x = "동체 혈액량 (mL)", y = "ICP (mmHg)", colour = "시간 (h)",
           title = "압력-용적 궤적 (백막의 지수적 경화)") + THEME
  })
  output$p_vocc <- renderPlot({
    attempt()$sim %>% select(time, VOCCo, GVENo) %>% pivot_longer(-time) %>%
      mutate(name = recode(name, VOCCo = "정맥폐쇄 관여도 VOCC",
                           GVENo = "정맥 누출 전도도 G_ven")) %>%
      ggplot(aes(time, value)) + geom_line(colour = "#00897B", linewidth = 1) +
      facet_wrap(~name, scales = "free_y") + labs(x = "시간 (h)", y = NULL) +
      THEME
  })

  ## ---- 6 dose-response ---------------------------------------------
  drtab <- reactive({
    b <- patient()
    do.call(rbind, lapply(c(0, 5, 10, 25, 50, 100, 150, 200), function(ds) {
      r <- ed_attempt(b, dose = ds, lead = input$lead)
      data.frame(dose = ds, Cu = r$CUpk, ratio = r$RATIO, cGMP = r$CGMPpk,
                 R = r$Rpk, ICP = r$ICPpk, RIG = r$RIGpk, IIEF = r$IIEF)
    }))
  })
  output$p_dr <- renderPlot({
    drtab() %>% select(dose, cGMP, R, ICP, IIEF) %>% pivot_longer(-dose) %>%
      ggplot(aes(dose, value)) +
      geom_line(colour = "#FB8C00", linewidth = 1) + geom_point() +
      facet_wrap(~name, scales = "free_y") +
      labs(x = sprintf("%s 용량 (mg)", input$drug), y = NULL) + THEME
  })
  output$t_dr <- renderTable(drtab(), digits = 3)

  ## ---- 7 duration window -------------------------------------------
  wintab <- reactive({
    b <- patient()
    lags <- c(0.5, 1, 2, 4, 6, 8, 12, 18, 24, 36, 48, 72)
    do.call(rbind, lapply(lags, function(L) {
      r <- ed_attempt(b, dose = input$dose, lead = L, window = L + 3)
      data.frame(lag = L, Cu = r$CUpk, ICP = r$ICPpk, RIG = r$RIGpk,
                 SEP3 = r$PS3, IIEF = r$IIEF)
    }))
  })
  output$p_win <- renderPlot({
    ggplot(wintab(), aes(lag, RIG)) +
      geom_line(colour = "#1E88E5", linewidth = 1) + geom_point() +
      geom_hline(yintercept = 60, linetype = 2, colour = "#D81B60") +
      scale_x_log10() +
      labs(x = "복용 후 경과 시간 (h, 로그)", y = "최고 강직도 (%)",
           title = "단일 용량 후 유효 창") + THEME
  })
  output$t_win <- renderTable(wintab(), digits = 3)

  ## ---- 8 endpoints --------------------------------------------------
  output$p_end <- renderPlot({
    b <- patient()
    r0 <- ed_attempt(b); r1 <- attempt()
    df <- data.frame(
      arm = rep(c("무치료", "치료"), each = 4),
      k = rep(c("SEP2", "SEP3", "IIEF-EF/30", "강직도/100"), 2),
      v = c(r0$PS2, r0$PS3, r0$IIEF / 30, r0$RIGpk / 100,
            r1$PS2, r1$PS3, r1$IIEF / 30, r1$RIGpk / 100))
    ggplot(df, aes(k, v, fill = arm)) +
      geom_col(position = "dodge", width = .7) +
      scale_fill_manual(values = PAL) + ylim(0, 1) +
      labs(x = NULL, y = "비율 / 정규화 점수", fill = NULL) + THEME
  })
  output$p_mcid <- renderPlot({
    b <- patient(); r0 <- ed_attempt(b); r1 <- attempt()
    d <- r1$IIEF - r0$IIEF
    sev <- if (r0$IIEF < 11) "중증" else if (r0$IIEF < 17) "중등도" else "경증"
    mcid <- c("중증" = 7, "중등도" = 5, "경증" = 2)[[sev]]
    ggplot(data.frame(x = c("ΔIIEF-EF", "MCID"), y = c(d, mcid)),
           aes(x, y, fill = x)) + geom_col(width = .55) +
      scale_fill_manual(values = PAL) +
      labs(x = NULL, y = "IIEF-EF 점수",
           title = sprintf("기저 중증도: %s (MCID %d점)", sev, mcid)) + THEME
  })
  output$t_end <- renderTable({
    b <- patient(); r0 <- ed_attempt(b); r1 <- attempt()
    data.frame(`구분` = c("무치료", "치료", "차이"),
               `IIEF-EF` = c(r0$IIEF, r1$IIEF, r1$IIEF - r0$IIEF),
               SEP2 = c(r0$PS2, r1$PS2, r1$PS2 - r0$PS2),
               SEP3 = c(r0$PS3, r1$PS3, r1$PS3 - r0$PS3),
               `최고 ICP` = c(r0$ICPpk, r1$ICPpk, r1$ICPpk - r0$ICPpk),
               check.names = FALSE)
  }, digits = 3)

  ## ---- 9 scenarios --------------------------------------------------
  scn <- reactive({
    b <- patient()
    rows <- list(
      c("무치료", 0, 0, 0, 0, 0),
      c("PDE5i 저용량", input$dose / 2, 0, 0, 0, 0),
      c("PDE5i 표준용량", input$dose, 0, 0, 0, 0),
      c("알프로스타딜 10 µg ICI", 0, 10, 0, 0, 0),
      c("알프로스타딜 20 µg ICI", 0, 20, 0, 0, 0),
      c("삼중혼합 (trimix)", 0, 10, 15, 0.5, 0),
      c("PDE5i + 니트로글리세린", input$dose, 0, 0, 0, 0.76))
    do.call(rbind, lapply(rows, function(z) {
      r <- ed_attempt(b, dose = as.numeric(z[2]), lead = input$lead,
                      pge = as.numeric(z[3]), pap = as.numeric(z[4]),
                      phen = as.numeric(z[5]), gtn = as.numeric(z[6]))
      data.frame(`시나리오` = z[1], ICP = r$ICPpk, `강직도` = r$RIGpk,
                 SEP3 = r$PS3, `IIEF-EF` = r$IIEF, `최저 MAP` = r$MAPmin,
                 check.names = FALSE)
    }))
  })
  output$p_scn <- renderPlot({
    scn() %>% mutate(`시나리오` = factor(`시나리오`, `시나리오`)) %>%
      ggplot(aes(`시나리오`, `IIEF-EF`, fill = ICP)) +
      geom_col(width = .7) + coord_flip() +
      scale_fill_viridis_c(option = "D") +
      labs(x = NULL, y = "IIEF-EF", fill = "최고 ICP") + THEME
  })
  output$t_scn <- renderTable(scn(), digits = 3)

  ## ---- 10 structural arm -------------------------------------------
  strun <- reactive({
    m <- patient()$mod %>% param(ATTEVERY = 120)
    as.data.frame(mrgsim(m, end = 24 * 30 * input$months, delta = 12,
                         hmax = 0.05))
  })
  output$p_struct <- renderPlot({
    strun() %>% mutate(day = time / 24) %>%
      select(day, SM, COL, SMIo, TGFB, PO2, LEN, NRV, ERFR) %>%
      pivot_longer(-day) %>%
      ggplot(aes(day, value)) + geom_line(colour = "#6D4C41", linewidth = 1) +
      facet_wrap(~name, scales = "free_y") +
      labs(x = "일", y = NULL, title = "느린 팔: 주–개월 시간상수") + THEME
  })
  output$t_struct <- renderTable({
    d <- strun()
    data.frame(`지표` = c("SMI 시작", "SMI 종료", "콜라겐 종료", "평활근 종료",
                        "음경 길이 (cm)", "NRV 종료", "동체 pO₂"),
               `값` = sprintf("%.3g", c(d$SMIo[1], tail(d$SMIo, 1),
                                      tail(d$COL, 1), tail(d$SM, 1),
                                      tail(d$LEN, 1), tail(d$NRV, 1),
                                      tail(d$PO2, 1))))
  })

  ## ---- 11 rehabilitation -------------------------------------------
  rehab <- reactive({
    tab <- ed_s9_rehab(ed_mod)
    tab
  })
  output$p_rehab <- renderPlot({
    rehab() %>% select(arm, SMI_onTx, SMI_wash, ICP_unassisted,
                       IIEF_unassisted) %>%
      pivot_longer(-arm) %>%
      ggplot(aes(arm, value, fill = arm)) + geom_col(width = .6) +
      facet_wrap(~name, scales = "free_y") +
      scale_fill_manual(values = PAL) +
      labs(x = NULL, y = NULL, fill = NULL,
           title = "9개월 치료 후 6주 약물중단 세척기") + THEME
  })
  output$t_rehab <- renderTable(rehab(), digits = 3)

  ## ---- 12 safety ----------------------------------------------------
  output$p_map <- renderPlot({
    b <- patient()
    arms <- list(`무치료` = list(0, 0), `니트로만` = list(0, 0.76),
                 `PDE5i만` = list(input$dose, 0),
                 `PDE5i + 니트로` = list(input$dose, 0.76))
    d <- do.call(rbind, lapply(names(arms), function(a) {
      r <- ed_attempt(b, dose = arms[[a]][[1]], lead = input$lead,
                      gtn = arms[[a]][[2]], window = input$lead + 3)
      data.frame(arm = a, time = r$sim$time, MAP = r$sim$MAP)
    }))
    ggplot(d, aes(time, MAP, colour = arm)) + geom_line(linewidth = 1) +
      scale_colour_manual(values = PAL) +
      labs(x = "시간 (h)", y = "평균 동맥압 (mmHg)", colour = NULL,
           title = "같은 두 인자 곱을 전신 혈관에서 읽은 결과") + THEME
  })
  output$p_iso <- renderPlot({
    d <- attempt()$sim
    d %>% select(time, P6INH, P11INH) %>% pivot_longer(-time) %>%
      mutate(name = recode(name, P6INH = "PDE6 점유 (%) — 시각",
                           P11INH = "PDE11 점유 (%) — 근육통")) %>%
      ggplot(aes(time, value, colour = name)) + geom_line(linewidth = 1) +
      scale_colour_manual(values = PAL) +
      labs(x = "시간 (h)", y = "억제율 (%)", colour = NULL) + THEME
  })
  output$t_safe <- renderTable({
    r <- attempt()
    data.frame(`지표` = c("최저 MAP (mmHg)", "MAP 변화 (mmHg)",
                        "최대 PDE6 점유 (%)", "최대 PDE11 점유 (%)",
                        "최대 전신 cGMP (배)"),
               `값` = sprintf("%.3g", c(r$MAPmin, min(r$sim$DMAP),
                                      max(r$sim$P6INH), max(r$sim$P11INH),
                                      max(r$sim$SCG))))
  })

  ## ---- 13 virtual population ---------------------------------------
  poptab <- reactive({
    q <- (seq_len(input$npop) - 0.5) / input$npop
    fn <- exp(input$sigma * qnorm(q))
    do.call(rbind, lapply(fn, function(f) {
      b <- ed_baseline(ed_mod, "_ui", input$drug, fnoi = f)
      r0 <- ed_attempt(b); r1 <- ed_attempt(b, dose = input$dose,
                                            lead = input$lead)
      data.frame(FNOI = f, IIEF0 = r0$IIEF, IIEF1 = r1$IIEF,
                 ICP0 = r0$ICPpk, ICP1 = r1$ICPpk,
                 SEP30 = r0$PS3, SEP31 = r1$PS3)
    }))
  })
  output$p_pop <- renderPlot({
    poptab() %>% select(FNOI, IIEF0, IIEF1) %>% pivot_longer(-FNOI) %>%
      mutate(name = recode(name, IIEF0 = "무치료", IIEF1 = "치료")) %>%
      ggplot(aes(FNOI, value, colour = name)) +
      geom_line(linewidth = 1) + geom_point(size = 1) + scale_x_log10() +
      scale_colour_manual(values = PAL) +
      labs(x = "개인별 NO 생성능 배수 FNOI (로그)", y = "IIEF-EF",
           colour = NULL,
           title = "문턱 판독 + 개체간 변이 = 집단 평균") + THEME
  })
  output$t_pop <- renderTable({
    d <- poptab()
    data.frame(`구분` = c("무치료", "치료"),
               `평균 IIEF-EF` = c(mean(d$IIEF0), mean(d$IIEF1)),
               `평균 SEP3` = c(mean(d$SEP30), mean(d$SEP31)),
               `반응자 비율 (SEP3≥0.5)` = c(mean(d$SEP30 >= .5),
                                            mean(d$SEP31 >= .5)),
               check.names = FALSE)
  }, digits = 3)

  ## ---- 14 biomarkers ------------------------------------------------
  output$p_bio <- renderPlot({
    strun() %>% mutate(day = time / 24) %>%
      select(day, ROS, ADMA, ECPL, OXD, SGCOX, TT, LH, HCT, PSV) %>%
      pivot_longer(-day) %>%
      ggplot(aes(day, value)) + geom_line(colour = "#546E7A", linewidth = 1) +
      facet_wrap(~name, scales = "free_y") + labs(x = "일", y = NULL) + THEME
  })
  output$t_bio <- renderTable({
    d <- strun()
    data.frame(`바이오마커` = c("ROS (배)", "ADMA (µmol/L)", "eNOS 결합",
                              "산화 손상", "산화 sGC", "테스토스테론",
                              "LH", "헤마토크릿 (%)",
                              "동체동맥 최고속도 대용치 (cm/s)"),
               `값` = sprintf("%.3g", c(tail(d$ROS, 1), tail(d$ADMA, 1),
                                      tail(d$ECPL, 1), tail(d$OXD, 1),
                                      tail(d$SGCOX, 1), tail(d$TT, 1),
                                      tail(d$LH, 1), tail(d$HCT, 1),
                                      max(d$PSV))))
  })
}

shinyApp(ui, server)
