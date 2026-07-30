## =============================================================================
##  dlb_shiny_app.R
##  Dementia with Lewy Bodies — interactive QSP dashboard
##
##  Run with:
##      setwd("dementia-with-lewy-bodies")
##      shiny::runApp("dlb_shiny_app.R")
##
##  Requires: shiny, ggplot2, dplyr, tidyr, mrgsolve  (and dlb_mrgsolve_model.R
##  in the same directory).
##
##  ---------------------------------------------------------------------------
##  WHAT THIS APP IS FOR
##  ---------------------------------------------------------------------------
##  Not "here are ten line plots of a simulation".  Each tab is built around one
##  claim the model makes, and is laid out so that the claim can be CHECKED --
##  usually by putting the thing that is supposed to explain the effect directly
##  next to the effect.
##
##    Tab 3  puts the three transducers side by side, with the receptor-density
##           equation for each, so the sign difference is visible rather than
##           asserted.
##    Tab 4  plots the attention cubic's potential landscape next to the
##           measured fluctuation score, so "fluctuation is a variance" is a
##           picture and not a sentence.
##    Tab 5  decomposes the hallucination drive into its three factors and shows
##           what happens to the product when you knock out one at a time.
##    Tab 6  overlays D2 occupancy, postsynaptic reserve, and the resulting
##           deficit for whichever phenotype you have selected, so neuroleptic
##           sensitivity can be traced from cause to consequence in one screen.
##    Tab 7  lets you switch ambroxol on at any day and watch the GBA1 loop.
##
##  Licence: see repository LICENSE.  EDUCATIONAL / RESEARCH USE ONLY.
## =============================================================================

library(shiny)
library(mrgsolve)
library(ggplot2)
library(dplyr)
library(tidyr)

source("dlb_mrgsolve_model.R")

## ----------------------------------------------------------------- theming --
qsp_theme <- theme_minimal(base_size = 12) +
  theme(panel.grid.minor = element_blank(),
        strip.text  = element_text(face = "bold"),
        plot.title  = element_text(face = "bold", size = 13),
        plot.subtitle = element_text(colour = "grey35", size = 10),
        legend.position = "bottom",
        legend.title = element_blank())

PAL <- c(DLB = "#1f4a7a", `PD-dementia` = "#b8860b", AD = "#7a2f4a",
         treated = "#1f7a4a", untreated = "#8a8a99")

yr <- function(d) d$time / 365.25

## Every plot in the app goes through this, so that the x axis, the theme and
## the "years since diagnosis" convention are identical everywhere.
trace_plot <- function(d, vars, labels = vars, title = "", subtitle = "",
                       ylab = "", hlines = NULL) {
  long <- d %>%
    select(time, scenario, all_of(vars)) %>%
    pivot_longer(all_of(vars), names_to = "var", values_to = "y") %>%
    mutate(var = factor(var, levels = vars, labels = labels),
           yrs = time / 365.25)
  p <- ggplot(long, aes(yrs, y, colour = scenario, linetype = scenario)) +
    geom_line(linewidth = 0.75) +
    facet_wrap(~var, scales = "free_y") +
    labs(title = title, subtitle = subtitle,
         x = "진단 후 경과 연수 (years since diagnosis)", y = ylab) +
    qsp_theme
  if (!is.null(hlines))
    p <- p + geom_hline(yintercept = hlines, linetype = "dotted",
                        colour = "grey50")
  p
}

## =============================================================================
##  UI
## =============================================================================
ui <- fluidPage(
  titlePanel("루이소체 치매 QSP 대시보드 · Dementia with Lewy Bodies QSP Dashboard"),
  tags$p(style = "color:#555;margin-top:-8px;",
         "66-compartment mrgsolve model · three-transmitter presynaptic/postsynaptic dissociation · ",
         tags$b("교육·연구용 (educational / research use only)")),

  sidebarLayout(
    sidebarPanel(
      width = 3,

      h4("① 환자 (Patient)"),
      selectInput("pheno", "표현형 (phenotype)",
                  choices = c("DLB" = 0, "PD-dementia" = 1, "AD" = 2),
                  selected = 0),
      helpText(tags$small(
        "표현형은 라벨이 아니라 ", tags$b("초기 병리 분포"),
        "입니다. 뇌간/변연계/신피질 α-시누클레인과 tau·아밀로이드의 시작값만 바뀌고, ",
        "66개 방정식은 완전히 동일합니다.")),
      sliderInput("gbaf", "GBA1 기능 대립유전자 용량 (1 = WT)",
                  min = 0.30, max = 1.00, value = 1.00, step = 0.05),
      sliderInput("apoe4", "APOE ε4 대립유전자 수", min = 0, max = 2, value = 0, step = 1),
      sliderInput("snca", "SNCA 유전자 용량", min = 0.8, max = 1.6, value = 1.0, step = 0.05),
      sliderInput("antich", "항콜린제 부담 (anticholinergic burden)",
                  min = 0, max = 3, value = 0, step = 1),

      hr(), h4("② 치료 (Treatment)"),
      selectInput("chei", "콜린에스터라제 억제제",
                  choices = c("없음 (none)" = "none",
                              "리바스티그민 6 mg BID 경구" = "riv_oral",
                              "리바스티그민 9.5 mg/24h 패치" = "riv_patch",
                              "도네페질 10 mg qd" = "don")),
      numericInput("chei_start", "ChEI 시작일 (day)", value = 180, min = 0, max = 2000, step = 30),
      checkboxInput("mem",  "메만틴 20 mg qd", FALSE),
      checkboxInput("pim",  "피마반세린 34 mg qd (day 365~)", FALSE),
      selectInput("apd", "항정신병약 (day 730~, 180일)",
                  choices = c("없음 (none)" = "none",
                              "리스페리돈 1 mg qd" = "risp",
                              "쿠에티아핀 50 mg qd" = "quet")),
      checkboxInput("ld",   "레보도파/카비도파 150 mg TID (day 365~)", FALSE),
      checkboxInput("zon",  "조노사마이드 25 mg qd (day 365~)", FALSE),
      checkboxInput("amb",  "앰브록솔 1.26 g/일", FALSE),
      numericInput("amb_start", "앰브록솔 시작일 (day)", value = 0, min = 0, max = 2500, step = 30),
      checkboxInput("mab",  "항 α-시누클레인 항체 4500 mg IV q4w", FALSE),
      checkboxInput("mel",  "멜라토닌/클로나제팜 (RBD)", FALSE),
      checkboxInput("drox", "드록시도파/미도드린 (기립성 저혈압)", FALSE),
      checkboxInput("mod",  "모다피닐 (주간 졸림)", FALSE),

      hr(), h4("③ 시뮬레이션 (Simulation)"),
      sliderInput("years", "시뮬레이션 기간 (년)", min = 2, max = 10, value = 8, step = 1),
      checkboxInput("show_ref", "무치료 대조군 함께 표시", TRUE),
      actionButton("go", "실행 (Run)", class = "btn-primary btn-block"),
      hr(),
      tags$small(tags$b("주의."), " 본 모델은 공개 문헌으로 구성된 반정량적 교육용 모델이며, ",
                 "임상 의사결정·처방·규제 제출에 사용해서는 안 됩니다.")
    ),

    mainPanel(
      width = 9,
      tabsetPanel(
        id = "tabs",

        ## ------------------------------------------------------- 1. profile
        tabPanel("1 · 환자 프로파일",
                 br(), h4("진단 시점의 환자 상태 (state at diagnosis)"),
                 tableOutput("tbl_profile"),
                 h4("핵심·지표 특징 (core & indicative features)"),
                 plotOutput("p_profile", height = 420),
                 helpText("MIBG H/M < 2.0, DaTSCAN 결합비 감소, EEG 우세주파수 5.6–7.9 Hz(pre-alpha)는 ",
                          "4차 합의기준의 지표 바이오마커입니다. 세 값 모두 모델이 계산한 것이며, ",
                          "따로 지정된 상수가 아닙니다.")),

        ## ------------------------------------------------------------ 2. PK
        tabPanel("2 · 약동학 (PK)",
                 br(), plotOutput("p_pk", height = 480),
                 h4("리바스티그민: 혈장 반감기와 약력학 반감기의 분리"),
                 plotOutput("p_carb", height = 320),
                 helpText("리바스티그민의 혈장 반감기는 약 1.5시간이지만, 카바밀화된 효소의 회복은 ",
                          "효소 재합성 속도(약 9시간)로 결정됩니다. 모델은 카바밀화 효소를 ",
                          "별도의 상태변수로 두며, 그래서 12시간 간격 투여가 성립합니다. ",
                          "패치는 같은 AChE 억제율을 훨씬 낮은 Cmax로 달성하고, ",
                          "위장관 부작용 지표는 Cmax로 구동되므로 패치에서 낮게 나옵니다.")),

        ## ------------------------------------------- 3. the three transducers
        tabPanel("3 · 세 개의 전달자 ★",
                 br(),
                 div(style = "background:#fffbe6;border-left:4px solid #b8860b;padding:10px;",
                     tags$b("이 탭이 모델의 핵심 주장입니다."), br(),
                     "세 상행 전달계는 모두 같은 3-인자 곱을 통과합니다 — ",
                     tags$code("전시냅스 용량 × 후시냅스 밀도 × (1 − 점유율)"), ". ",
                     "다른 것은 가운데 인자의 미분방정식뿐이며, 그 부호가 계마다 다릅니다."),
                 br(),
                 plotOutput("p_transducer", height = 300),
                 h4("후시냅스 수용체 밀도 — 부호가 다른 세 개의 방정식"),
                 plotOutput("p_receptor", height = 300),
                 tableOutput("tbl_transducer")),

        ## --------------------------------------------- 4. cognition & CAF
        tabPanel("4 · 인지와 변동 ★",
                 br(), plotOutput("p_cog", height = 340),
                 h4("주의 상태의 퍼텐셜 지형 (attention potential landscape)"),
                 plotOutput("p_potential", height = 340),
                 helpText("왼쪽: 선택한 시점에서의 3차 퍼텐셜. 우물이 두 개면 양안정 영역 안에 있고, ",
                          "이때 상태 전환이 일어나므로 임상에서 '인지 변동'으로 채점됩니다. ",
                          "오른쪽: 양안정 깊이(BISTAB)와 실제 CAF 점수. ",
                          "ChEI는 구동력을 올려 환자를 이 띠 밖으로 밀어내며, ",
                          "그래서 평균(MMSE)보다 분산(CAF)을 먼저, 더 크게 움직입니다.")),

        ## -------------------------------------------------- 5. hallucinations
        tabPanel("5 · 환시 (PAD 모형)",
                 br(), plotOutput("p_vh", height = 340),
                 h4("세 인자의 곱 분해 (factor decomposition)"),
                 plotOutput("p_pad", height = 340),
                 helpText("환시 구동력은 합이 아니라 ", tags$b("곱"), "입니다: ",
                          tags$code("(1−상향 증거 충실도) × (1−하향 주의 결속) × 5-HT2A 이득"), ". ",
                          "곱이기 때문에 한 인자만으로는 환시를 설명할 수도, 없앨 수도 없습니다 — ",
                          "피마반세린이 환시를 없애는 게 아니라 약 1/3만 줄이는 이유입니다.")),

        ## ---------------------------------- 6. motor & neuroleptic sensitivity
        tabPanel("6 · 운동과 신경이완제 민감성 ★",
                 br(), plotOutput("p_motor", height = 340),
                 h4("민감성의 인과 사슬 한 화면 (cause → consequence)"),
                 plotOutput("p_nsens", height = 360),
                 helpText("D2 점유율은 세 표현형에서 동일합니다. 다른 것은 ",
                          tags$b("후시냅스 예비능(RESERVE)"), " 뿐이고, ",
                          "약물이 만드는 결손(DEFICIT)은 예비능이 적을수록 증폭됩니다. ",
                          "무투약 상태에서는 점유율이 0이므로 결손도 0 — ",
                          "즉 민감성 반응은 질환 중증도의 이정표가 아니라 ", tags$b("약물 사건"), "입니다.")),

        ## ---------------------------------------------- 7. pathology & GBA1
        tabPanel("7 · 병리 진행과 GBA1",
                 br(), plotOutput("p_path", height = 340),
                 h4("GBA1 양성 되먹임 고리"),
                 plotOutput("p_gba", height = 340),
                 helpText("GCase↓ → GlcCer↑ → 올리고머 분해 지연 → 올리고머↑ → GCase 수송 차단 → GCase↓. ",
                          "기본 파라미터에서 이 고리의 이득은 1보다 작으므로 ",
                          tags$b("단안정(monostable)"), "입니다 — 앰브록솔을 끊으면 되돌아갑니다. ",
                          "따라서 '기한을 놓치면 끝'이 아니라 '일찍 걸수록 지렛대가 길다'가 ",
                          "이 모델이 실제로 지지하는 주장입니다. ",
                          "앰브록솔 시작일을 바꿔가며 확인해 보십시오.")),

        ## ------------------------------------------------------ 8. biomarkers
        tabPanel("8 · 바이오마커",
                 br(), plotOutput("p_bio", height = 480),
                 tableOutput("tbl_bio")),

        ## ------------------------------------------------ 9. scenario compare
        tabPanel("9 · 시나리오 비교",
                 br(),
                 helpText("모델 파일에 정의된 23개 시나리오를 모두 실행하여 비교합니다. ",
                          "처음 실행 시 30초 내외가 걸립니다."),
                 actionButton("go_all", "전체 시나리오 실행 (run all)", class = "btn-warning"),
                 br(), br(),
                 tableOutput("tbl_all"),
                 plotOutput("p_all", height = 460)),

        ## --------------------------------------------------------- 10. survival
        tabPanel("10 · 생존과 위험",
                 br(), plotOutput("p_surv", height = 400),
                 tableOutput("tbl_surv"),
                 helpText("위험함수는 인지 손실·운동 점수·낙상 지수에 대해 로그선형이고, ",
                          "신경이완제 민감성은 ", tags$b("곱셈적"), "으로 들어갑니다. ",
                          "항정신병약 노출 후 사망률이 2–3배가 되는 신호를 재현하는 것은 ",
                          "바로 이 곱셈 구조입니다.")),

        ## ------------------------------------------------------ 11. model card
        tabPanel("11 · 모델 카드",
                 br(), htmlOutput("model_card"))
      )
    )
  )
)

## =============================================================================
##  SERVER
## =============================================================================
server <- function(input, output, session) {

  ## ------------------------------------------------------------ dose builder
  build_dose <- function() {
    end <- input$years * 365
    ev_list <- list()
    add <- function(e) ev_list[[length(ev_list) + 1]] <<- e

    st <- input$chei_start
    if (input$chei == "riv_oral")  add(riv_oral(6,   st, max(end - st, 1)))
    if (input$chei == "riv_patch") add(riv_patch(9.5, st, max(end - st, 1)))
    if (input$chei == "don")       add(don_oral(10,  st, max(end - st, 1)))
    if (input$mem)  add(mem_oral(20, st, max(end - st, 1)))
    if (input$pim)  add(pim_oral(34, 365, max(end - 365, 1)))
    if (input$ld)   add(ld_oral(150, 365, max(end - 365, 1)))
    if (input$zon)  add(zon_oral(25, 365, max(end - 365, 1)))
    if (input$amb)  add(amb_oral(420, input$amb_start,
                                 max(end - input$amb_start, 1)))
    if (input$mab)  add(mab_iv(4500, 180, max(floor((end - 180)/28), 1)))
    if (input$apd == "risp") add(apd_oral(1,  730, 180))
    if (input$apd == "quet") add(apd_oral(50, 730, 180))

    if (!length(ev_list)) return(NULL)
    Reduce(c, ev_list)
  }

  build_param <- function() {
    p <- list(PHENO  = as.numeric(input$pheno),
              GBAF   = input$gbaf,
              APOE4  = input$apoe4,
              SNCADOSE = input$snca,
              ANTICH = input$antich,
              MELON  = as.numeric(input$mel),
              DROXON = as.numeric(input$drox),
              MODON  = as.numeric(input$mod))
    if (input$apd == "quet")
      p <- c(p, list(EC50D2 = 900, CLAPD = 1400, VAPD = 700,
                     APDA1 = 1.0, APDSED = 1.0))
    p
  }

  sim <- eventReactive(input$go, {
    end <- input$years * 365
    m   <- param(mod, build_param())
    d   <- build_dose()
    trt <- if (is.null(d)) mrgsim_df(m, end = end, delta = 1)
           else            mrgsim_df(m, events = d, end = end, delta = 1)
    trt$scenario <- "treated"
    if (isTRUE(input$show_ref)) {
      ref <- mrgsim_df(param(mod, build_param()), end = end, delta = 1)
      ref$scenario <- "untreated"
      rbind(ref, trt)
    } else trt
  }, ignoreNULL = FALSE)

  ph_label <- reactive(c("0" = "DLB", "1" = "PD-dementia", "2" = "AD")[input$pheno])

  ## ----------------------------------------------------------- 1. profile --
  output$tbl_profile <- renderTable({
    d <- sim(); d0 <- d[d$time == 0 & d$scenario == d$scenario[1], ][1, ]
    data.frame(
      `항목` = c("표현형", "MMSE", "MDS-UPDRS III", "CAF (인지 변동)",
               "NPI 총점", "환시 부담", "RBD 중증도", "기립성 SBP 강하 (mmHg)",
               "Epworth 졸림", "MIBG H/M", "DaTSCAN 결합비", "EEG 우세주파수 (Hz)",
               "신피질 섬유 부하", "변연계 섬유 부하", "뇌간 섬유 부하"),
      `값` = c(ph_label(),
             sprintf("%.1f", d0$MMSE), sprintf("%.1f", d0$MOT), sprintf("%.1f", d0$CAF),
             sprintf("%.1f", d0$NPI), sprintf("%.2f", d0$VHB), sprintf("%.1f", d0$RBDS),
             sprintf("%.1f", d0$AUTS), sprintf("%.1f", d0$EDSS), sprintf("%.2f", d0$MIBG),
             sprintf("%.2f", d0$DATSBR), sprintf("%.2f", d0$EEGF),
             sprintf("%.3f", d0$FIBN), sprintf("%.3f", d0$FIBL), sprintf("%.3f", d0$FIBB)),
      check.names = FALSE)
  })

  output$p_profile <- renderPlot({
    trace_plot(sim(), c("MIBG", "DATSBR", "EEGF", "RBDS", "AUTS", "EDSS"),
               c("MIBG H/M", "DaTSCAN 결합비", "EEG 우세주파수 (Hz)",
                 "RBD 중증도", "기립성 SBP 강하 (mmHg)", "Epworth"),
               title = "지표 바이오마커와 지지 특징의 시간 경과",
               subtitle = paste0("표현형: ", ph_label())) +
      scale_colour_manual(values = PAL)
  })

  ## ---------------------------------------------------------------- 2. PK --
  output$p_pk <- renderPlot({
    d <- sim(); d <- d[d$scenario == "treated", ]
    trace_plot(d, c("CRIV", "CDON", "CPIM", "CLD", "CAPD", "CAMBB"),
               c("리바스티그민 (ng/mL)", "도네페질 (ng/mL)", "피마반세린 (ng/mL)",
                 "레보도파 (ng/mL)", "항정신병약 (ng/mL)", "뇌 앰브록솔"),
               title = "약물 농도 시간 경과",
               subtitle = "출력 격자는 1일 간격이므로 반감기가 짧은 약물의 봉우리는 표시되지 않을 수 있습니다") +
      scale_colour_manual(values = PAL)
  })

  output$p_carb <- renderPlot({
    d <- sim(); d <- d[d$scenario == "treated", ]
    trace_plot(d, c("INHACHE", "INHBCHE", "CHDRIVE"),
               c("AChE 억제율", "BuChE 억제율", "콜린성 전달자 출력"),
               title = "카바밀화 효소 — 약력학 반감기는 효소 재합성이 결정한다",
               hlines = 0.6) + scale_colour_manual(values = PAL)
  })

  ## ------------------------------------------------------- 3. transducers --
  output$p_transducer <- renderPlot({
    trace_plot(sim(), c("CHDRIVE", "DADRIVE", "HT2ASIG"),
               c("콜린성 (ACh × M1)", "도파민성 (DA × D2 × 후시냅스)",
                 "세로토닌성 (5-HT × 5-HT2A)"),
               title = "세 전달자의 출력 — 같은 곱, 다른 궤적",
               subtitle = "세 계 모두 '전시냅스 × 후시냅스 × (1 − 점유율)'") +
      scale_colour_manual(values = PAL)
  })

  output$p_receptor <- renderPlot({
    trace_plot(sim(), c("M1R", "D2R", "HT2A"),
               c("M1/M4 (tau에만 반응 → DLB에서 보존)",
                 "D2 (변연계 α-syn이 상향조절을 막음)",
                 "5-HT2A (탈신경 + 신피질 섬유 → 상승)"),
               title = "후시냅스 수용체 밀도 — 부호가 서로 다른 세 개의 방정식",
               subtitle = "ChEI 우월성 · 레보도파 실패 · 신경이완제 민감성 · 피마반세린 효능은 모두 이 그림의 귀결") +
      scale_colour_manual(values = PAL)
  })

  output$tbl_transducer <- renderTable({
    d <- sim(); d <- d[d$scenario == tail(unique(d$scenario), 1), ]
    idx <- sapply(c(0, 1, 3, 5) * 365.25, function(t) which.min(abs(d$time - t)))
    data.frame(
      `연차` = c(0, 1, 3, 5),
      `ACh × M1` = round(d$CHDRIVE[idx], 3),
      `M1 밀도`  = round(d$M1R[idx], 3),
      `DA 신호`  = round(d$DADRIVE[idx], 3),
      `D2 밀도`  = round(d$D2R[idx], 3),
      `후시냅스 예비능` = round(d$RESERVE[idx], 3),
      `5-HT2A 신호` = round(d$HT2ASIG[idx], 3),
      `5-HT2A 밀도` = round(d$HT2A[idx], 3),
      check.names = FALSE)
  })

  ## --------------------------------------------------------- 4. cognition --
  output$p_cog <- renderPlot({
    trace_plot(sim(), c("MMSE", "CAF", "ATTM", "BISTAB"),
               c("MMSE", "CAF (인지 변동, 0–16)", "주의 상태 ATTM", "양안정 깊이 BISTAB"),
               title = "인지: 평균과 분산은 다른 것을 측정한다") +
      scale_colour_manual(values = PAL)
  })

  output$p_potential <- renderPlot({
    d <- sim()
    dt <- d[d$scenario == tail(unique(d$scenario), 1), ]
    pick <- c(0.5, 1.5, 3, 5)
    A <- seq(-1.6, 1.6, by = 0.01)
    alpha <- as.numeric(param(mod)$ALPHAB)
    off   <- as.numeric(param(mod)$ATTOFF)
    grid <- do.call(rbind, lapply(pick, function(t) {
      i <- which.min(abs(dt$time - t * 365.25))
      D <- dt$DRIVEA[i] - off
      data.frame(A = A,
                 V = A^4/4 - alpha*A^2/2 - D*A,
                 lab = sprintf("%.1f년  (구동력 %.3f)", t, dt$DRIVEA[i]))
    }))
    grid <- grid %>% group_by(lab) %>% mutate(V = V - min(V)) %>% ungroup()
    ggplot(grid, aes(A, V, colour = lab)) +
      geom_line(linewidth = 0.8) +
      labs(title = "주의 상태의 퍼텐셜 V(A) = A⁴/4 − αA²/2 − D·A",
           subtitle = "우물이 두 개인 동안 상태 전환이 가능하고, 그것이 임상적 '인지 변동'이다",
           x = "주의 상태 A", y = "퍼텐셜 (상대값)") +
      qsp_theme
  })

  ## ---------------------------------------------------- 5. hallucinations --
  output$p_vh <- renderPlot({
    trace_plot(sim(), c("VHB", "VHDRIVE", "HT2ASIG", "OCC2A"),
               c("환시 부담 (NPI-hall)", "환시 구동력 (세 인자의 곱)",
                 "5-HT2A 신호", "피마반세린 5-HT2A 점유율"),
               title = "환시") + scale_colour_manual(values = PAL)
  })

  output$p_pad <- renderPlot({
    d <- sim()
    long <- d %>%
      mutate(`상향 결손 (1−BU)` = 1 - BOTTOMUP,
             `하향 결손 (1−TD)` = 1 - TOPDOWN,
             `5-HT2A 이득`      = HT2ASIG,
             yrs = time/365.25) %>%
      select(yrs, scenario, `상향 결손 (1−BU)`, `하향 결손 (1−TD)`, `5-HT2A 이득`) %>%
      pivot_longer(-c(yrs, scenario))
    ggplot(long, aes(yrs, value, colour = scenario)) +
      geom_line(linewidth = 0.75) + facet_wrap(~name, scales = "free_y") +
      labs(title = "PAD 곱의 세 인자",
           subtitle = "세 인자가 동시에 커져야 환시가 나온다 — 그래서 단일 기전 설명도, 단일 기전 치료도 실패한다",
           x = "진단 후 경과 연수", y = "") +
      scale_colour_manual(values = PAL) + qsp_theme
  })

  ## ------------------------------------------------------------- 6. motor --
  output$p_motor <- renderPlot({
    trace_plot(sim(), c("MOT", "DADRIVE", "OCCD2", "NSENS"),
               c("MDS-UPDRS III", "도파민성 신호", "선조체 D2 점유율",
                 "신경이완제 민감성 지수"),
               title = "운동 증상") + scale_colour_manual(values = PAL)
  })

  output$p_nsens <- renderPlot({
    end <- input$years * 365
    d <- do.call(rbind, lapply(0:2, function(ph) {
      m <- param(mod, c(build_param(), list(PHENO = ph)))
      x <- mrgsim_df(m, events = if (input$apd == "quet") apd_oral(50, 730, 180)
                                 else apd_oral(1, 730, 180),
                     end = min(end, 1460), delta = 1)
      x$scenario <- c("DLB", "PD-dementia", "AD")[ph + 1]; x
    }))
    trace_plot(d, c("OCCD2", "RESERVE", "DEFICIT", "NSENS", "MOT", "SURV"),
               c("D2 점유율 (세 군 동일)", "후시냅스 예비능 (다르다)",
                 "약물 유발 결손 = 점유율 × 예비능 증폭", "민감성 지수",
                 "MDS-UPDRS III", "생존 확률"),
               title = "동일한 리스페리돈 노출, 세 개의 수용체 지형, 세 개의 결과",
               subtitle = "PHENO 외에는 어떤 것도 바뀌지 않았고, PHENO는 초기 병리 분포만 설정한다") +
      scale_colour_manual(values = PAL)
  })

  ## --------------------------------------------------------- 7. pathology --
  output$p_path <- renderPlot({
    trace_plot(sim(), c("FIBB", "FIBL", "FIBN", "SEED"),
               c("뇌간 섬유", "변연계 섬유", "신피질 섬유", "간질 seed 풀"),
               title = "지역별 병기 진행",
               subtitle = "하류 지역의 흡수는 상류 지역의 부하로 게이팅되어 순차적 병기가 나온다") +
      scale_colour_manual(values = PAL)
  })

  output$p_gba <- renderPlot({
    trace_plot(sim(), c("GCASE", "GLCCER", "OLIGTOT", "GCTGT"),
               c("GCase 활성", "글루코실세라마이드", "총 올리고머", "GCase 수송 목표"),
               title = "GBA1 양성 되먹임 고리",
               subtitle = "앰브록솔을 중단하면 되돌아온다 = 단안정. '기한'이 아니라 '지렛대 길이'의 문제") +
      scale_colour_manual(values = PAL)
  })

  ## -------------------------------------------------------- 8. biomarkers --
  output$p_bio <- renderPlot({
    trace_plot(sim(), c("MIBG", "DATSBR", "EEGF", "PTAU", "ABETA", "QTC"),
               c("MIBG H/M", "DaTSCAN 결합비", "EEG 우세주파수 (Hz)",
                 "tau 병리 (au)", "아밀로이드 (Centiloid)", "QTc 변화 (ms)"),
               title = "바이오마커") + scale_colour_manual(values = PAL)
  })

  output$tbl_bio <- renderTable({
    d <- sim(); d <- d[d$scenario == tail(unique(d$scenario), 1), ]
    idx <- sapply(c(0, 2, 4, 6) * 365.25, function(t) which.min(abs(d$time - t)))
    data.frame(`연차` = c(0, 2, 4, 6),
               `MIBG H/M` = round(d$MIBG[idx], 2),
               `DaTSCAN`  = round(d$DATSBR[idx], 2),
               `EEG Hz`   = round(d$EEGF[idx], 2),
               `p-tau`    = round(d$PTAU[idx], 3),
               `Centiloid`= round(d$ABETA[idx], 1),
               `QTc (ms)` = round(d$QTC[idx], 1),
               check.names = FALSE)
  })

  ## -------------------------------------------------- 9. all the scenarios --
  all_sc <- eventReactive(input$go_all, {
    withProgress(message = "23개 시나리오 실행 중...", value = 0, {
      run_all()
    })
  })

  output$tbl_all <- renderTable({ summarise_scn(all_sc()) })

  output$p_all <- renderPlot({
    d <- all_sc()
    ggplot(d, aes(time/365.25, MMSE, colour = scenario)) +
      geom_line(linewidth = 0.6) +
      labs(title = "모든 시나리오의 MMSE 궤적",
           x = "진단 후 경과 연수", y = "MMSE") +
      qsp_theme + theme(legend.position = "right",
                        legend.text = element_text(size = 7)) +
      guides(colour = guide_legend(ncol = 1))
  })

  ## --------------------------------------------------------- 10. survival --
  output$p_surv <- renderPlot({
    trace_plot(sim(), c("SURV", "HAZ", "FALLS", "CUMH"),
               c("생존 확률", "순간 위험", "낙상 지수", "누적 위험"),
               title = "생존", hlines = NULL) + scale_colour_manual(values = PAL)
  })

  output$tbl_surv <- renderTable({
    d <- sim()
    do.call(rbind, lapply(split(d, d$scenario), function(x) {
      i <- which(x$SURV <= 0.5)[1]
      data.frame(`시나리오` = x$scenario[1],
                 `중간생존 (년)` = if (is.na(i)) NA_real_ else round(x$time[i]/365.25, 2),
                 `5년 생존` = round(x$SURV[which.min(abs(x$time - 1825))], 3),
                 check.names = FALSE)
    }))
  })

  ## ------------------------------------------------------- 11. model card --
  output$model_card <- renderUI({
    HTML('
<h3>모델 카드 (Model card)</h3>
<table border="0" cellpadding="6">
<tr><td><b>구획 수</b></td><td>66 ODE</td></tr>
<tr><td><b>시간 단위</b></td><td>일 (days)</td></tr>
<tr><td><b>적분기</b></td><td>LSODA (mrgsolve), rtol 1e-8 / atol 1e-8</td></tr>
<tr><td><b>약물</b></td><td>리바스티그민(경구·패치), 도네페질, 갈란타민, 메만틴,
피마반세린(+AC-279), 레보도파/카비도파, 일반 항정신병약 구획, 조노사마이드,
앰브록솔, 항 α-시누클레인 항체, 드록시도파, 멜라토닌/클로나제팜, 모다피닐</td></tr>
<tr><td><b>표현형</b></td><td>DLB · PD-dementia · AD — <i>초기 병리 분포만</i> 다르고 방정식은 동일</td></tr>
<tr><td><b>개인간 변이</b></td><td>진단 시점의 nbM/PPN·SNc/LC 잔존량, 변연계·신피질 부하, 아밀로이드</td></tr>
</table>

<h4>모델이 <b>주장</b>하는 것</h4>
<ol>
<li>세 전달계가 하나의 전달자를 공유하며, 후시냅스 가소성의 <b>부호</b>만 다르다.
ChEI 우월성·레보도파 실패·신경이완제 민감성·피마반세린 효능은 각각의 규칙이 아니라
이 구조의 <b>결과</b>다.</li>
<li>인지 변동은 중증도가 아니라 <b>분산</b>이다 — 양안정 주의 상태의 전환 빈도.
따라서 ChEI는 MMSE보다 변동 지표를 비례적으로 더 크게 개선한다(모델 계산: 약 3.5배).</li>
<li>환시는 세 결손의 <b>곱</b>이다. 그래서 단일 기전으로 설명되지도, 없어지지도 않는다.</li>
<li>신경이완제 민감성은 <b>약물 사건</b>이지 중증도 이정표가 아니다 —
무투약 상태에서는 D2 점유율이 0이므로 결손도 0이다.</li>
</ol>

<h4>모델이 <b>주장하지 않는</b> 것</h4>
<ul>
<li><b>GBA1 고리는 양안정이 아니다.</b> 초안에서는 안장-절점 분기를 넣으려 했으나,
문헌에서 지지되는 범위의 파라미터로는 고리 이득이 1에 못 미쳤다. 앰브록솔을 중단하면
GCase가 되돌아온다(탭 7에서 직접 확인 가능). 따라서 "기한을 놓치면 끝"이 아니라
"일찍 시작할수록 지렛대가 길다"가 이 모델이 지지하는 주장이다.</li>
<li>신경 소실 속도상수들은 개별적으로 측정된 값이 아니라 자연경과 기준점에
동시에 맞춘 값이며 서로 상관되어 있다. 개별 값에 생물학적 의미를 두면 안 된다.</li>
<li>어떤 예측도 전향적으로 검증되지 않았다.</li>
</ul>

<p style="color:#a00;"><b>임상 사용 금지.</b> 교육 및 연구 목적의 반정량적 모델입니다.</p>
')
  })
}

shinyApp(ui, server)
