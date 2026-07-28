## =============================================================================
## CRPS QSP Shiny dashboard
## 복합부위통증증후군 QSP 모델 인터랙티브 대시보드
##
##   library(shiny); library(mrgsolve); library(ggplot2); library(dplyr)
##   shiny::runApp("crps_shiny_app.R")
##
## The app loads crps_mrgsolve_model.R from the same directory and exposes the
## model's structure rather than only its output: the left panel separates
## PATIENT covariates (trait, sympathetic tone, insult, cast) from TREATMENT
## choices (arms, start day, intensity), because in this model the interaction
## between those two blocks -- not the potency of any arm -- is what decides
## the outcome.
##
## Tabs (8):
##   1. 환자 프로파일          patient profile, attractor read-out, ring gain
##   2. 약물 PK                all drug concentrations on one time axis
##   3. 말초 노드 (PD)         SP/CGRP, cytokines, NGF, ROS, oedema, alpha1
##   4. 중추 노드 · 글리아 잠금 spinal, glial latch, descending tone, cortex
##   5. 임상 엔드포인트        NRS, CSS, ROM, temperature asymmetry
##   6. 치료 창 (window)       start-day scan -> t*, with the arm decomposition
##   7. 시나리오 비교          9 prebuilt scenarios side by side
##   8. 바이오마커 · 골        CTX-I, BMD, blister-fluid-like cytokine panel
## =============================================================================

library(shiny)
library(mrgsolve)
library(ggplot2)
library(dplyr)
library(tidyr)

MODEL_FILE <- "crps_mrgsolve_model.R"
mod <- mrgsolve::mread_cache("crps", MODEL_FILE)
E   <- mod@envir            # event builders + analysis functions
DAY <- 24

STATE_LABELS <- c(
  NP = "SP/CGRP 신경펩타이드", CYT = "사이토카인 (IL-6/TNF/IL-1b)",
  NGF = "NGF", EDEMA = "부종", ROS = "산화 스트레스 (ROS)",
  AAB = "자가항체 역가", ALPHA1 = "alpha1-수용체 상향조절",
  PSENS = "말초 감작", PERF = "관류 지수", HYPOX = "조직 저산소/산증",
  SSENS = "척수 중추 감작", GLIA = "교세포 활성 (잠금 변수)",
  DINH = "하행 억제 톤", CORTEX = "피질 지도 열화", DISUSE = "불용/회피",
  PAIN = "통증 NRS", ROM = "관절 운동범위", OC = "파골세포 활성",
  BMD = "국소 BMD", CTXI = "혈청 CTX-I")

theme_crps <- function() {
  theme_minimal(base_size = 13) +
    theme(panel.grid.minor = element_blank(),
          strip.text = element_text(face = "bold"),
          legend.position = "bottom")
}

## ---------------------------------------------------------------------------
## UI
## ---------------------------------------------------------------------------
ui <- fluidPage(
  titlePanel("복합부위통증증후군 (CRPS) QSP 시뮬레이터 — 말초 노드 · 행동-피질 고리 · 교세포 잠금"),
  sidebarLayout(
    sidebarPanel(
      width = 3,
      h4("환자 (Patient covariates)"),
      sliderInput("kfear", "공포-회피 이득 KFEAR (기질)",
                  min = 0.2, max = 1.6, value = 1.0, step = 0.05),
      sliderInput("symp", "교감신경 톤 SYMP_TONE",
                  min = 0.5, max = 2.0, value = 1.0, step = 0.1),
      sliderInput("inj", "유발사건 강도 INJ_AMP",
                  min = 0.1, max = 1.6, value = 1.0, step = 0.05),
      sliderInput("immob", "석고고정 기간 (일)",
                  min = 0, max = 90, value = 30, step = 5),
      hr(),
      h4("치료 (Treatment)"),
      sliderInput("start", "치료 시작일 (day)",
                  min = 3, max = 365, value = 7, step = 1),
      checkboxGroupInput("arms", "치료 arm",
                         choices = c("프레드니솔론 (경구 테이퍼)" = "pred",
                                     "NAC / DMSO (항산화)" = "nac",
                                     "재활 · GMI · 단계적 노출" = "rehab",
                                     "네리드로네이트 100mg IV x4" = "nerid",
                                     "케타민 100시간 주입" = "ket",
                                     "가바펜틴 600mg TID" = "gbp",
                                     "아미트립틸린 50mg" = "amt",
                                     "IVIG 0.5 g/kg x6" = "ivig",
                                     "척수자극 (SCS)" = "scs",
                                     "교감신경 차단" = "sb",
                                     "혈관확장제 (PDE5i 유사)" = "vaso"),
                         selected = c("pred", "nac", "rehab")),
      sliderInput("rehab_int", "재활 강도/순응도 (0-1)",
                  min = 0, max = 1, value = 1, step = 0.05),
      sliderInput("ketrate", "케타민 최대 주입속도 (mg/h)",
                  min = 5, max = 90, value = 22, step = 1),
      sliderInput("years", "시뮬레이션 기간 (년)",
                  min = 1, max = 5, value = 3, step = 1),
      hr(),
      helpText("치료를 켠 뒤 '치료 시작일'만 움직여 보십시오. 동일한 처방이 ",
               "90일 이전에는 완전 소실을, 95일 이후에는 무효를 보이는 것이 ",
               "이 모델의 핵심 결과입니다.")
    ),
    mainPanel(
      width = 9,
      tabsetPanel(
        id = "tabs",
        tabPanel("1. 환자 프로파일",
                 fluidRow(column(4, wellPanel(h4("3년 후 상태"), tableOutput("attractor"))),
                          column(8, plotOutput("profile_pain", height = "300px"))),
                 plotOutput("profile_ring", height = "260px")),
        tabPanel("2. 약물 PK", plotOutput("pk", height = "620px")),
        tabPanel("3. 말초 노드 (PD)", plotOutput("periph", height = "620px")),
        tabPanel("4. 중추 노드 · 교세포 잠금",
                 plotOutput("central", height = "420px"),
                 plotOutput("latch", height = "250px")),
        tabPanel("5. 임상 엔드포인트", plotOutput("endpoints", height = "620px")),
        tabPanel("6. 치료 창 (window)",
                 fluidRow(column(6, plotOutput("window", height = "380px")),
                          column(6, plotOutput("armscan", height = "380px"))),
                 tableOutput("window_tab")),
        tabPanel("7. 시나리오 비교",
                 plotOutput("scen_plot", height = "420px"),
                 tableOutput("scen_tab")),
        tabPanel("8. 바이오마커 · 골", plotOutput("bio", height = "620px"))
      )
    )
  )
)

## ---------------------------------------------------------------------------
## server
## ---------------------------------------------------------------------------
server <- function(input, output, session) {

  build_events <- function(arms, start_day, ketrate) {
    st <- start_day * DAY
    evs <- list()
    if ("pred"  %in% arms) evs <- c(evs, list(E$ev_prednisolone(st)))
    if ("nac"   %in% arms) evs <- c(evs, list(E$ev_nac(st, days = 180)))
    if ("nerid" %in% arms) evs <- c(evs, list(E$ev_neridronate(st)))
    if ("ket"   %in% arms) evs <- c(evs, list(E$ev_ketamine(st, rate_max = ketrate)))
    if ("gbp"   %in% arms) evs <- c(evs, list(E$ev_gabapentin(st, days = 180)))
    if ("amt"   %in% arms) evs <- c(evs, list(E$ev_amitriptyline(st, days = 180)))
    if ("ivig"  %in% arms) evs <- c(evs, list(E$ev_ivig(st)))
    if (!length(evs)) return(NULL)
    do.call(E$comb_ev, evs)
  }

  build_pars <- function(input) {
    st <- input$start * DAY
    p <- list(KFEAR = input$kfear, SYMP_TONE = input$symp,
              INJ_AMP = input$inj, IMMOB_DUR = input$immob * DAY)
    if ("rehab" %in% input$arms)
      p <- c(p, list(REHAB = input$rehab_int, REHAB_T0 = st,
                     REHAB_DUR = 180 * DAY))
    if ("scs"  %in% input$arms) p <- c(p, list(SCS_ON = 1, SCS_T0 = st))
    if ("sb"   %in% input$arms) p <- c(p, list(SYMPBLOCK = 1, SB_T0 = st))
    if ("vaso" %in% input$arms) p <- c(p, list(VASODIL = 1))
    p
  }

  simulate <- reactive({
    end <- input$years * 8760
    E$sim(mod, build_pars(input),
          build_events(input$arms, input$start, input$ketrate),
          end = end, delta = 6)
  })

  untreated <- reactive({
    end <- input$years * 8760
    E$sim(mod, list(KFEAR = input$kfear, SYMP_TONE = input$symp,
                    INJ_AMP = input$inj, IMMOB_DUR = input$immob * DAY),
          NULL, end = end, delta = 6)
  })

  long <- function(d, cols) {
    d %>% select(all_of(c("time", cols))) %>%
      mutate(day = time / 24) %>% select(-time) %>%
      pivot_longer(-day, names_to = "state", values_to = "value") %>%
      mutate(label = ifelse(state %in% names(STATE_LABELS),
                            paste0(state, " — ", STATE_LABELS[state]), state))
  }

  ## ---- tab 1 -------------------------------------------------------------
  output$attractor <- renderTable({
    f <- E$summarise_run(simulate())
    u <- E$summarise_run(untreated())
    tab <- data.frame(
      metric   = c("NRS", "CSS (0-16)", "ROM", "BMD", "GLIA", "CORTEX",
                   "DISUSE", "TEMP_ASYM (C)", "latched", "remission"),
      treated  = c(f$NRS, f$CSS, f$ROM, f$BMD, f$GLIA, f$CORTEX,
                   f$DISUSE, f$TEMP, f$LATCHED, f$REMISSION),
      untreated = c(u$NRS, u$CSS, u$ROM, u$BMD, u$GLIA, u$CORTEX,
                    u$DISUSE, u$TEMP, u$LATCHED, u$REMISSION))
    names(tab) <- c("\uc9c0\ud45c", "\uce58\ub8cc", "\ubb34\uce58\ub8cc")
    tab
  }, digits = 3)

  output$profile_pain <- renderPlot({
    d <- bind_rows(mutate(simulate(), arm = "치료"),
                   mutate(untreated(), arm = "무치료"))
    ggplot(d, aes(time / 24, PAIN, colour = arm)) +
      geom_line(linewidth = 0.9) +
      geom_vline(xintercept = input$start, linetype = 2, colour = "grey40") +
      labs(x = "일 (day)", y = "통증 NRS (0-10)", colour = NULL,
           title = "통증 경과 — 점선 = 치료 시작") +
      scale_colour_manual(values = c("치료" = "#1565c0", "무치료" = "#b71c1c")) +
      theme_crps()
  })

  output$profile_ring <- renderPlot({
    d <- simulate()
    ggplot(d, aes(time / 24, RING_GAIN)) +
      geom_hline(yintercept = 1, linetype = 2, colour = "#b71c1c") +
      geom_line(colour = "#6a1b9a", linewidth = 0.9) +
      labs(x = "일", y = "고리 이득 (ring gain)",
           title = "통증→불용→피질→통증 고리의 국소 이득 (>1 이면 자기유지 가능)") +
      theme_crps()
  })

  ## ---- tab 2 -------------------------------------------------------------
  output$pk <- renderPlot({
    d <- simulate() %>% mutate(day = time / 24) %>%
      select(day, `케타민 (ng/mL)` = KET_ng, `노르케타민 (ng/mL)` = NORK_ng,
             `프레드니솔론 (ng/mL)` = PRED_ng, `가바펜틴 (ng/mL)` = GBP_ng,
             `아미트립틸린 (ng/mL)` = AMT_ng, `NAC (ng/mL)` = NAC_ng,
             `IVIG (g/L)` = IVIG_gL, `네리드로네이트 혈장 (mg/L)` = NER_mgL) %>%
      pivot_longer(-day)
    ggplot(d, aes(day, value)) + geom_line(colour = "#1565c0") +
      facet_wrap(~name, scales = "free_y", ncol = 2) +
      labs(x = "일", y = NULL, title = "약물 농도-시간 프로파일") + theme_crps()
  })

  ## ---- tab 3 -------------------------------------------------------------
  output$periph <- renderPlot({
    ggplot(long(simulate(), c("NP", "CYT", "NGF", "EDEMA", "ROS", "AAB",
                              "ALPHA1", "PSENS", "HYPOX")),
           aes(day, value)) +
      geom_line(colour = "#c62828") +
      facet_wrap(~label, scales = "free_y", ncol = 3) +
      labs(x = "일", y = "활성 지수 (0-1)",
           title = "빠른 말초 노드 — 유발사건에 의해 구동되고 스스로 소멸한다") +
      theme_crps()
  })

  ## ---- tab 4 -------------------------------------------------------------
  output$central <- renderPlot({
    ggplot(long(simulate(), c("SSENS", "GLIA", "DINH", "CORTEX", "DISUSE")),
           aes(day, value)) +
      geom_line(colour = "#4527a0") +
      facet_wrap(~label, scales = "free_y", ncol = 3) +
      labs(x = "일", y = "지수", title = "느린 중추 고리 (수주-수개월 시상수)") +
      theme_crps()
  })

  output$latch <- renderPlot({
    g50 <- as.numeric(mrgsolve::param(mod)$GLIA50)
    d <- bind_rows(mutate(simulate(), arm = "치료"),
                   mutate(untreated(), arm = "무치료"))
    ggplot(d, aes(time / 24, GLIA, colour = arm)) +
      geom_hline(yintercept = g50, linetype = 2, colour = "#b71c1c") +
      geom_line(linewidth = 0.9) +
      annotate("text", x = Inf, y = g50, label = " GLIA50 (잠금 문턱)",
               hjust = 1.02, vjust = -0.6, size = 3.4, colour = "#b71c1c") +
      scale_colour_manual(values = c("치료" = "#1565c0", "무치료" = "#b71c1c")) +
      labs(x = "일", y = "GLIA", colour = NULL,
           title = "교세포 잠금 변수: 이 선을 넘으면 척수 감작이 구심 입력 없이 유지된다") +
      theme_crps()
  })

  ## ---- tab 5 -------------------------------------------------------------
  output$endpoints <- renderPlot({
    d <- simulate() %>% mutate(day = time / 24) %>%
      select(day, `통증 NRS` = PAIN, `CRPS 중증도 점수 (0-16)` = CSS,
             `관절 운동범위 (정상비)` = ROM,
             `온도 비대칭 (섭씨, +온난/-한랭)` = TEMP_ASYM,
             `부종 지수` = EDEMA, `활성 CRPS 플래그` = ACTIVE_CRPS) %>%
      pivot_longer(-day)
    ggplot(d, aes(day, value)) + geom_line(colour = "#00695c") +
      geom_vline(xintercept = input$start, linetype = 2, colour = "grey50") +
      facet_wrap(~name, scales = "free_y", ncol = 2) +
      labs(x = "일", y = NULL, title = "임상 엔드포인트") + theme_crps()
  })

  ## ---- tab 6 -------------------------------------------------------------
  window_res <- reactive({
    days <- c(3, 7, 14, 21, 30, 45, 60, 75, 90, 95, 105, 120, 180, 240, 365)
    base <- build_pars(input)
    arms <- intersect(input$arms, c("pred", "nac", "rehab"))
    if (!length(arms)) arms <- c("pred", "nac", "rehab")
    out <- lapply(days, function(dd) {
      st <- dd * DAY
      p <- base
      p$REHAB_T0 <- st
      if ("rehab" %in% arms) p$REHAB <- input$rehab_int else p$REHAB <- 0
      evs <- list()
      if ("pred" %in% arms) evs <- c(evs, list(E$ev_prednisolone(st)))
      if ("nac"  %in% arms) evs <- c(evs, list(E$ev_nac(st, days = 180)))
      ev <- if (length(evs)) do.call(E$comb_ev, evs) else NULL
      d <- E$sim(mod, p, ev, end = 3 * 8760, delta = 24)
      f <- E$final_row(d)
      data.frame(start_day = dd, NRS = f$PAIN, CSS = f$CSS,
                 GLIA = f$GLIA, BMD = f$BMD, ROM = f$ROM)
    })
    bind_rows(out)
  })

  output$window <- renderPlot({
    w <- window_res()
    ggplot(w, aes(start_day, NRS)) +
      geom_line(colour = "#2e7d32", linewidth = 1) +
      geom_point(colour = "#2e7d32") +
      labs(x = "치료 시작일 (day)", y = "3년 후 NRS",
           title = "치료 창: 동일 처방, 시작일만 다름") +
      theme_crps()
  })

  output$armscan <- renderPlot({
    days <- c(7, 30, 60, 120, 240)
    combos <- list(pred = "pred", nac = "nac", rehab = "rehab",
                   all = c("pred", "nac", "rehab"))
    res <- lapply(names(combos), function(nm) {
      sapply(days, function(dd) {
        p <- E$package(dd, combos[[nm]])
        d <- E$sim(mod, utils::modifyList(build_pars(input), p$pars),
                   p$events, end = 3 * 8760, delta = 24)
        E$final_row(d)$PAIN
      })
    })
    df <- data.frame(day = rep(days, length(combos)),
                     arm = rep(names(combos), each = length(days)),
                     NRS = unlist(res))
    ggplot(df, aes(day, NRS, colour = arm)) +
      geom_line(linewidth = 0.9) + geom_point() +
      labs(x = "치료 시작일", y = "3년 후 NRS", colour = "arm",
           title = "어느 arm이 창을 만드는가") + theme_crps()
  })

  output$window_tab <- renderTable(window_res(), digits = 3)

  ## ---- tab 7 -------------------------------------------------------------
  scen <- reactive(E$run_scenarios(mod, end = 3 * 8760))

  output$scen_plot <- renderPlot({
    s <- scen()
    d <- bind_rows(lapply(names(s), function(n)
      data.frame(day = s[[n]]$time / 24, NRS = s[[n]]$PAIN, scenario = n)))
    ggplot(d, aes(day, NRS, colour = scenario)) + geom_line(linewidth = 0.8) +
      labs(x = "일", y = "NRS", colour = NULL, title = "9개 시나리오") +
      theme_crps()
  })

  output$scen_tab <- renderTable(attr(scen(), "summary"), digits = 3)

  ## ---- tab 8 -------------------------------------------------------------
  output$bio <- renderPlot({
    d <- simulate() %>% mutate(day = time / 24) %>%
      select(day, `혈청 CTX-I (ng/mL)` = CTXI, `국소 BMD (대측비)` = BMD,
             `파골세포 활성` = OC, `사이토카인 지수 (수포액 대응)` = CYT,
             `자가항체 역가` = AAB, `산화 스트레스` = ROS) %>%
      pivot_longer(-day)
    ggplot(d, aes(day, value)) + geom_line(colour = "#6d4c41") +
      facet_wrap(~name, scales = "free_y", ncol = 2) +
      labs(x = "일", y = NULL, title = "바이오마커 및 골 축") + theme_crps()
  })
}

shinyApp(ui, server)
