## =============================================================================
##  INVASIVE PULMONARY ASPERGILLOSIS — QSP Shiny DASHBOARD
## =============================================================================
##  11 tabs over the 53-ODE model in `ipa_mrgsolve_model.R`.
##
##  The app is organised around one question the user should leave with:
##  WHICH LEVER ARE YOU ACTUALLY PULLING? Tab 4 lets you move the drug; tab 5
##  lets you move the marrow; tab 8 lets you move the clock. The mortality
##  readout is the same in all three, so the three can be compared in the only
##  currency that matters.
##
##  Run:  shiny::runApp("ipa_shiny_app.R")
##  Requires: shiny, mrgsolve, ggplot2, dplyr, tidyr, DT
## =============================================================================

library(shiny)
library(mrgsolve)
suppressMessages({library(ggplot2); library(dplyr); library(tidyr); library(DT)})

source("ipa_mrgsolve_model.R", local = TRUE)   # defines mod, rx_*, sim(), CYP2C19

THEME <- theme_bw(base_size = 12) +
  theme(panel.grid.minor = element_blank(),
        strip.background = element_rect(fill = "grey92", colour = NA),
        legend.position = "bottom")

PAL <- c("#2E86C1", "#E74C3C", "#28B463", "#B9770E", "#7D3C98", "#17A589",
         "#CB4335", "#5D6D7E")

## -----------------------------------------------------------------------------
## UI
## -----------------------------------------------------------------------------
ui <- fluidPage(
  titlePanel("침습성 폐 아스페르길루스증 (IPA) — QSP 대시보드"),
  tags$p(tags$b("한 문장: "),
         "트리아졸은 성장항에 곱해지고 호중구와 폴리엔만이 제거항에 더해진다. ",
         "따라서 호중구가 없는 숙주에서 아졸 단독으로는 균 부하가 결코 감소하지 않는다 ",
         "— 느려질 뿐이다."),
  hr(),
  sidebarLayout(
    sidebarPanel(
      width = 3,
      h4("숙주 (host)"),
      selectInput("phenotype", "숙주 표현형",
                  c("호중구감소증 (혈액암 관해유도)" = "neutropenic",
                    "스테로이드 숙주 (호중구 정상)"   = "steroid",
                    "이식 (타크로리무스 + 저용량 스테로이드)" = "transplant",
                    "면역정상 (대조)"                 = "normal")),
      sliderInput("anc_rec", "호중구 회복일 (day; 99 = 회복 없음)",
                  min = 3, max = 99, value = 10, step = 1),
      sliderInput("inoc", "폐포 도달 분생포자 (log10 CFUe)",
                  min = 2, max = 5, value = 3, step = 0.25),
      sliderInput("pred", "프레드니손 등가 (mg/day)", min = 0, max = 100, value = 0, step = 10),
      hr(),
      h4("항진균제 (antifungal)"),
      selectInput("drug", "치료",
                  c("없음" = "none", "보리코나졸" = "vrc", "이사부코나졸" = "isa",
                    "리포솜 암포테리신 B" = "amb", "아니둘라풍긴" = "ech",
                    "보리코나졸 + 아니둘라풍긴" = "vrc_ech",
                    "포사코나졸 (예방/구제)" = "pos")),
      sliderInput("start_d", "치료 시작일 (day)", min = 1, max = 14, value = 3, step = 1),
      sliderInput("vrc_dose", "보리코나졸 유지용량 (mg/kg q12h)",
                  min = 1, max = 9, value = 4, step = 0.2),
      sliderInput("amb_dose", "L-AmB 용량 (mg/kg/day)", min = 1, max = 10, value = 3, step = 1),
      checkboxInput("gcsf", "G-CSF 병용", FALSE),
      hr(),
      h4("균주 (isolate)"),
      selectInput("geno", "CYP2C19 유전형", names(CYP2C19), selected = "NM (*1/*1)"),
      sliderInput("mic", "보리코나졸 MIC (mg/L)", min = 0.125, max = 8, value = 0.5, step = 0.125),
      checkboxInput("tr34", "TR34/L98H 환경형 내성주", FALSE),
      hr(),
      sliderInput("horizon", "관찰 기간 (일)", min = 28, max = 120, value = 84, step = 7),
      actionButton("go", "시뮬레이션 실행", class = "btn-primary btn-block")
    ),
    mainPanel(
      width = 9,
      tabsetPanel(
        id = "tabs", type = "tabs",

        ## ------------------------------------------------------------ tab 1
        tabPanel(
          "1 · 환자 프로파일",
          br(),
          fluidRow(
            column(4, wellPanel(h4("숙주 상태"), tableOutput("tbl_host"))),
            column(4, wellPanel(h4("치료"),      tableOutput("tbl_rx"))),
            column(4, wellPanel(h4("결과 요약"), tableOutput("tbl_out")))),
          plotOutput("p_profile", height = "440px"),
          helpText("호중구(ANC)·병소 부피·관류·산소화를 한 화면에 둔 이유는, ",
                   "이 네 곡선이 사실상 서로의 원인이기 때문입니다.")),

        ## ------------------------------------------------------------ tab 2
        tabPanel(
          "2 · 약동학 (PK)",
          br(),
          plotOutput("p_pk", height = "330px"),
          plotOutput("p_site", height = "330px"),
          helpText("아래 그림의 점선(병소 내 농도 = ELF 농도 × 관류^γ)이 실선(ELF)에서 ",
                   "멀어지는 구간이 곧 혈관침습이 약을 차단하고 있는 구간입니다.")),

        ## ------------------------------------------------------------ tab 3
        tabPanel(
          "3 · CYP2C19 · TDM",
          br(),
          plotOutput("p_geno", height = "340px"),
          DTOutput("dt_geno"),
          helpText("같은 mg/kg 용량에서 유전형만으로 AUC가 5배 이상 벌어집니다. ",
                   "TDM은 이 분산을 좁히는 장치이지, 평균 노출을 올리는 장치가 아닙니다.")),

        ## ------------------------------------------------------------ tab 4
        tabPanel(
          "4 · 진균 부하 · 병소",
          br(),
          plotOutput("p_burden", height = "340px"),
          plotOutput("p_lesion", height = "300px"),
          helpText("굵은 수평 점선은 아졸 최대효과에서의 성장 하한 ",
                   "k_grow·(1 − Imax)에 해당하는 궤적입니다. ",
                   "호중구가 없으면 어떤 아졸 노출도 이 선 아래로 내려가지 못합니다.")),

        ## ------------------------------------------------------------ tab 5
        tabPanel(
          "5 · 바이오마커: GM · BDG · PCR",
          br(),
          plotOutput("p_bio", height = "330px"),
          plotOutput("p_gmdiss", height = "330px"),
          helpText("아래 그림이 이 모델의 임상적으로 가장 도발적인 예측입니다. ",
                   "1주차에 GM이 더 낮은 팔이 균 부하는 더 높을 수 있습니다 — ",
                   "GM은 재고(stock)가 아니라 유량(flux)의 지표이기 때문입니다.")),

        ## ------------------------------------------------------------ tab 6
        tabPanel(
          "6 · 임상 엔드포인트",
          br(),
          plotOutput("p_surv", height = "330px"),
          plotOutput("p_organ", height = "330px")),

        ## ------------------------------------------------------------ tab 7
        tabPanel(
          "7 · 시나리오 비교",
          br(),
          checkboxGroupInput(
            "cmp", "비교할 팔", inline = TRUE,
            choices = c("무치료" = "none", "보리코나졸" = "vrc", "이사부코나졸" = "isa",
                        "L-AmB 3" = "amb", "아니둘라풍긴" = "ech",
                        "VRC+아니둘라풍긴" = "vrc_ech"),
            selected = c("none", "vrc", "amb")),
          plotOutput("p_cmp", height = "380px"),
          DTOutput("dt_cmp")),

        ## ------------------------------------------------------------ tab 8
        tabPanel(
          "8 · 치료 지연 분기",
          br(),
          plotOutput("p_delay", height = "380px"),
          DTOutput("dt_delay"),
          helpText("관류 되먹임(γ)을 0으로 끄면 이 절벽이 사라집니다 — ",
                   "즉 늦은 치료가 실패하는 이유는 '균이 많아서'가 아니라 ",
                   "'약이 도달하지 못해서'라는 것이 모델의 주장입니다.")),

        ## ------------------------------------------------------------ tab 9
        tabPanel(
          "9 · 약물 지렛대 vs 숙주 지렛대",
          br(),
          plotOutput("p_lever", height = "420px"),
          DTOutput("dt_lever"),
          helpText("가로축을 고정하고 약을 바꾼 폭과, 약을 고정하고 골수를 바꾼 폭을 ",
                   "같은 단위(12주 사망률 %p)로 비교합니다.")),

        ## ----------------------------------------------------------- tab 10
        tabPanel(
          "10 · 내성 · MIC",
          br(),
          plotOutput("p_mic", height = "340px"),
          plotOutput("p_resist", height = "320px")),

        ## ----------------------------------------------------------- tab 11
        tabPanel(
          "11 · 독성 · 약물상호작용",
          br(),
          plotOutput("p_tox", height = "340px"),
          plotOutput("p_ddi", height = "320px"),
          helpText("타크로리무스 곡선은 아졸을 '항진균제'가 아니라 ",
                   "'CYP3A4 억제제'로 본 그림입니다."))
      )
    )
  ),
  hr(),
  tags$small("교육·연구용 QSP 모델입니다. 임상 의사결정에 사용하지 마십시오. ",
             "모든 수치는 ipa_reference_model.py 로 독립 재현·검증되었습니다.")
)

## -----------------------------------------------------------------------------
## SERVER
## -----------------------------------------------------------------------------
server <- function(input, output, session) {

  build <- reactive({
    ph  <- input$phenotype
    par <- list(
      INOCULUM   = 10^input$inoc,
      CYP2C19    = unname(CYP2C19[[input$geno]]),
      MIC_VRC    = input$mic,
      RESIST_FRAC = if (isTRUE(input$tr34)) 1 else 0,
      T_ANC_REC  = if (ph == "neutropenic") ifelse(input$anc_rec >= 99, 1e9, input$anc_rec * 24) else 0,
      ANC0       = switch(ph, neutropenic = 0.05, steroid = 4, transplant = 3, normal = 4))
    d <- switch(input$drug,
      none    = rx_none(),
      vrc     = rx_vrc(input$start_d, maint = input$vrc_dose),
      isa     = rx_isa(input$start_d),
      amb     = rx_amb(input$start_d, mgkg = input$amb_dose),
      ech     = rx_ech(input$start_d),
      vrc_ech = rbind(rx_vrc(input$start_d, maint = input$vrc_dose), rx_ech(input$start_d)),
      pos     = rx_pos(0, 90))
    if (input$pred > 0) d <- rbind(d, rx_steroid(input$pred, 60))
    if (ph == "steroid" && input$pred == 0)   d <- rbind(d, rx_steroid(60, 60))
    if (ph == "transplant") d <- rbind(d, rx_tac(), rx_steroid(20, 84))
    if (isTRUE(input$gcsf)) d <- rbind(d, rx_gcsf(input$start_d, 12))
    list(dose = d, par = par)
  })

  run <- eventReactive(input$go, {
    b <- build()
    sim(b$dose, b$par, end = input$horizon * 24, delta = 2)
  }, ignoreNULL = FALSE)

  day <- function(df) df$time / 24

  ## ---------------------------------------------------------------- tab 1
  output$tbl_host <- renderTable({
    o <- run()
    data.frame(항목 = c("표현형", "호중구 회복", "접종량 (log10)", "스테로이드"),
               값 = c(input$phenotype,
                      ifelse(input$anc_rec >= 99, "없음", paste0("day ", input$anc_rec)),
                      sprintf("%.2f", input$inoc),
                      paste0(input$pred, " mg/day")))
  })
  output$tbl_rx <- renderTable({
    data.frame(항목 = c("약제", "시작일", "VRC 용량", "CYP2C19", "MIC"),
               값 = c(input$drug, paste0("day ", input$start_d),
                      paste0(input$vrc_dose, " mg/kg q12h"), input$geno,
                      paste0(input$mic, " mg/L")))
  })
  output$tbl_out <- renderTable({
    o <- run()
    data.frame(항목 = c("최고 log10 부하", "최저 관류", "최대 GM", "12주 사망률(%)"),
               값 = sprintf("%.2f", c(max(o$LOGB), min(o$PERF), max(o$GMser),
                                      tail(o$MORT, 1))))
  })
  output$p_profile <- renderPlot({
    o <- run()
    d <- data.frame(day = day(o),
                    `ANC (10^9/L)` = o$ANC, `병소 부피 (mL)` = o$VLES,
                    `관류 지수` = o$PERF, `PaO2/FiO2` = o$PFR,
                    check.names = FALSE) %>%
      pivot_longer(-day)
    ggplot(d, aes(day, value)) + geom_line(linewidth = 0.9, colour = PAL[1]) +
      facet_wrap(~name, scales = "free_y") + THEME +
      labs(x = "일", y = NULL, title = "숙주 프로파일")
  })

  ## ---------------------------------------------------------------- tab 2
  output$p_pk <- renderPlot({
    o <- run()
    d <- data.frame(day = day(o), 보리코나졸 = o$CP_VRC, 이사부코나졸 = o$CP_ISA,
                    포사코나졸 = o$CP_POS, `L-AmB` = o$CP_AMB,
                    에키노칸딘 = o$CP_ECH, check.names = FALSE) %>%
      pivot_longer(-day) %>% filter(value > 1e-6)
    if (!nrow(d)) return(NULL)
    ggplot(d, aes(day, value, colour = name)) + geom_line(linewidth = 0.8) +
      scale_colour_manual(values = PAL) + THEME +
      labs(x = "일", y = "혈장 농도 (mg/L)", colour = NULL, title = "혈장 약동학")
  })
  output$p_site <- renderPlot({
    o <- run()
    d <- rbind(
      data.frame(day = day(o), conc = o$VRC_elf,          층 = "ELF (상피이장액)", drug = "VRC"),
      data.frame(day = day(o), conc = o$VRC_elf * o$PERF, 층 = "병소 내 (× 관류)",  drug = "VRC"),
      data.frame(day = day(o), conc = o$AMB_elf,          층 = "ELF (상피이장액)", drug = "L-AmB"),
      data.frame(day = day(o), conc = o$AMB_elf * o$PERF, 층 = "병소 내 (× 관류)",  drug = "L-AmB")) %>%
      filter(conc > 1e-6)
    if (!nrow(d)) return(NULL)
    ggplot(d, aes(day, conc, colour = drug, linetype = 층)) +
      geom_line(linewidth = 0.85) + scale_colour_manual(values = PAL) + THEME +
      labs(x = "일", y = "농도 (mg/L)", colour = NULL, linetype = NULL,
           title = "폐 대 병소: 관류가 만드는 간극")
  })

  ## ---------------------------------------------------------------- tab 3
  geno_tab <- reactive({
    bind_rows(lapply(names(CYP2C19), function(k) {
      o <- sim(rx_vrc(0, 14, maint = input$vrc_dose),
               list(CYP2C19 = unname(CYP2C19[[k]]), INOCULUM = 0),
               end = 14 * 24, delta = 0.25)
      w <- o[o$time >= 12 * 24 & o$time <= 13 * 24, ]
      data.frame(유전형 = k, Cmin = min(w$CP_VRC), Cmax = max(w$CP_VRC),
                 AUC24 = sum(diff(w$time) * (head(w$CP_VRC, -1) + tail(w$CP_VRC, -1)) / 2))
    }))
  })
  output$p_geno <- renderPlot({
    g <- geno_tab() %>% mutate(유전형 = factor(유전형, levels = names(CYP2C19)))
    ggplot(g, aes(유전형, AUC24, fill = 유전형)) +
      geom_col(width = 0.65) +
      geom_hline(yintercept = c(24, 132), linetype = 2, colour = "grey30") +
      annotate("text", x = 1, y = 138, label = "치료 범위 (trough 1–5.5 mg/L 상당)",
               hjust = 0, size = 3.4) +
      scale_fill_manual(values = PAL, guide = "none") + THEME +
      labs(x = NULL, y = "정상상태 AUC24 (mg·h/L)",
           title = sprintf("동일한 %.1f mg/kg q12h 용량에서의 유전형별 노출", input$vrc_dose))
  })
  output$dt_geno <- renderDT(datatable(geno_tab() %>% mutate(across(where(is.numeric), ~round(.x, 2))),
                                       options = list(dom = "t")))

  ## ---------------------------------------------------------------- tab 4
  output$p_burden <- renderPlot({
    o <- run()
    p <- as.list(mrgsolve::param(mod))
    fl <- p$KGROW * (1 - p$IMAX_AZOLE)
    t0 <- input$start_d * 24
    ref <- data.frame(day = day(o)) %>%
      mutate(y = o$LOGB[which.min(abs(o$time - t0))] +
               pmax(0, (time <- day * 24) - t0) * fl / log(10))
    ggplot() +
      geom_line(data = data.frame(day = day(o), y = o$LOGB), aes(day, y),
                linewidth = 1.0, colour = PAL[2]) +
      geom_line(data = ref, aes(day, y), linetype = 2, linewidth = 0.8, colour = "grey35") +
      annotate("text", x = max(day(o)) * 0.55, y = max(o$LOGB),
               label = "회색 점선 = 아졸 최대효과에서의 성장 하한", size = 3.6, hjust = 0) +
      THEME + labs(x = "일", y = "log10 균사 부하 (CFUe)", title = "진균 부하")
  })
  output$p_lesion <- renderPlot({
    o <- run()
    d <- data.frame(day = day(o), `병소 부피 (mL)` = o$VLES,
                    `경색 부피 (mL)` = o$NEC, `관류 지수` = o$PERF,
                    `혈관침습 지수` = o$ANG, check.names = FALSE) %>% pivot_longer(-day)
    ggplot(d, aes(day, value)) + geom_line(linewidth = 0.9, colour = PAL[4]) +
      facet_wrap(~name, scales = "free_y") + THEME + labs(x = "일", y = NULL)
  })

  ## ---------------------------------------------------------------- tab 5
  output$p_bio <- renderPlot({
    o <- run()
    d <- data.frame(day = day(o), `혈청 GM (ODI)` = o$GMser,
                    `BDG (pg/mL)` = o$BDG, `PCR (copies/mL)` = o$PCRs,
                    `log10 부하` = o$LOGB, check.names = FALSE) %>% pivot_longer(-day)
    ggplot(d, aes(day, value)) + geom_line(linewidth = 0.9, colour = PAL[3]) +
      facet_wrap(~name, scales = "free_y") + THEME + labs(x = "일", y = NULL)
  })
  gm_tab <- reactive({
    arms <- list(`무치료` = rx_none(), `보리코나졸` = rx_vrc(input$start_d),
                 `L-AmB 3 mg/kg` = rx_amb(input$start_d, mgkg = 3),
                 `아니둘라풍긴` = rx_ech(input$start_d))
    bind_rows(lapply(names(arms), function(a) {
      o <- sim(arms[[a]], list(T_ANC_REC = 336), end = 20 * 24, delta = 1)
      data.frame(day = o$time / 24, arm = a, GM = o$GMser, logB = o$LOGB)
    }))
  })
  output$p_gmdiss <- renderPlot({
    d <- gm_tab() %>% pivot_longer(c(GM, logB))
    d$name <- factor(d$name, c("GM", "logB"), c("혈청 GM 지수", "log10 균 부하"))
    ggplot(d, aes(day, value, colour = arm)) + geom_line(linewidth = 0.9) +
      facet_wrap(~name, scales = "free_y") +
      scale_colour_manual(values = PAL) + THEME +
      labs(x = "일", y = NULL, colour = NULL,
           title = "GM은 재고가 아니라 유량을 본다 — 1주차 순위가 뒤집힌다")
  })

  ## ---------------------------------------------------------------- tab 6
  output$p_surv <- renderPlot({
    o <- run()
    ggplot(data.frame(day = day(o), s = 100 * o$SURV), aes(day, s)) +
      geom_line(linewidth = 1.0, colour = PAL[5]) + ylim(0, 100) + THEME +
      labs(x = "일", y = "생존 확률 (%)", title = "모델 생존 곡선")
  })
  output$p_organ <- renderPlot({
    o <- run()
    d <- data.frame(day = day(o), `체온 (°C)` = o$TEMP, `CRP (mg/L)` = o$CRP,
                    `IL-6 (pg/mL)` = o$IL6, `PaO2/FiO2` = o$PFR,
                    check.names = FALSE) %>% pivot_longer(-day)
    ggplot(d, aes(day, value)) + geom_line(linewidth = 0.9, colour = PAL[6]) +
      facet_wrap(~name, scales = "free_y") + THEME + labs(x = "일", y = NULL)
  })

  ## ---------------------------------------------------------------- tab 7
  cmp_run <- reactive({
    req(length(input$cmp) > 0)
    b <- build()
    arms <- list(none = rx_none(), vrc = rx_vrc(input$start_d, maint = input$vrc_dose),
                 isa = rx_isa(input$start_d), amb = rx_amb(input$start_d, mgkg = input$amb_dose),
                 ech = rx_ech(input$start_d),
                 vrc_ech = rbind(rx_vrc(input$start_d, maint = input$vrc_dose),
                                 rx_ech(input$start_d)))
    bind_rows(lapply(input$cmp, function(k) {
      sim(arms[[k]], b$par, end = input$horizon * 24, delta = 4) %>%
        mutate(arm = k)
    }))
  })
  output$p_cmp <- renderPlot({
    d <- cmp_run() %>%
      transmute(day = time / 24, arm,
                `log10 부하` = LOGB, `혈청 GM` = GMser,
                `관류` = PERF, `사망률 (%)` = MORT) %>%
      pivot_longer(-c(day, arm))
    ggplot(d, aes(day, value, colour = arm)) + geom_line(linewidth = 0.9) +
      facet_wrap(~name, scales = "free_y") +
      scale_colour_manual(values = PAL) + THEME +
      labs(x = "일", y = NULL, colour = NULL)
  })
  output$dt_cmp <- renderDT({
    d <- cmp_run() %>% group_by(arm) %>%
      summarise(`최고 log10 부하` = round(max(LOGB), 2),
                `d42 log10 부하`  = round(LOGB[which.min(abs(time - 42 * 24))], 2),
                `최저 관류`       = round(min(PERF), 3),
                `최대 GM`         = round(max(GMser), 2),
                `12주 사망률 (%)` = round(tail(MORT, 1), 1), .groups = "drop")
    datatable(d, options = list(dom = "t"))
  })

  ## ---------------------------------------------------------------- tab 8
  delay_tab <- reactive({
    bind_rows(lapply(c(1:8, 10, 12, 14), function(dd) {
      a <- sim(rx_vrc(dd, maint = input$vrc_dose), list(T_ANC_REC = 336), end = 84 * 24, delta = 6)
      b <- sim(rx_vrc(dd, maint = input$vrc_dose),
               list(T_ANC_REC = 336, PERF_GAMMA = 0), end = 84 * 24, delta = 6)
      data.frame(`시작일` = dd,
                 `d84 log10 부하 (γ=1)` = round(tail(a$LOGB, 1), 2),
                 `d84 log10 부하 (γ=0)` = round(tail(b$LOGB, 1), 2),
                 `최저 관류` = round(min(a$PERF), 3),
                 `12주 사망률 (%)` = round(tail(a$MORT, 1), 1),
                 check.names = FALSE)
    }))
  })
  output$p_delay <- renderPlot({
    d <- delay_tab()
    dd <- data.frame(day = rep(d$`시작일`, 2),
                     y = c(d$`d84 log10 부하 (γ=1)`, d$`d84 log10 부하 (γ=0)`),
                     model = rep(c("관류 되먹임 있음 (γ=1)", "없음 (γ=0)"), each = nrow(d)))
    ggplot(dd, aes(day, y, colour = model)) +
      geom_line(linewidth = 1.0) + geom_point(size = 2) +
      scale_colour_manual(values = PAL[c(2, 1)]) + THEME +
      labs(x = "치료 시작일", y = "12주째 log10 균 부하", colour = NULL,
           title = "치료 지연 절벽은 관류 되먹임이 만든다")
  })
  output$dt_delay <- renderDT(datatable(delay_tab(), options = list(dom = "t")))

  ## ---------------------------------------------------------------- tab 9
  lever_tab <- reactive({
    arms <- list(무치료 = rx_none(), VRC = rx_vrc(input$start_d),
                 ISA = rx_isa(input$start_d), `L-AmB` = rx_amb(input$start_d, mgkg = 3))
    bind_rows(lapply(c(5, 10, 14, 21, 28, 99), function(r) {
      bind_rows(lapply(names(arms), function(a) {
        o <- sim(arms[[a]], list(T_ANC_REC = ifelse(r >= 99, 1e9, r * 24)),
                 end = 84 * 24, delta = 8)
        data.frame(`호중구 회복일` = ifelse(r >= 99, 99, r), 팔 = a,
                   `12주 사망률` = round(tail(o$MORT, 1), 1), check.names = FALSE)
      }))
    }))
  })
  output$p_lever <- renderPlot({
    d <- lever_tab()
    ggplot(d, aes(factor(`호중구 회복일`), `12주 사망률`, fill = 팔)) +
      geom_col(position = position_dodge(0.8), width = 0.75) +
      scale_fill_manual(values = PAL) + THEME +
      labs(x = "호중구 회복일 (99 = 회복 없음)", y = "12주 사망률 (%)", fill = NULL,
           title = "막대 사이(약)보다 막대 그룹 사이(골수)가 더 멀다")
  })
  output$dt_lever <- renderDT({
    d <- lever_tab() %>% pivot_wider(names_from = 팔, values_from = `12주 사망률`)
    datatable(d, options = list(dom = "t"))
  })

  ## --------------------------------------------------------------- tab 10
  output$p_mic <- renderPlot({
    d <- bind_rows(lapply(c(0.125, 0.25, 0.5, 1, 2, 4, 8), function(m) {
      o <- sim(rx_vrc(input$start_d, maint = input$vrc_dose),
               list(T_ANC_REC = 336, MIC_VRC = m), end = 84 * 24, delta = 8)
      data.frame(MIC = m, logB_d42 = o$LOGB[which.min(abs(o$time - 42 * 24))],
                 mort = tail(o$MORT, 1))
    }))
    ggplot(d, aes(MIC, mort)) + geom_line(linewidth = 1.0, colour = PAL[2]) +
      geom_point(size = 2.2) + scale_x_log10(breaks = c(0.125, 0.5, 2, 8)) +
      geom_vline(xintercept = 1, linetype = 2) +
      annotate("text", x = 1.1, y = min(d$mort), label = "EUCAST ECOFF", hjust = 0, size = 3.5) +
      THEME + labs(x = "보리코나졸 MIC (mg/L, log 축)", y = "12주 사망률 (%)",
                   title = "MIC 사다리")
  })
  output$p_resist <- renderPlot({
    arms <- list(무치료 = rx_none(), `VRC 단독` = rx_vrc(input$start_d),
                 `L-AmB` = rx_amb(input$start_d, mgkg = 3),
                 `VRC+에키노칸딘` = rbind(rx_vrc(input$start_d), rx_ech(input$start_d)))
    d <- bind_rows(lapply(names(arms), function(a) {
      o <- sim(arms[[a]], list(T_ANC_REC = 336), end = 84 * 24, delta = 8)
      data.frame(day = o$time / 24, arm = a, rfrac = pmax(o$RFRAC, 1e-12))
    }))
    ggplot(d, aes(day, rfrac, colour = arm)) + geom_line(linewidth = 0.9) +
      scale_y_log10() + scale_colour_manual(values = PAL) + THEME +
      labs(x = "일", y = "내성 아집단 비율 (log 축)", colour = NULL,
           title = "아졸 단독요법 하의 내성 선택")
  })

  ## --------------------------------------------------------------- tab 11
  output$p_tox <- renderPlot({
    o <- run()
    d <- data.frame(day = day(o), `ALT (U/L)` = o$ALT, `크레아티닌 (mg/dL)` = o$SCR,
                    `QTc 증가 (ms)` = o$QTC, `신경/시각 효과` = o$NEUROE,
                    check.names = FALSE) %>% pivot_longer(-day)
    ggplot(d, aes(day, value)) + geom_line(linewidth = 0.9, colour = PAL[7]) +
      facet_wrap(~name, scales = "free_y") + THEME + labs(x = "일", y = NULL)
  })
  output$p_ddi <- renderPlot({
    arms <- list(`타크로리무스 단독` = rx_tac(),
                 `+ 보리코나졸` = rbind(rx_tac(), rx_vrc(0)),
                 `+ 이사부코나졸` = rbind(rx_tac(), rx_isa(0)),
                 `+ 포사코나졸` = rbind(rx_tac(), rx_pos(0, 30)))
    d <- bind_rows(lapply(names(arms), function(a) {
      o <- sim(arms[[a]], list(T_ANC_REC = 0, ANC0 = 3), end = 20 * 24, delta = 1)
      data.frame(day = o$time / 24, arm = a, tac = o$CP_TAC)
    }))
    ggplot(d, aes(day, tac, colour = arm)) + geom_line(linewidth = 0.9) +
      geom_hline(yintercept = c(5, 15), linetype = 2, colour = "grey40") +
      scale_colour_manual(values = PAL) + THEME +
      labs(x = "일", y = "타크로리무스 (ng/mL)", colour = NULL,
           title = "회색 띠 = 통상 목표 농도 5–15 ng/mL")
  })
}

shinyApp(ui, server)
