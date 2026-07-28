## =============================================================================
##  Familial Chylomicronemia Syndrome (FCS) — QSP Shiny dashboard
##  fcs_shiny_app.R
## -----------------------------------------------------------------------------
##  Run with:
##      setwd("familial-chylomicronemia-syndrome")
##      shiny::runApp("fcs_shiny_app.R")
##
##  Requires: shiny, mrgsolve, ggplot2, DT  (dplyr/tidyr optional)
##  The app sources fcs_mrgsolve_model.R, so keep both files together.
##
##  The dashboard is organised around the model's central claim: TRL clearance
##  is a SUM of two limbs, the genotype multiplies the first by zero, and every
##  therapeutic decision in FCS is therefore a decision about the second limb
##  or about the input.  Each tab makes one part of that argument visible.
## =============================================================================

library(shiny)
library(mrgsolve)
library(ggplot2)

if(!exists("mod")) source("fcs_mrgsolve_model.R")

THEME <- theme_minimal(base_size = 12) +
  theme(panel.grid.minor = element_blank(),
        plot.title = element_text(face = "bold", size = 13),
        legend.position = "bottom")

PAL <- c(placebo = "#90a4ae", diet = "#7cb342", fibrate = "#546e7a",
         evinacumab = "#6d4c41", volanesorsen = "#00897b",
         olezarsen = "#43a047", plozasiran = "#0288d1")

## =============================================================================
##  UI
## =============================================================================
ui <- fluidPage(
  titlePanel(
    div(style = "line-height:1.25",
        strong("가족성 유미지질혈증 증후군 (Familial Chylomicronemia Syndrome, FCS) — QSP 대시보드"),
        div(style = "font-size:13px;color:#555;font-weight:normal",
            "포화된 청소 용량(saturated clearance capacity)으로서의 FCS · ",
            "CL_total = CL_LPL(genotype ≈ 0) + CL_independent(apoC-III, saturable)"))),

  sidebarLayout(
    sidebarPanel(
      width = 3,

      h4("① 환자 (Patient)"),
      selectInput("geno", "유전형 / Genotype",
                  choices = c("LPL-null FCS"            = "fcs_null",
                              "APOC2 결핍 FCS"          = "fcs_apoc2",
                              "GPIHBP1 결핍 FCS"        = "fcs_gpihbp1",
                              "잔존 LPL 5% (부분결핍)"  = "partial_05",
                              "잔존 LPL 20% (경증)"     = "partial_20",
                              "다인성 유미지질혈증 MCS" = "mcs",
                              "정상 대조 (healthy)"     = "healthy"),
                  selected = "fcs_null"),
      sliderInput("fat", "식이 장쇄지방 LCT (g/day)", 5, 100, 20, step = 5),
      sliderInput("fmct", "MCT 대체 비율 (fraction)", 0, 0.8, 0, step = 0.1),
      checkboxInput("alc",  "알코올 (alcohol)", FALSE),
      checkboxInput("est",  "경구 에스트로겐 (oral oestrogen)", FALSE),
      checkboxInput("preg", "임신 (pregnancy)", FALSE),

      hr(),
      h4("② 치료 (Therapy)"),
      selectInput("drug", "약제 / Drug",
                  choices = c("없음 (diet only)"        = "none",
                              "페노피브레이트 + ω-3"     = "fibrate",
                              "에비나쿠맙 15 mg/kg q4w"  = "evinacumab",
                              "볼라네소르센 300 mg qw"   = "volanesorsen",
                              "올레자르센 80 mg q4w"     = "olezarsen",
                              "플로자시란 25 mg q12w"    = "plozasiran"),
                  selected = "plozasiran"),
      conditionalPanel("input.drug == 'olezarsen'",
                       sliderInput("ole_mg", "올레자르센 용량 (mg)", 50, 80, 80, 10)),
      conditionalPanel("input.drug == 'plozasiran'",
                       sliderInput("plo_mg", "플로자시란 용량 (mg)", 10, 50, 25, 5)),
      sliderInput("days", "관찰기간 (days)", 90, 730, 365, step = 30),
      checkboxInput("binge", "30일마다 60 g 지방 폭식 (non-adherence)", FALSE),

      hr(),
      h4("③ 모델 파라미터"),
      sliderInput("vmax", "Vmax(limb 2, LPL-비의존) mg/h", 600, 2600, 1300, 100),
      sliderInput("imax", "apoC-III 최대 증폭 IMAX_C3", 0.5, 4, 2.2, 0.1),
      sliderInput("hill", "AP 위험 볼록도 HILL_AP", 1.0, 3.0, 1.7, 0.1),
      actionButton("go", "시뮬레이션 실행", class = "btn-primary btn-block"),

      hr(),
      div(style = "font-size:11px;color:#777",
          "본 모델은 교육·연구용입니다. 임상 의사결정에 사용할 수 없습니다.",
          br(), "Research and education only — not for clinical use.")
    ),

    mainPanel(
      width = 9,
      tabsetPanel(
        id = "tabs", type = "tabs",

        ## ---------------------------------------------------------------- 1
        tabPanel("① 환자 프로파일",
          br(),
          fluidRow(
            column(4, wellPanel(h4(textOutput("kpi_tg")),  p("평균 TG (mg/dL)"))),
            column(4, wellPanel(h4(textOutput("kpi_c3")),  p("혈장 apoC-III (mg/dL)"))),
            column(4, wellPanel(h4(textOutput("kpi_ap")),  p("췌장염 위험 (events/py)")))),
          fluidRow(
            column(4, wellPanel(h4(textOutput("kpi_tat")), p("TG>880 초과 일수 / 년"))),
            column(4, wellPanel(h4(textOutput("kpi_plt")), p("혈소판 최저치 (10⁹/L)"))),
            column(4, wellPanel(h4(textOutput("kpi_lpl")), p("LPL 청소율 (dL/h)")))),
          hr(),
          h4("환자 요약 / Patient summary"),
          verbatimTextOutput("profile_txt")),

        ## ---------------------------------------------------------------- 2
        tabPanel("② 지질 시계열 (PD)",
          br(), plotOutput("p_tg", height = "330px"),
          plotOutput("p_species", height = "300px"),
          helpText("점선 = 880 mg/dL (10 mmol/L) 급성췌장염 임계치. ",
                   "FCS에서는 진정한 공복 상태가 존재하지 않으며 식후 변동이 ",
                   "위험의 상당 부분을 차지한다.")),

        ## ---------------------------------------------------------------- 3
        tabPanel("③ 포화 절벽 (saturation cliff)",
          br(), plotOutput("p_sat", height = "420px"),
          DT::dataTableOutput("t_sat"),
          helpText("정상에서는 지방 섭취량과 TG가 거의 선형이지만 ",
                   "LPL-null FCS에서는 쌍곡선이 된다. '하루 20 g 미만'이 ",
                   "기울기가 아니라 절벽 가장자리인 이유.")),

        ## ---------------------------------------------------------------- 4
        tabPanel("④ 두 경로 분해 (limb decomposition)",
          br(), plotOutput("p_limb", height = "330px"),
          DT::dataTableOutput("t_limb"),
          helpText("limb1(LPL)은 모든 LPL-null 시나리오에서 정확히 0 mg/h이며 ",
                   "어떤 약도 이를 움직이지 못한다. 0 × 무엇 = 0.")),

        ## ---------------------------------------------------------------- 5
        tabPanel("⑤ 약물 PK",
          br(), plotOutput("p_pk", height = "330px"),
          plotOutput("p_target", height = "300px"),
          helpText("GalNAc 접합은 표적을 바꾸지 않고 전신 노출을 20-30배 ",
                   "줄였다. 볼라네소르센의 혈소판 신호는 약리학이 아니라 ",
                   "PK 문제였다.")),

        ## ---------------------------------------------------------------- 6
        tabPanel("⑥ 임상 엔드포인트",
          br(),
          fluidRow(column(6, plotOutput("p_haz",  height = "300px")),
                   column(6, plotOutput("p_prob", height = "300px"))),
          plotOutput("p_tat", height = "260px"),
          helpText("위험함수는 TG에 대해 볼록(Hill 1.7)하므로, 평균 TG의 ",
                   "% 변화보다 '임계치 초과 시간'의 감소가 실제 사건 감소를 ",
                   "더 잘 설명한다.")),

        ## ---------------------------------------------------------------- 7
        tabPanel("⑦ Jensen 격차 (식후 변동)",
          br(), plotOutput("p_jensen", height = "330px"),
          DT::dataTableOutput("t_jensen"),
          helpText("E[λ(TG)] > λ(E[TG]). 공복 TG 한 점으로 위험을 평가하면 ",
                   "체계적으로 과소평가한다.")),

        ## ---------------------------------------------------------------- 8
        tabPanel("⑧ 안전성",
          br(),
          fluidRow(column(6, plotOutput("p_plt", height = "300px")),
                   column(6, plotOutput("p_alt", height = "300px"))),
          DT::dataTableOutput("t_safety"),
          helpText("혈소판 감소는 전신 PS-ASO 조직 노출에 비례한다. ",
                   "GalNAc 접합체에서는 사실상 사라진다.")),

        ## ---------------------------------------------------------------- 9
        tabPanel("⑨ 시나리오 비교",
          br(), plotOutput("p_cmp", height = "360px"),
          DT::dataTableOutput("t_ledger"),
          helpText("dTG(%)와 AP 사건 감소(%)를 나란히 비교하라. ",
                   "후자가 항상 더 크다 — 볼록성의 배당금.")),

        ## ---------------------------------------------------------------- 10
        tabPanel("⑩ 바이오마커 · 보정",
          br(), plotOutput("p_bio", height = "330px"),
          h4("문헌 보정 표 / Calibration targets"),
          DT::dataTableOutput("t_calib"),
          helpText("각 목표치의 출처는 fcs_references.md 참조."))
      )
    )
  )
)

## =============================================================================
##  SERVER
## =============================================================================
server <- function(input, output, session) {

  pars <- reactive({
    p <- geno[[input$geno]]
    p$FMCT     <- input$fmct
    p$ALC      <- as.numeric(input$alc || isTRUE(p$ALC == 1))
    p$EST      <- as.numeric(input$est)
    p$PREG     <- as.numeric(input$preg)
    p$VMAX_IND <- input$vmax
    p$IMAX_C3  <- input$imax
    p$HILL_AP  <- input$hill
    if(input$drug == "fibrate") p$OM3_ON <- 1
    p
  })

  drug_events <- reactive({
    d <- input$days
    switch(input$drug,
      none         = NULL,
      fibrate      = fenofibrate(0, d),
      evinacumab   = evinacumab(0, ceiling(d/28)),
      volanesorsen = volanesorsen(0, ceiling(d/7)),
      olezarsen    = olezarsen(0, ceiling(d/28), mg = input$ole_mg),
      plozasiran   = plozasiran(0, ceiling(d/90), mg = input$plo_mg))
  })

  sim <- eventReactive(input$go, {
    e <- drug_events()
    if(input$binge) {
      b <- do.call(c, lapply(seq(30, input$days, by = 30),
                             function(k) fat_binge(60, k)))
      e <- if(is.null(e)) b else e + b
    }
    withProgress(message = "시뮬레이션 중...", value = 0.5,
      run_arm(pars(), events = e, days = input$days, fat = input$fat, delta = 2))
  }, ignoreNULL = FALSE)

  on_treat <- reactive({ d <- sim(); subset(d, day >= 0) })

  ## ------------------------------------------------------------------ KPIs
  output$kpi_tg  <- renderText(sprintf("%.0f", mean(on_treat()$TG)))
  output$kpi_c3  <- renderText(sprintf("%.1f", mean(on_treat()$APOC3)))
  output$kpi_ap  <- renderText({
    d <- on_treat(); sprintf("%.3f", (max(d$CUMHAZ)-min(d$CUMHAZ))/(max(d$day)/365))})
  output$kpi_tat <- renderText({
    d <- on_treat()
    sprintf("%.0f", (max(d$TAT880)-min(d$TAT880))/24/(max(d$day)/365))})
  output$kpi_plt <- renderText(sprintf("%.0f", min(on_treat()$PLT)))
  output$kpi_lpl <- renderText(sprintf("%.2f", mean(on_treat()$CLLPL)))

  output$profile_txt <- renderPrint({
    d <- on_treat(); b <- sim()
    base <- mean(b$TG[b$day > -30 & b$day <= 0])
    cat("유전형              :", input$geno, "\n")
    cat("식이 LCT            :", input$fat, "g/day  (MCT 대체", input$fmct, ")\n")
    cat("치료                :", input$drug, "\n")
    cat("치료 전 평균 TG     :", sprintf("%.0f mg/dL (%.1f mmol/L)",
                                          base, base/88.57), "\n")
    cat("치료 후 평균 TG     :", sprintf("%.0f mg/dL  (%+.1f %%)",
                                          mean(d$TG), 100*(mean(d$TG)/base-1)), "\n")
    cat("유미지질 분율       :", sprintf("%.0f %%", 100*mean(d$CM_C)/mean(d$TG)), "\n")
    cat("혈장 점도 지수      :", sprintf("%.2f", mean(d$VISCI)), "\n")
    cat("유발성 황색종 점수  :", sprintf("%.2f / 10", tail(d$XANTH, 1)), "\n")
    cat("간비종대 지수       :", sprintf("%.2f / 10", tail(d$HSM, 1)), "\n")
    cat("인지 피로 점수      :", sprintf("%.2f / 10", tail(d$FOG, 1)), "\n")
    cat("1년 내 췌장염 확률  :", sprintf("%.1f %%", 100*max(d$PROB_AP)), "\n")
  })

  ## ------------------------------------------------------------------ tab 2
  output$p_tg <- renderPlot({
    d <- sim()
    ggplot(d, aes(day, TG)) +
      geom_hline(yintercept = 880, linetype = 2, colour = "grey40") +
      geom_vline(xintercept = 0, linetype = 3, colour = "grey60") +
      geom_line(colour = "#d84315", linewidth = 0.6) +
      labs(title = "총 혈장 중성지방 (total plasma triglyceride)",
           x = "치료 시작 후 일수 (day)", y = "TG (mg/dL)") + THEME
  })

  output$p_species <- renderPlot({
    d <- sim()
    long <- rbind(
      data.frame(day = d$day, conc = d$CM_C,   species = "chylomicron-TG"),
      data.frame(day = d$day, conc = d$VLDL_C, species = "VLDL-TG"),
      data.frame(day = d$day, conc = d$REM_C,  species = "remnant-TG"))
    ggplot(long, aes(day, conc, colour = species)) +
      geom_line(linewidth = 0.6) +
      scale_colour_manual(values = c("chylomicron-TG" = "#0277bd",
                                     "VLDL-TG" = "#ef6c00",
                                     "remnant-TG" = "#7e57c2")) +
      labs(title = "지단백 종별 분해 (FCS에서는 유미지질이 지배적)",
           x = "day", y = "mg/dL", colour = NULL) + THEME
  })

  ## ------------------------------------------------------------------ tab 3
  sat_tab <- reactive({
    fats <- c(5,10,15,20,25,30,40,50,60,80,100)
    do.call(rbind, lapply(c("fcs_null","partial_05","partial_20","healthy"),
      function(g) do.call(rbind, lapply(fats, function(f) {
        d  <- run_arm(geno[[g]], days = 30, burn = 120, fat = f, delta = 4)
        ss <- subset(d, day > 18)
        data.frame(genotype = g, fat = f, TG = mean(ss$TG),
                   peak = max(ss$TG),
                   pct_above = 100*mean(ss$TG > 880))
      }))))
  })

  output$p_sat <- renderPlot({
    ggplot(sat_tab(), aes(fat, TG, colour = genotype)) +
      geom_hline(yintercept = 880, linetype = 2, colour = "grey40") +
      geom_line(linewidth = 0.9) + geom_point(size = 2) +
      scale_colour_manual(values = c(fcs_null = "#d84315", partial_05 = "#f9a825",
                                     partial_20 = "#43a047", healthy = "#1565c0")) +
      labs(title = "포화 절벽: 식이 지방 대비 정상상태 TG",
           x = "식이 장쇄지방 (g/day)", y = "정상상태 TG (mg/dL)",
           colour = NULL) + THEME
  })
  output$t_sat <- DT::renderDataTable(
    DT::datatable(sat_tab(), options = list(pageLength = 8), rownames = FALSE) |>
      DT::formatRound(c("TG","peak","pct_above"), 1))

  ## ------------------------------------------------------------------ tab 4
  output$p_limb <- renderPlot({
    d <- sim()
    long <- rbind(
      data.frame(day = d$day, flux = d$FLUX1, limb = "limb 1 · LPL"),
      data.frame(day = d$day, flux = d$FLUX2, limb = "limb 2 · LPL-independent"),
      data.frame(day = d$day, flux = d$FLUX3, limb = "residual scavenging"))
    ggplot(long, aes(day, flux, colour = limb)) +
      geom_line(linewidth = 0.7) +
      scale_colour_manual(values = c("limb 1 · LPL" = "#d81b60",
                                     "limb 2 · LPL-independent" = "#7e57c2",
                                     "residual scavenging" = "#90a4ae")) +
      labs(title = "청소 경로별 유량 (mg/h)", x = "day", y = "flux (mg/h)",
           colour = NULL) + THEME
  })
  output$t_limb <- DT::renderDataTable({
    d <- on_treat()
    tb <- data.frame(
      limb = c("limb 1 (LPL)", "limb 2 (LPL-independent)", "residual"),
      mean_flux_mg_h = c(mean(d$FLUX1), mean(d$FLUX2), mean(d$FLUX3)))
    tb$percent <- 100*tb$mean_flux_mg_h/sum(tb$mean_flux_mg_h)
    DT::datatable(tb, rownames = FALSE, options = list(dom = "t")) |>
      DT::formatRound(c("mean_flux_mg_h","percent"), 2)
  })

  ## ------------------------------------------------------------------ tab 5
  output$p_pk <- renderPlot({
    d <- sim()
    long <- rbind(
      data.frame(day = d$day, amt = d$VOL_LIV,  cmt = "volanesorsen · liver"),
      data.frame(day = d$day, amt = d$VOL_SYS,  cmt = "volanesorsen · systemic"),
      data.frame(day = d$day, amt = d$OLE_LIV,  cmt = "olezarsen · liver"),
      data.frame(day = d$day, amt = d$OLE_SYS,  cmt = "olezarsen · systemic"),
      data.frame(day = d$day, amt = d$PLO_RISC, cmt = "plozasiran · RISC"),
      data.frame(day = d$day, amt = d$CEVI2*10, cmt = "evinacumab · plasma x10"))
    long <- long[long$amt > 1e-6, ]
    if(!nrow(long)) return(ggplot() + labs(title = "약물 없음 (no drug)") + THEME)
    ggplot(long, aes(day, amt, colour = cmt)) + geom_line(linewidth = 0.7) +
      labs(title = "약물 조직/효과 구획 (mg 또는 AU)", x = "day", y = "amount",
           colour = NULL) + THEME
  })
  output$p_target <- renderPlot({
    d <- sim()
    ggplot(d, aes(day, APOC3)) + geom_line(colour = "#00897b", linewidth = 0.8) +
      labs(title = "표적 관여: 혈장 apoC-III", x = "day",
           y = "apoC-III (mg/dL)") + THEME
  })

  ## ------------------------------------------------------------------ tab 6
  output$p_haz <- renderPlot({
    d <- sim()
    ggplot(d, aes(day, HAZ_YR)) + geom_line(colour = "#bf360c", linewidth = 0.7) +
      labs(title = "급성 췌장염 순간 위험", x = "day",
           y = "events / patient-year") + THEME
  })
  output$p_prob <- renderPlot({
    d <- on_treat()
    ggplot(d, aes(day, 100*PROB_AP)) +
      geom_line(colour = "#5e35b1", linewidth = 0.8) +
      labs(title = "누적 췌장염 확률", x = "day", y = "%") + THEME
  })
  output$p_tat <- renderPlot({
    d <- on_treat()
    ggplot(d, aes(day, TAT880/24)) +
      geom_line(colour = "#4527a0", linewidth = 0.8) +
      labs(title = "TG > 880 mg/dL 누적 초과 일수 (time above threshold)",
           x = "day", y = "누적 일수") + THEME
  })

  ## ------------------------------------------------------------------ tab 7
  jen <- eventReactive(input$go, { FCS_jensen_gap() })
  output$t_jensen <- DT::renderDataTable(
    DT::datatable(jen(), rownames = FALSE, options = list(dom = "t")) |>
      DT::formatRound(2:9, 3))
  output$p_jensen <- renderPlot({
    j <- jen()
    long <- rbind(
      data.frame(fat = j$fat_g_day, h = j$haz_of_fasting, k = "λ(공복 TG)"),
      data.frame(fat = j$fat_g_day, h = j$haz_of_mean,    k = "λ(평균 TG)"),
      data.frame(fat = j$fat_g_day, h = j$mean_of_haz,    k = "평균 λ(TG) — 진짜 위험"))
    ggplot(long, aes(factor(fat), h, fill = k)) +
      geom_col(position = "dodge") +
      scale_fill_manual(values = c("#90a4ae", "#f9a825", "#bf360c")) +
      labs(title = "Jensen 격차: 볼록한 위험함수 + 식후 변동",
           x = "식이 지방 (g/day)", y = "events / patient-year", fill = NULL) +
      THEME
  })

  ## ------------------------------------------------------------------ tab 8
  output$p_plt <- renderPlot({
    d <- sim()
    ggplot(d, aes(day, PLT)) +
      geom_hline(yintercept = 100, linetype = 2, colour = "#d84315") +
      geom_hline(yintercept = 50,  linetype = 3, colour = "#b71c1c") +
      geom_line(colour = "#5e35b1", linewidth = 0.8) +
      labs(title = "혈소판 수", x = "day", y = "10⁹/L") + THEME
  })
  output$p_alt <- renderPlot({
    d <- sim()
    ggplot(d, aes(day, ALT)) + geom_line(colour = "#6d4c41", linewidth = 0.8) +
      labs(title = "혈청 ALT", x = "day", y = "U/L") + THEME
  })
  output$t_safety <- DT::renderDataTable({
    d <- on_treat()
    tb <- data.frame(
      metric = c("혈소판 최저치 (10⁹/L)", "혈소판 <100 인 시간 (%)",
                 "ALT 최고치 (U/L)", "전신 ASO 조직량 (mg)",
                 "간 ASO 조직량 (mg)"),
      value = c(min(d$PLT), 100*mean(d$PLT < 100), max(d$ALT),
                max(d$VOL_SYS + d$OLE_SYS), max(d$VOL_LIV + d$OLE_LIV)))
    DT::datatable(tb, rownames = FALSE, options = list(dom = "t")) |>
      DT::formatRound("value", 2)
  })

  ## ------------------------------------------------------------------ tab 9
  ledger <- eventReactive(input$go, { FCS_trial_ledger(days = input$days) })
  output$t_ledger <- DT::renderDataTable(
    DT::datatable(ledger(), rownames = FALSE,
                  options = list(pageLength = 8, scrollX = TRUE)) |>
      DT::formatRound(2:12, 2))
  output$p_cmp <- renderPlot({
    L <- ledger()
    long <- rbind(
      data.frame(arm = L$arm, v = -L$dTG_m12_pct,     k = "TG 감소 (%)"),
      data.frame(arm = L$arm, v = L$AP_reduction_pct, k = "췌장염 사건 감소 (%)"))
    ggplot(long, aes(reorder(arm, v), v, fill = k)) +
      geom_col(position = "dodge") + coord_flip() +
      scale_fill_manual(values = c("#0288d1", "#bf360c")) +
      labs(title = "% TG 변화 vs 실제 사건 감소 — 볼록성의 배당금",
           x = NULL, y = "%", fill = NULL) + THEME
  })

  ## ------------------------------------------------------------------ tab 10
  output$p_bio <- renderPlot({
    d <- sim()
    long <- rbind(
      data.frame(day = d$day, v = d$XANTH, k = "유발성 황색종"),
      data.frame(day = d$day, v = d$HSM,   k = "간비종대 지수"),
      data.frame(day = d$day, v = d$FOG,   k = "인지 피로"),
      data.frame(day = d$day, v = d$VISCI, k = "혈장 점도 지수"))
    ggplot(long, aes(day, v, colour = k)) + geom_line(linewidth = 0.7) +
      facet_wrap(~k, scales = "free_y") +
      labs(title = "비췌장 질환 부담 바이오마커", x = "day", y = NULL) +
      THEME + theme(legend.position = "none")
  })

  output$t_calib <- DT::renderDataTable({
    tb <- data.frame(
      target = c("정상 공복 TG",
                 "FCS · 10 g 지방/일",
                 "FCS · 20 g 지방/일",
                 "FCS · 60 g 지방/일",
                 "잔존 LPL 5%",
                 "볼라네소르센 apoC-III",
                 "볼라네소르센 TG (3개월)",
                 "볼라네소르센 혈소판 <100",
                 "올레자르센 TG (12개월, 공개연장)",
                 "올레자르센 췌장염 사건",
                 "플로자시란 TG (10개월)",
                 "플로자시란 췌장염 사건",
                 "피브레이트+ω-3 (LPL-null)",
                 "위약군 췌장염 발생률"),
      published = c("< 150 mg/dL", "~600-800 mg/dL", "1000-2000 mg/dL",
                    "> 5000 mg/dL", "300-500 mg/dL", "-76 ± 8 %", "-77 %",
                    "약 47 % 환자", "-73.7 %", "1건 vs 위약 11건", "약 -80 %",
                    "-83 %", "-10 ~ -20 %", "0.20-0.30 /patient-year"),
      source = c("일반 참고치", "Williams 2018 / 임상경험", "Brahm 2015",
                 "Gaudet 2014", "Rahalkar 2009", "Witztum NEJM 2019 (APPROACH)",
                 "Witztum NEJM 2019", "Witztum NEJM 2019",
                 "Stroes NEJM 2024 (Balance)", "Stroes NEJM 2024",
                 "Watts NEJM 2024/25 (PALISADE)", "Watts NEJM 2024/25",
                 "Falko 2018 / Gaudet 2014", "Balance·PALISADE 위약군"))
    DT::datatable(tb, rownames = FALSE, options = list(pageLength = 14, dom = "t"))
  })
}

shinyApp(ui, server)
