# =============================================================================
#  nhb_shiny_app.R
#  Neonatal hyperbilirubinaemia QSP model — interactive dashboard
#  신생아 고빌리루빈혈증 QSP 모델 · Shiny 대시보드
#  ---------------------------------------------------------------------------
#  Run:
#     shiny::runApp("nhb_shiny_app.R")
#  Requires nhb_mrgsolve_model.R in the same directory (it is sourced for the
#  compiled model object `mod`, the genotype table and the dosing helpers).
#
#  The app is built around the model's central claim: the number on the chart
#  (TSB) is not the number that injures (free bilirubin).  Tab 3 exists to make
#  that visible, and every other tab reports both.
#
#  NOT FOR CLINICAL USE.  Thresholds are analytic approximations of the AAP
#  2022 curves and must not be used to make care decisions.
# =============================================================================

library(shiny)
library(mrgsolve)
library(dplyr)
library(tidyr)
library(ggplot2)

source("nhb_mrgsolve_model.R", local = TRUE)   # provides mod, GENOTYPES, doses

THEME <- theme_bw(base_size = 12) +
  theme(panel.grid.minor = element_blank(),
        strip.background = element_rect(fill = "grey92"),
        legend.position = "bottom")

# -----------------------------------------------------------------------------
#  UI
# -----------------------------------------------------------------------------
ui <- fluidPage(
  titlePanel(
    div(
      h3("신생아 고빌리루빈혈증 QSP 시뮬레이터 · Neonatal Hyperbilirubinaemia QSP"),
      h5(em(paste("34-ODE mrgsolve model · haem catabolism → albumin binding →",
                  "UGT1A1 → enterohepatic shunt → photochemistry → BBB"))),
      h6(strong("교육·연구 목적 전용 — 임상 판단에 사용하지 마십시오 / NOT FOR CLINICAL USE"),
         style = "color:#b03a2e")
    )
  ),

  sidebarLayout(
    sidebarPanel(
      width = 3,

      h4("1 · 환자 프로파일 / Infant"),
      sliderInput("ga", "재태 주수 · gestational age (wk)", 28, 42, 40, 1),
      sliderInput("bw", "출생 체중 · birth weight (kg)", 0.8, 4.8, 3.40, 0.05),
      sliderInput("alb", "혈청 알부민 · albumin (g/dL)", 1.8, 4.5, 3.40, 0.1),
      sliderInput("hb", "제대 혈색소 · cord haemoglobin (g/dL)", 8, 22, 17, 0.5),
      selectInput("geno", "UGT1A1 유전형 · genotype",
                  choices = c("wild type (*1/*1)"        = "wild",
                              "*28 heterozygote"          = "UGT1A1_28_het",
                              "Gilbert (*28/*28)"         = "Gilbert",
                              "*6/*6 (East Asian, G71R)"  = "UGT1A1_6",
                              "Crigler-Najjar II"         = "CN2",
                              "Crigler-Najjar I"          = "CN1"),
                  selected = "wild"),

      h4("2 · 용혈 부하 / Haemolytic load"),
      sliderInput("abmat", "동종면역 항체 부하 · alloantibody load", 0, 0.6, 0, 0.02),
      checkboxInput("g6pd", "G6PD 결핍 · G6PD deficient", FALSE),
      conditionalPanel(
        "input.g6pd == true",
        sliderInput("oxstart", "산화 자극 시작 · oxidant challenge start (h)",
                    0, 168, 60, 6),
        sliderInput("oxdur", "지속 시간 · duration (h)", 0, 48, 24, 6)),
      sliderInput("cephal", "두혈종 혈액량 · extravasated Hb (g)", 0, 20, 0, 1),

      h4("3 · 결합·장벽 수식인자 / Binding & barrier"),
      sliderInput("facid", "산증 계수 · acidaemia factor on KA", 0.5, 1.0, 1.0, 0.05),
      sliderInput("fdisp", "치환 약물 계수 · displacer factor", 0.5, 1.0, 1.0, 0.05),
      sliderInput("fbbb", "BBB 투과성 배수 · BBB permeability", 1.0, 3.0, 1.0, 0.1),

      h4("4 · 수유 / Feeding"),
      radioButtons("breast", "수유 방식 · feeding",
                   c("모유 · human milk" = 1, "분유 · formula" = 0), 1, inline = TRUE),
      sliderInput("fiearly", "초기 섭취 적절도 · early intake adequacy",
                  0.1, 1.0, 1.0, 0.05),
      sliderInput("fisw", "개선 시점 · intake improves at (h)", 0, 240, 96, 12),

      h4("5 · 치료 / Therapy"),
      radioButtons("ptmode", "광선치료 방식 · phototherapy mode",
                   c("없음 · none"                       = "none",
                     "임계값 기반 · threshold-driven"     = "auto",
                     "고정 구간 · fixed window"           = "fixed"),
                   "auto"),
      conditionalPanel(
        "input.ptmode != 'none'",
        sliderInput("irr", "분광 조도 · irradiance (µW/cm²/nm)", 0, 60, 30, 1),
        sliderInput("fbsa", "노출 체표면적 비율 · exposed BSA", 0.2, 1.0, 0.80, 0.05)),
      conditionalPanel(
        "input.ptmode == 'fixed'",
        sliderInput("ptwin", "광선치료 구간 · window (h)", 0, 336, c(24, 96), 6)),
      conditionalPanel(
        "input.ptmode == 'auto'",
        sliderInput("ptduty", "1일 최대 시간 · max hours/day", 6, 24, 24, 2)),

      checkboxInput("ivig", "IVIG 1 g/kg", FALSE),
      conditionalPanel("input.ivig == true",
                       sliderInput("tivig", "투여 시각 · time (h)", 0, 168, 20, 2)),
      checkboxInput("et", "이중용적 교환수혈 · double-volume exchange", FALSE),
      conditionalPanel("input.et == true",
                       sliderInput("tet", "시작 · start (h)", 6, 240, 30, 2),
                       sliderInput("etdur", "소요 시간 · duration (h)", 1.5, 5, 3, 0.5)),
      checkboxInput("snmp", "스타노포르핀 · stannsoporfin 4.5 mg/kg IM", FALSE),
      checkboxInput("pheno", "페노바르비탈 · phenobarbital 5 mg/kg/d", FALSE),
      checkboxInput("udca", "UDCA 10 mg/kg q12h", FALSE),
      checkboxInput("agar", "경구 아가 · oral agar 250 mg/kg/d", FALSE),
      checkboxInput("aav", "AAV8-hUGT1A1 유전자 치료 · gene transfer", FALSE),
      conditionalPanel("input.aav == true",
                       sliderInput("taav", "투여 시각 · time (days)", 1, 90, 60, 1),
                       sliderInput("tgxmax", "발현 수준 · transgene activity (% adult)",
                                   1, 40, 10, 1)),

      hr(),
      sliderInput("tend", "시뮬레이션 기간 · horizon (days)", 3, 60, 14, 1),
      actionButton("run", "시뮬레이션 실행 · Simulate",
                   class = "btn-primary", width = "100%")
    ),

    mainPanel(
      width = 9,
      tabsetPanel(
        id = "tabs", type = "tabs",

        # ------------------------------------------------------------- TAB 1
        tabPanel(
          "① 빌리루빈 궤적 · Bilirubin & thresholds",
          br(),
          plotOutput("p_traj", height = "430px"),
          fluidRow(
            column(6, h5("판단 요약 · decision summary"),
                   tableOutput("t_decision")),
            column(6, h5("체중·성장 · weight"), plotOutput("p_wt", height = "220px"))
          ),
          helpText(paste("파선 = AAP 2022 광선치료 임계값 근사, 점선 = 교환수혈 임계값.",
                         "임계값은 재태주수·시간·위험인자의 함수입니다."))
        ),

        # ------------------------------------------------------------- TAB 2
        tabPanel(
          "② 유리 빌리루빈 · Free bilirubin",
          br(),
          fluidRow(
            column(6, plotOutput("p_bf", height = "300px")),
            column(6, plotOutput("p_bf_vs_tsb", height = "300px"))
          ),
          fluidRow(
            column(6, plotOutput("p_ba", height = "260px")),
            column(6, h5("등위험(iso-Bf) 곡선 · iso-risk contour"),
                   plotOutput("p_iso", height = "260px"))
          ),
          helpText(paste("이 탭이 모델의 핵심입니다. 측정되는 값(TSB)과 손상을 일으키는",
                         "값(Bf)의 관계는 포화 결합 등온선이며, 알부민이 낮거나 친화도가",
                         "떨어지면 같은 TSB가 전혀 다른 위험을 뜻합니다."))
        ),

        # ------------------------------------------------------------- TAB 3
        tabPanel(
          "③ 플럭스 균형 · Production vs clearance",
          br(),
          plotOutput("p_flux", height = "340px"),
          fluidRow(
            column(6, plotOutput("p_prod", height = "260px")),
            column(6, plotOutput("p_ugt", height = "260px"))
          ),
          helpText(paste("TSB는 두 플럭스의 차이의 적분이므로 그 자체로는 원인을",
                         "구별하지 못합니다. ETCOc(호기말 CO)는 생성 플럭스를 직접",
                         "읽어내므로 생성 병변과 제거 병변을 분리합니다."))
        ),

        # ------------------------------------------------------------- TAB 4
        tabPanel(
          "④ 광선치료 용량-반응 · Phototherapy",
          br(),
          fluidRow(
            column(7, plotOutput("p_ptgrid", height = "340px")),
            column(5, plotOutput("p_photoflux", height = "340px"))
          ),
          h5("광이성질체는 '빌리루빈'으로 측정됩니다 · photoisomers are assayed as bilirubin"),
          plotOutput("p_isomer", height = "300px"),
          helpText(paste("조도를 30 → 50 µW/cm²/nm로 올리는 것보다 노출 체표면적을",
                         "0.35 → 0.80으로 넓히는 것이 훨씬 효과적입니다(광학적 포화)."))
        ),

        # ------------------------------------------------------------- TAB 5
        tabPanel(
          "⑤ 교환수혈 · Exchange transfusion",
          br(),
          plotOutput("p_et", height = "420px"),
          tableOutput("t_et"),
          helpText(paste("교환수혈은 '빠른 광선치료'가 아닙니다. 공혈 혈장이 알부민을",
                         "재설정하고 모체 항체를 제거하므로 부하·결합능·생성 세 항에",
                         "동시에 작용하며, 그래서 Bf의 하강폭이 TSB의 하강폭보다 큽니다."))
        ),

        # ------------------------------------------------------------- TAB 6
        tabPanel(
          "⑥ 약물 PK/PD · Drug exposure",
          br(),
          fluidRow(
            column(6, plotOutput("p_pk1", height = "300px")),
            column(6, plotOutput("p_pk2", height = "300px"))
          ),
          plotOutput("p_pdeff", height = "280px"),
          helpText(paste("페노바르비탈은 UGT1A1 전사를 유도하므로 CN-I(GENO = 0)에서는",
                         "정의상 무효입니다. 존재하지 않는 유전자 산물은 유도할 수 없습니다."))
        ),

        # ------------------------------------------------------------- TAB 7
        tabPanel(
          "⑦ 신경독성 · Neurotoxicity",
          br(),
          fluidRow(
            column(6, plotOutput("p_brain", height = "300px")),
            column(6, plotOutput("p_inj", height = "300px"))
          ),
          fluidRow(
            column(6, plotOutput("p_abr", height = "260px")),
            column(6, h5("위험 요약 · risk summary"), tableOutput("t_risk"))
          ),
          helpText(paste("ABR 결손은 가역적으로, 누적 손상은 비가역적으로 모델링됩니다.",
                         "치료 후 청각 뇌간 반응이 정상화되는 임상 관찰과 대응합니다."))
        ),

        # ------------------------------------------------------------- TAB 8
        tabPanel(
          "⑧ 바이오마커 · Biomarkers",
          br(),
          fluidRow(
            column(6, plotOutput("p_etco", height = "290px")),
            column(6, plotOutput("p_hb", height = "290px"))
          ),
          fluidRow(
            column(6, plotOutput("p_tcb", height = "270px")),
            column(6, plotOutput("p_stool", height = "270px"))
          ),
          helpText(paste("ETCOc는 헴 분해 플럭스의 직접 지표(정상 ~1.4 ppm, 용혈 2-4 ppm),",
                         "TcB는 피부 구획을 반영하므로 광선치료 중에는 혈청과 해리됩니다."))
        ),

        # ------------------------------------------------------------- TAB 9
        tabPanel(
          "⑨ 시나리오 비교 · Scenario library",
          br(),
          checkboxGroupInput("scen", "비교할 시나리오 · scenarios",
                             choices = names(scenarios),
                             selected = c("S1", "S3", "S5", "S6", "S7"),
                             inline = TRUE),
          actionButton("runscen", "시나리오 실행 · Run library", class = "btn-info"),
          br(), br(),
          plotOutput("p_scen", height = "420px"),
          tableOutput("t_scen")
        ),

        # ------------------------------------------------------------ TAB 10
        tabPanel(
          "⑩ Crigler-Najjar 성장 분석 · Growth decomposition",
          br(),
          h5(paste("광선치료의 한계는 왜 아이가 자라면서 올라가는가 —",
                   "체표면적/체중 항은 생성량 감소로 상쇄된다")),
          plotOutput("p_cn", height = "380px"),
          tableOutput("t_cn"),
          helpText(paste("[G] 기하학만 / [G+P] + kg당 생성량 감소 / [G+P+O] + 피부 광학.",
                         "BSA/W는 출생→18세에 0.42배로 줄지만 kg당 빌리루빈 생성량도",
                         "0.45배로 줄어 서로 상쇄되므로, 기하학만으로는 임상에서 관찰되는",
                         "광선치료 실패를 설명할 수 없습니다."))
        ),

        # ------------------------------------------------------------ TAB 11
        tabPanel(
          "⑪ 모델 정보 · Model card",
          br(),
          h4("상태 변수 34개 · 34 ODEs"),
          verbatimTextOutput("t_cmt"),
          h4("현재 파라미터 · current parameter set"),
          verbatimTextOutput("t_par"),
          h4("알려진 편향 · known biases"),
          tags$ul(
            tags$li(paste("이중용적 교환수혈이 체내 총 부하의 40-45 %를 제거합니다",
                          "(고전적 교육 수치 ~25 %). 혈관외 구획을 단일 균질 구획으로",
                          "둔 결과이며, 같은 이유로 CN-I의 광선치료 한계값이 낙관적입니다.")),
            tags$li(paste("TAUUGT는 UGT1A1 단백 성숙과 OATP1B1·ligandin 성숙을 하나의",
                          "시간상수로 묶어 TSB 궤적에 보정한 값입니다.")),
            tags$li("AAP 2022 임계값은 해석적 근사이며 표값이 아닙니다."),
            tags$li(paste("피부 광학 항(FOPTMIN, TAUOPT)의 방향성은 확립되어 있으나",
                          "값은 예시입니다."))
          )
        )
      )
    )
  )
)

# -----------------------------------------------------------------------------
#  SERVER
# -----------------------------------------------------------------------------
server <- function(input, output, session) {

  # -- assemble the parameter set from the sidebar ---------------------------
  pset <- reactive({
    p <- list(
      GA = input$ga, W0 = input$bw, ALB0 = input$alb, ALBSET = max(input$alb, 3.0),
      HB0 = input$hb, GENO = GENOTYPES[[input$geno]],
      ABMAT0 = input$abmat, BLEX0 = input$cephal,
      G6PD = as.numeric(input$g6pd),
      FACID = input$facid, FDISP = input$fdisp, FBBB = input$fbbb,
      BREAST = as.numeric(input$breast),
      FIEARLY = input$fiearly, FISW = input$fisw, FILATE = 1.0,
      # gestational age drives RBC lifespan, binding affinity and the barrier
      LRBC  = 50 + 30*pmin(1, pmax(0, (input$ga - 32)/8)),
      FMATK = 0.70 + 0.30*pmin(1, pmax(0, (input$ga - 30)/8)),
      # AAP neurotoxicity risk factors, as the guideline defines them
      RF = as.numeric(input$ga < 38 || input$alb < 3.0 || input$abmat > 0 ||
                      input$g6pd || input$facid < 1.0)
    )
    if (input$g6pd) {
      p$OXSTART <- input$oxstart; p$OXDUR <- input$oxdur; p$OXSET <- 1.0
    }
    if (input$ptmode == "auto") {
      p$AUTOPT <- 1; p$IRRSET <- input$irr; p$FBSASET <- input$fbsa
      p$PTOFFM <- 2.0
    } else if (input$ptmode == "fixed") {
      p$AUTOPT <- 0; p$IRRSET <- input$irr; p$FBSASET <- input$fbsa
      p$PTSTART <- input$ptwin[1]; p$PTSTOP <- input$ptwin[2]
    }
    if (input$et) {
      p$ETSTART <- input$tet; p$ETDUR <- input$etdur; p$ETRATE <- 170/input$etdur
    }
    if (input$aav) {
      p$GTSTART <- input$taav*24; p$TGXMAX <- input$tgxmax/100
    }
    p
  })

  eset <- reactive({
    wt <- input$bw
    e <- NULL
    add <- function(a, b) if (is.null(a)) b else c(a, b)
    if (input$ivig)  e <- add(e, dose_ivig(wt, input$tivig))
    if (input$snmp)  e <- add(e, dose_snmp(wt, 24))
    if (input$pheno) e <- add(e, dose_pheno(wt, 0, input$tend))
    if (input$udca)  e <- add(e, dose_udca(wt, 0, input$tend))
    if (input$agar)  e <- add(e, dose_agar(wt, 0, input$tend))
    e
  })

  sim <- eventReactive(input$run, {
    m <- param(mod, pset())
    end <- input$tend*24
    e <- eset()
    out <- if (is.null(e)) mrgsim(m, end = end, delta = 0.5, hmax = 0.25)
           else mrgsim(m, events = e, end = end, delta = 0.5, hmax = 0.25)
    as_tibble(out)
  }, ignoreNULL = FALSE)

  # ------------------------------------------------------------------ TAB 1
  output$p_traj <- renderPlot({
    d <- sim()
    shade <- d %>% mutate(on = IRRAD > 0)
    ggplot(d, aes(time/24)) +
      geom_area(data = filter(shade, on), aes(y = Inf), fill = "#5dade2",
                alpha = 0.12) +
      geom_line(aes(y = THRESHPT, colour = "AAP phototherapy threshold"),
                linetype = 2, linewidth = 0.7) +
      geom_line(aes(y = THRESHET, colour = "AAP exchange threshold"),
                linetype = 3, linewidth = 0.7) +
      geom_line(aes(y = TSBOUT, colour = "TSB (reported by the assay)"),
                linewidth = 1.1) +
      geom_line(aes(y = UCBNAT, colour = "native (4Z,15Z) pigment"),
                linewidth = 0.8) +
      geom_line(aes(y = TCB, colour = "transcutaneous TcB"), linewidth = 0.6) +
      scale_colour_manual(values = c(
        "TSB (reported by the assay)" = "#c0392b",
        "native (4Z,15Z) pigment"     = "#8e44ad",
        "transcutaneous TcB"          = "#27ae60",
        "AAP phototherapy threshold"  = "grey35",
        "AAP exchange threshold"      = "grey10")) +
      labs(x = "postnatal age (days)", y = "bilirubin (mg/dL)", colour = NULL,
           title = "Total serum bilirubin against this infant's own thresholds",
           subtitle = "shaded = lamps on") + THEME
  })

  output$p_wt <- renderPlot({
    ggplot(sim(), aes(time/24, WTOUT)) + geom_line(linewidth = 1) +
      labs(x = "days", y = "weight (kg)") + THEME
  })

  output$t_decision <- renderTable({
    d <- sim(); dt <- d$time[2] - d$time[1]
    tibble(
      quantity = c("peak TSB (mg/dL)", "time of peak (h)",
                   "peak free bilirubin (nM)", "phototherapy hours",
                   "hours above phototherapy threshold",
                   "hours above exchange threshold",
                   "AUC above threshold (mg·h/dL)",
                   "exchange volume given (mL/kg)"),
      value = c(sprintf("%.2f", max(d$TSBOUT)),
                sprintf("%.0f", d$time[which.max(d$TSBOUT)]),
                sprintf("%.1f", max(d$BF)),
                sprintf("%.0f", max(d$PTH)),
                sprintf("%.0f", sum(d$OVERPT)*dt),
                sprintf("%.0f", sum(d$OVERET)*dt),
                sprintf("%.0f", max(d$AUCX)),
                sprintf("%.0f", max(d$ETV))))
  })

  # ------------------------------------------------------------------ TAB 2
  output$p_bf <- renderPlot({
    ggplot(sim(), aes(time/24, BF)) +
      geom_hline(yintercept = 30, linetype = 2, colour = "#c0392b") +
      annotate("text", x = 0.2, y = 32, hjust = 0, size = 3.4,
               label = "~30 nM: injury accrual begins in this model") +
      geom_line(linewidth = 1.1, colour = "#2471a3") +
      labs(x = "postnatal age (days)", y = "free bilirubin (nM)",
           title = "The species that crosses the barrier") + THEME
  })

  output$p_bf_vs_tsb <- renderPlot({
    ggplot(sim(), aes(TSBOUT, BF, colour = time/24)) +
      geom_path(linewidth = 1) +
      geom_hline(yintercept = 30, linetype = 2) +
      scale_colour_viridis_c(name = "days") +
      labs(x = "TSB, what is measured (mg/dL)", y = "Bf, what injures (nM)",
           title = "A convex map, not a proportionality") + THEME
  })

  output$p_ba <- renderPlot({
    ggplot(sim(), aes(time/24)) +
      geom_line(aes(y = BARATIO), linewidth = 1, colour = "#7d3c98") +
      geom_hline(yintercept = 0.8, linetype = 2) +
      labs(x = "days", y = "bilirubin / albumin molar ratio",
           title = "B/A: albumin-free, but NOT affinity-free") + THEME
  })

  output$p_iso <- renderPlot({
    grid <- expand.grid(tsb = seq(2, 30, 0.5),
                        alb = c(2.0, 2.5, 3.0, 3.5, 4.0))
    grid$bf <- free_bilirubin(grid$tsb, grid$alb,
                              fk = input$facid*input$fdisp)
    ggplot(grid, aes(tsb, bf, colour = factor(alb))) +
      geom_line(linewidth = 0.9) +
      geom_hline(yintercept = 30, linetype = 2) +
      coord_cartesian(ylim = c(0, 150)) +
      labs(x = "TSB (mg/dL)", y = "free bilirubin (nM)",
           colour = "albumin\n(g/dL)",
           title = "Iso-risk: the TSB that carries equal Bf") + THEME
  })

  # ------------------------------------------------------------------ TAB 3
  output$p_flux <- renderPlot({
    sim() %>%
      select(time, production = PROD24) %>%
      ggplot(aes(time/24, production)) + geom_line(linewidth = 1) +
      labs(x = "days", y = "bilirubin production (mg/kg/day)",
           title = "The production arm of the flux balance",
           subtitle = "term reference ~7.7 mg/kg/day, roughly twice the adult 3.8") +
      THEME
  })

  output$p_prod <- renderPlot({
    ggplot(sim(), aes(time/24, ETCOC)) + geom_line(linewidth = 1, colour = "#c0392b") +
      geom_hline(yintercept = 1.7, linetype = 2) +
      labs(x = "days", y = "ETCOc (ppm)",
           title = "End-tidal CO reads the production arm directly") + THEME
  })

  output$p_ugt <- renderPlot({
    ggplot(sim(), aes(time/24, UGTPCT)) + geom_line(linewidth = 1, colour = "#1e8449") +
      labs(x = "days", y = "UGT1A1 activity (% of adult)",
           title = "The clearance arm: ontogeny × genotype × induction") + THEME
  })

  # ------------------------------------------------------------------ TAB 4
  output$p_ptgrid <- renderPlot({
    grid <- expand.grid(IRRSET = c(0, 8, 15, 30, 50),
                        FBSASET = c(0.35, 0.60, 0.80, 1.00))
    res <- bind_rows(lapply(seq_len(nrow(grid)), function(i) {
      m <- param(mod, c(pset(), list(AUTOPT = 0, PTSTART = 24, PTSTOP = 1e9,
                                     IRRSET = grid$IRRSET[i],
                                     FBSASET = grid$FBSASET[i])))
      o <- as_tibble(mrgsim(m, end = 96, delta = 2, hmax = 0.25))
      tibble(irr = grid$IRRSET[i], fbsa = grid$FBSASET[i],
             change = 100*(o$TSBOUT[o$time == 48]/o$TSBOUT[o$time == 24] - 1))
    }))
    ggplot(res, aes(irr, change, colour = factor(fbsa))) +
      geom_line(linewidth = 1) + geom_point() +
      labs(x = "spectral irradiance (µW/cm²/nm)",
           y = "change in reported TSB over 24 h (%)",
           colour = "exposed\nBSA fraction",
           title = "Area beats irradiance once the skin saturates") + THEME
  })

  output$p_photoflux <- renderPlot({
    ggplot(sim(), aes(time/24, PHOTOFLX)) + geom_line(linewidth = 1) +
      labs(x = "days", y = "irreversible photochemical removal (mg/h)",
           title = "Lumirubin + E-isomer into bile",
           subtitle = "the only bilirubin exit that does not need UGT1A1") + THEME
  })

  output$p_isomer <- renderPlot({
    sim() %>%
      select(time, `reported TSB` = TSBOUT, `native pigment` = UCBNAT,
             photoisomers = ISOMER) %>%
      pivot_longer(-time) %>%
      ggplot(aes(time/24, value, colour = name)) + geom_line(linewidth = 1) +
      labs(x = "days", y = "mg/dL", colour = NULL,
           title = "Photoisomers are counted as bilirubin by the assay") + THEME
  })

  # ------------------------------------------------------------------ TAB 5
  output$p_et <- renderPlot({
    d <- sim()
    d %>% select(time, TSB = TSBOUT, `free bilirubin (nM)` = BF,
                 `albumin (g/dL)` = ALB, `haemoglobin (g/dL)` = HB,
                 `alloantibody load` = ABMAT) %>%
      pivot_longer(-time) %>%
      ggplot(aes(time, value)) + geom_line(linewidth = 0.9) +
      facet_wrap(~name, scales = "free_y", ncol = 2) +
      labs(x = "postnatal age (h)", y = NULL,
           title = "Exchange transfusion acts on load, capacity and production") +
      THEME
  })

  output$t_et <- renderTable({
    d <- sim()
    if (max(d$ETV) == 0) return(tibble(note = "no exchange transfusion in this run"))
    t0 <- input$tet; t1 <- input$tet + input$etdur
    pre  <- d$TSBOUT[which.min(abs(d$time - t0))]
    post <- d$TSBOUT[which.min(abs(d$time - t1))]
    reb  <- d$TSBOUT[which.min(abs(d$time - (t1 + 6)))]
    bfpre  <- d$BF[which.min(abs(d$time - t0))]
    bfpost <- d$BF[which.min(abs(d$time - t1))]
    tibble(quantity = c("pre-exchange TSB", "immediately post",
                        "removed (%)", "rebound at +6 h (% of pre)",
                        "pre-exchange Bf (nM)", "post-exchange Bf (nM)",
                        "Bf removed (%)"),
           value = sprintf("%.2f", c(pre, post, 100*(1-post/pre),
                                     100*reb/pre, bfpre, bfpost,
                                     100*(1-bfpost/bfpre))))
  })

  # ------------------------------------------------------------------ TAB 6
  output$p_pk1 <- renderPlot({
    sim() %>% select(time, `IVIG central (g)` = IGGC,
                     `stannsoporfin (µg/mL)` = CSNMP) %>%
      pivot_longer(-time) %>%
      ggplot(aes(time/24, value)) + geom_line(linewidth = 1) +
      facet_wrap(~name, scales = "free_y", ncol = 1) +
      labs(x = "days", y = NULL, title = "Production-arm drugs") + THEME
  })

  output$p_pk2 <- renderPlot({
    sim() %>% select(time, `phenobarbital (mg/L)` = CPB,
                     `UDCA (µmol/L)` = CUDCA) %>%
      pivot_longer(-time) %>%
      ggplot(aes(time/24, value)) + geom_line(linewidth = 1) +
      facet_wrap(~name, scales = "free_y", ncol = 1) +
      labs(x = "days", y = NULL, title = "Clearance-arm drugs") + THEME
  })

  output$p_pdeff <- renderPlot({
    sim() %>% select(time, `UGT1A1 (% adult)` = UGTPCT,
                     `transgene (fraction)` = TGX,
                     `luminal binder (mg)` = GBIND) %>%
      pivot_longer(-time) %>%
      ggplot(aes(time/24, value)) + geom_line(linewidth = 1) +
      facet_wrap(~name, scales = "free_y", ncol = 3) +
      labs(x = "days", y = NULL, title = "Pharmacodynamic states") + THEME
  })

  # ------------------------------------------------------------------ TAB 7
  output$p_brain <- renderPlot({
    ggplot(sim(), aes(time/24, BBR)) +
      geom_hline(yintercept = 35, linetype = 2, colour = "#c0392b") +
      geom_line(linewidth = 1.1, colour = "#943126") +
      labs(x = "days", y = "brain bilirubin (nM-equivalent)",
           title = "Basal-ganglia bilirubin",
           subtitle = "dashed = injury threshold") + THEME
  })

  output$p_inj <- renderPlot({
    sim() %>% select(time, `cumulative injury` = INJ,
                     `P(kernicterus)` = KERN) %>%
      pivot_longer(-time) %>%
      ggplot(aes(time/24, value, colour = name)) + geom_line(linewidth = 1.1) +
      coord_cartesian(ylim = c(0, 1)) +
      labs(x = "days", y = NULL, colour = NULL,
           title = "Irreversible accrual and its logistic read-out") + THEME
  })

  output$p_abr <- renderPlot({
    ggplot(sim(), aes(time/24, ABRD)) + geom_line(linewidth = 1) +
      labs(x = "days", y = "ABR deficit (relative)",
           title = "Auditory brainstem response — a reversible marker") + THEME
  })

  output$t_risk <- renderTable({
    d <- sim()
    tibble(quantity = c("peak brain bilirubin (nM-eq)",
                        "hours above the injury threshold",
                        "final cumulative injury index",
                        "P(kernicterus)",
                        "peak ABR deficit",
                        "AAP neurotoxicity risk factor flag"),
           value = c(sprintf("%.1f", max(d$BBR)),
                     sprintf("%.0f", sum(d$BBR > 35)*(d$time[2]-d$time[1])),
                     sprintf("%.3f", last(d$INJ)),
                     sprintf("%.3f", last(d$KERN)),
                     sprintf("%.3f", max(d$ABRD)),
                     ifelse(pset()$RF > 0, "present", "absent")))
  })

  # ------------------------------------------------------------------ TAB 8
  output$p_etco <- renderPlot({
    ggplot(sim(), aes(time/24, ETCOC)) + geom_line(linewidth = 1) +
      geom_hline(yintercept = c(1.3, 1.7), linetype = 3) +
      labs(x = "days", y = "ETCOc (ppm)",
           title = "End-tidal CO", subtitle = "dotted = normal band") + THEME
  })
  output$p_hb <- renderPlot({
    sim() %>% select(time, `haemoglobin (g/dL)` = HB,
                     `reticulocytes (%)` = RET) %>%
      pivot_longer(-time) %>%
      ggplot(aes(time/24, value)) + geom_line(linewidth = 1) +
      facet_wrap(~name, scales = "free_y", ncol = 1) +
      labs(x = "days", y = NULL, title = "Erythron") + THEME
  })
  output$p_tcb <- renderPlot({
    ggplot(sim(), aes(TSBOUT, TCB, colour = IRRAD > 0)) + geom_point(size = 0.7) +
      geom_abline(slope = 1, intercept = 0, linetype = 2) +
      labs(x = "TSB (mg/dL)", y = "TcB (mg/dL)", colour = "lamps on",
           title = "TcB dissociates from TSB under the lamps") + THEME
  })
  output$p_stool <- renderPlot({
    sim() %>% select(time, `cumulative output (mg)` = BSTL,
                     `gut unconjugated (mg)` = BGU) %>%
      pivot_longer(-time) %>%
      ggplot(aes(time/24, value)) + geom_line(linewidth = 1) +
      facet_wrap(~name, scales = "free_y", ncol = 1) +
      labs(x = "days", y = NULL, title = "The enterohepatic shunt") + THEME
  })

  # ------------------------------------------------------------------ TAB 9
  scen_sim <- eventReactive(input$runscen, {
    sel <- input$scen
    if (length(sel) == 0) return(NULL)
    bind_rows(lapply(sel, function(s) run_scenario(scenarios[[s]], end = 336)))
  })

  output$p_scen <- renderPlot({
    d <- scen_sim(); if (is.null(d)) return(NULL)
    ggplot(d, aes(time/24, TSBOUT, colour = scenario)) +
      geom_line(aes(y = THRESHPT), linetype = 2, colour = "grey55",
                show.legend = FALSE) +
      geom_line(linewidth = 1) +
      facet_wrap(~scenario, ncol = 2) +
      labs(x = "postnatal age (days)", y = "TSB (mg/dL)") +
      THEME + theme(legend.position = "none")
  })

  output$t_scen <- renderTable({
    d <- scen_sim(); if (is.null(d)) return(NULL)
    summarise_scenarios(d)
  })

  # ----------------------------------------------------------------- TAB 10
  cn <- eventReactive(input$run, { cn_ceiling(days = 21) })

  output$p_cn <- renderPlot({
    cn() %>%
      pivot_longer(c(TSB_G, TSB_GP, TSB_GPO)) %>%
      mutate(name = recode(name,
                           TSB_G   = "[G] geometry only",
                           TSB_GP  = "[G+P] + production ontogeny",
                           TSB_GPO = "[G+P+O] + skin optics")) %>%
      ggplot(aes(age_y, value, colour = name)) +
      geom_line(linewidth = 1.1) + geom_point() +
      labs(x = "age (years)", y = "phototherapy ceiling TSB (mg/dL)",
           colour = NULL,
           title = "Crigler-Najjar type I on 12 h/day intensive phototherapy",
           subtitle = paste("geometry alone predicts a rising ceiling; adding the",
                            "fall in production per kg cancels it")) + THEME
  })

  output$t_cn <- renderTable({ cn() })

  # ----------------------------------------------------------------- TAB 11
  output$t_cmt <- renderText(paste(names(init(mod)), collapse = " · "))
  output$t_par <- renderPrint({ str(pset()) })
}

shinyApp(ui, server)

# =============================================================================
#  END
# =============================================================================
