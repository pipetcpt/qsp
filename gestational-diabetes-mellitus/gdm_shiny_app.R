# =============================================================================
#  gdm_shiny_app.R
#  Interactive dashboard for the Gestational Diabetes Mellitus QSP model
#  임신성 당뇨병 QSP 모델 · Shiny 대시보드 (9 tabs)
#  ---------------------------------------------------------------------------
#  Run:  shiny::runApp("gdm_shiny_app.R")
#  Needs gdm_mrgsolve_model.R in the same directory (it is sourced for the
#  compiled model object `gdm` and the event/helper functions).
#
#  DESIGN INTENT
#  -------------
#  The app is built around the model's central claim: GDM is a race between the
#  placental contra-insulin clock and the beta-cell compensation clock. So the
#  primary control is NOT "how high is her glucose" — it is BCAP (beta-cell
#  adaptive capacity) and pre-gestational BMI. Glucose is an OUTPUT. Every tab
#  is arranged to make that causal direction visible, and the fetal-exposure tab
#  exists because that is the only axis on which the three drug options
#  genuinely differ.
# =============================================================================

library(shiny)
library(mrgsolve)
library(dplyr)
library(tidyr)
library(ggplot2)

# ---- load the model ---------------------------------------------------------
# Sourcing the model file compiles `gdm` and defines meal_events(), met_events(),
# insulin_events(), glyb_events(), combine_ev(), fasting_series(),
# titrate_insulin(), run_scenario(), at_delivery().
source("gdm_mrgsolve_model.R", local = FALSE, echo = FALSE)

THEME <- theme_bw(base_size = 12) +
  theme(panel.grid.minor = element_blank(),
        legend.position  = "bottom",
        plot.title       = element_text(face = "bold"))

TARGET_FPG <- 95    # mg/dL, ADA/ACOG fasting target in pregnancy
TARGET_1H  <- 140   # mg/dL, 1-hour post-prandial target

# =============================================================================
#  UI
# =============================================================================
ui <- fluidPage(
  titlePanel("임신성 당뇨병 (GDM) QSP 모델 — Gestational Diabetes Mellitus"),
  tags$p(style = "color:#555;",
         "태반 대항인슐린 구동 vs 베타세포 보상 — GDM 여부는 입력이 아니라 결과입니다.",
         tags$em("Whether GDM occurs is an output of beta-cell adaptive capacity, not an input.")),

  sidebarLayout(
    sidebarPanel(
      width = 3,

      h4("1. 환자 (Patient)"),
      sliderInput("bcap", "베타세포 적응능 BCAP",
                  min = 0.25, max = 1.20, value = 0.55, step = 0.05),
      helpText(tags$small("1.00 정상 · 0.85 분비결핍형 · 0.55 혼합형 · 0.30 중증")),
      sliderInput("bmi", "임신 전 BMI (kg/m²)",
                  min = 18, max = 42, value = 31, step = 0.5),
      sliderInput("age", "산모 연령 (세) — 위험도 표시용",
                  min = 18, max = 48, value = 34, step = 1),
      checkboxInput("lact", "산후 3개월 이상 수유 (lactation ≥3 mo)", FALSE),

      h4("2. 생활습관 (Lifestyle / MNT)"),
      sliderInput("cho", "일일 탄수화물 (g/day)",
                  min = 140, max = 280, value = 200, step = 5),
      checkboxInput("lowgi", "저혈당지수 식이 (low-GI, 흡수 완만화)", FALSE),
      sliderInput("ex", "운동 효과 EXEFF",
                  min = 0, max = 0.30, value = 0, step = 0.02),
      helpText(tags$small("0.20 ≈ 주 150분 이상 중등도 운동 (AMPK 경로, 인슐린 비의존)")),

      h4("3. 약물 (Pharmacotherapy)"),
      radioButtons("drug", NULL,
                   choices = c("없음 (MNT only)"        = "none",
                               "메트포민 (metformin)"   = "met",
                               "인슐린 (basal-bolus)"   = "ins",
                               "글리부라이드 (glyburide)" = "glyb",
                               "메트포민 + 인슐린 추가"  = "both"),
                   selected = "none"),
      sliderInput("tstart", "치료 시작 주수 (GW)",
                  min = 20, max = 36, value = 28, step = 1),
      conditionalPanel("input.drug == 'met' || input.drug == 'both'",
        sliderInput("metdose", "메트포민 1회 용량 (mg, BID)",
                    min = 250, max = 1250, value = 1000, step = 250),
        sliderInput("metadh", "복약 순응도 (GI 부작용 반영)",
                    min = 0.4, max = 1.0, value = 1.0, step = 0.05)),
      conditionalPanel("input.drug == 'ins' || input.drug == 'both'",
        checkboxInput("titrate", "공복혈당 95 mg/dL 목표로 자동 적정", TRUE),
        conditionalPanel("!input.titrate",
          sliderInput("insdose", "총 일일 인슐린 (U/day)",
                      min = 4, max = 140, value = 40, step = 2))),
      conditionalPanel("input.drug == 'glyb'",
        sliderInput("glybdose", "글리부라이드 1회 용량 (mg, BID)",
                    min = 1.25, max = 10, value = 5, step = 1.25)),

      h4("4. 분만 (Delivery)"),
      sliderInput("tdel", "분만 주수 (GW)",
                  min = 36, max = 41, value = 39, step = 0.5),

      hr(),
      actionButton("run", "시뮬레이션 실행 (Run)", class = "btn-primary",
                   width = "100%"),
      br(), br(),
      checkboxInput("addref", "정상 임신 대조군 함께 표시", TRUE),
      helpText(tags$small("첫 실행은 모델 컴파일로 20-40초 걸릴 수 있습니다."))
    ),

    mainPanel(
      width = 9,
      tabsetPanel(
        id = "tabs", type = "tabs",

        # ------------------------------------------------------------- TAB 1 --
        tabPanel(
          "1. 환자 프로파일",
          br(),
          fluidRow(
            column(6, h4("진단 요약 (Diagnostic summary)"),
                   tableOutput("tbl_profile")),
            column(6, h4("Powe 생리학적 아형 위치"),
                   plotOutput("plt_subtype", height = "300px"),
                   helpText(tags$small(
                     "Powe (Diabetes Care 2016;39:1052) 는 GDM을 인슐린저항 우세형·",
                     "분비결핍형·혼합형으로 구분했습니다. 이 모델에서 아형은 BCAP와 BMI로부터",
                     "자동으로 나타납니다.")))
          ),
          hr(),
          h4("임신 중 대사 적응 (Metabolic adaptation across gestation)"),
          plotOutput("plt_adapt", height = "330px")
        ),

        # ------------------------------------------------------------- TAB 2 --
        tabPanel(
          "2. 태반 내분비 축",
          br(),
          h4("태반 성장과 대항인슐린 호르몬 (Placental clock)"),
          plotOutput("plt_placenta", height = "420px"),
          helpText("hPL·프로게스테론·에스트라디올·TNF-α는 태반 질량을 따라가고, ",
                   "이들을 합산한 대항인슐린 지수(CID)가 인슐린 감수성을 결정합니다. ",
                   "이 곡선은 모든 임신에서 거의 동일합니다 — 다른 것은 베타세포 쪽입니다."),
          hr(),
          h4("아디포카인 (Adipokines)"),
          plotOutput("plt_adipo", height = "300px")
        ),

        # ------------------------------------------------------------- TAB 3 --
        tabPanel(
          "3. 모체 혈당 · 인슐린",
          br(),
          fluidRow(
            column(6, h4("공복혈당 궤적 (Fasting glucose)"),
                   plotOutput("plt_fpg", height = "310px")),
            column(6, h4("공복 인슐린 (Fasting insulin)"),
                   plotOutput("plt_ins", height = "310px"))
          ),
          hr(),
          h4("24시간 혈당 프로파일 (선택 주수) — CGM 유사 출력"),
          sliderInput("cgmweek", "표시 주수 (GW)", min = 12, max = 40,
                      value = 34, step = 1, width = "60%"),
          plotOutput("plt_cgm", height = "320px"),
          tableOutput("tbl_cgm")
        ),

        # ------------------------------------------------------------- TAB 4 --
        tabPanel(
          "4. 베타세포 보상",
          br(),
          h4("보상 vs 실패 (Compensation vs failure)"),
          plotOutput("plt_beta", height = "400px"),
          helpText("BCM = 베타세포 질량 x 기능, KG50E = 포도당 감수 곡선의 EC50 ",
                   "(좌측 이동 = 감수성 증가). DI = SI x BCM (disposition index)."),
          hr(),
          fluidRow(
            column(6, h4("Disposition index 궤적"),
                   plotOutput("plt_di", height = "300px")),
            column(6, h4("혈당 자극 인슐린 분비 곡선"),
                   plotOutput("plt_isr", height = "300px"))
          )
        ),

        # ------------------------------------------------------------- TAB 5 --
        tabPanel(
          "5. 약물 PK · 태아 노출",
          br(),
          h4("모체 대 태아 약물 농도 (Maternal vs fetal exposure)"),
          plotOutput("plt_pk", height = "360px"),
          tableOutput("tbl_pk"),
          hr(),
          tags$div(
            style = "background:#fff6f6; border-left:4px solid #c04040; padding:10px;",
            tags$b("이 탭이 존재하는 이유 (why this tab exists):"), br(),
            "세 약물은 모체 혈당 강하 효과에서는 서로 비슷합니다. 실제로 다른 것은 ",
            tags$b("태아 노출"), " 입니다 — 인슐린은 태반을 통과하지 않고(구조적으로 0), ",
            "메트포민은 자유 통과하여 umbilical:maternal ≈ 1.0, 글리부라이드는 ",
            "BCRP/ABCG2 유출에도 불구하고 cord:maternal ≈ 0.7 에 도달합니다."
          ),
          br(),
          h4("인슐린 요구량 (titration 결과)"),
          verbatimTextOutput("txt_titrate")
        ),

        # ------------------------------------------------------------- TAB 6 --
        tabPanel(
          "6. 태아 성장",
          br(),
          h4("태아 체중과 체성분 (Fetal weight and composition)"),
          plotOutput("plt_fetal", height = "360px"),
          helpText("회색 점선 = 기준(50백분위) 곡선. 과성장은 ", tags$b("비대칭적"),
                   "입니다 — 지방이 제지방보다 인슐린 탄성이 훨씬 큽니다."),
          hr(),
          fluidRow(
            column(6, h4("지방 vs 제지방 (Fat vs lean)"),
                   plotOutput("plt_comp", height = "300px")),
            column(6, h4("태아 혈당 · 인슐린 (Pedersen 연쇄)"),
                   plotOutput("plt_fetal_gi", height = "300px"))
          )
        ),

        # ------------------------------------------------------------- TAB 7 --
        tabPanel(
          "7. 임상 엔드포인트",
          br(),
          h4("분만 시점 엔드포인트 (Endpoints at delivery)"),
          tableOutput("tbl_endpoint"),
          plotOutput("plt_endpoint", height = "330px"),
          hr(),
          h4("보정 목표 대비 (against calibration targets)"),
          tableOutput("tbl_calib"),
          helpText("보정 근거: HAPO (NEJM 2008;358:1991) 연속적 혈당–LGA 기울기, ",
                   "Landon/MFMU 경증 GDM 치료 효과 (NEJM 2009;361:1339), ",
                   "ACHOIS (NEJM 2005;352:2477).")
        ),

        # ------------------------------------------------------------- TAB 8 --
        tabPanel(
          "8. 시나리오 비교",
          br(),
          h4("치료 전략 비교 (같은 환자 생리, 다른 전략)"),
          checkboxGroupInput(
            "cmp", NULL, inline = TRUE,
            choices = c("무치료" = "none", "MNT+운동" = "mnt",
                        "메트포민" = "met", "인슐린" = "ins",
                        "글리부라이드" = "glyb"),
            selected = c("none", "mnt", "ins")),
          actionButton("runcmp", "비교 실행 (Run comparison)", class = "btn-default"),
          br(), br(),
          plotOutput("plt_cmp", height = "330px"),
          tableOutput("tbl_cmp")
        ),

        # ------------------------------------------------------------- TAB 9 --
        tabPanel(
          "9. 바이오마커 · 산후",
          br(),
          fluidRow(
            column(6, h4("당화 지표 (HbA1c, 평균혈당)"),
                   plotOutput("plt_a1c", height = "300px")),
            column(6, h4("유리지방산 · 모체 지방량"),
                   plotOutput("plt_bio", height = "300px"))
          ),
          hr(),
          h4("산후 궤적과 5년 제2형 당뇨병 위험"),
          plotOutput("plt_pp", height = "330px"),
          tableOutput("tbl_pp"),
          helpText("분만 시 태반이 퇴축하면 CID가 붕괴하고 인슐린 감수성은 회복되지만, ",
                   "베타세포 결손은 회복되지 않습니다 — 산후 disposition index가 ",
                   "그 결손을 드러냅니다 (Bellamy Lancet 2009;373:1773, RR 7.43).")
        ),

        # ---------------------------------------------------------- TAB 10 --
        tabPanel(
          "10. HAPO 검증",
          br(),
          h4("HAPO 연속 기울기 재현 (continuous gradient, no threshold)"),
          actionButton("runhapo", "스윕 실행 (BCAP x BMI grid — 약 1-3분)",
                       class = "btn-default"),
          br(), br(),
          plotOutput("plt_hapo", height = "380px"),
          tableOutput("tbl_hapo"),
          helpText("관찰값: HAPO 범주 1 (FPG <75 mg/dL) LGA 5.3% → 범주 7 (≥100) 26.3%. ",
                   "모델이 재현해야 하는 것은 특정 절단점이 아니라 ", tags$b("역치 없는 연속성"),
                   " 입니다.")
        )
      )
    )
  )
)

# =============================================================================
#  SERVER
# =============================================================================
server <- function(input, output, session) {

  # ---- assemble parameters and events from the UI ---------------------------
  build <- function(drug = input$drug, bcap = input$bcap, bmi = input$bmi,
                    ex = input$ex, cho = input$cho, lowgi = input$lowgi) {
    tdel <- input$tdel*7
    pars <- list(BCAP = bcap, BMI = bmi, EXEFF = ex, CHOD = cho,
                 TDEL = tdel, LACT = as.numeric(input$lact),
                 KTR_G = if (lowgi) 16 else 26)
    end  <- max(tdel + 84, 364)
    ev   <- meal_events(cho, 56, end)
    titr <- NA_real_

    if (drug %in% c("met", "both")) {
      pars$METADH <- input$metadh
      ev <- combine_ev(ev, met_events(input$metdose, input$tstart, input$tdel))
    }
    if (drug %in% c("ins", "both")) {
      if (isTRUE(input$titrate)) {
        tt <- titrate_insulin(pars, TARGET_FPG, input$tstart, input$tdel, cho,
                              u_start = if (drug == "both") 12 else 20)
        titr <- tt$total_u
      } else {
        titr <- input$insdose
      }
      ev <- combine_ev(ev, insulin_events(titr, input$tstart, input$tdel))
    }
    if (drug == "glyb") {
      ev <- combine_ev(ev, glyb_events(input$glybdose, TRUE,
                                       input$tstart, input$tdel))
    }
    list(pars = pars, ev = ev, end = end, tdel = tdel, titr = titr)
  }

  simulate <- function(spec, delta = 0.05) {
    gdm |> param(spec$pars) |> data_set(spec$ev) |>
      mrgsim(start = 56, end = spec$end, delta = delta, hmax = 0.05) |>
      as.data.frame()
  }

  # ---- main reactive --------------------------------------------------------
  sim <- eventReactive(input$run, {
    withProgress(message = "시뮬레이션 중...", value = 0.3, {
      spec <- build()
      incProgress(0.4)
      d <- simulate(spec)
      d$arm <- "환자 (patient)"
      out <- list(spec = spec, main = d, titr = spec$titr)

      if (isTRUE(input$addref)) {
        rspec <- build(drug = "none", bcap = 1.00, bmi = 22, ex = 0,
                       cho = input$cho, lowgi = FALSE)
        r <- simulate(rspec)
        r$arm <- "정상 대조 (NGT)"
        out$ref <- r
      }
      incProgress(0.3)
      out
    })
  }, ignoreNULL = FALSE)

  allsim <- reactive({
    s <- sim()
    if (is.null(s$ref)) s$main else rbind(s$main, s$ref)
  })

  fast <- reactive({
    allsim() |>
      mutate(day = floor(time), frac = time - day) |>
      filter(frac > 0.20, frac < 0.29) |>
      group_by(arm, day) |>
      summarise(across(c(GLU, INS, MBG, HBA1C, FFA, BCM, SIP, DI, DIREF,
                         CID, KG50E, FATM, ADIPO, TNFA, HPL, PROG, E2,
                         PLAC, GLUF, INSF, BCF, FWT, FWREF, BWZ, FATPCT,
                         CMET, CMETF, CGLY, CGLYF, NHR, PLGA, PNHYP,
                         PSD, PCS, PPE, CIT2D, HYPOAUC),
                       mean), .groups = "drop") |>
      mutate(GW = day/7)
  })

  del <- reactive({
    s <- sim(); td <- s$spec$tdel
    allsim() |> group_by(arm) |> slice(which.min(abs(time - td))) |> ungroup()
  })

  # ========================================================== TAB 1 ==========
  output$tbl_profile <- renderTable({
    d <- del(); f <- fast()
    ft <- f |> filter(GW <= input$tdel) |> group_by(arm) |>
      summarise(FPG28 = mean(GLU[GW >= 27 & GW <= 29]),
                FPGterm = mean(tail(GLU, 7)),
                INSterm = mean(tail(INS, 7)), .groups = "drop")
    data.frame(
      항목 = c("BCAP (베타세포 적응능)", "임신 전 BMI", "산모 연령",
               "GW 28 공복혈당 (mg/dL)", "분만시 공복혈당 (mg/dL)",
               "분만시 공복인슐린 (µU/mL)", "인슐린 감수성 (임신전 대비 %)",
               "베타세포 질량x기능 BCM", "Disposition index (상대값)",
               "HbA1c (%)", "IADPSG 공복 기준(5.1 mmol/L=92) 초과"),
      값 = c(sprintf("%.2f", input$bcap), sprintf("%.1f", input$bmi),
             sprintf("%d", input$age),
             sprintf("%.1f", ft$FPG28[ft$arm == "환자 (patient)"]),
             sprintf("%.1f", ft$FPGterm[ft$arm == "환자 (patient)"]),
             sprintf("%.1f", ft$INSterm[ft$arm == "환자 (patient)"]),
             sprintf("%.0f%%", 100*d$MATSUDA[d$arm == "환자 (patient)"]),
             sprintf("%.2f", d$BCM[d$arm == "환자 (patient)"]),
             sprintf("%.2f", (d$DI/d$DIREF)[d$arm == "환자 (patient)"]),
             sprintf("%.2f", d$HBA1C[d$arm == "환자 (patient)"]),
             ifelse(ft$FPG28[ft$arm == "환자 (patient)"] > 92, "예 (YES)", "아니오 (no)")),
      check.names = FALSE)
  }, striped = TRUE, spacing = "xs")

  output$plt_subtype <- renderPlot({
    d <- del() |> filter(arm == "환자 (patient)")
    grid <- expand.grid(BCAP = seq(0.3, 1.15, 0.05), BMI = seq(20, 40, 1))
    grid$SI  <- 1/(1 + 1.40)*exp(-0.055*(grid$BMI - 22))
    grid$SEC <- grid$BCAP
    ggplot(grid, aes(SEC, SI)) +
      geom_raster(aes(fill = SEC*SI), alpha = 0.55) +
      scale_fill_gradient(low = "#c94040", high = "#3a7d44", guide = "none") +
      annotate("point", x = input$bcap,
               y = 1/(1 + 1.40)*exp(-0.055*(input$bmi - 22)),
               size = 5, shape = 21, fill = "white", stroke = 1.4) +
      annotate("text", x = 0.42, y = 0.40, label = "분비결핍형\n(deficient)",
               size = 3.4) +
      annotate("text", x = 1.05, y = 0.13, label = "인슐린저항형\n(resistant)",
               size = 3.4) +
      annotate("text", x = 0.45, y = 0.13, label = "혼합형\n(mixed)", size = 3.4) +
      labs(x = "분비 능력 (BCAP)", y = "말기 인슐린 감수성 (relative)",
           title = "생리학적 아형 좌표") + THEME
  })

  output$plt_adapt <- renderPlot({
    fast() |>
      select(arm, GW, `인슐린감수성 SI` = SIP, `베타세포 BCM` = BCM,
             `대항인슐린 CID` = CID, `공복혈당 FPG` = GLU) |>
      pivot_longer(-c(arm, GW)) |>
      ggplot(aes(GW, value, colour = arm)) +
      geom_line(linewidth = 0.8) +
      facet_wrap(~name, scales = "free_y") +
      labs(x = "임신 주수 (GW)", y = NULL, colour = NULL,
           title = "두 개의 시계: 태반 구동 vs 베타세포 보상") + THEME
  })

  # ========================================================== TAB 2 ==========
  output$plt_placenta <- renderPlot({
    fast() |>
      transmute(arm, GW,
                `태반 질량 (norm.)` = PLAC,
                `hPL (mg/L)` = HPL,
                `프로게스테론 (ng/mL)` = PROG,
                `에스트라디올 (ng/mL)` = E2/1000,
                `TNF-α (pg/mL)` = TNFA,
                `대항인슐린 지수 CID` = CID) |>
      pivot_longer(-c(arm, GW)) |>
      ggplot(aes(GW, value, colour = arm)) +
      geom_line(linewidth = 0.8) +
      geom_vline(xintercept = input$tdel, linetype = 3) +
      facet_wrap(~name, scales = "free_y", ncol = 3) +
      labs(x = "임신 주수 (GW)", y = NULL, colour = NULL,
           title = "태반 시계 — 모든 임신에서 거의 동일") + THEME
  })

  output$plt_adipo <- renderPlot({
    fast() |>
      transmute(arm, GW, `아디포넥틴 (µg/mL)` = ADIPO,
                `모체 지방량 (kg)` = FATM) |>
      pivot_longer(-c(arm, GW)) |>
      ggplot(aes(GW, value, colour = arm)) +
      geom_line(linewidth = 0.8) + facet_wrap(~name, scales = "free_y") +
      labs(x = "임신 주수 (GW)", y = NULL, colour = NULL) + THEME
  })

  # ========================================================== TAB 3 ==========
  output$plt_fpg <- renderPlot({
    ggplot(fast(), aes(GW, GLU, colour = arm)) +
      geom_line(linewidth = 0.9) +
      geom_hline(yintercept = TARGET_FPG, linetype = 2, colour = "#b03030") +
      geom_hline(yintercept = 92, linetype = 3, colour = "#3060b0") +
      geom_vline(xintercept = input$tdel, linetype = 3) +
      annotate("text", x = 14, y = TARGET_FPG + 3, size = 3,
               label = "치료 목표 95", colour = "#b03030") +
      annotate("text", x = 14, y = 89, size = 3,
               label = "IADPSG 92", colour = "#3060b0") +
      labs(x = "임신 주수 (GW)", y = "공복혈당 (mg/dL)", colour = NULL) + THEME
  })

  output$plt_ins <- renderPlot({
    ggplot(fast(), aes(GW, INS, colour = arm)) +
      geom_line(linewidth = 0.9) +
      geom_vline(xintercept = input$tdel, linetype = 3) +
      labs(x = "임신 주수 (GW)", y = "공복 인슐린 (µU/mL)", colour = NULL) + THEME
  })

  output$plt_cgm <- renderPlot({
    d0 <- floor(input$cgmweek*7)
    allsim() |> filter(time >= d0, time < d0 + 2) |>
      mutate(hour = (time - d0)*24) |>
      ggplot(aes(hour, GLU, colour = arm)) +
      geom_line(linewidth = 0.8) +
      geom_hline(yintercept = c(63, TARGET_1H), linetype = 2, colour = "grey40") +
      scale_x_continuous(breaks = seq(0, 48, 6)) +
      labs(x = sprintf("GW %d 이후 시간 (h)", input$cgmweek),
           y = "혈당 (mg/dL)", colour = NULL,
           title = "48시간 혈당 프로파일 (CGM 유사)",
           subtitle = "점선 = 63 및 140 mg/dL (TIR 범위)") + THEME
  })

  output$tbl_cgm <- renderTable({
    d0 <- floor(input$cgmweek*7)
    allsim() |> filter(time >= d0, time < d0 + 7) |>
      group_by(arm) |>
      summarise(`평균혈당 (mg/dL)` = round(mean(GLU), 1),
                `최고혈당` = round(max(GLU), 1),
                `최저혈당` = round(min(GLU), 1),
                `TIR 63-140 (%)` = round(100*mean(GLU >= 63 & GLU <= 140), 1),
                `>140 시간비율 (%)` = round(100*mean(GLU > 140), 1),
                `<70 시간비율 (%)` = round(100*mean(GLU < 70), 1),
                .groups = "drop")
  }, striped = TRUE, spacing = "xs")

  # ========================================================== TAB 4 ==========
  output$plt_beta <- renderPlot({
    fast() |>
      transmute(arm, GW, `베타세포 BCM` = BCM,
                `포도당 EC50 KG50E (mg/dL)` = KG50E,
                `Disposition index (rel.)` = DI/DIREF,
                `공복혈당 (mg/dL)` = GLU) |>
      pivot_longer(-c(arm, GW)) |>
      ggplot(aes(GW, value, colour = arm)) +
      geom_line(linewidth = 0.8) + facet_wrap(~name, scales = "free_y") +
      geom_vline(xintercept = input$tdel, linetype = 3) +
      labs(x = "임신 주수 (GW)", y = NULL, colour = NULL,
           title = "보상이 성공했는가, 실패했는가") + THEME
  })

  output$plt_di <- renderPlot({
    ggplot(fast(), aes(SIP, BCM, colour = arm)) +
      geom_path(linewidth = 0.9, arrow = arrow(length = unit(0.15, "cm"))) +
      labs(x = "인슐린 감수성 SI", y = "베타세포 BCM", colour = NULL,
           title = "쌍곡선 관계 위의 궤적",
           subtitle = "정상은 SI가 떨어질 때 BCM이 충분히 올라가 곡선을 따라 이동합니다") +
      THEME
  })

  output$plt_isr <- renderPlot({
    d <- del() |> filter(arm == "환자 (patient)")
    g <- seq(50, 250, 2)
    curves <- rbind(
      data.frame(GLU = g, arm = "분만시 (adapted)",
                 ISR = d$BCM*g^2/(d$KG50E^2 + g^2)),
      data.frame(GLU = g, arm = "임신 전 (non-adapted)",
                 ISR = 1.0*g^2/(130^2 + g^2)))
    ggplot(curves, aes(GLU, ISR, colour = arm)) +
      geom_line(linewidth = 0.9) +
      geom_vline(xintercept = d$GLU, linetype = 2) +
      labs(x = "혈당 (mg/dL)", y = "상대 인슐린 분비율", colour = NULL,
           title = "포도당 자극 인슐린 분비 곡선",
           subtitle = "임신은 곡선을 위로(질량) 그리고 왼쪽으로(감수성) 이동시킵니다") +
      THEME
  })

  # ========================================================== TAB 5 ==========
  output$plt_pk <- renderPlot({
    d <- allsim() |> filter(arm == "환자 (patient)",
                            time > input$tstart*7 - 1, time < input$tdel*7 + 1)
    long <- rbind(
      data.frame(GW = d$GWk, conc = d$CMET,  who = "모체 (maternal)", drug = "메트포민 (mg/L)"),
      data.frame(GW = d$GWk, conc = d$CMETF, who = "태아 (fetal)",    drug = "메트포민 (mg/L)"),
      data.frame(GW = d$GWk, conc = d$CGLY,  who = "모체 (maternal)", drug = "글리부라이드 (mg/L)"),
      data.frame(GW = d$GWk, conc = d$CGLYF, who = "태아 (fetal)",    drug = "글리부라이드 (mg/L)"))
    long <- long |> group_by(drug) |> filter(max(conc) > 1e-6) |> ungroup()
    if (nrow(long) == 0) {
      return(ggplot() + annotate("text", 0, 0, size = 5,
        label = "선택된 경구약이 없습니다.\n인슐린은 태반을 통과하지 않으므로 태아 노출은 구조적으로 0입니다.") +
        theme_void())
    }
    ggplot(long, aes(GW, conc, colour = who)) +
      geom_line(linewidth = 0.7) + facet_wrap(~drug, scales = "free_y") +
      labs(x = "임신 주수 (GW)", y = "혈중 농도", colour = NULL,
           title = "태아는 모체가 노출된 것을 함께 노출됩니다") + THEME
  })

  output$tbl_pk <- renderTable({
    d <- allsim() |> filter(arm == "환자 (patient)",
                            time > input$tdel*7 - 7, time <= input$tdel*7)
    data.frame(
      약물 = c("메트포민", "글리부라이드", "인슐린"),
      `모체 Cmax` = c(round(max(d$CMET), 3), round(max(d$CGLY), 4), NA),
      `모체 Cavg` = c(round(mean(d$CMET), 3), round(mean(d$CGLY), 4), NA),
      `태아 Cavg` = c(round(mean(d$CMETF), 3), round(mean(d$CGLYF), 4), 0),
      `태아:모체 비` = c(
        ifelse(mean(d$CMET) > 1e-6, round(mean(d$CMETF)/mean(d$CMET), 2), NA),
        ifelse(mean(d$CGLY) > 1e-8, round(mean(d$CGLYF)/mean(d$CGLY), 2), NA),
        0),
      `문헌 보고 비` = c("~1.0 (Vanky 2005)", "~0.7 (Hemauer 2010)", "0 (통과 안 함)"),
      check.names = FALSE)
  }, striped = TRUE, spacing = "xs")

  output$txt_titrate <- renderPrint({
    s <- sim()
    if (is.na(s$titr)) {
      cat("인슐린을 사용하지 않는 시나리오입니다.\n")
    } else {
      wt <- input$bmi*1.62^2
      cat(sprintf("총 일일 인슐린: %.0f U/day  (%.2f U/kg/day, 체중 %.0f kg 가정)\n",
                  s$titr, s$titr/wt, wt))
      cat(sprintf("기저:볼루스 = 50:50, 볼루스는 3식 분할\n"))
      cat("참고: 3분기 요구량은 통상 0.9-1.0 U/kg/day 로 상승합니다.\n")
      cat("이 값은 지정된 것이 아니라 공복혈당 95 mg/dL 목표에 도달하도록\n")
      cat("적정(titration)된 결과입니다 — 즉 환자 생리의 산물입니다.\n")
    }
  })

  # ========================================================== TAB 6 ==========
  output$plt_fetal <- renderPlot({
    f <- fast() |> filter(GW <= input$tdel)
    ggplot(f, aes(GW, FWT, colour = arm)) +
      geom_line(linewidth = 0.9) +
      geom_line(aes(y = FWREF), colour = "grey45", linetype = 2) +
      labs(x = "임신 주수 (GW)", y = "태아 체중 (g)", colour = NULL,
           title = "태아 성장 궤적",
           subtitle = "회색 점선 = 기준 50백분위") + THEME
  })

  output$plt_comp <- renderPlot({
    allsim() |> filter(time <= input$tdel*7) |>
      transmute(arm, GW = GWk, `지방 (g)` = FATF, `제지방 (g)` = LEANF) |>
      pivot_longer(-c(arm, GW)) |>
      ggplot(aes(GW, value, colour = arm, linetype = name)) +
      geom_line(linewidth = 0.8) +
      labs(x = "임신 주수 (GW)", y = "질량 (g)", colour = NULL, linetype = NULL,
           title = "비대칭 과성장: 지방이 탄성 구획") + THEME
  })

  output$plt_fetal_gi <- renderPlot({
    fast() |> filter(GW <= input$tdel) |>
      transmute(arm, GW, `태아 혈당 (mg/dL)` = GLUF,
                `태아 인슐린 (µU/mL)` = INSF,
                `태아 베타세포 BCF` = BCF,
                `태아 지방 비율 (%)` = FATPCT) |>
      pivot_longer(-c(arm, GW)) |>
      ggplot(aes(GW, value, colour = arm)) +
      geom_line(linewidth = 0.8) + facet_wrap(~name, scales = "free_y") +
      labs(x = "임신 주수 (GW)", y = NULL, colour = NULL,
           title = "Pedersen 연쇄") + THEME
  })

  # ========================================================== TAB 7 ==========
  ep_tbl <- reactive({
    del() |> transmute(
      arm,
      `출생체중 (g)` = round(FWT),
      `체중 z-score` = round(BWZ, 2),
      `태아 지방 (%)` = round(FATPCT, 1),
      `제대 인슐린 비 (NHR)` = round(NHR, 2),
      `LGA (%)` = round(100*PLGA, 1),
      `거대아 >4000g (%)` = round(100*PMACR, 1),
      `신생아 저혈당 (%)` = round(100*PNHYP, 1),
      `견부난산 (%)` = round(100*PSD, 1),
      `전자간증 (%)` = round(100*PPE, 1),
      `제왕절개 (%)` = round(100*PCS, 1),
      `모체 저혈당 노출 (mg/dL·d)` = round(HYPOAUC, 0),
      `5년 T2DM (%)` = round(100*CIT2D, 1))
  })

  output$tbl_endpoint <- renderTable(ep_tbl(), striped = TRUE, spacing = "xs")

  output$plt_endpoint <- renderPlot({
    del() |>
      transmute(arm, LGA = 100*PLGA, `거대아` = 100*PMACR,
                `신생아저혈당` = 100*PNHYP, `견부난산` = 100*PSD,
                `전자간증` = 100*PPE, `제왕절개` = 100*PCS) |>
      pivot_longer(-arm) |>
      ggplot(aes(reorder(name, -value), value, fill = arm)) +
      geom_col(position = "dodge") +
      labs(x = NULL, y = "확률 (%)", fill = NULL,
           title = "분만 시점 임상 엔드포인트") + THEME
  })

  output$tbl_calib <- renderTable({
    d <- del() |> filter(arm == "환자 (patient)")
    data.frame(
      엔드포인트 = c("LGA", "거대아 >4000 g", "견부난산", "전자간증/HDP", "제왕절개"),
      `모델 (%)` = c(round(100*d$PLGA, 1), round(100*d$PMACR, 1),
                     round(100*d$PSD, 1), round(100*d$PPE, 1),
                     round(100*d$PCS, 1)),
      `MFMU 무치료 (%)` = c(14.5, 14.3, 4.0, 13.6, 33.8),
      `MFMU 치료 (%)`   = c(7.1, 5.9, 1.5, 8.6, 26.9),
      출처 = rep("Landon NEJM 2009;361:1339", 5),
      check.names = FALSE)
  }, striped = TRUE, spacing = "xs")

  # ========================================================== TAB 8 ==========
  cmp <- eventReactive(input$runcmp, {
    withProgress(message = "전략 비교 중...", value = 0.1, {
      arms <- input$cmp
      out <- list()
      n <- max(1, length(arms))
      for (a in arms) {
        spec <- switch(
          a,
          none = build(drug = "none", ex = 0, cho = input$cho, lowgi = FALSE),
          mnt  = build(drug = "none", ex = 0.20, cho = 175, lowgi = TRUE),
          met  = build(drug = "met",  ex = 0.20, cho = 175, lowgi = TRUE),
          ins  = build(drug = "ins",  ex = 0.20, cho = 175, lowgi = TRUE),
          glyb = build(drug = "glyb", ex = 0.20, cho = 175, lowgi = TRUE))
        d <- simulate(spec, delta = 0.05)
        d$arm <- c(none = "무치료", mnt = "MNT+운동", met = "메트포민",
                   ins = "인슐린", glyb = "글리부라이드")[[a]]
        d$titr <- spec$titr
        out[[a]] <- d
        incProgress(1/n)
      }
      do.call(rbind, out)
    })
  })

  output$plt_cmp <- renderPlot({
    cmp() |>
      mutate(day = floor(time), frac = time - day) |>
      filter(frac > 0.20, frac < 0.29) |>
      group_by(arm, day) |> summarise(FPG = mean(GLU), .groups = "drop") |>
      ggplot(aes(day/7, FPG, colour = arm)) +
      geom_line(linewidth = 0.9) +
      geom_hline(yintercept = TARGET_FPG, linetype = 2) +
      labs(x = "임신 주수 (GW)", y = "공복혈당 (mg/dL)", colour = NULL,
           title = "같은 환자, 다른 전략") + THEME
  })

  output$tbl_cmp <- renderTable({
    td <- input$tdel*7
    cmp() |> group_by(arm) |> slice(which.min(abs(time - td))) |>
      transmute(전략 = arm,
                `공복혈당` = round(GLU, 1),
                `평균혈당` = round(MBG, 1),
                `HbA1c` = round(HBA1C, 2),
                `출생체중 (g)` = round(FWT),
                `z-score` = round(BWZ, 2),
                `태아지방 %` = round(FATPCT, 1),
                `LGA %` = round(100*PLGA, 1),
                `신생아저혈당 %` = round(100*PNHYP, 1),
                `견부난산 %` = round(100*PSD, 1),
                `태아 메트포민` = round(CMETF, 3),
                `태아 글리부라이드` = round(CGLYF, 4),
                `5년 T2DM %` = round(100*CIT2D, 1)) |>
      ungroup()
  }, striped = TRUE, spacing = "xs")

  # ========================================================== TAB 9 ==========
  output$plt_a1c <- renderPlot({
    fast() |> transmute(arm, GW, `평균혈당 MBG (mg/dL)` = MBG,
                        `HbA1c (%)` = HBA1C) |>
      pivot_longer(-c(arm, GW)) |>
      ggplot(aes(GW, value, colour = arm)) +
      geom_line(linewidth = 0.8) + facet_wrap(~name, scales = "free_y") +
      geom_vline(xintercept = input$tdel, linetype = 3) +
      labs(x = "임신 주수 (GW)", y = NULL, colour = NULL) + THEME
  })

  output$plt_bio <- renderPlot({
    fast() |> transmute(arm, GW, `유리지방산 (mmol/L)` = FFA,
                        `모체 지방량 (kg)` = FATM) |>
      pivot_longer(-c(arm, GW)) |>
      ggplot(aes(GW, value, colour = arm)) +
      geom_line(linewidth = 0.8) + facet_wrap(~name, scales = "free_y") +
      labs(x = "임신 주수 (GW)", y = NULL, colour = NULL) + THEME
  })

  output$plt_pp <- renderPlot({
    fast() |>
      transmute(arm, GW, `인슐린 감수성 SI` = SIP, `베타세포 BCM` = BCM,
                `Disposition index (rel.)` = DI/DIREF,
                `공복혈당 (mg/dL)` = GLU) |>
      pivot_longer(-c(arm, GW)) |>
      ggplot(aes(GW, value, colour = arm)) +
      geom_line(linewidth = 0.8) +
      geom_vline(xintercept = input$tdel, linetype = 2, colour = "#b03030") +
      facet_wrap(~name, scales = "free_y") +
      labs(x = "임신 주수 (GW) — 분만 이후는 산후 주수", y = NULL, colour = NULL,
           title = "분만(붉은 선)에서 태반이 사라진 뒤 남는 것",
           subtitle = "SI는 회복되지만 베타세포 결손은 회복되지 않습니다") + THEME
  })

  output$tbl_pp <- renderTable({
    td <- input$tdel*7
    allsim() |> group_by(arm) |>
      slice(which.min(abs(time - (td + 56)))) |>   # ~8 weeks post partum
      transmute(arm,
                `산후 공복혈당 (mg/dL)` = round(GLU, 1),
                `산후 인슐린 (µU/mL)` = round(INS, 1),
                `산후 SI (임신전 대비 %)` = round(100*MATSUDA, 0),
                `산후 BCM` = round(BCM, 2),
                `산후 DI (상대값)` = round(DI/DIREF, 2),
                `5년 T2DM 누적발생 (%)` = round(100*CIT2D, 1),
                `수유 반영` = ifelse(input$lact, "예", "아니오")) |>
      ungroup()
  }, striped = TRUE, spacing = "xs")

  # ========================================================= TAB 10 ==========
  hapo <- eventReactive(input$runhapo, {
    withProgress(message = "HAPO 스윕 실행 중 (BCAP x BMI)...", value = 0.1, {
      grid <- expand.grid(BCAP = seq(0.45, 1.15, 0.10), BMI = c(22, 27, 32, 37))
      res <- lapply(seq_len(nrow(grid)), function(i) {
        p <- list(BCAP = grid$BCAP[i], BMI = grid$BMI[i], TDEL = 273)
        s <- gdm |> param(p) |> data_set(meal_events(200, 56, 280)) |>
          mrgsim(start = 56, end = 280, delta = 0.25, hmax = 0.05) |>
          as.data.frame()
        f <- fasting_series(s, "GLU")
        d <- s[which.min(abs(s$time - 273)), ]
        incProgress(1/nrow(grid))
        data.frame(BCAP = grid$BCAP[i], BMI = grid$BMI[i],
                   FPG = mean(tail(f$GLU, 7)), MBG = d$MBG,
                   BW = d$FWT, BWZ = d$BWZ, FATPCT = d$FATPCT,
                   LGA = 100*d$PLGA, NEOHYPO = 100*d$PNHYP,
                   SD = 100*d$PSD, T2DM = 100*d$CIT2D)
      })
      do.call(rbind, res)
    })
  })

  output$plt_hapo <- renderPlot({
    h <- hapo()
    obs <- data.frame(FPG = c(72, 77, 81, 85, 89, 94, 102),
                      LGA = c(5.3, 7.9, 9.4, 11.5, 14.4, 19.3, 26.3))
    ggplot(h, aes(FPG, LGA)) +
      geom_point(aes(colour = factor(BMI)), size = 2.4) +
      geom_smooth(se = FALSE, method = "loess", formula = y ~ x,
                  colour = "grey30", linewidth = 0.8) +
      geom_point(data = obs, shape = 4, size = 3.6, stroke = 1.2,
                 colour = "#b03030") +
      geom_line(data = obs, colour = "#b03030", linetype = 2) +
      labs(x = "말기 공복혈당 (mg/dL)", y = "LGA 확률 (%)",
           colour = "임신 전 BMI",
           title = "HAPO 연속 기울기 — 역치는 없다",
           subtitle = "붉은 X/점선 = HAPO 관찰값 (범주 1-7); 점 = 모델") + THEME
  })

  output$tbl_hapo <- renderTable({
    hapo() |> arrange(FPG) |>
      transmute(BCAP, BMI, `공복혈당` = round(FPG, 1),
                `평균혈당` = round(MBG, 1), `출생체중` = round(BW),
                `z-score` = round(BWZ, 2), `태아지방 %` = round(FATPCT, 1),
                `LGA %` = round(LGA, 1), `신생아저혈당 %` = round(NEOHYPO, 1),
                `견부난산 %` = round(SD, 1), `5년 T2DM %` = round(T2DM, 1))
  }, striped = TRUE, spacing = "xs")
}

shinyApp(ui, server)

# =============================================================================
#  NOTE ON HONEST PRESENTATION
#  ---------------------------------------------------------------------------
#  Tab 7 shows the model's endpoint probabilities SIDE BY SIDE with the observed
#  MFMU trial rates rather than in isolation, and tab 10 overlays the actual
#  HAPO category means on the model sweep. That is deliberate: a QSP dashboard
#  that only shows its own output invites the reader to mistake calibration for
#  validation. Where the model disagrees with the trial numbers, the disagreement
#  should be visible on the same axes.
# =============================================================================
