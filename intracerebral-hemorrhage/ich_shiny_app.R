## =============================================================================
##  ich_shiny_app.R
##  Interactive dashboard for the spontaneous intracerebral haemorrhage (ICH)
##  QSP model.
##
##  Run with:
##      shiny::runApp("ich_shiny_app.R")
##  or from this directory:
##      R -e 'shiny::runApp(".")'
##
##  The app sources ich_mrgsolve_model.R for the compiled model, the scenario
##  library and the dosing helpers, so the two files can never drift apart.
##
##  EIGHT TABS
##   1. 환자 · 병변 (Patient & Lesion)      — covariates, the volume budget
##   2. 약물 PK (Drug PK)                   — every drug on one time axis
##   3. 지혈 (Haemostasis)                  — the clot-competence factor
##   4. 혈종 · 출혈 플럭스 (Bleeding)        — the product that drives dV/dt
##   5. 두개내압 · 관류 (ICP & Perfusion)    — the two appearances of MAP
##   6. 철 · 염증 (Iron & Inflammation)      — the delayed chemical channel
##   7. 임상 결과 (Clinical Endpoints)       — NIHSS, GCS, mRS, mortality
##   8. 시나리오 비교 (Scenario Comparison)  — all 12 arms + the U-shape surface
## =============================================================================

library(shiny)
library(mrgsolve)
library(dplyr)
library(tidyr)
library(ggplot2)

## ---------------------------------------------------------------------------
##  Load the model, scenarios and dosing helpers
## ---------------------------------------------------------------------------
.here <- tryCatch(dirname(sys.frame(1)$ofile), error = function(e) ".")
if (is.null(.here) || !nzchar(.here)) .here <- "."
source(file.path(.here, "ich_mrgsolve_model.R"), local = FALSE)

THEME <- theme_bw(base_size = 12) +
  theme(legend.position = "bottom",
        strip.background = element_rect(fill = "grey92"),
        plot.title = element_text(face = "bold", size = 12))

PAL <- c("#3060b0", "#c05050", "#2e8b57", "#8000a0", "#e08050",
         "#7fa8b8", "#b8845f", "#555555")

## helper: long-format selected columns
gather_vars <- function(d, vars) {
  d %>% select(time, all_of(vars)) %>%
    pivot_longer(-time, names_to = "variable", values_to = "value") %>%
    mutate(variable = factor(variable, levels = vars))
}

## ===========================================================================
##  UI
## ===========================================================================
ui <- fluidPage(
  titlePanel("자발성 뇌내출혈 (ICH) QSP 모델 — Spontaneous Intracerebral Haemorrhage"),
  tags$p(style = "color:#666; margin-top:-10px; font-size:90%;",
         "출혈 플럭스 = 구동압 x 지혈 무능 x 열린 출혈점 수. ",
         "결과 = 질량효과 채널 + 헴·철 채널. ",
         "MAP는 부호가 반대인 두 자리에 등장합니다 (CPP = MAP - ICP)."),

  sidebarLayout(
    sidebarPanel(
      width = 3,

      selectInput("scenario", "시나리오 (Scenario)",
                  choices = names(scenarios),
                  selected = "3_bp_intensive_140"),
      helpText(textOutput("scen_note")),
      tags$hr(),

      h4("환자 (Patient)"),
      sliderInput("VHEM0", "초기 혈종 용적 (mL)", 5, 90, 30, step = 1),
      sliderInput("AGE", "연령 (years)", 40, 95, 68, step = 1),
      sliderInput("WT", "체중 (kg)", 45, 130, 75, step = 1),
      selectInput("LOC", "위치 (Location)",
                  choices = c("심부 (deep)" = 0, "뇌엽 (lobar)" = 1,
                              "후두개와 (infratentorial)" = 2),
                  selected = 0),
      sliderInput("SVD", "소혈관질환 부담 (0-1)", 0, 1, 0.4, step = 0.05),
      sliderInput("FIVH", "뇌실내 확장 분율", 0, 0.6, 0.10, step = 0.02),
      sliderInput("SBPBASE", "무치료 SBP 설정점 (mmHg)", 130, 240, 185, step = 5),
      sliderInput("HPG", "합토글로빈 소거능 (Hp1-1=1.0, Hp2-2=0.6)",
                  0.4, 1.2, 1.0, step = 0.05),

      tags$hr(),
      h4("항응고 상태 (Antithrombotic on board)"),
      sliderInput("FII0", "프로트롬빈 활성 (%)  [20% ~ INR 3.0]",
                  10, 100, 100, step = 5),
      sliderInput("APIX0", "아픽사반 체내량 (mg)", 0, 8, 0, step = 0.2),
      sliderInput("DABI0", "다비가트란 체내량 (mg)", 0, 30, 0, step = 1),
      sliderInput("IASA", "아스피린 혈소판 기능 소실", 0, 0.6, 0, step = 0.05),
      sliderInput("IP2Y12", "P2Y12 억제제 기능 소실", 0, 0.6, 0, step = 0.05),

      tags$hr(),
      h4("개입 (Interventions)"),
      checkboxInput("EVDON", "뇌실외배액 (EVD)", FALSE),
      checkboxInput("SURGON", "혈종 제거술 (evacuation)", FALSE),
      conditionalPanel("input.SURGON",
        sliderInput("FEVAC", "제거 분율", 0.2, 0.95, 0.70, step = 0.05),
        sliderInput("TSURG", "수술 시각 (h)", 4, 96, 30, step = 2)),
      checkboxInput("DCON", "감압 두개절제술", FALSE),
      checkboxInput("INSON", "인슐린 프로토콜", FALSE),
      checkboxInput("COOLON", "발열 관리 / 냉각", FALSE),
      checkboxInput("VITK", "비타민 K 10 mg IV", FALSE),
      sliderInput("REHAB", "재활 강도 (0-1.5)", 0, 1.5, 1.0, step = 0.1),

      tags$hr(),
      sliderInput("tmax", "시뮬레이션 기간 (일)", 2, 90, 14, step = 1),
      actionButton("run", "재계산 (Run)", class = "btn-primary")
    ),

    mainPanel(
      width = 9,
      tabsetPanel(
        id = "tabs",

        ## ---- 1 --------------------------------------------------------
        tabPanel("1. 환자 · 병변",
          h4("두개내 용적 수지 (Monro-Kellie volume budget)"),
          p("ICP는 상태변수가 아니라 이 수지의 잔차입니다. CSF 예비량이 소진될 때까지는 곡선이 평탄하고, 소진된 뒤부터 급격해집니다."),
          plotOutput("p_budget", height = "300px"),
          fluidRow(
            column(6, h4("혈종 · 부종 · 뇌실내혈액"), plotOutput("p_vols", height = "260px")),
            column(6, h4("압력-용적 관계 (실측 궤적)"), plotOutput("p_pv", height = "260px"))
          ),
          h4("요약 (Key numbers)"), tableOutput("t_summary")
        ),

        ## ---- 2 --------------------------------------------------------
        tabPanel("2. 약물 PK",
          h4("투여된 약물의 농도 · 효과부위"),
          p("빈 패널은 그 시나리오에서 해당 약물이 투여되지 않았음을 뜻합니다."),
          plotOutput("p_pk", height = "420px"),
          fluidRow(
            column(6, h4("역전제 결합 (안덱사네트 / 이다루시주맙)"),
                      plotOutput("p_reversal", height = "280px")),
            column(6, h4("유리 항응고제 활성"),
                      plotOutput("p_anticoag", height = "280px"))
          ),
          p(strong("반동(rebound)은 별도 항으로 코딩되어 있지 않습니다."),
            " 안덱사네트 중단 후 anti-Xa 상승은 복합체 해리(KOFFA)와 말초구획 재분포만으로 창발합니다. KOFFA = 0 으로 두면 사라집니다.")
        ),

        ## ---- 3 --------------------------------------------------------
        tabPanel("3. 지혈",
          h4("클롯 경쟁력 — 출혈 플럭스의 두 번째 인자"),
          plotOutput("p_clot", height = "320px"),
          fluidRow(
            column(6, h4("응고 인자 · 혈소판"), plotOutput("p_coag", height = "280px")),
            column(6, h4("트롬빈 생성능 · 섬유소분해"), plotOutput("p_thr", height = "280px"))
          ),
          p("CLOT는 ", strong("두 팔"), "로 구성됩니다: 혈소판 팔(WPLT, VKA·DOAC에 영향받지 않음)과 응고 팔(WFIB). 이 분리 때문에 혈소판 수혈·DDAVP는 와파린 출혈에 효과가 없고, INR 3.0이 출혈을 3.5배가 아니라 약 2배로 늘립니다.")
        ),

        ## ---- 4 --------------------------------------------------------
        tabPanel("4. 혈종 · 출혈 플럭스",
          h4("dV/dt = KBLEED x NOPEN x (1 - CLOT) x max(0, MAP - PTISS)"),
          p("네 인자를 함께 표시합니다. 어느 하나라도 0이 되면 플럭스 전체가 사라집니다 — 이것이 단일 인자 개입의 효과가 작고 병용이 상승적인 이유입니다."),
          plotOutput("p_flux", height = "420px"),
          fluidRow(
            column(6, h4("혈종 용적 및 확장"), plotOutput("p_vhem", height = "280px")),
            column(6, h4("탐포네이드: 조직 대항압"), plotOutput("p_tamp", height = "280px"))
          )
        ),

        ## ---- 5 --------------------------------------------------------
        tabPanel("5. 두개내압 · 관류",
          h4("MAP의 두 얼굴"),
          p("같은 MAP가 출혈을 밀어내는 구동압이면서(낮을수록 좋음) CPP의 분자입니다(낮을수록 나쁨). U자 곡선은 어느 방정식에도 코딩되어 있지 않습니다."),
          plotOutput("p_icp", height = "380px"),
          fluidRow(
            column(6, h4("자동조절 · 상대 관류"), plotOutput("p_auto", height = "280px")),
            column(6, h4("누적 노출 (관측량)"), plotOutput("p_accum", height = "280px"))
          )
        ),

        ## ---- 6 --------------------------------------------------------
        tabPanel("6. 철 · 염증",
          h4("지연성 화학 채널: 적혈구 용해 → 헴 → HO-1 → 유리 Fe2+ → 지질과산화"),
          plotOutput("p_iron", height = "380px"),
          fluidRow(
            column(6, h4("신경염증"), plotOutput("p_inflam", height = "280px")),
            column(6, h4("조직 생존 (뉴런 · 백질)"), plotOutput("p_tissue", height = "280px"))
          ),
          p("데페록사민은 ", strong("VHEM으로 가는 간선이 하나도 없습니다"),
            " — 이 채널에만 작용합니다. 그래서 24시간 혈종 용적은 대조군과 동일하면서 90일 결과만 움직입니다 (i-DEF 패턴).")
        ),

        ## ---- 7 --------------------------------------------------------
        tabPanel("7. 임상 결과",
          h4("NIHSS · GCS · 90일 mRS · 사망 확률"),
          plotOutput("p_clin", height = "380px"),
          fluidRow(
            column(6, h4("ICH score 구성요소"), tableOutput("t_ichscore")),
            column(6, h4("결과 확률"), plotOutput("p_outcome", height = "260px"))
          ),
          p(strong("주의: "), "DNR/치료 중단에 의한 자기실현적 예언(Becker 2001)은 모든 ICH 시험 해석을 교란합니다. 모델은 이를 지도에 노드로 포함하지만 ODE로 구현하지는 않았습니다 — 모델의 사망률은 생물학적 궤적만 반영합니다.")
        ),

        ## ---- 8 --------------------------------------------------------
        tabPanel("8. 시나리오 비교",
          h4("12개 시나리오"),
          p("아래 표는 문헌 앵커와 직접 비교할 수 있도록 각 시험이 실제로 보고한 지표를 계산합니다. 계산에 30-90초 걸립니다."),
          actionButton("run_all", "12개 시나리오 모두 실행", class = "btn-warning"),
          tags$br(), tags$br(),
          tableOutput("t_all"),
          tags$hr(),
          h4("창발하는 U자 곡선 (Emergent U-shape)"),
          p("혈압을 얼마나 깊이 내릴지를 0에서 1까지 훑습니다. 혈종 확장은 단조 감소하지만 CPP<60 노출은 단조 증가하므로, 결과는 중간에서 최적이 됩니다."),
          actionButton("run_surf", "반응면 계산", class = "btn-warning"),
          plotOutput("p_surface", height = "340px"),
          tableOutput("t_surface")
        )
      )
    )
  )
)

## ===========================================================================
##  SERVER
## ===========================================================================
server <- function(input, output, session) {

  output$scen_note <- renderText({
    s <- scenarios[[input$scenario]]
    if (is.null(s)) "" else s$note
  })

  ## when a scenario is chosen, adopt its parameter overrides in the sidebar
  observeEvent(input$scenario, {
    p <- scenarios[[input$scenario]]$param
    if (!is.null(p$VHEM0))   updateSliderInput(session, "VHEM0",   value = p$VHEM0)
    if (!is.null(p$SBPBASE)) updateSliderInput(session, "SBPBASE", value = p$SBPBASE)
    if (!is.null(p$FII0))    updateSliderInput(session, "FII0",    value = p$FII0)
    if (!is.null(p$APIX0))   updateSliderInput(session, "APIX0",   value = p$APIX0)
    if (!is.null(p$FIVH))    updateSliderInput(session, "FIVH",    value = p$FIVH)
    updateCheckboxInput(session, "EVDON",  value = isTRUE(p$EVDON  == 1))
    updateCheckboxInput(session, "SURGON", value = isTRUE(p$SURGON == 1))
    updateCheckboxInput(session, "INSON",  value = isTRUE(p$INSON  == 1))
    updateCheckboxInput(session, "COOLON", value = isTRUE(p$COOLON == 1))
    updateCheckboxInput(session, "VITK",   value = isTRUE(p$VITK   == 1))
  })

  ## ---- the single simulation the first seven tabs share -----------------
  sim <- eventReactive(input$run, {
    s   <- scenarios[[input$scenario]]
    end <- input$tmax * 24

    pars <- c(
      s$param,
      list(VHEM0 = input$VHEM0, AGE = input$AGE, WT = input$WT,
           LOC = as.numeric(input$LOC), SVD = input$SVD, FIVH = input$FIVH,
           SBPBASE = input$SBPBASE, HPG = input$HPG,
           FII0 = input$FII0, APIX0 = input$APIX0, DABI0 = input$DABI0,
           IASA = input$IASA, IP2Y12 = input$IP2Y12,
           EVDON = as.numeric(input$EVDON), SURGON = as.numeric(input$SURGON),
           DCON = as.numeric(input$DCON), INSON = as.numeric(input$INSON),
           COOLON = as.numeric(input$COOLON), VITK = as.numeric(input$VITK),
           REHAB = input$REHAB)
    )
    ## sidebar wins over the scenario's own overrides for the shared covariates
    pars <- pars[!duplicated(names(pars), fromLast = TRUE)]
    if (input$SURGON) {
      pars$FEVAC <- input$FEVAC
      pars$TSURG <- input$TSURG
    }

    m <- param(ich, pars)
    withProgress(message = "적분 중 (integrating)...", value = 0.5, {
      out <- if (is.null(s$events)) mrgsim(m, end = end, delta = 0.25)
             else mrgsim(m, events = s$events, end = end, delta = 0.25)
    })
    as_tibble(out) %>% mutate(day = time / 24)
  }, ignoreNULL = FALSE)

  ## ================= TAB 1 : patient & lesion ==========================
  output$p_budget <- renderPlot({
    d <- sim()
    d %>% transmute(day,
                    `혈종 VHEM` = VHEM,
                    `조기 부종 OEDE` = OEDE,
                    `후기 부종 OEDL` = OEDL,
                    `뇌실내혈액 VIVH` = VIVH,
                    `CSF 변위 (음수 = 대상성)` = VCSF - 140) %>%
      pivot_longer(-day) %>%
      ggplot(aes(day, value, fill = name)) +
      geom_area(alpha = 0.85) +
      geom_line(data = transmute(d, day, value = VADD, name = "비대상 용적 VADD"),
                aes(day, value), inherit.aes = FALSE,
                colour = "black", linewidth = 1.1, linetype = "dashed") +
      scale_fill_manual(values = PAL) +
      labs(x = "일 (days)", y = "용적 (mL)", fill = NULL,
           subtitle = "점선 = 대상성으로 흡수되지 않은 잔여 용적 (ICP를 결정)") + THEME
  })

  output$p_vols <- renderPlot({
    gather_vars(sim(), c("VHEM", "OEDE", "OEDL", "VIVH")) %>%
      mutate(day = time / 24) %>%
      ggplot(aes(day, value, colour = variable)) +
      geom_line(linewidth = 0.8) + scale_colour_manual(values = PAL) +
      labs(x = "일", y = "mL", colour = NULL) + THEME
  })

  output$p_pv <- renderPlot({
    ggplot(sim(), aes(VADD, ICP, colour = day)) +
      geom_path(linewidth = 1) +
      scale_colour_viridis_c(option = "plasma") +
      labs(x = "비대상 용적 VADD (mL)", y = "ICP (mmHg)", colour = "일",
           subtitle = "CSF 예비량 소진 전/후로 기울기가 달라집니다") + THEME
  })

  output$t_summary <- renderTable({
    d <- sim(); e <- max(d$time)
    at <- function(v, t) d[[v]][which.min(abs(d$time - t))]
    tibble(
      `지표` = c("초기 혈종 (mL)", "24시간 혈종 (mL)", "확장 (mL)", "확장 (%)",
               "PHE 정점 (mL)", "PHE 정점 (일)", "ICP 정점 (mmHg)",
               "CPP 최저 (mmHg)", "CPP<60 누적 (h)", "ICP>20 누적 (h)",
               "정중선 변위 정점 (mm)", "NIHSS (24 h)", "NIHSS (최종)",
               "GCS 최저", "ICH score (24 h)", "P(mRS 0-2)", "P(사망)"),
      `값` = c(round(at("VHEM", 0), 1), round(at("VHEM", 24), 1),
             round(at("VHEM", 24) - at("VHEM", 0), 2),
             round(100 * (at("VHEM", 24) - at("VHEM", 0)) / at("VHEM", 0), 1),
             round(max(d$PHE), 1), round(d$time[which.max(d$PHE)] / 24, 2),
             round(max(d$ICP), 1), round(min(d$CPP), 1),
             round(at("TCPP", e), 1), round(at("TICP", e), 1),
             round(max(d$MIDLINE), 1),
             round(at("NIH", 24), 1), round(at("NIH", e), 1),
             round(min(d$GCSE), 1), round(at("ICHSC", 24), 0),
             round(at("PMRS02", e), 3), round(at("PMORT", e), 3))
    )
  }, digits = 3)

  ## ================= TAB 2 : PK ========================================
  output$p_pk <- renderPlot({
    d <- sim()
    d %>% transmute(day,
                    `니카르디핀 중앙 (mg)` = NICA1,
                    `니카르디핀 효과부위` = NICE,
                    `라베탈롤 (mg)` = LABE,
                    `TXA 중앙 (mg)` = TXAC,
                    `4F-PCC 잔여 (IU)` = PCCA,
                    `데페록사민 (mg)` = DFOC,
                    `만니톨 (g)` = MANN,
                    `뇌실내 알테플라제 (mg)` = IVTPA) %>%
      pivot_longer(-day) %>%
      ggplot(aes(day, value)) + geom_line(colour = PAL[1], linewidth = 0.8) +
      facet_wrap(~name, scales = "free_y", ncol = 4) +
      labs(x = "일", y = NULL) + THEME
  })

  output$p_reversal <- renderPlot({
    sim() %>% transmute(day,
                        `유리 안덱사네트 (NEQ)` = ANDXE,
                        `안덱사네트-아픽사반 복합체` = CPLXA,
                        `유리 이다루시주맙 (NEQ)` = IDAE,
                        `이다-다비 복합체` = CPLXD) %>%
      pivot_longer(-day) %>%
      ggplot(aes(day, value, colour = name)) + geom_line(linewidth = 0.8) +
      scale_colour_manual(values = PAL) +
      labs(x = "일", y = "mg-당량", colour = NULL) + THEME
  })

  output$p_anticoag <- renderPlot({
    sim() %>% transmute(day, `유리 아픽사반 anti-Xa (ng/mL)` = AXA,
                        `유리 다비가트란 (ng/mL)` = DTT, INR = INR) %>%
      pivot_longer(-day) %>%
      ggplot(aes(day, value)) + geom_line(colour = PAL[2], linewidth = 0.9) +
      facet_wrap(~name, scales = "free_y", ncol = 1) +
      labs(x = "일", y = NULL) + THEME
  })

  ## ================= TAB 3 : haemostasis ===============================
  output$p_clot <- renderPlot({
    sim() %>% transmute(day, `클롯 경쟁력 CLOT` = CLOT,
                        `열린 출혈점 NOPEN` = NOPEN,
                        `트롬빈 생성능 THRGEN` = THRGEN) %>%
      pivot_longer(-day) %>%
      ggplot(aes(day, value, colour = name)) + geom_line(linewidth = 1) +
      scale_colour_manual(values = PAL) +
      labs(x = "일", y = "상대값", colour = NULL) + THEME
  })

  output$p_coag <- renderPlot({
    sim() %>% transmute(day, `프로트롬빈 활성 (%)` = FII,
                        `피브리노겐 (mg/dL)` = FIB,
                        `혈소판 기능 (0-1)` = PLTF) %>%
      pivot_longer(-day) %>%
      ggplot(aes(day, value)) + geom_line(colour = PAL[3], linewidth = 0.9) +
      facet_wrap(~name, scales = "free_y", ncol = 1) +
      labs(x = "일", y = NULL) + THEME
  })

  output$p_thr <- renderPlot({
    sim() %>% transmute(day, `플라스민 활성` = PLAS,
                        `조직 트롬빈 THRT` = THRT,
                        `INR` = INR) %>%
      pivot_longer(-day) %>%
      ggplot(aes(day, value)) + geom_line(colour = PAL[4], linewidth = 0.9) +
      facet_wrap(~name, scales = "free_y", ncol = 1) +
      labs(x = "일", y = NULL) + THEME
  })

  ## ================= TAB 4 : bleeding flux =============================
  output$p_flux <- renderPlot({
    sim() %>% transmute(day,
                        `출혈 플럭스 (mL/h)` = BLEEDR,
                        `구동압 PDRIVE (mmHg)` = PDRIVE,
                        `지혈 무능 (1-CLOT)` = 1 - CLOT,
                        `열린 출혈점 NOPEN` = NOPEN) %>%
      pivot_longer(-day) %>%
      ggplot(aes(day, value)) +
      geom_line(colour = PAL[2], linewidth = 0.9) +
      facet_wrap(~name, scales = "free_y", ncol = 2) +
      labs(x = "일", y = NULL,
           subtitle = "위 왼쪽 = 세 인자의 곱") + THEME
  })

  output$p_vhem <- renderPlot({
    d <- sim(); v0 <- d$VHEM[1]
    ggplot(d, aes(day, VHEM)) +
      geom_line(colour = PAL[2], linewidth = 1) +
      geom_hline(yintercept = v0, linetype = "dotted") +
      geom_hline(yintercept = v0 * 1.33, linetype = "dashed", colour = "grey40") +
      annotate("text", x = max(d$day) * 0.7, y = v0 * 1.36,
               label = "확장 정의 (+33%)", size = 3.2, colour = "grey30") +
      labs(x = "일", y = "혈종 용적 (mL)") + THEME
  })

  output$p_tamp <- renderPlot({
    sim() %>% transmute(day, `MAP` = MAP, `조직 대항압 PTISS` = PTISS,
                        `ICP` = ICP) %>%
      pivot_longer(-day) %>%
      ggplot(aes(day, value, colour = name)) + geom_line(linewidth = 0.9) +
      scale_colour_manual(values = PAL) +
      labs(x = "일", y = "mmHg", colour = NULL,
           subtitle = "MAP와 PTISS가 만나면 출혈이 스스로 멈춥니다") + THEME
  })

  ## ================= TAB 5 : ICP & perfusion ===========================
  output$p_icp <- renderPlot({
    sim() %>% transmute(day, `SBP` = SBP, `MAP` = MAP, `ICP` = ICP,
                        `CPP` = CPP) %>%
      pivot_longer(-day) %>%
      ggplot(aes(day, value, colour = name)) + geom_line(linewidth = 0.9) +
      geom_hline(yintercept = 60, linetype = "dashed", colour = "grey50") +
      geom_hline(yintercept = 20, linetype = "dotted", colour = "grey50") +
      scale_colour_manual(values = PAL) +
      labs(x = "일", y = "mmHg", colour = NULL,
           subtitle = "파선 = CPP 60 임계, 점선 = ICP 20 임계") + THEME
  })

  output$p_auto <- renderPlot({
    sim() %>% transmute(day, `자동조절 능력 AUTOR` = AUTOR,
                        `상대 관류 CBFP` = CBFP,
                        `허혈 스트레스 ISCH` = ISCH) %>%
      pivot_longer(-day) %>%
      ggplot(aes(day, value, colour = name)) + geom_line(linewidth = 0.9) +
      scale_colour_manual(values = PAL) +
      labs(x = "일", y = "상대값", colour = NULL) + THEME
  })

  output$p_accum <- renderPlot({
    sim() %>% transmute(day, `SBP>140 AUC` = ASBP, `CPP<60 시간` = TCPP,
                        `ICP>20 시간` = TICP, `Fe2+ AUC` = AFE,
                        `IL-6 AUC` = AIL6) %>%
      pivot_longer(-day) %>%
      ggplot(aes(day, value)) + geom_line(colour = PAL[6], linewidth = 0.9) +
      facet_wrap(~name, scales = "free_y", ncol = 3) +
      labs(x = "일", y = NULL,
           subtitle = "이들은 파라미터가 아니라 적분 결과(관측량)입니다") + THEME
  })

  ## ================= TAB 6 : iron & inflammation =======================
  output$p_iron <- renderPlot({
    sim() %>% transmute(day, `유리 헴 HEME` = HEME, `유리 Fe2+ FEII` = FEII,
                        `페리틴 격리능 FERR` = FERR,
                        `지질과산화 LPO` = LPO,
                        `페록사민 (배설)` = FERX) %>%
      pivot_longer(-day) %>%
      ggplot(aes(day, value, colour = name)) + geom_line(linewidth = 0.9) +
      scale_colour_manual(values = PAL) +
      labs(x = "일", y = "상대값", colour = NULL) + THEME
  })

  output$p_inflam <- renderPlot({
    sim() %>% transmute(day, `전염증 미세아교세포 MG1` = MG1,
                        `수복형 MG2` = MG2, `호중구 NEU` = NEU,
                        `IL-6` = IL6, `IL-10` = IL10,
                        `MMP-9` = MMP9, `BBB 투과성` = BBBP) %>%
      pivot_longer(-day) %>%
      ggplot(aes(day, value, colour = name)) + geom_line(linewidth = 0.8) +
      scale_colour_manual(values = c(PAL, "#999999")) +
      labs(x = "일", y = "상대값", colour = NULL) + THEME
  })

  output$p_tissue <- renderPlot({
    sim() %>% transmute(day, `생존 뉴런 분율 NEUR` = NEUR,
                        `백질 보전 WMI` = WMI,
                        `가소성 예비 PLAST` = PLAST,
                        `기계적 변형 STRAIN` = STRAIN) %>%
      pivot_longer(-day) %>%
      ggplot(aes(day, value, colour = name)) + geom_line(linewidth = 0.9) +
      scale_colour_manual(values = PAL) +
      labs(x = "일", y = "0-1 (STRAIN은 지수)", colour = NULL) + THEME
  })

  ## ================= TAB 7 : clinical ==================================
  output$p_clin <- renderPlot({
    sim() %>% transmute(day, `NIHSS` = NIH, `추정 GCS` = GCSE,
                        `체온 (degC)` = TEMP, `혈당 (mmol/L)` = GLU,
                        `정중선 변위 (mm)` = MIDLINE,
                        `ICH score` = ICHSC) %>%
      pivot_longer(-day) %>%
      ggplot(aes(day, value)) + geom_line(colour = PAL[1], linewidth = 0.9) +
      facet_wrap(~name, scales = "free_y", ncol = 3) +
      labs(x = "일", y = NULL) + THEME
  })

  output$t_ichscore <- renderTable({
    d <- sim(); at <- function(v, t) d[[v]][which.min(abs(d$time - t))]
    g <- at("GCSE", 24); v <- at("VHEM", 24); iv <- at("VIVH", 24)
    tibble(
      `구성요소` = c("GCS", "용적 >= 30 mL", "뇌실내출혈", "후두개와", "연령 >= 80", "합계"),
      `값` = c(sprintf("%.1f", g), sprintf("%.1f mL", v), sprintf("%.1f mL", iv),
             ifelse(as.numeric(input$LOC) > 1.5, "예", "아니오"),
             sprintf("%d세", input$AGE), ""),
      `점수` = c(ifelse(g <= 4, 2, ifelse(g <= 12, 1, 0)),
               ifelse(v >= 30, 1, 0), ifelse(iv > 0.5, 1, 0),
               ifelse(as.numeric(input$LOC) > 1.5, 1, 0),
               ifelse(input$AGE >= 80, 1, 0),
               at("ICHSC", 24))
    )
  })

  output$p_outcome <- renderPlot({
    sim() %>% transmute(day, `P(mRS 0-2)` = PMRS02, `P(사망)` = PMORT,
                        `효용 가중 mRS` = UWMRS) %>%
      pivot_longer(-day) %>%
      ggplot(aes(day, value, colour = name)) + geom_line(linewidth = 1) +
      scale_colour_manual(values = PAL) + ylim(0, 1) +
      labs(x = "일", y = "확률 / 효용", colour = NULL) + THEME
  })

  ## ================= TAB 8 : scenario comparison =======================
  all_sim <- eventReactive(input$run_all, {
    withProgress(message = "12개 시나리오 실행 중...", value = 0.3, {
      summarise_scenarios()
    })
  })

  output$t_all <- renderTable({
    req(all_sim())
    all_sim() %>%
      transmute(`시나리오` = scenario,
                `확장 (mL)` = expansion_mL,
                `확장 (%)` = expansion_pct,
                `SBP 최저` = SBP_nadir_mmHg,
                `PHE 정점 (mL)` = PHE_peak_mL,
                `ICP 정점` = ICP_peak_mmHg,
                `CPP<60 (h)` = hours_CPP_lt60,
                `Fe AUC (14d)` = AUC_Fe2_14d,
                `NIHSS 90d` = NIHSS_90d,
                `P(mRS 0-2)` = P_mRS02_90d,
                `P(사망)` = P_death)
  }, digits = 3)

  surf <- eventReactive(input$run_surf, {
    withProgress(message = "반응면 계산 중...", value = 0.3, {
      bp_response_surface(speeds = c(1, 4), end = 1440)
    })
  })

  output$p_surface <- renderPlot({
    req(surf())
    s <- surf()
    best <- s %>% group_by(ttt) %>% slice_max(P_mRS02, n = 1)
    ggplot(s, aes(SBP_nadir, P_mRS02, colour = factor(ttt))) +
      geom_line(linewidth = 1) + geom_point(size = 2) +
      geom_point(data = best, size = 5, shape = 21, stroke = 1.4,
                 fill = NA, colour = "black") +
      scale_x_reverse() +
      scale_colour_manual(values = PAL[c(1, 2)]) +
      labs(x = "달성 SBP 최저치 (mmHg) — 오른쪽으로 갈수록 깊은 강압",
           y = "P(mRS 0-2) at 90 d", colour = "목표 도달 시간 (h)",
           subtitle = "검은 원 = 최적점. U자 곡선은 코딩된 것이 아니라 창발한 것입니다.") +
      THEME
  })

  output$t_surface <- renderTable({
    req(surf())
    surf() %>% filter(ttt == 1) %>%
      transmute(`강압 강도` = intensity, `SBP 최저` = SBP_nadir,
                `CPP 최저` = CPP_min, `확장 (mL)` = expansion,
                `CPP<60 (h)` = hrs_CPP60, `뉴런 생존` = NEUR_end,
                `NIHSS` = NIHSS_end, `P(mRS 0-2)` = P_mRS02)
  }, digits = 3)
}

shinyApp(ui, server)
