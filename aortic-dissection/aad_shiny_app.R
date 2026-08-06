# =============================================================================
#  aad_shiny_app.R
#  Acute Type B Aortic Dissection — QSP simulator dashboard
#  급성 B형 대동맥 박리 — QSP 시뮬레이터 대시보드
# =============================================================================
#
#  The app exists to make one comparison hard to avoid: set any two arms to the
#  SAME systolic pressure and look at what else differs.  Tab 3 does nothing but
#  that.  Every other tab is there to show which mechanism carried the
#  difference.
#
#  Tabs
#    1  환자·병변 프로파일        Patient & lesion profile
#    2  약동학·수용체 점유율      Pharmacokinetics & receptor occupancy
#    3  압력 vs 충격량 (핵심)     Pressure versus impulse — THE dissociation
#    4  두 내강 분압기            The two-lumen divider
#    5  벽 응력과 파열 위험       Wall stress & rupture hazard
#    6  가강 혈전과 역설          False-lumen thrombus & the paradox
#    7  분지 관류와 레닌 함정     Branch perfusion & the renin trap
#    8  척수와 측부순환           Spinal cord & the collateral network
#    9  시나리오 비교             Scenario comparison
#   10  임상 종말점·바이오마커    Clinical endpoints & biomarkers
#
#  Run with:
#    shiny::runApp("aortic-dissection/aad_shiny_app.R")
# =============================================================================

library(shiny)
library(mrgsolve)
library(dplyr)
library(tidyr)
library(ggplot2)

# -----------------------------------------------------------------------------
#  model
# -----------------------------------------------------------------------------
MODEL_DIR <- if (dir.exists("aortic-dissection")) "aortic-dissection" else "."
mod <- mread_cache("aad_mrgsolve_model", MODEL_DIR)

THEME <- theme_minimal(base_size = 12) +
  theme(panel.grid.minor = element_blank(),
        strip.text = element_text(face = "bold"),
        legend.position = "bottom")

# The 22 arms of the reference model, so the app and the paper agree.
SCENARIOS <- list(
  "S01 미치료 (untreated)"                        = list(),
  "S02 니카르디핀 단독 (15 mg/h)"                 = list(VD_AGENT = 1, VD_MAX = 15),
  "S03 에스몰롤 단독 (혈압 적정)"                 = list(BB_AGENT = 1, BB_MAX = 1350, BB_ON_SBP = 1),
  "S04 지침 순서: 에스몰롤 → 니카르디핀"          = list(BB_AGENT = 1, BB_MAX = 900, VD_AGENT = 1, VD_MAX = 15),
  "S05 니트로프루시드 단독"                       = list(VD_AGENT = 3, VD_MAX = 13.5),
  "S06 니트로프루시드 + 에스몰롤 (Wheat)"         = list(VD_AGENT = 3, VD_MAX = 13.5, BB_AGENT = 1, BB_MAX = 900),
  "S07 라베탈롤 단독"                             = list(BB_AGENT = 2, BB_MAX = 120, BB_ON_SBP = 1),
  "S08 클레비디핀 + 에스몰롤"                     = list(VD_AGENT = 2, VD_MAX = 21, BB_AGENT = 1, BB_MAX = 900),
  "S09 진통만 (펜타닐 75 µg/h)"                   = list(R_FEN = 0.075),
  "S10 완전 프로토콜 (진통+β+CCB)"                = list(R_FEN = 0.075, BB_AGENT = 1, BB_MAX = 900,
                                                        VD_AGENT = 1, VD_MAX = 15),
  "S21 니카르디핀 초과용량 (압력 일치)"           = list(VD_AGENT = 1, VD_MAX = 30),
  "S11 공격적 목표 SBP 95 + 신동맥 관류장애"      = list(VD_AGENT = 1, VD_MAX = 32, BB_AGENT = 1, BB_MAX = 900,
                                                        SBP_SET = 95, OBST_REN = 0.45, INV_REN = 1, INV_MES = 1),
  "S12 신동맥 관류장애, 표준 목표"                = list(VD_AGENT = 1, VD_MAX = 32, BB_AGENT = 1, BB_MAX = 900,
                                                        OBST_REN = 0.45, INV_REN = 1, INV_MES = 1),
  "S13 고용량 SNP + 신장애 (시안화물)"            = list(VD_AGENT = 3, VD_MAX = 45, BB_AGENT = 1, BB_MAX = 900,
                                                        OBST_REN = 0.45, INV_REN = 1, SBP_SET = 100),
  "S14 만성 5년: 가강 개통"                       = list(GRE0 = 2.40, D_MET = 100, D_LOS = 50, .end = 43800),
  "S15 만성 5년: 부분 혈전화 (맹관)"              = list(GRE0 = 0.05, D_MET = 100, D_LOS = 50, .end = 43800),
  "S16 만성 5년: 완전 혈전화 (작은 입구)"         = list(GRE0 = 0.05, AEN0 = 9, D_MET = 100, D_LOS = 50, .end = 43800),
  "S17 만성 5년: TEVAR (14일)"                    = list(GRE0 = 0.05, TEVAR_T = 336, D_MET = 100, D_LOS = 50, .end = 43800),
  "S22 만성 5년: TEVAR + 관류압 증강"             = list(GRE0 = 0.05, TEVAR_T = 336, D_MET = 50, D_LOS = 25,
                                                        SBP_SET = 135, .end = 43800),
  "S18 마르판형 + 로사르탄 (5년)"                 = list(GEN_MMP = 2.2, GRE0 = 0.05, D_MET = 100, D_LOS = 50, .end = 43800),
  "S19 마르판형, β차단만 (5년)"                   = list(GEN_MMP = 2.2, GRE0 = 0.05, D_MET = 100, D_LOS = 0, .end = 43800),
  "S20 만성 5년: 순응도 25%"                      = list(GRE0 = 0.05, D_MET = 100, D_LOS = 50, ADHERE = 0.25, .end = 43800)
)

simulate_arm <- function(pars, end = NULL) {
  e <- if (!is.null(pars$.end)) pars$.end else (end %||% 72)
  pars$.end <- NULL
  d <- if (e > 800) 24 else 0.1
  mod %>% param(pars) %>% mrgsim(end = e, delta = d) %>% as_tibble()
}
`%||%` <- function(a, b) if (is.null(a)) b else a

# -----------------------------------------------------------------------------
#  UI
# -----------------------------------------------------------------------------
ui <- fluidPage(
  titlePanel("급성 B형 대동맥 박리 QSP 시뮬레이터 — Acute Type B Aortic Dissection"),
  p(tags$b("한 개의 곱:"),
    "ASI = (가강 수축기 응력 / 기준) × (dP/dt max / 기준).",
    "압력은 첫 항에만, 교감신경은 둘째 항에만 들어간다.",
    "혈관확장제는 첫 항을 내리고 둘째 항을 올린다."),
  sidebarLayout(
    sidebarPanel(
      width = 3,
      h4("환자와 병변"),
      sliderInput("AEN0", "입구 파열구 면적 A_en (mm²)", 4, 120, 42, step = 1),
      sliderInput("GRE0", "재유입 전도도 G_re (0 = 맹관)", 0, 3, 1, step = 0.05),
      sliderInput("DFL0", "가강 직경 (mm)", 8, 40, 17, step = 1),
      sliderInput("DTL0", "진강 직경 (mm)", 10, 32, 20, step = 1),
      sliderInput("ELN", "탄성판 완전성 (0-1)", 0.3, 1, 0.82, step = 0.01),
      sliderInput("GEN_MMP", "유전적 MMP-9 배수 (마르판 2.2)", 1, 3, 1, step = 0.1),
      sliderInput("BARO_REF", "압수용체 작동점 (mmHg)", 85, 125, 106, step = 1),
      checkboxInput("INV_REN", "신동맥이 박리 구역에서 기시", FALSE),
      sliderInput("OBST_REN", "정적 신동맥 폐색 (0-1)", 0, 0.8, 0, step = 0.05),

      hr(), h4("급성기 약물"),
      selectInput("BB", "β차단제",
                  c("없음" = "0", "에스몰롤" = "1", "라베탈롤" = "2"), "1"),
      sliderInput("BB_MAX", "β차단제 최대 주입 (mg/h)", 0, 1350, 900, step = 25),
      checkboxInput("BB_ON_SBP", "β차단제를 혈압으로 적정 (기본은 심박수)", FALSE),
      selectInput("VD", "혈관확장제",
                  c("없음" = "0", "니카르디핀" = "1", "클레비디핀" = "2",
                    "니트로프루시드" = "3"), "1"),
      sliderInput("VD_MAX", "혈관확장제 최대 주입 (mg/h)", 0, 45, 15, step = 1),
      sliderInput("R_FEN", "펜타닐 (mg/h)", 0, 0.2, 0, step = 0.005),
      sliderInput("SBP_SET", "수축기 목표 (mmHg)", 85, 140, 110, step = 1),

      hr(), h4("만성기·중재"),
      sliderInput("D_MET", "메토프롤롤 (mg / 12 h)", 0, 200, 0, step = 25),
      sliderInput("D_LOS", "로사르탄 (mg / 12 h)", 0, 100, 0, step = 25),
      sliderInput("ADHERE", "복약 순응도", 0, 1, 1, step = 0.05),
      numericInput("TEVAR_T", "TEVAR 시점 (h, −1 = 안 함)", -1, min = -1, max = 8760),
      sliderInput("TEVAR_COVER", "분절동맥 피복 분율", 0, 0.7, 0.30, step = 0.05),
      sliderInput("END", "시뮬레이션 기간 (h)", 24, 43800, 72, step = 24),
      actionButton("go", "시뮬레이션 실행", class = "btn-primary")
    ),

    mainPanel(
      width = 9,
      tabsetPanel(
        id = "tabs",

        # ---- 1 -------------------------------------------------------------
        tabPanel(
          "1 · 환자·병변 프로파일",
          h4("병변 형상과 그 형상이 이미 정해 놓은 것"),
          tableOutput("profile"),
          p(tags$b("읽는 법:"), "가강 외벽은 중막 외측 + 외막뿐이므로 두께가 전체 벽의 40%",
            "(1.90 mm → 0.76 mm)다. 같은 압력에서 가강의 원주응력은 진강의 약 2.5배이며,",
            "이 비대칭은 약물이 아니라 해부학이 결정한다."),
          plotOutput("p_geom", height = "330px"),
          plotOutput("p_stressbar", height = "300px")
        ),

        # ---- 2 -------------------------------------------------------------
        tabPanel(
          "2 · 약동학·수용체 점유율",
          fluidRow(column(6, plotOutput("p_pk", height = "340px")),
                   column(6, plotOutput("p_occ", height = "340px"))),
          plotOutput("p_rate", height = "290px"),
          p(tags$b("적정은 적분 제어기다."), "주입 속도는 목표와의 오차를 적분해 결정되며,",
            "포화되면 표에 그대로 드러난다 — 니카르디핀은 표시 상한 15 mg/h에서",
            "수축기 121 mmHg를 넘어 내려가지 못한다."),
          tableOutput("t_exposure")
        ),

        # ---- 3 -------------------------------------------------------------
        tabPanel(
          "3 · 압력 vs 충격량 (핵심)",
          h4("같은 압력에서 무엇이 다른가"),
          fluidRow(column(4, selectInput("cmpA", "비교 A", names(SCENARIOS),
                                         "S21 니카르디핀 초과용량 (압력 일치)")),
                   column(4, selectInput("cmpB", "비교 B", names(SCENARIOS),
                                         "S04 지침 순서: 에스몰롤 → 니카르디핀")),
                   column(4, selectInput("cmpC", "기준 C", names(SCENARIOS),
                                         "S01 미치료 (untreated)"))),
          plotOutput("p_dissoc", height = "430px"),
          tableOutput("t_dissoc"),
          p(tags$b("모델이 말하는 것:"), "표시 상한의 니카르디핀은 수축기압을 33 mmHg 내리고",
            "dP/dt max는 0.3% 바꾼다. 같은 압력에서 판막 응력-충격량은 2.8배,",
            "72시간 파열구 성장은 14배, 72시간 사망률은 9배 차이가 난다.",
            "혈관확장제 단독군의 사망률은 무치료군보다 높다.")
        ),

        # ---- 4 -------------------------------------------------------------
        tabPanel(
          "4 · 두 내강 분압기",
          h4("전도도 두 개가 표현형 세 개를 만든다"),
          plotOutput("p_divider", height = "380px"),
          plotOutput("p_gre_sweep", height = "330px"),
          p(tags$b("맹관의 산수:"), "G_re → 0이면 유출이 없으므로 압력강하도 없고",
            "(가강 평균압 = 전신 평균압), 맥압은 τ_FL로 저역통과 필터되어 사라지므로",
            "가강 확장기압이 전신 평균압까지 올라간다 — 진강 확장기압보다 높다.",
            "원위 정체구역이 응고해 재유입구를 닫는 동안 입구는 그대로 열려 있다.")
        ),

        # ---- 5 -------------------------------------------------------------
        tabPanel(
          "5 · 벽 응력과 파열 위험",
          plotOutput("p_sigma", height = "350px"),
          fluidRow(column(6, plotOutput("p_hazard", height = "320px")),
                   column(6, plotOutput("p_wall", height = "320px"))),
          p(tags$b("파열 위험은 σ/S의 지수함수다:"),
            "h = 1.10e−4 · exp((σ/S − 1)/0.155) /h.",
            "σ/S = 0.40에서 연 2%, 0.67에서 연 11%로 보정했다.")
        ),

        # ---- 6 -------------------------------------------------------------
        tabPanel(
          "6 · 가강 혈전과 역설",
          h4("혈전은 벽을 선형으로 지지하고 지수로 얇게 만든다"),
          plotOutput("p_paradox", height = "360px"),
          p(tags$b("해석해:"), "성장 ∝ (1 + A·THR)^n · (1 − s·THR)에서",
            "dg/dTHR = 0 → THR* = (nA − s)/((1+n)sA) = 0.426.",
            "수치 최적값 0.425. 최악의 가강은 절반만 응고한 가강이다."),
          plotOutput("p_thr_time", height = "330px")
        ),

        # ---- 7 -------------------------------------------------------------
        tabPanel(
          "7 · 분지 관류와 레닌 함정",
          plotOutput("p_perf", height = "350px"),
          plotOutput("p_raas", height = "330px"),
          p(tags$b("함정:"), "혈압을 내리면 신장이 더 허혈해지고, 레닌이 안지오텐신 II로",
            "저항을 되돌린다. 같은 처방·같은 목표에서 SBP 141 대 113이 되며",
            "혈관확장제는 이미 상한에 있다. 해법은 4번째 혈관확장제가 아니라 신장을 여는 것이다."),
          tableOutput("t_cn")
        ),

        # ---- 8 -------------------------------------------------------------
        tabPanel(
          "8 · 척수와 측부순환",
          plotOutput("p_cord", height = "350px"),
          plotOutput("p_collat", height = "320px"),
          p(tags$b("측부순환망은 시계다:"), "분절동맥 손실은 즉시 일어나지만 측부순환 모집은",
            "약 5일이 걸린다. 하지마비 위험은 그 창(window) 안에 있고, 단계적 피복과",
            "관류압 증강은 그 창을 겨냥한 처치다 — 퇴원 처방이 아니다.")
        ),

        # ---- 9 -------------------------------------------------------------
        tabPanel(
          "9 · 시나리오 비교",
          checkboxGroupInput("multi", "비교할 시나리오",
                             choices = names(SCENARIOS),
                             selected = c("S01 미치료 (untreated)",
                                          "S02 니카르디핀 단독 (15 mg/h)",
                                          "S04 지침 순서: 에스몰롤 → 니카르디핀",
                                          "S07 라베탈롤 단독",
                                          "S21 니카르디핀 초과용량 (압력 일치)"),
                             inline = TRUE),
          plotOutput("p_multi", height = "460px"),
          tableOutput("t_multi")
        ),

        # ---- 10 ------------------------------------------------------------
        tabPanel(
          "10 · 임상 종말점·바이오마커",
          plotOutput("p_end", height = "380px"),
          fluidRow(column(6, plotOutput("p_bio", height = "320px")),
                   column(6, plotOutput("p_growth", height = "320px"))),
          tableOutput("t_end"),
          p(tags$b("주의:"), "이 모델은 교육·가설생성용이며 환자 진료에 쓰도록 검증되지 않았다.",
            "로사르탄의 무효과 예측과 관류압 증강의 순효과는 모델의 가장 취약한 예측이며",
            "aad_mrgsolve_model.R의 CALIBRATION NOTES에 그대로 표시해 두었다.")
        )
      )
    )
  )
)

# -----------------------------------------------------------------------------
#  server
# -----------------------------------------------------------------------------
server <- function(input, output, session) {

  user_pars <- reactive({
    list(AEN0 = input$AEN0, AEN_REF = 42, GRE0 = input$GRE0,
         DFL0 = input$DFL0, DTL0 = input$DTL0, GEN_MMP = input$GEN_MMP,
         BARO_REF = input$BARO_REF,
         INV_REN = as.numeric(input$INV_REN), OBST_REN = input$OBST_REN,
         BB_AGENT = as.numeric(input$BB), BB_MAX = input$BB_MAX,
         BB_ON_SBP = as.numeric(input$BB_ON_SBP),
         VD_AGENT = as.numeric(input$VD), VD_MAX = input$VD_MAX,
         R_FEN = input$R_FEN, SBP_SET = input$SBP_SET,
         D_MET = input$D_MET, D_LOS = input$D_LOS, ADHERE = input$ADHERE,
         TEVAR_T = input$TEVAR_T, TEVAR_COVER = input$TEVAR_COVER)
  })

  sim <- eventReactive(input$go, {
    p <- user_pars()
    d <- if (input$END > 800) 24 else 0.1
    m <- mod %>% param(p)
    m <- m %>% init(ELN = input$ELN, DFL = input$DFL0, DTL = input$DTL0,
                    AEN = input$AEN0)
    m %>% mrgsim(end = input$END, delta = d) %>% as_tibble()
  }, ignoreNULL = FALSE)

  tail_mean <- function(d, col, hours = 6) {
    tt <- max(d$time) - hours
    mean(d[[col]][d$time >= tt], na.rm = TRUE)
  }

  # ---- 1 -------------------------------------------------------------------
  output$profile <- renderTable({
    d <- sim()
    h_fl <- 1.90 * 0.40
    data.frame(
      항목 = c("가강 외벽 두께 h_FL (mm)", "진강 벽 두께 (mm)",
               "두께 비 (진강/가강)", "제시 시 총 직경 D_ao (mm)",
               "가강 분율", "가강 수축기 응력 σ_FL (kPa)",
               "진강 수축기 응력 σ_TL (kPa)", "응력 비 (가강/진강)",
               "τ_FL (s)", "가강 확장기압 (mmHg)", "ASI (제시 시)"),
      값 = c(sprintf("%.2f", h_fl), sprintf("%.2f", 1.90 * 0.60),
             sprintf("%.2f", 0.60 / 0.40),
             sprintf("%.1f", d$D_AORTA[1]),
             sprintf("%.2f", d$FL_FRACTION[1]),
             sprintf("%.1f", d$SIGMA_FL_SYS[1]),
             sprintf("%.1f", d$SIGMA_TL[1]),
             sprintf("%.2f", d$SIGMA_FL_SYS[1] / d$SIGMA_TL[1]),
             sprintf("%.3f", d$TAU_FL[1]),
             sprintf("%.1f", d$PFL_DIASTOLIC[1]),
             sprintf("%.2f", d$ASI[1])))
  }, striped = TRUE)

  output$p_geom <- renderPlot({
    d <- sim()
    d %>% select(time, `총 직경 D_ao` = D_AORTA, `가강 D_FL` = DFL,
                 `진강 D_TL` = DTL) %>%
      pivot_longer(-time) %>%
      ggplot(aes(time, value, colour = name)) +
      geom_line(linewidth = 0.9) +
      geom_hline(yintercept = 55, linetype = "dashed", colour = "grey40") +
      annotate("text", x = 0, y = 56.5, hjust = 0, size = 3.2,
               label = "중재 역치 55 mm") +
      labs(x = "시간 (h)", y = "직경 (mm)", colour = NULL,
           title = "대동맥 형상의 시간 경과") + THEME
  })

  output$p_stressbar <- renderPlot({
    d <- sim()
    tibble(구획 = c("가강 (h 0.76 mm)", "진강 (h 1.14 mm)", "내막 판막"),
           응력 = c(tail_mean(d, "SIGMA_FL_SYS"), tail_mean(d, "SIGMA_TL"),
                    tail_mean(d, "SIGMA_FLAP"))) %>%
      ggplot(aes(reorder(구획, 응력), 응력, fill = 구획)) +
      geom_col(width = 0.6, show.legend = FALSE) +
      geom_hline(yintercept = 118, linetype = "dashed") +
      annotate("text", x = 0.6, y = 126, hjust = 0, size = 3.2,
               label = "정상 대동맥 118 kPa") +
      coord_flip() +
      labs(x = NULL, y = "원주 응력 (kPa)",
           title = "같은 압력, 다른 두께 — 라플라스가 만드는 비대칭") + THEME
  })

  # ---- 2 -------------------------------------------------------------------
  output$p_pk <- renderPlot({
    d <- sim()
    d %>% select(time, 에스몰롤 = A_ESM, 니카르디핀 = A_NICC,
                 니트로프루시드 = A_SNP, 라베탈롤 = A_LABC,
                 클레비디핀 = A_CLV) %>%
      pivot_longer(-time) %>% filter(value > 1e-9) %>%
      ggplot(aes(time, value, colour = name)) + geom_line(linewidth = 0.9) +
      scale_y_log10() +
      labs(x = "시간 (h)", y = "체내 약물량 (mg, log)", colour = NULL,
           title = "약동학") + THEME
  })

  output$p_occ <- renderPlot({
    d <- sim()
    d %>% select(time, `β1 점유` = OCC_BETA1, `DHP 점유` = OCC_CCB,
                 `NO 공여` = OCC_NOD, `AT1 차단` = OCC_AT1) %>%
      pivot_longer(-time) %>%
      ggplot(aes(time, value, colour = name)) + geom_line(linewidth = 0.9) +
      ylim(0, 1) +
      labs(x = "시간 (h)", y = "점유율 / 차단율", colour = NULL,
           title = "수용체 점유율 (경쟁적 가산)") + THEME
  })

  output$p_rate <- renderPlot({
    d <- sim()
    d %>% select(time, `혈관확장제 주입 (mg/h)` = RATE_VD,
                 `β차단제 주입 (mg/h)` = RATE_BB) %>%
      pivot_longer(-time) %>%
      ggplot(aes(time, value, colour = name)) + geom_line(linewidth = 0.9) +
      facet_wrap(~name, scales = "free_y") +
      labs(x = "시간 (h)", y = NULL, colour = NULL,
           title = "폐쇄회로 적정 — 간호사의 적정은 적분 제어기다") +
      THEME + theme(legend.position = "none")
  })

  output$t_exposure <- renderTable({
    d <- sim()
    data.frame(
      약물 = c("에스몰롤", "니카르디핀", "니트로프루시드", "라베탈롤",
               "클레비디핀", "시안화물", "티오시아네이트"),
      `정상상태 농도_mg_per_L` = sprintf(
        "%.4f", c(tail_mean(d, "A_ESM") / 255, tail_mean(d, "A_NICC") / 45,
                  tail_mean(d, "A_SNP") / 15, tail_mean(d, "A_LABC") / 420,
                  tail_mean(d, "A_CLV") / 12, tail_mean(d, "CYANIDE"),
                  tail_mean(d, "THIOCYANATE"))))
  }, striped = TRUE)

  # ---- 3 -------------------------------------------------------------------
  cmp <- reactive({
    keys <- c(input$cmpA, input$cmpB, input$cmpC)
    bind_rows(lapply(keys, function(k)
      simulate_arm(SCENARIOS[[k]]) %>% mutate(arm = k)))
  })

  output$p_dissoc <- renderPlot({
    cmp() %>%
      select(time, arm, `수축기압 (mmHg)` = SBP, `dP/dt max (mmHg/s)` = DPDT,
             `ASI (곱)` = ASI, `판막 응력-충격량` = ASI_FLAP) %>%
      pivot_longer(c(-time, -arm)) %>%
      ggplot(aes(time, value, colour = arm)) + geom_line(linewidth = 0.9) +
      facet_wrap(~name, scales = "free_y") +
      labs(x = "시간 (h)", y = NULL, colour = NULL,
           title = "압력은 내려가는데 충격량은 그대로다") + THEME
  })

  output$t_dissoc <- renderTable({
    cmp() %>% group_by(arm) %>%
      filter(time >= max(time) - 6) %>%
      summarise(`SBP` = mean(SBP), `HR` = mean(HR), `dP/dt` = mean(DPDT),
                `ASI` = mean(ASI), `ASI_flap` = mean(ASI_FLAP),
                `σ_FL` = mean(SIGMA_FL_SYS), .groups = "drop") %>%
      left_join(cmp() %>% group_by(arm) %>%
                  summarise(`파열구 성장 %` = 100 * (last(AEN) / first(AEN) - 1),
                            `사망률 %` = last(MORTALITY), .groups = "drop"),
                by = "arm")
  }, digits = 3, striped = TRUE)

  # ---- 4 -------------------------------------------------------------------
  output$p_divider <- renderPlot({
    d <- sim()
    d %>% select(time, `진강 수축기압` = SBP, `진강 확장기압` = DBP,
                 `가강 평균압` = PFL, `가강 수축기압` = PFL_SYSTOLIC,
                 `가강 확장기압` = PFL_DIASTOLIC) %>%
      pivot_longer(-time) %>%
      ggplot(aes(time, value, colour = name)) + geom_line(linewidth = 0.9) +
      labs(x = "시간 (h)", y = "압력 (mmHg)", colour = NULL,
           title = "두 내강의 압력 — 확장기압이 판별 변수다") + THEME
  })

  output$p_gre_sweep <- renderPlot({
    base <- user_pars()
    gres <- c(0.02, 0.05, 0.1, 0.25, 0.5, 1, 2, 3)
    res <- bind_rows(lapply(gres, function(g) {
      p <- base; p$GRE0 <- g
      d <- mod %>% param(p) %>% mrgsim(end = 24, delta = 1) %>% as_tibble()
      tibble(G_re = g,
             `가강 확장기압` = tail(d$PFL_DIASTOLIC, 1),
             `가강 평균압` = tail(d$PFL, 1),
             `가강 유량` = tail(d$FLOW_FALSE, 1) * 100,
             `진강 확장기압` = tail(d$DBP, 1))
    }))
    res %>% pivot_longer(-G_re) %>%
      ggplot(aes(G_re, value, colour = name)) +
      geom_line(linewidth = 0.9) + geom_point() + scale_x_log10() +
      labs(x = "재유입 전도도 G_re (log)", y = "mmHg  /  유량 ×100",
           colour = NULL,
           title = "재유입구를 닫으면 가강 확장기압이 전신 평균압으로 올라간다") +
      THEME
  })

  # ---- 5 -------------------------------------------------------------------
  output$p_sigma <- renderPlot({
    d <- sim()
    d %>% select(time, `가강 수축기 응력` = SIGMA_FL_SYS,
                 `가강 평균 응력` = SIGMA_FL, `진강 응력` = SIGMA_TL,
                 `판막 응력` = SIGMA_FLAP, `벽 강도 S` = WALL_STRENGTH) %>%
      pivot_longer(-time) %>%
      ggplot(aes(time, value, colour = name)) + geom_line(linewidth = 0.9) +
      labs(x = "시간 (h)", y = "응력 / 강도 (kPa)", colour = NULL,
           title = "벽 응력과 벽 강도가 만나는 곳이 파열이다") + THEME
  })

  output$p_hazard <- renderPlot({
    d <- sim()
    d %>% mutate(`연간 파열 위험 (%)` =
                   100 * (1 - exp(-8760 * 1.10e-4 *
                                    exp(pmin((STRESS_RATIO - 1) / 0.155, 30))))) %>%
      select(time, `σ/S` = STRESS_RATIO, `연간 파열 위험 (%)`) %>%
      pivot_longer(-time) %>%
      ggplot(aes(time, value, colour = name)) + geom_line(linewidth = 0.9) +
      facet_wrap(~name, scales = "free_y", ncol = 1) +
      labs(x = "시간 (h)", y = NULL, colour = NULL, title = "파열 위험") +
      THEME + theme(legend.position = "none")
  })

  output$p_wall <- renderPlot({
    d <- sim()
    d %>% select(time, `가강 외벽 두께 (mm)` = HFL, `엘라스틴 ELN` = ELN,
                 `외막 콜라겐 COL` = COL) %>%
      pivot_longer(-time) %>%
      ggplot(aes(time, value, colour = name)) + geom_line(linewidth = 0.9) +
      facet_wrap(~name, scales = "free_y", ncol = 1) +
      labs(x = "시간 (h)", y = NULL, colour = NULL,
           title = "벽 재형성 — 콜라겐은 유일한 음성 되먹임") +
      THEME + theme(legend.position = "none")
  })

  # ---- 6 -------------------------------------------------------------------
  output$p_paradox <- renderPlot({
    n <- 1.9; s <- 0.85; A <- 1.0
    thr <- seq(0, 1, by = 0.005)
    rel <- (1 + A * thr)^n * (1 - s * thr)
    star <- (n * A - s) / ((1 + n) * s * A)
    tibble(THR = thr, `상대 성장률` = rel / rel[1]) %>%
      ggplot(aes(THR, `상대 성장률`)) +
      geom_line(linewidth = 1.1, colour = "#b8860b") +
      geom_vline(xintercept = star, linetype = "dashed") +
      geom_hline(yintercept = 1, linetype = "dotted") +
      annotate("text", x = star + 0.02, y = 0.75, hjust = 0, size = 3.6,
               label = sprintf("해석해 THR* = %.3f", star)) +
      labs(x = "가강 혈전 분율 THR", y = "성장률 (혈전 0 대비)",
           title = "부분 혈전화 역설: 지지는 선형, 얇아짐은 지수 1.9") + THEME
  })

  output$p_thr_time <- renderPlot({
    d <- sim()
    d %>% select(time, `가강 본체 혈전` = THRFL, `원위 혈전` = THRD,
                 `입구 전도도` = G_ENTRY, `재유입 전도도` = G_REENTRY,
                 `가강 유량` = FLOW_FALSE) %>%
      pivot_longer(-time) %>%
      ggplot(aes(time, value, colour = name)) + geom_line(linewidth = 0.9) +
      labs(x = "시간 (h)", y = "분율 / 전도도 / 유량", colour = NULL,
           title = "혈전화의 시간 경과 — 봉인은 0.70을 넘어야 시작된다") + THEME
  })

  # ---- 7 -------------------------------------------------------------------
  output$p_perf <- renderPlot({
    d <- sim()
    d %>% select(time, `신장 관류` = QREN, `장간막 관류` = Q_MESENTERIC,
                 `척수 관류` = Q_SPINAL, `진강 압박` = TL_COLLAPSE) %>%
      pivot_longer(-time) %>%
      ggplot(aes(time, value, colour = name)) + geom_line(linewidth = 0.9) +
      labs(x = "시간 (h)", y = "상대 관류 / 압박 지수", colour = NULL,
           title = "분지 관류: 정적 폐색과 동적 폐색") + THEME
  })

  output$p_raas <- renderPlot({
    d <- sim()
    d %>% select(time, `레닌 활성도` = PRA, `안지오텐신 II` = ANGII,
                 `알도스테론` = ALDO, `SVR` = SVR, `SBP` = SBP) %>%
      pivot_longer(-time) %>%
      ggplot(aes(time, value, colour = name)) + geom_line(linewidth = 0.9) +
      facet_wrap(~name, scales = "free_y") +
      labs(x = "시간 (h)", y = NULL, colour = NULL,
           title = "레닌 함정 — 혈압을 내리면 신장이 되돌린다") +
      THEME + theme(legend.position = "none")
  })

  output$t_cn <- renderTable({
    d <- sim()
    data.frame(
      항목 = c("시안화물 최고 (mg/L)", "독성 역치 (mg/L)",
               "티오시아네이트 (mg/L)", "젖산 (mmol/L)",
               "크레아티닌 (mg/dL)", "eGFR"),
      값 = sprintf("%.2f", c(max(d$CYANIDE), 1.0, tail_mean(d, "THIOCYANATE"),
                            tail_mean(d, "LACT"), tail(d$CREAT, 1),
                            tail_mean(d, "EGFR"))))
  }, striped = TRUE)

  # ---- 8 -------------------------------------------------------------------
  output$p_cord <- renderPlot({
    d <- sim()
    d %>% select(time, `척수 관류 지수` = Q_SPINAL,
                 `허혈 부담 SCI` = SCI, `하지마비 확률` = PARAPLEGIA_PROB) %>%
      pivot_longer(-time) %>%
      ggplot(aes(time, value, colour = name)) + geom_line(linewidth = 0.9) +
      facet_wrap(~name, scales = "free_y", ncol = 1) +
      labs(x = "시간 (h)", y = NULL, colour = NULL,
           title = "척수 허혈은 회복 가능한 부담이며 창(window)이 있다") +
      THEME + theme(legend.position = "none")
  })

  output$p_collat <- renderPlot({
    d <- sim()
    d %>% select(time, `측부순환망 COLLAT` = COLLAT, `평균동맥압` = MAP) %>%
      pivot_longer(-time) %>%
      ggplot(aes(time, value, colour = name)) + geom_line(linewidth = 0.9) +
      facet_wrap(~name, scales = "free_y", ncol = 1) +
      labs(x = "시간 (h)", y = NULL, colour = NULL,
           title = "측부순환 모집 (기저 0.55, τ ≈ 5 일)") +
      THEME + theme(legend.position = "none")
  })

  # ---- 9 -------------------------------------------------------------------
  multi <- reactive({
    req(length(input$multi) > 0)
    bind_rows(lapply(input$multi, function(k)
      simulate_arm(SCENARIOS[[k]]) %>% mutate(arm = k)))
  })

  output$p_multi <- renderPlot({
    multi() %>%
      select(time, arm, SBP, HR, DPDT, ASI, SIGMA_FL_SYS, D_AORTA,
             STRESS_RATIO, MORTALITY) %>%
      pivot_longer(c(-time, -arm)) %>%
      ggplot(aes(time, value, colour = arm)) + geom_line(linewidth = 0.8) +
      facet_wrap(~name, scales = "free_y") +
      labs(x = "시간 (h)", y = NULL, colour = NULL,
           title = "시나리오 비교") +
      THEME + theme(legend.text = element_text(size = 8))
  })

  output$t_multi <- renderTable({
    multi() %>% group_by(arm) %>%
      summarise(`기간 (h)` = max(time),
                SBP = mean(SBP[time >= max(time) - 6]),
                HR = mean(HR[time >= max(time) - 6]),
                `dP/dt` = mean(DPDT[time >= max(time) - 6]),
                ASI = mean(ASI[time >= max(time) - 6]),
                `D_ao (mm)` = last(D_AORTA),
                `σ/S` = mean(STRESS_RATIO[time >= max(time) - 6]),
                `목표 내 %` = 100 * last(TIT) / max(time),
                `사망률 %` = last(MORTALITY), .groups = "drop")
  }, digits = 3, striped = TRUE)

  # ---- 10 ------------------------------------------------------------------
  output$p_end <- renderPlot({
    d <- sim()
    d %>% select(time, `통증 0-10` = PAIN, `목표 내 누적시간 (h)` = TIT,
                 `누적 위험` = CUMH, `생존 확률` = SURV) %>%
      pivot_longer(-time) %>%
      ggplot(aes(time, value, colour = name)) + geom_line(linewidth = 0.9) +
      facet_wrap(~name, scales = "free_y") +
      labs(x = "시간 (h)", y = NULL, colour = NULL,
           title = "임상 종말점") +
      THEME + theme(legend.position = "none")
  })

  output$p_bio <- renderPlot({
    d <- sim()
    d %>% select(time, `IL-6` = IL6, `MMP-9` = MMP9, `TGF-β` = TGFB,
                 `D-이합체` = DDIM, `젖산` = LACT) %>%
      pivot_longer(-time) %>%
      ggplot(aes(time, value, colour = name)) + geom_line(linewidth = 0.9) +
      labs(x = "시간 (h)", y = "상대값 / 농도", colour = NULL,
           title = "바이오마커") + THEME
  })

  output$p_growth <- renderPlot({
    d <- sim()
    yrs <- max(d$time) / 8760
    d %>% mutate(`성장 (mm)` = D_AORTA - first(D_AORTA)) %>%
      ggplot(aes(time / 8760, `성장 (mm)`)) +
      geom_line(linewidth = 1, colour = "#b8860b") +
      labs(x = "시간 (년)", y = "직경 증가 (mm)",
           title = sprintf("평균 성장률 %.2f mm/년",
                           (last(d$D_AORTA) - first(d$D_AORTA)) / max(yrs, 1e-9))) +
      THEME
  })

  output$t_end <- renderTable({
    d <- sim()
    hit55 <- d$time[which(d$D_AORTA >= 55)][1]
    data.frame(
      종말점 = c("최종 수축기압 (mmHg)", "최종 심박수 (bpm)",
                 "dP/dt max (mmHg/s)", "ASI", "판막 응력-충격량",
                 "파열구 성장 (%)", "총 직경 (mm)", "성장률 (mm/년)",
                 "55 mm 도달 (년)", "σ/S", "목표 내 시간 (%)",
                 "누적 사망률 (%)", "하지마비 확률"),
      값 = c(sprintf("%.1f", tail_mean(d, "SBP")),
             sprintf("%.1f", tail_mean(d, "HR")),
             sprintf("%.0f", tail_mean(d, "DPDT")),
             sprintf("%.3f", tail_mean(d, "ASI")),
             sprintf("%.3f", tail_mean(d, "ASI_FLAP")),
             sprintf("%.2f", 100 * (tail(d$AEN, 1) / d$AEN[1] - 1)),
             sprintf("%.1f", tail(d$D_AORTA, 1)),
             sprintf("%.2f", (tail(d$D_AORTA, 1) - d$D_AORTA[1]) /
                       max(max(d$time) / 8760, 1e-9)),
             ifelse(is.na(hit55), "도달 안 함", sprintf("%.2f", hit55 / 8760)),
             sprintf("%.3f", tail_mean(d, "STRESS_RATIO")),
             sprintf("%.1f", 100 * tail(d$TIT, 1) / max(d$time)),
             sprintf("%.2f", tail(d$MORTALITY, 1)),
             sprintf("%.3f", max(d$PARAPLEGIA_PROB))))
  }, striped = TRUE)
}

shinyApp(ui, server)
