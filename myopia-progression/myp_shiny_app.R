## =============================================================================
##  myp_shiny_app.R
##  Interactive dashboard for the progressive-myopia QSP model
##
##  10 tabs:
##    1  환자 프로파일        patient profile, baseline optics, risk stratum
##    2  성장 궤적            axial length and refraction over time
##    3  광학 분해            the three optical elements, and why SER misleads
##    4  맥락막 vs 공막       the two actuators that share one measurement
##    5  아트로핀 PK          ocular and systemic exposure
##    6  용량-반응            efficacy vs side effect, and the ratio
##    7  시나리오 비교        up to six arms side by side
##    8  리바운드             stopping, tapering, receptor up-regulation
##    9  공막 생물학          the matrix cascade behind the creep term
##   10  임상 엔드포인트      lifetime risk, CARE, NNT
##
##  Run:  shiny::runApp("myp_shiny_app.R")
##  Needs: shiny, mrgsolve, dplyr, tidyr, ggplot2, DT
## =============================================================================

library(shiny)
library(mrgsolve)
library(dplyr)
library(tidyr)
library(ggplot2)
library(DT)

`%||%` <- function(a, b) if (is.null(a)) b else a

MOD <- mread_cache("myp", "myp_mrgsolve_model.R")

## --- the baseline axial length that reproduces a prescribed refraction -----
## The optics block is inverted by bisection, exactly as solve_AL0() does in
## the Python reference model, so that t = 0 reproduces the requested SER.
ser_of <- function(AL, ACD = 3.60, LT = 3.45, PL = 22.60, CR = 7.80,
                   naq = 1.336, nvit = 1.336, vtx = 0.012) {
  PCORN <- (naq - 1) / (CR / 1000)
  d1 <- ACD + LT / 2
  d2 <- pmax(AL - d1, 5)
  L2 <- 1000 * nvit / d2 - PL
  L1 <- L2 / (1 + (d1 / 1000 / naq) * L2)
  FC <- L1 - PCORN
  FC / (1 + vtx * FC)
}
solve_AL0 <- function(SER0, ...) {
  lo <- 18; hi <- 34
  for (i in 1:200) {
    mid <- (lo + hi) / 2
    if (ser_of(mid, ...) > SER0) lo <- mid else hi <- mid
  }
  (lo + hi) / 2
}

DEVICE <- c("단일시 안경 (single-vision)"      = 0.15,
            "누진다초점 (PAL)"                 = -0.30,
            "MiSight 이중초점 콘택트렌즈"      = -0.90,
            "DIMS 안경"                        = -1.20,
            "HAL 렌즐릿 안경"                  = -2.20,
            "저교정 +0.75 D (해로움)"          = 0.90)

## simulate one arm; `days` in days
run_arm <- function(age0, ser0, npar, grs, ethn, outd, neard,
                    atro, atrostop, atrotap, trtpd, okon, rlrl,
                    fuiris, adhfix, trtstart, days) {
  al0 <- solve_AL0(ser0)
  pars <- list(AGE0 = age0, SER0 = ser0, AL0 = al0, NPAR = npar, GRS = grs,
               ETHN = ethn, OUTD = outd, NEARD = neard,
               ATROPCT = atro, ATROSTOP = atrostop, ATROTAP = atrotap,
               TRTPD = trtpd, OKON = okon, RLRLON = rlrl,
               FUIRISR = fuiris, ADHFIX = adhfix, TRTSTART = trtstart,
               OPTSTOP = if (atrostop < 1e8 && atro == 0) atrostop else 1e9)
  ev_drops <- if (atro > 0)
    ev(time = 0, amt = atro * 10 * 30 / 289.37 * 1000, cmt = "ATEAR",
       ii = 1, addl = max(1, floor(min(days, atrostop))) - 1) else ev()
  MOD %>% param(pars) %>%
    mrgsim(events = ev_drops, end = days, delta = 7, hmax = 1) %>%
    as_tibble() %>% mutate(yr = time / 365)
}

lifetime_vi <- function(AL) 1 / (1 + exp(-(-2.060 + 0.950 * (AL - 26))))

THEME <- theme_minimal(base_size = 12) +
  theme(panel.grid.minor = element_blank(),
        legend.position = "bottom", legend.title = element_blank(),
        strip.text = element_text(face = "bold"))

## =============================================================================
##  UI
## =============================================================================
ui <- fluidPage(
  titlePanel("소아 근시 진행 QSP 모델 — Progressive Myopia QSP Dashboard"),
  sidebarLayout(
    sidebarPanel(
      width = 3,
      h4("환자 (Patient)"),
      sliderInput("age0", "시작 연령 (yr)", 6, 15, 9.7, 0.1),
      sliderInput("ser0", "기저 굴절 SER (D)", -8, 1, -3, 0.25),
      selectInput("npar", "근시 부모 수", c(0, 1, 2), selected = 2),
      sliderInput("grs", "다유전자 위험점수", 0, 1, 0.70, 0.05),
      radioButtons("ethn", "인종", c("동아시아" = 1, "유럽" = 0), inline = TRUE),
      sliderInput("outd", "야외활동 (h/day)", 0.2, 4, 1.0, 0.1),
      sliderInput("neard", "근거리 작업 (h/day)", 0.5, 8, 3.0, 0.5),
      radioButtons("iris", "홍채 색조", c("어두움" = 1, "중간" = 2, "밝음" = 4),
                   inline = TRUE),
      hr(),
      h4("치료 (Treatment)"),
      selectInput("atro", "아트로핀 농도 (% w/v)",
                  c(0, 0.01, 0.025, 0.05, 0.1, 0.5, 1.0), selected = 0),
      selectInput("dev", "광학 처방", names(DEVICE), selected = names(DEVICE)[1]),
      checkboxInput("okon", "오르토케라톨로지 (ortho-K)", FALSE),
      checkboxInput("rlrl", "반복 저수준 적색광 (RLRL)", FALSE),
      sliderInput("trtstart", "치료 시작 (yr)", 0, 5, 0, 0.5),
      sliderInput("stopyr", "치료 중단 (yr; 10 = 중단 없음)", 0.5, 10, 10, 0.5),
      sliderInput("tapyr", "감량 기간 (yr)", 0, 2, 0, 0.25),
      sliderInput("adh", "순응도 (0 = 모델이 계산)", 0, 1, 0, 0.05),
      hr(),
      sliderInput("years", "추적 기간 (yr)", 1, 11, 8, 1),
      helpText("모든 곡선은 측정 안축장(각막→RPE)을 기준으로 하며, ",
               "굴절값은 상태변수가 아니라 광학식의 출력입니다.")
    ),
    mainPanel(
      width = 9,
      tabsetPanel(
        id = "tabs",
        # ------------------------------------------------------------------
        tabPanel("1. 환자 프로파일",
          br(),
          fluidRow(
            column(6, h4("기저 광학 (baseline optics)"), tableOutput("t_base")),
            column(6, h4("위험 계층 (risk stratum)"), tableOutput("t_risk"))),
          hr(),
          plotOutput("p_profile", height = "330px"),
          helpText("AL/CR 비가 3.0을 넘으면 굴절값이 아직 정상이어도 ",
                   "근시로 분류됩니다. 안축장은 굴절값보다 먼저 움직입니다.")),
        # ------------------------------------------------------------------
        tabPanel("2. 성장 궤적",
          br(), plotOutput("p_traj", height = "560px"),
          hr(), h4("연도별 진행 (annual progression)"),
          DTOutput("t_traj")),
        # ------------------------------------------------------------------
        tabPanel("3. 광학 분해",
          br(),
          h4("굴절값은 세 광학 요소의 차이일 뿐"),
          plotOutput("p_optics", height = "400px"),
          hr(),
          h4("dSER/dAL — 유도된 값이며 가정한 상수가 아니다"),
          plotOutput("p_dsda", height = "300px"),
          helpText("각막면에서는 안축장이 길어질수록 민감도가 31% 떨어지지만, ",
                   "안경면 환산이 이를 거의 상쇄하여 7%만 남습니다. ",
                   "그래서 '2.7 D/mm' 경험식은 안경 굴절에서만 잘 맞습니다.")),
        # ------------------------------------------------------------------
        tabPanel("4. 맥락막 vs 공막",
          br(),
          h4("하나의 측정값을 공유하는 두 개의 작동기"),
          plotOutput("p_chor", height = "400px"),
          hr(),
          h4("측정된 1년 효과 중 가역적(맥락막) 성분의 비율"),
          tableOutput("t_decomp"),
          helpText("생체계측기는 각막에서 RPE까지를 재기 때문에 ",
                   "맥락막이 두꺼워지면 안축장이 짧아진 것처럼 보입니다. ",
                   "이 성분은 치료를 멈추면 되돌아갑니다.")),
        # ------------------------------------------------------------------
        tabPanel("5. 아트로핀 PK",
          br(), plotOutput("p_pk", height = "420px"),
          hr(), h4("전신 노출 (systemic exposure)"),
          tableOutput("t_pk"),
          helpText("혈장 C_max는 0.01%에서 약 8 pg/mL, 1%에서 약 765 pg/mL로, ",
                   "1% 점안 후 측정된 300–900 pg/mL 범위와 일치합니다.")),
        # ------------------------------------------------------------------
        tabPanel("6. 용량-반응",
          br(),
          h4("효능과 부작용은 서로 다른 용량-반응 곡선을 가진다"),
          plotOutput("p_dr", height = "420px"),
          hr(), DTOutput("t_dr"),
          helpText("축성 반응의 Hill 기울기는 1.25이고 ED50은 약 0.05%로, ",
                   "LAMP가 권고한 농도와 사실상 같습니다. ",
                   "효능/산동 비는 0.01%에서 가장 나쁩니다.")),
        # ------------------------------------------------------------------
        tabPanel("7. 시나리오 비교",
          br(),
          checkboxGroupInput("arms", "비교할 치료 (최대 6개)",
            choices = c("무치료", "아트로핀 0.01%", "아트로핀 0.05%",
                        "아트로핀 1%", "DIMS 안경", "HAL 안경",
                        "MiSight", "오르토케라톨로지", "RLRL",
                        "Ortho-K + 아트로핀 0.01%"),
            selected = c("무치료", "아트로핀 0.05%", "DIMS 안경",
                         "오르토케라톨로지"), inline = TRUE),
          plotOutput("p_cmp", height = "480px"),
          hr(), DTOutput("t_cmp")),
        # ------------------------------------------------------------------
        tabPanel("8. 리바운드",
          br(),
          h4("중단 후: 맥락막의 빠른 붕괴 + 수용체 상향조절의 느린 과민성"),
          plotOutput("p_reb", height = "440px"),
          hr(), h4("중단 방식별 세척기 진행 (washout-year progression)"),
          tableOutput("t_reb"),
          helpText("리바운드는 두 성분의 합입니다. 수용체 상향조절은 ",
                   "τ = 90일로 감쇠하므로, 감량(taper)이 리바운드를 줄입니다.")),
        # ------------------------------------------------------------------
        tabPanel("9. 공막 생물학",
          br(), plotOutput("p_scl", height = "520px"),
          helpText("모든 결합은 거듭제곱 법칙이므로 어떤 변수도 음수가 될 수 ",
                   "없습니다. MMP-2:TIMP-2 균형이 기질 회전 방향을 정하는 ",
                   "단일 스칼라이며, 공막 크리프가 최종 공통 경로입니다.")),
        # ------------------------------------------------------------------
        tabPanel("10. 임상 엔드포인트",
          br(),
          fluidRow(
            column(6, h4("누적 절대 감소량 (CARE)"),
                   plotOutput("p_care", height = "320px")),
            column(6, h4("평생 시각장애 위험"),
                   plotOutput("p_vi", height = "320px"))),
          hr(),
          h4("1 mm의 가격은 일정하지 않다"),
          tableOutput("t_mm"),
          helpText("안축장 23 mm에서 1 mm는 평생 시각장애 위험 1.1%p를 ",
                   "뜻하지만 28 mm에서는 22.8%p입니다 — 20배 차이. ",
                   "따라서 '% 감소'는 잘못된 화폐 단위입니다."))
      )
    )
  )
)

## =============================================================================
##  SERVER
## =============================================================================
server <- function(input, output, session) {

  days <- reactive(input$years * 365)
  stopd <- reactive(if (input$stopyr >= 10) 1e9 else input$stopyr * 365)

  base_args <- reactive(list(
    age0 = input$age0, ser0 = input$ser0, npar = as.numeric(input$npar),
    grs = input$grs, ethn = as.numeric(input$ethn), outd = input$outd,
    neard = input$neard, fuiris = as.numeric(input$iris),
    adhfix = if (input$adh > 0) input$adh else -1,
    trtstart = input$trtstart * 365, days = days()))

  sim_one <- function(atro = 0, trtpd = 0.15, okon = 0, rlrl = 0,
                      atrostop = 1e9, atrotap = 0) {
    a <- base_args()
    do.call(run_arm, c(a, list(atro = atro, atrostop = atrostop,
                               atrotap = atrotap, trtpd = trtpd,
                               okon = okon, rlrl = rlrl)))
  }

  cur <- reactive(sim_one(atro = as.numeric(input$atro),
                          trtpd = DEVICE[[input$dev]],
                          okon = as.numeric(input$okon),
                          rlrl = as.numeric(input$rlrl),
                          atrostop = stopd(),
                          atrotap = input$tapyr * 365))
  ctl <- reactive(sim_one())

  ## ---- tab 1 ------------------------------------------------------------
  output$t_base <- renderTable({
    al0 <- solve_AL0(input$ser0)
    data.frame(
      항목 = c("기저 굴절 SER (D)", "기저 안축장 (mm)", "각막 곡률반경 (mm)",
               "각막 굴절력 (D)", "수정체 굴절력 (D)", "AL/CR 비",
               "dSER/dAL (안경면, D/mm)"),
      값 = c(sprintf("%.2f", input$ser0), sprintf("%.3f", al0), "7.80",
             sprintf("%.2f", 0.336 / 0.00780), "22.60",
             sprintf("%.3f", al0 / 7.80),
             sprintf("%.3f", (ser_of(al0 + 0.01) - ser_of(al0)) / 0.01)))
  }, striped = TRUE)

  output$t_risk <- renderTable({
    d <- cur(); e <- tail(d, 1)
    data.frame(
      항목 = c("최종 안축장 (mm)", "최종 SER (D)", "고도근시 도달 (≤ −6 D)",
               "안축장 26 mm 도달", "평생 시각장애 위험",
               "맥락막 두께 (µm)", "후방 공막 두께 (µm)"),
      값 = c(sprintf("%.3f", e$AL), sprintf("%.2f", e$SERTR),
             ifelse(any(d$HIGHMYO > 0), "예", "아니오"),
             ifelse(any(d$AL26 > 0), "예", "아니오"),
             sprintf("%.1f%%", 100 * lifetime_vi(e$AL)),
             sprintf("%.0f", e$CHT), sprintf("%.0f", e$SCT)))
  }, striped = TRUE)

  output$p_profile <- renderPlot({
    d <- cur()
    d %>% select(yr, AL, ALCR) %>%
      pivot_longer(-yr) %>%
      mutate(name = recode(name, AL = "측정 안축장 (mm)",
                           ALCR = "AL/CR 비")) %>%
      ggplot(aes(yr, value)) +
      geom_line(linewidth = 1.1, colour = "#2f6fb5") +
      geom_hline(data = data.frame(name = "AL/CR 비", y = 3.0),
                 aes(yintercept = y), linetype = 2, colour = "#b8564a") +
      facet_wrap(~name, scales = "free_y") +
      labs(x = "추적 기간 (yr)", y = NULL) + THEME
  })

  ## ---- tab 2 ------------------------------------------------------------
  output$p_traj <- renderPlot({
    bind_rows(mutate(ctl(), arm = "무치료 (untreated)"),
              mutate(cur(), arm = "현재 설정 (this setting)")) %>%
      select(yr, arm, AL, SERTR, CHT, CREEPS) %>%
      pivot_longer(c(AL, SERTR, CHT, CREEPS)) %>%
      mutate(name = recode(name,
        AL = "측정 안축장 (mm)", SERTR = "굴절값 SER (D)",
        CHT = "맥락막 두께 (µm)", CREEPS = "공막 크리프 (상대값)")) %>%
      ggplot(aes(yr, value, colour = arm)) +
      geom_line(linewidth = 1.1) +
      facet_wrap(~name, scales = "free_y") +
      scale_colour_manual(values = c("#8a8f98", "#2f6fb5")) +
      labs(x = "추적 기간 (yr)", y = NULL) + THEME
  })

  output$t_traj <- renderDT({
    d <- cur(); c0 <- ctl()
    yrs <- seq_len(input$years)
    pick <- function(df, y) df[which.min(abs(df$yr - y)), ]
    tibble(
      `연도` = yrs,
      `연령` = sapply(yrs, function(y) round(pick(d, y)$AGEy, 1)),
      `안축장 (mm)` = sapply(yrs, function(y) round(pick(d, y)$AL, 3)),
      `연간 안축장 (mm/yr)` = sapply(yrs, function(y)
        round(pick(d, y)$AL - pick(d, y - 1)$AL, 3)),
      `SER (D)` = sapply(yrs, function(y) round(pick(d, y)$SERTR, 2)),
      `무치료 안축장 (mm)` = sapply(yrs, function(y) round(pick(c0, y)$AL, 3)),
      `CARE (mm)` = sapply(yrs, function(y)
        round(pick(c0, y)$AL - pick(d, y)$AL, 3))
    ) %>% datatable(rownames = FALSE, options = list(dom = "t", pageLength = 12))
  })

  ## ---- tab 3 ------------------------------------------------------------
  output$p_optics <- renderPlot({
    d <- cur()
    al0 <- solve_AL0(input$ser0)
    ## contribution of each element to the refractive change, in dioptres
    d %>% transmute(
      yr,
      `안축장 기여 (element 1)` = ser_of(AL, 3.60, 3.45, 22.60, 7.80) -
        ser_of(al0, 3.60, 3.45, 22.60, 7.80),
      `수정체 기여 (element 3)` = ser_of(al0, 3.60, 3.45, PLENS, 7.80) -
        ser_of(al0, 3.60, 3.45, 22.60, 7.80),
      `실제 SER 변화 (net)` = SERTR - input$ser0) %>%
      pivot_longer(-yr) %>%
      ggplot(aes(yr, value, colour = name)) +
      geom_hline(yintercept = 0, colour = "grey70") +
      geom_line(linewidth = 1.1) +
      scale_colour_manual(values = c("#2f6fb5", "#b8564a", "#1f1f1f")) +
      labs(x = "추적 기간 (yr)", y = "굴절 변화 (D)") + THEME
  })

  output$p_dsda <- renderPlot({
    al <- seq(22.5, 31, 0.1)
    spec <- (ser_of(al + 0.01) - ser_of(al)) / 0.01
    ## the corneal-plane sensitivity, obtained by undoing the vertex conversion
    fc <- function(s) s / (1 - 0.012 * s)
    corn <- (fc(ser_of(al + 0.01)) - fc(ser_of(al))) / 0.01
    tibble(al, `안경면 (spectacle plane)` = spec,
           `각막면 (corneal plane)` = corn) %>%
      pivot_longer(-al) %>%
      ggplot(aes(al, value, colour = name)) +
      geom_line(linewidth = 1.1) +
      scale_colour_manual(values = c("#b8564a", "#2f6fb5")) +
      labs(x = "안축장 (mm)", y = "dSER/dAL (D/mm)") + THEME
  })

  ## ---- tab 4 ------------------------------------------------------------
  output$p_chor <- renderPlot({
    d <- cur()
    al0 <- solve_AL0(input$ser0)
    d %>% transmute(yr,
      `공막 안축장 (permanent)` = al0 + AXCUM,
      `측정 안축장 (what a biometer reads)` = AL,
      `맥락막 두께 (µm, 우측 축 아님)` = CHT / 1000 + al0) %>%
      pivot_longer(-yr) %>%
      ggplot(aes(yr, value, colour = name)) +
      geom_line(linewidth = 1.1) +
      scale_colour_manual(values = c("#c9955f", "#1f1f1f", "#2f6fb5")) +
      labs(x = "추적 기간 (yr)", y = "mm") + THEME
  })

  output$t_decomp <- renderTable({
    a <- base_args(); a$days <- 365
    c0 <- do.call(run_arm, c(a, list(atro = 0, atrostop = 1e9, atrotap = 0,
                                     trtpd = 0.15, okon = 0, rlrl = 0)))
    al0 <- solve_AL0(input$ser0)
    arms <- list(`아트로핀 0.01%` = list(atro = 0.01),
                 `아트로핀 0.05%` = list(atro = 0.05),
                 `아트로핀 1%`    = list(atro = 1.00),
                 `DIMS 안경`      = list(trtpd = -1.20),
                 `오르토케라톨로지` = list(okon = 1),
                 `RLRL`           = list(rlrl = 1))
    rows <- lapply(names(arms), function(nm) {
      z <- arms[[nm]]
      d <- do.call(run_arm, c(a, list(
        atro = z$atro %||% 0, atrostop = 1e9, atrotap = 0,
        trtpd = z$trtpd %||% 0.15, okon = z$okon %||% 0, rlrl = z$rlrl %||% 0)))
      dm <- (tail(c0, 1)$AL - al0) - (tail(d, 1)$AL - al0)
      ds <- tail(c0, 1)$AXCUM - tail(d, 1)$AXCUM
      data.frame(치료 = nm,
                 `측정 감소 (mm)` = round(dm, 4),
                 `공막 감소 (mm)` = round(ds, 4),
                 `맥락막 성분 (mm)` = round(dm - ds, 4),
                 `가역 비율` = sprintf("%.1f%%", 100 * (dm - ds) / dm),
                 check.names = FALSE)
    })
    bind_rows(rows)
  }, striped = TRUE)

  ## ---- tab 5 ------------------------------------------------------------
  output$p_pk <- renderPlot({
    if (as.numeric(input$atro) == 0)
      return(ggplot() + annotate("text", 0, 0, label = "아트로핀 농도를 선택하세요") +
               theme_void())
    a <- base_args(); a$days <- 14
    d <- do.call(run_arm, c(a, list(atro = as.numeric(input$atro),
                                    atrostop = 1e9, atrotap = 0, trtpd = 0.15,
                                    okon = 0, rlrl = 0)))
    d %>% select(time, CPLA_PG, CCHOR_NM, CSCL_NM) %>%
      pivot_longer(-time) %>%
      mutate(name = recode(name, CPLA_PG = "혈장 (pg/mL)",
                           CCHOR_NM = "맥락막 (nM)", CSCL_NM = "공막 (nM)")) %>%
      ggplot(aes(time, value)) + geom_line(linewidth = 1) +
      facet_wrap(~name, scales = "free_y", ncol = 1) +
      labs(x = "시간 (day)", y = NULL) + THEME
  })

  output$t_pk <- renderTable({
    a <- base_args(); a$days <- 30
    rows <- lapply(c(0.01, 0.025, 0.05, 0.1, 0.5, 1.0), function(pc) {
      d <- do.call(run_arm, c(a, list(atro = pc, atrostop = 1e9, atrotap = 0,
                                      trtpd = 0.15, okon = 0, rlrl = 0)))
      data.frame(`농도 (%)` = pc,
                 `혈장 C_max (pg/mL)` = round(max(d$CPLA_PG), 1),
                 `동공 (mm)` = round(tail(d, 1)$PUPD, 2),
                 `조절력 (D)` = round(tail(d, 1)$AAMP, 2),
                 `구동력 억제 EATR` = round(tail(d, 1)$EATRo, 3),
                 check.names = FALSE)
    })
    bind_rows(rows)
  }, striped = TRUE)

  ## ---- tab 6 ------------------------------------------------------------
  dr_tab <- reactive({
    a <- base_args(); a$days <- 365
    al0 <- solve_AL0(input$ser0)
    c0 <- do.call(run_arm, c(a, list(atro = 0, atrostop = 1e9, atrotap = 0,
                                     trtpd = 0.15, okon = 0, rlrl = 0)))
    base <- tail(c0, 1)$AL - al0
    doses <- c(0.005, 0.01, 0.025, 0.05, 0.1, 0.25, 0.5, 1.0)
    bind_rows(lapply(doses, function(pc) {
      d <- tail(do.call(run_arm, c(a, list(atro = pc, atrostop = 1e9,
        atrotap = 0, trtpd = 0.15, okon = 0, rlrl = 0))), 1)
      eff <- 1 - (d$AL - al0) / base
      myd <- (d$PUPD - 4.60) / 3.80
      tibble(dose = pc, efficacy = eff, mydriasis = myd,
             accom_loss = 1 - d$AAMP / 13.40,
             ratio = eff / pmax(myd, 1e-6), adherence = d$ADH)
    }))
  })

  output$p_dr <- renderPlot({
    dr_tab() %>%
      select(dose, `축성 효능 (efficacy)` = efficacy,
             `산동 (mydriasis)` = mydriasis,
             `조절력 손실 (accommodation loss)` = accom_loss,
             `효능/산동 비 (ratio)` = ratio) %>%
      pivot_longer(-dose) %>%
      ggplot(aes(dose, value, colour = name)) +
      geom_line(linewidth = 1.1) + geom_point(size = 1.8) +
      scale_x_log10(breaks = c(0.005, 0.01, 0.025, 0.05, 0.1, 0.25, 0.5, 1)) +
      geom_vline(xintercept = 0.05, linetype = 2, colour = "grey60") +
      labs(x = "아트로핀 농도 (% w/v, log)", y = NULL) + THEME
  })

  output$t_dr <- renderDT({
    dr_tab() %>%
      transmute(`농도 (%)` = dose,
                `축성 효능` = round(efficacy, 3),
                `산동 분율` = round(mydriasis, 3),
                `조절력 손실 분율` = round(accom_loss, 3),
                `효능/산동 비` = round(ratio, 2),
                `순응도` = round(adherence, 2)) %>%
      datatable(rownames = FALSE, options = list(dom = "t"))
  })

  ## ---- tab 7 ------------------------------------------------------------
  ARMS <- list(
    `무치료`                   = list(),
    `아트로핀 0.01%`           = list(atro = 0.01),
    `아트로핀 0.05%`           = list(atro = 0.05),
    `아트로핀 1%`              = list(atro = 1.00),
    `DIMS 안경`                = list(trtpd = -1.20),
    `HAL 안경`                 = list(trtpd = -2.20),
    `MiSight`                  = list(trtpd = -0.90),
    `오르토케라톨로지`         = list(okon = 1, trtpd = 0),
    `RLRL`                     = list(rlrl = 1),
    `Ortho-K + 아트로핀 0.01%` = list(okon = 1, trtpd = 0, atro = 0.01))

  cmp <- reactive({
    sel <- head(input$arms, 6)
    bind_rows(lapply(sel, function(nm) {
      z <- ARMS[[nm]]
      mutate(sim_one(atro = z$atro %||% 0, trtpd = z$trtpd %||% 0.15,
                     okon = z$okon %||% 0, rlrl = z$rlrl %||% 0), arm = nm)
    }))
  })

  output$p_cmp <- renderPlot({
    cmp() %>% select(yr, arm, AL, SERTR, CHT, VIRISK) %>%
      pivot_longer(c(AL, SERTR, CHT, VIRISK)) %>%
      mutate(name = recode(name, AL = "측정 안축장 (mm)",
                           SERTR = "굴절값 SER (D)",
                           CHT = "맥락막 두께 (µm)",
                           VIRISK = "평생 시각장애 위험")) %>%
      ggplot(aes(yr, value, colour = arm)) + geom_line(linewidth = 1) +
      facet_wrap(~name, scales = "free_y") +
      labs(x = "추적 기간 (yr)", y = NULL) + THEME
  })

  output$t_cmp <- renderDT({
    al0 <- solve_AL0(input$ser0)
    cmp() %>% group_by(arm) %>% slice_tail(n = 1) %>% ungroup() %>%
      transmute(`치료` = arm,
                `안축장 (mm)` = round(AL, 3),
                `총 신장 (mm)` = round(AL - al0, 3),
                `SER (D)` = round(SERTR, 2),
                `평생 시각장애 위험` = sprintf("%.1f%%", 100 * VIRISK),
                `동공 (mm)` = round(PUPD, 2),
                `조절력 (D)` = round(AAMP, 2)) %>%
      arrange(`안축장 (mm)`) %>%
      datatable(rownames = FALSE, options = list(dom = "t"))
  })

  ## ---- tab 8 ------------------------------------------------------------
  output$p_reb <- renderPlot({
    a <- base_args(); a$days <- 3 * 365
    mk <- function(nm, pc, tap) mutate(do.call(run_arm, c(a, list(
      atro = pc, atrostop = 2 * 365, atrotap = tap, trtpd = 0.15,
      okon = 0, rlrl = 0))), arm = nm)
    bind_rows(
      mutate(do.call(run_arm, c(a, list(atro = 0, atrostop = 1e9, atrotap = 0,
        trtpd = 0.15, okon = 0, rlrl = 0))), arm = "무치료"),
      mk("0.01% → 중단", 0.01, 0),
      mk("0.5% → 급중단", 0.50, 0),
      mk("0.5% → 1년 감량", 0.50, 365)) %>%
      select(yr, arm, AL, SERTR, CHT) %>%
      pivot_longer(c(AL, SERTR, CHT)) %>%
      mutate(name = recode(name, AL = "측정 안축장 (mm)",
                           SERTR = "굴절값 SER (D)",
                           CHT = "맥락막 두께 (µm)")) %>%
      ggplot(aes(yr, value, colour = arm)) +
      geom_vline(xintercept = 2, linetype = 2, colour = "grey60") +
      geom_line(linewidth = 1) +
      facet_wrap(~name, scales = "free_y", ncol = 1) +
      labs(x = "추적 기간 (yr)", y = NULL) + THEME
  })

  output$t_reb <- renderTable({
    a <- base_args(); a$days <- 3 * 365
    pick <- function(df, y) df[which.min(abs(df$yr - y)), ]
    rows <- lapply(list(c(0.01, 0), c(0.10, 0), c(0.50, 0),
                        c(0.50, 90), c(0.50, 180), c(0.50, 365)),
      function(z) {
        d <- do.call(run_arm, c(a, list(atro = z[1], atrostop = 2 * 365,
          atrotap = z[2], trtpd = 0.15, okon = 0, rlrl = 0)))
        data.frame(`농도 (%)` = z[1], `감량 (day)` = z[2],
                   `세척기 SER 변화 (D)` =
                     round(pick(d, 3)$SERTR - pick(d, 2)$SERTR, 3),
                   `세척기 안축장 (mm)` =
                     round(pick(d, 3)$AL - pick(d, 2)$AL, 3),
                   check.names = FALSE)
      })
    bind_rows(rows)
  }, striped = TRUE)

  ## ---- tab 9 ------------------------------------------------------------
  output$p_scl <- renderPlot({
    d <- cur()
    d %>% select(yr, MTRo, CREEPS, SCT, CHOX, PSTAPH, LACQ) %>%
      pivot_longer(-yr) %>%
      mutate(name = recode(name,
        MTRo = "MMP-2 : TIMP-2 균형", CREEPS = "공막 크리프 (상대값)",
        SCT = "후방 공막 두께 (µm)", CHOX = "맥락막 저산소 지수",
        PSTAPH = "후방 포도상 지수", LACQ = "래커 균열 / 위축 지수")) %>%
      ggplot(aes(yr, value)) + geom_line(linewidth = 1, colour = "#a05a95") +
      facet_wrap(~name, scales = "free_y") +
      labs(x = "추적 기간 (yr)", y = NULL) + THEME
  })

  ## ---- tab 10 -----------------------------------------------------------
  output$p_care <- renderPlot({
    c0 <- ctl(); d <- cur()
    pick <- function(df, y) df[which.min(abs(df$yr - y)), ]
    yrs <- seq_len(input$years)
    tibble(yr = yrs,
           CARE = sapply(yrs, function(y) pick(c0, y)$AL - pick(d, y)$AL)) %>%
      mutate(increment = CARE - lag(CARE, default = 0)) %>%
      pivot_longer(-yr) %>%
      mutate(name = recode(name, CARE = "누적 (CARE, mm)",
                           increment = "연간 증가분 (mm)")) %>%
      ggplot(aes(yr, value, fill = name)) +
      geom_col(position = "dodge") +
      scale_fill_manual(values = c("#2f6fb5", "#c9955f")) +
      labs(x = "연도", y = "mm") + THEME
  })

  output$p_vi <- renderPlot({
    bind_rows(mutate(ctl(), arm = "무치료"), mutate(cur(), arm = "현재 설정")) %>%
      ggplot(aes(yr, 100 * VIRISK, colour = arm)) +
      geom_line(linewidth = 1.1) +
      scale_colour_manual(values = c("#8a8f98", "#2f6fb5")) +
      labs(x = "추적 기간 (yr)", y = "평생 시각장애 위험 (%)") + THEME
  })

  output$t_mm <- renderTable({
    al <- 23:30
    data.frame(`안축장 (mm)` = al,
               `평생 위험` = sprintf("%.1f%%", 100 * lifetime_vi(al)),
               `+1 mm 후` = sprintf("%.1f%%", 100 * lifetime_vi(al + 1)),
               `증가 (%p)` = round(100 * (lifetime_vi(al + 1) -
                                            lifetime_vi(al)), 2),
               check.names = FALSE)
  }, striped = TRUE)
}

shinyApp(ui, server)
