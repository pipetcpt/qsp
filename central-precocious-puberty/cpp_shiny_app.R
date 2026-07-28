# =============================================================================
#  cpp_shiny_app.R
#  Central Precocious Puberty (CPP) — interactive QSP dashboard
#  ---------------------------------------------------------------------------
#  성조숙증(중추성) QSP 대시보드 · 10 tabs
#
#  Run with:
#      shiny::runApp("cpp_shiny_app.R")
#  (the file sources cpp_mrgsolve_model.R from the same directory, which
#   compiles the 44-ODE mrgsolve model)
#
#  DESIGN PRINCIPLE OF THIS DASHBOARD
#  ----------------------------------
#  The app is built around the one thing the model exists to show: the SIGN of
#  the treatment effect on adult height is not a property of the drug, it is a
#  property of how much growth-plate reserve remains.  Every tab therefore
#  displays the treated trajectory AGAINST ITS OWN UNTREATED COUNTERFACTUAL
#  (same phenotype, no drug), because a treated trajectory on its own is
#  uninterpretable — and against the NORMAL reference girl/boy, because that is
#  the height the family is actually hoping for.
#
#  A second deliberate choice: the "Monitoring" tab plots growth velocity and
#  dBA/dCA side by side and labels which one is trustworthy.  Growth velocity
#  FALLS on effective therapy; a dashboard that shows it without that context
#  invites exactly the wrong clinical conclusion.
# =============================================================================

library(shiny)
library(mrgsolve)
library(dplyr)
library(ggplot2)
library(tidyr)

source("cpp_mrgsolve_model.R", local = TRUE, chdir = TRUE)

THEME <- theme_bw(base_size = 12) +
  theme(panel.grid.minor = element_blank(),
        strip.background = element_rect(fill = "grey92"),
        legend.position = "bottom")

PAL <- c(untreated = "#c0392b", treated = "#2471a3", normal = "#7f8c8d",
         alt = "#e67e22", alt2 = "#16a085")

# =============================================================================
#  UI
# =============================================================================
ui <- fluidPage(
  titlePanel("중추성 성조숙증 QSP 모델 · Central Precocious Puberty QSP Dashboard"),
  p(strong("모델의 핵심 명제:"),
    "에스트라디올은 성인키 적분식에 ", strong("서로 반대 부호로 두 번"),
    " 들어간다 — (+) GH/IGF-1 증폭을 통해 성장속도를 올리고,",
    "(−) 성장판 예비능을 소모해 성장 기간을 끝낸다.",
    "GnRH 작용제는 두 팔을 동시에 없애므로, 치료가 이득이 되는지는",
    strong("약이 아니라 남은 성장판 예비능"), "이 결정한다."),
  hr(),
  sidebarLayout(
    sidebarPanel(
      width = 3,
      h4("① 환자 (Phenotype)"),
      radioButtons("sex", "성별 / Sex", c("여아 Girl" = "F", "남아 Boy" = "M"),
                   inline = TRUE),
      sliderInput("mk0", "MKRN3 브레이크 잔존 (MK0) — 사춘기 시작은 결과물",
                  min = 0.20, max = 1.00, value = 0.40, step = 0.02),
      helpText("MK0 0.40 → 유방발육 6.5세 · 1.00 → 10.3세 · 0.66 → 8.6세(서행형)"),
      sliderInput("ht0", "만 5세 신장 (cm)", 96, 120, 108, step = 0.5),
      sliderInput("ba0", "만 5세 골연령 (yr)", 4.0, 7.0, 5.0, step = 0.1),
      checkboxInput("mas", "GnRH-비의존성 (McCune-Albright, GNAS)", FALSE),
      checkboxInput("lesion", "시상하부 병변/하마르토마 (ectopic drive)", FALSE),

      h4("② 치료 (Therapy)"),
      selectInput("drug", "제형 / Formulation", c(
        "없음 (untreated)"                  = "none",
        "류프롤리드 3.75 mg IM q28d"        = "leup1m",
        "류프롤리드 7.5 mg IM q28d"         = "leup1mH",
        "류프롤리드 11.25 mg IM q12wk"      = "leup3m",
        "류프롤리드 30 mg IM q12wk"         = "leup3mH",
        "트립토렐린 11.25 mg IM q12wk"      = "trip3m",
        "트립토렐린 22.5 mg IM q24wk"       = "trip6m",
        "히스트렐린 50 mg 임플란트 q12mo"   = "hist",
        "나파렐린 비강분무 1800 µg/day"     = "naf",
        "GnRH 길항제 18 mg q28d"            = "antag"),
        selected = "leup3m"),
      sliderInput("start", "치료 시작 연령 (yr)", 5.5, 12.0, 7.4, step = 0.1),
      sliderInput("dur", "치료 기간 (yr)", 0.5, 8.0, 5.0, step = 0.5),
      sliderInput("late", "실제 투여 간격 지연 배수 (adherence proxy)",
                  1.0, 2.0, 1.0, step = 0.05),
      helpText("1.0 = 정확히 지킴, 1.5 = 28일 제형을 42일마다 맞음"),

      h4("③ 병용 (Add-on)"),
      checkboxInput("ai", "아로마타제 억제제 (아나스트로졸 1 mg/day)", FALSE),
      checkboxInput("ai_potent", "  └ 고효능 AI (레트로졸형)", FALSE),
      checkboxInput("gh", "rhGH 0.043 mg/kg/day", FALSE),
      checkboxInput("tam", "타목시펜 20 mg/day", FALSE),
      checkboxInput("cavd", "칼슘 + 비타민 D", FALSE),

      h4("④ 시뮬레이션"),
      sliderInput("years", "시뮬레이션 기간 (yr from age 5)", 10, 18, 16.5,
                  step = 0.5),
      actionButton("go", "시뮬레이션 실행", class = "btn-primary btn-block")
    ),

    mainPanel(
      width = 9,
      tabsetPanel(
        id = "tabs",
        # ---------------------------------------------------------------
        tabPanel(
          "1 · 환자 프로파일",
          h4("사춘기 시작 시점은 입력이 아니라 결과물이다"),
          p("MKRN3 브레이크(MK0) 하나만 바꾸면 유방발육·초경·골연령·최종키가",
            "모두 함께 움직인다. 아래 표는 지금 설정된 환자와, 같은 환자의",
            "무치료 대조군, 그리고 정상 기준 아동을 나란히 보여준다."),
          tableOutput("tbl_profile"),
          hr(),
          h4("MK0 스윕 — 표현형 스펙트럼"),
          plotOutput("plt_mk_sweep", height = "330px"),
          helpText("서행성 변이형(MK0 0.60-0.76)은 조기 유방발육으로 내원하지만",
                   "치료 이득이 매우 작다 — 과잉치료 문제의 정량적 형태.")
        ),
        # ---------------------------------------------------------------
        tabPanel(
          "2 · PK · 수용체",
          h4("작용제는 수용체를 막지 않는다 — 박동성 코드를 파괴한다"),
          plotOutput("plt_pk", height = "620px"),
          helpText("fa = 작용제 점유율, fx = 길항제 점유율, RS = 감작 수용체 분율.",
                   "작용제는 AINT = 1.6의 내인활성 때문에 첫 투여 시 자극이",
                   "기저보다 올라간다(flare) — 억제는 RS가 붕괴한 뒤에 온다.",
                   "길항제는 flare가 없지만 탈감작도 없으므로 점유율을 계속",
                   "유지해야 한다.")
        ),
        # ---------------------------------------------------------------
        tabPanel(
          "3 · HPG 축",
          h4("KNDy 박동발생기 → 박동 빈도 → LH/FSH → 성선 스테로이드"),
          plotOutput("plt_hpg", height = "620px"),
          helpText("PULS(박동/24h)는 명시적 상태변수다. 지속적 작용제 노출은",
                   "빈도 정보를 0으로 만든다 — 이것이 치료 원리다.")
        ),
        # ---------------------------------------------------------------
        tabPanel(
          "4 · 성장 · 골연령",
          h4("두 개의 시계: 신장 시계와 골 시계"),
          plotOutput("plt_growth", height = "640px"),
          helpText("골연령은 '같은 나이 정상 아동 대비 성숙 속도'로 정의된다",
                   "(Greulich-Pyle 아틀라스의 정의 그대로). 그래서 정상 아동은",
                   "정의상 ΔBA/ΔCA = 1.0이고, 억제된 아동은 CA 10-12세에서",
                   "0.6-0.7로 '정지'처럼 보인다 — 모델에 정지 항은 없다.")
        ),
        # ---------------------------------------------------------------
        tabPanel(
          "5 · 임상 엔드포인트",
          h4("최종 성인키와, 무치료 대조군 대비 이득"),
          tableOutput("tbl_endpoints"),
          hr(),
          plotOutput("plt_height", height = "380px"),
          helpText("BP-PAH = Bayley-Pinneau 예측 성인키(임상에서 계산하는 값),",
                   "true final = 모델이 실제로 적분한 최종키. 두 값의 차이가",
                   "곧 임상 예측식의 편향이다.")
        ),
        # ---------------------------------------------------------------
        tabPanel(
          "6 · 사인 플립 (핵심 결과)",
          h4("치료 이득의 부호는 시작 시점의 골연령이 정한다"),
          p("아래 곡선의 어디에도 '8세 이전에 치료' 규칙은 코딩되어 있지 않다.",
            "규칙은 (+)팔과 (−)팔의 경쟁에서 나오는 ", strong("출력"), "이다."),
          plotOutput("plt_signflip", height = "420px"),
          tableOutput("tbl_signflip")
        ),
        # ---------------------------------------------------------------
        tabPanel(
          "7 · 모니터링 (함정)",
          h4("성장속도는 신뢰할 수 없는 지표다"),
          p("효과적인 억제는 성장속도를 정상 사춘기전 값보다도 아래로",
            "떨어뜨린다(성 스테로이드와 그 GH/IGF-1 증폭을 동시에 잃기 때문).",
            strong("성장속도 저하는 치료 실패가 아니라 치료 성공의 신호다.")),
          plotOutput("plt_monitor", height = "560px"),
          hr(),
          tableOutput("tbl_monitor")
        ),
        # ---------------------------------------------------------------
        tabPanel(
          "8 · 제형 비교",
          h4("결정 변수는 효능이 아니라 최저농도 구간의 커버리지다"),
          plotOutput("plt_forms", height = "420px"),
          tableOutput("tbl_forms"),
          helpText("정확히 투여된 데포 제형들은 서로 거의 구분되지 않는다.",
                   "차이는 (a) 비강분무, (b) 지연된 주사, (c) 저용량에서 생긴다.")
        ),
        # ---------------------------------------------------------------
        tabPanel(
          "9 · 바이오마커",
          h4("표적 조직 반응 — 자궁 부피 · Tanner 유방 · 자궁내막"),
          plotOutput("plt_biomarker", height = "620px"),
          helpText("음모(부신 안드로겐)는 GnRH 작용제로 억제되지 않는다 —",
                   "치료 중에도 진행하므로 보호자 상담에서 반드시 언급해야 한다.")
        ),
        # ---------------------------------------------------------------
        tabPanel(
          "10 · 안전성 · 삶의 질",
          h4("치료의 비용: 골밀도 · 체질량 · 열성 홍조 · 심리사회 지표"),
          plotOutput("plt_safety", height = "620px"),
          helpText("골밀도 Z-점수의 '하강'은 대부분 진행된 기저치가 정상으로",
                   "회귀하는 것이며, 최대 골질량은 회복된다.")
        )
      )
    )
  )
)

# =============================================================================
#  SERVER
# =============================================================================
server <- function(input, output, session) {

  base_pars <- reactive({
    male <- input$sex == "M"
    p <- if (male) pheno_boy(input$mk0) else pheno_girl(input$mk0)
    p$HT00 <- input$ht0
    p$BA00 <- input$ba0
    p$CAVD <- as.numeric(input$cavd)
    if (input$mas) { p$AUTSET <- 0.42; p$KSEC <- 0.22 }
    if (input$lesion) p$KNDLES <- 0.25
    if (input$ai_potent) { p$IMAXAI <- 0.992; p$IC50AI <- 2.2 }
    p
  })

  build_events <- reactive({
    a0 <- input$start; a1 <- input$start + input$dur
    L <- input$late
    ev <- switch(input$drug,
      none    = NULL,
      leup1m  = ev_leup(3.75, 28 * L, a0, a1),
      leup1mH = ev_leup(7.5, 28 * L, a0, a1),
      leup3m  = ev_leup(11.25, 84 * L, a0, a1),
      leup3mH = ev_leup(30, 84 * L, a0, a1),
      trip3m  = ev_trip(11.25, 84 * L, a0, a1),
      trip6m  = ev_trip(22.5, 168 * L, a0, a1),
      hist    = ev_hist(a0, a1),
      naf     = ev_naf(1800, a0, a1),
      antag   = ev_antag(18, 28 * L, a0, a1))
    if (input$ai)  ev <- c(ev, ev_ai(if (input$ai_potent) 2.5 else 1.0, a0, a1))
    if (input$gh)  ev <- c(ev, ev_gh(0.043 * 30, a0, a1))
    if (input$tam) ev <- c(ev, ev_tam(20, a0, a1))
    ev
  })

  sims <- eventReactive(input$go, {
    p <- base_pars()
    yrs <- input$years
    dT <- sim_cpp(p, build_events(), years = yrs, delta = 0.5)
    dU <- sim_cpp(p, NULL, years = yrs, delta = 0.5)
    pN <- p; pN$MK0 <- MK_NORMAL
    pN$AUTSET <- 0; pN$KSEC <- 0; pN$KNDLES <- 0
    dN <- sim_cpp(pN, NULL, years = yrs, delta = 0.5)
    list(treated = dT, untreated = dU, normal = dN, male = input$sex == "M")
  }, ignoreNULL = FALSE)

  long3 <- function(cols) {
    s <- sims()
    bind_rows(
      mutate(s$treated[, c("CA_out", cols)], arm = "treated"),
      mutate(s$untreated[, c("CA_out", cols)], arm = "untreated"),
      mutate(s$normal[, c("CA_out", cols)], arm = "normal")) %>%
      pivot_longer(all_of(cols), names_to = "var", values_to = "value")
  }

  panel3 <- function(cols, labs, logy = FALSE) {
    d <- long3(cols)
    d$var <- factor(d$var, levels = cols, labels = labs)
    g <- ggplot(d, aes(CA_out, value, colour = arm)) +
      geom_line(linewidth = 0.7) +
      facet_wrap(~var, scales = "free_y") +
      scale_colour_manual(values = PAL, name = NULL) +
      labs(x = "역연령 Chronological age (yr)", y = NULL) + THEME
    if (logy) g <- g + scale_y_log10()
    g
  }

  # ---- Tab 1 -----------------------------------------------------------
  output$tbl_profile <- renderTable({
    s <- sims()
    f <- function(d) c(
      `유방/사춘기 시작 (yr)` = thelarche(d),
      `초경 (yr)` = menarche(d),
      `골성숙 완료 (yr)` = fusion_age(d),
      `최종 성인키 (cm)` = adult_height(d),
      `최대 성장속도 (cm/yr)` = max(d$GVo),
      `최대 E2 (pg/mL)` = max(d$E2),
      `CA 8세 골연령 (yr)` = at_age(d, 8, "BA"),
      `ΔBA/ΔCA 7-9세` = ba_ca_ratio(d, 7, 9),
      `누적 E2 (pg/mL·yr)` = tail(d$CUME2, 1),
      `성인 BMD Z` = tail(d$BMDZ, 1))
    m <- cbind(`현재 설정 (치료)` = f(s$treated),
               `같은 환자 · 무치료` = f(s$untreated),
               `정상 기준 아동` = f(s$normal))
    as.data.frame(round(m, 2))
  }, rownames = TRUE, digits = 2)

  output$plt_mk_sweep <- renderPlot({
    male <- input$sex == "M"
    mks <- seq(0.28, 0.92, by = 0.08)
    res <- do.call(rbind, lapply(mks, function(mk) {
      p <- if (male) pheno_boy(mk) else pheno_girl(mk)
      p$HT00 <- input$ht0; p$BA00 <- input$ba0
      dU <- sim_cpp(p, NULL, years = 17, delta = 2)
      a <- min(max(ifelse(is.na(thelarche(dU)), 12, thelarche(dU)) + 0.5, 6.3), 11.5)
      dT <- sim_cpp(p, ev_leup(11.25, 84, a, a + 5), years = 17, delta = 2)
      data.frame(MK0 = mk, onset = thelarche(dU),
                 noRx = adult_height(dU), Rx = adult_height(dT))
    }))
    res$gain <- res$Rx - res$noRx
    ggplot(res, aes(onset)) +
      geom_hline(yintercept = 0, linetype = 3) +
      geom_line(aes(y = gain), colour = PAL["treated"], linewidth = 1) +
      geom_point(aes(y = gain), colour = PAL["treated"], size = 2.4) +
      geom_line(aes(y = (noRx - min(res$noRx))), colour = PAL["untreated"],
                linetype = 2) +
      annotate("text", x = max(res$onset, na.rm = TRUE), y = 0.4,
               label = "치료 이득 = 0", hjust = 1, size = 3.4) +
      labs(x = "무치료 시 사춘기 시작 연령 (yr)",
           y = "성인키 이득 (cm, 실선) / 무치료 최종키 상대값 (점선)") + THEME
  })

  # ---- Tab 2 -----------------------------------------------------------
  output$plt_pk <- renderPlot({
    s <- sims()
    cols <- c("CPTOT", "CANTo", "FAo", "FXo", "RS", "So")
    labs <- c("작용제 농도 (ng/mL)", "길항제 농도 (ng/mL)",
              "작용제 점유율 fa", "길항제 점유율 fx",
              "감작 수용체 RS", "성선자극세포 자극 S")
    d <- bind_rows(mutate(s$treated[, c("CA_out", cols)], arm = "treated"),
                   mutate(s$untreated[, c("CA_out", cols)], arm = "untreated")) %>%
      pivot_longer(all_of(cols), names_to = "var", values_to = "value")
    d$var <- factor(d$var, levels = cols, labels = labs)
    ggplot(d, aes(CA_out, value, colour = arm)) +
      geom_line(linewidth = 0.6) +
      facet_wrap(~var, scales = "free_y", ncol = 2) +
      scale_colour_manual(values = PAL, name = NULL) +
      labs(x = "역연령 (yr)", y = NULL) + THEME
  })

  # ---- Tab 3 -----------------------------------------------------------
  output$plt_hpg <- renderPlot({
    panel3(c("KND", "PULS", "LH", "FSH", "E2", "INHB"),
           c("KNDy 구동 (-)", "GnRH 박동 빈도 (pulses/24h)",
             "LH (IU/L)", "FSH (IU/L)", "에스트라디올 (pg/mL)",
             "인히빈 B (pg/mL)"))
  })

  # ---- Tab 4 -----------------------------------------------------------
  output$plt_growth <- renderPlot({
    s <- sims()
    cols <- c("HT", "GVo", "IGF1", "GPRES", "BA", "DBADCA")
    labs <- c("신장 (cm)", "성장속도 (cm/yr)", "IGF-1 (ng/mL)",
              "성장판 예비능 GPRES (-)", "골연령 (yr)", "ΔBA/ΔCA (-)")
    d <- long3(cols); d$var <- factor(d$var, levels = cols, labels = labs)
    ref <- data.frame(var = factor("골연령 (yr)", levels = labs),
                      x = range(s$treated$CA_out))
    ggplot(d, aes(CA_out, value, colour = arm)) +
      geom_line(linewidth = 0.7) +
      geom_abline(data = ref, aes(slope = 1, intercept = 0), linetype = 3,
                  inherit.aes = FALSE) +
      facet_wrap(~var, scales = "free_y", ncol = 2) +
      scale_colour_manual(values = PAL, name = NULL) +
      labs(x = "역연령 (yr)", y = NULL,
           caption = "골연령 패널의 점선은 BA = CA (정상 기준선)") + THEME
  })

  # ---- Tab 5 -----------------------------------------------------------
  output$tbl_endpoints <- renderTable({
    s <- sims()
    male <- s$male
    f <- function(d) c(
      `최종 성인키 (cm)` = adult_height(d),
      `골성숙 완료 연령 (yr)` = fusion_age(d),
      `초경 연령 (yr)` = menarche(d),
      `최종 Tanner 유방 단계` = tail(d$BST, 1),
      `치료 종료 시 BP-PAH (cm)` = pah_bp(at_age(d, min(input$start + input$dur,
                                                       max(d$CA_out)), "HT"),
                                          at_age(d, min(input$start + input$dur,
                                                        max(d$CA_out)), "BA"),
                                          male),
      `최고 심리사회 지표 (0-10)` = max(d$QOL),
      `최고 열성홍조 지표 (0-10)` = max(d$HF),
      `BMD Z 최저치` = min(d$BMDZ))
    m <- cbind(`치료` = f(s$treated), `무치료` = f(s$untreated),
               `정상 기준` = f(s$normal))
    out <- as.data.frame(round(m, 2))
    out["성인키 이득 (cm)", ] <- c(
      round(adult_height(s$treated) - adult_height(s$untreated), 2), 0,
      round(adult_height(s$normal) - adult_height(s$untreated), 2))
    out
  }, rownames = TRUE, digits = 2)

  output$plt_height <- renderPlot({
    s <- sims()
    d <- long3(c("HT"))
    bp <- bind_rows(
      data.frame(CA_out = s$treated$CA_out, arm = "treated",
                 value = pah_bp(s$treated$HT, s$treated$BA, s$male)),
      data.frame(CA_out = s$untreated$CA_out, arm = "untreated",
                 value = pah_bp(s$untreated$HT, s$untreated$BA, s$male)))
    ggplot() +
      geom_line(data = d, aes(CA_out, value, colour = arm), linewidth = 0.8) +
      geom_line(data = bp, aes(CA_out, value, colour = arm), linetype = 2) +
      scale_colour_manual(values = PAL, name = NULL) +
      labs(x = "역연령 (yr)", y = "신장 (실선) / BP 예측 성인키 (점선), cm") +
      THEME
  })

  # ---- Tab 6 -----------------------------------------------------------
  signflip <- reactive({
    p <- base_pars()
    dU <- sim_cpp(p, NULL, years = 17, delta = 2)
    hU <- adult_height(dU)
    starts <- seq(6.4, 11.2, by = 0.4)
    do.call(rbind, lapply(starts, function(a) {
      dT <- sim_cpp(p, ev_leup(11.25, 84, a, a + 6), years = 17, delta = 2)
      data.frame(start_CA = a, start_BA = at_age(dU, a, "BA"),
                 GPRES = at_age(dU, a, "GPRES"),
                 gain = adult_height(dT) - hU,
                 GV_on_Rx = gv_mean(dT, a + 0.3, min(a + 2, 16.4)))
    }))
  })

  output$plt_signflip <- renderPlot({
    r <- signflip()
    ggplot(r, aes(start_BA, gain)) +
      geom_hline(yintercept = 0, linetype = 3) +
      geom_line(linewidth = 1, colour = PAL["treated"]) +
      geom_point(size = 2.4, colour = PAL["treated"]) +
      geom_vline(xintercept = 12, linetype = 2, colour = PAL["untreated"]) +
      annotate("text", x = 12.05, y = max(r$gain), hjust = 0, size = 3.6,
               label = "임상 합의: 골연령 12세 이후 이득 없음") +
      labs(x = "치료 시작 시점의 골연령 (yr)",
           y = "무치료 대비 성인키 이득 (cm)") + THEME
  })

  output$tbl_signflip <- renderTable({ signflip() }, digits = 2)

  # ---- Tab 7 -----------------------------------------------------------
  output$plt_monitor <- renderPlot({
    s <- sims()
    cols <- c("GVo", "DBADCA", "LH", "UTV")
    labs <- c("성장속도 (cm/yr) — 믿을 수 없음",
              "ΔBA/ΔCA — 신뢰할 수 있음", "기저 LH (IU/L)",
              "자궁 부피 (mL)")
    d <- long3(cols); d$var <- factor(d$var, levels = cols, labels = labs)
    ggplot(d, aes(CA_out, value, colour = arm)) +
      geom_line(linewidth = 0.7) +
      facet_wrap(~var, scales = "free_y", ncol = 2) +
      scale_colour_manual(values = PAL, name = NULL) +
      labs(x = "역연령 (yr)", y = NULL) + THEME
  })

  output$tbl_monitor <- renderTable({
    s <- sims()
    lo <- input$start + 0.35; hi <- input$start + input$dur
    f <- function(d, col) mean(d[[col]][d$CA_out >= lo & d$CA_out <= hi])
    data.frame(
      surrogate = c("기저 LH (IU/L)", "E2 (pg/mL)", "자궁 부피 (mL)",
                    "성장속도 (cm/yr)", "ΔBA/ΔCA", "Tanner 유방", "IGF-1"),
      treated = c(f(s$treated, "LH"), f(s$treated, "E2"), f(s$treated, "UTV"),
                  f(s$treated, "GVo"), f(s$treated, "DBADCA"),
                  f(s$treated, "BST"), f(s$treated, "IGF1")),
      untreated = c(f(s$untreated, "LH"), f(s$untreated, "E2"),
                    f(s$untreated, "UTV"), f(s$untreated, "GVo"),
                    f(s$untreated, "DBADCA"), f(s$untreated, "BST"),
                    f(s$untreated, "IGF1")),
      check.names = FALSE)
  }, digits = 3)

  # ---- Tab 8 -----------------------------------------------------------
  forms <- reactive({
    p <- base_pars()
    a0 <- input$start; a1 <- input$start + input$dur
    regs <- list(
      "류프롤리드 3.75 q28d"   = ev_leup(3.75, 28, a0, a1),
      "류프롤리드 11.25 q12wk" = ev_leup(11.25, 84, a0, a1),
      "트립토렐린 22.5 q24wk"  = ev_trip(22.5, 168, a0, a1),
      "히스트렐린 임플란트"    = ev_hist(a0, a1),
      "나파렐린 비강분무"      = ev_naf(1800, a0, a1),
      "3.75 mg 42일마다(지연)" = ev_leup(3.75, 42, a0, a1),
      "1.875 mg q28d(저용량)"  = ev_leup(1.875, 28, a0, a1),
      "GnRH 길항제 18 mg"      = ev_antag(18, 28, a0, a1))
    hU <- adult_height(sim_cpp(p, NULL, years = 16.5, delta = 2))
    list(hU = hU, tab = do.call(rbind, lapply(names(regs), function(nm) {
      d <- sim_cpp(p, regs[[nm]], years = 16.5, delta = 1)
      data.frame(formulation = nm,
                 mean_Cp = mean(d$CPTOT[d$CA_out >= a0 + 0.25 & d$CA_out <= a1]),
                 mean_fa = mean(d$FAo[d$CA_out >= a0 + 0.25 & d$CA_out <= a1]),
                 pct_LH_esc = pct_time(d, a0 + 0.25, a1, "LH", 0.5),
                 pct_E2_esc = pct_time(d, a0 + 0.25, a1, "E2", 10),
                 adult_height = adult_height(d),
                 gain = adult_height(d) - hU)
    })))
  })

  output$plt_forms <- renderPlot({
    r <- forms()$tab
    ggplot(r, aes(reorder(formulation, gain), gain, fill = pct_E2_esc)) +
      geom_col() + coord_flip() +
      scale_fill_gradient(low = "#2471a3", high = "#c0392b",
                          name = "E2 > 10 pg/mL 인 시간 (%)") +
      labs(x = NULL, y = "무치료 대비 성인키 이득 (cm)") + THEME
  })

  output$tbl_forms <- renderTable({ forms()$tab }, digits = 3)

  # ---- Tab 9 -----------------------------------------------------------
  output$plt_biomarker <- renderPlot({
    panel3(c("UTV", "BST", "ENDO", "DHEAS", "TESTO", "FOL"),
           c("자궁 부피 (mL)", "Tanner 유방 단계 (1-5)",
             "자궁내막 상태 (-)", "DHEAS (µg/dL) — 억제되지 않음",
             "테스토스테론 (ng/dL)", "난포 기능 질량 (-)"))
  })

  # ---- Tab 10 ----------------------------------------------------------
  output$plt_safety <- renderPlot({
    panel3(c("BMDZ", "BMIZ", "HF", "QOL", "GH", "CUME2"),
           c("요추 BMD Z-점수", "BMI Z-점수", "열성홍조 지표 (0-10)",
             "심리사회 고통 지표 (0-10)", "GH 박동 진폭 (단위)",
             "누적 E2 노출 (pg/mL·yr)"))
  })
}

shinyApp(ui, server)
