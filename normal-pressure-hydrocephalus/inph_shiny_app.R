## =============================================================================
##  inph_shiny_app.R
##  Interactive dashboard for the iNPH QSP model
##  정상압 수두증 (iNPH) QSP 모델 — 인터랙티브 대시보드
## =============================================================================
##
##  13 tabs:
##    1  환자 프로파일        patient phenotype and the hydraulic baseline
##    2  CSF 역학             Davson / Marmarou — pressure, compliance, amplitude
##    3  션트 수력학          valve + gravitational unit + posture, resolved
##    4  밸브 적정 지도       the titration map (the central clinical result)
##    5  압력-용적 곡선       P-V curve, pulse amplitude, RAP
##    6  뇌실 형태            Evans index, callosal angle, DESH, two-clock creep
##    7  백질 손상            recoverable vs permanent white matter
##    8  임상 삼징            gait / cognition / continence, fast + slow arms
##    9  진단 검사            tap test and external lumbar drain simulator
##   10  약물 PK/PD           acetazolamide (saturable RBC binding) and others
##   11  CSF 바이오마커       the dilution/flux decomposition
##   12  시나리오 비교        all 21 scenarios side by side
##   13  안전성               subdural collection, headache, acidosis
##
##  Run:  shiny::runApp("inph_shiny_app.R")
##  Requires inph_mrgsolve_model.R in the same directory.
## =============================================================================

library(shiny)
library(mrgsolve)
library(dplyr)
library(tidyr)
library(ggplot2)

source("inph_mrgsolve_model.R", local = FALSE)

THEME <- theme_bw(base_size = 12) +
  theme(legend.position = "bottom",
        strip.background = element_rect(fill = "#eef3f8"),
        plot.title = element_text(face = "bold"))

## Colour scale used consistently across tabs
PAL <- c("#2c6fb5", "#c0504d", "#4e8b4e", "#c98a2e", "#7b4fa8", "#2f8177",
         "#b0518a", "#8a8a3f")

fmt <- function(x, d = 2) formatC(x, format = "f", digits = d)

# --------------------------------------------------------------------------- #
#  UI
# --------------------------------------------------------------------------- #
ui <- fluidPage(
  titlePanel(paste("정상압 수두증 (Idiopathic Normal Pressure Hydrocephalus)",
                   "QSP 대시보드")),
  tags$p(style = "color:#555;margin-top:-8px",
         paste("평균 뇌압은 정의상 정상이다. 병태생리의 실제 변수는",
               "유출저항(Rout) · 순응도(C) · 박동성(AMP)이며,",
               "치료의 안전역은 생물학이 아니라 정수압(hydrostatics)이 결정한다.")),

  sidebarLayout(
    sidebarPanel(
      width = 3,
      h4("환자 표현형 (Phenotype)"),
      sliderInput("Rout_init", "CSF 유출저항 Rout [mmHg/(mL/min)]",
                  6, 26, 19, 0.5),
      ## E1 is defined against the P-V reference pressure P0_marm = -2 mmHg;
      ## healthy ~0.156, iNPH ~0.222, post-decompression floor 0.162.
      sliderInput("E1_init", "탄성계수 E1 [1/mL] (순응도 반비례)",
                  0.13, 0.34, 0.222, 0.005),
      sliderInput("Vv0", "뇌실 용적 Vv [mL]", 30, 160, 95, 5),
      sliderInput("WM_init", "백질 기능 보전도 WMint", 0.35, 0.95, 0.68, 0.01),
      sliderInput("WMperm_init", "영구 축삭 손실 WMperm", 0.0, 0.45, 0.10, 0.01),
      sliderInput("APOE", "APOE 계수 (1.0 = e3/e3, 2.2 = e4)",
                  1.0, 2.5, 1.0, 0.1),
      sliderInput("atrophy", "피질 위축 (경막하 공간 크기)", 0, 1, 0.55, 0.05),
      sliderInput("f_up", "하루 중 직립 비율 f_up", 0.2, 0.9, 0.60, 0.05),
      sliderInput("comorb", "비수두증성 보행 장애 부담", 0, 0.5, 0, 0.05),

      hr(),
      h4("션트 하드웨어 (Hardware)"),
      checkboxInput("shunt", "션트 삽입 (VP shunt)", FALSE),
      sliderInput("Popen_cm", "밸브 개방압 Popen [cmH2O]", 0, 30, 10, 1),
      sliderInput("Ggrav_cm", "중력식 보조기 [cmH2O] (직립 시에만 작동)",
                  0, 40, 0, 5),
      sliderInput("asd_eff", "막형 anti-siphon 효율", 0, 1, 0, 0.05),
      sliderInput("Rsh", "션트 저항 Rsh [mmHg/(mL/min)]", 1, 40, 2.5, 0.5),
      sliderInput("Hcol_cm", "정수압 컬럼 높이 [cmH2O]", 20, 60, 45, 1),
      checkboxInput("etv", "ETV (제3뇌실 천공술)", FALSE),

      hr(),
      h4("약물 (Drugs)"),
      selectInput("az_dose", "Acetazolamide (BID)",
                  c("없음" = "0", "250 mg" = "250", "500 mg" = "500",
                    "1000 mg" = "1000"), "0"),
      checkboxInput("melatonin", "Melatonin 2 mg 취침 전", FALSE),
      checkboxInput("solifenacin", "Solifenacin 5 mg", FALSE),
      checkboxInput("donepezil", "Donepezil 10 mg", FALSE),
      checkboxInput("antithrombotic", "항혈전제 복용", FALSE),

      hr(),
      sliderInput("tend", "시뮬레이션 기간 [개월]", 3, 60, 24, 3),
      actionButton("go", "시뮬레이션 실행", class = "btn-primary"),
      tags$p(style = "font-size:11px;color:#777;margin-top:10px",
             paste("교육/연구 목적의 모델입니다. 임상 의사결정에",
                   "사용할 수 없습니다."))
    ),

    mainPanel(
      width = 9,
      tabsetPanel(
        id = "tabs", type = "tabs",

        tabPanel("1 환자 프로파일",
                 h4("기저 수력학 (closed-form flow balance)"),
                 tableOutput("tbl_baseline"),
                 h4("건강한 75세 대조군과의 비교"),
                 tableOutput("tbl_vs_control"),
                 plotOutput("plt_profile", height = "320px"),
                 tags$small(paste("Rout는 2배 이상 상승하지만 평균 뇌압은",
                                  "7-15 mmHg 정상 범위를 벗어나지 않는다.",
                                  "박동압(AMP)은 그 손실을 증폭한다."))),

        tabPanel("2 CSF 역학",
                 plotOutput("plt_hydro", height = "560px"),
                 tags$small(paste("Davson: ICP = Pss + If·Rout.",
                                  "Marmarou: C = 1/(E1·(P−P0)).",
                                  "AMP = ΔVp/C."))),

        tabPanel("3 션트 수력학",
                 h4("자세별 분해 (supine vs upright)"),
                 tableOutput("tbl_posture"),
                 plotOutput("plt_posture", height = "380px"),
                 tags$small(paste("직립 시 뇌실-복강 정수압 컬럼이 구동압에",
                                  "그대로 더해진다. 중력식 보조기는 직립에서만",
                                  "개방압을 올려 이 항을 상쇄한다."))),

        tabPanel("4 밸브 적정 지도",
                 h4("개방압 스윕 — 중력식 보조기 유무별"),
                 actionButton("run_map", "적정 지도 계산 (수십 초 소요)"),
                 plotOutput("plt_map", height = "560px"),
                 tableOutput("tbl_map")),

        tabPanel("5 압력-용적 곡선",
                 plotOutput("plt_pv", height = "420px"),
                 plotOutput("plt_amp", height = "300px")),

        tabPanel("6 뇌실 형태",
                 plotOutput("plt_morph", height = "520px"),
                 tags$small(paste("성공적인 션트 후에도 Evans index는 거의",
                                  "변하지 않는다 — 가소성(plastic) 확장이",
                                  "되돌아오지 않기 때문이다."))),

        tabPanel("7 백질 손상",
                 plotOutput("plt_wm", height = "520px"),
                 tags$small(paste("WMint는 회복 가능한 풀, WMperm은",
                                  "비가역 풀이다. 션트는 전자를 되돌리고",
                                  "후자의 시계를 멈출 뿐이다."))),

        tabPanel("8 임상 삼징",
                 plotOutput("plt_triad", height = "560px"),
                 tableOutput("tbl_triad")),

        tabPanel("9 진단 검사",
                 h4("Tap test / 지속 배액 시뮬레이터"),
                 fluidRow(
                   column(4, sliderInput("tap_vol", "천자 제거량 [mL]",
                                         0, 60, 40, 5)),
                   column(4, sliderInput("eld_h", "지속 배액 기간 [시간]",
                                         0, 120, 72, 12)),
                   column(4, sliderInput("thr", "판정 역치 [m/s]",
                                         0.02, 0.20, 0.10, 0.01))),
                 plotOutput("plt_dx", height = "440px"),
                 tableOutput("tbl_dx"),
                 tags$small(paste("τ = Rout·C ≈ 6분은 선형화 값이며 40 mL",
                                  "천자에는 적용되지 않는다. 정맥동압 아래로",
                                  "내려가면 흡수가 0이 되므로 회복은 생성률",
                                  "(0.35 mL/min)에 제한되어 약 114분이 걸린다.",
                                  "빠른 성분은 비대칭이다 — 충전 29분, 방전",
                                  "21.6시간."))),

        tabPanel("10 약물 PK/PD",
                 plotOutput("plt_pk", height = "440px"),
                 tableOutput("tbl_pk"),
                 tags$small(paste("Acetazolamide는 적혈구 탄산탈수효소에",
                                  "포화 결합한다 — 용량에 따라 겉보기",
                                  "반감기가 길어지는 이유. PD는 EC50이 낮아",
                                  "250 mg BID 이상에서 평탄역이다."))),

        tabPanel("11 CSF 바이오마커",
                 plotOutput("plt_bio", height = "480px"),
                 tableOutput("tbl_bio"),
                 tags$small(paste("c = J(brain→CSF) / [(If+Qsh)/Vcsf + kdeg].",
                                  "iNPH에서는 분자(글림파틱 유출)와 분모의",
                                  "회전율이 동시에 감소하므로 Aβ42와 p-tau가",
                                  "모두 낮게 측정된다."))),

        tabPanel("12 시나리오 비교",
                 actionButton("run_scen", "21개 시나리오 실행 (수 분 소요)"),
                 tableOutput("tbl_scen"),
                 plotOutput("plt_scen", height = "460px")),

        tabPanel("13 안전성",
                 plotOutput("plt_safety", height = "520px"),
                 tableOutput("tbl_safety"),
                 tags$small(paste("경막하 수액낭의 발생률은 이 모델에서",
                                  "서열척도(ordinal)이다 — k_haz가 보정되지",
                                  "않았기 때문이며, 절대 발생률로 읽어서는",
                                  "안 된다.")))
      )
    )
  )
)

# --------------------------------------------------------------------------- #
#  Server
# --------------------------------------------------------------------------- #
server <- function(input, output, session) {

  ## ---- assemble the parameter overrides from the sidebar ---------------- ##
  overrides <- reactive({
    list(Rout_init = input$Rout_init, atrophy = input$atrophy,
         f_up = input$f_up, comorb = input$comorb, APOE = input$APOE,
         shunt = as.numeric(input$shunt), Popen_cm = input$Popen_cm,
         Ggrav_cm = input$Ggrav_cm, asd_eff = input$asd_eff,
         Rsh = input$Rsh, Hcol_cm = input$Hcol_cm,
         etv = as.numeric(input$etv),
         antithrombotic = as.numeric(input$antithrombotic))
  })

  ## The phenotype sliders that are INITIAL CONDITIONS rather than parameters
  ## have to be pushed into the init vector after inph_init() has built it.
  patient_now <- reactive({
    p <- as.list(param(mod))
    p[names(overrides())] <- overrides()
    y <- inph_init(p, healthy = FALSE)
    y["Vv"] <- input$Vv0
    y["Vplast"] <- p$f_plastic * (input$Vv0 - p$Vv_norm)
    y["E1"] <- input$E1_init
    y["WMint"] <- input$WM_init
    y["WMperm"] <- input$WMperm_init
    ## re-seat the clinical scores and the references on the edited state
    yl <- as.list(y)
    hy <- inph_hydro(yl, p)
    wex <- max(0, yl$Wpv - p$Wpv_norm)
    CBF <- p$CBF0 * (yl$Autoreg / p$Autoreg_norm) /
      (1 + p$kCBF_W * wex + p$kCBF_P * hy$dPtmP)
    perf <- 0.55 + 0.45 * min(1, CBF / p$CBF0)
    y["Gslow"] <- max(0.05, p$G_max * yl$WMint^p$G_wm_pow * perf *
                        (1 - 0.25 * yl$Ab_plq) * (1 - p$comorb))
    taub <- min(0.4, 0.4 * max(0, yl$Tau_isf - 1) / 3)
    y["Cslow"] <- max(5, p$MMSE_max - p$kcog_wm * (1 - yl$WMint) -
                        p$kcog_ad * min(1, 0.6 * yl$Ab_plq + taub))
    y["Urin"] <- max(0, min(p$Ur_max, p$kur_wm * (1 - yl$WMint)))
    p$AMP_ref <- hy$AMPday
    p$Pday_ref <- hy$Pday
    p$E1_floor_eff <- min(p$E1_floor, yl$E1)
    list(p = p, init = y, hyd = hy, CBF = CBF)
  })

  ## ---- drug events ------------------------------------------------------ ##
  drug_events <- reactive({
    tend <- input$tend * 30.4
    e <- NULL
    add <- function(a, b) if (is.null(a)) b else c(a, b)
    azd <- as.numeric(input$az_dose)
    if (azd > 0)
      e <- add(e, dose_ev("Aaz_g", azd, bioav = 0.9, ii = 0.5,
                          addl = ceiling(tend / 0.5) - 1))
    if (input$melatonin)
      e <- add(e, dose_ev("Ame_g", 2, bioav = 0.15, ii = 1,
                          addl = ceiling(tend) - 1))
    if (input$solifenacin)
      e <- add(e, dose_ev("Aso_g", 5, bioav = 0.9, ii = 1,
                          addl = ceiling(tend) - 1))
    if (input$donepezil)
      e <- add(e, dose_ev("Ado_g", 10, bioav = 1.0, ii = 1,
                          addl = ceiling(tend) - 1))
    e
  })

  ## ---- the main simulation --------------------------------------------- ##
  sim <- eventReactive(input$go, {
    pt <- patient_now()
    tend <- input$tend * 30.4
    m <- mod %>% param(pt$p) %>% init(pt$init)
    ev0 <- drug_events()
    out <- if (is.null(ev0)) mrgsim(m, end = tend, delta = 1, hmax = 0.25)
           else mrgsim(m, data = ev0, end = tend, delta = 1, hmax = 0.25)
    as_tibble(out)
  }, ignoreNULL = FALSE)

  ## ---- reference: the untreated / healthy comparators ------------------ ##
  control <- reactive({
    p <- as.list(param(mod))
    r <- inph_refs(p, healthy = TRUE)
    p[names(r)] <- r
    y <- inph_init(p, healthy = TRUE)
    hy <- inph_hydro(as.list(y), p)
    list(p = p, init = y, hyd = hy)
  })

  ## ===================== TAB 1 — patient profile ======================== ##
  output$tbl_baseline <- renderTable({
    pt <- patient_now(); h <- pt$hyd
    tibble(
      항목 = c("Rout [mmHg/(mL/min)]", "E1 [1/mL]",
               "평균 ICP 누운 자세 [mmHg]", "평균 ICP 직립 [mmHg]",
               "일평균 ICP [mmHg]", "순응도 C (누운 자세) [mL/mmHg]",
               "ICP 박동압 AMP [mmHg]", "박동성 transmantle 구배 [mmHg]",
               "평균 transmantle 구배 [mmHg]",
               "션트 유량 [mL/day]", "CSF 생성 [mL/day]",
               "뇌실주위 백질 관류 [mL/100g/min]",
               "τ = Rout·C (선형화) [min]",
               "40 mL 천자 회복 최소 시간 [min]"),
      값 = c(fmt(pt$init[["Rout"]], 1), fmt(pt$init[["E1"]], 3),
             fmt(h$Ps), fmt(h$Pu), fmt(h$Pday), fmt(h$Csup, 3),
             fmt(h$AMP), fmt(h$dPtmP, 3), fmt(h$dPtmM, 3),
             fmt(h$Qsh, 0), fmt(h$Ifd, 0), fmt(pt$CBF, 1),
             ## Rout [mmHg/(mL/min)] x C [mL/mmHg] is already in minutes
             fmt(pt$init[["Rout"]] * h$Csup, 1),
             ## but the tap-test recovery is PRODUCTION-limited, not
             ## resistance-limited: below sinus pressure nothing is absorbed
             fmt(40 / pt$p$If0, 0))
    )
  }, striped = TRUE)

  output$tbl_vs_control <- renderTable({
    pt <- patient_now(); ct <- control()
    a <- ct$hyd; b <- pt$hyd
    rat <- function(x, y) if (abs(x) > 1e-9) paste0(fmt(y / x), "x") else "n/a"
    tibble(
      항목 = c("Rout", "평균 ICP (누운 자세)", "박동압 AMP",
               "박동성 transmantle 구배", "Evans index"),
      `건강 75세` = c(fmt(ct$init[["Rout"]], 1), fmt(a$Ps), fmt(a$AMP),
                     fmt(a$dPtmP, 3),
                     fmt(0.20 + 0.0035 * (ct$init[["Vv"]] - 25), 3)),
      `환자` = c(fmt(pt$init[["Rout"]], 1), fmt(b$Ps), fmt(b$AMP),
                fmt(b$dPtmP, 3),
                fmt(0.20 + 0.0035 * (pt$init[["Vv"]] - 25), 3)),
      `비` = c(rat(ct$init[["Rout"]], pt$init[["Rout"]]),
              rat(a$Ps, b$Ps), rat(a$AMP, b$AMP),
              rat(a$dPtmP, b$dPtmP),
              rat(0.20 + 0.0035 * (ct$init[["Vv"]] - 25),
                  0.20 + 0.0035 * (pt$init[["Vv"]] - 25)))
    )
  }, striped = TRUE)

  output$plt_profile <- renderPlot({
    ## mean ICP and pulse amplitude as functions of Rout, at this patient's E1
    pt <- patient_now(); p <- pt$p
    grid <- lapply(seq(6, 26, 0.5), function(R) {
      y <- as.list(pt$init); y$Rout <- R
      h <- inph_hydro(y, p)
      tibble(Rout = R, `평균 ICP (누운 자세)` = h$Ps, `박동압 AMP` = h$AMP)
    }) %>% bind_rows() %>% pivot_longer(-Rout)
    ggplot(grid, aes(Rout, value, colour = name)) +
      geom_line(linewidth = 1) +
      geom_vline(xintercept = pt$init[["Rout"]], linetype = 2,
                 colour = "#555555") +
      annotate("rect", xmin = -Inf, xmax = Inf, ymin = 7, ymax = 15,
               alpha = 0.08, fill = "#2c6fb5") +
      annotate("text", x = 7.5, y = 15.6, hjust = 0, size = 3.4,
               label = "정상 ICP 범위 7-15 mmHg") +
      geom_hline(yintercept = 4, linetype = 3, colour = "#c0504d") +
      annotate("text", x = 24, y = 4.4, hjust = 1, size = 3.4,
               colour = "#c0504d", label = "AMP 4 mmHg 역치") +
      scale_colour_manual(values = PAL) +
      labs(x = "CSF 유출저항 Rout [mmHg/(mL/min)]", y = "mmHg", colour = NULL,
           title = "질환은 자기 이름이 가리키는 변수에 거의 나타나지 않는다") +
      THEME
  })

  ## ===================== TAB 2 — CSF hydrodynamics ====================== ##
  output$plt_hydro <- renderPlot({
    d <- sim()
    d %>% select(time, ICP_sup, ICP_up, ICP_day, AMP, Cspine, dPtm_pulse,
                 dPtm_mean, Qsh_day, If_day) %>%
      pivot_longer(-time) %>%
      mutate(name = factor(name, levels = c(
        "ICP_sup", "ICP_up", "ICP_day", "AMP", "Cspine", "dPtm_pulse",
        "dPtm_mean", "Qsh_day", "If_day"))) %>%
      ggplot(aes(time / 30.4, value)) +
      geom_line(colour = PAL[1], linewidth = 0.9) +
      facet_wrap(~name, scales = "free_y", ncol = 3) +
      labs(x = "개월", y = NULL, title = "뇌척수액 역학 시간 경과") + THEME
  })

  ## ===================== TAB 3 — shunt hydraulics ======================= ##
  posture_table <- reactive({
    pt <- patient_now(); p <- pt$p; y <- as.list(pt$init)
    cm <- function(x) x / 1.36
    rows <- lapply(c(FALSE, TRUE), function(up) {
      pp <- p; pp$f_up <- if (up) 1 else 0
      h <- inph_hydro(y, pp)
      tibble(자세 = if (up) "직립 (upright)" else "누운 자세 (supine)",
             `ICP [mmHg]` = fmt(if (up) h$Pu else h$Ps),
             `정수압 컬럼 [mmHg]` = fmt(if (up)
               cm(p$Hcol_cm) * (1 - p$asd_eff) else 0),
             `유효 개방압 [mmHg]` = fmt(cm(p$Popen_cm) +
                                         if (up) cm(p$Ggrav_cm) else 0),
             `말단 압력 [mmHg]` = fmt(cm(if (up) p$Pd_up_cm else p$Pd_sup_cm)),
             `션트 유량 [mL/day]` = fmt(h$Qsh, 0),
             `생성량 대비` = fmt(h$Qsh / h$Ifd))
    })
    bind_rows(rows)
  })

  output$tbl_posture <- renderTable(posture_table(), striped = TRUE)

  output$plt_posture <- renderPlot({
    pt <- patient_now(); p <- pt$p; y <- as.list(pt$init)
    grid <- lapply(seq(0, 30, 1), function(Po) {
      bind_rows(lapply(c(0, 30), function(G) {
        pp <- p; pp$shunt <- 1; pp$Popen_cm <- Po; pp$Ggrav_cm <- G
        h <- inph_hydro(y, pp)
        tibble(Popen = Po,
               arm = if (G > 0) "중력식 보조기 30 cmH2O" else "보호 없음",
               `일평균 ICP` = h$Pday, `직립 ICP` = h$Pu,
               `배액량 mL/day` = h$Qsh)
      }))
    }) %>% bind_rows() %>% pivot_longer(-c(Popen, arm))
    ggplot(grid, aes(Popen, value, colour = arm)) +
      geom_line(linewidth = 1) +
      geom_hline(yintercept = 0, linetype = 3) +
      facet_wrap(~name, scales = "free_y") +
      scale_colour_manual(values = PAL[c(2, 3)]) +
      labs(x = "밸브 개방압 [cmH2O]", y = NULL, colour = NULL,
           title = "정수압이 안전역을 결정한다") + THEME
  })

  ## ===================== TAB 4 — titration map ========================== ##
  map_data <- eventReactive(input$run_map, {
    withProgress(message = "적정 지도 계산 중...", {
      pt <- patient_now()
      tend <- min(input$tend * 30.4, 730)
      rows <- list()
      grid <- expand.grid(G = c(0, 30),
                          Po = c(2, 4, 6, 8, 10, 12, 14, 16, 20, 24, 28))
      for (i in seq_len(nrow(grid))) {
        incProgress(1 / nrow(grid))
        p <- pt$p; p$shunt <- 1
        p$Popen_cm <- grid$Po[i]; p$Ggrav_cm <- grid$G[i]
        m <- mod %>% param(p) %>% init(pt$init)
        o <- as_tibble(mrgsim(m, end = tend, delta = 30, hmax = 0.25))
        rows[[i]] <- tibble(
          Ggrav = grid$G[i], Popen = grid$Po[i],
          ICP_day = last(o$ICP_day), ICP_up = last(o$ICP_up),
          AMP = last(o$AMP), dGait = last(o$gait) - first(o$gait),
          Vsdh = last(o$Vsdh), SDH_pct = 100 * last(o$SDH_inc),
          headache = last(o$Headx))
      }
      bind_rows(rows) %>%
        mutate(utility = dGait - 1.1 * SDH_pct / 100 - 0.03 * headache,
               arm = ifelse(Ggrav > 0, "중력식 보조기 30 cmH2O",
                            "보호 없음"))
    })
  })

  output$plt_map <- renderPlot({
    tm <- map_data()
    tm %>% select(arm, Popen, dGait, SDH_pct, Vsdh, utility) %>%
      pivot_longer(-c(arm, Popen)) %>%
      ggplot(aes(Popen, value, colour = arm)) +
      geom_line(linewidth = 1) + geom_point(size = 1.8) +
      facet_wrap(~name, scales = "free_y", ncol = 2) +
      scale_colour_manual(values = PAL[c(3, 2)]) +
      labs(x = "밸브 개방압 [cmH2O]", y = NULL, colour = NULL,
           title = "밸브 적정 지도 — 중력식 보조기는 최적점 자체를 이동시킨다") +
      THEME
  })

  output$tbl_map <- renderTable({
    tm <- map_data()
    tm %>% group_by(arm) %>%
      slice_max(utility, n = 1) %>%
      transmute(arm, `최적 개방압 [cmH2O]` = Popen,
                `보행 이득 [m/s]` = fmt(dGait, 3),
                `경막하 발생률 [%]` = fmt(SDH_pct, 1),
                `수액낭 [mL]` = fmt(Vsdh, 1),
                `일평균 ICP [mmHg]` = fmt(ICP_day)) %>% ungroup()
  }, striped = TRUE)

  ## ===================== TAB 5 — P-V curve ============================== ##
  output$plt_pv <- renderPlot({
    pt <- patient_now(); ct <- control()
    mk <- function(E1, lab) {
      P <- seq(-4, 40, 0.25)
      tibble(P = P, C = pmin(1 / (E1 * pmax(P - pt$p$P0_marm, 0.5)),
                             pt$p$C_max), arm = lab)
    }
    d <- bind_rows(mk(pt$init[["E1"]], sprintf("환자 E1 = %.3f",
                                               pt$init[["E1"]])),
                   mk(ct$init[["E1"]], sprintf("건강 대조 E1 = %.3f",
                                               ct$init[["E1"]])))
    ggplot(d, aes(P, C, colour = arm)) + geom_line(linewidth = 1) +
      geom_vline(xintercept = pt$hyd$Ps, linetype = 2) +
      scale_colour_manual(values = PAL[c(1, 3)]) +
      labs(x = "ICP [mmHg]", y = "순응도 C [mL/mmHg]", colour = NULL,
           title = "Marmarou 지수 압력-용적 곡선",
           subtitle = paste("점선 = 환자의 현재 작동점.",
                            "낮은 순응도가 같은 동맥 박동을 더 큰",
                            "압력 파형으로 바꾼다.")) + THEME
  })

  output$plt_amp <- renderPlot({
    d <- sim()
    d %>% select(time, AMP, AMP_day, Cspine) %>% pivot_longer(-time) %>%
      ggplot(aes(time / 30.4, value, colour = name)) +
      geom_line(linewidth = 0.9) +
      scale_colour_manual(values = PAL) +
      labs(x = "개월", y = NULL, colour = NULL,
           title = "박동압과 순응도의 시간 경과") + THEME
  })

  ## ===================== TAB 6 — morphology ============================= ##
  output$plt_morph <- renderPlot({
    d <- sim()
    d %>% select(time, Vv, Vsas, Vplast, EvansIdx, CallAngle, dPtm_pulse) %>%
      pivot_longer(-time) %>%
      ggplot(aes(time / 30.4, value)) +
      geom_line(colour = PAL[4], linewidth = 0.9) +
      facet_wrap(~name, scales = "free_y") +
      labs(x = "개월", y = NULL,
           title = "뇌실 형태 — 가소성 확장은 되돌아오지 않는다") + THEME
  })

  ## ===================== TAB 7 — white matter =========================== ##
  output$plt_wm <- renderPlot({
    d <- sim()
    d %>% select(time, WMint, WMperm, Myel, Wpv, AQ, CBF_pv, Astro, Micro) %>%
      pivot_longer(-time) %>%
      ggplot(aes(time / 30.4, value)) +
      geom_line(colour = PAL[3], linewidth = 0.9) +
      facet_wrap(~name, scales = "free_y", ncol = 4) +
      labs(x = "개월", y = NULL,
           title = "두 개의 시계: 회복 가능한 풀과 비가역 풀") + THEME
  })

  ## ===================== TAB 8 — clinical triad ========================= ##
  output$plt_triad <- renderPlot({
    d <- sim()
    d %>% select(time, gait, MMSE, urin, iNPHGS, Gfast, Gslow) %>%
      pivot_longer(-time) %>%
      ggplot(aes(time / 30.4, value)) +
      geom_line(colour = PAL[1], linewidth = 0.9) +
      facet_wrap(~name, scales = "free_y", ncol = 3) +
      labs(x = "개월", y = NULL,
           title = "Hakim 삼징 — 빠른 성분(압력 결합)과 느린 성분(조직 회복)") +
      THEME
  })

  output$tbl_triad <- renderTable({
    d <- sim()
    tibble(항목 = c("보행 속도 [m/s]", "MMSE", "요절박 지수", "iNPHGS"),
           시작 = c(fmt(first(d$gait), 3), fmt(first(d$MMSE), 1),
                   fmt(first(d$urin)), fmt(first(d$iNPHGS), 0)),
           종료 = c(fmt(last(d$gait), 3), fmt(last(d$MMSE), 1),
                   fmt(last(d$urin)), fmt(last(d$iNPHGS), 0)),
           변화 = c(fmt(last(d$gait) - first(d$gait), 3),
                   fmt(last(d$MMSE) - first(d$MMSE), 1),
                   fmt(last(d$urin) - first(d$urin)),
                   fmt(last(d$iNPHGS) - first(d$iNPHGS), 0)))
  }, striped = TRUE)

  ## ===================== TAB 9 — diagnostics ============================ ##
  dx_sim <- reactive({
    pt <- patient_now()
    m <- mod %>% param(pt$p) %>% init(pt$init)
    tap <- as_tibble(mrgsim(m, data = tap_ev(0.5, input$tap_vol),
                            end = 8, delta = 0.01, hmax = 0.02)) %>%
      mutate(arm = sprintf("Tap test %g mL", input$tap_vol))
    ph <- pt$p; ph$eld_rate <- 10 / 60
    me <- mod %>% param(ph) %>% init(pt$init)
    ## drain for input$eld_h hours, then stop
    tsw <- input$eld_h / 24
    o1 <- as_tibble(mrgsim(me, end = max(tsw, 0.01), delta = 0.01, hmax = 0.02))
    y2 <- unlist(o1[nrow(o1), names(pt$init)])
    m2 <- mod %>% param(pt$p) %>% init(y2)
    o2 <- as_tibble(mrgsim(m2, end = 8 - tsw, delta = 0.01, hmax = 0.02)) %>%
      mutate(time = time + tsw)
    eld <- bind_rows(o1, o2) %>%
      mutate(arm = sprintf("ELD 10 mL/h x %g h", input$eld_h))
    bind_rows(tap, eld)
  })

  output$plt_dx <- renderPlot({
    dx_sim() %>% select(time, arm, ICP_sup, AMP, gait) %>%
      pivot_longer(-c(time, arm)) %>%
      ggplot(aes(time * 24, value, colour = arm)) +
      geom_line(linewidth = 0.9) +
      facet_wrap(~name, scales = "free_y", ncol = 1) +
      scale_colour_manual(values = PAL[c(1, 2)]) +
      labs(x = "시간 [h]", y = NULL, colour = NULL,
           title = "진단적 섭동: 부피 계단(천자) vs 유량 계단(지속 배액)") +
      THEME
  })

  output$tbl_dx <- renderTable({
    d <- dx_sim()
    d %>% group_by(arm) %>%
      summarise(`기저 보행` = fmt(first(gait), 3),
                `최대 보행` = fmt(max(gait), 3),
                `최대 이득` = fmt(max(gait) - first(gait), 3),
                `역치 통과` = ifelse(max(gait) - first(gait) >= input$thr,
                                    "양성", "음성"),
                `최저 ICP` = fmt(min(ICP_sup)), .groups = "drop")
  }, striped = TRUE)

  ## ===================== TAB 10 — PK/PD ================================= ##
  output$plt_pk <- renderPlot({
    pt <- patient_now()
    bind_rows(lapply(c(125, 250, 500, 1000), function(dz) {
      m <- mod %>% param(pt$p) %>% init(pt$init)
      as_tibble(mrgsim(m, data = dose_ev("Aaz_g", dz, bioav = 0.9, ii = 0.5,
                                         addl = 39),
                       end = 20, delta = 0.05, hmax = 0.02)) %>%
        mutate(dose = sprintf("%g mg BID", dz))
    })) %>%
      select(time, dose, AZ_free, AZ_eff, AZ_rbcfrac, If_day, ICP_sup) %>%
      pivot_longer(-c(time, dose)) %>%
      ggplot(aes(time, value, colour = dose)) +
      geom_line(linewidth = 0.8) +
      facet_wrap(~name, scales = "free_y") +
      scale_colour_manual(values = PAL) +
      labs(x = "일", y = NULL, colour = NULL,
           title = "Acetazolamide: 포화 적혈구 결합 + PD 평탄역") + THEME
  })

  output$tbl_pk <- renderTable({
    pt <- patient_now()
    bind_rows(lapply(c(0, 125, 250, 500, 1000), function(dz) {
      m <- mod %>% param(pt$p) %>% init(pt$init)
      o <- if (dz == 0) as_tibble(mrgsim(m, end = 20, delta = 1, hmax = 0.25))
           else as_tibble(mrgsim(m, data = dose_ev("Aaz_g", dz, bioav = 0.9,
                                                   ii = 0.5, addl = 39),
                                 end = 20, delta = 1, hmax = 0.02))
      tibble(`용량 (BID)` = if (dz == 0) "없음" else sprintf("%g mg", dz),
             `유리 농도 [mg/L]` = fmt(last(o$AZ_free), 4),
             `CSF 생성 억제 [%]` = fmt(100 * last(o$AZ_eff), 1),
             `적혈구 결합 분율` = fmt(last(o$AZ_rbcfrac), 3),
             `CSF 생성 [mL/day]` = fmt(last(o$If_day), 0),
             `ICP 누운 자세 [mmHg]` = fmt(last(o$ICP_sup)),
             `HCO3 [mmol/L]` = fmt(last(o$HCO3), 1))
    }))
  }, striped = TRUE)

  ## ===================== TAB 11 — biomarkers ============================ ##
  output$plt_bio <- renderPlot({
    d <- sim()
    d %>% select(time, CSF_Ab42, CSF_pTau, CSF_NfL, CSF_LRG, turnover, AQ) %>%
      pivot_longer(-time) %>%
      ggplot(aes(time / 30.4, value)) +
      geom_line(colour = PAL[8], linewidth = 0.9) +
      facet_wrap(~name, scales = "free_y", ncol = 3) +
      labs(x = "개월", y = NULL,
           title = "CSF 바이오마커 — 희석과 유출 감소의 합성") + THEME
  })

  output$tbl_bio <- renderTable({
    d <- sim(); ct <- control()
    pc <- as.list(ct$p); yc <- as.list(ct$init)
    Vc <- yc$Vv + yc$Vsas
    koutc <- pc$If0 * 1440 / Vc + pc$kdeg_csf
    glyc <- yc$AQ / (1 + 0.35 * yc$Wpv)
    tibble(
      항목 = c("CSF 회전율 [/day]", "CSF Aβ42 [pg/mL]", "CSF p-tau [pg/mL]",
               "CSF NfL [pg/mL]", "CSF LRG [a.u.]", "AQP4 극성"),
      `건강 대조` = c(fmt(pc$If0 * 1440 / Vc, 3),
                     fmt(pc$kflux_ab * glyc * yc$Ab_isf / koutc / Vc, 0),
                     fmt(pc$kflux_tau * glyc * yc$Tau_isf / koutc / Vc, 1),
                     "-", "-", fmt(yc$AQ, 3)),
      시작 = c(fmt(first(d$turnover), 3), fmt(first(d$CSF_Ab42), 0),
              fmt(first(d$CSF_pTau), 1), fmt(first(d$CSF_NfL), 0),
              fmt(first(d$CSF_LRG)), "-"),
      종료 = c(fmt(last(d$turnover), 3), fmt(last(d$CSF_Ab42), 0),
              fmt(last(d$CSF_pTau), 1), fmt(last(d$CSF_NfL), 0),
              fmt(last(d$CSF_LRG)), "-"))
  }, striped = TRUE)

  ## ===================== TAB 12 — scenarios ============================= ##
  scen_data <- eventReactive(input$run_scen, {
    withProgress(message = "21개 시나리오 실행 중...", {
      res <- lapply(names(scen), function(nm) {
        incProgress(1 / length(scen), detail = nm)
        scen[[nm]]()
      })
      bind_rows(res)
    })
  })

  output$tbl_scen <- renderTable({
    summarise_scenarios(scen_data()) %>%
      mutate(across(where(is.numeric), ~ round(.x, 3)))
  }, striped = TRUE)

  output$plt_scen <- renderPlot({
    scen_data() %>%
      filter(time <= 730) %>%
      ggplot(aes(time / 30.4, gait, colour = scenario)) +
      geom_line(linewidth = 0.7, show.legend = FALSE) +
      facet_wrap(~scenario, ncol = 5) +
      labs(x = "개월", y = "보행 속도 [m/s]",
           title = "21개 시나리오의 보행 궤적") + THEME
  })

  ## ===================== TAB 13 — safety ================================ ##
  output$plt_safety <- renderPlot({
    d <- sim()
    d %>% mutate(SDH_pct = 100 * SDH_inc) %>%
      select(time, Vsdh, SDH_pct, Headx, HCO3, Kser, Occl) %>%
      pivot_longer(-time) %>%
      ggplot(aes(time / 30.4, value)) +
      geom_line(colour = PAL[2], linewidth = 0.9) +
      facet_wrap(~name, scales = "free_y", ncol = 3) +
      labs(x = "개월", y = NULL,
           title = "안전성: 과배액 · 저압 두통 · 대사성 산증 · 폐색") + THEME
  })

  output$tbl_safety <- renderTable({
    d <- sim()
    tibble(항목 = c("경막하 수액낭 [mL]", "증상성 경막하 누적 발생률 [%]",
                    "저압 두통 지수", "혈청 HCO3 [mmol/L]",
                    "혈청 K [mmol/L]", "션트 폐색 분율",
                    "누적 배액량 [mL]"),
           종료값 = c(fmt(last(d$Vsdh), 1), fmt(100 * last(d$SDH_inc), 1),
                     fmt(last(d$Headx)), fmt(last(d$HCO3), 1),
                     fmt(last(d$Kser)), fmt(last(d$Occl), 3),
                     fmt(last(d$Vdr), 0)))
  }, striped = TRUE)
}

shinyApp(ui, server)
