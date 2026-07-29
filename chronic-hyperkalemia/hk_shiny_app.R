## =============================================================================
##  CHRONIC HYPERKALAEMIA -- QSP DASHBOARD (Shiny)
##  Companion to hk_mrgsolve_model.R
## =============================================================================
##
##  RUN
##      install.packages(c("shiny","mrgsolve","dplyr","tidyr","ggplot2","DT","gridExtra"))
##      shiny::runApp("hk_shiny_app.R")
##
##  DESIGN NOTE -- WHY THIS APP IS LAID OUT THE WAY IT IS
##  ------------------------------------------------------
##  Almost every hyperkalaemia tool ever built plots serum potassium and stops.
##  Serum potassium is the LEAST informative quantity in this disease, because
##  it is a ratio of two things that move on different timescales and for
##  different reasons.  So every tab here is built to show the reader something
##  the potassium number itself hides:
##
##    Tab 1  Patient          what you are simulating
##    Tab 2  The two pools    K_total vs LAMrel -- the decomposition
##    Tab 3  Kidney           the kaliuresis reserve (FE_K), the thing that
##                            was moving for years while serum K sat still
##    Tab 4  Drug PK          exposures, MR occupancy, binder capture fraction
##    Tab 5  Acute rescue     what each emergency drug actually does, with a
##                            running "mmol removed" counter next to it
##    Tab 6  The dilemma      RAASi dose vs potassium, clinician in the loop
##    Tab 7  Scenarios        side-by-side long-term comparison
##    Tab 8  Population       who a binder rescues, and who it cannot
##    Tab 9  Validation       every trial the model was checked against,
##                            including the ones it disagrees with
##
## =============================================================================

library(shiny)
library(mrgsolve)
library(dplyr)
library(tidyr)
library(ggplot2)
library(DT)

mod <- mread_cache("hk", "hk_mrgsolve_model.R")

RUN_IN <- 250

thm <- theme_bw(base_size = 12) +
  theme(panel.grid.minor = element_blank(),
        strip.background = element_rect(fill = "grey92", colour = NA),
        legend.position = "bottom")

## Reference bands used on every potassium plot
k_bands <- function() list(
  annotate("rect", xmin = -Inf, xmax = Inf, ymin = 3.5, ymax = 5.0,
           fill = "#2f7f4f", alpha = 0.07),
  geom_hline(yintercept = 5.5, linetype = "22", colour = "#c07000"),
  geom_hline(yintercept = 6.0, linetype = "22", colour = "#b03030"),
  geom_hline(yintercept = 3.5, linetype = "22", colour = "#1f6bb8"))

sim_ss <- function(pars, tend = RUN_IN) {
  mod %>% param(pars) %>%
    mrgsim(end = tend, delta = 1, atol = 1e-8, rtol = 1e-8) %>%
    as.data.frame()
}

# =============================================================================
#  UI
# =============================================================================
ui <- fluidPage(
  titlePanel("만성 고칼륨혈증 QSP 대시보드 — Chronic Hyperkalaemia QSP Dashboard"),
  tags$p(style = "color:#555;margin-top:-8px",
    HTML("serum K is a RATIO, not a pool: &nbsp;<code>K_total = Ce&middot;V_ECF +
          Ci0&middot;LAMrel&middot;(Ce/Ce0)<sup>&alpha;</sup>&middot;V_ICF</code>
          &nbsp;&mdash;&nbsp; every tab below is built to show what that number hides")),

  sidebarLayout(
    sidebarPanel(
      width = 3,
      h4("환자 (Patient)"),
      sliderInput("gfr", "eGFR (mL/min/1.73)", 8, 100, 25, step = 1),
      sliderInput("intake", "식이 칼륨 Dietary K (mmol/day)", 20, 180, 80, step = 5),
      sliderInput("hco3", "혈청 HCO3 (mmol/L) — 0 = model-derived",
                  0, 28, 0, step = 1),
      sliderInput("bw", "체중 Body weight (kg)", 40, 130, 70, step = 1),

      hr(), h4("만성 약물 (Chronic therapy)"),
      sliderInput("ace", "ACE inhibitor (mg/day, lisinopril-equivalent)",
                  0, 40, 20, step = 5),
      sliderInput("mra", "MRA (mg/day, spironolactone-equivalent)",
                  0, 50, 25, step = 5),
      radioButtons("binder", "칼륨 결합제 K binder",
                   c("none", "patiromer", "SZC"), inline = TRUE),
      conditionalPanel("input.binder == 'patiromer'",
        sliderInput("pat", "patiromer (g/day)", 0, 33.6, 16.8, step = 4.2)),
      conditionalPanel("input.binder == 'SZC'",
        sliderInput("szc", "SZC (g/day)", 0, 30, 10, step = 5)),
      sliderInput("fur", "furosemide (mg/day)", 0, 160, 0, step = 20),
      sliderInput("bic", "oral alkali (mmol/day)", 0, 140, 0, step = 10),
      checkboxInput("sglt2i", "SGLT2 inhibitor", FALSE),

      hr(), h4("장기 시뮬레이션 (Long run)"),
      checkboxInput("titrate", "임상의 감량 루프 (clinician titration loop)", TRUE),
      checkboxInput("progress", "eGFR 진행 (let eGFR decline)", TRUE),
      sliderInput("years", "years", 1, 10, 5, step = 1),

      hr(),
      helpText(HTML("<b>읽는 법.</b> 혈청 K 열이 아니라 <b>RAASi 용량 열</b>과
        <b>제거된 mmol</b> 열을 보십시오. 거의 모든 처방은 결국 K 를 5.5 아래로
        내립니다 — 차이는 <i>그 대가로 무엇을 잃었는가</i> 입니다."))
    ),

    mainPanel(
      width = 9,
      tabsetPanel(
        id = "tabs",

        tabPanel("1. 환자 Patient",
          br(), fluidRow(
            column(4, wellPanel(h4("정상 상태 (steady state)"), tableOutput("tbl_ss"))),
            column(8, plotOutput("p_runin", height = "380px"))),
          helpText("모델을 250일 동안 돌려 정상 상태에 도달시킨 뒤의 값입니다.
                    ALDO·RASDN 이 안정되는 데 며칠이 걸린다는 점에 주목하십시오 —
                    MRA 로 인한 고칼륨혈증이 즉시가 아니라 며칠 뒤에 나타나는
                    이유가 바로 이것입니다.")),

        tabPanel("2. 두 저장고 Pools",
          br(), fluidRow(
            column(6, plotOutput("p_pool", height = "330px")),
            column(6, plotOutput("p_buffer", height = "330px"))),
          hr(),
          fluidRow(column(12, wellPanel(
            h4("같은 숫자, 반대의 저장고 (same number, opposite pool)"),
            tableOutput("tbl_two")))),
          helpText(HTML("교환환율(mmol 당 mmol/L)은 <b>적합된 값이 아니라</b>
            &alpha; 로부터 유도됩니다. 만성 ~224 mmol/(mmol/L) 은 고전적
            결핍 노모그램(혈청 K 3.0 → 200-400 mmol 결핍)을 재현하고,
            급성 ~66 은 정맥 KCl 40 mmol 이 위험한 이유를 설명합니다."))),

        tabPanel("3. 신장 예비능 Kidney reserve",
          br(), plotOutput("p_fek", height = "420px"),
          hr(), fluidRow(
            column(6, plotOutput("p_thresh", height = "300px")),
            column(6, tableOutput("tbl_thresh"))),
          helpText("FE_K 곡선이 이 질환의 전부입니다. 혈청 K 는 eGFR 이
                    떨어져도 오랫동안 평평하고, 그동안 움직이고 있었던 것은
                    분획배설률입니다. 혈청 K 는 이미 소진된 예비능의
                    '후행 지표'입니다.")),

        tabPanel("4. 약물 PK/PD",
          br(), fluidRow(
            column(6, plotOutput("p_pk", height = "300px")),
            column(6, plotOutput("p_occ", height = "300px"))),
          hr(), plotOutput("p_binder_dr", height = "320px"),
          helpText("결합제의 phi(포획 분율)에는 구조적 상한이 있습니다:
                    어떤 용량으로도 phi_max x f_abs x 식이섭취량 이상을 제거할 수
                    없습니다. 이 상한 아래로 균형이 무너진 환자에게 투석은
                    선호가 아니라 질량보존의 요구입니다.")),

        tabPanel("5. 급성 처치 Acute",
          br(),
          fluidRow(column(3, numericInput("acuteK", "시작 혈청 K", 6.8, 5.5, 9, 0.1)),
                   column(3, numericInput("insU", "insulin (U)", 10, 0, 20, 1)),
                   column(3, numericInput("dex", "dextrose (g)", 25, 0, 50, 5)),
                   column(3, numericInput("salb", "salbutamol neb (mg)", 20, 0, 20, 5))),
          plotOutput("p_acute", height = "330px"),
          plotOutput("p_removed", height = "260px"),
          helpText(HTML("<b>위 두 그림을 반드시 함께 보십시오.</b> 위 그림에서
            내려간 선이 아래 그림에서 0 에 머물러 있다면, 그 약은 칼륨을
            <i>옮겼을 뿐 제거하지 않았습니다</i> — 반드시 되돌아옵니다."))),

        tabPanel("6. RAASi 딜레마",
          br(), plotOutput("p_dilemma", height = "440px"),
          hr(), tableOutput("tbl_dilemma"),
          helpText("칼슘·인슐린이 아니라 이 탭이 이 질환의 실제 임상 문제입니다.
                    감량 루프를 끄고 켜 보면, 고칼륨혈증의 해악 대부분이
                    부정맥이 아니라 '중단된 신보호 치료'를 통해 발생함을
                    볼 수 있습니다.")),

        tabPanel("7. 시나리오 비교 Scenarios",
          br(), checkboxGroupInput("scen", "비교할 처방",
            choices = c("no binder", "patiromer 16.8 g", "SZC 10 g",
                        "low-K diet 50", "furosemide 40", "SGLT2i",
                        "alkali 70", "no MRA"),
            selected = c("no binder", "patiromer 16.8 g", "low-K diet 50", "no MRA"),
            inline = TRUE),
          plotOutput("p_scen", height = "460px"),
          DTOutput("dt_scen")),

        tabPanel("8. 가상 인구 Population",
          br(), fluidRow(
            column(4, sliderInput("pop_n", "격자 해상도 (grid per axis)", 4, 12, 8)),
            column(8, helpText("eGFR x 식이 K x HCO3 격자. 결합제가 실제로
                                구조해내는 환자는 '역치를 결합제의 상한보다
                                적게 초과한' 사람들이며, 그 밴드는 계산 가능한
                                양입니다."))),
          plotOutput("p_pop", height = "440px"),
          tableOutput("tbl_pop")),

        tabPanel("9. 검증 Validation",
          br(), h4("모델이 맞춘 것과 틀린 것 (what the model got right, and wrong)"),
          DTOutput("dt_valid"),
          hr(),
          h4("모델이 다루지 못하는 것 (out of scope)"),
          tags$ul(
            tags$li("AKI on CKD — 실제 중증 고칼륨혈증의 최다 원인"),
            tags$li("횡문근융해·종양용해 — ICF 저장고로부터의 방출"),
            tags$li("가성 고칼륨혈증 (용혈, 백혈구/혈소판 증가)"),
            tags$li("투석 — 간헐적 초고청소율 제거항"),
            tags$li("운동 유발 일과성 상승"),
            tags$li("결과층(hazard)은 관찰연구 기반이며 인과적이지 않음 —
                     '숫자를 낮추는 것'의 이득을 과대평가할 가능성이 큼")))
      )
    )
  )
)

# =============================================================================
#  SERVER
# =============================================================================
server <- function(input, output, session) {

  base_pars <- reactive({
    p <- list(GFR_SET = input$gfr, INTAKE = input$intake, BW = input$bw,
              ACE_TGT = input$ace, MRA_TGT = input$mra,
              FUR_RATE = input$fur, BIC_RATE = input$bic,
              SGLT2I = as.numeric(input$sglt2i),
              HCO3_SET = input$hco3,
              PAT_RATE = if (input$binder == "patiromer") input$pat else 0,
              SZC_RATE = if (input$binder == "SZC") input$szc else 0)
    p
  })

  run_in <- reactive(sim_ss(base_pars()))
  ss     <- reactive(dplyr::filter(run_in(), time == max(time)))

  # ---- Tab 1 ---------------------------------------------------------------
  output$tbl_ss <- renderTable({
    s <- ss()
    data.frame(
      quantity = c("serum K (mmol/L)", "total body K (mmol)", "LAMrel",
                   "urine K (mmol/day)", "FE_K (%)", "colonic K (mmol/day)",
                   "aldosterone (ng/dL)", "MR occupancy", "HCO3 (mmol/L)",
                   "arterial pH", "binder capture phi", "QRS (ms)",
                   "composite hazard"),
      value = c(sprintf("%.2f", s$K), sprintf("%.0f", s$KTOT),
                sprintf("%.3f", s$LAMrel), sprintf("%.1f", s$UK),
                sprintf("%.1f", s$FEK), sprintf("%.1f", s$COLK_EX),
                sprintf("%.1f", s$ALDO), sprintf("%.2f", s$MROCC),
                sprintf("%.1f", s$HCO3ser), sprintf("%.3f", s$PHA),
                sprintf("%.3f", s$PHI), sprintf("%.0f", s$QRS),
                sprintf("%.2f", s$HZTOT)))
  }, colnames = FALSE)

  output$p_runin <- renderPlot({
    run_in() %>%
      select(time, K, ALDO, RASDN, UK, FEK) %>%
      pivot_longer(-time) %>%
      ggplot(aes(time, value)) +
      geom_line(linewidth = 0.9, colour = "#1f6bb8") +
      facet_wrap(~name, scales = "free_y") +
      labs(x = "days of run-in", y = NULL,
           title = "정상 상태까지의 접근 — ASDN 수송체 풀(RASDN)이 가장 느립니다") +
      thm
  })

  # ---- Tab 2 ---------------------------------------------------------------
  output$p_pool <- renderPlot({
    s <- ss()
    data.frame(pool = factor(c("ECF", "fast ICF", "slow ICF (muscle)"),
                             levels = c("ECF", "fast ICF", "slow ICF (muscle)")),
               mmol = c(s$KE, s$KIF, s$KIS)) %>%
      ggplot(aes(pool, mmol, fill = pool)) +
      geom_col(width = 0.6) +
      geom_text(aes(label = sprintf("%.0f\n(%.1f%%)", mmol, 100*mmol/s$KTOT)),
                vjust = -0.2, size = 3.6) +
      scale_fill_manual(values = c("#1f6bb8", "#7fb0dd", "#c9dcef"), guide = "none") +
      expand_limits(y = max(s$KIS)*1.15) +
      labs(title = "체내 칼륨 분포", subtitle = "혈청이 보고 있는 것은 왼쪽 막대뿐입니다",
           y = "mmol", x = NULL) + thm
  })

  output$p_buffer <- renderPlot({
    s <- ss()
    kk <- seq(2, 8, by = 0.05)
    a <- 0.25; ci0 <- 140; ce0 <- 4.2
    v_e <- 0.20*input$bw; v_i <- 0.36*input$bw; v_f <- 0.25*v_i
    ktot <- kk*v_e + ci0*(kk/ce0)^a*v_i
    kref <- ce0*v_e + ci0*v_i
    data.frame(K = kk, delta = ktot - kref) %>%
      ggplot(aes(K, delta)) +
      geom_hline(yintercept = 0, colour = "grey60") +
      geom_line(linewidth = 1.1, colour = "#a02a72") +
      geom_vline(xintercept = s$K, linetype = "22") +
      annotate("point", x = 3.0,
               y = 3.0*v_e + ci0*(3.0/ce0)^a*v_i - kref, size = 3, colour = "#1f6bb8") +
      annotate("text", x = 3.0, y = 3.0*v_e + ci0*(3.0/ce0)^a*v_i - kref,
               label = "  K 3.0: 고전 노모그램 200-400 mmol 결핍",
               hjust = 0, size = 3.4, colour = "#1f6bb8") +
      labs(title = "혈청 K 대 전신 칼륨 과부족",
           subtitle = "이 곡선의 기울기가 '한 단위의 검사 수치'가 몇 mmol 인지입니다",
           x = "serum K (mmol/L)", y = "whole-body K, deviation from normal (mmol)") +
      thm
  })

  output$tbl_two <- renderTable({
    s <- ss()
    data.frame(
      patient = c("A — 질량수지형 (CKD, 결합제가 유일한 답)",
                  "B — 분배형 (DKA, 인슐린이 답이고 결합제는 해롭다)"),
      `serum K` = c(sprintf("%.2f", s$K), "4.7"),
      LAMrel = c(sprintf("%.3f", s$LAMrel), "0.86"),
      `whole-body K` = c(sprintf("%+.0f mmol", s$KTOT - (4.2*0.20*input$bw + 140*0.36*input$bw)),
                         "-400 mmol"),
      `치료 후` = c("결합제로 서서히 하강, 되돌아오지 않음",
                     "인슐린으로 K 2.7 까지 하강 — 그것이 심정지"),
      check.names = FALSE)
  })

  # ---- Tab 3 ---------------------------------------------------------------
  reserve <- reactive({
    gfrs <- c(90, 75, 60, 45, 35, 25, 20, 15, 12, 10)
    arms <- list("none" = list(ACE_TGT = 0, MRA_TGT = 0),
                 "ACEi" = list(ACE_TGT = 20, MRA_TGT = 0),
                 "ACEi+MRA" = list(ACE_TGT = 20, MRA_TGT = 25),
                 "ACEi+MRA+binder" = list(ACE_TGT = 20, MRA_TGT = 25, PAT_RATE = 16.8))
    do.call(rbind, lapply(names(arms), function(a)
      do.call(rbind, lapply(gfrs, function(g) {
        r <- dplyr::filter(sim_ss(c(base_pars()[c("INTAKE", "BW")], arms[[a]],
                                    list(GFR_SET = g))), time == max(time))
        data.frame(arm = a, eGFR = g, K = r$K, FEK = r$FEK, UK = r$UK,
                   SCAP = r$SCAP, COL = r$COLK_EX)
      }))))
  })

  output$p_fek <- renderPlot({
    d <- reserve()
    p1 <- ggplot(d, aes(eGFR, K, colour = arm)) + k_bands() +
      geom_line(linewidth = 1.1) + geom_point(size = 1.6) +
      scale_x_reverse() +
      labs(title = "혈청 K 는 예비능이 소진될 때까지 평평하다",
           y = "serum K (mmol/L)", colour = NULL) + thm
    p2 <- ggplot(dplyr::filter(d, arm == "none"), aes(eGFR, FEK)) +
      geom_line(linewidth = 1.1, colour = "#2f7f8c") + geom_point(size = 1.6) +
      scale_x_reverse() +
      labs(title = "…그동안 움직이고 있던 것은 이쪽이다 (칼륨 분획배설률)",
           y = "FE_K (%)") + thm
    gridExtra::grid.arrange(p1, p2, ncol = 1)
  })

  thresholds <- reactive({
    arms <- list("none" = list(ACE_TGT = 0, MRA_TGT = 0),
                 "ACEi" = list(ACE_TGT = 20, MRA_TGT = 0),
                 "ACEi+MRA" = list(ACE_TGT = 20, MRA_TGT = 25),
                 "ACEi+MRA+binder" = list(ACE_TGT = 20, MRA_TGT = 25, PAT_RATE = 16.8))
    data.frame(arm = names(arms), threshold = vapply(arms, function(a) {
      lo <- 5; hi <- 120
      for (i in 1:30) {
        mid <- (lo + hi)/2
        k <- dplyr::filter(sim_ss(c(base_pars()[c("INTAKE","BW")], a,
                                    list(GFR_SET = mid))), time == max(time))$K
        if (k > 5.5) lo <- mid else hi <- mid
      }
      (lo + hi)/2
    }, numeric(1)))
  })

  output$p_thresh <- renderPlot({
    ggplot(thresholds(), aes(reorder(arm, threshold), threshold, fill = arm)) +
      geom_col(width = 0.6) + coord_flip() +
      geom_text(aes(label = sprintf("%.1f", threshold)), hjust = -0.15, size = 4) +
      scale_fill_brewer(palette = "Blues", guide = "none") +
      expand_limits(y = max(thresholds()$threshold)*1.2) +
      labs(title = "K 5.5 를 넘는 eGFR 역치",
           subtitle = "MRA 는 K 를 올리는 것이 아니라 이 역치를 올린다",
           x = NULL, y = "eGFR (mL/min/1.73)") + thm
  })
  output$tbl_thresh <- renderTable(thresholds(), digits = 1)

  # ---- Tab 4 ---------------------------------------------------------------
  output$p_pk <- renderPlot({
    run_in() %>% select(time, CACE, CMRA, PATC, SZCP) %>%
      pivot_longer(-time) %>% dplyr::filter(time <= 20) %>%
      ggplot(aes(time, value)) + geom_line(linewidth = 0.9, colour = "#7040a0") +
      facet_wrap(~name, scales = "free_y") +
      labs(title = "약물 노출 (첫 20일)", x = "days", y = NULL) + thm
  })

  output$p_occ <- renderPlot({
    run_in() %>% select(time, MROCC, RASDN, PHI) %>%
      pivot_longer(-time) %>% dplyr::filter(time <= 30) %>%
      ggplot(aes(time, value, colour = name)) + geom_line(linewidth = 1.0) +
      labs(title = "MR 점유율 · ASDN 수송체 풀 · 결합제 포획분율",
           subtitle = "RASDN 의 t1/2 1.5일이 MRA 고칼륨혈증의 지연을 만든다",
           x = "days", y = NULL, colour = NULL) + thm
  })

  output$p_binder_dr <- renderPlot({
    doses <- seq(0, 34, by = 2)
    d <- do.call(rbind, lapply(c("patiromer", "SZC"), function(drug)
      do.call(rbind, lapply(doses, function(dz) {
        pp <- base_pars(); pp$PAT_RATE <- 0; pp$SZC_RATE <- 0
        if (drug == "patiromer") pp$PAT_RATE <- dz else pp$SZC_RATE <- dz
        r <- dplyr::filter(sim_ss(pp), time == max(time))
        data.frame(drug = drug, dose = dz, K = r$K, phi = r$PHI)
      }))))
    ggplot(d, aes(dose, K, colour = drug)) + k_bands() +
      geom_line(linewidth = 1.1) +
      labs(title = "결합제 용량-반응 (정상 상태)",
           subtitle = sprintf("식이 %d mmol/day 에서의 구조적 상한에 주목", input$intake),
           x = "dose (g/day)", y = "serum K (mmol/L)", colour = NULL) + thm
  })

  # ---- Tab 5 ---------------------------------------------------------------
  acute <- reactive({
    p <- base_pars()
    s <- ss()
    dK <- (input$acuteK - s$K)*s$BUF_CHR
    v_e <- 0.20*input$bw
    i0 <- list(KE = s$KE + dK*v_e/s$BUF_CHR,
               KIF = s$KIF + dK*(1 - v_e/s$BUF_CHR)*0.25,
               KIS = s$KIS + dK*(1 - v_e/s$BUF_CHR)*0.75)
    run <- function(lab, ev_ = NULL, extra = list()) {
      m <- mod %>% param(c(p, extra)) %>% init(i0)
      o <- if (is.null(ev_)) mrgsim(m, end = 2, delta = 1/288)
           else mrgsim(m, ev_, end = 2, delta = 1/288)
      as.data.frame(o) %>% mutate(arm = lab)
    }
    rbind(
      run("무처치 none"),
      run("calcium 1 g", ev(amt = 0.25, cmt = "CAE")),
      run(sprintf("insulin %g U + D50 %g g", input$insU, input$dex),
          ev(amt = input$insU*1000/12, cmt = "INS") +
          ev(amt = input$dex/180.15*1000/16, cmt = "GLU")),
      run(sprintf("salbutamol %g mg", input$salb),
          ev(amt = input$salb*0.20/150, cmt = "B2C")),
      run("SZC 10 g TID", NULL, list(SZC_RATE = 30)))
  })

  output$p_acute <- renderPlot({
    ggplot(acute(), aes(time*24, K, colour = arm)) + k_bands() +
      geom_line(linewidth = 1.0) +
      labs(title = "혈청 칼륨 — 위쪽만 보면 인슐린이 최고의 약처럼 보인다",
           x = "hours", y = "serum K (mmol/L)", colour = NULL) + thm
  })

  output$p_removed <- renderPlot({
    ggplot(acute(), aes(time*24, KREM, colour = arm)) +
      geom_line(linewidth = 1.0) +
      labs(title = "…실제로 몸에서 제거된 칼륨 (mmol)",
           subtitle = "칼슘·인슐린·살부타몰은 정확히 0 에 머문다",
           x = "hours", y = "cumulative K removed (mmol)", colour = NULL) + thm
  })

  # ---- Tab 6/7 -------------------------------------------------------------
  SCEN <- list("no binder"        = list(),
               "patiromer 16.8 g" = list(PAT_RATE = 16.8),
               "SZC 10 g"         = list(SZC_RATE = 10),
               "low-K diet 50"    = list(INTAKE = 50),
               "furosemide 40"    = list(FUR_RATE = 40),
               "SGLT2i"           = list(SGLT2I = 1),
               "alkali 70"        = list(BIC_RATE = 70),
               "no MRA"           = list(MRA_TGT = 0))

  long_run <- reactive({
    sel <- if (input$tabs == "6. RAASi 딜레마") names(SCEN) else input$scen
    do.call(rbind, lapply(sel, function(a) {
      pp <- base_pars()
      pp$PAT_RATE <- 0; pp$SZC_RATE <- 0
      pp <- modifyList(pp, SCEN[[a]])
      pp$TITRATE <- as.numeric(input$titrate)
      pp$PROGRESS <- as.numeric(input$progress)
      o <- as.data.frame(mrgsim(param(mod, pp), end = input$years*365, delta = 7,
                                atol = 1e-8, rtol = 1e-8))
      o$arm <- a; o
    }))
  })

  output$p_dilemma <- renderPlot({
    long_run() %>% select(time, arm, K, RDo, GFR, HZTOT) %>%
      pivot_longer(c(K, RDo, GFR, HZTOT)) %>%
      mutate(name = recode(name, K = "serum K (mmol/L)",
                           RDo = "prescribed RAASi dose (fraction)",
                           GFR = "eGFR (mL/min/1.73)",
                           HZTOT = "composite hazard ratio")) %>%
      ggplot(aes(time/365, value, colour = arm)) +
      geom_line(linewidth = 0.95) +
      facet_wrap(~name, scales = "free_y") +
      labs(title = "RAASi–칼륨 딜레마",
           subtitle = "K 패널이 아니라 RAASi 용량 패널을 보십시오",
           x = "years", y = NULL, colour = NULL) + thm
  })

  output$tbl_dilemma <- renderTable({
    long_run() %>% group_by(arm) %>% dplyr::filter(time == max(time)) %>%
      transmute(arm, `K (5y)` = K, `RAASi dose` = RDo, `eGFR (5y)` = GFR,
                `days K>5.5` = T55, `mean hazard` = CHAZ/time,
                `K removed (mmol)` = KREM) %>% ungroup()
  }, digits = 2)

  output$p_scen <- renderPlot({
    long_run() %>%
      ggplot(aes(time/365, K, colour = arm)) + k_bands() +
      geom_line(linewidth = 1.0) +
      labs(x = "years", y = "serum K (mmol/L)", colour = NULL,
           title = "시나리오 비교") + thm
  })
  output$dt_scen <- renderDT({
    long_run() %>% group_by(arm) %>% dplyr::filter(time == max(time)) %>%
      transmute(arm, K = round(K, 2), RAASi = round(RDo, 2),
                eGFR = round(GFR, 1), days_K_over_5.5 = round(T55),
                mean_hazard = round(CHAZ/time, 3),
                K_removed_mmol = round(KREM)) %>% ungroup()
  }, options = list(dom = "t", pageLength = 10))

  # ---- Tab 8 ---------------------------------------------------------------
  output$p_pop <- renderPlot({
    n <- input$pop_n
    grid <- expand.grid(GFR_SET = seq(15, 60, length.out = n),
                        INTAKE = seq(40, 140, length.out = n),
                        HCO3_SET = seq(16, 26, length.out = 4))
    res <- do.call(rbind, lapply(seq_len(nrow(grid)), function(i) {
      pp <- as.list(grid[i, ])
      a <- dplyr::filter(sim_ss(c(pp, list(ACE_TGT = 20, MRA_TGT = 25,
                                           BW = input$bw))), time == max(time))$K
      b <- dplyr::filter(sim_ss(c(pp, list(ACE_TGT = 20, MRA_TGT = 25,
                                           BW = input$bw, PAT_RATE = 16.8))),
                         time == max(time))$K
      data.frame(grid[i, ], K_no = a, K_bind = b)
    }))
    res %>% mutate(status = case_when(
        K_no <= 5.5 ~ "결합제 불필요 (already < 5.5)",
        K_bind <= 5.5 ~ "결합제가 구조 (rescued)",
        TRUE ~ "결합제로도 부족 (still > 5.5)")) %>%
      ggplot(aes(GFR_SET, INTAKE, fill = status)) +
      geom_tile(colour = "white") +
      facet_wrap(~sprintf("HCO3 %.0f", HCO3_SET)) +
      scale_fill_manual(values = c("결합제 불필요 (already < 5.5)" = "#cfe6d2",
                                   "결합제가 구조 (rescued)" = "#1f6bb8",
                                   "결합제로도 부족 (still > 5.5)" = "#b03030")) +
      labs(title = "결합제가 실제로 구조하는 환자는 누구인가",
           x = "eGFR (mL/min/1.73)", y = "dietary K (mmol/day)", fill = NULL) + thm
  })

  output$tbl_pop <- renderTable({
    data.frame(note = c(
      "결합제의 이득은 균일하지 않습니다.",
      "역치를 '결합제의 상한보다 적게' 초과한 환자에게 집중됩니다.",
      "그 상한은 phi_max x f_abs x 식이섭취량 이며 계산 가능합니다.",
      "이 상한을 넘어선 환자에게 투석은 선호가 아니라 질량보존의 요구입니다."))
  }, colnames = FALSE)

  # ---- Tab 9 ---------------------------------------------------------------
  output$dt_valid <- renderDT({
    data.frame(
      target = c("eGFR 60 정상상태 K", "eGFR 45 정상상태 K", "eGFR 30 정상상태 K",
                 "혈청 K 3.0 의 전신 결핍", "OPAL-HK 파티로머 16.8 g 4주",
                 "HARMONIZE SZC 5/10/15 g 28일", "HARMONIZE SZC 10 g TID 48시간",
                 "RALES 스피로노락톤 25 mg", "FIDELIO 피네레논 20 mg",
                 "인슐린 10 U + 포도당", "살부타몰 20 mg 분무",
                 "ACE 억제제 단독", "중탄산염 HCO3 18→24",
                 "인슐린/포도당 후 저혈당 발생률"),
      observed = c("4.40", "4.55", "4.75", "-200 ~ -400 mmol", "-1.01 mmol/L",
                   "4.8 / 4.5 / 4.4", "약 -1.1 mmol/L", "+0.30 mmol/L",
                   "+0.23 mmol/L", "-0.6 ~ -1.0 (30-60분)", "-0.6 ~ -1.0",
                   "+0.1 ~ +0.4", "약 -0.2 ~ -0.3", "15-20%"),
      model = c("4.47", "4.61", "4.80", "-302 mmol", "-0.95 mmol/L",
                "4.89 / 4.63 / 4.49", "-0.95 mmol/L", "+0.30 mmol/L",
                "+0.10 mmol/L", "-0.85 (약 50분)", "-0.74 (2-4시간)",
                "+0.21 ~ +0.23", "-0.28 (전량 신장 경로)", "0% — 재현 못함"),
      verdict = c("예측(held-out)", "예측(held-out)", "예측(held-out)",
                  "예측 — α 로부터 유도, 적합 아님",
                  "일치", "일치 (일관되게 +0.1 높음)", "일치", "적합 앵커",
                  "불일치 — 가정한 MR 부하가 낮음 (역산 시 72% 필요)",
                  "일치", "일치", "일치", "적합 앵커",
                  "불일치 — 구조적 한계, 모델 신뢰 불가"),
      check.names = FALSE)
  }, options = list(dom = "t", pageLength = 20))
}

shinyApp(ui, server)
