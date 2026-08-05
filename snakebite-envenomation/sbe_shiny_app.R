## =============================================================================
##  sbe_shiny_app.R
##  Snakebite envenoming + antivenom — interactive QSP dashboard
##  독사 교상 + 항독소 — 인터랙티브 QSP 대시보드
##
##  15 tabs.  The organising idea of the whole app is the model's thesis:
##
##        ANTIVENOM BINDS.  IT DOES NOT UNDO.
##
##  so the app is laid out as two clocks.  Tabs 2-3 are the VENOM clock (hours),
##  where antivenom acts.  Tabs 4-9 are the SUBSTRATE clock (days), where every
##  clinical endpoint is read.  Tab 12 puts the four claims side by side so the
##  user can try to break them.
##
##  Run with:  shiny::runApp("sbe_shiny_app.R")
##  Requires:  shiny, mrgsolve, dplyr, tidyr, ggplot2, DT
##  Depends on: sbe_mrgsolve_model.R  (in the same directory)
## =============================================================================

library(shiny)
library(mrgsolve)
library(dplyr)
library(tidyr)
library(ggplot2)
library(DT)

source("sbe_mrgsolve_model.R")

THEME <- theme_bw(base_size = 12) +
  theme(legend.position = "bottom", legend.title = element_blank(),
        plot.title = element_text(face = "bold", size = 12),
        plot.subtitle = element_text(size = 10, colour = "grey30"),
        strip.background = element_rect(fill = "grey92"))

# a shaded band for a clinically meaningful threshold
band <- function(ymin, ymax, fill = "#D94A4A", alpha = 0.10) {
  annotate("rect", xmin = -Inf, xmax = Inf, ymin = ymin, ymax = ymax,
           fill = fill, alpha = alpha)
}

# =============================================================================
#  UI
# =============================================================================
ui <- fluidPage(
  titlePanel(
    div(
      h3("독사 교상(뱀독 중독) + 항독소 QSP 모델 · Snakebite Envenoming & Antivenom"),
      h5(em("50-ODE model · 4 toxin classes × 3 antivenom fragment formats × 8 snake archetypes")),
      h5(strong("THESIS: antivenom binds; it does not undo."),
         " A venom clock in hours drives a substrate clock in days; every endpoint is read off the slow one.")
    )
  ),

  sidebarLayout(
    sidebarPanel(
      width = 3,

      h4("1 · 교상 (the bite)"),
      selectInput("snake", "종 (species)", choices = names(SNAKES),
                  selected = "Daboia russelii (Sri Lanka)"),
      sliderInput("venom_mult", "주입 독량 배수 (venom dose multiplier)",
                  min = 0.1, max = 3.5, value = 1.0, step = 0.1),
      helpText("A bite is not a calibrated injection. Between-subject CV is ~55%; ",
               "20-50% of viperid bites inject nothing at all."),

      hr(),
      h4("2 · 항독소 (antivenom)"),
      selectInput("product", "제품 형식 (fragment format)",
                  choices = c("none", names(AV)), selected = "IgG"),
      sliderInput("vials", "바이알 수 (vials)", min = 0, max = 40, value = 10, step = 1),
      sliderInput("t_av", "투여 시각 (h post-bite)", min = 0.5, max = 24,
                  value = 4, step = 0.5),
      sliderInput("av_infuse", "주입 시간 (infusion duration, h)",
                  min = 0.25, max = 4, value = 1, step = 0.25),
      helpText("Acute reactions are driven by the infusion RATE, not the total dose."),
      checkboxInput("maint", "유지 용량 (maintenance: 2 vials q6h × 3)", FALSE),

      hr(),
      h4("3 · 보조 치료 (adjuncts)"),
      checkboxInput("neo", "네오스티그민 0.5 mg IM q4h × 8", FALSE),
      checkboxInput("vares", "바레스플라딥 500 mg PO q12h × 8 (from 0.5 h)", FALSE),
      checkboxInput("txa", "트라넥사믹산 1 g IV q8h × 3", FALSE),
      checkboxInput("cryo", "동결침전제제 (cryoprecipitate, +0.9 g/L × 2)", FALSE),
      checkboxInput("fluid", "정질용액 0.30 L/h, 1-8 h", TRUE),

      hr(),
      h4("4 · 지지 치료 가용성 (what the hospital has)"),
      checkboxInput("icu", "인공환기 가능 (mechanical ventilation)", TRUE),
      checkboxInput("dial", "신대체요법 가능 (dialysis)", TRUE),
      helpText(strong("These two switches change mortality more than any antivenom "),
               strong("choice in the panel above."), " That is not a modelling artefact."),

      hr(),
      sliderInput("tmax", "관찰 기간 (h)", min = 72, max = 336, value = 336, step = 24),
      actionButton("run", "시뮬레이션 실행 (simulate)", class = "btn-primary"),
      br(), br(),
      wellPanel(
        h5(strong("화학량 (stoichiometry)")),
        verbatimTextOutput("stoich_box")
      )
    ),

    mainPanel(
      width = 9,
      tabsetPanel(
        id = "tabs", type = "pills",

        # ---- 1 ------------------------------------------------------------
        tabPanel(
          "1 · 프로파일",
          h4("교상 프로파일 및 즉시 요약 (bite profile and immediate summary)"),
          fluidRow(
            column(6, tableOutput("profile")),
            column(6, tableOutput("keynum"))
          ),
          hr(),
          h4("네 가지 종말점을 한 화면에 (four endpoints, one screen)"),
          plotOutput("overview", height = "480px"),
          helpText("Fibrinogen and the neuromuscular safety factor live on the ",
                   "substrate clock; free venom and antivenom live on the venom clock. ",
                   "The gap between them is the model.")
        ),

        # ---- 2 ------------------------------------------------------------
        tabPanel(
          "2 · 독액 PK",
          h4("독액 약동학 — 두 개의 저장고와 하나의 전신 구획"),
          plotOutput("pk_venom", height = "420px"),
          hr(),
          h4("혈장 유리 독소, 계열별 (free plasma toxin by class)"),
          plotOutput("pk_class", height = "340px"),
          hr(),
          h4("항원 검사가 보는 것 vs 환자가 느끼는 것"),
          plotOutput("pk_antigen", height = "300px"),
          helpText(strong("A venom antigen EIA measures free AND antivenom-bound toxin."),
                   " Antigen can be high while the patient is completely protected — ",
                   "which is why 'venom still detectable' is not an indication to re-dose.")
        ),

        # ---- 3 ------------------------------------------------------------
        tabPanel(
          "3 · 항독소 PK · 재발 부등식",
          h4("항독소 약동학"),
          plotOutput("pk_av", height = "360px"),
          hr(),
          h4("재발 부등식 (the recurrence inequality)"),
          tableOutput("halflives"),
          helpText("Recurrence is a ratio of two half-lives. If the fragment's ",
                   "terminal half-life is shorter than the venom INPUT half-life, ",
                   "free venom must reappear — and no larger front-loaded dose ",
                   "prevents it, only re-dosing or a longer-lived fragment."),
          hr(),
          h4("바이알 수는 화학량을 결정하지 않는다 (vials are not a margin)"),
          plotOutput("stoich_plot", height = "320px")
        ),

        # ---- 4 ------------------------------------------------------------
        tabPanel(
          "4 · 응고",
          h4("피브리노겐 — 두 개의 시계가 만나는 곳"),
          plotOutput("coag_fg", height = "360px"),
          helpText(strong("The falling limb is the venom's. The rising limb is the liver's."),
                   " The hepatic synthesis term in this model contains no antivenom ",
                   "variable at all: dose and timing move the DEPTH of the nadir, ",
                   "not the slope out of it."),
          hr(),
          h4("혈소판 · 제X인자 · D-dimer"),
          plotOutput("coag_other", height = "380px"),
          hr(),
          h4("20분 전혈응고검사 (20WBCT) — 유리관 하나로 하는 게이트"),
          plotOutput("coag_wbct", height = "220px")
        ),

        # ---- 5 ------------------------------------------------------------
        tabPanel(
          "5 · 신경근",
          h4("두 가지 마비, 같은 임상 소견, 정반대의 약리학"),
          plotOutput("nmj_sf", height = "380px"),
          helpText("SF = SF₀ · TERM_eff · R_eff — a PRODUCT. Zero in either factor ",
                   "is zero overall, and the two factors answer opposite questions."),
          hr(),
          fluidRow(
            column(6, h5("후시냅스 점유율 (postsynaptic occupancy, BR)"),
                   plotOutput("nmj_br", height = "300px"),
                   helpText("Reversible: decays on its own 35 h clock once free ",
                            "toxin is gone, and can be out-competed by raising ",
                            "acetylcholine.")),
            column(6, h5("전시냅스 종말 완전성 (presynaptic terminal, TERM)"),
                   plotOutput("nmj_term", height = "300px"),
                   helpText("Irreversible: the only recovery term is regeneration, ",
                            "t½ 5 days. Neither antivenom nor neostigmine appears ",
                            "anywhere in this equation."))
          ),
          hr(),
          h4("단일호흡계수 및 인공환기 필요 시간"),
          plotOutput("nmj_sbc", height = "280px")
        ),

        # ---- 6 ------------------------------------------------------------
        tabPanel(
          "6 · 국소 손상",
          h4("항독소가 도달할 수 없는 구획"),
          plotOutput("local_nec", height = "380px"),
          helpText("k_b,loc / k_b0 = 0.038 — a 27-fold penalty at the bite site, ",
                   "and a further ×4 in the deep depot. Intravenous antibody ",
                   "controls the plasma compartment and leaves the limb."),
          hr(),
          h4("부종 · 구획압 · CK"),
          plotOutput("local_cp", height = "380px"),
          helpText("Fasciotomy threshold 30 mmHg (dashed). CK release is driven by ",
                   "the necrosis STOCK, not its flux, which is why it peaks at ",
                   "12-48 h rather than at 4 h.")
        ),

        # ---- 7 ------------------------------------------------------------
        tabPanel(
          "7 · 신장",
          h4("신장은 수치가 아니라 적분이다 (the kidney is an integral)"),
          plotOutput("renal", height = "420px"),
          helpText(strong("Creatinine peaks on day 3-5, long after venom is undetectable."),
                   " Nothing observable at presentation predicts this; the ",
                   "time-integral of fibrin flux and pigment casts does."),
          hr(),
          h4("영구 손실의 천장 (the ceiling recovery cannot exceed)"),
          plotOutput("renal_scar", height = "300px")
        ),

        # ---- 8 ------------------------------------------------------------
        tabPanel(
          "8 · 순환",
          h4("혈액보다 혈장이 먼저 빠져나간다"),
          plotOutput("haemo", height = "420px"),
          helpText("The Russell's viper paradox: haematocrit RISES while the ",
                   "patient bleeds, because the glycocalyx is shed and plasma ",
                   "leaves the circulation before red cells do."),
          hr(),
          h4("염증 (inflammation)"),
          plotOutput("infl", height = "300px")
        ),

        # ---- 9 ------------------------------------------------------------
        tabPanel(
          "9 · 임상 종말점",
          h4("누적 위험 (cumulative hazards)"),
          plotOutput("endpoints", height = "380px"),
          hr(),
          h4("종말점 요약 (endpoint summary)"),
          tableOutput("endpoint_tbl"),
          hr(),
          h4("무엇이 실제로 환자를 살리는가 (what actually saves the patient)"),
          plotOutput("support_compare", height = "320px"),
          helpText("Same venom, same antivenom, same timing — the only difference ",
                   "is whether the hospital has a ventilator and a dialysis machine.")
        ),

        # ---- 10 -----------------------------------------------------------
        tabPanel(
          "10 · 항독소 이상반응",
          h4("장부의 비용 쪽 (the cost side of the ledger)"),
          plotOutput("ae", height = "400px"),
          helpText("The anaphylactoid term is driven by the INFUSION RATE of ",
                   "foreign protein — so slowing the infusion is a real ",
                   "intervention, and giving more vials faster is a real risk. ",
                   "Serum sickness peaks on day 5-12, i.e. after discharge, ",
                   "which is why every trial under-reports it."),
          hr(),
          h4("주입 속도 민감도 (infusion-rate sensitivity)"),
          plotOutput("ae_rate", height = "320px")
        ),

        # ---- 11 -----------------------------------------------------------
        tabPanel(
          "11 · 시나리오 비교",
          h4("27개 사전 정의 시나리오 (27 predefined scenarios)"),
          actionButton("run_all", "전체 시나리오 실행 (run all 27)", class = "btn-warning"),
          helpText("This runs 27 full 14-day simulations; allow a few seconds."),
          br(),
          DT::dataTableOutput("scen_tbl"),
          hr(),
          h4("선택 시나리오 겹쳐 보기 (overlay selected scenarios)"),
          selectInput("ov_var", "변수",
                      choices = c("Fibrinogen" = "FG", "Free venom" = "CVFREE",
                                  "Antivenom" = "CA_OUT", "Safety factor" = "SFOUT",
                                  "Creatinine" = "SCR", "Platelets" = "PLT",
                                  "Necrosis" = "NEC", "CK" = "CK",
                                  "Serum sickness" = "SS"),
                      selected = "FG"),
          plotOutput("scen_overlay", height = "420px")
        ),

        # ---- 12 -----------------------------------------------------------
        tabPanel(
          "12 · 네 가지 주장",
          h4("이 탭은 모델을 깨뜨려 보기 위한 것입니다 (try to break the claims)"),
          br(),
          h4("주장 1 · 항독소는 나디르의 깊이를 바꾸고, 회복 기울기는 바꾸지 않는다"),
          plotOutput("claim1", height = "360px"),
          tableOutput("claim1_tbl"),
          hr(),
          h4("주장 2 · 재발은 절편의 성질이며 바이알 수의 성질이 아니다"),
          plotOutput("claim2", height = "360px"),
          tableOutput("claim2_tbl"),
          hr(),
          h4("주장 3 · 같은 마비, 정반대의 약리학"),
          plotOutput("claim3", height = "360px"),
          tableOutput("claim3_tbl"),
          hr(),
          h4("주장 4 · 신장은 적분이다"),
          plotOutput("claim4", height = "340px"),
          tableOutput("claim4_tbl")
        ),

        # ---- 13 -----------------------------------------------------------
        tabPanel(
          "13 · 바이오마커 패널",
          h4("병상 검사 패널 (the panel a clinician actually orders)"),
          plotOutput("biomarker", height = "620px"),
          hr(),
          h4("시각별 값 (values at selected times)"),
          tableOutput("biomarker_tbl")
        ),

        # ---- 14 -----------------------------------------------------------
        tabPanel(
          "14 · 가상 집단",
          h4("임상시험의 한 팔은 혼합물이다 (a trial arm is a mixture)"),
          helpText("Every endpoint in this model is a THRESHOLD, and the mean of a ",
                   "threshold is not the threshold of a mean. A single median ",
                   "patient therefore cannot reproduce a trial's event rate."),
          sliderInput("npop", "피험자 수 (subjects per arm)", min = 20, max = 300,
                      value = 60, step = 20),
          actionButton("run_pop", "집단 시뮬레이션 (run population)", class = "btn-warning"),
          br(), br(),
          DT::dataTableOutput("pop_tbl"),
          hr(),
          h4("분포 (distributions)"),
          plotOutput("pop_plot", height = "400px")
        ),

        # ---- 15 -----------------------------------------------------------
        tabPanel(
          "15 · 모델 문서",
          h4("이 모델이 주장하는 것과 주장하지 않는 것"),
          verbatimTextOutput("docs")
        )
      )
    )
  )
)

# =============================================================================
#  SERVER
# =============================================================================
server <- function(input, output, session) {

  args_from_input <- function(icu = NULL, dial = NULL) {
    a <- list(
      snake = input$snake, venom_mult = input$venom_mult,
      product = if (input$product == "none") NULL else input$product,
      vials = input$vials, t_av = input$t_av, av_infuse = input$av_infuse,
      icu = if (is.null(icu)) input$icu else icu,
      dialysis = if (is.null(dial)) input$dial else dial,
      tmax = input$tmax
    )
    if (input$maint) a$maintenance <- c(input$t_av + 6, 6, 2, 3)
    if (input$neo)   a$neostigmine <- c(max(input$t_av, 1), 4, 8)
    if (input$vares) a$varespladib <- c(0.5, 12, 8)
    if (input$txa)   a$txa <- c(max(input$t_av, 1), 8, 3)
    if (input$cryo)  a$cryo <- data.frame(time = c(input$t_av + 0.5, input$t_av + 6),
                                          dFG = c(0.9, 0.9))
    if (input$fluid) a$fluids <- c(1, 8, 0.30)
    a
  }

  sim <- eventReactive(input$run, {
    withProgress(message = "integrating 50 ODEs...", {
      do.call(run_sbe, args_from_input())
    })
  }, ignoreNULL = FALSE)

  summ <- reactive({
    summarise_sbe(sim(), input$snake, input$venom_mult,
                  if (input$product == "none") NULL else input$product,
                  input$vials, input$t_av, "current")
  })

  # ---- sidebar stoichiometry box -----------------------------------------
  output$stoich_box <- renderText({
    need <- as.numeric(stoich_need(input$snake, input$venom_mult))
    if (need <= 0) return("Dry bite: nothing to neutralise.\nAll of the risk, none of the benefit.")
    mg <- if (input$product == "none") 0 else input$vials * AV[[input$product]]$mgNE_vial
    paste0(
      sprintf("need   %6.1f mgNE\n", need),
      sprintf("given  %6.1f mgNE\n", mg),
      sprintf("margin %6.2fx\n", mg / need),
      if (mg / need < 1) ">>> UNDER-DOSED <<<" else "adequate at t=0"
    )
  })

  # ---- 1 : profile --------------------------------------------------------
  output$profile <- renderTable({
    s <- SNAKES[[input$snake]]
    data.frame(
      항목 = c("종 (species)", "한글명", "주입 독량 (mg)", "SVMP / SVSP / PLA2 / 3FTx",
               "전시냅스 PLA2 분율 f_pre", "국소 괴사 강도 f_loc", "항독소",
               "바이알 · mgNE", "투여 지연 (h)"),
      값 = c(input$snake, s$ko,
             sprintf("%.0f", s$VDOSE * input$venom_mult),
             sprintf("%.2f / %.2f / %.2f / %.2f", s$F_SVMP, s$F_SVSP, s$F_PLA2, s$F_TFTX),
             sprintf("%.2f", s$F_PRE), sprintf("%.2f", s$F_LOC),
             if (input$product == "none") "none" else AV[[input$product]]$label,
             if (input$product == "none") "-" else
               sprintf("%d vials = %.0f mgNE", input$vials,
                       input$vials * AV[[input$product]]$mgNE_vial),
             sprintf("%.1f", input$t_av))
    )
  })

  output$keynum <- renderTable({
    s <- summ()
    data.frame(
      종말점 = c("피브리노겐 나디르 (g/L)", "응고불능 시간 (h)",
                 "안전계수 나디르", "인공환기 시간 (h)",
                 "국소 괴사 (%)", "CK 최고 (U/L)",
                 "크레아티닌 최고 (mg/dL) / 시각",
                 "유리 독액 AUC (mg·h/L)", "재발 (venom rebound)",
                 "P(전신출혈)", "P(사망)"),
      값 = c(sprintf("%.2f", s$FG_nadir), sprintf("%.1f", s$hrs_unclottable),
             sprintf("%.2f", s$SF_nadir), sprintf("%.0f", s$hrs_ventilated),
             sprintf("%.1f", 100 * s$NEC_final), sprintf("%.0f", s$CK_peak),
             sprintf("%.2f at %.0f h", s$SCR_peak, s$t_SCR_peak),
             sprintf("%.2f", s$VAUC),
             if (s$recurrence) sprintf("YES at %.0f h", s$recur_t) else "no",
             sprintf("%.1f %%", 100 * s$p_bleed),
             sprintf("%.1f %%", 100 * s$p_death))
    )
  })

  output$overview <- renderPlot({
    d <- sim()
    d %>%
      transmute(time,
                `free venom (mg/L)` = CVFREE,
                `antivenom (mgNE/L)` = CA_OUT,
                `fibrinogen (g/L)` = FG,
                `safety factor` = SFOUT) %>%
      pivot_longer(-time) %>%
      mutate(clock = ifelse(grepl("venom|antivenom", name),
                            "VENOM CLOCK (hours) — where antivenom acts",
                            "SUBSTRATE CLOCK (days) — where endpoints are read")) %>%
      ggplot(aes(time, value, colour = name)) +
      geom_line(linewidth = 0.9) +
      facet_wrap(~clock + name, scales = "free_y", ncol = 2) +
      labs(x = "hours post-bite", y = NULL) + THEME +
      theme(legend.position = "none")
  })

  # ---- 2 : venom PK -------------------------------------------------------
  output$pk_venom <- renderPlot({
    sim() %>%
      transmute(time,
                `fast depot (mg)` = D_SVMP + D_SVSP + D_PLA2 + D_TFTX,
                `slow / sequestered depot (mg)` = P_SVMP + P_SVSP + P_PLA2 + P_TFTX,
                `plasma free toxin (mg)` = V_SVMP + V_SVSP + V_PLA2 + V_TFTX,
                `antivenom-toxin complex (mg)` = C_SVMP + C_SVSP + C_PLA2 + C_TFTX) %>%
      pivot_longer(-time) %>%
      ggplot(aes(time, value, colour = name)) +
      geom_line(linewidth = 0.9) +
      labs(x = "hours post-bite", y = "toxin mass (mg)",
           title = "The venom leaves the body faster than it leaves the bite",
           subtitle = "the slow depot sets the apparent terminal half-life (flip-flop) — the number antivenom must outlast") +
      THEME
  })

  output$pk_class <- renderPlot({
    sim() %>%
      transmute(time,
                SVMP = V_SVMP / 3.5, SVSP = V_SVSP / 3.5,
                PLA2 = V_PLA2 / 4.0, `3FTx` = V_TFTX / 8.0) %>%
      pivot_longer(-time) %>%
      filter(value > 1e-9) %>%
      ggplot(aes(time, value, colour = name)) +
      geom_line(linewidth = 0.9) + scale_y_log10() +
      labs(x = "hours post-bite", y = "free plasma concentration (mg/L, log)",
           title = "Four toxin classes, four different absorption and clearance clocks") +
      THEME
  })

  output$pk_antigen <- renderPlot({
    sim() %>%
      transmute(time,
                `free venom (what harms)` = CVFREE,
                `total antigen (what the EIA sees)` = CVTOTAL) %>%
      pivot_longer(-time) %>%
      filter(value > 1e-9) %>%
      ggplot(aes(time, value, colour = name)) +
      geom_line(linewidth = 0.9) + scale_y_log10() +
      labs(x = "hours post-bite", y = "mg/L (log scale)",
           title = "Detectable is not the same as dangerous") +
      THEME
  })

  # ---- 3 : antivenom PK ---------------------------------------------------
  output$pk_av <- renderPlot({
    sim() %>%
      transmute(time,
                `central (mgNE)` = A_C, `peripheral (mgNE)` = A_P,
                `bite-site tissue (mgNE)` = A_T) %>%
      pivot_longer(-time) %>%
      ggplot(aes(time, value, colour = name)) +
      geom_line(linewidth = 0.9) +
      labs(x = "hours post-bite", y = "antivenom (mgNE)",
           title = "Antivenom disposition",
           subtitle = "note how little ever reaches the bite site — that is why local necrosis is not preventable intravenously") +
      THEME
  })

  output$halflives <- renderTable({
    h <- analyse_halflives()
    venom_term <- log(2) / (0.0120 + 0.0035)   # SVMP slow depot, from the model
    h %>%
      mutate(`t½ (h)` = round(t_half_beta_h, 1),
             `ratio vs venom input (44.7 h)` = round(t_half_beta_h / venom_term, 2),
             verdict = ifelse(grepl("venom", entity), "—",
                              ifelse(t_half_beta_h / venom_term >= 1,
                                     "outlasts the tail", "LOSES the race"))) %>%
      select(entity, `t½ (h)`, `ratio vs venom input (44.7 h)`, verdict)
  })

  output$stoich_plot <- renderPlot({
    analyse_stoichiometry() %>%
      mutate(snake = reorder(snake, mgNE_need)) %>%
      ggplot(aes(snake, mgNE_need)) +
      geom_col(fill = "#4A78B8") +
      geom_hline(yintercept = 10 * AV$IgG$mgNE_vial, linetype = 2, colour = "#D94A4A") +
      geom_text(aes(label = sprintf("%.0f mgNE\n(%.1f IgG vials)", mgNE_need, IgG_vials)),
                hjust = -0.05, size = 3) +
      coord_flip() + expand_limits(y = 110) +
      labs(x = NULL, y = "mgNE required to neutralise a median bite",
           title = "Naja carries LESS venom protein than Daboia and needs MORE antivenom",
           subtitle = "dashed line = the standard 10-vial Indian ASV dose (60 mgNE). 55% of cobra venom is the class antivenom covers worst (ε = 0.35).") +
      THEME
  })

  # ---- 4 : coagulation ----------------------------------------------------
  output$coag_fg <- renderPlot({
    ggplot(sim(), aes(time, FG)) +
      band(0, 0.5) +
      geom_hline(yintercept = 0.5, linetype = 2, colour = "#D94A4A") +
      geom_line(linewidth = 1) +
      annotate("text", x = Inf, y = 0.5, hjust = 1.05, vjust = -0.6, size = 3.2,
               colour = "#D94A4A", label = "20WBCT positive below 0.5 g/L") +
      labs(x = "hours post-bite", y = "fibrinogen (g/L)",
           title = "The falling limb is the venom's; the rising limb is the liver's",
           subtitle = "hepatic synthesis ceiling ≈ 1.55 g/L/day, and no antivenom variable appears in that term") +
      THEME
  })

  output$coag_other <- renderPlot({
    sim() %>%
      transmute(time, `platelets (10⁹/L)` = PLT,
                `factor X activity` = FX, `D-dimer (µg/mL)` = XDP) %>%
      pivot_longer(-time) %>%
      ggplot(aes(time, value)) + geom_line(linewidth = 0.9, colour = "#4A78B8") +
      facet_wrap(~name, scales = "free_y") +
      labs(x = "hours post-bite", y = NULL) + THEME
  })

  output$coag_wbct <- renderPlot({
    ggplot(sim(), aes(time, WBCT20)) +
      geom_area(fill = "#D94A4A", alpha = 0.5) +
      scale_y_continuous(breaks = c(0, 1), labels = c("clots", "UNCLOTTABLE")) +
      labs(x = "hours post-bite", y = NULL,
           title = "20-minute whole blood clotting test over time",
           subtitle = "a glass tube, no reagent, no electricity — and it is the gate between a dry bite and an envenoming") +
      THEME
  })

  # ---- 5 : neuromuscular --------------------------------------------------
  output$nmj_sf <- renderPlot({
    ggplot(sim(), aes(time, SFOUT)) +
      band(0, 1.3) + band(1.3, 2.2, fill = "#E8A33D") +
      geom_hline(yintercept = c(1.3, 2.2), linetype = 2) +
      geom_line(linewidth = 1) +
      annotate("text", x = Inf, y = 2.2, hjust = 1.05, vjust = -0.6, size = 3.2,
               label = "ptosis below 2.2") +
      annotate("text", x = Inf, y = 1.3, hjust = 1.05, vjust = -0.6, size = 3.2,
               colour = "#D94A4A", label = "ventilation below 1.3") +
      labs(x = "hours post-bite", y = "neuromuscular safety factor",
           title = "SF = SF₀ · TERM_eff · R_eff — a product, not a sum") +
      THEME
  })

  output$nmj_br <- renderPlot({
    ggplot(sim(), aes(time, BR)) + geom_line(linewidth = 1, colour = "#B8544A") +
      ylim(0, 1) + labs(x = "h", y = "fraction of nAChR occupied by 3FTx") + THEME
  })

  output$nmj_term <- renderPlot({
    ggplot(sim(), aes(time, TERM)) + geom_line(linewidth = 1, colour = "#4A78B8") +
      ylim(0, 1) + labs(x = "h", y = "presynaptic terminal integrity") + THEME
  })

  output$nmj_sbc <- renderPlot({
    sim() %>%
      transmute(time, `single-breath count` = SBC,
                `ptosis (1 = present)` = PTOSIS,
                `ventilation required (1 = yes)` = VENT) %>%
      pivot_longer(-time) %>%
      ggplot(aes(time, value)) + geom_line(linewidth = 0.9, colour = "#4A78B8") +
      facet_wrap(~name, scales = "free_y", ncol = 1) +
      labs(x = "hours post-bite", y = NULL) + THEME
  })

  # ---- 6 : local ----------------------------------------------------------
  output$local_nec <- renderPlot({
    sim() %>%
      transmute(time, `myonecrosis fraction` = NEC,
                `antivenom in bite-site tissue (mgNE)` = A_T,
                `antivenom in plasma (mgNE)` = A_C) %>%
      pivot_longer(-time) %>%
      ggplot(aes(time, value)) + geom_line(linewidth = 0.9, colour = "#C97A4A") +
      facet_wrap(~name, scales = "free_y", ncol = 1) +
      labs(x = "hours post-bite", y = NULL) + THEME
  })

  output$local_cp <- renderPlot({
    d <- sim()
    p1 <- d %>% transmute(time, `oedema (L)` = EDEMA, `CK (U/L)` = CK,
                          `myoglobin (µg/mL)` = MB) %>%
      pivot_longer(-time) %>%
      ggplot(aes(time, value)) + geom_line(linewidth = 0.9, colour = "#C97A4A") +
      facet_wrap(~name, scales = "free_y", ncol = 3) +
      labs(x = "h", y = NULL) + THEME
    p2 <- ggplot(d, aes(time, CPRESS)) +
      geom_hline(yintercept = 30, linetype = 2, colour = "#D94A4A") +
      geom_line(linewidth = 1) +
      labs(x = "h", y = "compartment pressure (mmHg)") + THEME
    gridExtra_or_stack(p1, p2)
  })

  # gridExtra is optional: fall back to printing the first panel only
  gridExtra_or_stack <- function(p1, p2) {
    if (requireNamespace("gridExtra", quietly = TRUE)) {
      gridExtra::grid.arrange(p1, p2, ncol = 1, heights = c(1, 1))
    } else {
      print(p1)
    }
  }

  # ---- 7 : renal ----------------------------------------------------------
  output$renal <- renderPlot({
    sim() %>%
      transmute(time, `serum creatinine (mg/dL)` = SCR,
                `eGFR (mL/min/1.73m²)` = EGFR,
                `pigment cast burden (AU)` = CAST,
                `free venom (mg/L)` = CVFREE) %>%
      pivot_longer(-time) %>%
      ggplot(aes(time, value)) + geom_line(linewidth = 0.9, colour = "#4A9B5C") +
      facet_wrap(~name, scales = "free_y") +
      labs(x = "hours post-bite", y = NULL,
           title = "Creatinine peaks days after the venom has gone") + THEME
  })

  output$renal_scar <- renderPlot({
    sim() %>%
      transmute(time, `functional nephron fraction` = NEPH,
                `recovery ceiling (1 − scar)` = 1 - FIBR) %>%
      pivot_longer(-time) %>%
      ggplot(aes(time, value, colour = name)) + geom_line(linewidth = 1) +
      ylim(0, 1) + labs(x = "hours post-bite", y = NULL,
                        title = "Permanent scar sets the ceiling recovery cannot exceed") +
      THEME
  })

  # ---- 8 : haemodynamics --------------------------------------------------
  output$haemo <- renderPlot({
    sim() %>%
      transmute(time, `plasma volume (L)` = PV, `MAP (mmHg)` = MAPOUT,
                `haematocrit (%)` = HCT, `glycocalyx integrity` = GLX) %>%
      pivot_longer(-time) %>%
      ggplot(aes(time, value)) + geom_line(linewidth = 0.9, colour = "#3D8FA8") +
      facet_wrap(~name, scales = "free_y") +
      labs(x = "hours post-bite", y = NULL) + THEME
  })

  output$infl <- renderPlot({
    sim() %>% transmute(time, `IL-6 (pg/mL)` = IL6, `TNF-α (pg/mL)` = TNF) %>%
      pivot_longer(-time) %>%
      ggplot(aes(time, value)) + geom_line(linewidth = 0.9, colour = "#7A7A7A") +
      facet_wrap(~name, scales = "free_y") + labs(x = "h", y = NULL) + THEME
  })

  # ---- 9 : endpoints ------------------------------------------------------
  output$endpoints <- renderPlot({
    sim() %>%
      transmute(time, `P(systemic bleeding)` = 100 * P_BLEED,
                `P(death)` = 100 * P_DEATH) %>%
      pivot_longer(-time) %>%
      ggplot(aes(time, value, colour = name)) + geom_line(linewidth = 1) +
      labs(x = "hours post-bite", y = "cumulative probability (%)") + THEME
  })

  output$endpoint_tbl <- renderTable({
    s <- summ()
    data.frame(
      종말점 = c("응고불능 총 시간 (h)", "안검하수 총 시간 (h)",
                 "인공환기 총 시간 (h)", "AKI 최고 크레아티닌 (mg/dL)",
                 "영구 신반흔 (%)", "국소 괴사 (%)",
                 "혈청병 최고 (day)", "P(전신출혈) %", "P(사망) %"),
      값 = c(round(s$hrs_unclottable, 1), round(s$hrs_ptosis, 0),
             round(s$hrs_ventilated, 0), round(s$SCR_peak, 2),
             round(100 * s$FIBR_final, 1), round(100 * s$NEC_final, 1),
             round(s$t_SS_peak / 24, 1), round(100 * s$p_bleed, 1),
             round(100 * s$p_death, 1))
    )
  })

  output$support_compare <- renderPlot({
    combos <- list(c(TRUE, TRUE), c(TRUE, FALSE), c(FALSE, TRUE), c(FALSE, FALSE))
    lbl <- c("ventilator + dialysis", "ventilator only", "dialysis only", "neither")
    bind_rows(lapply(seq_along(combos), function(i) {
      d <- do.call(run_sbe, args_from_input(icu = combos[[i]][1], dial = combos[[i]][2]))
      tibble(setting = lbl[i], p_death = 100 * tail(d$P_DEATH, 1))
    })) %>%
      mutate(setting = factor(setting, levels = lbl)) %>%
      ggplot(aes(setting, p_death, fill = setting)) +
      geom_col() + geom_text(aes(label = sprintf("%.1f %%", p_death)), vjust = -0.4) +
      labs(x = NULL, y = "modelled P(death) %",
           title = "Supportive care versus antivenom choice") +
      THEME + theme(legend.position = "none")
  })

  # ---- 10 : adverse events -----------------------------------------------
  output$ae <- renderPlot({
    sim() %>%
      transmute(time, `anaphylactoid activation (MCA)` = MCA,
                `circulating immune complexes (IC)` = IC,
                `serum sickness score (SS)` = SS) %>%
      pivot_longer(-time) %>%
      ggplot(aes(time, value)) + geom_line(linewidth = 0.9, colour = "#A8705C") +
      facet_wrap(~name, scales = "free_y", ncol = 1) +
      labs(x = "hours post-bite", y = NULL) + THEME
  })

  output$ae_rate <- renderPlot({
    durs <- c(0.25, 0.5, 1, 2, 4)
    bind_rows(lapply(durs, function(dd) {
      a <- args_from_input(); a$av_infuse <- dd; a$tmax <- 72
      d <- do.call(run_sbe, a)
      tibble(`infusion duration (h)` = dd, `peak anaphylactoid index` = max(d$MCA))
    })) %>%
      ggplot(aes(`infusion duration (h)`, `peak anaphylactoid index`)) +
      geom_hline(yintercept = 1, linetype = 2, colour = "#D94A4A") +
      geom_line(linewidth = 1) + geom_point(size = 2.5) +
      labs(title = "The same total dose, five infusion rates",
           subtitle = "dashed line = threshold for a clinical reaction. Nothing about the vial count changed.") +
      THEME
  })

  # ---- 11 : scenarios -----------------------------------------------------
  allres <- eventReactive(input$run_all, {
    withProgress(message = "running 27 scenarios...", run_all())
  })

  output$scen_tbl <- DT::renderDataTable({
    bind_rows(lapply(allres(), `[[`, "sum")) %>%
      transmute(scenario = name, product, vials, `t_av (h)` = t_av,
                `margin` = round(molar_margin, 2),
                `FG nadir` = round(FG_nadir, 2),
                `unclot h` = round(hrs_unclottable, 1),
                `recur` = ifelse(recurrence, sprintf("%.0f h", recur_t), "-"),
                `SF nadir` = round(SF_nadir, 2),
                `vent h` = round(hrs_ventilated, 0),
                `NEC %` = round(100 * NEC_final, 1),
                `SCr peak` = round(SCR_peak, 2),
                `P(bleed) %` = round(100 * p_bleed, 1),
                `P(death) %` = round(100 * p_death, 1)) %>%
      DT::datatable(options = list(pageLength = 27, scrollX = TRUE), rownames = FALSE)
  })

  output$scen_overlay <- renderPlot({
    r <- allres()
    k <- names(r)[1:8]
    bind_rows(lapply(k, function(n) r[[n]]$sim %>% mutate(scenario = n))) %>%
      rename(y = !!input$ov_var) %>%
      ggplot(aes(time, y, colour = scenario)) + geom_line(linewidth = 0.8) +
      labs(x = "hours post-bite", y = input$ov_var) + THEME +
      guides(colour = guide_legend(ncol = 2))
  })

  # ---- 12 : the four claims ----------------------------------------------
  output$claim1 <- renderPlot({ plot_claim1(allres()) })
  output$claim1_tbl <- renderTable({ analyse_depth_vs_slope(allres()) })
  output$claim2 <- renderPlot({ plot_claim2(allres()) })
  output$claim2_tbl <- renderTable({ analyse_recurrence(allres()) })
  output$claim3 <- renderPlot({ plot_claim3(allres()) })
  output$claim3_tbl <- renderTable({ analyse_two_paralyses(allres()) })
  output$claim4 <- renderPlot({
    r <- allres()
    bind_rows(lapply(r, `[[`, "sum")) %>%
      filter(grepl("^0[1-5] Russell", name)) %>%
      ggplot(aes(VAUC, SCR_peak, label = name)) +
      geom_point(size = 3, colour = "#4A9B5C") +
      geom_text(hjust = -0.05, size = 3) + expand_limits(x = 130) +
      labs(x = "free venom AUC (mg·h/L) — the integral antivenom reduces",
           y = "peak serum creatinine (mg/dL) — read days later",
           title = "The kidney reads the integral, not the level") + THEME
  })
  output$claim4_tbl <- renderTable({ analyse_kidney_integral(allres()) })

  # ---- 13 : biomarkers ----------------------------------------------------
  output$biomarker <- renderPlot({
    sim() %>%
      transmute(time,
                `Fibrinogen (g/L)` = FG, `Platelets (10⁹/L)` = PLT,
                `D-dimer (µg/mL)` = XDP, `Factor X` = FX,
                `CK (U/L)` = CK, `Myoglobin (µg/mL)` = MB,
                `Creatinine (mg/dL)` = SCR, `eGFR` = EGFR,
                `Haematocrit (%)` = HCT, `MAP (mmHg)` = MAPOUT,
                `IL-6 (pg/mL)` = IL6, `Venom antigen (mg/L)` = CVTOTAL) %>%
      pivot_longer(-time) %>%
      ggplot(aes(time, value)) + geom_line(linewidth = 0.85, colour = "#3F5D8C") +
      facet_wrap(~name, scales = "free_y", ncol = 3) +
      labs(x = "hours post-bite", y = NULL) + THEME
  })

  output$biomarker_tbl <- renderTable({
    d <- sim()
    ts <- c(0, 2, 6, 12, 24, 48, 72, 120, 168, 240, 336)
    ts <- ts[ts <= max(d$time)]
    pick <- function(v) sapply(ts, function(t) approx(d$time, v, t)$y)
    data.frame(
      `h` = ts,
      Fibrinogen = round(pick(d$FG), 2),
      Platelets = round(pick(d$PLT), 0),
      `D-dimer` = round(pick(d$XDP), 1),
      CK = round(pick(d$CK), 0),
      Creatinine = round(pick(d$SCR), 2),
      `Safety factor` = round(pick(d$SFOUT), 2),
      `Free venom` = signif(pick(d$CVFREE), 3),
      Antivenom = signif(pick(d$CA_OUT), 3),
      check.names = FALSE
    )
  })

  # ---- 14 : population ---------------------------------------------------
  pop <- eventReactive(input$run_pop, {
    withProgress(message = "running virtual population...",
                 virtual_population(n = input$npop))
  })

  output$pop_tbl <- DT::renderDataTable({
    pop() %>% mutate(across(where(is.numeric), ~round(.x, 1))) %>%
      DT::datatable(options = list(pageLength = 10, scrollX = TRUE), rownames = FALSE)
  })

  output$pop_plot <- renderPlot({
    pop() %>%
      select(arm, pct_unclottable, pct_recurrence, pct_AKI3,
             pct_ventilated, mean_p_death, pct_underdosed) %>%
      pivot_longer(-arm) %>%
      ggplot(aes(arm, value, fill = name)) +
      geom_col(position = "dodge") + coord_flip() +
      labs(x = NULL, y = "%",
           title = "Event rates are properties of a MIXTURE, not of a median patient") +
      THEME
  })

  # ---- 15 : documentation ------------------------------------------------
  output$docs <- renderText(paste(collapse = "\n", c(
    "이 모델이 주장하는 것 (WHAT THIS MODEL CLAIMS)",
    "==============================================",
    "",
    "  ANTIVENOM BINDS.  IT DOES NOT UNDO.",
    "",
    "  clock 1 (hours)  the VENOM clock     — antivenom acts here",
    "  clock 2 (days)   the SUBSTRATE clock — resynthesis / regeneration",
    "",
    "  Every clinical endpoint is a read-out of clock 2, driven by the",
    "  time-integral of clock 1.  Four consequences follow arithmetically:",
    "",
    "  1  Antivenom timing sets the DEPTH of the nadir, not the SLOPE of",
    "     recovery.  The hepatic fibrinogen synthesis term contains no",
    "     antivenom variable at all.",
    "  2  Recurrence is a ratio of two half-lives.  When the fragment's",
    "     terminal half-life is shorter than the venom INPUT half-life, free",
    "     venom must reappear, and no larger front-loaded dose prevents it.",
    "  3  Presynaptic and postsynaptic paralysis look identical and answer",
    "     opposite questions.  Neither antivenom nor neostigmine appears in",
    "     the dTERM/dt equation.",
    "  4  Acute kidney injury is an integral, not a level.  Creatinine peaks",
    "     on day 3-5, long after venom is undetectable.",
    "",
    "이 모델이 주장하지 않는 것 (WHAT IT DOES NOT CLAIM)",
    "===================================================",
    "",
    "  * The ε vector (SVMP 1.20 / SVSP 1.00 / PLA2 0.60 / 3FTx 0.35) is a",
    "    MODEL PARAMETER, not a measured potency.  The direction is well",
    "    supported by antivenomics; the four numbers are not.",
    "  * The slow-depot compartmentalisation (f_slow, ka_s) is a structural",
    "    hypothesis.  Prolonged venom antigenaemia is observed; '35% of the",
    "    dose sequestered with a 45 h release half-life' is not.",
    "  * The model reproduces recurrent venom ANTIGENAEMIA after a short-lived",
    "    fragment from absorption kinetics alone.  It does NOT reproduce severe",
    "    recurrent HYPOFIBRINOGENAEMIA: that needs a late venom flux ~17x",
    "    larger than a sequestered depot can supply, and neither a larger",
    "    sequestered fraction nor complex dissociation gets there.  This is",
    "    reported as a negative result rather than tuned away.",
    "  * All hazard-function baselines are calibrated, not measured.",
    "",
    "검증 (VERIFICATION)",
    "===================",
    "",
    "  The build environment had no R runtime, so every equation was",
    "  independently re-implemented in Python (sbe_reference_model.py) and",
    "  actually integrated.  Integration exposed and fixed eight real defects,",
    "  each marked NOTE(defect n) in both files — including one that was",
    "  mechanistic rather than numerical (neostigmine written as a free-standing",
    "  multiplier on the safety factor, which made it rescue a presynaptic",
    "  krait bite as efficiently as a postsynaptic cobra bite, i.e. destroyed",
    "  claim 3 while looking perfectly plausible).",
    "",
    "  All 145 PubMed references were machine-verified against NCBI",
    "  E-utilities; none was written from memory.",
    "",
    "면책 (DISCLAIMER)",
    "=================",
    "",
    "  Educational and research use only.  Not validated for clinical",
    "  decision-making, prescribing, or regulatory submission.  Snakebite is",
    "  a medical emergency: follow local WHO/national treatment guidelines and",
    "  the antivenom manufacturer's label, not this model."
  )))
}

shinyApp(ui, server)
