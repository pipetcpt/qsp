## ===========================================================================
##  osa_shiny_app.R
##  Osteosarcoma QSP model — interactive dashboard
##  골육종 QSP 모델 인터랙티브 대시보드
##
##  The app is organised around the one thing the model is for: making the
##  exposure BUDGET and its three organ CEILINGS visible at the same time as
##  the cure probability they buy.  Every tab is a different face of
##
##      P(cure) = exp( -lambda0 * [1 - exp(-n0 * e^-K)] )
##      overall survival = P(cure) x (1 - TRM)
##
##  Tab 1  환자 프로파일     patient, tumour and regimen builder
##  Tab 2  약동학 (PK)       methotrexate / doxorubicin / cisplatin profiles
##  Tab 3  신장 피드백 루프  the bifurcation: C_tub vs S(pH), critical pH
##  Tab 4  종양 · PD         viable vs necrotic primary, log-kill accrual
##  Tab 5  임상 엔드포인트   Huvos, cure, TRM, survival, LVEF, hearing, eGFR
##  Tab 6  시나리오 비교     all 20 arms, ranked by SURVIVAL not by kill
##  Tab 7  바이오마커 · 독성 ALP, ANC, mucositis, the go/no-go gate log
##  Tab 8  가상 집단         the detectability wall, and where the budget belongs
##
##  Run:  shiny::runApp("osa_shiny_app.R")
##  Requires osa_mrgsolve_model.R in the same directory.
## ===========================================================================

library(shiny)
library(mrgsolve)
library(ggplot2)
library(dplyr)
library(tidyr)

source("osa_mrgsolve_model.R", local = TRUE)

THEME <- theme_minimal(base_size = 12) +
  theme(panel.grid.minor = element_blank(),
        strip.background = element_rect(fill = "#f0f0f0", colour = NA),
        legend.position = "bottom")

CEIL <- "#c00000"      # colour of every organ ceiling in the app
GOOD <- "#1f7a4d"      # colour of every protective term

## ===========================================================================
##  UI
## ===========================================================================
ui <- fluidPage(
  titlePanel(
    div(
      h3("골육종 QSP 모델 · Osteosarcoma Quantitative Systems Pharmacology"),
      p(style = "color:#555; font-size:13px; margin-top:-6px",
        HTML("Cure is a Poisson bet on lesions nobody can see, paid for out of an
             exposure budget with three hard organ ceilings.<br/>
             <b>P(cure) = exp(&minus;&lambda;&#8320;&middot;[1 &minus; exp(&minus;n&#8320;e<sup>&minus;K</sup>)])</b>
             &nbsp;&middot;&nbsp; &lambda;&#8320; = 1.80 occult lesions
             &nbsp;&middot;&nbsp; methotrexate is cleared through the organ it destroys"))
    )
  ),

  sidebarLayout(
    sidebarPanel(
      width = 3,
      h4("① 환자 · 종양 (Patient & tumour)"),
      sliderInput("bsa", "BSA (m²)", 1.0, 2.2, 1.7, 0.05),
      sliderInput("prim0", "Primary tumour volume at diagnosis (mL)",
                  30, 600, 150, 10),
      sliderInput("res0", "Intrinsic resistant fraction RES₀",
                  0, 0.45, 0.02, 0.01),
      helpText(HTML("RES₀ is what 'poor responder' <i>means</i>: every kill
                     term carries (1 &minus; RES).")),
      sliderInput("lambda0", "λ₀ — occult lesions at diagnosis (Poisson mean)",
                  0, 15, 1.80, 0.1),
      helpText(HTML("λ₀ = 1.80 is fixed by one number: cure after amputation
                     alone is exp(&minus;λ₀) = 0.165.")),
      sliderInput("kid0", "Starting kidney function (fraction)",
                  0.4, 1.0, 1.0, 0.05),

      hr(),
      h4("② 지지요법 (Supportive care — the real lever)"),
      sliderInput("ph", "Urine pH", 5.5, 8.0, 7.5, 0.1),
      sliderInput("hydr", "Hydration (L/m²/day)", 0.75, 6.0, 3.0, 0.25),
      helpText(HTML("S(pH) = 0.86&middot;10<sup>0.682(pH&minus;5)</sup> mM.
                     One pH unit is worth <b>4.81×</b> the fluid, or a
                     4.81-fold dose reduction.")),
      checkboxInput("gcsf", "G-CSF support", FALSE),

      hr(),
      h4("③ 레지멘 (Regimen)"),
      selectInput("scenario", "Preset scenario",
                  choices = names(SCENARIOS), selected = "S02_MAP"),
      sliderInput("mtxdose", "Methotrexate (g/m² per course)", 0, 18, 12, 1),
      sliderInput("doxdose", "Doxorubicin (mg/m² per course)", 0, 120, 75, 5),
      sliderInput("cisdose", "Cisplatin (mg/m² per course)", 0, 150, 120, 10),
      checkboxInput("dexra", "Dexrazoxane 10:1 (raises the cardiac ceiling)",
                    FALSE),
      checkboxInput("ie", "Add ifosfamide/etoposide (MAPIE escalation)", FALSE),
      numericInput("glucat", "Glucarpidase at (h after each MTX; NA = none)",
                   NA, min = 6, max = 96, step = 6),

      hr(),
      h4("④ 집단 (Population)"),
      sliderInput("cv", "CV of achieved log-kill between patients",
                  0.02, 0.65, 0.25, 0.01),
      helpText(HTML("This slider, not the dose sliders, is where the forty-year
                     plateau lives. Watch tab 8."))
    ),

    mainPanel(
      width = 9,
      tabsetPanel(
        id = "tabs", type = "tabs",

        ## ---------------------------------------------------------------
        tabPanel(
          "① 환자 프로파일",
          br(),
          fluidRow(
            column(4, wellPanel(h4("Exposure budget delivered"),
                                tableOutput("budget"))),
            column(4, wellPanel(h4("Organ ceilings reached"),
                                tableOutput("ceilings"))),
            column(4, wellPanel(h4("Outcome"), tableOutput("outcome")))
          ),
          h4("Where the log-kill came from"),
          plotOutput("attribution", height = "280px"),
          helpText(HTML("Twelve methotrexate courses and six doxorubicin
            courses are worth comparable <i>totals</i>, but the per-course
            value differs several-fold — and each is limited by a different
            organ, so they are not interchangeable currency."))
        ),

        ## ---------------------------------------------------------------
        tabPanel(
          "② 약동학 (PK)",
          br(),
          h4("High-dose methotrexate — a single course"),
          plotOutput("pk_mtx", height = "320px"),
          helpText(HTML("Dashed lines are the protocol monitoring thresholds:
            <b>C24 &lt; 10 µM, C48 &lt; 1 µM, C72 &lt; 0.1 µM</b>. Crossing any
            of them is what triggers rescue escalation.")),
          h4("Doxorubicin, doxorubicinol and cisplatin over the whole protocol"),
          plotOutput("pk_other", height = "300px")
        ),

        ## ---------------------------------------------------------------
        tabPanel(
          "③ 신장 피드백 루프",
          br(),
          h4("The loop: C_tub = CL_ren(KID)·C_p / UF(KID) versus S(pH)"),
          plotOutput("loop_ts", height = "300px"),
          helpText(HTML("Injury lowers <i>both</i> the clearance and the urine
            flow, which <i>raises</i> the concentration that caused it. Positive
            feedback, so there is a critical pH rather than a risk factor.")),
          fluidRow(
            column(6,
                   h4("The bifurcation"),
                   plotOutput("crit_ph", height = "300px")),
            column(6,
                   h4("The exchange rate"),
                   plotOutput("exchange", height = "300px"),
                   helpText(HTML("C_tub scales as dose/flow; S scales as
                     10<sup>0.682·pH</sup>. Any factor f on C_tub is worth
                     log₁₀(f)/0.682 pH units — so doubling the fluid and
                     halving the dose are both worth 0.44 pH units.")))
          )
        ),

        ## ---------------------------------------------------------------
        tabPanel(
          "④ 종양 · PD",
          br(),
          h4("Primary tumour: viable versus necrotic volume"),
          plotOutput("tumour", height = "300px"),
          helpText(HTML("Osteosarcoma volume barely shrinks — the kill converts
            viable tumour into <i>retained</i> necrotic osteoid, which is
            exactly what the pathologist grades at week 11.")),
          fluidRow(
            column(6, h4("Micrometastatic log-kill accrual"),
                   plotOutput("logkill", height = "280px")),
            column(6, h4("Resistant fraction under selection"),
                   plotOutput("resist", height = "280px"))
          )
        ),

        ## ---------------------------------------------------------------
        tabPanel(
          "⑤ 임상 엔드포인트",
          br(),
          fluidRow(
            column(3, wellPanel(h4("Huvos"), h2(textOutput("ep_huvos")),
                                p("necrosis % at week 11; ≥90% = good"))),
            column(3, wellPanel(h4("Cure"), h2(textOutput("ep_cure")),
                                p("population cure fraction"))),
            column(3, wellPanel(h4("TRM"), h2(textOutput("ep_trm")),
                                p("treatment-related mortality"))),
            column(3, wellPanel(style = "background:#eef7f0",
                                h4("Survival"), h2(textOutput("ep_surv")),
                                p("cure × (1 − TRM) — read THIS one")))
          ),
          h4("The three irreversible ceilings over time"),
          plotOutput("ceil_ts", height = "320px"),
          helpText(HTML("Heart, cochlea and tubule are cumulative and
            essentially irreversible. They are the walls of the box, and
            standard MAP already stands against all three.")),
          h4("The double exponential"),
          plotOutput("double_exp", height = "280px"),
          helpText(HTML("Survival of one lesion is 1 &minus; exp(&minus;n₀e<sup>&minus;K</sup>);
            cure needs <i>every</i> lesion to fail. That is why the same
            log-kill means something completely different at λ₀ = 1.8 and at
            λ₀ = 15."))
        ),

        ## ---------------------------------------------------------------
        tabPanel(
          "⑥ 시나리오 비교",
          br(),
          h4("All 20 arms — ranked by SURVIVAL, not by log-kill"),
          plotOutput("scen_plot", height = "460px"),
          helpText(HTML("Several arms buy log-kill by damaging the kidney. The
            kill they buy is real, and they still lose. Under-hydration buys
            +0.4 log₁₀ and pays 12% treatment-related mortality, for survival
            0.595 against MAP's 0.596 — the regimen sits on an
            <b>indifference curve</b>, so moves along it read as null however
            large they are. Only dexrazoxane moves the wall.")),
          h4("Full table"),
          tableOutput("scen_table")
        ),

        ## ---------------------------------------------------------------
        tabPanel(
          "⑦ 바이오마커 · 독성",
          br(),
          fluidRow(
            column(6, h4("Neutrophils (Friberg transit model)"),
                   plotOutput("anc", height = "260px")),
            column(6, h4("Mucositis grade"),
                   plotOutput("muc", height = "260px"))
          ),
          fluidRow(
            column(6, h4("Alkaline phosphatase (× ULN)"),
                   plotOutput("alp", height = "260px")),
            column(6, h4("Osteoclast drive and matrix TGF-β"),
                   plotOutput("bone", height = "260px"))
          ),
          h4("The go/no-go gate log — where escalation spends the budget"),
          tableOutput("gatelog"),
          helpText(HTML("MTX needs ANC ≥ 0.75, mucositis &lt; 3 and kidney
            ≥ 0.60; AP and IE need ANC ≥ 1.0; doxorubicin needs LVEF ≥ 50%.
            A deferred course returns two weeks later at &minus;25% dose, and is
            dropped after that. This is the mechanism, not an assumption, by
            which MAPIE delivers less of the backbone that was doing the
            killing."))
        ),

        ## ---------------------------------------------------------------
        tabPanel(
          "⑧ 가상 집단",
          br(),
          h4("The detectability wall"),
          plotOutput("detect", height = "320px"),
          helpText(HTML("At the real spread of achieved log-kill (CV ≈ 25%) a
            realistic +0.2 log intensification is worth about two points of
            survival and needs <b>2764 patients per arm</b>. EURAMOS-1
            randomised 618. The trial was not underpowered by accident; it was
            underpowered by the shape of the dose-response surface it was
            built on top of.")),
          fluidRow(
            column(6, h4("What the spread does to the value of +0.5 log"),
                   plotOutput("spread", height = "300px")),
            column(6, h4("Where the budget belongs — the shoulder"),
                   plotOutput("tailplot", height = "300px"),
                   helpText(HTML("Same total extra log-kill, different
                     recipient. P(cure) saturates at <i>both</i> ends, so the
                     extra log is wasted on patients already cured
                     <i>and</i> on patients no achievable dose can rescue. The
                     second quartile — the near-misses at roughly 85-89%
                     necrosis — is worth <b>2.4×</b> uniform intensification,
                     and the worst quartile is worth <i>less</i> than uniform.")))
          ),
          h4("Marginal value of an extra log, by decile of achieved kill"),
          plotOutput("popdist", height = "300px")
        )
      )
    )
  )
)

## ===========================================================================
##  SERVER
## ===========================================================================
server <- function(input, output, session) {

  ## Build the regimen from the preset plus the manual overrides -----------
  reg <- reactive({
    r <- SCENARIOS[[input$scenario]]
    r$mtx_dose    <- input$mtxdose
    r$dox_dose    <- input$doxdose
    r$cis_dose    <- input$cisdose
    r$dexrazoxane <- input$dexra
    if (input$ie && is.null(r$ifo_weeks)) r$ifo_weeks <- c(14, 18, 23, 28, 32)
    if (!is.na(input$glucat)) r$gluc_at <- input$glucat
    r$p <- c(r$p, list(
      BSA      = input$bsa,
      PRIM0    = input$prim0,
      RES0     = input$res0,
      LAMBDA0  = input$lambda0,
      URINE_PH = input$ph,
      UF0      = 0.21 * input$hydr / 3.0,
      GCSF     = if (input$gcsf) 1.5 else 0.0))
    r
  })

  run <- reactive({
    withProgress(message = "Integrating 55 ODEs over 40 weeks...", {
      simulate_osa(reg())
    })
  })

  fin <- reactive({ s <- run()$sim; s[nrow(s), ] })

  cure <- reactive({
    lam <- rpois(POP_N, input$lambda0)
    pop_cure(fin()$KNATS, cv = input$cv, lam = lam)
  })
  surv <- reactive(cure() * (1 - fin()$PTRM))

  ## ---- Tab 1 -----------------------------------------------------------
  output$budget <- renderTable({
    g <- run()$given; p <- run()$planned
    data.frame(
      agent = c("Methotrexate", "Doxorubicin", "Cisplatin", "Ifos/Etop"),
      given = sprintf("%d / %d", g, p),
      exposure = sprintf("%.1f", c(fin()$EMTX, fin()$EDOX, fin()$ECIS,
                                   fin()$EIFO)))
  }, striped = TRUE)

  output$ceilings <- renderTable({
    s <- run()$sim
    data.frame(
      organ = c("Heart (LVEF %)", "Cochlea (dB shift)",
                "Tubule (eGFR nadir %)", "Marrow (ANC nadir)",
                "Cumulative doxorubicin"),
      value = c(sprintf("%.1f", fin()$LVEFo),
                sprintf("%.1f", fin()$HLo),
                sprintf("%.0f", 100 * min(s$eGFRo) / 120),
                sprintf("%.2f", min(s$CIRC)),
                sprintf("%.0f mg/m2", fin()$CUMDOX)))
  }, striped = TRUE)

  output$outcome <- renderTable({
    data.frame(
      metric = c("Huvos necrosis %", "log10 kill (micromets)",
                 "Cure fraction", "TRM", "Overall survival",
                 "Courses on schedule"),
      value = c(sprintf("%.1f", run()$huvos),
                sprintf("%.2f", fin()$KLOG10),
                sprintf("%.3f", cure()),
                sprintf("%.1f%%", 100 * fin()$PTRM),
                sprintf("%.3f", surv()),
                sprintf("%.0f%%", 100 * run()$on_schedule)))
  }, striped = TRUE)

  output$attribution <- renderPlot({
    f <- fin()
    d <- data.frame(
      agent = factor(c("Methotrexate", "Doxorubicin", "Cisplatin",
                       "Ifos/Etoposide", "Immune (lung)"),
                     levels = c("Methotrexate", "Doxorubicin", "Cisplatin",
                                "Ifos/Etoposide", "Immune (lung)")),
      log10 = c(f$KI_MTX, f$KI_DOX, f$KI_CIS, f$KI_IE, f$KI_IMM) / log(10))
    ggplot(d, aes(agent, log10, fill = agent)) +
      geom_col(width = 0.65, show.legend = FALSE) +
      geom_text(aes(label = sprintf("%.2f", log10)), vjust = -0.4, size = 4) +
      labs(x = NULL, y = expression(log[10]~"cell kill delivered")) + THEME
  })

  ## ---- Tab 2 -----------------------------------------------------------
  output$pk_mtx <- renderPlot({
    r <- mtx_single_course(pH = input$ph, dose_g_m2 = input$mtxdose,
                           uf_mult = input$hydr / 3, kid0 = input$kid0,
                           gluc_at = input$glucat, hours = 120)
    d <- data.frame(t = c(24, 48, 72),
                    thr = c(10, 1, 0.1),
                    obs = c(r$C24, r$C48, r$C72))
    ggplot(d, aes(t, obs)) +
      geom_hline(aes(yintercept = thr), linetype = "dashed",
                 colour = CEIL) +
      geom_point(size = 4, colour = "#1f4d63") +
      geom_text(aes(label = sprintf("%.2f µM", obs)), vjust = -1, size = 4) +
      scale_y_log10() +
      labs(x = "hours after the start of the 4-h infusion",
           y = "plasma methotrexate (µM, log scale)",
           subtitle = sprintf(
             "Cmax %.0f µM · supersaturation ratio %.2f · %s · %s",
             r$Cmax, r$SS_max,
             if (r$delayed) "DELAYED ELIMINATION" else "elimination normal",
             if (r$aki) "AKI" else "no AKI")) + THEME
  })

  output$pk_other <- renderPlot({
    run()$sim %>%
      select(time, Doxorubicin = CDOXo, Cisplatin = CCISo,
             Etoposide = CETOo) %>%
      pivot_longer(-time) %>%
      ggplot(aes(time / 168, value, colour = name)) +
      geom_line(linewidth = 0.5) +
      facet_wrap(~name, scales = "free_y", ncol = 3) +
      labs(x = "week", y = "mg/L", colour = NULL) + THEME
  })

  ## ---- Tab 3 -----------------------------------------------------------
  output$loop_ts <- renderPlot({
    run()$sim %>%
      filter(time <= 6 * 168) %>%
      select(time, `tubular C_tub (mM)` = CTUBo,
             `solubility S(pH) (mM)` = SOLo) %>%
      pivot_longer(-time) %>%
      ggplot(aes(time / 168, value, colour = name)) +
      geom_line(linewidth = 0.7) +
      scale_colour_manual(values = c("tubular C_tub (mM)" = CEIL,
                                     "solubility S(pH) (mM)" = GOOD)) +
      labs(x = "week", y = "mM", colour = NULL,
           subtitle = "Crystal forms wherever the red line rises above the green one") +
      THEME
  })

  output$crit_ph <- renderPlot({
    phs <- seq(5.8, 8.0, by = 0.1)
    d <- do.call(rbind, lapply(phs, function(p) {
      r <- mtx_single_course(pH = p, dose_g_m2 = input$mtxdose,
                             uf_mult = input$hydr / 3, hours = 168)
      data.frame(pH = p, SS = r$SS_max, kid = 100 * r$kid_final)
    }))
    ggplot(d, aes(pH, SS)) +
      geom_hline(yintercept = 1, linetype = "dashed", colour = CEIL) +
      geom_line(linewidth = 0.8, colour = "#1f4d63") +
      annotate("text", x = 6.1, y = 1.15, label = "supersaturated above 1",
               colour = CEIL, size = 3.4, hjust = 0) +
      scale_y_log10() +
      labs(x = "urine pH", y = "peak supersaturation ratio C_tub / S",
           subtitle = "The critical pH is where this curve crosses 1") + THEME
  })

  output$exchange <- renderPlot({
    d <- data.frame(
      factor_on_Ctub = c(0.25, 0.5, 1, 2, 4),
      label = c("¼ the dose / 4× fluid", "½ the dose / 2× fluid",
                "reference", "2× the dose / ½ fluid",
                "4× the dose / ¼ fluid"))
    d$pH_equivalent <- -log10(d$factor_on_Ctub) / 0.682
    ggplot(d, aes(reorder(label, pH_equivalent), pH_equivalent)) +
      geom_col(fill = "#d8b44a", width = 0.6) +
      geom_text(aes(label = sprintf("%+.2f pH", pH_equivalent)),
                hjust = ifelse(d$pH_equivalent >= 0, -0.1, 1.1), size = 3.6) +
      coord_flip() +
      labs(x = NULL, y = "equivalent shift in urine pH") + THEME
  })

  ## ---- Tab 4 -----------------------------------------------------------
  output$tumour <- renderPlot({
    run()$sim %>%
      select(time, Viable = PRIM, Necrotic = NEC) %>%
      pivot_longer(-time) %>%
      ggplot(aes(time / 168, value, fill = name)) +
      geom_area(alpha = 0.75) +
      geom_vline(xintercept = 10, linetype = "dashed") +
      annotate("text", x = 10.2, y = Inf, label = "surgery, week 11",
               vjust = 1.5, hjust = 0, size = 3.4) +
      scale_fill_manual(values = c(Viable = "#e0b681", Necrotic = "#9a9a9a")) +
      labs(x = "week", y = "volume (mL)", fill = NULL) + THEME
  })

  output$logkill <- renderPlot({
    run()$sim %>%
      ggplot(aes(time / 168, KLOG10)) +
      geom_line(linewidth = 0.8, colour = "#2d2d6b") +
      labs(x = "week", y = expression(log[10]~"kill, micrometastatic pool"),
           subtitle = sprintf("MAP needs 6.6 log10 to reach cure 0.60; this run delivers %.2f",
                              fin()$KLOG10)) + THEME
  })

  output$resist <- renderPlot({
    run()$sim %>%
      ggplot(aes(time / 168, RES)) +
      geom_line(linewidth = 0.8, colour = CEIL) +
      labs(x = "week", y = "resistant clone fraction",
           subtitle = "Selection under cytotoxic pressure caps the deliverable log-kill") +
      THEME
  })

  ## ---- Tab 5 -----------------------------------------------------------
  output$ep_huvos <- renderText(sprintf("%.1f%%", run()$huvos))
  output$ep_cure  <- renderText(sprintf("%.3f", cure()))
  output$ep_trm   <- renderText(sprintf("%.1f%%", 100 * fin()$PTRM))
  output$ep_surv  <- renderText(sprintf("%.3f", surv()))

  output$ceil_ts <- renderPlot({
    run()$sim %>%
      transmute(week = time / 168,
                `Heart — LVEF (%)` = LVEFo,
                `Cochlea — hearing loss (dB)` = HLo,
                `Tubule — eGFR (mL/min)` = eGFRo) %>%
      pivot_longer(-week) %>%
      ggplot(aes(week, value)) +
      geom_line(linewidth = 0.8, colour = CEIL) +
      facet_wrap(~name, scales = "free_y", ncol = 3) +
      labs(x = "week", y = NULL) + THEME
  })

  output$double_exp <- renderPlot({
    K <- seq(0, 20, by = 0.1)
    d <- expand.grid(K = K, lambda0 = c(0.5, 1.8, 6, 15)) %>%
      mutate(cure = exp(-lambda0 * (1 - exp(-1e6 * exp(-K)))),
             lambda0 = factor(lambda0))
    ggplot(d, aes(K / log(10), cure, colour = lambda0)) +
      geom_line(linewidth = 0.8) +
      geom_vline(xintercept = fin()$KLOG10, linetype = "dashed") +
      annotate("text", x = fin()$KLOG10 + 0.1, y = 0.05,
               label = "this run", hjust = 0, size = 3.4) +
      labs(x = expression(log[10]~"cell kill K"), y = "P(cure)",
           colour = expression(lambda[0])) + THEME
  })

  ## ---- Tab 6 -----------------------------------------------------------
  scen_res <- reactive({
    withProgress(message = "Running all 20 scenarios...", run_all())
  })

  output$scen_plot <- renderPlot({
    scen_res() %>%
      select(scenario, `log10 kill` = log10_kill, `cure` = cure,
             `survival` = survival) %>%
      pivot_longer(-scenario) %>%
      mutate(name = factor(name,
                           levels = c("log10 kill", "cure", "survival"))) %>%
      ggplot(aes(reorder(scenario,
                         ave(value, name, FUN = function(x) x)),
                 value, fill = name)) +
      geom_col(show.legend = FALSE) +
      facet_wrap(~name, scales = "free_x", ncol = 3) +
      coord_flip() +
      labs(x = NULL, y = NULL) + THEME
  })

  output$scen_table <- renderTable({
    scen_res() %>%
      transmute(scenario, Huvos = huvos_pct, log10 = log10_kill,
                cure, `TRM %` = 100 * TRM, survival,
                `EFS HR` = EFS_HR, `on sched %` = on_sched,
                `eGFR nadir %` = eGFR_nadir, LVEF, dB = hearing_dB,
                `ANC nadir` = ANC_nadir)
  }, digits = 3, striped = TRUE)

  ## ---- Tab 7 -----------------------------------------------------------
  output$anc <- renderPlot({
    ggplot(run()$sim, aes(time / 168, CIRC)) +
      geom_hline(yintercept = c(0.5, 1.0), linetype = "dashed",
                 colour = c(CEIL, "#888888")) +
      geom_line(linewidth = 0.6, colour = "#a86f49") +
      labs(x = "week", y = expression(ANC~(10^9/L)),
           subtitle = "dashed: grade-4 threshold 0.5 and the AP/IE gate 1.0") +
      THEME
  })

  output$muc <- renderPlot({
    ggplot(run()$sim, aes(time / 168, MUCo)) +
      geom_hline(yintercept = 3, linetype = "dashed", colour = CEIL) +
      geom_line(linewidth = 0.6, colour = "#a86f49") +
      ylim(0, 4) +
      labs(x = "week", y = "mucositis grade",
           subtitle = "dashed: the gate — grade 3 defers every course") + THEME
  })

  output$alp <- renderPlot({
    ggplot(run()$sim, aes(time / 168, ALP)) +
      geom_line(linewidth = 0.6, colour = "#b08a20") +
      labs(x = "week", y = "ALP (× ULN)") + THEME
  })

  output$bone <- renderPlot({
    run()$sim %>%
      select(time, `osteoclast activity` = OCL,
             `matrix TGF-β` = TGFM) %>%
      pivot_longer(-time) %>%
      ggplot(aes(time / 168, value, colour = name)) +
      geom_line(linewidth = 0.6) +
      labs(x = "week", y = NULL, colour = NULL,
           subtitle = "TGF-β multiplies GROWTH, not survival — which is why bone-directed drugs move this panel and not the cure") +
      THEME
  })

  output$gatelog <- renderTable({
    r <- run()
    data.frame(
      agent   = names(r$planned),
      planned = as.integer(r$planned),
      given   = as.integer(r$given),
      `on schedule` = as.integer(r$on_time),
      omitted = as.integer(r$omitted),
      check.names = FALSE)
  }, striped = TRUE)

  ## ---- Tab 8 -----------------------------------------------------------
  output$detect <- renderPlot({
    d <- detectability(fin()$KNATS, cv = input$cv)
    ggplot(d, aes(delta_log10, n_per_arm)) +
      geom_hline(yintercept = 618, linetype = "dashed", colour = CEIL) +
      annotate("text", x = 0.42, y = 700, hjust = 1,
               label = "EURAMOS-1 randomised 618 per arm",
               colour = CEIL, size = 3.6) +
      geom_line(linewidth = 0.8, colour = "#2d2d6b") +
      geom_point(size = 3, colour = "#2d2d6b") +
      geom_text(aes(label = sprintf("%.0f", n_per_arm)), vjust = -0.8,
                size = 3.4) +
      scale_y_log10() +
      labs(x = expression(Delta~log[10]~"cell kill from an intensification"),
           y = "patients per arm for 80% power") + THEME
  })

  output$spread <- renderPlot({
    Kbar <- fin()$KNATS
    d <- do.call(rbind, lapply(c(0.02, 0.05, 0.10, 0.15, 0.25, 0.35, 0.50,
                                 0.65), function(cv) {
      data.frame(cv = cv,
                 gain = pop_cure(Kbar + 0.5 * log(10), cv) - pop_cure(Kbar, cv))
    }))
    ggplot(d, aes(cv, gain)) +
      geom_col(fill = "#8f8fc4", width = 0.04) +
      geom_vline(xintercept = input$cv, linetype = "dashed", colour = CEIL) +
      labs(x = "between-patient CV of achieved log-kill",
           y = "absolute survival gain from +0.5 log kill",
           subtitle = "Same drug, same half-log. Only the spread changes.") +
      THEME
  })

  output$tailplot <- renderPlot({
    d <- budget_targeting(fin()$KNATS, cv = input$cv)
    d$hi <- d$gain == max(d$gain)
    ggplot(d, aes(reorder(who, gain), gain, fill = hi)) +
      geom_col(width = 0.6, show.legend = FALSE) +
      scale_fill_manual(values = c("FALSE" = "#8f8fc4", "TRUE" = GOOD)) +
      geom_text(aes(label = sprintf("%+.3f", gain)), hjust = -0.15,
                size = 3.8) +
      coord_flip() + ylim(0, max(d$gain) * 1.25) +
      labs(x = NULL, y = "absolute gain in cure fraction") + THEME
  })

  output$popdist <- renderPlot({
    d <- marginal_by_decile(fin()$KNATS, cv = input$cv)
    d$mid <- (d$lo_log10 + d$hi_log10) / 2
    d$hi  <- d$marginal_value == max(d$marginal_value)
    ggplot(d, aes(factor(decile), marginal_value, fill = hi)) +
      geom_col(width = 0.75, show.legend = FALSE) +
      scale_fill_manual(values = c("FALSE" = "#8f8fc4", "TRUE" = GOOD)) +
      geom_text(aes(label = sprintf("%.2f-%.2f", lo_log10, hi_log10)),
                vjust = -0.5, size = 3, colour = "#555555") +
      labs(x = "decile of achieved log-kill (labels: log10 range)",
           y = "marginal value of +1 log per unit mass",
           subtitle = "Zero at both ends. The shoulder is where the budget belongs.") +
      THEME
  })
}

shinyApp(ui, server)
