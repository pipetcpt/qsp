##  Idiopathic Intracranial Hypertension — QSP explorer (Shiny)
##  ============================================================================
##  Companion dashboard for iih_mrgsolve_model.R.
##
##  DESIGN INTENT
##  -------------
##  Most PK/PD dashboards let you turn the dose up and watch the effect grow.
##  This one is built to show you where that stops working. The organising
##  quantity is the VENOUS FLOOR
##
##        ICP_floor = P_cv + G(ICP_floor)
##
##  the pressure below which no secretion-blocking drug can go at any dose.
##  Every pressure plot draws the floor alongside ICP, and the headroom
##  between them (DRUG_ROOM) is displayed as a first-class read-out. The
##  "Composition" tab exists because that headroom depends on a split —
##  venous gradient versus outflow resistance — that two decades of trials
##  did not measure: you can hold the opening pressure at exactly the value
##  written in a patient's notes and watch the achievable drug effect vary
##  more than two-fold.
##
##  Tabs
##    1  Patient & phenotype   — build the virtual patient, see the floor
##    2  Composition & bound   — the headline result: same ICP, different ceiling
##    3  Drug PK/PD            — exposure, CA inhibition, tolerability ceiling
##    4  Pressure & the loop   — ICP, P_sss, gradient, loop gain, amplification
##    5  Optic nerve & vision  — swelling vs axon loss, Frisen, PMD, OCT traps
##    6  Symptoms              — headache with sensitisation, tinnitus, diplopia
##    7  Scenario comparison   — the 14 prebuilt therapy scenarios side by side
##    8  Diagnostic studies    — CSF infusion test and lumbar puncture dynamics
##    9  Trial reproduction    — model vs IIHTT / IIH:WT / exenatide / stenting
##
##  Run:  shiny::runApp("iih_shiny_app.R")
##  Requires: shiny, mrgsolve, ggplot2, dplyr, tidyr  (+ iih_mrgsolve_model.R)
##  ============================================================================

library(shiny)
library(mrgsolve)
library(ggplot2)
library(dplyr)
library(tidyr)

source("iih_mrgsolve_model.R")

REMISSION_CM <- 25.0        # Friedman 2013 diagnostic threshold

## Published anchors, shown in the trial-reproduction tab.
ANCHORS <- tibble::tribble(
  ~study,                   ~endpoint,                  ~observed, ~unit,
  "IIHTT (PMID 24756514)",  "ICP, between-arm",             -59.9, "mmH2O",
  "IIHTT",                  "PMD, between-arm",              0.71, "dB",
  "IIHTT",                  "Frisen, between-arm",          -0.70, "grade",
  "IIHTT",                  "weight, between-arm",          -4.05, "kg",
  "IIHTT OCT (26198807)",   "RNFL, acetazolamide arm",     -175.0, "um",
  "IIH:WT (33900360)",      "ICP at 12 mo vs control",       -6.0, "cmH2O",
  "IIH:WT",                 "weight at 12 mo vs control",   -21.4, "kg",
  "VLCD (20610512)",        "ICP after -15.7 kg",            -8.0, "cmH2O",
  "Exenatide (36907221)",   "ICP at 2.5 h",                  -5.7, "cmH2O",
  "Stenting (29871989)",    "CSF-OP at 3 mo",               -16.8, "cmH2O",
  "Stenting manometry (29922401)", "sagittal sinus pressure", -8.1, "mmHg"
)

## ============================================================================
## UI
## ============================================================================
ui <- fluidPage(
  titlePanel("Idiopathic Intracranial Hypertension — QSP model explorer"),
  tags$p(style = "color:#555;margin-top:-8px;",
         tags$b("ICP = I_f x R_out + P_sss."),
         " Every drug in IIH acts on I_f. P_sss is an additive floor none of",
         " them touches. This app is organised around that bound."),
  sidebarLayout(
    sidebarPanel(
      width = 3,
      selectInput("pheno", "Phenotype (each pinned to its published baseline pressure)",
                  choices = c("IIHTT-like (mild)" = "iihtt",
                              "Pure outflow resistance" = "resistive",
                              "Stenosis-dominated" = "stenotic",
                              "Fulminant" = "fulminant",
                              "Healthy control" = "normal"),
                  selected = "iihtt"),
      hr(),
      h4("Patient"),
      sliderInput("W0", "Body weight (kg)", 60, 170, 107.7, 0.5),
      sliderInput("GMAX", "Max trans-stenotic gradient (mmHg)", 0, 30, 13, 0.5),
      sliderInput("PCV0", "Post-stenotic venous pressure (mmHg)", 3, 16, 8, 0.5),
      sliderInput("ROUT0", "CSF outflow resistance (mmHg/(mL/min))",
                  5, 80, 27.9, 0.1),
      helpText("Normal outflow resistance is 6-10. The three sliders above",
               "decide how much of the pressure is FLOOR and how much is",
               "resistive — the split that decides whether a drug can work."),
      hr(),
      h4("Therapy"),
      sliderInput("acz", "Acetazolamide (mg/day, 4 divided doses)",
                  0, 4000, 0, 250),
      sliderInput("tpm", "Topiramate (mg/day)", 0, 400, 0, 25),
      checkboxInput("glp", "Exenatide 10 ug SC bid", FALSE),
      sliderInput("diet", "Diet: weight set-point shift (kg)", -30, 0, 0, 1),
      checkboxInput("bari", "Bariatric surgery at day 30", FALSE),
      checkboxInput("stent", "Venous sinus stenting at day 14", FALSE),
      sliderInput("stent_eff", "  fraction of gradient abolished",
                  0, 1, 0.90, 0.05),
      checkboxInput("shunt", "CSF shunt at day 14", FALSE),
      checkboxInput("onsf", "Optic nerve sheath fenestration", FALSE),
      hr(),
      sliderInput("end", "Simulation horizon (days)", 30, 730, 365, 5),
      actionButton("go", "Simulate", class = "btn-primary")
    ),
    mainPanel(
      width = 9,
      tabsetPanel(
        id = "tabs",
        tabPanel("1 Patient & phenotype",
                 br(), fluidRow(
                   column(4, wellPanel(h4("Baseline"), tableOutput("t_base"))),
                   column(8, plotOutput("p_decomp", height = "330px"))),
                 htmlOutput("txt_floor")),
        tabPanel("2 Composition & bound",
                 br(),
                 htmlOutput("txt_bound"),
                 numericInput("fix_icp", "Hold the measured opening pressure at (cmH2O)",
                              37, 26, 60, 1),
                 plotOutput("p_ceiling", height = "360px"),
                 tableOutput("t_ceiling")),
        tabPanel("3 Drug PK/PD",
                 br(), plotOutput("p_pk", height = "300px"),
                 plotOutput("p_pd", height = "260px"),
                 htmlOutput("txt_pd")),
        tabPanel("4 Pressure & the loop",
                 br(), plotOutput("p_icp", height = "330px"),
                 plotOutput("p_loop", height = "280px"),
                 htmlOutput("txt_loop")),
        tabPanel("5 Optic nerve & vision",
                 br(), plotOutput("p_eye", height = "330px"),
                 plotOutput("p_oct", height = "280px"),
                 htmlOutput("txt_oct")),
        tabPanel("6 Symptoms",
                 br(), plotOutput("p_sym", height = "330px"),
                 htmlOutput("txt_sym")),
        tabPanel("7 Scenario comparison",
                 br(), checkboxGroupInput(
                   "scn", "Scenarios", inline = TRUE,
                   choices = c("01_natural_history", "02_IIHTT_placebo_diet",
                               "03_IIHTT_acetazolamide", "04_acetazolamide_4g",
                               "06_very_low_calorie_diet", "07_bariatric_surgery",
                               "08_exenatide", "09_venous_sinus_stenting",
                               "11_maximal_drug_vs_floor", "12_csf_shunt",
                               "14_acetazolamide_withdrawal"),
                   selected = c("01_natural_history", "03_IIHTT_acetazolamide",
                                "09_venous_sinus_stenting",
                                "11_maximal_drug_vs_floor")),
                 plotOutput("p_scn", height = "420px"),
                 tableOutput("t_scn")),
        tabPanel("8 Diagnostic studies",
                 br(), h4("CSF infusion test — what it actually measures"),
                 plotOutput("p_inf", height = "300px"),
                 tableOutput("t_inf"),
                 h4("Therapeutic lumbar puncture — and how fast CSF re-forms"),
                 plotOutput("p_lp", height = "260px"),
                 htmlOutput("txt_diag")),
        tabPanel("9 Trial reproduction",
                 br(), tableOutput("t_anchor"),
                 htmlOutput("txt_anchor"))
      )
    )
  )
)

## ============================================================================
## SERVER
## ============================================================================
server <- function(input, output, session) {

  ## keep the sliders in step with the selected phenotype
  observeEvent(input$pheno, {
    p <- IIH_phenotype(input$pheno)
    updateSliderInput(session, "W0", value = p$W0)
    updateSliderInput(session, "GMAX", value = p$GMAX)
    updateSliderInput(session, "PCV0", value = p$PCV0)
    updateSliderInput(session, "ROUT0", value = p$ROUT0)
  })

  pars <- reactive({
    p <- IIH_phenotype(input$pheno)
    p$W0 <- input$W0; p$GMAX <- input$GMAX
    p$PCV0 <- input$PCV0; p$ROUT0 <- input$ROUT0
    p$D_DIET <- input$diet
    p$ONSF <- as.numeric(input$onsf)
    if (input$bari)  { p$D_BARI <- -30; p$T_BARI <- 30 }
    if (input$shunt) { p$SHUNT_ON <- 1; p$T_SHUNT <- 14 }
    p
  })

  evts <- reactive({
    e <- NULL
    add <- function(a, b) if (is.null(a)) b else c(a, b)
    if (input$acz > 0) e <- add(e, acz_ev(input$acz, input$end))
    if (input$tpm > 0) e <- add(e, tpm_ev(input$tpm, input$end))
    if (input$glp)     e <- add(e, glp_ev(10, input$end))
    if (input$stent)   e <- add(e, stent_ev(input$stent_eff, 14))
    e
  })

  sim <- eventReactive(input$go, {
    m <- param(iih, pars())
    e <- evts()
    out <- if (is.null(e)) mrgsim(m, end = input$end, delta = 1)
           else mrgsim(m, events = e, end = input$end, delta = 1)
    as.data.frame(out)
  }, ignoreNULL = FALSE)

  base_row <- reactive(sim()[1, ])

  ## ---- tab 1 --------------------------------------------------------------
  output$t_base <- renderTable({
    b <- base_row()
    data.frame(
      quantity = c("ICP (cmH2O)", "sinus pressure P_sss (cmH2O)",
                   "trans-stenotic gradient (mmHg)", "loop gain gamma",
                   "amplification 1/(1-gamma)", "VENOUS FLOOR (cmH2O)",
                   "max fall any I_f drug can give (cmH2O)",
                   "remission reachable by drug?",
                   "papilloedema (Frisen)", "PMD (dB)"),
      value = c(sprintf("%.1f", b$ICP_cm), sprintf("%.1f", b$PSSS_cm),
                sprintf("%.1f", b$GRAD), sprintf("%.3f", b$GAMMA),
                sprintf("%.2f", b$AMPLIF), sprintf("%.1f", b$FLOOR_cm),
                sprintf("%.1f", b$DRUG_ROOM),
                ifelse(b$FLOOR_OK > 0.5, "yes", "NO — structural"),
                sprintf("%.2f", b$FRISEN), sprintf("%.2f", b$PMD)))
  }, colnames = FALSE)

  output$p_decomp <- renderPlot({
    b <- base_row()
    resistive <- b$ICP_cm - b$PSSS_cm
    d <- data.frame(
      component = factor(c("P_cv (abdominal/venous)",
                           "G(ICP) trans-stenotic gradient",
                           "I_f x R_out (resistive)"),
                         levels = c("I_f x R_out (resistive)",
                                    "G(ICP) trans-stenotic gradient",
                                    "P_cv (abdominal/venous)")),
      cm = c(b$PSSS_cm - b$GRAD * 1.35951, b$GRAD * 1.35951, resistive))
    ggplot(d, aes("patient", cm, fill = component)) +
      geom_col(width = 0.45) +
      geom_hline(yintercept = REMISSION_CM, linetype = 2, colour = "red3") +
      annotate("text", x = 1.35, y = REMISSION_CM + 1.2,
               label = "remission threshold 25 cmH2O", colour = "red3",
               size = 3.4, hjust = 0.5) +
      scale_fill_manual(values = c("I_f x R_out (resistive)" = "#8fd6b4",
                                   "G(ICP) trans-stenotic gradient" = "#ffb066",
                                   "P_cv (abdominal/venous)" = "#f4a0bd")) +
      labs(title = "What the opening pressure is made of",
           subtitle = paste("Only the green slice is available to a drug.",
                            "The two lower slices are the floor."),
           x = NULL, y = "cmH2O", fill = NULL) +
      coord_flip() + theme_minimal(base_size = 12)
  })

  output$txt_floor <- renderUI({
    b <- base_row()
    HTML(sprintf(
      "<p style='margin-top:8px'>The floor for this patient is <b>%.1f cmH2O</b>.
       A drug that abolished CSF secretion entirely would take the pressure from
       %.1f to %.1f cmH2O — a fall of %.1f cmH2O and no more. Remission
       (&le; 25 cmH2O) is therefore <b>%s</b> by any secretion-blocking drug at
       any dose.%s</p>",
      b$FLOOR_cm, b$ICP_cm, b$FLOOR_cm, b$DRUG_ROOM,
      ifelse(b$FLOOR_OK > 0.5, "reachable in principle",
             "UNREACHABLE — the failure is structural, not a dosing problem"),
      ifelse(b$FLOOR_OK > 0.5, "",
             " Only stenting (which moves the floor), weight loss (which lowers
              its base) or CSF diversion (which bypasses the equation) can help
              this patient.")))
  })

  ## ---- tab 2: the headline calculation ------------------------------------
  ceiling_tab <- reactive({
    icp_cm <- input$fix_icp
    grid <- seq(0, 26, by = 2)
    rout_for <- function(gmax) {
      f <- function(rout) {
        p <- modifyList(pars(), list(GMAX = gmax, ROUT0 = rout))
        as.numeric(mrgsim(param(iih, p), end = 0, delta = 1)$ICP_cm[1]) - icp_cm
      }
      tryCatch(stats::uniroot(f, c(0.2, 500))$root, error = function(e) NA_real_)
    }
    do.call(rbind, lapply(grid, function(g) {
      r <- rout_for(g)
      if (is.na(r)) return(NULL)
      p <- modifyList(pars(), list(GMAX = g, ROUT0 = r))
      s <- mrgsim(param(iih, p), end = 0, delta = 1)
      data.frame(gradient = as.numeric(s$GRAD[1]), R_out = r,
                 ICP = as.numeric(s$ICP_cm[1]),
                 floor = as.numeric(s$FLOOR_cm[1]),
                 room = as.numeric(s$DRUG_ROOM[1]),
                 reachable = as.numeric(s$FLOOR_OK[1]) > 0.5)
    }))
  })

  output$p_ceiling <- renderPlot({
    d <- ceiling_tab()
    req(nrow(d) > 1)
    ggplot(d, aes(gradient)) +
      geom_ribbon(aes(ymin = floor, ymax = ICP), fill = "#8fd6b4", alpha = 0.55) +
      geom_line(aes(y = ICP), linewidth = 1.1) +
      geom_line(aes(y = floor), colour = "#d2691e", linewidth = 1.2) +
      geom_hline(yintercept = REMISSION_CM, linetype = 2, colour = "red3") +
      annotate("text", x = max(d$gradient) * 0.62, y = REMISSION_CM + 1.1,
               label = "remission 25 cmH2O", colour = "red3", size = 3.4) +
      labs(title = paste0("Every point here is the SAME measured pressure (",
                          input$fix_icp, " cmH2O)"),
           subtitle = paste("black = measured ICP; orange = venous floor;",
                            "green band = all the room a drug has"),
           x = "trans-stenotic gradient (mmHg) — not measured in any IIH drug trial",
           y = "cmH2O") +
      theme_minimal(base_size = 12)
  })

  output$t_ceiling <- renderTable({
    d <- ceiling_tab()
    req(nrow(d) > 0)
    d %>%
      filter(gradient %% 4 < 2) %>%
      transmute(`gradient (mmHg)` = round(gradient, 1),
                `R_out` = round(R_out, 1),
                `floor (cmH2O)` = round(floor, 1),
                `max drug fall (cmH2O)` = round(room, 1),
                `remission reachable` = ifelse(reachable, "yes", "NO"))
  })

  output$txt_bound <- renderUI(HTML(
    "<p>This is the central result of the model. Hold the opening pressure
     fixed at a value you might read in a patient's notes, and vary only the
     <i>composition</i> of that pressure. The green band is the entire
     therapeutic room available to acetazolamide, topiramate, a GLP-1 agonist,
     a diuretic, or any combination of them &mdash; and it shrinks by more than
     two-fold across patients who are clinically indistinguishable without
     venography. Where the orange floor crosses the red line, remission is
     unreachable by drugs at any dose.</p>"))

  ## ---- tab 3 --------------------------------------------------------------
  output$p_pk <- renderPlot({
    d <- sim()
    d %>% select(time, `acetazolamide (mg/L)` = CACZ,
                 `serum HCO3 (mmol/L)` = HCO3) %>%
      pivot_longer(-time) %>%
      ggplot(aes(time, value)) + geom_line(linewidth = 0.9) +
      facet_wrap(~name, scales = "free_y") +
      labs(title = "Exposure and the tolerability signal that caps the dose",
           x = "day", y = NULL) + theme_minimal(base_size = 12)
  })

  output$p_pd <- renderPlot({
    d <- sim()
    d %>% select(time, `CA inhibition effect (EACZ)` = EACZ,
                 `GLP-1R effect (EGLP)` = EGLP,
                 `total I_f suppression` = EPS_TOT) %>%
      pivot_longer(-time) %>%
      ggplot(aes(time, value, colour = name)) + geom_line(linewidth = 0.9) +
      geom_hline(yintercept = 0.55, linetype = 3) +
      annotate("text", x = Inf, y = 0.57, hjust = 1.05, size = 3.3,
               label = "carbonic-anhydrase ceiling ~0.55") +
      labs(title = "Fractional suppression of CSF formation",
           x = "day", y = "fraction", colour = NULL) +
      theme_minimal(base_size = 12) + theme(legend.position = "bottom")
  })

  output$txt_pd <- renderUI(HTML(
    "<p>CSF flow does not fall until choroid-plexus carbonic anhydrase is
     almost completely inhibited, which is why the exposure&ndash;response here
     is steep (Hill 3) and why IIHTT titrated towards 4 g/day. The practical
     ceiling is set by tolerability rather than pharmacology: bicarbonate falls,
     paraesthesia and fatigue follow, adherence drops, and 48% of the IIHTT arm
     never reached the target dose.</p>"))

  ## ---- tab 4 --------------------------------------------------------------
  output$p_icp <- renderPlot({
    d <- sim()
    ggplot(d, aes(time)) +
      geom_ribbon(aes(ymin = FLOOR_cm, ymax = pmax(ICP_cm, FLOOR_cm)),
                  fill = "#8fd6b4", alpha = 0.4) +
      geom_line(aes(y = ICP_cm, colour = "ICP"), linewidth = 1.1) +
      geom_line(aes(y = PSSS_cm, colour = "sinus pressure P_sss"),
                linewidth = 0.9) +
      geom_line(aes(y = FLOOR_cm, colour = "venous floor"), linewidth = 1.1) +
      geom_hline(yintercept = REMISSION_CM, linetype = 2, colour = "red3") +
      scale_colour_manual(values = c("ICP" = "black",
                                     "sinus pressure P_sss" = "#3b7dd8",
                                     "venous floor" = "#d2691e")) +
      labs(title = "Pressure, its floor, and the room in between",
           x = "day", y = "cmH2O", colour = NULL) +
      theme_minimal(base_size = 12) + theme(legend.position = "bottom")
  })

  output$p_loop <- renderPlot({
    d <- sim()
    d %>% select(time, `trans-stenotic gradient (mmHg)` = GRAD,
                 `loop gain gamma` = GAMMA,
                 `amplification 1/(1-gamma)` = AMPLIF,
                 `drug headroom (cmH2O)` = DRUG_ROOM) %>%
      pivot_longer(-time) %>%
      ggplot(aes(time, value)) + geom_line(linewidth = 0.9) +
      facet_wrap(~name, scales = "free_y") +
      labs(x = "day", y = NULL,
           title = "The loop: gain, amplification, and what is left for a drug") +
      theme_minimal(base_size = 12)
  })

  output$txt_loop <- renderUI(HTML(
    "<p>Because the sinus is collapsible, ICP feeds back on its own outflow
     pressure. Any input &mdash; a drug, a kilogram, an infusion &mdash; is
     amplified by 1/(1-gamma). The same factor inflates the apparent outflow
     resistance measured by a CSF infusion test (tab 8), so a resistance
     measured in IIH is a property of the loop rather than of the tissue.</p>"))

  ## ---- tab 5 --------------------------------------------------------------
  output$p_eye <- renderPlot({
    d <- sim()
    d %>% select(time, `translaminar gradient (mmHg)` = TLPG_out,
                 `papilloedema swelling (um)` = EDEMA,
                 `Frisen grade` = FRISEN, `PMD (dB)` = PMD) %>%
      pivot_longer(-time) %>%
      ggplot(aes(time, value)) + geom_line(linewidth = 0.9) +
      facet_wrap(~name, scales = "free_y") +
      labs(title = "Pressure becomes vision through a threshold and an integral",
           x = "day", y = NULL) + theme_minimal(base_size = 12)
  })

  output$p_oct <- renderPlot({
    d <- sim()
    d %>% select(time, `OCT RNFL (um)` = RNFL, `OCT RGCL (um)` = RGCL,
                 `axon loss (%)` = AXLOSS,
                 `cumulative injurious exposure (mmHg.d)` = AUCTL) %>%
      pivot_longer(-time) %>%
      ggplot(aes(time, value)) + geom_line(linewidth = 0.9, colour = "#2f6fa8") +
      facet_wrap(~name, scales = "free_y") +
      labs(title = "The two OCT traps, and the quantity that actually predicts loss",
           x = "day", y = NULL) + theme_minimal(base_size = 12)
  })

  output$txt_oct <- renderUI(HTML(
    "<p>RNFL thickness is axons <i>plus</i> swelling <i>minus</i> atrophy, so a
     falling RNFL is either recovery or blindness and the number alone cannot
     tell you which. The IIHTT's own ganglion-cell numbers show the trap: the
     acetazolamide arm 'thinned' more than placebo (3.6 vs 2.1 um), which reads
     as more damage but is oedema contamination resolving faster in the
     better-treated arm. Permanent loss tracks the <i>integral</i> of
     translaminar gradient above threshold, not today's pressure &mdash; which
     is why time-to-treatment matters more than the last increment of dose.</p>"))

  ## ---- tab 6 --------------------------------------------------------------
  output$p_sym <- renderPlot({
    d <- sim()
    d %>% select(time, `HIT-6 headache impact` = HIT6,
                 `central sensitisation` = SENS,
                 `pulsatile tinnitus (0-10)` = TINN,
                 `ICP (cmH2O)` = ICP_cm) %>%
      pivot_longer(-time) %>%
      ggplot(aes(time, value)) + geom_line(linewidth = 0.9, colour = "#7b52ab") +
      facet_wrap(~name, scales = "free_y") +
      labs(x = "day", y = NULL,
           title = "Symptoms: one tracks the gradient, one outlives the pressure") +
      theme_minimal(base_size = 12)
  })

  output$txt_sym <- renderUI(HTML(
    "<p>Pulsatile tinnitus tracks the trans-stenotic gradient, which is why
     stenting abolishes it almost immediately. Headache does not track pressure
     nearly as well: once central sensitisation and analgesic overuse are
     established they persist after the pressure is controlled, which is the
     commonest source of apparent treatment failure in clinic. That term is
     structurally motivated and quantitatively unanchored &mdash; treat it as a
     hypothesis, not a prediction.</p>"))

  ## ---- tab 7 --------------------------------------------------------------
  scns <- reactive(IIH_simulate_scenarios(end = input$end))

  output$p_scn <- renderPlot({
    req(length(input$scn) > 0)
    S <- scns()
    d <- bind_rows(lapply(input$scn, function(n)
      as.data.frame(S[[n]]) %>% mutate(scenario = n)))
    d %>% select(time, scenario, `ICP (cmH2O)` = ICP_cm,
                 `venous floor (cmH2O)` = FLOOR_cm,
                 `PMD (dB)` = PMD, `Frisen` = FRISEN) %>%
      pivot_longer(c(-time, -scenario)) %>%
      ggplot(aes(time, value, colour = scenario)) +
      geom_line(linewidth = 0.9) +
      facet_wrap(~name, scales = "free_y") +
      labs(x = "day", y = NULL, colour = NULL) +
      theme_minimal(base_size = 12) + theme(legend.position = "bottom")
  })

  output$t_scn <- renderTable({
    req(length(input$scn) > 0)
    S <- scns()
    do.call(rbind, lapply(input$scn, function(n) {
      d <- as.data.frame(S[[n]]); e <- d[nrow(d), ]
      data.frame(scenario = n,
                 `ICP end (cmH2O)` = round(e$ICP_cm, 1),
                 `floor (cmH2O)` = round(e$FLOOR_cm, 1),
                 `headroom left` = round(e$DRUG_ROOM, 1),
                 `PMD (dB)` = round(e$PMD, 2),
                 `axon loss (%)` = round(e$AXLOSS, 1),
                 remission = ifelse(e$REMISSION > 0.5, "yes", "no"),
                 check.names = FALSE)
    }))
  })

  ## ---- tab 8 --------------------------------------------------------------
  inf_runs <- reactive(IIH_infusion_study(rate = 1.5, minutes = 60))

  output$p_inf <- renderPlot({
    R <- inf_runs()
    d <- bind_rows(lapply(names(R), function(n)
      as.data.frame(R[[n]]) %>% mutate(phenotype = n, minutes = time * 1440)))
    d %>% select(minutes, phenotype, ICP, PSSS) %>%
      pivot_longer(c(ICP, PSSS)) %>%
      ggplot(aes(minutes, value, colour = phenotype, linetype = name)) +
      geom_line(linewidth = 0.9) +
      labs(title = "CSF infusion at 1.5 mL/min, with the sinus trace alongside",
           x = "minutes", y = "mmHg", colour = NULL, linetype = NULL) +
      theme_minimal(base_size = 12) + theme(legend.position = "bottom")
  })

  output$t_inf <- renderTable({
    R <- inf_runs()
    do.call(rbind, lapply(names(R), function(n) {
      d <- as.data.frame(R[[n]])
      p <- IIH_phenotype(n)
      d0 <- d[1, ]; d1 <- d[nrow(d), ]
      app <- (d1$ICP - d0$ICP) / 1.5
      data.frame(phenotype = n,
                 `true R_out` = round(p$ROUT0, 1),
                 `apparent R_out` = round(app, 1),
                 `ratio` = round(app / p$ROUT0, 2),
                 `1/(1-gamma)` = round(d0$AMPLIF, 2),
                 `dP_sss/dICP` = round((d1$PSSS - d0$PSSS) /
                                         max(d1$ICP - d0$ICP, 1e-9), 2),
                 check.names = FALSE)
    }))
  })

  output$p_lp <- renderPlot({
    d <- as.data.frame(IIH_lumbar_puncture())
    ggplot(d, aes(time * 24, ICP_cm)) + geom_line(linewidth = 1) +
      geom_hline(yintercept = REMISSION_CM, linetype = 2, colour = "red3") +
      labs(title = "25 mL removed over 10 minutes",
           subtitle = "the relaxation time constant is ~10 minutes, which is why an LP is a test and not a treatment",
           x = "hours", y = "ICP (cmH2O)") + theme_minimal(base_size = 12)
  })

  output$txt_diag <- renderUI(HTML(
    "<p>The table above is the model's practical recommendation. A classical
     infusion test returns R_out/(1-gamma), not R_out, so in a stenosis-dominated
     patient it over-states the resistance by exactly the amplification factor
     &mdash; and any drug effect predicted from that resistance is over-stated
     with it. Adding a simultaneous sinus pressure trace gives dP_sss/dICP,
     which identifies gamma directly and separates drug potency from loop gain.
     Nothing about ICP measured alone can do that.</p>"))

  ## ---- tab 9 --------------------------------------------------------------
  output$t_anchor <- renderTable(ANCHORS)

  output$txt_anchor <- renderUI(HTML(
    "<p>The model is calibrated on the IIHTT <i>between-arm</i> differences,
     because that is the randomisation-protected estimand. It reproduces the
     between-arm pressure difference to about 5%, and it does <b>not</b>
     reproduce either arm's absolute fall: it is short by 49.7 mmH2O in the
     acetazolamide arm and 46.8 mmH2O in the placebo arm. Those two shortfalls
     being nearly equal is the point &mdash; an arm-independent offset is what a
     non-treatment effect looks like, and enrolment in IIHTT required a raised
     opening pressure, which is the classic setup for regression to the mean.
     Reading the trial at face value instead requires the drug to suppress CSF
     formation by 53%, at or above the ceiling for carbonic anhydrase
     inhibition.</p>
     <p>Two things the model gets wrong, kept in view rather than tuned away:
     it over-predicts the between-arm visual difference (1.31 vs 0.71 dB) and
     the weight-mediated share of the visual benefit (~16% vs the trial's 4%);
     and the exenatide anchor forces the GLP-1R secretory ceiling up to roughly
     acetazolamide's, with the 10 ug dose already saturated.</p>"))
}

shinyApp(ui, server)
