##  Bronchopulmonary Dysplasia (BPD) — QSP Shiny Dashboard
##  ============================================================================
##  Companion interface to bpd_mrgsolve_model.R.
##
##  The app is organised around ONE question: at 36 weeks postmenstrual age you
##  will be handed a number — the FiO2 and the mode of support this infant still
##  needs. That number is the definition of BPD, but it is NOT the deficit. The
##  deficit is the alveolar surface the developmental programme was scheduled to
##  build and did not, and this dashboard shows both side by side so the gap
##  between them is visible.
##
##  Run:
##      shiny::runApp("bpd_shiny_app.R")
##  Requires: shiny, mrgsolve, dplyr, tidyr, ggplot2, DT
##  ============================================================================

library(shiny)
library(mrgsolve)
library(dplyr)
library(tidyr)
library(ggplot2)
library(DT)

## Load the model and the scenario library from the sibling file. `local=TRUE`
## with a fresh environment keeps the app from clobbering a user's workspace.
.bpd_env <- new.env()
source("bpd_mrgsolve_model.R", local = .bpd_env)
mod           <- get("mod", envir = .bpd_env)
BPD_reference <- get("BPD_reference", envir = .bpd_env)
BPD_scenarios <- get("BPD_scenarios", envir = .bpd_env)
BPD_run       <- get("BPD_run", envir = .bpd_env)
BPD_endpoints <- get("BPD_endpoints", envir = .bpd_env)
BPD_window_remaining <- get("BPD_window_remaining", envir = .bpd_env)

THEME <- theme_minimal(base_size = 12) +
  theme(panel.grid.minor = element_blank(),
        plot.title = element_text(face = "bold", size = 12),
        legend.position = "bottom", legend.title = element_blank())

GRADE_LAB <- c("0 = no BPD", "1 = <=2 L/min", "2 = >2 L/min or NIPPV",
               "3 = invasive MV")

# ============================================================================
#  UI
# ============================================================================
ui <- fluidPage(
  titlePanel(
    div(
      h3("Bronchopulmonary Dysplasia — QSP Dashboard"),
      p(em(paste("Alveolarisation is a scheduled programme. Injury does not",
                 "destroy alveoli, it arrests the clock while the window",
                 "S_dev(PMA) runs out — so TIMING beats POTENCY.")))
    )
  ),

  sidebarLayout(
    sidebarPanel(
      width = 3,

      h4("1. The infant"),
      sliderInput("GA", "Gestational age at birth (weeks)",
                  min = 22, max = 34, value = 25, step = 0.5),
      numericInput("BW", "Birth weight (kg); 0 = GA default",
                   value = 0.70, min = 0, max = 3, step = 0.05),
      checkboxInput("IUGR", "IUGR / SGA (anti-angiogenic intrauterine milieu)",
                    FALSE),
      checkboxInput("CHORIO", "Histologic chorioamnionitis (HIT #1)", TRUE),
      checkboxInput("ANTESTER", "Complete antenatal corticosteroid course",
                    TRUE),
      checkboxInput("UREA0", "Ureaplasma airway colonisation at birth", FALSE),
      sliderInput("GENRISK", "Genetic risk index", 0, 1, 0, step = 0.1),

      hr(),
      h4("2. Respiratory strategy (the exposure term)"),
      checkboxInput("LISA", "LISA / nCPAP instead of delivery-room intubation",
                    FALSE),
      checkboxInput("VTV", "Volume-targeted ventilation", FALSE),
      checkboxInput("PHC", "Permissive hypercapnia", FALSE),
      radioButtons("SPO2HI", "SpO2 target band",
                   choices = c("91-95% (high)" = 1, "85-89% (low)" = 0),
                   selected = 1, inline = TRUE),
      sliderInput("PDA", "Days with a haemodynamically significant PDA",
                  min = 0, max = 30, value = c(1, 14)),
      sliderInput("SEP_T0", "Day of a late-onset sepsis episode (0 = none)",
                  min = 0, max = 56, value = 0, step = 1),

      hr(),
      h4("3. Drugs (the transduction and programme terms)"),
      checkboxInput("CAF", "Caffeine citrate", TRUE),
      conditionalPanel(
        "input.CAF",
        sliderInput("CAF_T0", "Caffeine start day", 1, 28, 1, step = 1),
        sliderInput("CAF_MD", "Caffeine citrate maintenance (mg/kg/day)",
                    2.5, 20, 5, step = 2.5)),
      selectInput("STEROID", "Corticosteroid schedule",
                  choices = c("none",
                              "PREMILOC hydrocortisone (d1-10, 8.5 mg/kg)",
                              "DART dexamethasone (10-day taper, 0.89 mg/kg)",
                              "early high-dose dexamethasone (0.5 mg/kg/d)",
                              "intratracheal budesonide + surfactant"),
                  selected = "none"),
      conditionalPanel(
        "input.STEROID != 'none'",
        sliderInput("STER_T0", "Steroid start day", 1, 35, 1, step = 1)),
      checkboxInput("VITA", "Vitamin A 5000 IU IM 3x/week x 4 weeks (Tyson)",
                    FALSE),
      checkboxInput("AZI", "Azithromycin 10 mg/kg/day x 7 days", FALSE),
      checkboxInput("FURO", "Furosemide (from day 10)", FALSE),
      checkboxInput("SIL", "Sildenafil 3 mg/kg/day (from day 28)", FALSE),
      sliderInput("INO", "Inhaled NO (ppm)", 0, 20, 0, step = 5),
      checkboxInput("NUTR", "Aggressive early protein/energy", FALSE),
      checkboxInput("MSC", "Intratracheal MSC (day 7)", FALSE),
      checkboxInput("RHIGF", "rhIGF-1 / rhIGFBP-3", FALSE),

      hr(),
      sliderInput("END_PMA", "Simulate to PMA (weeks)", 36, 60, 44, step = 2),
      helpText(HTML(paste(
        "<b>The 36-week vertical line</b> on every plot is the definitional",
        "clock stop. Everything to the right of it is the part of the story",
        "the trial endpoint cannot see.")))
    ),

    mainPanel(
      width = 9,
      tabsetPanel(
        id = "tabs", type = "tabs",

        # ---- 1 -------------------------------------------------------------
        tabPanel(
          "1. Endpoint card",
          br(),
          fluidRow(
            column(4, wellPanel(h4("At 36 weeks PMA"),
                                tableOutput("card36"))),
            column(4, wellPanel(h4("The developmental deficit"),
                                tableOutput("cardDeficit"))),
            column(4, wellPanel(h4("The price paid"),
                                tableOutput("cardCost")))
          ),
          hr(),
          h4("Grade trajectory — the support this infant needs, day by day"),
          plotOutput("gradePlot", height = "230px"),
          h4("What the chart records (FiO2) vs what was lost (LOSTW)"),
          plotOutput("fio2VsLostw", height = "300px"),
          helpText(paste("These two curves are the whole argument. FiO2 has a",
                         "reserve zone: it does not move until the surface",
                         "deficit is substantial, and it can be improved by",
                         "drugs that change no structure at all. LOSTW has no",
                         "reserve and no way back."))
        ),

        # ---- 2 -------------------------------------------------------------
        tabPanel(
          "2. The window",
          br(),
          h4("S_dev(PMA): the septation drive, and when it is spent"),
          plotOutput("windowPlot", height = "300px"),
          h4("Alveolar surface built vs the gestation-matched ideal"),
          plotOutput("alvPlot", height = "300px"),
          h4("How much of the programme is still ahead at birth?"),
          DT::dataTableOutput("windowTable"),
          helpText(paste("A 24-week infant must complete almost the whole",
                         "alveolarisation programme ex utero, and a large part",
                         "of it before the 36-week clock even stops. That",
                         "geometry, not a fitted risk coefficient, is where",
                         "the gestational gradient of BPD comes from."))
        ),

        # ---- 3 -------------------------------------------------------------
        tabPanel(
          "3. Drug PK",
          br(),
          h4("Concentrations (real neonatal PK, mg/L unless stated)"),
          plotOutput("pkPlot", height = "420px"),
          h4("Glucocorticoid receptor occupancy — one receptor, three ligands"),
          plotOutput("grPlot", height = "260px"),
          helpText(paste("Dexamethasone, hydrocortisone and budesonide are",
                         "converted to a single dexamethasone-equivalent",
                         "concentration before the occupancy equation, so the",
                         "comparison between schedules is a comparison of",
                         "exposure and timing rather than of units. The LUNG",
                         "and SYSTEMIC curves separate only for intratracheal",
                         "budesonide — which is the point of giving it that",
                         "way."))
        ),

        # ---- 4 -------------------------------------------------------------
        tabPanel(
          "4. The growth gate G",
          br(),
          h4("G and its five components"),
          plotOutput("gatePlot", height = "340px"),
          h4("What is holding the gate shut?"),
          plotOutput("gateBlockPlot", height = "280px"),
          helpText(paste("G is a weighted geometric mean, so it is dominated",
                         "by whichever component is lowest — a therapy that",
                         "doubles an already-adequate component buys almost",
                         "nothing. Use this tab to find the binding constraint",
                         "before choosing a drug on the left."))
        ),

        # ---- 5 -------------------------------------------------------------
        tabPanel(
          "5. Injury and inflammation",
          br(),
          h4("Where the injury drive comes from, decomposed"),
          plotOutput("injPlot", height = "320px"),
          h4("The inflammatory cascade"),
          plotOutput("inflamPlot", height = "300px"),
          helpText(paste("The escalation loop is visible here: a rising",
                         "support requirement raises the mechanical and",
                         "oxidative terms, which raise NF-kB and IL-1beta,",
                         "which close the gate, which worsens the lung, which",
                         "raises the support requirement."))
        ),

        # ---- 6 -------------------------------------------------------------
        tabPanel(
          "6. Support and physiology",
          br(),
          h4("Respiratory support"),
          plotOutput("supPlot", height = "300px"),
          h4("Lung water, septal thickening, airway smooth muscle"),
          plotOutput("structPlot", height = "280px"),
          h4("Weight — which sets every mg/kg dose and every clearance"),
          plotOutput("wtPlot", height = "240px")
        ),

        # ---- 7 -------------------------------------------------------------
        tabPanel(
          "7. BPD-PH and the right ventricle",
          br(),
          h4("Pulmonary vascular resistance and its three terms"),
          plotOutput("phPlot", height = "320px"),
          h4("Microvascular density and RV hypertrophy"),
          plotOutput("rvPlot", height = "280px"),
          helpText(HTML(paste(
            "PVR is modelled as <b>(area term) x (tone term) x (remodelling",
            "term)</b>. A PDE5 inhibitor or inhaled NO enters only the tone",
            "term, which is why they lower the number without changing the",
            "disease. The only durable way down is to build the capillary",
            "bed — and that has to happen while the window is open.")))
        ),

        # ---- 8 -------------------------------------------------------------
        tabPanel(
          "8. Scenario comparison",
          br(),
          checkboxGroupInput(
            "scen", "Prebuilt scenarios to compare",
            choices = names(BPD_scenarios),
            selected = names(BPD_scenarios)[c(1, 2, 4, 7, 8, 16)],
            inline = FALSE),
          actionButton("runScen", "Run selected scenarios",
                       class = "btn-primary"),
          br(), br(),
          DT::dataTableOutput("scenTable"),
          br(),
          plotOutput("scenBar", height = "340px")
        ),

        # ---- 9 -------------------------------------------------------------
        tabPanel(
          "9. Timing vs potency",
          br(),
          p(paste("An idealised transduction-class gate opener, given at",
                  "different start days and two different potencies. The",
                  "question is how late a perfect drug can arrive and still",
                  "beat a mediocre one given on day 1.")),
          fluidRow(
            column(4, sliderInput("eff_weak", "Weak effect (fraction)",
                                  0.1, 0.9, 0.4, step = 0.1)),
            column(4, sliderInput("eff_strong", "Strong effect (fraction)",
                                  0.5, 1.0, 1.0, step = 0.1)),
            column(4, br(), actionButton("runTiming", "Run the sweep",
                                         class = "btn-primary"))
          ),
          plotOutput("timingPlot", height = "360px"),
          DT::dataTableOutput("timingTable")
        ),

        # ---- 10 ------------------------------------------------------------
        tabPanel(
          "10. Trade-offs",
          br(),
          h4("Net utility = survival without BPD  -  w x neurodevelopmental cost"),
          plotOutput("utilPlot", height = "360px"),
          sliderInput("wndi", "Weight on neurodevelopmental cost (w)",
                      0, 2, 0.45, step = 0.05, width = "60%"),
          helpText(HTML(paste(
            "<b>This weight is a value judgement, not a measurement.</b> The",
            "model can tell you how much lung a steroid buys and how much",
            "dexamethasone-equivalent it costs; it cannot tell you the",
            "exchange rate between a ventilator day and a point of",
            "developmental quotient. Move the slider and watch the ranking",
            "change — that instability IS the honest answer, and it is why",
            "early high-dose dexamethasone was abandoned despite reducing",
            "BPD."))),
          hr(),
          h4("The three classes of therapy, and what each can and cannot do"),
          tableOutput("classTable")
        )
      )
    )
  )
)

# ============================================================================
#  SERVER
# ============================================================================
server <- function(input, output, session) {

  ## ---- assemble the parameter set from the sidebar -------------------------
  pars <- reactive({
    p <- BPD_reference
    p$GA       <- input$GA
    p$BW       <- input$BW
    p$IUGR     <- as.numeric(input$IUGR)
    p$CHORIO   <- as.numeric(input$CHORIO)
    p$ANTESTER <- as.numeric(input$ANTESTER)
    p$UREA0    <- as.numeric(input$UREA0)
    p$GENRISK  <- input$GENRISK
    p$LISA     <- as.numeric(input$LISA)
    p$VTV      <- as.numeric(input$VTV)
    p$PHC      <- as.numeric(input$PHC)
    p$HYPOCAP  <- as.numeric(!input$VTV)
    p$SPO2HI   <- as.numeric(input$SPO2HI)

    p$PDA_T0 <- if (input$PDA[2] > input$PDA[1]) max(input$PDA[1], 0.5) else 999
    p$PDA_T1 <- if (input$PDA[2] > input$PDA[1]) input$PDA[2] else 0
    p$SEP_T0 <- if (input$SEP_T0 > 0) input$SEP_T0 else 999

    if (input$CAF) { p$CAF_T0 <- input$CAF_T0; p$CAF_MD <- input$CAF_MD }

    p$HC_PREM <- 0; p$DEX_DART <- 0
    if (input$STEROID == "PREMILOC hydrocortisone (d1-10, 8.5 mg/kg)") {
      p$HC_PREM <- 1; p$HC_T0 <- input$STER_T0
    } else if (input$STEROID == "DART dexamethasone (10-day taper, 0.89 mg/kg)") {
      p$DEX_DART <- 1; p$DEX_T0 <- input$STER_T0
    } else if (input$STEROID == "early high-dose dexamethasone (0.5 mg/kg/d)") {
      p$DEX_RATE <- 0.5; p$DEX_T0 <- input$STER_T0
      p$DEX_T1 <- input$STER_T0 + 6
    } else if (input$STEROID == "intratracheal budesonide + surfactant") {
      p$BUD_DOSE <- 0.25; p$BUD_T0 <- input$STER_T0; p$BUD_N <- 4
    }

    if (input$VITA)  { p$VITA_T0 <- 1 }
    if (input$AZI)   { p$AZI_RATE <- 10; p$AZI_T0 <- 1; p$AZI_T1 <- 7 }
    if (input$FURO)  { p$FURO <- 1; p$FURO_T0 <- 10 }
    if (input$SIL)   { p$SIL_RATE <- 3; p$SIL_T0 <- 28; p$SIL_T1 <- 200 }
    if (input$INO > 0) p$INO_PPM <- input$INO
    if (input$NUTR)  p$NUTR <- 1
    if (input$MSC)   { p$MSC <- 1; p$MSC_T0 <- 7 }
    if (input$RHIGF) p$RHIGF <- 1
    p
  })

  sim <- reactive({
    withProgress(message = "Integrating 37 ODEs...", value = 0.5, {
      BPD_run(pars(), end_pma = input$END_PMA, dt = 0.25)
    })
  })

  ep36 <- reactive(BPD_endpoints(sim(), ga = input$GA))

  vline36 <- function() {
    geom_vline(xintercept = 36, linetype = 3, colour = "grey35")
  }

  ## ---- tab 1: the endpoint card -------------------------------------------
  output$card36 <- renderTable({
    e <- ep36()
    data.frame(
      metric = c("BPD grade", "FiO2 required", "on invasive MV (fraction)",
                 "cumulative ventilator days", "cumulative O2 days",
                 "P(survive to 36 wk)"),
      value = c(GRADE_LAB[e$grade + 1], sprintf("%.2f", e$FiO2_36),
                sprintf("%.2f", e$MV_frac_36), sprintf("%.0f", e$vent_days),
                sprintf("%.0f", e$O2_days), sprintf("%.3f", e$P_survive)))
  }, colnames = FALSE)

  output$cardDeficit <- renderTable({
    e <- ep36()
    data.frame(
      metric = c("alveolar surface (36-wk normal = 1.0)",
                 "% of gestation-matched ideal",
                 "LOSTW (surface never built)",
                 "microvascular density",
                 "PVR (healthy = 1.0)"),
      value = c(sprintf("%.3f", e$ALV36), sprintf("%.1f%%", e$ALV_pct_ideal),
                sprintf("%.3f", e$LOSTW), sprintf("%.3f", e$CAP36),
                sprintf("%.2f", e$PVR36)))
  }, colnames = FALSE)

  output$cardCost <- renderTable({
    e <- ep36()
    data.frame(
      metric = c("cumulative dexamethasone-equivalent (mg/kg)",
                 "neurodevelopmental cost (au)",
                 "weight at 36 wk (kg)",
                 "RV hypertrophy index",
                 "net utility"),
      value = c(sprintf("%.3f", e$DEXEQ_mgkg), sprintf("%.3f", e$NDI_cost),
                sprintf("%.2f", e$weight36), sprintf("%.2f", e$RVH36),
                sprintf("%.3f", e$net_utility)))
  }, colnames = FALSE)

  output$gradePlot <- renderPlot({
    s <- sim()
    ggplot(s, aes(PMAout, GRADE)) +
      geom_step(linewidth = 1.1, colour = "#b5651d") +
      scale_y_continuous(breaks = 0:3, limits = c(0, 3),
                         labels = c("none", "1", "2", "3")) +
      vline36() +
      labs(x = "PMA (weeks)", y = "support grade",
           title = "Mode of support (the 2019 NICHD grading, evaluated continuously)") +
      THEME
  })

  output$fio2VsLostw <- renderPlot({
    s <- sim()
    d <- s %>%
      transmute(PMAout,
                `FiO2 (what the chart records)` = FIO2,
                `LOSTW (what was lost)` = LOSTW) %>%
      pivot_longer(-PMAout)
    ggplot(d, aes(PMAout, value, colour = name)) +
      geom_line(linewidth = 1.1) + vline36() +
      facet_wrap(~name, scales = "free_y") +
      scale_colour_manual(values = c("#2c7fb8", "#d95f0e")) +
      labs(x = "PMA (weeks)", y = NULL) + THEME +
      theme(legend.position = "none")
  })

  ## ---- tab 2: the window ---------------------------------------------------
  output$windowPlot <- renderPlot({
    s <- sim()
    ggplot(s, aes(PMAout, SDEVo)) +
      geom_area(fill = "#a8e6a8", alpha = 0.6) +
      geom_line(linewidth = 1.1, colour = "#2b7a2b") +
      geom_vline(xintercept = input$GA, linetype = 2, colour = "#b5651d") +
      annotate("text", x = input$GA, y = 0.95, hjust = -0.05, size = 3.4,
               label = "birth", colour = "#b5651d") +
      vline36() +
      annotate("text", x = 36, y = 0.15, hjust = -0.05, size = 3.4,
               label = "36 wk: the clock stops", colour = "grey30") +
      labs(x = "PMA (weeks)", y = "S_dev",
           title = "The septation window is at its peak while the infant is in the NICU") +
      THEME
  })

  output$alvPlot <- renderPlot({
    s <- sim()
    d <- s %>%
      transmute(PMAout, built = ALV, ideal = ALV_ideal,
                `capillary bed` = CAPD) %>%
      pivot_longer(-PMAout)
    ggplot(d, aes(PMAout, value, colour = name, linetype = name)) +
      geom_line(linewidth = 1.1) + vline36() +
      scale_colour_manual(values = c("#7570b3", "#1b7837", "grey35")) +
      scale_linetype_manual(values = c(1, 1, 2)) +
      labs(x = "PMA (weeks)", y = "normalised to healthy 36-wk lung = 1.0",
           title = "The gap between the two solid lines is the disease") +
      THEME
  })

  output$windowTable <- DT::renderDataTable({
    DT::datatable(BPD_window_remaining(), options = list(dom = "t"),
                  rownames = FALSE)
  })

  ## ---- tab 3: PK -----------------------------------------------------------
  output$pkPlot <- renderPlot({
    s <- sim()
    d <- s %>%
      transmute(PMAout,
                `caffeine (mg/L)` = Ccaf_o,
                `dexamethasone (mg/L)` = Cdex_o,
                `hydrocortisone (mg/L)` = Chc_o,
                `budesonide, lung (mg/L)` = Cbud_o,
                `plasma retinol (umol/L)` = Cret_o,
                `sildenafil (mg/L)` = Csil_o,
                `azithromycin (mg/L)` = Cazi_o) %>%
      pivot_longer(-PMAout) %>%
      group_by(name) %>% filter(max(value) > 1e-8) %>% ungroup()
    if (!nrow(d)) {
      return(ggplot() + annotate("text", 0, 0, label = "No drug on board") +
               theme_void())
    }
    ggplot(d, aes(PMAout, value)) +
      geom_line(linewidth = 1.05, colour = "#2c7fb8") +
      facet_wrap(~name, scales = "free_y") + vline36() +
      labs(x = "PMA (weeks)", y = NULL) + THEME
  })

  output$grPlot <- renderPlot({
    s <- sim()
    d <- s %>% transmute(PMAout, lung = GRo_l, systemic = GRo_s) %>%
      pivot_longer(-PMAout)
    ggplot(d, aes(PMAout, value, colour = name)) +
      geom_line(linewidth = 1.1) + vline36() +
      scale_colour_manual(values = c("#1b7837", "#d95f0e")) +
      ylim(0, 1) +
      labs(x = "PMA (weeks)", y = "GR occupancy",
           title = "Lung occupancy buys the benefit; systemic occupancy pays the bill") +
      THEME
  })

  ## ---- tab 4: the gate -----------------------------------------------------
  output$gatePlot <- renderPlot({
    s <- sim()
    d <- s %>%
      transmute(PMAout, VEGF, `NO / cGMP` = NOB, `IGF-1` = IGF1,
                retinoate = RASIG, progenitors = EPC) %>%
      pivot_longer(-PMAout)
    ggplot(d, aes(PMAout, value, colour = name)) +
      geom_line(linewidth = 1.05) + vline36() + ylim(0, NA) +
      labs(x = "PMA (weeks)", y = "normalised (1.0 = healthy)",
           title = "G is a geometric mean — the lowest component is the binding constraint") +
      THEME
  })

  output$gateBlockPlot <- renderPlot({
    s <- sim()
    d <- s %>%
      transmute(PMAout,
                `IL-1beta` = IL1B, `TGF-beta1` = TGFB, ROS = ROS,
                `neutrophil elastase` = NEUT) %>%
      pivot_longer(-PMAout)
    ggplot(d, aes(PMAout, value, colour = name)) +
      geom_line(linewidth = 1.05) + vline36() +
      labs(x = "PMA (weeks)", y = "activity (au)",
           title = "The brakes on the gate") + THEME
  })

  ## ---- tab 5: injury -------------------------------------------------------
  output$injPlot <- renderPlot({
    s <- sim()
    d <- s %>%
      transmute(PMAout,
                `oxygenation demand index` = Dmdo,
                `hypoxaemia burden` = hbo,
                `FiO2` = FIO2,
                `invasive MV fraction` = MVF) %>%
      pivot_longer(-PMAout)
    ggplot(d, aes(PMAout, value, colour = name)) +
      geom_line(linewidth = 1.05) + facet_wrap(~name, scales = "free_y") +
      vline36() + labs(x = "PMA (weeks)", y = NULL,
                       title = "The exposure terms") + THEME +
      theme(legend.position = "none")
  })

  output$inflamPlot <- renderPlot({
    s <- sim()
    d <- s %>% transmute(PMAout, `NF-kB` = NFKB, `IL-1beta` = IL1B,
                         `neutrophil/elastase` = NEUT, ROS = ROS,
                         `TGF-beta1` = TGFB) %>%
      pivot_longer(-PMAout)
    ggplot(d, aes(PMAout, value, colour = name)) +
      geom_line(linewidth = 1.05) + vline36() +
      labs(x = "PMA (weeks)", y = "activity (au)") + THEME
  })

  ## ---- tab 6: support ------------------------------------------------------
  output$supPlot <- renderPlot({
    s <- sim()
    d <- s %>% transmute(PMAout, `support intensity` = SUP,
                         `invasive MV fraction` = MVF, FiO2 = FIO2) %>%
      pivot_longer(-PMAout)
    ggplot(d, aes(PMAout, value, colour = name)) +
      geom_line(linewidth = 1.1) + vline36() + ylim(0, 1) +
      labs(x = "PMA (weeks)", y = NULL) + THEME
  })

  output$structPlot <- renderPlot({
    s <- sim()
    d <- s %>% transmute(PMAout, `lung water` = LW,
                         `septal thickening` = SEPT,
                         `airway smooth muscle` = ASM) %>%
      pivot_longer(-PMAout)
    ggplot(d, aes(PMAout, value, colour = name)) +
      geom_line(linewidth = 1.05) + vline36() +
      labs(x = "PMA (weeks)", y = "au") + THEME
  })

  output$wtPlot <- renderPlot({
    ggplot(sim(), aes(PMAout, WT)) +
      geom_line(linewidth = 1.1, colour = "#7570b3") + vline36() +
      labs(x = "PMA (weeks)", y = "weight (kg)",
           title = "Weight sets the mg/kg dose AND scales every clearance") +
      THEME
  })

  ## ---- tab 7: PH -----------------------------------------------------------
  output$phPlot <- renderPlot({
    s <- sim()
    d <- s %>% transmute(PMAout, PVR = PVR,
                         `remodelling (VREM)` = VREM,
                         `hypoxaemia burden` = hbo) %>%
      pivot_longer(-PMAout)
    ggplot(d, aes(PMAout, value, colour = name)) +
      geom_line(linewidth = 1.05) + facet_wrap(~name, scales = "free_y") +
      vline36() + labs(x = "PMA (weeks)", y = NULL) + THEME +
      theme(legend.position = "none")
  })

  output$rvPlot <- renderPlot({
    s <- sim()
    d <- s %>% transmute(PMAout, `microvascular density` = CAPD,
                         `RV hypertrophy` = RVHY) %>%
      pivot_longer(-PMAout)
    ggplot(d, aes(PMAout, value, colour = name)) +
      geom_line(linewidth = 1.1) + vline36() +
      labs(x = "PMA (weeks)", y = NULL) + THEME
  })

  ## ---- tab 8: scenarios ----------------------------------------------------
  scenRes <- eventReactive(input$runScen, {
    req(length(input$scen) > 0)
    withProgress(message = "Running scenarios...", value = 0, {
      n <- length(input$scen)
      out <- lapply(seq_along(input$scen), function(i) {
        incProgress(1 / n, detail = input$scen[i])
        p <- BPD_scenarios[[input$scen[i]]]
        s <- BPD_run(p, end_pma = max(input$END_PMA, 36.5))
        e <- BPD_endpoints(s, ga = if (!is.null(p$GA)) p$GA else 25)
        bind_cols(tibble(scenario = input$scen[i]), e)
      })
      bind_rows(out)
    })
  })

  output$scenTable <- DT::renderDataTable({
    DT::datatable(scenRes(), options = list(pageLength = 20, scrollX = TRUE),
                  rownames = FALSE)
  })

  output$scenBar <- renderPlot({
    d <- scenRes() %>%
      transmute(scenario,
                `alveolar surface at 36 wk` = ALV36,
                `LOSTW (never built)` = LOSTW,
                `ventilator days / 100` = vent_days / 100,
                `NDI cost` = NDI_cost) %>%
      pivot_longer(-scenario)
    ggplot(d, aes(reorder(scenario, value), value, fill = name)) +
      geom_col(position = "dodge") + coord_flip() +
      labs(x = NULL, y = NULL,
           title = "Scenario comparison at 36 weeks PMA") + THEME
  })

  ## ---- tab 9: timing -------------------------------------------------------
  timingRes <- eventReactive(input$runTiming, {
    days <- c(1, 3, 5, 7, 10, 14, 21, 28)
    base <- pars()
    withProgress(message = "Sweeping start days...", value = 0, {
      rows <- lapply(days, function(d) {
        incProgress(1 / length(days), detail = paste("day", d))
        w <- BPD_run(modifyList(base, list(IL1RA = 1, IL1RA_E = input$eff_weak,
                                           IL1RA_T0 = d)), end_pma = 36.05)
        st <- BPD_run(modifyList(base, list(IL1RA = 1,
                                            IL1RA_E = input$eff_strong,
                                            IL1RA_T0 = d)), end_pma = 36.05)
        data.frame(start_day = d,
                   weak = tail(w$ALV, 1),
                   strong = tail(st$ALV, 1))
      })
      bind_rows(rows)
    })
  })

  output$timingPlot <- renderPlot({
    d <- timingRes()
    bar <- d$weak[d$start_day == 1]
    dl <- pivot_longer(d, -start_day, names_to = "potency",
                       values_to = "ALV36")
    ggplot(dl, aes(start_day, ALV36, colour = potency)) +
      geom_line(linewidth = 1.1) + geom_point(size = 2.2) +
      geom_hline(yintercept = bar, linetype = 2, colour = "#d95f0e") +
      annotate("text", x = max(d$start_day), y = bar, vjust = -0.6, hjust = 1,
               size = 3.6, colour = "#d95f0e",
               label = "the bar to beat: weak drug, day 1") +
      scale_colour_manual(values = c("#2c7fb8", "#1b7837")) +
      labs(x = "start day of therapy",
           y = "alveolar surface at 36 weeks PMA",
           title = paste("Where the green line crosses below the dashed line,",
                         "a perfect late drug has lost to a mediocre early one")) +
      THEME
  })

  output$timingTable <- DT::renderDataTable({
    d <- timingRes()
    bar <- d$weak[d$start_day == 1]
    d$strong_still_wins <- d$strong > bar
    DT::datatable(d, options = list(dom = "t"), rownames = FALSE) %>%
      DT::formatRound(c("weak", "strong"), 4)
  })

  ## ---- tab 10: trade-offs --------------------------------------------------
  output$utilPlot <- renderPlot({
    arms <- c("S00 reference (no adjunct)",
              "S01 caffeine (CAP dosing)",
              "S06 PREMILOC early low-dose hydrocortisone (d1-10)",
              "S07 DART late low-dose dexamethasone (d14-23)",
              "S08 early HIGH-dose dexamethasone (0.5 mg/kg/d, d1-7)",
              "S09 intratracheal budesonide + surfactant",
              "S15 EARLY bundle (all of it, from day 1)")
    arms <- intersect(arms, names(BPD_scenarios))
    withProgress(message = "Running the trade-off panel...", value = 0, {
      d <- bind_rows(lapply(arms, function(a) {
        incProgress(1 / length(arms), detail = a)
        p <- BPD_scenarios[[a]]
        e <- BPD_endpoints(BPD_run(p, end_pma = 36.05), ga = 25)
        data.frame(arm = a, benefit = e$P_survive * (e$grade < 1),
                   ndi = e$NDI_cost, alv = e$ALV36)
      }))
    })
    d$net <- d$benefit - input$wndi * d$ndi
    d$label <- sub("^S[0-9]+ ", "", d$arm)
    ggplot(d, aes(reorder(label, net), net, fill = net > 0)) +
      geom_col() + coord_flip() +
      geom_hline(yintercept = 0, colour = "grey25") +
      scale_fill_manual(values = c("TRUE" = "#1b7837", "FALSE" = "#b2182b"),
                        guide = "none") +
      labs(x = NULL, y = sprintf("net utility at w = %.2f", input$wndi),
           title = "Move the weight and the ranking moves. That instability is the answer.") +
      THEME
  })

  output$classTable <- renderTable({
    data.frame(
      class = c("EXPOSURE reducers", "TRANSDUCTION suppressors",
                "PROGRAMME supporters", "none of the three"),
      `acts on` = c("the injury input (vili, oxtox, PDA, Ureaplasma)",
                    "injury -> gate transduction (NF-kB, IL-1beta, TGF-beta1)",
                    "the gate G itself (VEGF, NO, IGF-1, retinoate)",
                    "a downstream number only"),
      examples = c(paste("caffeine (drive arm), LISA/nCPAP,",
                         "volume-targeted ventilation, permissive",
                         "hypercapnia, PDA closure, SpO2 band"),
                   paste("hydrocortisone, dexamethasone, intratracheal",
                         "budesonide, azithromycin if Ureaplasma+, IL-1Ra"),
                   paste("vitamin A, rhIGF-1/BP3, MSC,",
                         "the direct A2A arm of caffeine"),
                   "furosemide, bronchodilator, sildenafil, inhaled NO"),
      `ceiling set by` = c("the ventilator/oxygen-mediated share of injury",
                           paste("the steroid-responsive fraction, minus the",
                                 "antiproliferative cost on septation itself"),
                           "how much of the window is left when you start",
                           "nothing - the 36-week deficit is unchanged"),
      check.names = FALSE)
  })
}

shinyApp(ui, server)
