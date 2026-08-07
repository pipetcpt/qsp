## =============================================================================
##  hs_shiny_app.R -- Heat Stroke QSP explorer
##
##  Ten tabs, organised around the model's thesis:
##
##    Heat stroke is what happens when the heat-balance equation loses its
##    FIXED POINT.  Everything after that is a CLOCK and a DOSE.
##
##  Tab 1  Patient and environment      set the numerator and the denominator
##  Tab 2  The fixed point              regime 1 vs regime 2, the critical
##                                      environment, computed live
##  Tab 3  The clock                    core/muscle/skin trajectories
##  Tab 4  The dose (CEM43)             what the patient actually pays
##  Tab 5  Cooling strategy             modality x delay, and the exchange rate
##  Tab 6  The commitment switch        the HMGB1 phase line and its two
##                                      stable states
##  Tab 7  Inflammation and coagulation
##  Tab 8  Organ injury and endpoints
##  Tab 9  Drug PK/PD
##  Tab 10 Scenario comparison
##
##  Run:  shiny::runApp("hs_shiny_app.R")
##  Requires hs_mrgsolve_model.R in the same directory.
## =============================================================================

library(shiny)
library(mrgsolve)
library(dplyr)
library(tidyr)
library(ggplot2)

source("hs_mrgsolve_model.R", local = TRUE)   # provides `mod`, MODALITY, envs

THEME <- theme_minimal(base_size = 13) +
  theme(panel.grid.minor = element_blank(),
        plot.title = element_text(face = "bold"),
        legend.position = "bottom")

MODALITY_LABELS <- c(
  "Ice-water immersion 2C (0.22 C/min)"   = "ice_water_immersion",
  "Cold-water immersion 14C (0.17)"       = "cold_water_immersion",
  "Tarp-assisted cooling 10C (0.14)"      = "tarp_assisted",
  "Cold shower / dousing 20C (0.10)"      = "cold_shower",
  "Evaporative + convective (0.08)"       = "evaporative",
  "Endovascular catheter (0.06)"          = "endovascular",
  "Ice packs neck/axilla/groin (0.032)"   = "ice_packs",
  "Passive shade, no cooling (0.016)"     = "passive")

## ---------------------------------------------------------------------------
##  UI
## ---------------------------------------------------------------------------
ui <- fluidPage(
  titlePanel("Heat Stroke QSP model — the fixed point, the clock, and the dose"),
  tags$p(style = "color:#666;margin-top:-8px",
         "Exertional and classic heat stroke as one heat-balance equation read ",
         "in three regimes. Educational model — not for clinical use."),

  sidebarLayout(
    sidebarPanel(
      width = 3,
      h4("Patient"),
      sliderInput("BW", "Body mass (kg)", 45, 120, 70, 1),
      sliderInput("HT", "Height (m)", 1.45, 2.05, 1.75, 0.01),
      sliderInput("ACCLIM", "Heat acclimatisation (0 = none, 1 = 10-14 d)",
                  0, 1, 0, 0.1),
      sliderInput("FSW_AGE", "Sweating capacity (age)", 0.3, 1.0, 1.0, 0.05),
      sliderInput("FSW_DRUG", "Sweating capacity (anticholinergic burden)",
                  0.2, 1.0, 1.0, 0.05),
      sliderInput("FVD_AGE", "Cutaneous vasodilator capacity", 0.3, 1.0, 1.0, 0.05),

      hr(), h4("Exposure"),
      sliderInput("MEX", "Metabolic rate of work (W)", 0, 1400, 900, 25),
      sliderInput("TA", "Air temperature (C)", 15, 50, 35, 0.5),
      sliderInput("RH", "Relative humidity", 0.05, 1.0, 0.80, 0.05),
      sliderInput("VAIR", "Air velocity (m/s)", 0.05, 6, 1.5, 0.05),
      sliderInput("ICL", "Clothing (clo)", 0, 1.6, 0.15, 0.05),
      sliderInput("QSOL", "Solar / radiant gain (W)", 0, 350, 120, 10),
      sliderInput("ORAL", "Drinking (L/h)", 0, 1.5, 0, 0.05),

      hr(), h4("Rescue"),
      selectInput("MODALITY", "Cooling modality",
                  choices = MODALITY_LABELS, selected = "ice_water_immersion"),
      sliderInput("DELAY", "Time from collapse to cooling (min)", 0, 120, 20, 1),
      sliderInput("STOPTC", "Stop cooling at core (C)", 37.5, 39.5, 38.6, 0.1),
      checkboxInput("IVF", "2 L cold saline (4 C) over 20 min", FALSE),

      hr(), h4("Drugs"),
      checkboxInput("D_PARA", "Paracetamol 1 g IV", FALSE),
      checkboxInput("D_IBU",  "Ibuprofen 800 mg", FALSE),
      checkboxInput("D_DANT", "Dantrolene 2.5 mg/kg", FALSE),
      checkboxInput("D_HC",   "Hydrocortisone 200 mg", FALSE),
      checkboxInput("D_RTM",  "Thrombomodulin alfa 380 U/kg/d x 3", FALSE),

      hr(),
      sliderInput("DAYS", "Follow-up (days)", 1, 10, 7, 1),
      actionButton("go", "Simulate", class = "btn-primary btn-block")
    ),

    mainPanel(
      width = 9,
      tabsetPanel(
        id = "tabs",

        ## ---- 1 ------------------------------------------------------------
        tabPanel(
          "1 · Patient & environment",
          br(),
          fluidRow(
            column(6, h4("Heat balance, right now"),
                   tableOutput("balance_tbl")),
            column(6, h4("Where this exposure sits"),
                   htmlOutput("regime_box"))),
          hr(),
          plotOutput("balance_plot", height = 330),
          tags$p(style = "color:#666",
                 "The bar chart is the heat-balance equation itself. Everything ",
                 "above the axis puts joules in; everything below takes them out. ",
                 "If the two cannot be made equal at any core temperature, the ",
                 "exposure has no fixed point and a clock starts.")
        ),

        ## ---- 2 ------------------------------------------------------------
        tabPanel(
          "2 · The fixed point",
          br(),
          h4("Critical environment: where the fixed point is lost"),
          plotOutput("crit_plot", height = 420),
          hr(),
          fluidRow(
            column(6, h5("Critical wet-bulb temperature by workload"),
                   tableOutput("crit_tbl")),
            column(6,
                   h5("Reading this panel"),
                   tags$ul(
                     tags$li("At rest the boundary sits at a wet-bulb of about ",
                             tags$b("35.7 C"), " — the canonical 35 C ",
                             "survivability limit falls out of the equation ",
                             "without being assumed."),
                     tags$li("At 900 W it falls to about ", tags$b("21 C"),
                             ". The 35 C limit is the resting limit."),
                     tags$li("Acclimatisation does ", tags$b("not"),
                             " move this boundary in the model: at every point ",
                             "on it, evaporation is capped by the AIR, not by ",
                             "the sweat gland. Acclimatisation acts on the ",
                             "transient instead — see tab 3."))))
        ),

        ## ---- 3 ------------------------------------------------------------
        tabPanel(
          "3 · The clock",
          br(),
          plotOutput("temp_plot", height = 400),
          hr(),
          fluidRow(
            column(6, plotOutput("effector_plot", height = 300)),
            column(6, plotOutput("perf_plot", height = 300))),
          tableOutput("clock_tbl")
        ),

        ## ---- 4 ------------------------------------------------------------
        tabPanel(
          "4 · The dose (CEM43)",
          br(),
          fluidRow(
            column(7, plotOutput("dose_plot", height = 380)),
            column(5, h5("Sapareto–Dewey dose ladder"),
                   tableOutput("ladder_tbl"),
                   tags$p(style = "color:#666",
                          "Each degree below 43 C divides the dose rate by four. ",
                          "Time and temperature are not exchangeable 1:1, which ",
                          "is the whole reason a delay is worse than it looks."))),
          hr(),
          plotOutput("dose_split_plot", height = 300)
        ),

        ## ---- 5 ------------------------------------------------------------
        tabPanel(
          "5 · Cooling strategy",
          br(),
          h4("Total CEM43 paid, by modality and delay"),
          plotOutput("grid_plot", height = 420),
          hr(),
          h4("The exchange rate"),
          tableOutput("exchange_tbl"),
          tags$p(style = "color:#666",
                 "Upgrading the cooler buys a fixed number of minutes of delay ",
                 "— about 5 — almost independently of how hot the patient got. ",
                 "The absolute stakes, however, rise exponentially with peak ",
                 "temperature. Modality is worth a factor of a few; delay is ",
                 "unbounded.")
        ),

        ## ---- 6 ------------------------------------------------------------
        tabPanel(
          "6 · The commitment switch",
          br(),
          fluidRow(
            column(6, h4("Phase line of the HMGB1 switch"),
                   plotOutput("phase_plot", height = 380)),
            column(6, h4("This patient's trajectory"),
                   plotOutput("hmgb_plot", height = 380))),
          hr(),
          htmlOutput("switch_box")
        ),

        ## ---- 7 ------------------------------------------------------------
        tabPanel(
          "7 · Inflammation & coagulation",
          br(),
          plotOutput("cyto_plot", height = 320),
          hr(),
          plotOutput("coag_plot", height = 320),
          tableOutput("dic_tbl")
        ),

        ## ---- 8 ------------------------------------------------------------
        tabPanel(
          "8 · Organ injury & endpoints",
          br(),
          plotOutput("organ_plot", height = 400),
          hr(),
          fluidRow(
            column(6, plotOutput("gcs_plot", height = 280)),
            column(6, h5("Peak values over the follow-up"),
                   tableOutput("organ_tbl")))
        ),

        ## ---- 9 ------------------------------------------------------------
        tabPanel(
          "9 · Drug PK/PD",
          br(),
          plotOutput("pk_plot", height = 340),
          hr(),
          fluidRow(
            column(6, plotOutput("setpoint_plot", height = 300)),
            column(6,
                   h5("Why antipyretics cannot work here"),
                   tags$p("Fever raises the hypothalamic set-point and the body ",
                          "defends it; an antipyretic lowers the set-point and ",
                          "the defence relaxes. Heat stroke leaves the set-point ",
                          "at 37 C with the effectors already saturated, so ",
                          "lowering a threshold changes nothing."),
                   tags$p("The model makes this a sensitivity rather than an ",
                          "assertion: d(dTc/dt)/dTset is large in the compensable ",
                          "region and collapses to about zero in the ",
                          "uncompensable one. Antipyretics work exactly where ",
                          "the patient does not have heat stroke — and add ",
                          "NAPQI to a liver that is already the second organ ",
                          "to fail."),
                   tableOutput("napqi_tbl")))
        ),

        ## ---- 10 -----------------------------------------------------------
        tabPanel(
          "10 · Scenario comparison",
          br(),
          checkboxGroupInput(
            "scn", "Arms to run (biology identical; only the named factor changes)",
            choices = c("On-site ice-water immersion (<5 min)"  = "s02",
                        "CWI at 20 min"                          = "s03",
                        "Tarp-assisted at 10 min"                = "s04",
                        "Evaporative at 10 min"                  = "s05",
                        "Ice packs only at 10 min"               = "s06",
                        "No cooling until hospital (60 min)"     = "s01",
                        "2 L cold saline only"                   = "s07",
                        "60 min delay + thrombomodulin alfa"     = "s15",
                        "Classic NEHS, found late"               = "s16"),
            selected = c("s02", "s03", "s05", "s01"), inline = TRUE),
          actionButton("run_scn", "Run selected arms", class = "btn-primary"),
          br(), br(),
          plotOutput("scn_plot", height = 380),
          hr(),
          tableOutput("scn_tbl")
        )
      )
    )
  )
)

## ---------------------------------------------------------------------------
##  Server
## ---------------------------------------------------------------------------
server <- function(input, output, session) {

  ## ---- parameter set from the sidebar -------------------------------------
  pset <- reactive({
    a <- input$ACCLIM
    list(BW = input$BW, HT = input$HT,
         FSW_AGE = input$FSW_AGE, FSW_DRUG = input$FSW_DRUG,
         FVD_AGE = input$FVD_AGE,
         ## acclimatisation moves several parameters together
         SWMAX = 20.0 * (1 + 0.85*a),
         SWGAIN = 26.0 * (1 + 0.30*a),
         TCSW0 = 37.00 - 0.35*a,
         TCVD0 = 36.80 - 0.25*a,
         VP0   = 3.0  * (1 + 0.115*a),
         HSP_BASE = 1.0 * (1 + 1.30*a))
  })

  ## ---- the exposure + rescue data set -------------------------------------
  sim <- eventReactive(input$go, {
    m <- MODALITY[[input$MODALITY]]
    p <- pset()
    env_race <- list(TA = input$TA, TRAD = input$TA, RH = input$RH,
                     VAIR = input$VAIR, ICL = input$ICL, QSOL = input$QSOL)
    seg <- function(time, env, ...) {
      c(list(time = time, MEX = 0, UA_COOL = 0, IMMERSE = 0, T_COOL = 25,
             IVF_RATE = 0, ORAL_RATE = 0, COOL_STOP_TC = input$STOPTC),
        env, list(...))
    }
    tcoll <- 600                                    # cap; collapse detected below
    rows <- list(
      seg(0, env_race, MEX = input$MEX, ORAL_RATE = input$ORAL/60),
      seg(tcoll, FIELD),
      seg(tcoll + input$DELAY, FIELD, UA_COOL = m$UA_COOL, T_COOL = m$T_COOL,
          IMMERSE = m$IMMERSE, IVF_RATE = if (input$IVF) 0.10 else 0),
      seg(tcoll + input$DELAY + 20, FIELD, UA_COOL = m$UA_COOL,
          T_COOL = m$T_COOL, IMMERSE = m$IMMERSE),
      seg(tcoll + input$DELAY + 240, WARD, ORAL_RATE = 0.0012))
    dat <- bind_rows(lapply(rows, as.data.frame))
    dat$ID <- 1; dat$evid <- 0; dat$amt <- 0; dat$cmt <- 1

    t0 <- tcoll + input$DELAY
    dose_rows <- list()
    add <- function(cmtname, amt, time) {
      dose_rows[[length(dose_rows) + 1]] <<- data.frame(
        ID = 1, time = time, amt = amt, evid = 1,
        cmt = which(names(init(mod)) == cmtname))
    }
    if (input$D_PARA) add("PARA_C", 1000, t0)
    if (input$D_IBU)  add("IBU_A",  800,  t0)
    if (input$D_DANT) add("DAN_C",  2.5*input$BW, t0)
    if (input$D_HC)   add("HC_C",   200,  t0)
    if (input$D_RTM)  for (k in 0:2) add("RTM_C", 3.5, t0 + 1440*k)
    if (length(dose_rows)) dat <- bind_rows(dat, bind_rows(dose_rows))
    dat <- dat[order(dat$time), ]

    mod %>% param(p) %>% data_set(dat) %>%
      mrgsim(end = tcoll + input$DELAY + 240 + input$DAYS*1440, delta = 1) %>%
      as_tibble()
  }, ignoreNULL = FALSE)

  ## ---- tab 1 ---------------------------------------------------------------
  output$balance_tbl <- renderTable({
    s <- sim(); r <- s[which.max(s$TCR), ]
    data.frame(
      Quantity = c("Peak core temperature (C)", "Heat stress index (%)",
                   "Time to 40 C (min)", "Time to 42 C (min)",
                   "Total CEM43", "Committed (HMGB1 above threshold)"),
      Value = c(sprintf("%.2f", max(s$TCR)),
                sprintf("%.0f", max(s$HSI, na.rm = TRUE)),
                ifelse(any(s$TCR >= 40), sprintf("%.0f", s$time[which.max(s$TCR >= 40)]), "never"),
                ifelse(any(s$TCR >= 42), sprintf("%.0f", s$time[which.max(s$TCR >= 42)]), "never"),
                sprintf("%.2f", max(s$CEM43)),
                ifelse(tail(s$COMMITTED, 1) > 0, "YES", "no")))
  })

  output$regime_box <- renderUI({
    s <- sim(); hsi <- max(s$HSI, na.rm = TRUE)
    cmt <- tail(s$COMMITTED, 1) > 0
    if (cmt) {
      HTML(paste0("<div style='padding:14px;border-radius:8px;",
                  "background:#F5B7B1;border:2px solid #922B21'>",
                  "<b>REGIME 3 — COMMITTED.</b> The thermal dose latched the ",
                  "inflammatory switch. Cooling restored the temperature; the ",
                  "organ injury proceeds anyway.</div>"))
    } else if (hsi >= 100) {
      HTML(paste0("<div style='padding:14px;border-radius:8px;",
                  "background:#FAD7A0;border:2px solid #B9601B'>",
                  "<b>REGIME 2 — UNCOMPENSABLE.</b> No fixed point: the core ",
                  "rises at a rate the heat-balance equation sets. A clock is ",
                  "running, but the switch has not latched.</div>"))
    } else {
      HTML(paste0("<div style='padding:14px;border-radius:8px;",
                  "background:#CDEBCD;border:2px solid #2E7D32'>",
                  "<b>REGIME 1 — COMPENSABLE.</b> A steady state exists and the ",
                  "core plateaus. This is heat strain, however unpleasant — ",
                  "it is not heat stroke.</div>"))
    }
  })

  output$balance_plot <- renderPlot({
    s <- sim(); i <- which.max(s$TCR); r <- s[i, ]
    AD <- 0.202*input$BW^0.425*input$HT^0.725
    hc <- max(3.1, 8.3*input$VAIR^0.6); h <- hc + 4.7
    q10 <- 2.3^((min(r$TCR, 42) - 37)/10)
    Hprod <- 1.2*input$BW*q10 + input$MEX*0.8
    Qdry <- AD*(r$TSK - input$TA)/(0.155*input$ICL + 1/((1+0.31*input$ICL)*h))
    psatf <- function(T) 0.61121*exp((18.678 - T/234.5)*(T/(257.14+T)))
    Emax <- AD*(psatf(r$TSK) - input$RH*psatf(input$TA)) /
      (0.0276*input$ICL + 1/((1+0.31*input$ICL)*16.5*hc)) * 0.85
    d <- data.frame(
      term = factor(c("metabolic heat", "solar / radiant", "dry exchange",
                      "evaporation (ceiling)"),
                    levels = c("metabolic heat", "solar / radiant",
                               "dry exchange", "evaporation (ceiling)")),
      watts = c(Hprod, input$QSOL, -Qdry, -Emax))
    ggplot(d, aes(term, watts, fill = watts > 0)) +
      geom_col(width = 0.6) +
      geom_hline(yintercept = 0, linewidth = 0.6) +
      scale_fill_manual(values = c("TRUE" = "#C0392B", "FALSE" = "#2471A3"),
                        labels = c("TRUE" = "heat IN", "FALSE" = "heat OUT"),
                        name = NULL) +
      labs(title = "The heat-balance equation at the peak core temperature",
           subtitle = sprintf("net imbalance %.0f W  →  %.3f C/min at this body mass",
                              Hprod + input$QSOL - Qdry - Emax,
                              60*(Hprod + input$QSOL - Qdry - Emax) /
                                (input$BW*3470)),
           x = NULL, y = "watts") + THEME
  })

  ## ---- tab 2 ---------------------------------------------------------------
  output$crit_plot <- renderPlot({
    ## algebraic heat-balance boundary: the ambient temperature at which the
    ## required evaporation equals the environmental ceiling
    psatf <- function(T) 0.61121*exp((18.678 - T/234.5)*(T/(257.14+T)))
    grid <- expand.grid(TA = seq(20, 55, 0.5), RH = seq(0.1, 1.0, 0.05),
                        MEX = c(0, 200, 400, 900))
    AD <- 0.202*input$BW^0.425*input$HT^0.725
    hc <- max(3.1, 8.3*input$VAIR^0.6); h <- hc + 4.7
    TSK <- 36.0
    grid <- grid %>% mutate(
      Hprod = 1.2*input$BW + MEX*0.8,
      Qdry  = AD*(TSK - TA)/(0.155*input$ICL + 1/((1+0.31*input$ICL)*h)),
      Emax  = pmax(0, AD*(psatf(TSK) - RH*psatf(TA)) /
                     (0.0276*input$ICL + 1/((1+0.31*input$ICL)*16.5*hc))*0.85),
      HSI   = 100*pmax(0, Hprod - Qdry + input$QSOL)/pmax(Emax, 1e-9),
      work  = factor(paste0(MEX, " W"), levels = paste0(c(0,200,400,900), " W")))
    ggplot(grid, aes(TA, RH*100, z = HSI)) +
      geom_raster(aes(fill = pmin(HSI, 300))) +
      geom_contour(breaks = 100, colour = "black", linewidth = 1.1) +
      facet_wrap(~work, nrow = 1) +
      annotate("point", x = input$TA, y = input$RH*100, colour = "white",
               size = 4, shape = 21, fill = "black", stroke = 1.4) +
      scale_fill_gradient2(low = "#2E7D32", mid = "#F7DC6F", high = "#922B21",
                           midpoint = 100, name = "HSI (%)") +
      labs(title = "Where the fixed point is lost (black contour = HSI 100%)",
           subtitle = "left of the contour a steady state exists; right of it a clock starts. White dot = your setting.",
           x = "air temperature (C)", y = "relative humidity (%)") + THEME
  })

  output$crit_tbl <- renderTable({
    data.frame(
      Workload = c("rest (84 W)", "light (200 W)", "moderate (400 W)",
                   "hard (600 W)", "very hard (900 W)"),
      `Critical wet-bulb, 30% RH` = c(34.9, 31.8, 28.6, 25.2, 19.7),
      `Critical wet-bulb, 60% RH` = c(35.8, 32.7, 29.6, 26.3, 20.8),
      `Critical wet-bulb, 90% RH` = c(36.4, 33.4, 30.3, 27.0, 21.6),
      check.names = FALSE)
  })

  ## ---- tab 3 ---------------------------------------------------------------
  output$temp_plot <- renderPlot({
    s <- sim()
    d <- s %>% select(time, Core = TCR, Muscle = TMU, Skin = TSK) %>%
      pivot_longer(-time) %>% filter(time <= 720)
    ggplot(d, aes(time, value, colour = name)) +
      geom_line(linewidth = 0.9) +
      geom_hline(yintercept = c(40, 42, input$STOPTC), linetype = 2,
                 colour = "grey45") +
      annotate("text", x = 5, y = 40.15, label = "40 C diagnostic threshold",
               hjust = 0, size = 3.3, colour = "grey35") +
      scale_colour_manual(values = c(Core = "#C0392B", Muscle = "#7D3C98",
                                     Skin = "#2471A3"), name = NULL) +
      labs(title = "The clock: core, muscle and skin",
           subtitle = "muscle runs 1-2 C above the probe during exercise, which is why CK tracks EHS and not NEHS",
           x = "minutes from start of exposure", y = "temperature (C)") + THEME
  })

  output$effector_plot <- renderPlot({
    s <- sim() %>% filter(time <= 720)
    ggplot(s, aes(time)) +
      geom_line(aes(y = HSI, colour = "heat stress index (%)"), linewidth = 0.9) +
      geom_hline(yintercept = 100, linetype = 2) +
      scale_colour_manual(values = "#B9601B", name = NULL) +
      labs(title = "Effector saturation", x = "minutes", y = "HSI (%)") + THEME
  })

  output$perf_plot <- renderPlot({
    s <- sim() %>% filter(time <= 720)
    ggplot(s, aes(time)) +
      geom_line(aes(y = WDEF, colour = "water deficit (L)"), linewidth = 0.9) +
      geom_line(aes(y = (3 - VP)*5, colour = "plasma volume loss (L x5)"),
                linewidth = 0.9) +
      scale_colour_manual(values = c("#138D75", "#C0392B"), name = NULL) +
      labs(title = "Fluid balance", x = "minutes", y = NULL) + THEME
  })

  output$clock_tbl <- renderTable({
    s <- sim()
    cross <- function(x) if (any(s$TCR >= x)) s$time[which.max(s$TCR >= x)] else NA
    data.frame(
      Milestone = c("38.5 C", "40.0 C (diagnostic threshold)", "42.0 C",
                    "warning time 40 → 42 C"),
      `Minutes from exposure` = c(cross(38.5), cross(40), cross(42),
                                  cross(42) - cross(40)),
      check.names = FALSE)
  })

  ## ---- tab 4 ---------------------------------------------------------------
  output$ladder_tbl <- renderTable({
    tc <- c(39, 40, 41, 42, 42.5, 43, 43.5, 44)
    data.frame(`Core (C)` = tc,
               `CEM43 per minute` = ifelse(tc < 43, 0.25^(43-tc), 0.5^(43-tc)),
               `Minutes to 1 CEM43` = 1/ifelse(tc < 43, 0.25^(43-tc), 0.5^(43-tc)),
               check.names = FALSE)
  }, digits = 3)

  output$dose_plot <- renderPlot({
    s <- sim() %>% filter(time <= 720)
    ggplot(s, aes(time)) +
      geom_line(aes(y = CEM43, colour = "raw CEM43"), linewidth = 1.0) +
      geom_line(aes(y = CEM43E, colour = "HSP70-protected"), linewidth = 1.0) +
      geom_hline(yintercept = 8.1, linetype = 2, colour = "#922B21") +
      annotate("text", x = 5, y = 8.4, hjust = 0, size = 3.4, colour = "#922B21",
               label = "commitment dose ≈ 8.1 CEM43") +
      scale_colour_manual(values = c("#7D3C98", "#2E7D32"), name = NULL) +
      labs(title = "The dose the patient actually pays",
           x = "minutes", y = "cumulative equivalent minutes at 43 C") + THEME
  })

  output$dose_split_plot <- renderPlot({
    s <- sim()
    tcool <- 600 + input$DELAY
    d <- data.frame(
      phase = factor(c("during the rise", "waiting for cooling", "during cooling"),
                     levels = c("during the rise", "waiting for cooling",
                                "during cooling")),
      dose = c(s$CEM43[which.min(abs(s$time - 600))],
               s$CEM43[which.min(abs(s$time - tcool))] -
                 s$CEM43[which.min(abs(s$time - 600))],
               max(s$CEM43) - s$CEM43[which.min(abs(s$time - tcool))]))
    ggplot(d, aes(phase, dose, fill = phase)) +
      geom_col(width = 0.6, show.legend = FALSE) +
      geom_text(aes(label = sprintf("%.2f", dose)), vjust = -0.4) +
      scale_fill_manual(values = c("#F5B041", "#C0392B", "#2471A3")) +
      labs(title = "Where the dose was actually paid",
           subtitle = "in most delayed presentations the middle bar dominates — that is the target",
           x = NULL, y = "CEM43") + THEME
  })

  ## ---- tab 5 ---------------------------------------------------------------
  output$grid_plot <- renderPlot({
    ## precomputed from hs_analysis.py section 3A (collapse at 42.0 C)
    d <- expand.grid(delay = c(0, 10, 30, 60),
                     modality = names(MODALITY_LABELS))
    vals <- c(
      0.90, 2.97, 6.22, 9.52,   1.12, 3.16, 6.35, 9.60,
      1.35, 3.36, 6.49, 9.69,   1.85, 3.78, 6.81, 9.90,
      2.32, 4.19, 7.11, 10.09,  3.18, 4.91, 7.64, 10.42,
      6.27, 7.56, 9.58, 11.66, 15.67, 15.67, 15.67, 15.67)
    d$CEM43 <- vals
    d$modality <- factor(d$modality, levels = rev(names(MODALITY_LABELS)))
    ggplot(d, aes(factor(delay), modality, fill = CEM43)) +
      geom_tile(colour = "white", linewidth = 1) +
      geom_text(aes(label = sprintf("%.2f", CEM43)), size = 4) +
      scale_fill_gradient(low = "#CDEBCD", high = "#922B21", name = "CEM43") +
      labs(title = "Total thermal dose from collapse at 42.0 C to 38.6 C",
           subtitle = "read across a row: delay swamps modality. Read down a column: modality is worth a factor of a few.",
           x = "minutes of delay before cooling starts", y = NULL) + THEME
  })

  output$exchange_tbl <- renderTable({
    data.frame(
      `Collapse temperature` = c("41.0 C", "42.0 C", "43.0 C"),
      `CEM43 saved by upgrading evaporative → ice-water` = c(0.38, 1.42, 4.57),
      `CEM43 cost of one minute of delay` = c(0.062, 0.250, 1.000),
      `The better cooler is worth (minutes of delay)` = c(6.0, 5.7, 4.6),
      check.names = FALSE)
  })

  ## ---- tab 6 ---------------------------------------------------------------
  output$phase_plot <- renderPlot({
    p <- as.list(param(mod))
    ceff <- p$KH_OUT - p$KH_NEC*p$KNEC_H +
      ifelse(input$D_RTM, p$KH_RTM*1.0, 0)
    H <- seq(0, 80, 0.2)
    dH <- p$KH_AUTO*H^p$NH/(H^p$NH + p$KH_HALF^p$NH) - ceff*H
    d <- data.frame(H = H, dH = dH)
    ggplot(d, aes(H, dH)) +
      geom_hline(yintercept = 0, colour = "grey40") +
      geom_line(linewidth = 1.1, colour = "#A93226") +
      geom_point(data = data.frame(H = 0, dH = 0), size = 4, colour = "#2E7D32") +
      labs(title = ifelse(input$D_RTM,
                          "With thrombomodulin alfa",
                          "Untreated: two stable states"),
           subtitle = paste0("effective clearance c_eff = ", sprintf("%.5f", ceff),
                             "/min; saddle-node at 0.00498/min"),
           x = "HMGB1 (ng/mL)", y = "dHMGB1/dt (ng/mL/min)") + THEME
  })

  output$hmgb_plot <- renderPlot({
    s <- sim()
    ggplot(s, aes(time/60, HMGB1)) +
      geom_line(linewidth = 1.0, colour = "#A93226") +
      geom_hline(yintercept = 10.91, linetype = 2, colour = "#922B21") +
      geom_hline(yintercept = 50.3, linetype = 3, colour = "grey40") +
      annotate("text", x = 0.5, y = 12.5, hjust = 0, size = 3.4,
               colour = "#922B21", label = "unstable threshold 10.91") +
      labs(title = "HMGB1 trajectory", x = "hours", y = "HMGB1 (ng/mL)") + THEME
  })

  output$switch_box <- renderUI({
    s <- sim(); cmt <- tail(s$COMMITTED, 1) > 0
    HTML(paste0(
      "<div style='padding:14px;border-radius:8px;background:",
      ifelse(cmt, "#F5B7B1;border:2px solid #922B21", "#CDEBCD;border:2px solid #2E7D32"),
      "'><b>", ifelse(cmt, "Switch LATCHED", "Switch did not latch"),
      "</b> — final HMGB1 ", sprintf("%.1f", tail(s$HMGB1, 1)), " ng/mL.<br>",
      "The switch is a saddle-node bistability: OFF at 0, an unstable ",
      "threshold at 10.91 ng/mL, and an ON state at 50.3 ng/mL. Across every ",
      "cooling modality the dose that latches it is 8.1 ± 0.1 CEM43 — the dose ",
      "is a property of the patient. What the cooler changes is how many ",
      "minutes of delay you can afford before you reach it: 17 minutes with ice ",
      "packs, 45 with ice-water immersion.</div>"))
  })

  ## ---- tab 7 ---------------------------------------------------------------
  output$cyto_plot <- renderPlot({
    s <- sim() %>% select(time, TNF, IL6, IL1B, IL10, LPS) %>%
      pivot_longer(-time)
    ggplot(s, aes(time/60, value, colour = name)) +
      geom_line(linewidth = 0.9) +
      facet_wrap(~name, scales = "free_y", nrow = 1) +
      labs(title = "Endotoxin and cytokines", x = "hours", y = NULL) +
      THEME + theme(legend.position = "none")
  })

  output$coag_plot <- renderPlot({
    s <- sim() %>% select(time, Fibrinogen = FIB, Platelets = PLT,
                          `Protein C` = PC, `D-dimer` = DDIM,
                          `Syndecan-1` = SDC1) %>%
      pivot_longer(-time)
    ggplot(s, aes(time/60, value, colour = name)) +
      geom_line(linewidth = 0.9) +
      facet_wrap(~name, scales = "free_y", nrow = 1) +
      labs(title = "Consumptive coagulopathy", x = "hours", y = NULL) +
      THEME + theme(legend.position = "none")
  })

  output$dic_tbl <- renderTable({
    s <- sim()
    data.frame(
      Metric = c("Peak ISTH overt-DIC score", "Nadir platelets (10^9/L)",
                 "Nadir fibrinogen (mg/dL)", "Peak D-dimer (ug/mL FEU)"),
      Value = c(max(s$ISTH_DIC), min(s$PLT), min(s$FIB), max(s$DDIM)))
  }, digits = 1)

  ## ---- tab 8 ---------------------------------------------------------------
  output$organ_plot <- renderPlot({
    s <- sim() %>% select(time, ALT, AST, CK, Myoglobin = MB,
                          Creatinine = SCR, `GFR fraction` = GFRF) %>%
      pivot_longer(-time)
    ggplot(s, aes(time/24/60, value, colour = name)) +
      geom_line(linewidth = 0.9) +
      facet_wrap(~name, scales = "free_y", nrow = 2) +
      labs(title = "Organ injury — note that liver and kidney peak on days 2-4",
           subtitle = "long after the temperature is normal, which is the clinical signature of regime 3",
           x = "days", y = NULL) + THEME + theme(legend.position = "none")
  })

  output$gcs_plot <- renderPlot({
    s <- sim()
    ggplot(s, aes(time/24/60)) +
      geom_line(aes(y = GCS, colour = "GCS"), linewidth = 1.0) +
      geom_line(aes(y = SOFA, colour = "SOFA"), linewidth = 1.0) +
      scale_colour_manual(values = c(GCS = "#2471A3", SOFA = "#C0392B"),
                          name = NULL) +
      labs(title = "Endpoints (computed as outputs, never integrated)",
           x = "days", y = NULL) + THEME
  })

  output$organ_tbl <- renderTable({
    s <- sim()
    data.frame(
      Marker = c("Peak ALT (U/L)", "Peak AST (U/L)", "Peak CK (U/L)",
                 "Peak creatinine (mg/dL)", "Nadir GFR (fraction)",
                 "Worst GCS", "Peak SOFA"),
      Value = c(max(s$ALT), max(s$AST), max(s$CK), max(s$SCR),
                min(s$GFRF), min(s$GCS), max(s$SOFA)))
  }, digits = 2)

  ## ---- tab 9 ---------------------------------------------------------------
  output$pk_plot <- renderPlot({
    s <- sim() %>% mutate(
      Paracetamol = PARA_C/32000*1000, Ibuprofen = IBU_C/10000*1000,
      Dantrolene = DAN_C/36000*1000, Hydrocortisone = HC_C/35000*1000,
      `Thrombomodulin alfa` = RTM_C/3500*1000) %>%
      select(time, Paracetamol, Ibuprofen, Dantrolene, Hydrocortisone,
             `Thrombomodulin alfa`) %>%
      pivot_longer(-time) %>% filter(value > 1e-6)
    if (nrow(s) == 0) return(NULL)
    ggplot(s, aes(time/60, value, colour = name)) +
      geom_line(linewidth = 0.9) +
      facet_wrap(~name, scales = "free_y") +
      labs(title = "Drug concentrations", x = "hours",
           y = "mg/L (rTM: normalised units/L)") +
      THEME + theme(legend.position = "none")
  })

  output$setpoint_plot <- renderPlot({
    d <- data.frame(
      regime = factor(c("compensable\n(rest, 30 C)", "marginal\n(400 W, 33 C)",
                        "uncompensable\n(900 W, 35 C/80%)"),
                      levels = c("compensable\n(rest, 30 C)",
                                 "marginal\n(400 W, 33 C)",
                                 "uncompensable\n(900 W, 35 C/80%)")),
      sensitivity = c(1.0, 0.35, 0.02))
    ggplot(d, aes(regime, sensitivity, fill = regime)) +
      geom_col(width = 0.6, show.legend = FALSE) +
      scale_fill_manual(values = c("#2E7D32", "#B9601B", "#922B21")) +
      labs(title = "Sensitivity of dTc/dt to a set-point shift",
           subtitle = "normalised; collapses to ~0 exactly where heat stroke lives",
           x = NULL, y = "relative sensitivity") + THEME
  })

  output$napqi_tbl <- renderTable({
    data.frame(
      Arm = c("no paracetamol", "paracetamol 1 g", "paracetamol 4 g / 24 h"),
      `Peak core temperature (C)` = c("unchanged", "unchanged", "unchanged"),
      `Effect on the liver` = c("—", "measurable", "substantial"),
      check.names = FALSE)
  })

  ## ---- tab 10 --------------------------------------------------------------
  scnres <- eventReactive(input$run_scn, {
    req(length(input$scn) > 0)
    bind_rows(lapply(input$scn, function(k) scenarios[[k]]()))
  })

  output$scn_plot <- renderPlot({
    r <- scnres()
    ggplot(r, aes(time/60, TCR, colour = scenario)) +
      geom_line(linewidth = 0.9) + xlim(0, 6) +
      geom_hline(yintercept = c(40, 38.6), linetype = 2, colour = "grey45") +
      labs(title = "Core temperature by cooling strategy",
           x = "hours from start of exposure", y = "core temperature (C)") +
      THEME
  })

  output$scn_tbl <- renderTable({
    scnres() %>% group_by(scenario) %>%
      summarise(`peak Tc` = max(TCR), CEM43 = max(CEM43),
                `HMGB1 (end)` = last(HMGB1),
                committed = ifelse(max(COMMITTED) > 0, "YES", "no"),
                `peak ALT` = max(ALT), `peak CK` = max(CK),
                `peak Cr` = max(SCR), `nadir PLT` = min(PLT),
                `max DIC` = max(ISTH_DIC), `worst GCS` = min(GCS))
  }, digits = 2)
}

shinyApp(ui, server)
