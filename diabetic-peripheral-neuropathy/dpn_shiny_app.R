## =====================================================================
## Diabetic Peripheral Neuropathy (DPN) — Shiny Dashboard
## ---------------------------------------------------------------------
## Author: Claude Code Routine — QSP Disease Model Library
## File:   dpn_shiny_app.R
##
## 8 tabs:
##   1) Patient profile (baseline glycaemia, comorbid, genotype)
##   2) Drug PK (plasma concentration-time)
##   3) Mechanistic PD (sorbitol, AGE, ROS, NGF, VASA, IENFD, NCV)
##   4) Pain & sensitisation (NRS, peripheral / central, descending)
##   5) Clinical endpoints (MNSI, TCNS, BPI, Norfolk-QoL, sleep)
##   6) Scenario comparison (vs no-treatment baseline)
##   7) Biomarker panel (HbA1c, IENFD, NCV, sudomotor)
##   8) Outcomes & footulcer hazard (cumulative incidence)
##
## Run:  shiny::runApp("dpn_shiny_app.R")
## Requires: shiny, shinydashboard, mrgsolve, dplyr, tidyr, ggplot2,
##           DT, scales, plotly
## =====================================================================

suppressPackageStartupMessages({
  library(shiny)
  library(shinydashboard)
  library(mrgsolve)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(DT)
  library(scales)
})

mod <- mread("dpn_mrgsolve_model.R")

scenario_events <- function(scn, sim_days) {
  switch(scn,
    "Untreated" = ev(amt = 0, time = 0),
    "Pregabalin 300 mg BID" =
      ev(amt = 150, ii = 12, addl = floor(sim_days/0.5), cmt = "GUT_PG"),
    "Duloxetine 60 mg QD" =
      ev(amt = 60,  ii = 24, addl = sim_days,            cmt = "GUT_DLX"),
    "Pregabalin + Duloxetine" =
      c(ev(amt = 150, ii = 12, addl = floor(sim_days/0.5), cmt = "GUT_PG"),
        ev(amt =  60, ii = 24, addl = sim_days,            cmt = "GUT_DLX")),
    "α-Lipoic acid 600 mg IV→PO" =
      c(ev(amt = 600, ii = 24, addl = 20,            cmt = "PLAS_ALA"),
        ev(amt = 600, ii = 24, addl = sim_days - 21, cmt = "GUT_ALA", time = 21)),
    "Epalrestat 150 mg QD" =
      ev(amt = 150, ii = 24, addl = sim_days, cmt = "GUT_EP"),
    "Capsaicin 8% patch q90d" =
      ev(amt = 1, time = seq(0, sim_days, by = 90), cmt = "CAP_EFF"),
    "Combination + intensive glucose (HbA1c 6.5)" =
      c(ev(amt = 150, ii = 12, addl = floor(sim_days/0.5), cmt = "GUT_PG"),
        ev(amt = 60,  ii = 24, addl = sim_days,            cmt = "GUT_DLX"),
        ev(amt = 600, ii = 24, addl = sim_days,            cmt = "GUT_ALA"))
  )
}

run_sim <- function(scn, patient, sim_days) {
  m <- mod
  pars <- list(
    HbA1c0      = patient$HbA1c0,
    HbA1c_tgt   = patient$HbA1c_tgt,
    DM_DURATION = patient$DM_DURATION,
    AGE_PT      = patient$AGE_PT,
    BMI         = patient$BMI,
    HTN_FLAG    = patient$HTN_FLAG,
    SMOKE_FLAG  = patient$SMOKE_FLAG,
    EGFR        = patient$EGFR,
    GLUC_VAR    = patient$GLUC_VAR,
    DLX_CYP2D6  = patient$CYP2D6
  )
  if (grepl("Combination \\+ intensive", scn)) pars$HbA1c_tgt <- 6.5
  m <- update(m, param = pars)

  init <- list(
    HbA1c      = patient$HbA1c0,
    GSH        = 1.0,
    NGF        = 1.0,
    VASA       = max(0.4, 1.0 - patient$DM_DURATION * 0.03),
    IENFD      = max(2.0, 12 - patient$DM_DURATION * 0.4),
    NCV        = max(30, 50 - patient$DM_DURATION * 0.7),
    DESC_TONE  = 1.0
  )
  m <- init(m, init)

  ev_x <- scenario_events(scn, sim_days)
  out  <- mrgsim(m, events = ev_x, end = sim_days, delta = 1) %>% as.data.frame()
  out$scenario <- scn
  out
}

ui <- dashboardPage(
  dashboardHeader(title = "DPN QSP Dashboard"),
  dashboardSidebar(
    sidebarMenu(
      menuItem("1. Patient Profile",      tabName = "patient",  icon = icon("user")),
      menuItem("2. Drug PK",              tabName = "pk",       icon = icon("flask")),
      menuItem("3. Mechanistic PD",       tabName = "pd",       icon = icon("dna")),
      menuItem("4. Pain & Sensitisation", tabName = "pain",     icon = icon("bolt")),
      menuItem("5. Clinical Endpoints",   tabName = "clin",     icon = icon("stethoscope")),
      menuItem("6. Scenario Comparison",  tabName = "compare",  icon = icon("balance-scale")),
      menuItem("7. Biomarker Panel",      tabName = "biom",     icon = icon("vial")),
      menuItem("8. Outcomes & Hazard",    tabName = "outcome",  icon = icon("notes-medical"))
    ),
    hr(),
    h4("Patient", style = "color:#fff;padding-left:15px;"),
    sliderInput("HbA1c0",      "Baseline HbA1c (%)",   6.0, 12.0, 8.5, 0.1),
    sliderInput("HbA1c_tgt",   "Treatment target (%)", 5.5, 9.5,  7.0, 0.1),
    sliderInput("DM_DURATION", "Diabetes duration (y)", 0, 30, 10, 1),
    sliderInput("AGE_PT",      "Age (y)", 30, 85, 60),
    sliderInput("BMI",         "BMI (kg/m²)", 18, 45, 30, 0.5),
    sliderInput("EGFR",        "eGFR (mL/min/1.73m²)", 15, 120, 75, 5),
    sliderInput("GLUC_VAR",    "Glycaemic variability×", 0.5, 2.5, 1.0, 0.1),
    selectInput("CYP2D6",      "CYP2D6 phenotype",
                choices = c("Poor (0.4)"=0.4, "Extensive (1.0)"=1.0, "Ultra (2.0)"=2.0),
                selected = 1.0),
    checkboxInput("HTN_FLAG",  "Hypertension", TRUE),
    checkboxInput("SMOKE_FLAG","Current smoker", FALSE),
    hr(),
    h4("Simulation", style = "color:#fff;padding-left:15px;"),
    sliderInput("sim_days", "Days", 30, 1095, 365, 30),
    selectInput("scn", "Scenario",
                choices = c("Untreated",
                            "Pregabalin 300 mg BID",
                            "Duloxetine 60 mg QD",
                            "Pregabalin + Duloxetine",
                            "α-Lipoic acid 600 mg IV→PO",
                            "Epalrestat 150 mg QD",
                            "Capsaicin 8% patch q90d",
                            "Combination + intensive glucose (HbA1c 6.5)"),
                selected = "Pregabalin 300 mg BID"),
    actionButton("run", "Run simulation", icon = icon("play"),
                 style = "color:#fff;background:#117864;margin:10px;")
  ),
  dashboardBody(
    tabItems(
      tabItem("patient",
        fluidRow(
          box(width = 12, status = "primary", title = "Patient & disease profile",
              tableOutput("ptbl"),
              tags$small("Initial nerve compartments scale with diabetes duration.
              Adjust sliders, then press 'Run simulation'."))
        )
      ),
      tabItem("pk",
        fluidRow(
          box(width = 12, status = "info", title = "Plasma concentrations",
              plotOutput("pk_plot", height = 480))
        )
      ),
      tabItem("pd",
        fluidRow(
          box(width = 6, title = "Polyol & AGE",            plotOutput("pd_age")),
          box(width = 6, title = "ROS & GSH",               plotOutput("pd_ros")),
          box(width = 6, title = "Vasa nervorum / Hypoxia", plotOutput("pd_vasa")),
          box(width = 6, title = "NGF & IENFD",             plotOutput("pd_nerve"))
        )
      ),
      tabItem("pain",
        fluidRow(
          box(width = 6, title = "Pain NRS (0-10)",     plotOutput("pn_nrs")),
          box(width = 6, title = "Peripheral / Central sensitisation", plotOutput("pn_sens")),
          box(width = 6, title = "Descending tone (DLX response)",     plotOutput("pn_desc")),
          box(width = 6, title = "Sleep interference",                 plotOutput("pn_sleep"))
        )
      ),
      tabItem("clin",
        fluidRow(
          box(width = 6, title = "MNSI",          plotOutput("cl_mnsi")),
          box(width = 6, title = "TCNS",          plotOutput("cl_tcns")),
          box(width = 6, title = "BPI interference", plotOutput("cl_bpi")),
          box(width = 6, title = "Norfolk QoL",   plotOutput("cl_qol"))
        )
      ),
      tabItem("compare",
        fluidRow(
          box(width = 12, title = "All scenarios — NRS, IENFD, NCV, Norfolk-QoL",
              plotOutput("cmp_plot", height = 580))
        )
      ),
      tabItem("biom",
        fluidRow(
          box(width = 6, title = "HbA1c trajectory", plotOutput("bm_hba1c")),
          box(width = 6, title = "Sural NCV",        plotOutput("bm_ncv")),
          box(width = 6, title = "IENFD",            plotOutput("bm_ienfd")),
          box(width = 6, title = "AGE burden",       plotOutput("bm_age"))
        )
      ),
      tabItem("outcome",
        fluidRow(
          box(width = 12, title = "Foot-ulcer cumulative incidence",
              plotOutput("ou_hazard", height = 360)),
          box(width = 12, title = "Endpoint summary table", DTOutput("ou_tbl"))
        )
      )
    )
  )
)

server <- function(input, output, session) {
  patient <- reactive({
    list(HbA1c0 = input$HbA1c0, HbA1c_tgt = input$HbA1c_tgt,
         DM_DURATION = input$DM_DURATION, AGE_PT = input$AGE_PT,
         BMI = input$BMI, EGFR = input$EGFR,
         HTN_FLAG = as.integer(input$HTN_FLAG),
         SMOKE_FLAG = as.integer(input$SMOKE_FLAG),
         GLUC_VAR = input$GLUC_VAR,
         CYP2D6 = as.numeric(input$CYP2D6))
  })

  sim_main <- eventReactive(input$run, {
    run_sim(input$scn, patient(), input$sim_days)
  }, ignoreNULL = FALSE)

  sim_all <- eventReactive(input$run, {
    scns <- c("Untreated","Pregabalin 300 mg BID","Duloxetine 60 mg QD",
              "Pregabalin + Duloxetine","α-Lipoic acid 600 mg IV→PO",
              "Epalrestat 150 mg QD","Capsaicin 8% patch q90d",
              "Combination + intensive glucose (HbA1c 6.5)")
    do.call(rbind, lapply(scns, function(s) run_sim(s, patient(), input$sim_days)))
  }, ignoreNULL = FALSE)

  output$ptbl <- renderTable({
    p <- patient()
    data.frame(Variable = names(p), Value = as.character(unlist(p)))
  })

  base_theme <- theme_minimal(base_size = 12)

  output$pk_plot <- renderPlot({
    s <- sim_main()
    long <- s %>% select(time, Cp_PG, Cp_DLX, Cp_ALA, Cp_EP) %>%
      pivot_longer(-time, names_to = "drug", values_to = "Cp")
    ggplot(long, aes(time, Cp, color = drug)) + geom_line(size = 0.8) +
      labs(x = "Days", y = "Plasma concentration (mg/L)") + base_theme
  })

  pd_plot <- function(s, vars, ylab) {
    s %>% select(time, all_of(vars)) %>%
      pivot_longer(-time) %>%
      ggplot(aes(time, value, color = name)) + geom_line(size = 0.8) +
      labs(x = "Days", y = ylab) + base_theme
  }
  output$pd_age   <- renderPlot(pd_plot(sim_main(), c("SORB","AGE"), "Burden (a.u.)"))
  output$pd_ros   <- renderPlot(pd_plot(sim_main(), c("ROS","GSH"),  "Level (a.u.)"))
  output$pd_vasa  <- renderPlot(pd_plot(sim_main(), c("VASA","hypoxia"), "Index"))
  output$pd_nerve <- renderPlot(pd_plot(sim_main(), c("NGF","IENFD"), "NGF / IENFD"))

  output$pn_nrs   <- renderPlot(pd_plot(sim_main(), c("NRS"), "NRS (0-10)") +
                                  ylim(0,10))
  output$pn_sens  <- renderPlot(pd_plot(sim_main(), c("PERIPH_S","CENTRAL_S"), "a.u."))
  output$pn_desc  <- renderPlot(pd_plot(sim_main(), c("DESC_TONE"), "Descending tone"))
  output$pn_sleep <- renderPlot(pd_plot(sim_main(), c("SLEEP_INTERF"), "Sleep interf."))

  output$cl_mnsi <- renderPlot(pd_plot(sim_main(), c("MNSI"), "MNSI score"))
  output$cl_tcns <- renderPlot(pd_plot(sim_main(), c("TCNS"), "TCNS score"))
  output$cl_bpi  <- renderPlot(pd_plot(sim_main(), c("BPI_INT"), "BPI interference"))
  output$cl_qol  <- renderPlot(pd_plot(sim_main(), c("Norfolk_QoL"), "Norfolk QoL (0-100)"))

  output$bm_hba1c <- renderPlot(pd_plot(sim_main(), "HbA1c", "HbA1c (%)"))
  output$bm_ncv   <- renderPlot(pd_plot(sim_main(), "NCV", "NCV (m/s)"))
  output$bm_ienfd <- renderPlot(pd_plot(sim_main(), "IENFD", "IENFD (fibres/mm)"))
  output$bm_age   <- renderPlot(pd_plot(sim_main(), "AGE", "AGE burden"))

  output$ou_hazard <- renderPlot({
    s <- sim_main()
    ggplot(s, aes(time, DFU_INC)) + geom_line(size = 0.9, color = "#C0392B") +
      labs(x = "Days", y = "Cumulative foot-ulcer incidence") +
      base_theme + scale_y_continuous(labels = percent)
  })

  output$cmp_plot <- renderPlot({
    s <- sim_all() %>%
      select(scenario, time, NRS, IENFD, NCV, Norfolk_QoL) %>%
      pivot_longer(c(NRS, IENFD, NCV, Norfolk_QoL))
    ggplot(s, aes(time, value, color = scenario)) +
      geom_line(size = 0.8) + facet_wrap(~ name, scales = "free_y") +
      base_theme + theme(legend.position = "bottom")
  })

  output$ou_tbl <- renderDT({
    s <- sim_all() %>%
      group_by(scenario) %>%
      filter(time == max(time)) %>%
      summarise(
        NRS    = round(NRS, 2),
        MNSI   = round(MNSI, 2),
        TCNS   = round(TCNS, 2),
        IENFD  = round(IENFD, 2),
        NCV    = round(NCV, 2),
        QoL    = round(Norfolk_QoL, 1),
        DFU_inc_pct = scales::percent(DFU_INC, 0.1)
      )
    datatable(s, options = list(dom = 't', pageLength = 10), rownames = FALSE)
  })
}

shinyApp(ui, server)
