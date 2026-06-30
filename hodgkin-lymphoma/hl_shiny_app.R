# ===========================================================
# Hodgkin Lymphoma QSP - Shiny dashboard
# Eight tabs:
#   1) Patient profile     - stage / IPS / histology inputs
#   2) Drug PK             - concentration-time curves for chosen regimen
#   3) Tumor & TARC        - HRS mass, MTV, serum TARC
#   4) Immune dynamics     - Effector / Exhausted T cells, PD-1 occupancy, Treg
#   5) Hematology          - Friberg ANC profile, neutropenia nadir
#   6) Toxicity            - cardiotox, pulmonary, neuropathy, irAE indices
#   7) Endpoints           - Deauville surrogate, PFS hazard, regimen compare
#   8) Biomarkers          - TARC, IL-6, IL-13, sCD30, lab panel
# ===========================================================

suppressPackageStartupMessages({
  library(shiny)
  library(mrgsolve)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
})

source("hl_mrgsolve_model.R", local = TRUE)   # exposes mod_hl & ev_* / run_scenario

# Helper to translate UI choice to event builder
build_events <- function(regimen, cycles = 6){
  switch(regimen,
         "ABVD"        = ev_ABVD(cycles),
         "AVD"         = ev_AVD(cycles),
         "BV_AVD"      = ev_BV_AVD(cycles),
         "N_AVD"       = ev_N_AVD(cycles),
         "escBEACOPP"  = ev_escBEACOPP(cycles),
         "BV_mono"     = ev_BV_mono(cycles),
         "NIVO_mono"   = ev_NIVO_mono(cycles),
         "BV_NIVO"     = ev_BV_NIVO(cycles))
}

ui <- fluidPage(
  titlePanel("Hodgkin Lymphoma QSP Explorer"),
  sidebarLayout(
    sidebarPanel(
      width = 3,
      h4("Patient profile"),
      selectInput("histology", "Histology",
                  choices = c("Nodular sclerosis (NS)", "Mixed cellularity (MC)",
                              "Lymphocyte-rich (LR)", "Lymphocyte-depleted (LD)",
                              "NLPHL")),
      selectInput("stage", "Lugano stage",
                  choices = c("IA", "IIA early-favorable", "IIB early-unfavorable",
                              "III", "IV"), selected = "III"),
      sliderInput("IPS", "IPS score (advanced)", min = 0, max = 7, value = 2, step = 1),
      sliderInput("baseline_MTV", "Baseline MTV (mL)",
                  min = 100, max = 4000, value = 960, step = 50),
      checkboxInput("EBV", "EBV-positive HRS", value = FALSE),
      checkboxInput("Bsymp", "B-symptoms", value = TRUE),
      hr(),
      h4("Therapy"),
      selectInput("regimen", "Regimen",
                  choices = c("ABVD", "AVD", "BV_AVD", "N_AVD",
                              "escBEACOPP", "BV_mono", "NIVO_mono", "BV_NIVO"),
                  selected = "BV_AVD"),
      sliderInput("cycles", "Cycles", min = 2, max = 8, value = 6, step = 1),
      checkboxInput("GCSF", "G-CSF support", value = TRUE),
      sliderInput("tend", "Simulation horizon (days)",
                  min = 90, max = 730, value = 365, step = 30),
      hr(),
      actionButton("run", "Run simulation", class = "btn-primary")
    ),
    mainPanel(
      tabsetPanel(
        tabPanel("1. Patient",
                 h4("Risk synthesis"),
                 verbatimTextOutput("riskTxt"),
                 plotOutput("riskPlot", height = 280)),
        tabPanel("2. Drug PK",
                 plotOutput("pkPlot", height = 480)),
        tabPanel("3. Tumor / TARC",
                 plotOutput("tumorPlot", height = 280),
                 plotOutput("tarcPlot",  height = 280)),
        tabPanel("4. Immune",
                 plotOutput("immunePlot",  height = 280),
                 plotOutput("pdl1Plot",    height = 280)),
        tabPanel("5. Hematology",
                 plotOutput("ancPlot",     height = 280),
                 verbatimTextOutput("nadirText")),
        tabPanel("6. Toxicity",
                 plotOutput("toxPlot",     height = 480),
                 tableOutput("toxSummary")),
        tabPanel("7. Endpoints",
                 plotOutput("deauvillePlot", height = 240),
                 plotOutput("hazPlot",       height = 240),
                 tableOutput("endpointTab")),
        tabPanel("8. Biomarkers",
                 plotOutput("biomarkerPlot", height = 480))
      )
    )
  )
)

server <- function(input, output, session){
  sim <- eventReactive(input$run, {
    ev_obj <- build_events(input$regimen, input$cycles)
    mod_run <- mod_hl
    if (!is.null(input$GCSF) && input$GCSF){
      mod_run <- mod_run %>% param(GCSF_factor = 0.5)
    }
    mod_run %>%
      param(TumorW0 = input$baseline_MTV / 1200) %>%
      ev(ev_obj) %>%
      mrgsim(end = input$tend * 24, delta = 6) %>%
      as.data.frame()
  })

  output$riskTxt <- renderPrint({
    cat("Stage:        ", input$stage,        "\n")
    cat("Histology:    ", input$histology,    "\n")
    cat("IPS score:    ", input$IPS,          "\n")
    cat("EBV+:         ", input$EBV,          "\n")
    cat("B-symptoms:   ", input$Bsymp,        "\n")
    cat("Baseline MTV: ", input$baseline_MTV, " mL\n")
    cat("Regimen:      ", input$regimen, " x ", input$cycles, " cycles\n", sep="")
  })

  output$riskPlot <- renderPlot({
    df <- data.frame(
      Factor = c("IPS", "MTV (norm)", "Stage", "EBV", "B-symp"),
      Score  = c(input$IPS/7,
                 input$baseline_MTV/4000,
                 match(input$stage, c("IA","IIA early-favorable","IIB early-unfavorable","III","IV"))/5,
                 as.numeric(input$EBV),
                 as.numeric(input$Bsymp))
    )
    ggplot(df, aes(Factor, Score, fill=Factor)) +
      geom_col(width=0.6) + ylim(0,1) + theme_minimal() +
      theme(legend.position="none") +
      labs(y = "Normalized risk", title = "Patient risk profile")
  })

  output$pkPlot <- renderPlot({
    d <- sim()
    d_long <- d %>%
      select(time, Cdox, Cnivo_out) %>%
      mutate(BV = NA, MMAE = NA) %>%
      pivot_longer(-time, names_to="Drug", values_to="Conc")
    # add BV and MMAE from raw if available
    if ("BV_C" %in% names(d))   d_long$Conc[d_long$Drug=="BV"]   <- d$BV_C/6
    if ("MMAE_C" %in% names(d)) d_long$Conc[d_long$Drug=="MMAE"] <- d$MMAE_C/100
    ggplot(d_long, aes(time/24, Conc, color=Drug)) +
      geom_line(linewidth=0.7) + theme_bw() +
      labs(x="Time (days)", y="Concentration (mg/L)", title="Drug PK profiles")
  })

  output$tumorPlot <- renderPlot({
    d <- sim()
    ggplot(d, aes(time/24, MTV)) + geom_line(color="#d62728", linewidth=0.8) +
      theme_bw() + labs(x="Days", y="MTV (mL)", title="Metabolic tumor volume (HRS surrogate)")
  })
  output$tarcPlot <- renderPlot({
    d <- sim()
    ggplot(d, aes(time/24, TARC)) + geom_line(color="#1f77b4", linewidth=0.8) +
      theme_bw() + labs(x="Days", y="Serum TARC (pg/mL)", title="CCL17 (TARC) kinetics")
  })

  output$immunePlot <- renderPlot({
    d <- sim()
    d_long <- d %>%
      select(time, T_eff, T_exh, TREG) %>%
      pivot_longer(-time, names_to="Pool", values_to="Level")
    ggplot(d_long, aes(time/24, Level, color=Pool)) +
      geom_line(linewidth=0.8) + theme_bw() +
      labs(x="Days", y="Cell pool (AU)", title="T cell compartments")
  })
  output$pdl1Plot <- renderPlot({
    d <- sim()
    ggplot(d, aes(time/24, PDL1)) + geom_line(linewidth=0.8, color="#9467bd") +
      theme_bw() + labs(x="Days", y="PD-L1 (AU)", title="HRS PD-L1 expression")
  })

  output$ancPlot <- renderPlot({
    d <- sim()
    ggplot(d, aes(time/24, ANC)) +
      geom_line(linewidth=0.8) +
      geom_hline(yintercept = 0.5, linetype=2, color="red") +
      geom_hline(yintercept = 1.0, linetype=3, color="orange") +
      theme_bw() + labs(x="Days", y="ANC (10^9/L)", title="Neutrophil kinetics")
  })
  output$nadirText <- renderPrint({
    d <- sim()
    nadir <- min(d$ANC)
    fn    <- mean(d$ANC < 0.5)
    cat(sprintf("ANC nadir: %.2f x10^9/L\n", nadir))
    cat(sprintf("Time spent <0.5 (Grade 4): %.1f%% of horizon\n", 100*fn))
  })

  output$toxPlot <- renderPlot({
    d <- sim()
    d_long <- d %>%
      select(time, CARDIO, LUNG, NEURO, irAE) %>%
      pivot_longer(-time, names_to="Domain", values_to="Index")
    ggplot(d_long, aes(time/24, Index, color=Domain)) +
      geom_line(linewidth=0.8) + theme_bw() + ylim(0,1) +
      labs(x="Days", y="Toxicity index (0-1)", title="Cumulative toxicity")
  })
  output$toxSummary <- renderTable({
    d <- sim()
    data.frame(
      Domain = c("Cardiotoxicity", "Pulmonary (BLM)",
                 "Peripheral neuropathy", "Immune-related AE"),
      Final  = round(c(tail(d$CARDIO,1), tail(d$LUNG,1),
                       tail(d$NEURO,1), tail(d$irAE,1)), 3)
    )
  })

  output$deauvillePlot <- renderPlot({
    d <- sim()
    ggplot(d, aes(time/24, Deauville)) +
      geom_step(linewidth=0.7) + theme_bw() + ylim(1,5) +
      geom_hline(yintercept=3, linetype=2, color="red") +
      labs(x="Days", y="Deauville (1-5)", title="Surrogate Deauville score")
  })
  output$hazPlot <- renderPlot({
    d <- sim()
    ggplot(d, aes(time/24, PFS_surv)) +
      geom_line(linewidth=0.8) + theme_bw() + ylim(0,1) +
      labs(x="Days", y="PFS surrogate", title="Cumulative PFS (model)")
  })
  output$endpointTab <- renderTable({
    d <- sim()
    data.frame(
      Endpoint = c("End MTV (mL)", "End TARC (pg/mL)", "Deauville",
                   "PFS surrogate (1y)", "pCR achieved (HRS<0.02)"),
      Value    = c(round(tail(d$MTV,1),1),
                   round(tail(d$TARC,1),1),
                   tail(d$Deauville,1),
                   round(exp(-d$HAZ[which.min(abs(d$time-365*24))]),3),
                   ifelse(tail(d$pCR,1) > 0.5, "Yes", "No"))
    )
  })

  output$biomarkerPlot <- renderPlot({
    d <- sim()
    d_long <- d %>%
      select(time, TARC, IL6, IL13) %>%
      pivot_longer(-time, names_to="Marker", values_to="Value")
    ggplot(d_long, aes(time/24, Value, color=Marker)) +
      geom_line(linewidth=0.7) + theme_bw() +
      facet_wrap(~Marker, scales="free_y", ncol=1) +
      labs(x="Days", y="Concentration", title="Serum biomarker kinetics")
  })
}

shinyApp(ui, server)
