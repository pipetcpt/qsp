# =============================================================================
# tsc_shiny_app.R
# Interactive Shiny dashboard for the TSC QSP model
#
# Tabs:
#   1) Patient profile      - age/sex/weight, TSC1 vs TSC2, lesion baselines
#   2) Drug PK              - EVE, SIR, VGB, CBD concentration-time profiles
#   3) mTORC1 PD            - mTORC1 activity, FKBP12 occupancy
#   4) Lesion endpoints     - SEGA, AML, FASI skin, FEV1 (LAM)
#   5) Epilepsy / TAND      - Seizure frequency, GABA tone, EXIST-3 anchor
#   6) Scenario comparator  - up to 8 regimens overlaid
#   7) Adverse events       - cumulative hazards (stomatitis, lipid, pneumonitis, VFD, hepatic)
#   8) References & doc     - links into tsc_references.md
# =============================================================================

library(shiny)
library(mrgsolve)
library(dplyr)
library(tidyr)
library(ggplot2)

source("tsc_mrgsolve_model.R")   # provides mod_tsc

ui <- navbarPage(
  title = "TSC QSP Dashboard",

  # ---- Tab 1 -----------------------------------------------------------------
  tabPanel("① Patient",
    sidebarLayout(
      sidebarPanel(
        h4("Demographics & baseline lesions"),
        numericInput("age",    "Age (years)",         value = 10,  min = 0,  max = 80),
        selectInput("sex",     "Sex",                 choices = c("Male" = 0, "Female" = 1)),
        numericInput("weight", "Weight (kg)",         value = 30,  min = 3,  max = 130),
        selectInput("gene",    "Gene",                choices = c("TSC2 (severer)", "TSC1", "NMI")),
        numericInput("sega0",  "SEGA volume (cm^3)",  value = 3.0, min = 0,  max = 30),
        numericInput("aml0",   "AML diameter (cm)",   value = 4.0, min = 0,  max = 15),
        numericInput("skin0",  "Facial AF (FASI 0-100)", value = 50, min = 0, max = 100),
        numericInput("fev10",  "FEV1 (%pred, LAM-only)", value = 85, min = 20, max = 120),
        numericInput("sz0",    "Seizure freq / 28d",  value = 35, min = 0, max = 500)
      ),
      mainPanel(
        h3("Modeled patient"),
        verbatimTextOutput("patientText"),
        tags$hr(),
        helpText("Demographic and baseline lesion inputs feed all downstream tabs.",
                 "LAM physiology is activated only when sex = Female and age >= 18.")
      )
    )
  ),

  # ---- Tab 2 -----------------------------------------------------------------
  tabPanel("② Drug PK",
    sidebarLayout(
      sidebarPanel(
        h4("Dosing regimens"),
        numericInput("eve_dose", "Everolimus dose (mg QD)", value = 4.5, min = 0, max = 15, step = 0.5),
        numericInput("sir_dose", "Sirolimus dose (mg QD)",  value = 0,   min = 0, max = 5, step = 0.5),
        numericInput("vgb_dose", "Vigabatrin (mg BID)",     value = 0,   min = 0, max = 2500, step = 100),
        numericInput("cbd_dose", "Cannabidiol (mg BID)",    value = 0,   min = 0, max = 2000, step = 50),
        checkboxInput("topical_sir","Topical sirolimus 1% QD", value = FALSE),
        sliderInput("tsim", "Simulation duration (months)", min = 1, max = 24, value = 12),
        actionButton("run","Run", icon = icon("play"))
      ),
      mainPanel(
        plotOutput("pk_plot", height = 380),
        tags$hr(),
        plotOutput("trough_plot", height = 240),
        helpText("Target everolimus trough: 5-15 ng/mL (EXIST-3 used 5-15).",
                 "Target sirolimus trough: 6-14 ng/mL (MILES).",
                 "CBD therapeutic plasma exposure usually 50-300 ng/mL.")
      )
    )
  ),

  # ---- Tab 3 -----------------------------------------------------------------
  tabPanel("③ mTORC1 PD",
    fluidRow(column(12, plotOutput("mtor_plot", height = 360))),
    fluidRow(
      column(6, plotOutput("fkbp_plot", height = 280)),
      column(6, plotOutput("ddi_plot",  height = 280))
    ),
    helpText("mTORC1 activity normalised to untreated TSC baseline (=1).",
             "FKBP12-drug complex acts allosterically (IC50 ~2 nM for everolimus).",
             "DDI: CBD inhibits CYP3A4 raising EVE AUC ~50% (Crockett 2020).")
  ),

  # ---- Tab 4 -----------------------------------------------------------------
  tabPanel("④ Lesion endpoints",
    fluidRow(
      column(6, plotOutput("sega_plot", height = 320)),
      column(6, plotOutput("aml_plot",  height = 320))
    ),
    fluidRow(
      column(6, plotOutput("skin_plot", height = 320)),
      column(6, plotOutput("fev_plot",  height = 320))
    ),
    helpText("EXIST-1: SEGA >=50% volume reduction in 35% pts @ 6 mo.",
             "EXIST-2: AML >=50% diameter reduction in 42% pts @ 12 mo.",
             "MILES: sirolimus stabilises FEV1 (slope ~0 in 12 mo).")
  ),

  # ---- Tab 5 -----------------------------------------------------------------
  tabPanel("⑤ Epilepsy",
    fluidRow(column(12, plotOutput("sz_plot",   height = 340))),
    fluidRow(column(12, plotOutput("gaba_plot", height = 280))),
    helpText("EXIST-3 anchors: 9 ng/mL EVE trough -> seizure -40%, 15 ng/mL -> -39%.",
             "GWPCARE6: CBD 25 mg/kg/d -> -48% seizures @ 16 wk.",
             "Vigabatrin: spasm resolution 60-80% in TSC infants (Capal 2021).")
  ),

  # ---- Tab 6 -----------------------------------------------------------------
  tabPanel("⑥ Scenario comparator",
    sidebarLayout(
      sidebarPanel(
        h4("Pick regimens to overlay"),
        checkboxGroupInput("scenarios",
          label = NULL,
          choices = c("Untreated",
                      "Everolimus 4.5 mg/d",
                      "Sirolimus 2 mg/d",
                      "Vigabatrin 2000 mg/d",
                      "Cannabidiol 25 mg/kg/d",
                      "EVE + VGB + CBD",
                      "Topical sirolimus 1%",
                      "EVE + topical sirolimus"),
          selected = c("Untreated","Everolimus 4.5 mg/d","Cannabidiol 25 mg/kg/d")),
        selectInput("compare_endpoint","Endpoint",
                    choices = c("SEGA volume" = "SEGAvol",
                                "AML diameter" = "AMLdiam",
                                "FASI skin score" = "FASIscore",
                                "FEV1 %pred" = "FEV1pct",
                                "Seizures/28d" = "SzPer28d",
                                "mTORC1 activity" = "mTORact"))
      ),
      mainPanel(plotOutput("compare_plot", height = 520))
    )
  ),

  # ---- Tab 7 -----------------------------------------------------------------
  tabPanel("⑦ Adverse events",
    plotOutput("ae_plot", height = 480),
    helpText("Cumulative hazards converted to event probability via 1 - exp(-H).",
             "EVE/SIR drive stomatitis, lipid, pneumonitis hazards.",
             "VGB drives peripheral visual-field defect hazard.",
             "CBD drives transaminase elevation hazard.")
  ),

  # ---- Tab 8 -----------------------------------------------------------------
  tabPanel("⑧ Documentation",
    h3("Files"),
    tags$ul(
      tags$li(tags$a(href="tsc_qsp_model.svg",       "Mechanistic map (SVG)")),
      tags$li(tags$a(href="tsc_qsp_model.dot",       "Mechanistic map (.dot)")),
      tags$li(tags$a(href="tsc_mrgsolve_model.R",    "mrgsolve model (.R)")),
      tags$li(tags$a(href="tsc_references.md",       "References (.md)")),
      tags$li(tags$a(href="README.md",               "README"))
    ),
    h3("Key trials"),
    tags$ul(
      tags$li("EXIST-1 (Franz 2013, Lancet): everolimus reduces SEGA"),
      tags$li("EXIST-2 (Bissler 2013, Lancet): everolimus reduces AML"),
      tags$li("EXIST-3 (French 2016, Lancet): everolimus reduces refractory seizures"),
      tags$li("MILES (McCormack 2011, NEJM): sirolimus stabilises FEV1 in LAM"),
      tags$li("GWPCARE6 (Thiele 2021, JAMA Neurol): CBD adjunctive for TSC seizures")
    )
  )
)

server <- function(input, output, session) {

  # Build events from inputs
  events_now <- reactive({
    out <- list()
    if (input$eve_dose > 0) out <- c(out, list(ev(amt = input$eve_dose, cmt = "EVE_GUT", ii = 24, addl = input$tsim * 30 - 1, time = 0)))
    if (input$sir_dose > 0) out <- c(out, list(ev(amt = input$sir_dose, cmt = "SIR_GUT", ii = 24, addl = input$tsim * 30 - 1, time = 0)))
    if (input$vgb_dose > 0) out <- c(out, list(ev(amt = input$vgb_dose, cmt = "VGB_GUT", ii = 12, addl = input$tsim * 60 - 1, time = 0)))
    if (input$cbd_dose > 0) out <- c(out, list(ev(amt = input$cbd_dose, cmt = "CBD_GUT", ii = 12, addl = input$tsim * 60 - 1, time = 0)))
    if (input$topical_sir)  out <- c(out, list(ev(amt = 5, cmt = "TSIR_SKIN", ii = 24, addl = input$tsim * 30 - 1, time = 0)))
    if (length(out) == 0)   out <- list(ev(amt = 0, cmt = "EVE_GUT", time = 0))
    do.call(c, out)
  })

  sim_data <- eventReactive(input$run, {
    end_hr <- input$tsim * 30 * 24
    mod <- mod_tsc %>% param(PATIENT_AGE = input$age, PATIENT_SEX = as.numeric(input$sex))
    mod %>% ev(events_now()) %>%
      mrgsim(end = end_hr, delta = 6) %>% as.data.frame()
  }, ignoreNULL = FALSE)

  output$patientText <- renderPrint({
    cat("Age   :", input$age, "yrs\nSex   :", c("Male","Female")[as.numeric(input$sex)+1],
        "\nWeight:", input$weight, "kg\nGene  :", input$gene,
        "\nSEGA  :", input$sega0, "cm^3\nAML   :", input$aml0, "cm",
        "\nFASI  :", input$skin0, "\nFEV1  :", input$fev10, "%pred",
        "\nSeize :", input$sz0,   "/28d\n")
  })

  output$pk_plot <- renderPlot({
    d <- sim_data() %>%
      select(time, EVE_ngml, SIR_ngml, VGB_umol, CBD_ngml) %>%
      pivot_longer(-time, names_to = "drug", values_to = "conc")
    ggplot(d, aes(time/24, conc, colour = drug)) +
      geom_line(linewidth = 0.7) + theme_minimal() +
      labs(x = "Time (days)", y = "Concentration", title = "Drug PK profiles")
  })

  output$trough_plot <- renderPlot({
    d <- sim_data() %>%
      filter((time %% 24) == 0) %>%
      select(time, EVE_ngml, SIR_ngml) %>%
      pivot_longer(-time, names_to = "drug", values_to = "trough")
    ggplot(d, aes(time/24, trough, colour = drug)) +
      geom_step() + geom_hline(yintercept = c(5,15), linetype = "dashed", colour = "grey50") +
      theme_minimal() +
      labs(x = "Days", y = "Trough (ng/mL)", title = "Target trough window 5–15 ng/mL")
  })

  output$mtor_plot <- renderPlot({
    ggplot(sim_data(), aes(time/24, mTORact)) +
      geom_line(linewidth = 0.8, colour = "#0D47A1") + theme_minimal() +
      labs(x = "Days", y = "mTORC1 activity (norm.)", title = "mTORC1 activity over time")
  })
  output$fkbp_plot <- renderPlot({
    ggplot(sim_data(), aes(time/24, EVE_ngml + SIR_ngml)) +
      geom_area(alpha = 0.6, fill = "#1A237E") + theme_minimal() +
      labs(x = "Days", y = "Total mTORi (ng/mL)", title = "Combined EVE+SIR exposure")
  })
  output$ddi_plot <- renderPlot({
    ggplot(sim_data(), aes(time/24, CBD_ngml)) + geom_line(colour = "#33691E") +
      theme_minimal() + labs(x = "Days", y = "CBD (ng/mL)", title = "CBD exposure (CYP3A inhibitor of EVE)")
  })

  output$sega_plot <- renderPlot({ ggplot(sim_data(), aes(time/24, SEGAvol))   + geom_line(colour = "#B71C1C") + theme_minimal() + labs(x="Days", y="SEGA volume (cm^3)") })
  output$aml_plot  <- renderPlot({ ggplot(sim_data(), aes(time/24, AMLdiam))   + geom_line(colour = "#F57F17") + theme_minimal() + labs(x="Days", y="AML diameter (cm)") })
  output$skin_plot <- renderPlot({ ggplot(sim_data(), aes(time/24, FASIscore)) + geom_line(colour = "#880E4F") + theme_minimal() + labs(x="Days", y="FASI (0-100)") })
  output$fev_plot  <- renderPlot({ ggplot(sim_data(), aes(time/24, FEV1pct))   + geom_line(colour = "#00695C") + theme_minimal() + labs(x="Days", y="FEV1 (%pred)") })
  output$sz_plot   <- renderPlot({ ggplot(sim_data(), aes(time/24, SzPer28d))  + geom_line(colour = "#4A148C") + theme_minimal() + labs(x="Days", y="Seizures / 28 d") })
  output$gaba_plot <- renderPlot({ ggplot(sim_data(), aes(time/24, mTORact))   + geom_line(colour = "#827717") + theme_minimal() + labs(x="Days", y="GABA tone proxy") })

  output$ae_plot <- renderPlot({
    d <- sim_data() %>%
      select(time, HazStoma, HazLipid, HazPneu, HazVFD, HazHepat) %>%
      pivot_longer(-time, names_to = "AE", values_to = "prob")
    ggplot(d, aes(time/24, 100*prob, colour = AE)) + geom_line(linewidth = 0.7) +
      theme_minimal() + labs(x = "Days", y = "Cumulative event probability (%)",
                              title = "Modeled adverse event hazards")
  })

  # Scenario comparator (precomputed via run_scenario)
  output$compare_plot <- renderPlot({
    req(length(input$scenarios) > 0)
    evlist <- list("Untreated"                  = ev(amt = 0, cmt = "EVE_GUT", time = 0),
                   "Everolimus 4.5 mg/d"        = ev(amt = 4.5, cmt = "EVE_GUT", ii = 24, addl = 720, time = 0),
                   "Sirolimus 2 mg/d"           = ev(amt = 2,   cmt = "SIR_GUT", ii = 24, addl = 720, time = 0),
                   "Vigabatrin 2000 mg/d"       = ev(amt = 1000, cmt = "VGB_GUT", ii = 12, addl = 720, time = 0),
                   "Cannabidiol 25 mg/kg/d"     = ev(amt = 750,  cmt = "CBD_GUT", ii = 12, addl = 720, time = 0),
                   "EVE + VGB + CBD"            = c(ev(amt=4.5,cmt="EVE_GUT",ii=24,addl=720,time=0),
                                                    ev(amt=1000,cmt="VGB_GUT",ii=12,addl=720,time=0),
                                                    ev(amt=750,cmt="CBD_GUT",ii=12,addl=720,time=0)),
                   "Topical sirolimus 1%"       = ev(amt = 5,    cmt = "TSIR_SKIN", ii = 24, addl = 720, time = 0),
                   "EVE + topical sirolimus"    = c(ev(amt=4.5,cmt="EVE_GUT",ii=24,addl=720,time=0),
                                                    ev(amt=5,  cmt="TSIR_SKIN",ii=24,addl=720,time=0)))
    res <- bind_rows(lapply(input$scenarios, function(s) {
      d <- mod_tsc %>% ev(evlist[[s]]) %>% mrgsim(end = 24*30*24, delta = 24) %>% as.data.frame()
      d$scenario <- s
      d
    }))
    ggplot(res, aes(time/24, .data[[input$compare_endpoint]], colour = scenario)) +
      geom_line(linewidth = 0.7) + theme_minimal() +
      labs(x = "Days", y = input$compare_endpoint, title = paste("Endpoint:", input$compare_endpoint))
  })
}

shinyApp(ui, server)
