# =============================================================================
# Cardiogenic Shock — Shiny dashboard for the QSP / mrgsolve model
# -----------------------------------------------------------------------------
# Interactive comparison of inotrope / vasopressor / MCS strategies with
# SCAI-stage stratification, hemodynamic trajectories, MOF biomarkers,
# mortality survival curves and sensitivity sweeps.
# =============================================================================

library(shiny)
library(mrgsolve)
library(dplyr)
library(tibble)
library(tidyr)
library(ggplot2)
library(scales)
library(DT)

source("cs_mrgsolve_model.R")  # loads `build_cs()` + `run_cs_scenarios()`

CS_MOD <- tryCatch(build_cs(), error = function(e) NULL)

# ----------------------------------------------------------------- UI ---------
ui <- fluidPage(
  titlePanel("Cardiogenic Shock — QSP Dashboard"),
  tags$style(HTML(".well{background-color:#fafafa;}")),

  sidebarLayout(
    sidebarPanel(
      width = 3,
      h4("Patient profile"),
      sliderInput("LVEF0", "Baseline LVEF",
                  min = 0.10, max = 0.45, value = 0.22, step = 0.01),
      sliderInput("MAP0",  "Baseline MAP (mmHg)",
                  min = 45,   max = 80,   value = 62,   step = 1),
      sliderInput("LAC0",  "Baseline lactate (mmol/L)",
                  min = 1.5,  max = 12,   value = 4.0,  step = 0.1),
      sliderInput("CO0",   "Baseline CO (L/min)",
                  min = 1.0,  max = 3.5,  value = 2.0,  step = 0.1),
      sliderInput("PCWP0", "Baseline PCWP (mmHg)",
                  min = 12,   max = 35,   value = 24,   step = 1),
      sliderInput("BW",    "Body weight (kg)",
                  min = 50,   max = 110,  value = 75,   step = 1),

      hr(),
      h4("Drug doses"),
      sliderInput("NE_inf",   "Norepinephrine (µg/min)",  0,  0.5,  0.10, step = 0.01),
      sliderInput("EPI_inf",  "Epinephrine    (µg/min)",  0,  0.5,  0.00, step = 0.01),
      sliderInput("DOBU_inf", "Dobutamine (µg/kg/min)",   0,   20,    5,   step = 0.5),
      sliderInput("DA_inf",   "Dopamine   (µg/kg/min)",   0,   20,    0,   step = 0.5),
      sliderInput("MIL_inf",  "Milrinone  (µg/kg/min)",   0,  0.75,  0.00, step = 0.05),
      sliderInput("LEVO_inf", "Levosimendan (µg/kg/min)", 0,  0.20,  0.00, step = 0.01),
      sliderInput("VAS_inf",  "Vasopressin (U/min)",      0,  0.06,  0.00, step = 0.005),

      hr(),
      h4("Mechanical circulatory support"),
      checkboxInput("IABP_ON",    "IABP",    FALSE),
      checkboxInput("Impella_ON", "Impella CP", FALSE),
      checkboxInput("ECMO_ON",    "VA-ECMO", FALSE),

      hr(),
      sliderInput("Tend", "Simulation horizon (hours)", 24, 168, 72, step = 6),
      actionButton("go", "Run simulation", class = "btn-primary btn-block")
    ),

    mainPanel(
      tabsetPanel(
        id = "tabs",
        tabPanel("① Patient profile",
                 br(),
                 h4("SCAI staging snapshot"),
                 verbatimTextOutput("scai_text"),
                 plotOutput("scai_plot", height = 320),
                 helpText("SCAI Stage 1-A · 2-B(beginning) · 3-C(classic) · 4-D(deteriorating) · 5-E(extremis)")),
        tabPanel("② Drug PK",
                 br(),
                 plotOutput("pk_plot", height = 460),
                 helpText("Plasma concentrations of administered agents.")),
        tabPanel("③ Hemodynamic PD",
                 br(),
                 plotOutput("hd_plot", height = 520),
                 helpText("MAP / CO / SVR / PCWP / HR / SvO2 over time.")),
        tabPanel("④ Clinical endpoints",
                 br(),
                 plotOutput("ep_plot", height = 460),
                 plotOutput("survival_plot", height = 280),
                 helpText("Lactate clearance, urine output, creatinine, ALT and survival probability.")),
        tabPanel("⑤ Scenario comparison",
                 br(),
                 plotOutput("scen_panel", height = 720),
                 DT::dataTableOutput("scen_table"),
                 helpText("Pre-specified strategies (NE+Dobu vs NE+Mil vs Levosimendan vs MCS arms).")),
        tabPanel("⑥ Biomarkers / inflammation",
                 br(),
                 plotOutput("bio_plot", height = 520),
                 helpText("NO surge (vasoplegia), TNF-α, DAMP, RAAS, SNS dynamics."))
      )
    )
  )
)

# --------------------------------------------------------------- Server -------
server <- function(input, output, session) {
  simData <- eventReactive(input$go, {
    req(CS_MOD)
    pars <- list(
      LVEF0 = input$LVEF0, MAP0 = input$MAP0, LAC0 = input$LAC0,
      CO0 = input$CO0, PCWP0 = input$PCWP0, BW = input$BW,
      NE_inf = input$NE_inf, EPI_inf = input$EPI_inf,
      DOBU_inf = input$DOBU_inf, DA_inf = input$DA_inf,
      MIL_inf = input$MIL_inf, LEVO_inf = input$LEVO_inf,
      VAS_inf = input$VAS_inf,
      IABP_ON = as.integer(input$IABP_ON),
      Impella_ON = as.integer(input$Impella_ON),
      ECMO_ON = as.integer(input$ECMO_ON)
    )
    CS_MOD %>% param(pars) %>%
      mrgsim(end = input$Tend, delta = 0.2) %>% as_tibble()
  })

  scenData <- eventReactive(input$go, {
    req(CS_MOD)
    run_cs_scenarios(end = input$Tend)
  })

  output$scai_text <- renderPrint({
    d <- simData(); if (nrow(d) == 0) return(invisible())
    last <- tail(d, 1)
    cat(sprintf("Final SCAI proxy: %d   (LAC %.2f · MAP %.1f · CO %.2f)\n",
                last$SCAI, last$LAC, last$MAP, last$COx))
    cat(sprintf("Estimated survival probability (Cox-like): %.1f%%\n",
                100*last$SurvP))
  })
  output$scai_plot <- renderPlot({
    d <- simData(); req(nrow(d) > 0)
    ggplot(d, aes(time, SCAI)) +
      geom_step(linewidth = 1.0, colour = "#B71C1C") +
      scale_y_continuous(breaks = 1:5, labels = c("A", "B", "C", "D", "E")) +
      labs(x = "Hours", y = "SCAI stage") + theme_bw(base_size = 13)
  })

  output$pk_plot <- renderPlot({
    d <- simData(); req(nrow(d) > 0)
    pk <- d %>% select(time, NE_amt, EPI_amt, DOBU_amt, MIL_amt, LEVO_amt, OR_amt, VAS_amt) %>%
      pivot_longer(-time, names_to = "Drug", values_to = "Amount")
    ggplot(pk, aes(time, Amount, colour = Drug)) +
      geom_line(linewidth = 0.7) +
      facet_wrap(~Drug, scales = "free_y") +
      labs(x = "Hours", y = "Plasma amount (µg or U)") + theme_bw(base_size = 12)
  })

  output$hd_plot <- renderPlot({
    d <- simData(); req(nrow(d) > 0)
    hd <- d %>% select(time, MAP, COx, SVR, PCWP, HR, SVO2) %>%
      pivot_longer(-time, names_to = "var", values_to = "val")
    ggplot(hd, aes(time, val, colour = var)) +
      geom_line(linewidth = 0.7) +
      facet_wrap(~var, scales = "free_y") +
      labs(x = "Hours", y = NULL) + theme_bw(base_size = 12) +
      theme(legend.position = "none")
  })

  output$ep_plot <- renderPlot({
    d <- simData(); req(nrow(d) > 0)
    ep <- d %>% select(time, LAC, UO, CR, ALT) %>%
      pivot_longer(-time, names_to = "var", values_to = "val")
    ggplot(ep, aes(time, val, colour = var)) +
      geom_line(linewidth = 0.7) +
      facet_wrap(~var, scales = "free_y") +
      labs(x = "Hours", y = NULL) + theme_bw(base_size = 12) +
      theme(legend.position = "none")
  })
  output$survival_plot <- renderPlot({
    d <- simData(); req(nrow(d) > 0)
    ggplot(d, aes(time, SurvP)) +
      geom_line(linewidth = 0.8, colour = "#1B5E20") +
      scale_y_continuous(labels = percent_format(accuracy = 1)) +
      labs(x = "Hours", y = "Survival probability") + theme_bw(base_size = 13)
  })

  output$scen_panel <- renderPlot({
    s <- scenData(); req(nrow(s) > 0)
    long <- s %>% select(time, Scenario, MAP, COx, LAC, SurvP) %>%
      pivot_longer(c(MAP, COx, LAC, SurvP), names_to = "var", values_to = "val")
    ggplot(long, aes(time, val, colour = Scenario)) +
      geom_line(linewidth = 0.7) +
      facet_wrap(~var, scales = "free_y") +
      labs(x = "Hours", y = NULL) + theme_bw(base_size = 12)
  })
  output$scen_table <- DT::renderDataTable({
    s <- scenData(); req(nrow(s) > 0)
    s %>% group_by(Scenario) %>%
      summarise(
        MAP_end   = round(last(MAP), 1),
        CO_end    = round(last(COx), 2),
        Lac_end   = round(last(LAC), 2),
        SCAI_end  = last(SCAI),
        Surv_end  = sprintf("%.1f%%", 100*last(SurvP)),
        .groups   = "drop"
      ) %>%
      DT::datatable(options = list(dom = "t"))
  })

  output$bio_plot <- renderPlot({
    d <- simData(); req(nrow(d) > 0)
    bio <- d %>% select(time, NOex, TNFa) %>%
      pivot_longer(-time, names_to = "var", values_to = "val")
    ggplot(bio, aes(time, val, colour = var)) +
      geom_line(linewidth = 0.7) +
      facet_wrap(~var, scales = "free_y") +
      labs(x = "Hours", y = NULL) + theme_bw(base_size = 12) +
      theme(legend.position = "none")
  })
}

if (interactive()) shinyApp(ui, server)
