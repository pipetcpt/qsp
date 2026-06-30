# =============================================================================
# Restless Legs Syndrome (RLS) — QSP Shiny App skeleton
#   8 tabs covering: (1) Patient Profile, (2) Drug PK, (3) Brain Iron / Ferritin,
#                    (4) Dopamine·Adenosine·α2δ tone, (5) IRLS / PLMS endpoints,
#                    (6) Augmentation hazard, (7) Scenario comparison,
#                    (8) Biomarkers & QoL surrogates
#   Requires the model object `mod` from rls_mrgsolve_model.R
# =============================================================================

library(shiny)
library(mrgsolve)
library(dplyr)
library(tidyr)
library(ggplot2)
library(DT)

# Source the compiled model (must be in same dir)
source("rls_mrgsolve_model.R")

# UI --------------------------------------------------------------------------
ui <- fluidPage(
  titlePanel("Restless Legs Syndrome — QSP Model Dashboard"),

  sidebarLayout(
    sidebarPanel(
      h4("Patient profile"),
      numericInput("WT",   "Weight (kg)",   value = 70,  min = 30,  max = 150),
      numericInput("eGFR", "eGFR (mL/min)", value = 90,  min = 15,  max = 130),
      selectInput("SEX",   "Sex",
                  choices = c("Male"=0, "Female"=1)),
      selectInput("PREG",  "Pregnancy (3rd trimester)?",
                  choices = c("No"=0, "Yes"=1)),
      selectInput("ESRD",  "ESRD / hemodialysis?",
                  choices = c("No"=0, "Yes"=1)),
      selectInput("SSRI",  "SSRI / SNRI trigger?",
                  choices = c("No"=0, "Yes"=1)),
      numericInput("CAFFEINE_mg", "Caffeine intake (mg/day)",
                   value = 0, min = 0, max = 1000, step = 50),
      numericInput("Ferritin0", "Baseline serum ferritin (µg/L)",
                   value = 35, min = 5, max = 300, step = 5),
      hr(),

      h4("Treatment"),
      selectInput("drug", "Choose pharmacotherapy",
                  choices = c("Untreated",
                              "Pramipexole 0.25 mg QHS",
                              "Pramipexole 0.5 mg QHS",
                              "Ropinirole 2 mg QHS",
                              "Rotigotine patch 2 mg/24h",
                              "Rotigotine patch 3 mg/24h",
                              "Gabapentin enacarbil 600 mg",
                              "Pregabalin 300 mg QHS",
                              "Oxycodone/Naloxone PR 5/2.5 mg BID",
                              "IV FCM 1000 mg single",
                              "Pramipexole 0.5 + IV FCM combo")),
      numericInput("days", "Simulation days", value = 90, min = 14, max = 365),
      actionButton("run", "Run simulation", class = "btn-primary"),
      br(), br(),
      helpText("Each scenario assumes nightly worst-symptom phase.")
    ),

    mainPanel(
      tabsetPanel(
        id = "tabs",
        tabPanel("1. Patient profile",
                 h4("Inputs & risk summary"),
                 tableOutput("profile"),
                 plotOutput("riskBar")),

        tabPanel("2. Drug PK",
                 h4("Plasma concentration vs time"),
                 plotOutput("pkPlot"),
                 helpText("Plasma curves for the selected drug + relevant comparators.")),

        tabPanel("3. Iron / Ferritin",
                 h4("Serum ferritin and brain iron index over time"),
                 plotOutput("ironPlot"),
                 helpText("Brain iron 100 = healthy reference; RLS typically 60-80.")),

        tabPanel("4. Network tones (DA / A1 / α2δ / MOR)",
                 h4("Pharmacodynamic effects on each pathway"),
                 plotOutput("tonePlot")),

        tabPanel("5. IRLS / PLMS / Sleep",
                 h4("Symptom dynamics"),
                 plotOutput("irlsPlot"),
                 plotOutput("plmsPlot"),
                 plotOutput("slpPlot")),

        tabPanel("6. Augmentation",
                 h4("Augmentation index over time"),
                 plotOutput("augPlot"),
                 helpText("Calibrated to Allen 2014: ~7.7%/yr augmentation on pramipexole 0.5 mg.")),

        tabPanel("7. Scenario comparison",
                 h4("Side-by-side IRLS trajectories"),
                 plotOutput("comparePlot"),
                 DTOutput("compareTable")),

        tabPanel("8. Biomarkers & QoL",
                 h4("Composite endpoints"),
                 plotOutput("biomarkerPlot"),
                 helpText("CGI-I, RLS-QLI surrogate, ICD hazard, opioid AE index."))
      )
    )
  )
)

# Server ----------------------------------------------------------------------
server <- function(input, output, session) {

  scenario_events <- function(label) {
    switch(label,
      "Untreated"                            = ev(amt=0, evid=0),
      "Pramipexole 0.25 mg QHS"              = ev(amt=0.25,  cmt="GUT_pra", ii=24, addl=input$days-1),
      "Pramipexole 0.5 mg QHS"               = ev(amt=0.5,   cmt="GUT_pra", ii=24, addl=input$days-1),
      "Ropinirole 2 mg QHS"                  = ev(amt=2.0,   cmt="GUT_rop", ii=24, addl=input$days-1),
      "Rotigotine patch 2 mg/24h"            = ev(amt=2.0,   cmt="PATCH_rot", ii=24, addl=input$days-1, rate=2.0/24),
      "Rotigotine patch 3 mg/24h"            = ev(amt=3.0,   cmt="PATCH_rot", ii=24, addl=input$days-1, rate=3.0/24),
      "Gabapentin enacarbil 600 mg"          = ev(amt=600,   cmt="GUT_gab", ii=24, addl=input$days-1),
      "Pregabalin 300 mg QHS"                = ev(amt=300,   cmt="GUT_pre", ii=24, addl=input$days-1),
      "Oxycodone/Naloxone PR 5/2.5 mg BID"   = ev(amt=5,     cmt="GUT_oxy", ii=12, addl=input$days*2-1),
      "IV FCM 1000 mg single"                = ev(amt=1000,  cmt="CEN_fcm", rate=1000/0.25),
      "Pramipexole 0.5 + IV FCM combo"       = c(ev(amt=0.5, cmt="GUT_pra", ii=24, addl=input$days-1),
                                                  ev(amt=1000, cmt="CEN_fcm", rate=1000/0.25))
    )
  }

  run_sim <- eventReactive(input$run, {
    pars <- list(
      WT = input$WT, eGFR = input$eGFR,
      SEX = as.numeric(input$SEX), PREG = as.numeric(input$PREG),
      ESRD = as.numeric(input$ESRD), SSRI = as.numeric(input$SSRI),
      CAFFEINE_mg = input$CAFFEINE_mg,
      Ferritin0 = input$Ferritin0
    )
    mod %>%
      param(pars) %>%
      ev(scenario_events(input$drug)) %>%
      mrgsim(end = input$days*24, delta = 2) %>%
      as_tibble() %>%
      mutate(day = time/24, scenario = input$drug)
  })

  # ---- Tab 1
  output$profile <- renderTable({
    data.frame(
      Parameter = c("Weight (kg)", "eGFR", "Sex", "Pregnancy",
                    "ESRD", "SSRI/SNRI", "Caffeine (mg/d)",
                    "Baseline ferritin (µg/L)", "Drug"),
      Value     = c(input$WT, input$eGFR, ifelse(input$SEX=="1","F","M"),
                    ifelse(input$PREG=="1","Yes","No"),
                    ifelse(input$ESRD=="1","Yes","No"),
                    ifelse(input$SSRI=="1","Yes","No"),
                    input$CAFFEINE_mg, input$Ferritin0, input$drug)
    )
  })

  output$riskBar <- renderPlot({
    sim <- run_sim()
    end_state <- tail(sim, 1)
    risk <- data.frame(
      Domain = c("IRLS","PLMS","Augmentation","ICD hazard","Constipation"),
      Value  = c(end_state$IRLS, end_state$PLMS, end_state$AugIndex*40,
                 end_state$ICD_haz*40, end_state$Constip*40)
    )
    ggplot(risk, aes(Domain, Value, fill=Domain)) +
      geom_col() + theme_minimal() +
      labs(y = "Index (scaled)", title="End-of-treatment risk snapshot")
  })

  # ---- Tab 2
  output$pkPlot <- renderPlot({
    sim <- run_sim() %>% pivot_longer(starts_with("CP_"), names_to="drug_cp", values_to="conc")
    ggplot(sim, aes(day, conc, color=drug_cp)) +
      geom_line() + theme_minimal() +
      labs(x="Day", y="Plasma conc (drug-specific units)", title="Drug PK")
  })

  # ---- Tab 3
  output$ironPlot <- renderPlot({
    sim <- run_sim() %>%
      select(day, FERRITIN, BRAINFE) %>%
      pivot_longer(-day)
    ggplot(sim, aes(day, value, color=name)) +
      geom_line(size=1) + theme_minimal() +
      labs(x="Day", y="Iron compartment", title="Ferritin & brain iron")
  })

  # ---- Tab 4
  output$tonePlot <- renderPlot({
    sim <- run_sim() %>%
      select(day, DA_eff, a2d_eff, MOR_eff, DAtone_eff, A1tone_eff) %>%
      pivot_longer(-day)
    ggplot(sim, aes(day, value, color=name)) +
      geom_line() + theme_minimal() +
      labs(x="Day", y="Tone / efficacy (0-1)", title="DA / α2δ / MOR / A1 tones")
  })

  # ---- Tab 5
  output$irlsPlot <- renderPlot({
    ggplot(run_sim(), aes(day, IRLS)) +
      geom_line(color="#B91C1C", size=1) + ylim(0,40) +
      geom_hline(yintercept=c(10,20,30), linetype="dotted", color="grey50") +
      annotate("text", x=2, y=c(10,20,30), label=c("mild","mod","severe"), hjust=0) +
      theme_minimal() + labs(x="Day", y="IRLS (0-40)", title="IRLS over time")
  })
  output$plmsPlot <- renderPlot({
    ggplot(run_sim(), aes(day, PLMS)) +
      geom_line(color="#0F766E", size=1) +
      geom_hline(yintercept=15, linetype="dotted") +
      theme_minimal() + labs(x="Day", y="PLMS index (events/h)", title="PLMS")
  })
  output$slpPlot <- renderPlot({
    ggplot(run_sim(), aes(day, SleepEff)) +
      geom_line(color="#5B21B6", size=1) + ylim(0,100) +
      theme_minimal() + labs(x="Day", y="Sleep efficiency surrogate (%)",
                              title="Sleep efficiency")
  })

  # ---- Tab 6
  output$augPlot <- renderPlot({
    ggplot(run_sim(), aes(day, AugIndex)) +
      geom_line(color="#92400E", size=1) +
      geom_hline(yintercept=0.3, linetype="dotted") +
      theme_minimal() + labs(x="Day", y="Augmentation index",
                              title="Augmentation accumulation")
  })

  # ---- Tab 7
  output$comparePlot <- renderPlot({
    # Run multiple comparator scenarios
    scenarios <- c("Untreated", "Pramipexole 0.5 mg QHS",
                   "Rotigotine patch 2 mg/24h",
                   "Gabapentin enacarbil 600 mg",
                   "Pregabalin 300 mg QHS",
                   "Oxycodone/Naloxone PR 5/2.5 mg BID")
    bind_rows(lapply(scenarios, function(s) {
      mod %>% ev(scenario_events(s)) %>%
        mrgsim(end = input$days*24, delta = 4) %>%
        as_tibble() %>% mutate(day=time/24, scenario=s)
    })) -> df
    ggplot(df, aes(day, IRLS, color=scenario)) +
      geom_line() + ylim(0,40) + theme_minimal() +
      labs(title="IRLS trajectories — head-to-head")
  })

  output$compareTable <- renderDT({
    scenarios <- c("Untreated", "Pramipexole 0.5 mg QHS",
                   "Rotigotine patch 2 mg/24h",
                   "Gabapentin enacarbil 600 mg",
                   "Pregabalin 300 mg QHS",
                   "Oxycodone/Naloxone PR 5/2.5 mg BID")
    tbl <- bind_rows(lapply(scenarios, function(s) {
      mod %>% ev(scenario_events(s)) %>%
        mrgsim(end = input$days*24, delta = 24) %>%
        as_tibble() %>% tail(1) %>%
        transmute(scenario=s, IRLS=round(IRLS,1),
                  PLMS=round(PLMS,1), AugIndex=round(AugIndex,3),
                  SleepEff=round(SleepEff,1))
    }))
    datatable(tbl, options=list(dom="t"))
  })

  # ---- Tab 8
  output$biomarkerPlot <- renderPlot({
    sim <- run_sim() %>% select(day, CGI_I, AugIndex, ICD_haz, Constip) %>%
      pivot_longer(-day)
    ggplot(sim, aes(day, value, color=name)) +
      geom_line() + theme_minimal() +
      labs(x="Day", y="Index", title="QoL / AE surrogates")
  })
}

shinyApp(ui, server)
