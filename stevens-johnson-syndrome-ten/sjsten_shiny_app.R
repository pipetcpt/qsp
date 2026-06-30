## =========================================================================
## SJS/TEN QSP Dashboard — Shiny app
## Tabs:
##   1) Patient profile        2) Drug PK
##   3) Immune drivers (PD)    4) Clinical endpoints (BSA, SCORTEN, mortality)
##   5) Scenario comparison    6) Biomarkers          7) References
## Run:  shiny::runApp("sjsten_shiny_app.R")
## =========================================================================

library(shiny)
library(mrgsolve)
library(dplyr)
library(tidyr)
library(ggplot2)
library(DT)

# -- Load mrgsolve model from same directory
mod <- mread("sjsten_mrgsolve_model.R")

scenario_choices <- c(
  "Supportive only"        = "Supportive",
  "IVIG 1 g/kg/d × 4"      = "IVIG",
  "Cyclosporine 3 mg/kg"   = "CSA",
  "Etanercept 25 mg SC"    = "ETAN",
  "Methylpred 1 mg/kg"     = "PRED",
  "CSA + Etanercept"       = "CSA_ETAN",
  "JAK inhibitor (tofa)"   = "JAKI"
)

apply_scenario <- function(m, key) {
  pars <- list(
    Supportive = list(),
    IVIG       = list(SCEN_IVIG=1),
    CSA        = list(SCEN_CSA=1),
    ETAN       = list(SCEN_ETAN=1),
    PRED       = list(SCEN_PRED=1),
    CSA_ETAN   = list(SCEN_CSA=1, SCEN_ETAN=1),
    JAKI       = list(SCEN_JAKI=1)
  )[[key]]
  param(m, pars)
}

build_events <- function(key, dose_drug=400) {
  e <- ev(amt=dose_drug, ii=12, addl=2, cmt="A_drug_dep")
  if (key %in% c("ETAN","CSA_ETAN"))
    e <- c(e, ev(amt=25, ii=3,   addl=1, cmt="ETAN_dep"))
  if (key %in% c("CSA","CSA_ETAN"))
    e <- c(e, ev(amt=100, ii=0.5, addl=20, cmt="CSA_dep"))
  if (key == "PRED")
    e <- c(e, ev(amt=60, ii=1, addl=5, cmt="PRED"))
  if (key == "JAKI")
    e <- c(e, ev(amt=5, ii=0.5, addl=20, cmt="JAKI_dep"))
  e
}

sim_scenario <- function(key, covs, dose_drug) {
  m <- apply_scenario(mod, key) %>% param(covs)
  e <- build_events(key, dose_drug)
  m %>% ev(e) %>% mrgsim(end=30, delta=0.1) %>% as_tibble() %>%
    mutate(scenario = key)
}

# -- UI -------------------------------------------------------------------
ui <- fluidPage(
  titlePanel("SJS/TEN QSP Dashboard — drug → CD8⁺ CTL → epidermal detachment"),
  sidebarLayout(
    sidebarPanel(
      width = 3,
      h4("Patient"),
      sliderInput("AGE", "Age (y)", 5, 90, 45),
      sliderInput("WT",  "Weight (kg)", 20, 130, 65),
      checkboxInput("HLA", "HLA risk allele present (e.g. B*15:02)", TRUE),
      h4("Culprit drug"),
      selectInput("culprit", "Drug",
                  c("Carbamazepine"="CBZ","Allopurinol"="ALL","Sulfamethoxazole"="SMX",
                    "Lamotrigine"="LTG","Phenytoin"="PHT","Nevirapine"="NVP")),
      sliderInput("dose_drug", "Drug dose (mg q12h)", 100, 1200, 400, step=50),
      checkboxInput("withdraw", "Drug withdrawn on day 0", TRUE),
      h4("Scenarios"),
      checkboxGroupInput("scens", "Run",
                         choices = scenario_choices,
                         selected = c("Supportive","CSA","ETAN","CSA_ETAN")),
      actionButton("run", "Run simulation", class="btn-primary")
    ),
    mainPanel(
      width = 9,
      tabsetPanel(
        tabPanel("1. Patient profile",
                 plotOutput("plot_drug"),
                 verbatimTextOutput("patient_summary")),
        tabPanel("2. Drug PK",
                 plotOutput("plot_drugPK"),
                 plotOutput("plot_biopk")),
        tabPanel("3. Immune drivers (PD)",
                 plotOutput("plot_Tact"),
                 plotOutput("plot_cytokines")),
        tabPanel("4. Clinical endpoints",
                 plotOutput("plot_bsa"),
                 plotOutput("plot_scorten"),
                 plotOutput("plot_mortality")),
        tabPanel("5. Scenario comparison",
                 plotOutput("plot_compare_bsa"),
                 plotOutput("plot_compare_mort"),
                 DT::dataTableOutput("tbl_outcomes")),
        tabPanel("6. Biomarkers",
                 plotOutput("plot_gnly"),
                 plotOutput("plot_fasl"),
                 plotOutput("plot_il15"),
                 plotOutput("plot_hmgb")),
        tabPanel("7. References",
                 includeMarkdown("sjsten_references.md"))
      )
    )
  )
)

# -- Server ---------------------------------------------------------------
server <- function(input, output, session) {
  results <- eventReactive(input$run, {
    req(length(input$scens) > 0)
    covs <- list(
      AGE       = input$AGE,
      WT        = input$WT,
      HLA_RISK  = as.integer(input$HLA),
      SCEN_WD   = as.integer(input$withdraw)
    )
    bind_rows(lapply(input$scens, function(k)
      sim_scenario(k, covs, input$dose_drug)))
  }, ignoreNULL = FALSE)

  output$patient_summary <- renderText({
    paste0("Age: ", input$AGE, "  Weight: ", input$WT, " kg\n",
           "Culprit drug: ", input$culprit, "  Dose: ", input$dose_drug, " mg q12h\n",
           "HLA risk allele: ", ifelse(input$HLA,"Present","Absent"),
           "    Drug withdrawn: ", ifelse(input$withdraw,"Yes","No"),"\n",
           "Scenarios: ", paste(input$scens, collapse=", "))
  })

  pl <- function(df, y, lab) {
    ggplot(df, aes(time, .data[[y]], color=scenario)) +
      geom_line(linewidth=1) +
      labs(x="Day", y=lab) + theme_minimal()
  }

  output$plot_drug         <- renderPlot(pl(results(), "Cdrug",     "Culprit drug (mg/L)"))
  output$plot_drugPK       <- renderPlot(pl(results(), "Cdrug",     "Culprit drug PK"))
  output$plot_biopk        <- renderPlot({
    df <- results() %>%
      pivot_longer(c(IVIG, ETAN, INFL, CSA, PRED, JAKI), names_to="drug", values_to="conc")
    ggplot(df, aes(time, conc, color=drug)) + geom_line() + facet_wrap(~scenario) +
      labs(x="Day", y="Concentration") + theme_minimal()
  })
  output$plot_Tact         <- renderPlot(pl(results(), "T_act",     "Activated CTL clones"))
  output$plot_cytokines    <- renderPlot({
    df <- results() %>%
      pivot_longer(c(TNF, IFNg, IL15), names_to="cyt", values_to="val")
    ggplot(df, aes(time, val, color=scenario)) + geom_line() + facet_wrap(~cyt, scales="free_y") +
      labs(x="Day", y="pg/mL") + theme_minimal()
  })
  output$plot_bsa          <- renderPlot(pl(results(), "BSA_loss",  "BSA detachment (%)"))
  output$plot_scorten      <- renderPlot(pl(results(), "SCORTEN",   "SCORTEN"))
  output$plot_mortality    <- renderPlot(pl(results(), "PredMort",  "Predicted mortality"))
  output$plot_compare_bsa  <- renderPlot(pl(results(), "BSA_loss",  "BSA detachment (%)"))
  output$plot_compare_mort <- renderPlot(pl(results(), "PredMort",  "Predicted mortality"))
  output$plot_gnly         <- renderPlot(pl(results(), "GNLY",      "Granulysin (ng/mL)"))
  output$plot_fasl         <- renderPlot(pl(results(), "sFasL",     "sFasL (pg/mL)"))
  output$plot_il15         <- renderPlot(pl(results(), "IL15",      "IL-15 (pg/mL)"))
  output$plot_hmgb         <- renderPlot(pl(results(), "HMGB1",     "HMGB1 (ng/mL)"))

  output$tbl_outcomes <- DT::renderDataTable({
    results() %>% group_by(scenario) %>% summarise(
      "Max BSA (%)"        = round(max(BSA_loss),1),
      "Max SCORTEN"        = max(SCORTEN),
      "Day-14 mortality"   = round(approx(time, PredMort, xout=14, rule=2)$y, 3),
      "Re-epi day (50%)"   = round(approx(Re_epi, time, xout=50, rule=2)$y, 1),
      .groups="drop"
    )
  })
}

shinyApp(ui, server)
