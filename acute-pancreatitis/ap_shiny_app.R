# =============================================================================
# acute-pancreatitis/ap_shiny_app.R
# -----------------------------------------------------------------------------
# Interactive dashboard for the Acute Pancreatitis (AP) QSP model
# 8 tabs:
#   1. Patient Profile (etiology · genetics · severity drivers)
#   2. Drug PK (8 agents)
#   3. Trypsin / Cytokine cascade
#   4. Acinar Necrosis & DAMPs
#   5. Organ Failure (SOFA sub-scores, lung/kidney/liver/hemo/CNS)
#   6. Severity & Survival (BISAP / SOFA / mortality hazard)
#   7. Scenario Comparison (10 protocols)
#   8. Biomarker Heat-map (lipase/CRP/IL-6/TG)
# Launch with:  shiny::runApp("acute-pancreatitis/ap_shiny_app.R")
# =============================================================================

suppressPackageStartupMessages({
  library(shiny)
  library(mrgsolve)
  library(dplyr)
  library(ggplot2)
  library(tidyr)
  library(DT)
})

source("ap_mrgsolve_model.R", local = TRUE)

ui <- fluidPage(
  titlePanel("Acute Pancreatitis QSP — Trypsin Activation, SIRS, MODS, Therapy"),
  sidebarLayout(
    sidebarPanel(
      width = 3,
      h4("Patient & Etiology"),
      sliderInput("WT",  "Body weight (kg)", 40, 130, 70),
      sliderInput("AGE", "Age (years)",      18, 95,  55),
      sliderInput("BMI", "BMI (kg/m²)",      16, 50,  27),
      selectInput("ETIO", "Etiology",
                  choices = c("Gallstone (1)" = 1, "Alcohol (2)" = 2,
                              "Hypertriglyceridemia (3)" = 3,
                              "Post-ERCP (4)" = 4, "Other (5)" = 5),
                  selected = 1),
      sliderInput("TG0", "Baseline TG (mg/dL)", 100, 3000, 250),
      checkboxInput("PRSS1_FLAG",  "PRSS1 R122H carrier", FALSE),
      checkboxInput("SPINK1_FLAG", "SPINK1 N34S carrier", FALSE),
      hr(),
      h4("Resuscitation & Support"),
      sliderInput("LR_RATE", "Lactated Ringer's (mL/kg/h)", 1, 12, 5),
      checkboxInput("ENteral", "Early enteral nutrition", TRUE),
      hr(),
      h4("Drug regimen (toggle to add)"),
      checkboxGroupInput("drugs", NULL,
        choices = c("Indomethacin 100 mg PR ×1 (PEP prophylaxis)" = "ind",
                    "Octreotide 100 µg SC q8h" = "oct",
                    "Gabexate 600 mg IV q6h"   = "gab",
                    "Nafamostat 1.67 mg/h IV"  = "naf",
                    "Ulinastatin 200 000 U q8h" = "uli",
                    "Meropenem 1 g IV q8h (infected necrosis)" = "mer",
                    "Anakinra 100 mg SC q24h"  = "akr",
                    "Fentanyl PCA 50 µg q1h"   = "fen"),
        selected = c("fen")),
      hr(),
      sliderInput("TEND", "Simulation horizon (h)", 24, 720, 336, step = 24),
      actionButton("run", "Simulate", class = "btn-primary"),
      width = 3
    ),
    mainPanel(
      width = 9,
      tabsetPanel(
        id = "tabs",
        tabPanel("1. Patient Profile",
                 plotOutput("plot_profile", height = "350px"),
                 DTOutput("tab_profile")),
        tabPanel("2. Drug PK",
                 plotOutput("plot_pk", height = "500px")),
        tabPanel("3. Trypsin / Cytokines",
                 plotOutput("plot_trypsin", height = "260px"),
                 plotOutput("plot_cyto",    height = "320px")),
        tabPanel("4. Necrosis & DAMPs",
                 plotOutput("plot_necrosis", height = "300px"),
                 plotOutput("plot_damps",    height = "280px")),
        tabPanel("5. Organ Failure",
                 plotOutput("plot_sofa", height = "260px"),
                 plotOutput("plot_organs", height = "360px")),
        tabPanel("6. Severity & Survival",
                 plotOutput("plot_severity", height = "260px"),
                 plotOutput("plot_survival", height = "260px"),
                 DTOutput("tab_score")),
        tabPanel("7. Scenario Comparison",
                 plotOutput("plot_scen_sofa", height = "260px"),
                 plotOutput("plot_scen_surv", height = "260px"),
                 plotOutput("plot_scen_nec",  height = "260px"),
                 DTOutput("tab_scen_summary")),
        tabPanel("8. Biomarker Heat-map",
                 plotOutput("plot_heatmap", height = "500px"))
      )
    )
  )
)

server <- function(input, output, session) {

  build_dose <- reactive({
    dose <- ev(time = 0, amt = 0, cmt = "C_FEN", ID = 1)[-1, ]
    if ("ind" %in% input$drugs) dose <- rbind(dose,
        ev(time = 0, amt = 100, cmt = "DEPOT_IND", ID = 1))
    if ("oct" %in% input$drugs) dose <- rbind(dose,
        ev(time = seq(0, input$TEND - 8, by = 8), amt = 0.1, cmt = "DEPOT_OCT", ID = 1))
    if ("gab" %in% input$drugs) dose <- rbind(dose,
        ev(time = seq(0, input$TEND - 6, by = 6), amt = 600, cmt = "C_GAB", ID = 1))
    if ("naf" %in% input$drugs) dose <- rbind(dose,
        ev(time = seq(0, input$TEND - 1, by = 1), amt = 1.67, cmt = "C_NAF", ID = 1))
    if ("uli" %in% input$drugs) dose <- rbind(dose,
        ev(time = seq(0, input$TEND - 8, by = 8), amt = 2e5, cmt = "C_ULI", ID = 1))
    if ("mer" %in% input$drugs) dose <- rbind(dose,
        ev(time = seq(0, input$TEND - 8, by = 8), amt = 1000, cmt = "C_MER", ID = 1))
    if ("akr" %in% input$drugs) dose <- rbind(dose,
        ev(time = seq(0, input$TEND - 24, by = 24), amt = 100, cmt = "DEPOT_AKR", ID = 1))
    if ("fen" %in% input$drugs) dose <- rbind(dose,
        ev(time = seq(0, input$TEND - 1, by = 1), amt = 0.05, cmt = "C_FEN", ID = 1))
    dose
  })

  sim <- eventReactive(input$run, {
    p <- list(
      WT = input$WT, AGE = input$AGE, BMI = input$BMI,
      ETIO = as.numeric(input$ETIO), TG0 = input$TG0,
      PRSS1_FLAG  = as.numeric(input$PRSS1_FLAG),
      SPINK1_FLAG = as.numeric(input$SPINK1_FLAG),
      LR_RATE = input$LR_RATE,
      ENteral = as.numeric(input$ENteral)
    )
    mod <- update(ap_model, param = p)
    d <- build_dose()
    if (nrow(d) > 0) {
      out <- mod %>% data_set(d) %>% mrgsim(end = input$TEND, delta = 1)
    } else {
      out <- mod %>% mrgsim(end = input$TEND, delta = 1)
    }
    as.data.frame(out)
  }, ignoreNULL = FALSE)

  output$plot_profile <- renderPlot({
    df <- sim()
    df %>% select(time, MAP, Cr, BIL, PF) %>%
      gather("var", "val", -time) %>%
      ggplot(aes(time/24, val, colour = var)) + geom_line(linewidth = 1) +
      facet_wrap(~var, scales = "free_y") +
      labs(x = "Days", y = "", title = "Hemodynamics & organ markers") +
      theme_minimal(base_size = 13)
  })

  output$tab_profile <- renderDT({
    df <- sim()
    f <- df %>% filter(time == max(time))
    tibble(
      Variable  = c("Final SOFA","Necrosis %","BISAP","Survival Prob","Pain VAS","CRP"),
      Value     = c(round(f$SOFA, 1), round(f$NEC, 1), round(f$BISAP, 1),
                    round(f$SURVPROB, 3), round(f$VAS, 1), round(f$CRP, 1))
    )
  }, options = list(dom = "t"))

  output$plot_pk <- renderPlot({
    df <- sim()
    df %>% select(time, IND_C, OCT_C, GAB_C, NAF_C, ULI_C, MER_C, AKR_C, FEN_C) %>%
      gather("drug", "C", -time) %>%
      ggplot(aes(time/24, C, colour = drug)) + geom_line(linewidth = 0.9) +
      facet_wrap(~drug, scales = "free_y") +
      labs(x = "Days", y = "Concentration", title = "Drug PK time courses") +
      theme_minimal(base_size = 13)
  })

  output$plot_trypsin <- renderPlot({
    df <- sim()
    df %>% select(time, TRYP) %>% ggplot(aes(time/24, TRYP)) +
      geom_line(linewidth = 1, colour = "firebrick") +
      labs(x = "Days", y = "Active trypsin (a.u.)",
           title = "Active trypsin time course") +
      theme_minimal(base_size = 13)
  })

  output$plot_cyto <- renderPlot({
    df <- sim()
    df %>% select(time, TNF, IL1, IL6, IL8, CRP) %>%
      gather("cyt", "v", -time) %>%
      ggplot(aes(time/24, v, colour = cyt)) + geom_line(linewidth = 0.9) +
      facet_wrap(~cyt, scales = "free_y") +
      labs(x = "Days", y = "Level", title = "Cytokine & acute-phase response") +
      theme_minimal(base_size = 13)
  })

  output$plot_necrosis <- renderPlot({
    df <- sim()
    ggplot(df, aes(time/24, NEC)) + geom_line(linewidth = 1, colour = "darkred") +
      labs(x = "Days", y = "% acinar necrosis", title = "Pancreatic necrosis") +
      theme_minimal(base_size = 13)
  })

  output$plot_damps <- renderPlot({
    df <- sim()
    df %>% select(time, PERM, GUT, BT) %>%
      gather("v", "val", -time) %>%
      ggplot(aes(time/24, val, colour = v)) + geom_line(linewidth = 0.9) +
      facet_wrap(~v, scales = "free_y") +
      labs(x = "Days", y = "Score", title = "Vascular leak · gut barrier · BT") +
      theme_minimal(base_size = 13)
  })

  output$plot_sofa <- renderPlot({
    df <- sim()
    ggplot(df, aes(time/24, SOFA)) + geom_line(linewidth = 1, colour = "purple") +
      labs(x = "Days", y = "SOFA score", title = "Composite SOFA score") +
      theme_minimal(base_size = 13)
  })

  output$plot_organs <- renderPlot({
    df <- sim()
    df %>% select(time, PF, Cr, BIL, MAP, GCS) %>%
      gather("v", "val", -time) %>%
      ggplot(aes(time/24, val, colour = v)) + geom_line(linewidth = 0.9) +
      facet_wrap(~v, scales = "free_y") +
      labs(x = "Days", y = "", title = "Organ sub-scores") +
      theme_minimal(base_size = 13)
  })

  output$plot_severity <- renderPlot({
    df <- sim()
    df %>% select(time, BISAP, SOFA) %>%
      gather("v", "val", -time) %>%
      ggplot(aes(time/24, val, colour = v)) + geom_line(linewidth = 1) +
      labs(x = "Days", y = "Score", title = "BISAP & SOFA") +
      theme_minimal(base_size = 13)
  })

  output$plot_survival <- renderPlot({
    df <- sim()
    ggplot(df, aes(time/24, SURVPROB)) +
      geom_line(linewidth = 1, colour = "forestgreen") +
      ylim(0, 1) +
      labs(x = "Days", y = "Survival probability",
           title = "Survival (1 - cumulative mortality hazard)") +
      theme_minimal(base_size = 13)
  })

  output$tab_score <- renderDT({
    df <- sim()
    last <- df %>% filter(time == max(time))
    tibble(
      `Final SOFA` = round(last$SOFA, 1),
      `BISAP`       = round(last$BISAP, 1),
      `Necrosis %`  = round(last$NEC, 1),
      `Surv prob`   = round(last$SURVPROB, 3),
      `CRP (mg/L)`  = round(last$CRP, 1),
      `Cr (mg/dL)`  = round(last$Cr, 2),
      `PF ratio`    = round(last$PF, 0),
      `MAP`         = round(last$MAP, 0)
    )
  }, options = list(dom = "t"))

  scen_run <- eventReactive(input$run, {
    run_all_scenarios()
  }, ignoreNULL = FALSE)

  output$plot_scen_sofa <- renderPlot({
    df <- scen_run()
    ggplot(df, aes(time/24, SOFA, colour = scenario)) +
      geom_line(linewidth = 0.9) +
      labs(x = "Days", y = "SOFA", title = "SOFA across 10 scenarios") +
      theme_minimal(base_size = 13)
  })

  output$plot_scen_surv <- renderPlot({
    df <- scen_run()
    ggplot(df, aes(time/24, SURVPROB, colour = scenario)) +
      geom_line(linewidth = 0.9) + ylim(0, 1) +
      labs(x = "Days", y = "Survival prob",
           title = "Survival across 10 scenarios") +
      theme_minimal(base_size = 13)
  })

  output$plot_scen_nec <- renderPlot({
    df <- scen_run()
    ggplot(df, aes(time/24, NEC, colour = scenario)) +
      geom_line(linewidth = 0.9) +
      labs(x = "Days", y = "% Necrosis",
           title = "Necrosis trajectory across 10 scenarios") +
      theme_minimal(base_size = 13)
  })

  output$tab_scen_summary <- renderDT({
    df <- scen_run()
    df %>% group_by(scenario) %>%
      filter(time == max(time)) %>%
      summarise(
        SOFA_final  = round(SOFA, 1),
        NEC_pct     = round(NEC, 1),
        BISAP       = round(BISAP, 1),
        Surv_prob   = round(SURVPROB, 3),
        VAS_final   = round(VAS, 1)
      ) %>% arrange(SOFA_final)
  }, options = list(pageLength = 10, dom = "t"))

  output$plot_heatmap <- renderPlot({
    df <- scen_run()
    h <- df %>% group_by(scenario) %>%
      summarise(
        TRYP_peak  = max(TRYP),
        TNF_peak   = max(TNF),
        IL6_peak   = max(IL6),
        CRP_peak   = max(CRP),
        Necrosis   = max(NEC),
        SOFA_peak  = max(SOFA),
        BT_peak    = max(BT),
        VAS_peak   = max(VAS)
      ) %>% gather("biomarker", "val", -scenario)
    ggplot(h, aes(biomarker, scenario, fill = val)) +
      geom_tile() +
      scale_fill_viridis_c(option = "plasma") +
      labs(x = NULL, y = NULL, title = "Biomarker peak heat-map (scenarios)") +
      theme_minimal(base_size = 13) +
      theme(axis.text.x = element_text(angle = 35, hjust = 1))
  })
}

shinyApp(ui, server)
