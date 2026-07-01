# =============================================================================
# prurigo-nodularis/pn_shiny_app.R
# -----------------------------------------------------------------------------
# Interactive dashboard for the Prurigo Nodularis (PN) QSP model
# 8 tabs:
#   1. Patient Profile (subtype flags · baseline severity)
#   2. Drug PK (dupilumab · nemolizumab · gabapentin · nalbuphine ER · abrocitinib)
#   3. Barrier / Alarmin / Th2-Th17 Cascade
#   4. Nerve Sensitization (peripheral · central · IENFD · opioid tone)
#   5. Itch-Scratch-Fibrosis Cycle (scratch behavior · nodule burden)
#   6. Clinical Endpoints (WI-NRS · PN-IGA · sleep · DLQI)
#   7. Scenario Comparison (8 protocols)
#   8. Biomarker Summary Table
# Launch with:  shiny::runApp("prurigo-nodularis/pn_shiny_app.R")
# =============================================================================

suppressPackageStartupMessages({
  library(shiny)
  library(mrgsolve)
  library(dplyr)
  library(ggplot2)
  library(tidyr)
  library(DT)
})

source("pn_mrgsolve_model.R", local = TRUE)

ui <- fluidPage(
  titlePanel("Prurigo Nodularis QSP — Neuroimmune Itch-Scratch Cycle & Therapy"),
  sidebarLayout(
    sidebarPanel(
      width = 3,
      h4("Patient Profile"),
      sliderInput("WT",  "Body weight (kg)", 45, 130, 75),
      sliderInput("AGE", "Age (years)",      18, 90,  58),
      sliderInput("BASELINE_WINRS", "Baseline Worst-Itch NRS", 4, 10, 8.5, step = 0.5),
      sliderInput("NODULE0", "Baseline active nodule count", 5, 60, 25),
      sliderInput("DURATION_YR", "Disease duration (years)", 0.5, 15, 4, step = 0.5),
      checkboxInput("ATOPIC_FLAG", "Atopic diathesis subtype", FALSE),
      checkboxInput("CKD_FLAG", "CKD-associated pruritus overlap", FALSE),
      hr(),
      h4("Drug regimen (toggle to add)"),
      checkboxGroupInput("drugs", NULL,
        choices = c("Dupilumab 600mg LD / 300mg Q2W SC" = "dupi",
                    "Nemolizumab 60mg LD / 30mg Q4W SC"  = "nemo",
                    "Gabapentin 300mg TID oral"          = "gaba",
                    "Nalbuphine ER 162mg BID oral"        = "nal",
                    "Abrocitinib 200mg QD oral (off-label)" = "jaki"),
        selected = character(0)),
      checkboxInput("TCS_ON", "Topical corticosteroid (barrier repair)", TRUE),
      checkboxInput("KOR_AGONIST_ON", "Adjunct peripheral KOR agonist (CKD-aP)", FALSE),
      hr(),
      sliderInput("TEND_WK", "Simulation horizon (weeks)", 4, 52, 24, step = 4),
      actionButton("run", "Simulate", class = "btn-primary")
    ),
    mainPanel(
      width = 9,
      tabsetPanel(
        id = "tabs",
        tabPanel("1. Patient Profile",
                 plotOutput("plot_profile", height = "320px"),
                 DTOutput("tab_profile")),
        tabPanel("2. Drug PK",
                 plotOutput("plot_pk", height = "500px")),
        tabPanel("3. Barrier / Alarmin / Th2-Th17",
                 plotOutput("plot_barrier", height = "260px"),
                 plotOutput("plot_cyto",    height = "320px")),
        tabPanel("4. Nerve Sensitization",
                 plotOutput("plot_periph", height = "300px"),
                 plotOutput("plot_central", height = "280px")),
        tabPanel("5. Itch-Scratch-Fibrosis Cycle",
                 plotOutput("plot_scratch", height = "280px"),
                 plotOutput("plot_nodule",  height = "280px")),
        tabPanel("6. Clinical Endpoints",
                 plotOutput("plot_clin", height = "420px"),
                 DTOutput("tab_clin")),
        tabPanel("7. Scenario Comparison",
                 plotOutput("plot_scen_winrs", height = "260px"),
                 plotOutput("plot_scen_iga",   height = "260px"),
                 plotOutput("plot_scen_nod",   height = "260px"),
                 DTOutput("tab_scen_summary")),
        tabPanel("8. Biomarker Summary",
                 plotOutput("plot_heatmap", height = "500px"))
      )
    )
  )
)

server <- function(input, output, session) {

  build_dose <- reactive({
    end_h <- input$TEND_WK * 168
    dose <- data.frame(time = numeric(0), amt = numeric(0), cmt = character(0),
                        evid = numeric(0), ID = numeric(0))
    if ("dupi" %in% input$drugs) {
      dose <- rbind(dose, data.frame(time = 0, amt = 600, cmt = "DEPOT_DUPI", evid = 1, ID = 1))
      if (end_h > 336) dose <- rbind(dose, data.frame(
        time = seq(336, end_h, by = 336), amt = 300, cmt = "DEPOT_DUPI", evid = 1, ID = 1))
    }
    if ("nemo" %in% input$drugs) {
      dose <- rbind(dose, data.frame(time = 0, amt = 60, cmt = "DEPOT_NEMO", evid = 1, ID = 1))
      if (end_h > 672) dose <- rbind(dose, data.frame(
        time = seq(672, end_h, by = 672), amt = 30, cmt = "DEPOT_NEMO", evid = 1, ID = 1))
    }
    if ("gaba" %in% input$drugs) dose <- rbind(dose, data.frame(
      time = seq(0, max(end_h - 8, 0), by = 8), amt = 300, cmt = "GUT_GABA", evid = 1, ID = 1))
    if ("nal" %in% input$drugs) dose <- rbind(dose, data.frame(
      time = seq(0, max(end_h - 12, 0), by = 12), amt = 162, cmt = "GUT_NAL", evid = 1, ID = 1))
    if ("jaki" %in% input$drugs) dose <- rbind(dose, data.frame(
      time = seq(0, max(end_h - 24, 0), by = 24), amt = 200, cmt = "GUT_JAKI", evid = 1, ID = 1))
    dose[order(dose$time), ]
  })

  sim <- eventReactive(input$run, {
    p <- list(
      WT = input$WT, AGE = input$AGE,
      BASELINE_WINRS = input$BASELINE_WINRS, NODULE0 = input$NODULE0,
      DURATION_YR = input$DURATION_YR,
      ATOPIC_FLAG = as.numeric(input$ATOPIC_FLAG),
      CKD_FLAG = as.numeric(input$CKD_FLAG),
      TCS_REPAIR = if (input$TCS_ON) 0.15 else 0,
      JAKI_ON = as.numeric("jaki" %in% input$drugs),
      KOR_AGONIST_ON = as.numeric(input$KOR_AGONIST_ON)
    )
    mod <- update(pn_model, param = p)
    d <- build_dose()
    end_h <- input$TEND_WK * 168
    if (nrow(d) > 0) {
      out <- mod %>% data_set(d) %>% mrgsim(end = end_h, delta = 4)
    } else {
      out <- mod %>% mrgsim(end = end_h, delta = 4)
    }
    as.data.frame(out)
  }, ignoreNULL = FALSE)

  output$plot_profile <- renderPlot({
    df <- sim()
    df %>% select(time, WINRS, NODULE) %>%
      gather("var", "val", -time) %>%
      ggplot(aes(time/168, val, colour = var)) + geom_line(linewidth = 1) +
      facet_wrap(~var, scales = "free_y") +
      labs(x = "Weeks", y = "", title = "Baseline itch & nodule trajectory") +
      theme_minimal(base_size = 13)
  })

  output$tab_profile <- renderDT({
    df <- sim(); f <- df %>% filter(time == max(time)) %>% slice(1)
    tibble(
      Variable = c("Final WI-NRS", "Final PN-IGA", "Final Nodule Count",
                   "Final DLQI", "Responder (40% + IGA<=1)"),
      Value = c(round(f$WINRS, 1), round(f$IGA, 1), round(f$NODULE, 1),
                round(f$DLQI, 1), ifelse(f$RESPONDER > 0.5, "Yes", "No"))
    )
  }, options = list(dom = "t"))

  output$plot_pk <- renderPlot({
    df <- sim()
    df %>% select(time, CP_DUPI, CP_NEMO, CP_GABA, CP_NAL, CP_JAKI) %>%
      gather("drug", "C", -time) %>%
      ggplot(aes(time/168, C, colour = drug)) + geom_line(linewidth = 0.9) +
      facet_wrap(~drug, scales = "free_y") +
      labs(x = "Weeks", y = "Concentration", title = "Drug PK time courses") +
      theme_minimal(base_size = 13)
  })

  output$plot_barrier <- renderPlot({
    df <- sim()
    ggplot(df, aes(time/168, BARRIER)) + geom_line(linewidth = 1, colour = "chocolate") +
      labs(x = "Weeks", y = "Barrier dysfunction (0-1)", title = "Epidermal barrier dysfunction") +
      theme_minimal(base_size = 13)
  })

  output$plot_cyto <- renderPlot({
    df <- sim()
    df %>% select(time, TSLP, TH2, TH17, IL31) %>%
      gather("cyt", "v", -time) %>%
      ggplot(aes(time/168, v, colour = cyt)) + geom_line(linewidth = 0.9) +
      facet_wrap(~cyt, scales = "free_y") +
      labs(x = "Weeks", y = "Level (a.u.)", title = "Alarmin / Th2 / Th17 / IL-31 cascade") +
      theme_minimal(base_size = 13)
  })

  output$plot_periph <- renderPlot({
    df <- sim()
    df %>% select(time, PSENS, IENFD) %>%
      gather("v", "val", -time) %>%
      ggplot(aes(time/168, val, colour = v)) + geom_line(linewidth = 0.9) +
      facet_wrap(~v, scales = "free_y") +
      labs(x = "Weeks", y = "", title = "Peripheral sensitization & nerve fiber density") +
      theme_minimal(base_size = 13)
  })

  output$plot_central <- renderPlot({
    df <- sim()
    df %>% select(time, CSENS, OPIOID) %>%
      gather("v", "val", -time) %>%
      ggplot(aes(time/168, val, colour = v)) + geom_line(linewidth = 0.9) +
      facet_wrap(~v, scales = "free_y") +
      labs(x = "Weeks", y = "", title = "Central sensitization & opioid (MOR:KOR) tone") +
      theme_minimal(base_size = 13)
  })

  output$plot_scratch <- renderPlot({
    df <- sim()
    ggplot(df, aes(time/168, SCRATCH)) + geom_line(linewidth = 1, colour = "firebrick") +
      labs(x = "Weeks", y = "Scratch intensity (0-10)", title = "Scratch behavior") +
      theme_minimal(base_size = 13)
  })

  output$plot_nodule <- renderPlot({
    df <- sim()
    ggplot(df, aes(time/168, NODULE)) + geom_line(linewidth = 1, colour = "darkred") +
      labs(x = "Weeks", y = "Active nodule count", title = "Fibrotic nodule burden") +
      theme_minimal(base_size = 13)
  })

  output$plot_clin <- renderPlot({
    df <- sim()
    df %>% select(time, WINRS, IGA, SLEEP, DLQI) %>%
      gather("v", "val", -time) %>%
      ggplot(aes(time/168, val, colour = v)) + geom_line(linewidth = 0.9) +
      facet_wrap(~v, scales = "free_y") +
      labs(x = "Weeks", y = "", title = "Clinical endpoints") +
      theme_minimal(base_size = 13)
  })

  output$tab_clin <- renderDT({
    df <- sim(); last <- df %>% filter(time == max(time)) %>% slice(1)
    tibble(
      `WI-NRS`  = round(last$WINRS, 1),
      `% improve` = round(last$WINRS_PCT_IMPROVE, 0),
      `PN-IGA`  = round(last$IGA, 1),
      `Nodules` = round(last$NODULE, 1),
      `Sleep`   = round(last$SLEEP, 1),
      `DLQI`    = round(last$DLQI, 1)
    )
  }, options = list(dom = "t"))

  scen_run <- eventReactive(input$run, {
    end_h <- input$TEND_WK * 168
    lapply(scenarios, function(s) {
      run_scenario(s$label, dosing = s$dosing, params = s$params, end_t = end_h)
    }) %>% bind_rows()
  }, ignoreNULL = FALSE)

  output$plot_scen_winrs <- renderPlot({
    scen_run() %>% ggplot(aes(time/168, WINRS, colour = scenario)) +
      geom_line(linewidth = 0.9) +
      labs(x = "Weeks", y = "Worst-Itch NRS", title = "Scenario comparison — WI-NRS") +
      theme_minimal(base_size = 11) + theme(legend.position = "bottom")
  })

  output$plot_scen_iga <- renderPlot({
    scen_run() %>% ggplot(aes(time/168, IGA, colour = scenario)) +
      geom_line(linewidth = 0.9) +
      labs(x = "Weeks", y = "PN-IGA", title = "Scenario comparison — PN-IGA") +
      theme_minimal(base_size = 11) + theme(legend.position = "none")
  })

  output$plot_scen_nod <- renderPlot({
    scen_run() %>% ggplot(aes(time/168, NODULE, colour = scenario)) +
      geom_line(linewidth = 0.9) +
      labs(x = "Weeks", y = "Nodule count", title = "Scenario comparison — Nodule burden") +
      theme_minimal(base_size = 11) + theme(legend.position = "none")
  })

  output$tab_scen_summary <- renderDT({
    scen_run() %>% group_by(scenario) %>% filter(time == max(time)) %>%
      slice(1) %>% ungroup() %>%
      transmute(Scenario = scenario,
                `Final WI-NRS` = round(WINRS, 1),
                `% Improve` = round(WINRS_PCT_IMPROVE, 0),
                `Final IGA` = round(IGA, 1),
                `Final Nodules` = round(NODULE, 1),
                `Responder` = ifelse(RESPONDER > 0.5, "Yes", "No"))
  }, options = list(dom = "t", pageLength = 10))

  output$plot_heatmap <- renderPlot({
    df <- sim() %>% filter(time %% 168 < 4) %>%
      select(time, BARRIER, TSLP, TH2, TH17, IL31, PSENS, CSENS, OPIOID,
             SCRATCH, NODULE, WINRS, IGA, SLEEP, DLQI) %>%
      mutate(week = round(time/168)) %>% select(-time) %>%
      gather("biomarker", "value", -week) %>%
      group_by(biomarker) %>% mutate(value_norm = value / max(value, 1e-6)) %>% ungroup()
    ggplot(df, aes(x = week, y = biomarker, fill = value_norm)) +
      geom_tile() + scale_fill_viridis_c(name = "Normalized\nlevel") +
      labs(x = "Weeks", y = "", title = "Biomarker heat-map (normalized to peak)") +
      theme_minimal(base_size = 13)
  })
}

shinyApp(ui, server)
