## ============================================================================
## SIADH QSP — Shiny Dashboard
## ----------------------------------------------------------------------------
## 7 tabs: Patient profile · Drug PK · Renal V2R/AQP2 PD · Sodium & osmolality ·
##         Neuro symptoms & ODS risk · Scenario comparison · Biomarkers
## Dependencies: shiny, shinydashboard, mrgsolve, dplyr, tidyr, ggplot2, DT
## ----------------------------------------------------------------------------
suppressPackageStartupMessages({
  library(shiny)
  library(shinydashboard)
  library(mrgsolve)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(DT)
})

## ---------- Lazy model load ----------
MODEL_PATH <- "siadh_mrgsolve_model.R"
get_model <- function() {
  if (!exists(".SIADH_MOD", envir = .GlobalEnv) || is.null(.GlobalEnv$.SIADH_MOD)) {
    assign(".SIADH_MOD", mread_cache("siadh", project = "."), envir = .GlobalEnv)
  }
  .GlobalEnv$.SIADH_MOD
}

## ---------- Scenario builders ----------
build_events <- function(scenario, horizon) {
  if (scenario == "Untreated") {
    ev(amt = 0, cmt = "TOLVA_GUT")
  } else if (scenario == "Fluid restriction alone") {
    ev(amt = 0, cmt = "TOLVA_GUT")
  } else if (scenario == "Tolvaptan 15mg PO QD") {
    ev(amt = 1.0, cmt = "TOLVA_GUT", ii = 1, addl = horizon - 1)
  } else if (scenario == "Conivaptan IV load + infusion") {
    ev(amt = 1.4, cmt = "CONIV_CENT", ii = 1, addl = horizon - 1)
  } else if (scenario == "Demeclocycline 900mg/day") {
    ev(amt = 1.1, cmt = "DEMEC_GUT", ii = 1, addl = horizon - 1)
  } else if (scenario == "Hypertonic saline 3% (guideline-limited)") {
    ev(amt = 0.9, cmt = "HSAL_INF", ii = 0.5, addl = 2 * min(horizon, 3) - 1)
  } else if (scenario == "Overly rapid correction (unsafe example)") {
    ev(amt = 3.0, cmt = "HSAL_INF", ii = 0.25, addl = 4 * min(horizon, 3) - 1)
  } else {
    ev(amt = 0, cmt = "TOLVA_GUT")
  }
}

run_sim <- function(scenario, horizon, params) {
  mod <- get_model()
  ev_set <- build_events(scenario, horizon)
  par <- list(
    SEV            = params$sev,
    ECTOPIC        = ifelse(params$ectopic, 1, 0),
    CHRONIC        = ifelse(params$chronic, 1, 0),
    FLUID_RESTRICT = ifelse(scenario == "Fluid restriction alone" || params$fluid_restrict, 1, 0),
    LOOP_SALT      = ifelse(params$loop_salt, 1, 0),
    SGLT2I         = ifelse(params$sglt2i, 1, 0),
    DDAVP_CLAMP    = ifelse(params$ddavp_clamp, 1, 0),
    BASE_NA        = params$base_na
  )
  mod %>% param(par) %>% mrgsim(events = ev_set, end = horizon, delta = 0.25) %>%
    as_tibble() %>% mutate(scenario = scenario)
}

SCENARIO_LIST <- c(
  "Untreated", "Fluid restriction alone", "Tolvaptan 15mg PO QD",
  "Conivaptan IV load + infusion", "Demeclocycline 900mg/day",
  "Hypertonic saline 3% (guideline-limited)", "Overly rapid correction (unsafe example)"
)

## ---------- UI ----------
ui <- dashboardPage(
  dashboardHeader(title = "SIADH QSP"),
  dashboardSidebar(
    sidebarMenu(
      menuItem("1. Patient profile",        tabName = "patient", icon = icon("user")),
      menuItem("2. Drug PK",                tabName = "pk",      icon = icon("pills")),
      menuItem("3. Renal V2R / AQP2 PD",    tabName = "pd",      icon = icon("kidneys")),
      menuItem("4. Sodium & osmolality",    tabName = "na",      icon = icon("droplet")),
      menuItem("5. Neuro symptoms / ODS",   tabName = "neuro",   icon = icon("brain")),
      menuItem("6. Scenario comparison",    tabName = "compare", icon = icon("chart-line")),
      menuItem("7. Biomarkers",             tabName = "bio",     icon = icon("vial"))
    ),
    hr(),
    selectInput("scenario", "Scenario:", SCENARIO_LIST, selected = "Tolvaptan 15mg PO QD"),
    sliderInput("horizon", "Simulation horizon (days):", 3, 90, 14, step = 1),
    sliderInput("sev", "SIADH severity (0-1):", 0, 1, 0.6, step = 0.05),
    checkboxInput("ectopic", "Ectopic/malignancy-associated (e.g. SCLC)", FALSE),
    checkboxInput("chronic", "Chronic (>=48h, cerebral adaptation established)", TRUE),
    sliderInput("base_na", "Starting serum Na (mEq/L):", 105, 134, 122),
    checkboxInput("fluid_restrict", "+ Fluid restriction (adjunct)", FALSE),
    checkboxInput("loop_salt", "+ Loop diuretic & salt tablets", FALSE),
    checkboxInput("sglt2i", "+ SGLT2 inhibitor (emerging adjunct)", FALSE),
    checkboxInput("ddavp_clamp", "DDAVP clamp (overcorrection-prevention)", FALSE),
    actionButton("run", "Run simulation", icon = icon("play"),
                 style = "color:#fff;background:#1f6feb")
  ),
  dashboardBody(
    tabItems(
      tabItem("patient",
        fluidRow(
          box(width = 12, title = "Patient profile summary", status = "primary",
              solidHeader = TRUE,
              DTOutput("patient_table"),
              br(),
              p(strong("Note:"), "Ectopic/malignancy-associated SIADH (e.g. small-cell",
                "lung cancer) raises the effective non-osmotic AVP-drive ceiling and is",
                "typically less responsive to fluid restriction alone.")
          )
        ),
        fluidRow(
          box(width = 12, title = "About this dashboard", status = "info",
              p("This dashboard runs the mrgsolve QSP model for SIADH, linking",
                "non-osmotic/autonomous AVP drive, renal V2-receptor/aquaporin-2",
                "functional expression, water/sodium balance, cerebral osmotic",
                "adaptation, and correction-rate-dependent osmotic-demyelination",
                "syndrome (ODS) risk to a composite neuro-symptom score."),
              p("Pick a treatment scenario in the left panel, adjust patient severity",
                "and chronicity, then press ", strong("Run simulation"), " to update plots.")
          )
        )
      ),

      tabItem("pk",
        fluidRow(
          box(width = 12, title = "Drug exposure / PD-effect trajectories (selected scenario)",
              status = "primary", solidHeader = TRUE, plotOutput("pk_plot", height = 480))
        )
      ),

      tabItem("pd",
        fluidRow(
          box(width = 6, title = "AQP2 functional expression", plotOutput("aqp2_plot", height = 360)),
          box(width = 6, title = "Free-water clearance & total body water", plotOutput("fwcl_tbw_plot", height = 360))
        )
      ),

      tabItem("na",
        fluidRow(
          box(width = 6, title = "Serum sodium & osmolality", plotOutput("na_osm_plot", height = 360)),
          box(width = 6, title = "Urine osmolality & urine sodium", plotOutput("uosm_una_plot", height = 360))
        ),
        fluidRow(
          box(width = 12, title = "24h / 48h sodium correction rate vs guideline limit (8/24h, 18/48h)",
              plotOutput("correction_rate_plot", height = 360))
        )
      ),

      tabItem("neuro",
        fluidRow(
          box(width = 6, title = "Brain organic-osmolyte adaptation index", plotOutput("brainosm_plot", height = 360)),
          box(width = 6, title = "Neuro-symptom severity score", plotOutput("symptom_plot", height = 360))
        ),
        fluidRow(
          box(width = 12, title = "Cumulative osmotic-demyelination-syndrome (ODS) risk score",
              status = "warning", solidHeader = TRUE, plotOutput("ods_plot", height = 360))
        )
      ),

      tabItem("compare",
        fluidRow(
          box(width = 12, title = "Scenario comparison panel",
              status = "warning", solidHeader = TRUE,
              p("Runs all seven built-in scenarios with the current patient profile;",
                "press the button below."),
              actionButton("run_all", "Run all scenarios", icon = icon("rocket"),
                           style = "color:#fff;background:#0f5132"),
              br(), br(),
              plotOutput("compare_plot", height = 640)
          )
        )
      ),

      tabItem("bio",
        fluidRow(
          box(width = 12, title = "Endpoint summary table (Day-0, Day-2, Day-end)",
              status = "info", DTOutput("endpoint_table"))
        )
      )
    )
  )
)

## ---------- Server ----------
server <- function(input, output, session) {

  results <- reactiveVal(NULL)
  all_results <- reactiveVal(NULL)

  observeEvent(input$run, {
    showNotification("Running mrgsolve simulation…", type = "message", duration = 1)
    p <- list(
      sev = input$sev, ectopic = input$ectopic, chronic = input$chronic,
      base_na = input$base_na, fluid_restrict = input$fluid_restrict,
      loop_salt = input$loop_salt, sglt2i = input$sglt2i, ddavp_clamp = input$ddavp_clamp
    )
    results(run_sim(input$scenario, input$horizon, p))
  }, ignoreNULL = FALSE)

  observeEvent(input$run_all, {
    showNotification("Running 7 scenarios…", type = "message", duration = 1)
    p <- list(
      sev = input$sev, ectopic = input$ectopic, chronic = input$chronic,
      base_na = input$base_na, fluid_restrict = input$fluid_restrict,
      loop_salt = input$loop_salt, sglt2i = input$sglt2i, ddavp_clamp = input$ddavp_clamp
    )
    out <- lapply(SCENARIO_LIST, function(sc) run_sim(sc, input$horizon, p))
    all_results(bind_rows(out))
  })

  # --- Patient table ---
  output$patient_table <- renderDT({
    tibble(
      Field = c("SIADH severity", "Ectopic/malignancy-associated", "Chronic (>=48h)",
                "Starting serum Na (mEq/L)", "Fluid restriction", "Loop diuretic+salt",
                "SGLT2 inhibitor", "DDAVP clamp", "Scenario", "Horizon (d)"),
      Value = c(input$sev, input$ectopic, input$chronic, input$base_na,
                input$fluid_restrict, input$loop_salt, input$sglt2i,
                input$ddavp_clamp, input$scenario, input$horizon)
    ) %>% datatable(rownames = FALSE, options = list(dom = "t"))
  })

  # --- Drug PK ---
  output$pk_plot <- renderPlot({
    df <- results(); req(df)
    df %>% select(time, Tolvaptan_effect, Conivaptan_effect, Demeclocycline_effect, Urea_effect) %>%
      pivot_longer(-time) %>%
      ggplot(aes(time, value, colour = name)) +
        geom_line(linewidth = 0.8) +
        labs(x = "Time (days)", y = "Fractional drug effect (0-1)", colour = "Drug effect") +
        theme_minimal(base_size = 13)
  })

  # --- PD ---
  output$aqp2_plot <- renderPlot({
    df <- results(); req(df)
    ggplot(df, aes(time, AQP2_expression)) +
      geom_line(colour = "#2a9d8f", linewidth = 1) +
      labs(x = "Day", y = "AQP2 functional expression (a.u., 0-100)") +
      theme_minimal(base_size = 13)
  })
  output$fwcl_tbw_plot <- renderPlot({
    df <- results(); req(df)
    df %>% select(time, Free_water_clearance, Total_body_water) %>%
      pivot_longer(-time) %>%
      ggplot(aes(time, value, colour = name)) +
        geom_line(linewidth = 0.8) +
        facet_wrap(~name, scales = "free_y") +
        theme_minimal(base_size = 12) + theme(legend.position = "none")
  })

  # --- Sodium / osmolality ---
  output$na_osm_plot <- renderPlot({
    df <- results(); req(df)
    df %>% select(time, Serum_Na, Serum_osmolality) %>%
      pivot_longer(-time) %>%
      ggplot(aes(time, value, colour = name)) +
        geom_line(linewidth = 0.8) +
        facet_wrap(~name, scales = "free_y") +
        theme_minimal(base_size = 12) + theme(legend.position = "none")
  })
  output$uosm_una_plot <- renderPlot({
    df <- results(); req(df)
    df %>% select(time, Urine_osmolality, Urine_Na) %>%
      pivot_longer(-time) %>%
      ggplot(aes(time, value, colour = name)) +
        geom_line(linewidth = 0.8) +
        facet_wrap(~name, scales = "free_y") +
        theme_minimal(base_size = 12) + theme(legend.position = "none")
  })
  output$correction_rate_plot <- renderPlot({
    df <- results(); req(df)
    df %>% select(time, Na_correction_rate_24h, Na_correction_rate_48h) %>%
      pivot_longer(-time) %>%
      ggplot(aes(time, value, colour = name)) +
        geom_line(linewidth = 0.9) +
        geom_hline(yintercept = 8, linetype = "dashed", colour = "#c0392b") +
        geom_hline(yintercept = 18, linetype = "dotted", colour = "#8a6d3b") +
        labs(x = "Day", y = "Correction over window (mEq/L)", colour = NULL,
             caption = "Dashed = 8 mEq/L/24h limit · Dotted = 18 mEq/L/48h limit") +
        theme_minimal(base_size = 13)
  })

  # --- Neuro / ODS ---
  output$brainosm_plot <- renderPlot({
    df <- results(); req(df)
    ggplot(df, aes(time, Brain_osmolyte_adaptation)) +
      geom_line(colour = "#264653", linewidth = 1) +
      labs(x = "Day", y = "Brain organic-osmolyte index (a.u., lower=more adapted/depleted)") +
      theme_minimal(base_size = 13)
  })
  output$symptom_plot <- renderPlot({
    df <- results(); req(df)
    ggplot(df, aes(time, Neuro_symptom_score)) +
      geom_line(colour = "#9b1c1c", linewidth = 1) +
      scale_y_continuous(limits = c(0, 10)) +
      labs(x = "Day", y = "Neuro-symptom severity (0-10)") +
      theme_minimal(base_size = 13)
  })
  output$ods_plot <- renderPlot({
    df <- results(); req(df)
    ggplot(df, aes(time, ODS_risk_score)) +
      geom_line(colour = "#8a6d3b", linewidth = 1) +
      scale_y_continuous(limits = c(0, 100)) +
      labs(x = "Day", y = "Cumulative ODS risk score (0-100)") +
      theme_minimal(base_size = 13)
  })

  # --- Compare ---
  output$compare_plot <- renderPlot({
    df <- all_results(); req(df)
    df %>% select(time, scenario, Serum_Na, Na_correction_rate_24h, ODS_risk_score, Neuro_symptom_score) %>%
      pivot_longer(-c(time, scenario)) %>%
      ggplot(aes(time, value, colour = scenario)) +
        geom_line(linewidth = 0.6, alpha = 0.85) +
        facet_wrap(~name, scales = "free_y") +
        labs(x = "Day", y = NULL, colour = NULL) +
        theme_minimal(base_size = 12) +
        theme(legend.position = "bottom")
  })

  output$endpoint_table <- renderDT({
    df <- results(); req(df)
    tp <- c(0, min(2, max(df$time)), max(df$time))
    df %>% filter(time %in% tp) %>%
      select(Day = time, Serum_Na, Serum_osmolality, Urine_osmolality, Urine_Na,
             Na_correction_rate_24h, ODS_risk_score, Neuro_symptom_score) %>%
      mutate(across(-Day, ~ round(.x, 2))) %>%
      datatable(rownames = FALSE,
                options = list(dom = "t",
                               columnDefs = list(list(className='dt-center', targets="_all"))))
  })
}

shinyApp(ui, server)
