## ============================================================================
## Postural Orthostatic Tachycardia Syndrome (POTS) QSP — Shiny Dashboard
## ----------------------------------------------------------------------------
## 8 tabs: Patient · Drug PK · Autonomic/Volume PD · Hemodynamics & Stand Test ·
##         Deconditioning/Cerebral flow · Symptom & QoL · Scenario comparison ·
##         References
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
MODEL_PATH <- "pots_mrgsolve_model.R"
get_model <- function() {
  if (!exists(".POTS_MOD", envir = .GlobalEnv) || is.null(.GlobalEnv$.POTS_MOD)) {
    assign(".POTS_MOD", mread_cache("pots", project = "."), envir = .GlobalEnv)
  }
  .GlobalEnv$.POTS_MOD
}

## ---------- Scenario builders ----------
build_events <- function(scenario, horizon_h) {
  n_days <- ceiling(horizon_h / 24)
  if (scenario == "Untreated natural history") {
    ev(amt = 0, cmt = "CEN_MID")
  } else if (scenario == "Salt/fluid loading + compression") {
    ev(amt = 0, cmt = "CEN_MID")   # non-pharm flags set via param()
  } else if (scenario == "Propranolol 20 mg TID") {
    ev(amt = 20, cmt = "GUT_PROP", ii = 8/24, addl = 3 * n_days)
  } else if (scenario == "Midodrine 10 mg TID") {
    ev(amt = 10, cmt = "GUT_MID", ii = 8/24, addl = 3 * n_days)
  } else if (scenario == "Fludrocortisone 0.1 mg QD + salt") {
    ev(amt = 0.1, cmt = "GUT_FLUDRO", ii = 1, addl = n_days)
  } else if (scenario == "Ivabradine 5 mg BID") {
    ev(amt = 5, cmt = "GUT_IVA", ii = 12/24, addl = 2 * n_days)
  } else if (scenario == "Pyridostigmine 60 mg TID") {
    ev(amt = 60, cmt = "GUT_PYR", ii = 8/24, addl = 3 * n_days)
  } else if (scenario == "Droxidopa 100 mg TID") {
    ev(amt = 100, cmt = "GUT_DROX", ii = 8/24, addl = 3 * n_days)
  } else if (scenario == "Combination (BB+midodrine+compression+exercise)") {
    seq(ev(amt = 10, cmt = "GUT_PROP", ii = 8/24, addl = 3 * n_days),
        ev(amt = 5,  cmt = "GUT_MID",  ii = 8/24, addl = 3 * n_days))
  } else if (scenario == "Exercise training alone") {
    ev(amt = 0, cmt = "CEN_MID")
  } else {
    ev(amt = 0, cmt = "CEN_MID")
  }
}

flags_for_scenario <- function(scenario) {
  salt <- as.integer(scenario %in% c("Salt/fluid loading + compression",
                                      "Fludrocortisone 0.1 mg QD + salt"))
  compression <- as.integer(scenario %in% c("Salt/fluid loading + compression",
                                             "Combination (BB+midodrine+compression+exercise)"))
  exercise <- as.integer(scenario %in% c("Exercise training alone",
                                          "Combination (BB+midodrine+compression+exercise)"))
  list(salt = salt, compression = compression, exercise = exercise)
}

SCENARIO_LIST <- c(
  "Untreated natural history", "Salt/fluid loading + compression",
  "Propranolol 20 mg TID", "Midodrine 10 mg TID",
  "Fludrocortisone 0.1 mg QD + salt", "Ivabradine 5 mg BID",
  "Pyridostigmine 60 mg TID", "Droxidopa 100 mg TID",
  "Combination (BB+midodrine+compression+exercise)", "Exercise training alone"
)

run_sim <- function(scenario, horizon_h, params) {
  mod <- get_model()
  ev_set <- build_events(scenario, horizon_h)
  fl <- flags_for_scenario(scenario)
  par <- list(
    AGE          = params$age,
    WT           = params$wt,
    SUBTYPE      = as.numeric(params$subtype),
    BASE_HRSUP   = params$base_hrsup,
    BASE_COMPASS = params$base_compass,
    BASE_QOL     = params$base_qol,
    SALT_LOAD    = fl$salt,
    COMPRESSION  = fl$compression,
    EXERCISE     = fl$exercise
  )
  mod %>% param(par) %>% mrgsim(events = ev_set, end = horizon_h, delta = 1) %>%
    as_tibble() %>% mutate(scenario = scenario)
}

## ---------- UI ----------
ui <- dashboardPage(
  dashboardHeader(title = "POTS QSP"),
  dashboardSidebar(
    sidebarMenu(
      menuItem("1. Patient profile",        tabName = "patient", icon = icon("user")),
      menuItem("2. Drug PK",                tabName = "pk",      icon = icon("syringe")),
      menuItem("3. Autonomic / Volume PD",  tabName = "pd",      icon = icon("heartbeat")),
      menuItem("4. Hemodynamics / Stand test", tabName = "hemo", icon = icon("chart-line")),
      menuItem("5. Deconditioning / CBF",   tabName = "cbf",     icon = icon("brain")),
      menuItem("6. Symptom & QoL",          tabName = "qol",     icon = icon("smile")),
      menuItem("7. Scenario comparison",    tabName = "compare", icon = icon("layer-group")),
      menuItem("8. References",             tabName = "refs",    icon = icon("book"))
    ),
    hr(),
    selectInput("scenario", "Scenario:", SCENARIO_LIST, selected = "Midodrine 10 mg TID"),
    sliderInput("horizon_d", "Simulation horizon (days):", 3, 84, 14, step = 1),
    sliderInput("age",       "Age (years):",               12, 60, 28),
    sliderInput("wt",        "Weight (kg):",                40,120, 65),
    selectInput("subtype",   "POTS subtype:",
                choices = c("Mixed/idiopathic"=0, "Neuropathic"=1, "Hyperadrenergic"=2,
                            "Hypovolemic"=3, "Autoimmune/post-viral"=4), selected = 0),
    sliderInput("base_hrsup",   "Baseline supine HR (bpm):",   55, 95, 72),
    sliderInput("base_compass", "Baseline COMPASS-31 score:",   0,100, 45),
    sliderInput("base_qol",     "Baseline QoL score:",          0,100, 55),
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
              p(strong("Subtype legend:"),
                "Neuropathic = partial distal sympathetic denervation/impaired ",
                "venoconstriction; Hyperadrenergic = standing NE >600 pg/mL, excess ",
                "alpha1/beta1 drive; Hypovolemic = ~10-15% reduced plasma volume with ",
                "RAAS paradox; Autoimmune/post-viral = adrenergic/muscarinic ",
                "autoantibodies or post-COVID onset (subtypes commonly overlap).")
          )
        ),
        fluidRow(
          box(width = 12, title = "About this dashboard", status = "info",
              p("This dashboard runs the mrgsolve QSP model for postural ",
                "orthostatic tachycardia syndrome (POTS). Pick a treatment ",
                "scenario and subtype in the left panel, adjust the patient ",
                "profile, then press ", strong("Run simulation"), " to update plots."),
              p("Each tab focuses on a layer of the model: PK, autonomic/plasma-",
                "volume pharmacodynamics, standing heart-rate/hemodynamics ",
                "(the 10-min active-stand diagnostic test), cerebral blood-flow ",
                "and deconditioning, symptom/quality-of-life scores, scenario ",
                "comparison, and references.")
          )
        )
      ),

      tabItem("pk",
        fluidRow(
          box(width = 12, title = "Plasma drug exposure (selected scenario)",
              status = "primary", solidHeader = TRUE, plotOutput("pk_plot", height = 480))
        )
      ),

      tabItem("pd",
        fluidRow(
          box(width = 6, title = "Relative plasma volume & TPR indices", plotOutput("pv_tpr_plot", 360)),
          box(width = 6, title = "Standing norepinephrine & baroreflex sensitivity", plotOutput("ne_brs_plot", 360))
        )
      ),

      tabItem("hemo",
        fluidRow(
          box(width = 6, title = "Supine vs. standing heart rate (bpm)", plotOutput("hr_plot", 360)),
          box(width = 6, title = "ΔHR at 10-min stand (diagnostic threshold 30/40 bpm)", plotOutput("dhr_plot", 360))
        )
      ),

      tabItem("cbf",
        fluidRow(
          box(width = 6, title = "Cerebral blood-flow index (relative)", plotOutput("cbf_plot", 360)),
          box(width = 6, title = "Deconditioning index (0-1)", plotOutput("decon_plot", 360))
        )
      ),

      tabItem("qol",
        fluidRow(
          box(width = 6, title = "COMPASS-31 autonomic-symptom score", plotOutput("compass_plot", 360)),
          box(width = 6, title = "Quality-of-life score", plotOutput("qol_plot", 360))
        )
      ),

      tabItem("compare",
        fluidRow(
          box(width = 12, title = "Scenario comparison panel",
              status = "warning", solidHeader = TRUE,
              p("Runs all ten built-in scenarios with the current patient profile;",
                "press the button below."),
              actionButton("run_all", "Run all scenarios", icon = icon("rocket"),
                           style = "color:#fff;background:#0f5132"),
              br(), br(),
              plotOutput("compare_plot", height = 600)
          )
        ),
        fluidRow(
          box(width = 12, title = "Endpoint summary table (Day-0, Day-7, Day-end)",
              status = "info", DTOutput("endpoint_table"))
        )
      ),

      tabItem("refs",
        fluidRow(
          box(width = 12, title = "Key references", status = "primary", solidHeader = TRUE,
              p("See ", code("pots_references.md"), " in this directory for the full",
                "curated list (30+ PubMed-linked references)."),
              tags$ul(
                tags$li("Sheldon RS, et al. 2015 Heart Rhythm — POTS consensus statement (diagnostic criteria)."),
                tags$li("Raj SR, et al. 2009 Circulation — low-dose propranolol RCT."),
                tags$li("Fu Q, et al. 2010/2011 Circulation / Heart Rhythm — exercise training reverses deconditioning."),
                tags$li("Raj SR, et al. 2005 Circulation — reduced blood volume, RAAS paradox in POTS."),
                tags$li("Fedorowski A. 2022 Nat Rev Cardiol — POTS mechanisms & subtypes review."),
                tags$li("Vernino S, et al. 2021 Auton Neurosci — consensus statement incl. autoimmune/post-viral POTS.")
              )
          )
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
      age = input$age, wt = input$wt, subtype = input$subtype,
      base_hrsup = input$base_hrsup, base_compass = input$base_compass,
      base_qol = input$base_qol
    )
    results(run_sim(input$scenario, input$horizon_d * 24, p))
  }, ignoreNULL = FALSE)

  observeEvent(input$run_all, {
    showNotification("Running 10 scenarios…", type = "message", duration = 1)
    p <- list(
      age = input$age, wt = input$wt, subtype = input$subtype,
      base_hrsup = input$base_hrsup, base_compass = input$base_compass,
      base_qol = input$base_qol
    )
    out <- lapply(SCENARIO_LIST, function(sc) run_sim(sc, input$horizon_d * 24, p))
    all_results(bind_rows(out))
  })

  # --- Patient table ---
  output$patient_table <- renderDT({
    subtype_lbl <- c("0"="Mixed/idiopathic","1"="Neuropathic","2"="Hyperadrenergic",
                      "3"="Hypovolemic","4"="Autoimmune/post-viral")[[as.character(input$subtype)]]
    tibble(
      Field = c("Age", "Weight", "Subtype", "Baseline supine HR",
                "Baseline COMPASS-31", "Baseline QoL", "Scenario", "Horizon (d)"),
      Value = c(input$age, input$wt, subtype_lbl, input$base_hrsup,
                input$base_compass, input$base_qol, input$scenario, input$horizon_d)
    ) %>% datatable(rownames = FALSE, options = list(dom = "t"))
  })

  # --- PK ---
  output$pk_plot <- renderPlot({
    df <- results(); req(df)
    df %>% mutate(day = time/24) %>%
      select(day, conc_midodrine, conc_desglymidodrine, conc_propranolol,
             conc_fludrocortisone, conc_ivabradine, conc_pyridostigmine, conc_droxidopa) %>%
      pivot_longer(-day) %>%
      ggplot(aes(day, value, colour = name)) +
        geom_line(linewidth = 0.7) +
        scale_y_continuous(trans = "log1p") +
        labs(x = "Time (days)", y = "Plasma conc (mg/L)", colour = "Drug") +
        theme_minimal(base_size = 13)
  })

  # --- Autonomic / Volume PD ---
  output$pv_tpr_plot <- renderPlot({
    df <- results(); req(df)
    df %>% mutate(day = time/24) %>%
      select(day, PlasmaVolume_idx, TPR_idx) %>%
      pivot_longer(-day) %>%
      ggplot(aes(day, value, colour = name)) +
        geom_line(linewidth = 0.8) +
        geom_hline(yintercept = 1.0, lty = 2, colour = "grey50") +
        facet_wrap(~name, scales = "free_y") +
        theme_minimal(base_size = 12) + theme(legend.position = "none")
  })
  output$ne_brs_plot <- renderPlot({
    df <- results(); req(df)
    df %>% mutate(day = time/24) %>%
      select(day, StandingNE_pgmL, Baroreflex_idx) %>%
      pivot_longer(-day) %>%
      ggplot(aes(day, value, colour = name)) +
        geom_line(linewidth = 0.8) +
        facet_wrap(~name, scales = "free_y") +
        theme_minimal(base_size = 12) + theme(legend.position = "none")
  })

  # --- Hemodynamics ---
  output$hr_plot <- renderPlot({
    df <- results(); req(df)
    df %>% mutate(day = time/24) %>%
      select(day, HR_supine_bpm, HR_standing_bpm) %>%
      pivot_longer(-day) %>%
      ggplot(aes(day, value, colour = name)) +
        geom_line(linewidth = 1) +
        labs(x = "Day", y = "Heart rate (bpm)", colour = "") +
        theme_minimal(base_size = 13)
  })
  output$dhr_plot <- renderPlot({
    df <- results(); req(df)
    ggplot(df, aes(time/24, DeltaHR_bpm)) +
      geom_line(colour = "#cc3344", linewidth = 1) +
      geom_hline(yintercept = 30, lty = 2, colour = "grey50") +
      labs(x = "Day", y = "ΔHR (standing - supine, bpm)",
           caption = "Dashed line = adult POTS diagnostic threshold (≥30 bpm)") +
      theme_minimal(base_size = 13)
  })

  # --- Deconditioning / CBF ---
  output$cbf_plot <- renderPlot({
    df <- results(); req(df)
    ggplot(df, aes(time/24, CBF_idx)) +
      geom_line(colour = "#264653", linewidth = 1) +
      geom_hline(yintercept = 1.0, lty = 2, colour = "grey50") +
      labs(x = "Day", y = "Relative cerebral blood-flow index") +
      theme_minimal(base_size = 13)
  })
  output$decon_plot <- renderPlot({
    df <- results(); req(df)
    ggplot(df, aes(time/24, Deconditioning_idx)) +
      geom_line(colour = "#6d4c41", linewidth = 1) +
      scale_y_continuous(limits = c(0, 1)) +
      labs(x = "Day", y = "Deconditioning index (0-1)") +
      theme_minimal(base_size = 13)
  })

  # --- Symptom / QoL ---
  output$compass_plot <- renderPlot({
    df <- results(); req(df)
    ggplot(df, aes(time/24, COMPASS31_score)) +
      geom_line(colour = "#9b1c1c", linewidth = 1) +
      scale_y_continuous(limits = c(0, 100)) +
      labs(x = "Day", y = "COMPASS-31 score (0-100)") +
      theme_minimal(base_size = 13)
  })
  output$qol_plot <- renderPlot({
    df <- results(); req(df)
    ggplot(df, aes(time/24, QoL_score)) +
      geom_line(colour = "#2a9d8f", linewidth = 1) +
      scale_y_continuous(limits = c(0, 100)) +
      labs(x = "Day", y = "Quality-of-life score (0-100, higher=better)") +
      theme_minimal(base_size = 13)
  })

  # --- Compare ---
  output$compare_plot <- renderPlot({
    df <- all_results(); req(df)
    df %>% mutate(day = time/24) %>%
      ggplot(aes(day, DeltaHR_bpm, colour = scenario)) +
        geom_line(linewidth = 0.8) +
        geom_hline(yintercept = 30, lty = 2, colour = "grey50") +
        labs(x = "Day", y = "ΔHR (bpm)", colour = "Scenario") +
        theme_minimal(base_size = 13) +
        theme(legend.position = "bottom")
  })
  output$endpoint_table <- renderDT({
    df <- all_results(); req(df)
    horizon <- max(df$time)
    day7 <- min(7*24, horizon)
    df %>% filter(time %in% c(0, day7, horizon)) %>%
      mutate(Day = round(time/24)) %>%
      select(scenario, Day, DeltaHR_bpm, COMPASS31_score, QoL_score) %>%
      distinct(scenario, Day, .keep_all = TRUE) %>%
      arrange(scenario, Day) %>%
      datatable(rownames = FALSE, options = list(pageLength = 30))
  })
}

shinyApp(ui, server)
