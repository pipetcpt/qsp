## ============================================================================
## Opioid Use Disorder (OUD) QSP — Shiny Dashboard
## ----------------------------------------------------------------------------
## 8 tabs: Patient · Drug PK · MOR Occupancy/Tolerance PD · Withdrawal & Craving ·
##         Respiratory/Overdose Risk · Clinical Endpoints (QTc/Retention) ·
##         Scenario comparison · References
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
MODEL_PATH <- "oud_mrgsolve_model.R"
get_model <- function() {
  if (!exists(".OUD_MOD", envir = .GlobalEnv) || is.null(.GlobalEnv$.OUD_MOD)) {
    assign(".OUD_MOD", mread_cache("oud", project = "."), envir = .GlobalEnv)
  }
  .GlobalEnv$.OUD_MOD
}

## ---------- Scenario builders ----------
SCENARIO_LIST <- c(
  "Untreated chronic fentanyl use",
  "Fentanyl overdose + naloxone rescue",
  "Methadone maintenance (100 mg QD)",
  "Buprenorphine/naloxone maintenance (16 mg QD)",
  "Precipitated withdrawal (early buprenorphine induction)",
  "Naltrexone XR-IM depot (380 mg, post-detox)",
  "Untreated cessation + lofexidine",
  "Fentanyl + benzodiazepine co-use",
  "High-dose methadone (140 mg QD, QTc risk)",
  "Relapse after buprenorphine discontinuation"
)

build_events <- function(scenario, horizon_h) {
  nd <- ceiling(horizon_h / 24)
  if (scenario == "Untreated chronic fentanyl use") {
    ev(amt = 1.5, cmt = "GUT_FENT", ii = 5, addl = ceiling(horizon_h/5))
  } else if (scenario == "Fentanyl overdose + naloxone rescue") {
    seq(ev(amt = 3.0, cmt = "CEN_FENT", time = 0),
        ev(amt = 4.0, cmt = "GUT_NLX",  time = 5/60))
  } else if (scenario == "Methadone maintenance (100 mg QD)") {
    ev(amt = 100, cmt = "GUT_METH", ii = 24, addl = nd)
  } else if (scenario == "Buprenorphine/naloxone maintenance (16 mg QD)") {
    ev(amt = 16, cmt = "GUT_BUP", ii = 24, addl = nd)
  } else if (scenario == "Precipitated withdrawal (early buprenorphine induction)") {
    seq(ev(amt = 1.5, cmt = "GUT_FENT", ii = 5, addl = 6),
        ev(amt = 8, cmt = "GUT_BUP", time = 2))
  } else if (scenario == "Naltrexone XR-IM depot (380 mg, post-detox)") {
    ev(amt = 380, cmt = "DEPOT_NTX_XR", time = 0)
  } else if (scenario == "Untreated cessation + lofexidine") {
    ev(amt = 0.6, cmt = "GUT_LOF", ii = 8, addl = 3*nd)
  } else if (scenario == "Fentanyl + benzodiazepine co-use") {
    ev(amt = 1.5, cmt = "GUT_FENT", ii = 5, addl = ceiling(horizon_h/5))
  } else if (scenario == "High-dose methadone (140 mg QD, QTc risk)") {
    ev(amt = 140, cmt = "GUT_METH", ii = 24, addl = nd)
  } else if (scenario == "Relapse after buprenorphine discontinuation") {
    seq(ev(amt = 16, cmt = "GUT_BUP", ii = 24, addl = 30),
        ev(amt = 1.5, cmt = "GUT_FENT", time = 35*24))
  } else {
    ev(amt = 0, cmt = "CEN_FENT")
  }
}

flags_for_scenario <- function(scenario) {
  list(
    benzo    = as.integer(scenario == "Fentanyl + benzodiazepine co-use"),
    xylazine = 0,
    counseling = 0
  )
}

run_sim <- function(scenario, horizon_h, params) {
  mod <- get_model()
  ev_set <- build_events(scenario, horizon_h)
  fl <- flags_for_scenario(scenario)
  par <- list(
    AGE          = params$age,
    WT           = params$wt,
    YEARS_USE    = params$years_use,
    BASE_TOL     = params$base_tol,
    BASE_COWS    = params$base_cows,
    BASE_CRAVE   = params$base_crave,
    BASE_RETAIN  = params$base_retain,
    BENZO_COUSE  = fl$benzo,
    XYLAZINE_ADULT = fl$xylazine,
    COUNSELING   = as.integer(params$counseling)
  )
  mod %>% param(par) %>% mrgsim(events = ev_set, end = horizon_h, delta = 0.1) %>%
    as_tibble() %>% mutate(scenario = scenario)
}

## ---------- UI ----------
ui <- dashboardPage(
  dashboardHeader(title = "OUD QSP"),
  dashboardSidebar(
    sidebarMenu(
      menuItem("1. Patient profile",          tabName = "patient", icon = icon("user")),
      menuItem("2. Drug PK",                  tabName = "pk",      icon = icon("syringe")),
      menuItem("3. MOR occupancy/tolerance",  tabName = "pd",      icon = icon("dna")),
      menuItem("4. Withdrawal & craving",     tabName = "wd",      icon = icon("bolt")),
      menuItem("5. Respiratory/overdose risk",tabName = "resp",    icon = icon("lungs")),
      menuItem("6. Clinical endpoints",       tabName = "clin",    icon = icon("heart-pulse")),
      menuItem("7. Scenario comparison",      tabName = "compare", icon = icon("layer-group")),
      menuItem("8. References",               tabName = "refs",    icon = icon("book"))
    ),
    hr(),
    selectInput("scenario", "Scenario:", SCENARIO_LIST,
                selected = "Buprenorphine/naloxone maintenance (16 mg QD)"),
    sliderInput("horizon_d", "Simulation horizon (days):", 1, 60, 14, step = 1),
    sliderInput("age",       "Age (years):", 18, 70, 34),
    sliderInput("wt",        "Weight (kg):",  40,120, 72),
    sliderInput("years_use", "Years of regular opioid use:", 0, 30, 5),
    sliderInput("base_tol",   "Baseline tolerance index (0-1):", 0, 0.95, 0.55, step = 0.05),
    sliderInput("base_cows",  "Baseline inter-dose COWS (0-48):", 0, 20, 2),
    sliderInput("base_crave", "Baseline craving VAS (0-100):", 0,100, 30),
    checkboxInput("counseling", "Counseling / contingency management engaged", FALSE),
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
              p(strong("Model note:"),
                "MOR occupancy is computed from a competitive multi-ligand Emax ",
                "equation across illicit fentanyl, methadone, buprenorphine, ",
                "naloxone, and naltrexone plasma concentrations weighted by ",
                "relative affinity (Ki) and intrinsic efficacy; buprenorphine's ",
                "partial-agonist ceiling is enforced above ~55% occupancy.")
          )
        ),
        fluidRow(
          box(width = 12, title = "About this dashboard", status = "info",
              p("This dashboard runs the mrgsolve QSP model for opioid use ",
                "disorder (OUD). Pick a treatment/use scenario in the left ",
                "panel, adjust the patient profile, then press ",
                strong("Run simulation"), " to update plots."),
              p("Each tab focuses on a layer of the model: PK, MOR occupancy/",
                "tolerance pharmacodynamics, withdrawal (COWS)/craving, ",
                "respiratory-depression/overdose risk, clinical endpoints ",
                "(methadone QTc, treatment retention), scenario comparison, ",
                "and references.")
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
          box(width = 6, title = "Total MOR occupancy (fraction 0-1)", plotOutput("occ_plot", height = 360)),
          box(width = 6, title = "Tolerance index & phasic dopamine/euphoria", plotOutput("tol_da_plot", height = 360))
        )
      ),

      tabItem("wd",
        fluidRow(
          box(width = 6, title = "Locus-coeruleus tone & COWS withdrawal score", plotOutput("wd_plot", height = 360)),
          box(width = 6, title = "Craving VAS (0-100)", plotOutput("crave_plot", height = 360))
        )
      ),

      tabItem("resp",
        fluidRow(
          box(width = 12, title = "Respiratory-drive index (overdose threshold at 0.25, dashed)",
              status = "danger", solidHeader = TRUE, plotOutput("resp_plot", height = 420)),
          box(width = 12, title = "Overdose-flag events over time", DTOutput("overdose_table"))
        )
      ),

      tabItem("clin",
        fluidRow(
          box(width = 6, title = "Methadone-associated QTc (ms)", plotOutput("qtc_plot", height = 360)),
          box(width = 6, title = "Treatment-retention index (0-1)", plotOutput("retain_plot", height = 360))
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
          box(width = 12, title = "Endpoint summary table (start / mid / end of horizon)",
              status = "info", DTOutput("endpoint_table"))
        )
      ),

      tabItem("refs",
        fluidRow(
          box(width = 12, title = "Key references", status = "primary", solidHeader = TRUE,
              p("See ", code("oud_references.md"), " in this directory for the full",
                "curated list (30+ PubMed-linked references)."),
              tags$ul(
                tags$li("Volkow ND, et al. 2016 N Engl J Med — neurobiology of addiction, opponent-process model."),
                tags$li("Koob GF, Volkow ND. 2016 Lancet Psychiatry — anti-reward/dynorphin-KOR system."),
                tags$li("Dahan A, et al. 2005/2010 Anesthesiology/Br J Anaesth — opioid respiratory depression PK/PD."),
                tags$li("Walsh SL, et al. 1994 Clin Pharmacol Ther — buprenorphine dose-response/ceiling."),
                tags$li("Krantz MJ, et al. 2009 Ann Intern Med — methadone QTc prolongation."),
                tags$li("Mattick RP, et al. 2014 Cochrane — buprenorphine/methadone maintenance treatment."),
                tags$li("Sordo L, et al. 2017 BMJ — mortality risk reduction with opioid agonist treatment.")
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
    p <- list(age = input$age, wt = input$wt, years_use = input$years_use,
              base_tol = input$base_tol, base_cows = input$base_cows,
              base_crave = input$base_crave, base_retain = 0.20,
              counseling = input$counseling)
    results(run_sim(input$scenario, input$horizon_d * 24, p))
  }, ignoreNULL = FALSE)

  observeEvent(input$run_all, {
    showNotification("Running 10 scenarios…", type = "message", duration = 1)
    p <- list(age = input$age, wt = input$wt, years_use = input$years_use,
              base_tol = input$base_tol, base_cows = input$base_cows,
              base_crave = input$base_crave, base_retain = 0.20,
              counseling = input$counseling)
    out <- lapply(SCENARIO_LIST, function(sc) run_sim(sc, input$horizon_d * 24, p))
    all_results(bind_rows(out))
  })

  # --- Patient table ---
  output$patient_table <- renderDT({
    tibble(
      Field = c("Age", "Weight", "Years of use", "Baseline tolerance",
                "Baseline COWS", "Baseline craving VAS", "Counseling engaged",
                "Scenario", "Horizon (d)"),
      Value = c(input$age, input$wt, input$years_use, input$base_tol,
                input$base_cows, input$base_crave, input$counseling,
                input$scenario, input$horizon_d)
    ) %>% datatable(rownames = FALSE, options = list(dom = "t"))
  })

  # --- PK ---
  output$pk_plot <- renderPlot({
    df <- results(); req(df)
    df %>% mutate(hour = time) %>%
      select(hour, conc_fentanyl, conc_methadone, conc_buprenorphine,
             conc_naloxone, conc_naltrexone, conc_lofexidine) %>%
      pivot_longer(-hour) %>%
      ggplot(aes(hour, value, colour = name)) +
        geom_line(linewidth = 0.7) +
        scale_y_continuous(trans = "log1p") +
        labs(x = "Time (hours)", y = "Plasma conc (mg/L)", colour = "Drug") +
        theme_minimal(base_size = 13)
  })

  # --- MOR occupancy/tolerance ---
  output$occ_plot <- renderPlot({
    df <- results(); req(df)
    ggplot(df, aes(time, MOR_occupancy_total)) +
      geom_line(colour = "#264653", linewidth = 1) +
      scale_y_continuous(limits = c(0, 1)) +
      labs(x = "Hour", y = "Total MOR occupancy (fraction)") +
      theme_minimal(base_size = 13)
  })
  output$tol_da_plot <- renderPlot({
    df <- results(); req(df)
    df %>% select(time, Tolerance_idx, Dopamine_idx) %>%
      pivot_longer(-time) %>%
      ggplot(aes(time, value, colour = name)) +
        geom_line(linewidth = 0.8) +
        facet_wrap(~name, scales = "free_y") +
        theme_minimal(base_size = 12) + theme(legend.position = "none")
  })

  # --- Withdrawal & craving ---
  output$wd_plot <- renderPlot({
    df <- results(); req(df)
    df %>% select(time, LC_tone_idx, COWS_score) %>%
      pivot_longer(-time) %>%
      ggplot(aes(time, value, colour = name)) +
        geom_line(linewidth = 0.8) +
        facet_wrap(~name, scales = "free_y") +
        theme_minimal(base_size = 12) + theme(legend.position = "none")
  })
  output$crave_plot <- renderPlot({
    df <- results(); req(df)
    ggplot(df, aes(time, Craving_VAS)) +
      geom_line(colour = "#9b1c1c", linewidth = 1) +
      scale_y_continuous(limits = c(0, 100)) +
      labs(x = "Hour", y = "Craving VAS (0-100)") +
      theme_minimal(base_size = 13)
  })

  # --- Respiratory / overdose ---
  output$resp_plot <- renderPlot({
    df <- results(); req(df)
    ggplot(df, aes(time, Respiratory_idx)) +
      geom_line(colour = "#cc3344", linewidth = 1) +
      geom_hline(yintercept = 0.25, lty = 2, colour = "grey30") +
      scale_y_continuous(limits = c(0, 1)) +
      labs(x = "Hour", y = "Relative respiratory-drive index",
           caption = "Dashed line = overdose/apnea threshold (0.25)") +
      theme_minimal(base_size = 13)
  })
  output$overdose_table <- renderDT({
    df <- results(); req(df)
    df %>% filter(Overdose_flag == 1) %>%
      select(time, Respiratory_idx, MOR_occupancy_total) %>%
      datatable(rownames = FALSE, options = list(pageLength = 10))
  })

  # --- Clinical endpoints ---
  output$qtc_plot <- renderPlot({
    df <- results(); req(df)
    ggplot(df, aes(time, QTc_ms)) +
      geom_line(colour = "#6d4c41", linewidth = 1) +
      geom_hline(yintercept = 500, lty = 2, colour = "grey30") +
      labs(x = "Hour", y = "QTc (ms)",
           caption = "Dashed line = commonly used torsades-risk threshold (500 ms)") +
      theme_minimal(base_size = 13)
  })
  output$retain_plot <- renderPlot({
    df <- results(); req(df)
    ggplot(df, aes(time, Retention_idx)) +
      geom_line(colour = "#2a9d8f", linewidth = 1) +
      scale_y_continuous(limits = c(0, 1)) +
      labs(x = "Hour", y = "Treatment-retention index (0-1)") +
      theme_minimal(base_size = 13)
  })

  # --- Compare ---
  output$compare_plot <- renderPlot({
    df <- all_results(); req(df)
    ggplot(df, aes(time, Respiratory_idx, colour = scenario)) +
      geom_line(linewidth = 0.8) +
      geom_hline(yintercept = 0.25, lty = 2, colour = "grey50") +
      labs(x = "Hour", y = "Respiratory-drive index", colour = "Scenario") +
      theme_minimal(base_size = 13) +
      theme(legend.position = "bottom")
  })
  output$endpoint_table <- renderDT({
    df <- all_results(); req(df)
    horizon <- max(df$time)
    mid <- horizon/2
    df %>% filter(time %in% c(0, mid, horizon) | abs(time - mid) < 0.05) %>%
      mutate(Hour = round(time, 1)) %>%
      select(scenario, Hour, MOR_occupancy_total, COWS_score, Craving_VAS,
             Respiratory_idx, QTc_ms, Retention_idx) %>%
      distinct(scenario, Hour, .keep_all = TRUE) %>%
      arrange(scenario, Hour) %>%
      datatable(rownames = FALSE, options = list(pageLength = 30))
  })
}

shinyApp(ui, server)
