## ============================================================================
## Congenital Long QT Syndrome (LQTS) QSP — Shiny Dashboard
## ----------------------------------------------------------------------------
## 8 tabs: Patient/Genotype Profile · Drug PK · Ion Channel/PD · QTc & ECG
##         Surrogate · TdP/Arrhythmia Risk · Clinical Endpoints (Syncope/SCD) ·
##         Scenario Comparison · Biomarkers/Risk Stratification (Schwartz score)
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
MODEL_PATH <- "lqts_mrgsolve_model.R"
get_model <- function() {
  if (!exists(".LQTS_MOD", envir = .GlobalEnv) || is.null(.GlobalEnv$.LQTS_MOD)) {
    assign(".LQTS_MOD", mread_cache("lqts", project = "."), envir = .GlobalEnv)
  }
  .GlobalEnv$.LQTS_MOD
}

## ---------- Scenario builders ----------
build_events <- function(scenario, horizon_h) {
  n_days <- ceiling(horizon_h / 24)
  if (scenario == "Untreated (natural history)") {
    ev(amt = 0, cmt = "CENT_PROP")
  } else if (scenario == "Propranolol 40 mg TID (~2 mg/kg/day)") {
    ev(amt = 40, cmt = "GUT_PROP", ii = 8/24, addl = 3 * n_days)
  } else if (scenario == "Nadolol 60 mg QD (~1 mg/kg/day)") {
    ev(amt = 60, cmt = "GUT_NAD", ii = 1, addl = n_days)
  } else if (scenario == "Mexiletine 200 mg TID (LQT3-targeted)") {
    ev(amt = 200, cmt = "GUT_MEX", ii = 8/24, addl = 3 * n_days)
  } else if (scenario == "Mexiletine + Propranolol (combination)") {
    seq(ev(amt = 200, cmt = "GUT_MEX",  ii = 8/24, addl = 3 * n_days),
        ev(amt = 40,  cmt = "GUT_PROP", ii = 8/24, addl = 3 * n_days))
  } else if (scenario == "Beta-blocker + K+/spironolactone (LQT2)") {
    seq(ev(amt = 40,  cmt = "GUT_PROP",  ii = 8/24, addl = 3 * n_days),
        ev(amt = 600, cmt = "GUT_KCL",   ii = 8/24, addl = 3 * n_days),
        ev(amt = 25,  cmt = "GUT_SPIRO", ii = 1,    addl = n_days))
  } else if (scenario == "Added QT-prolonging drug X (acquired-on-congenital)") {
    ev(amt = 500, cmt = "GUT_QTX", ii = 1, addl = n_days)
  } else if (scenario == "LCSD + beta-blocker + ICD (high-risk)") {
    ev(amt = 40, cmt = "GUT_PROP", ii = 8/24, addl = 3 * n_days)
  } else {
    ev(amt = 0, cmt = "CENT_PROP")
  }
}

SCENARIO_LIST <- c(
  "Untreated (natural history)",
  "Propranolol 40 mg TID (~2 mg/kg/day)",
  "Nadolol 60 mg QD (~1 mg/kg/day)",
  "Mexiletine 200 mg TID (LQT3-targeted)",
  "Mexiletine + Propranolol (combination)",
  "Beta-blocker + K+/spironolactone (LQT2)",
  "Added QT-prolonging drug X (acquired-on-congenital)",
  "LCSD + beta-blocker + ICD (high-risk)"
)

run_sim <- function(scenario, horizon_h, params) {
  mod <- get_model()
  ev_set <- build_events(scenario, horizon_h)
  lcsd_flag  <- as.integer(scenario == "LCSD + beta-blocker + ICD (high-risk)")
  icd_flag   <- as.integer(scenario == "LCSD + beta-blocker + ICD (high-risk)")
  par <- list(
    AGE           = params$age,
    WT            = params$wt,
    GENOTYPE      = as.numeric(params$genotype),
    FEMALE        = as.integer(params$female),
    BASE_QTC      = params$base_qtc,
    PRIOR_SYNCOPE = as.integer(params$prior_syncope),
    TRIGGER       = as.numeric(params$trigger),
    SERUM_K_BASE  = params$serum_k,
    LCSD_DONE     = lcsd_flag,
    ICD_PRESENT   = icd_flag
  )
  mod %>% param(par) %>% mrgsim(events = ev_set, end = horizon_h, delta = 6) %>%
    as_tibble() %>% mutate(scenario = scenario)
}

## Schwartz risk score approximation (simplified, for illustrative biomarker tab)
schwartz_score <- function(qtc, twave_notch, torsades_hist, syncope, congenital_deaf,
                            resting_hr_low, female, family_lqts, family_scd) {
  score <- 0
  score <- score + if (qtc >= 480) 3 else if (qtc >= 460) 2 else if (qtc >= 450 && !female) 1 else 0
  score <- score + if (torsades_hist) 2 else 0
  score <- score + if (twave_notch) 1 else 0
  score <- score + if (resting_hr_low) 0.5 else 0
  score <- score + if (syncope) 2 else 0
  score <- score + if (congenital_deaf) 0.5 else 0
  score <- score + if (family_lqts) 1 else 0
  score <- score + if (family_scd) 0.5 else 0
  score
}

## ---------- UI ----------
ui <- dashboardPage(
  dashboardHeader(title = "LQTS QSP"),
  dashboardSidebar(
    sidebarMenu(
      menuItem("1. Patient/Genotype Profile", tabName = "patient", icon = icon("user")),
      menuItem("2. Drug PK",                  tabName = "pk",      icon = icon("syringe")),
      menuItem("3. Ion Channel / PD",         tabName = "pd",      icon = icon("bolt")),
      menuItem("4. QTc & ECG Surrogate",      tabName = "qtc",     icon = icon("heartbeat")),
      menuItem("5. TdP / Arrhythmia Risk",    tabName = "tdp",     icon = icon("exclamation-triangle")),
      menuItem("6. Clinical Endpoints",       tabName = "clin",    icon = icon("hospital")),
      menuItem("7. Scenario Comparison",      tabName = "compare", icon = icon("layer-group")),
      menuItem("8. Biomarkers / Risk Strat.", tabName = "biomk",   icon = icon("clipboard-list"))
    ),
    hr(),
    selectInput("scenario", "Scenario:", SCENARIO_LIST, selected = "Propranolol 40 mg TID (~2 mg/kg/day)"),
    sliderInput("horizon_d", "Simulation horizon (days):", 7, 365, 90, step = 1),
    selectInput("genotype", "Genotype:",
                choices = c("LQT1 (KCNQ1/IKs)" = 1, "LQT2 (KCNH2-hERG/IKr)" = 2, "LQT3 (SCN5A/late-INa)" = 3),
                selected = 1),
    selectInput("trigger", "Dominant trigger:",
                choices = c("Rest/baseline" = 0, "Exercise (classic LQT1)" = 1,
                            "Auditory/emotion (classic LQT2)" = 2, "Sleep/bradycardia (classic LQT3)" = 3),
                selected = 1),
    checkboxInput("female", "Female", value = FALSE),
    checkboxInput("prior_syncope", "Prior syncope / aborted cardiac arrest", value = FALSE),
    sliderInput("age", "Age (years):", 1, 70, 25),
    sliderInput("wt",  "Weight (kg):", 10, 120, 65),
    sliderInput("base_qtc", "Baseline (untreated) QTc (ms):", 420, 520, 470),
    sliderInput("serum_k", "Baseline serum K+ (mEq/L):", 3.0, 5.5, 4.0, step = 0.1),
    actionButton("run", "Run simulation", icon = icon("play"),
                 style = "color:#fff;background:#1f6feb")
  ),
  dashboardBody(
    tabItems(
      tabItem("patient",
        fluidRow(
          box(width = 12, title = "Patient / genotype profile summary", status = "primary",
              solidHeader = TRUE,
              DTOutput("patient_table"),
              br(),
              p(strong("Genotype legend:"),
                "LQT1 (KCNQ1/IKs loss-of-function) — exercise/swimming triggers, ",
                "best beta-blocker response. LQT2 (KCNH2-hERG/IKr loss-of-function) — ",
                "auditory/emotional triggers, female excess risk post-puberty. ",
                "LQT3 (SCN5A gain-of-function, persistent late-INa) — rest/sleep/",
                "bradycardia triggers, genotype-targeted mexiletine therapy.")
          )
        ),
        fluidRow(
          box(width = 12, title = "About this dashboard", status = "info",
              p("This dashboard runs the mrgsolve QSP model for congenital Long QT ",
                "Syndrome (LQTS), linking genotype-specific ion-channel conductance ",
                "deficits (IKs/IKr/late-INa) to a lumped repolarization-reserve/QTc ",
                "surrogate, EAD-substrate generation, autonomic triggering, and ",
                "cumulative Torsades de Pointes (TdP) risk."),
              p("Pick a genotype, dominant trigger, and treatment scenario in the ",
                "left panel, adjust the patient profile, then press ",
                strong("Run simulation"), " to update plots across all tabs.")
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
          box(width = 6, title = "Channel conductance indices (GKs / GKr / GNa-late)",
              plotOutput("channel_plot", height = 380)),
          box(width = 6, title = "Sympathetic drive index & serum K+",
              plotOutput("symp_k_plot", height = 380))
        ),
        fluidRow(
          box(width = 12, title = "qNet-like torsadogenic-risk metric (CiPA-style)",
              plotOutput("qnet_plot", height = 320))
        )
      ),

      tabItem("qtc",
        fluidRow(
          box(width = 12, title = "QTc surrogate over time (ECG proxy)",
              status = "primary", solidHeader = TRUE, plotOutput("qtc_plot", height = 480))
        )
      ),

      tabItem("tdp",
        fluidRow(
          box(width = 6, title = "EAD / arrhythmogenic substrate index",
              plotOutput("ead_plot", height = 380)),
          box(width = 6, title = "Cumulative TdP-risk hazard",
              plotOutput("hazard_plot", height = 380))
        ),
        fluidRow(
          box(width = 12, title = "Cumulative TdP event probability",
              plotOutput("tdp_prob_plot", height = 380))
        )
      ),

      tabItem("clin",
        fluidRow(
          box(width = 12, title = "Clinical endpoint summary (syncope / TdP / SCD proxy)",
              status = "warning", solidHeader = TRUE, DTOutput("clin_table"))
        ),
        fluidRow(
          box(width = 12, title = "Interpretation notes", status = "info",
              p("TdP_event_probability = 1 - exp(-cumulative hazard), a Poisson-",
                "process approximation of the probability of at least one TdP ",
                "episode over the simulated horizon. A subset of sustained TdP/VF ",
                "events are assumed fatal (sudden cardiac death) unless an ICD is ",
                "present and successfully rescues the patient (ICD_RESCUE_FRAC)."))
        )
      ),

      tabItem("compare",
        fluidRow(
          box(width = 12, title = "Scenario comparison panel",
              status = "warning", solidHeader = TRUE,
              p("Runs all eight built-in treatment scenarios with the current ",
                "patient profile; press the button below."),
              actionButton("run_all", "Run all scenarios", icon = icon("rocket"),
                           style = "color:#fff;background:#0f5132"),
              br(), br(),
              plotOutput("compare_qtc_plot", height = 380),
              br(),
              plotOutput("compare_risk_plot", height = 380)
          )
        ),
        fluidRow(
          box(width = 12, title = "Endpoint summary table (Day-0, Day-30, Day-end)",
              status = "info", DTOutput("endpoint_table"))
        )
      ),

      tabItem("biomk",
        fluidRow(
          box(width = 6, title = "Schwartz risk-score calculator (illustrative)",
              status = "primary", solidHeader = TRUE,
              checkboxInput("sw_notch", "Notched T-wave (low amplitude)", FALSE),
              checkboxInput("sw_tdp_hist", "Documented Torsades de Pointes", FALSE),
              checkboxInput("sw_low_hr", "Resting bradycardia for age", FALSE),
              checkboxInput("sw_deaf", "Congenital deafness (Jervell-Lange-Nielsen)", FALSE),
              checkboxInput("sw_fam_lqts", "Family member with definite LQTS", FALSE),
              checkboxInput("sw_fam_scd", "Unexplained SCD <30y in immediate family", FALSE),
              verbatimTextOutput("schwartz_result")
          ),
          box(width = 6, title = "Risk stratification biomarkers",
              status = "primary", solidHeader = TRUE,
              DTOutput("risk_biomarker_table"))
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
      age = input$age, wt = input$wt, genotype = input$genotype,
      female = input$female, base_qtc = input$base_qtc,
      prior_syncope = input$prior_syncope, trigger = input$trigger,
      serum_k = input$serum_k
    )
    results(run_sim(input$scenario, input$horizon_d * 24, p))
  }, ignoreNULL = FALSE)

  observeEvent(input$run_all, {
    showNotification("Running 8 scenarios…", type = "message", duration = 1)
    p <- list(
      age = input$age, wt = input$wt, genotype = input$genotype,
      female = input$female, base_qtc = input$base_qtc,
      prior_syncope = input$prior_syncope, trigger = input$trigger,
      serum_k = input$serum_k
    )
    out <- lapply(SCENARIO_LIST, function(sc) run_sim(sc, input$horizon_d * 24, p))
    all_results(bind_rows(out))
  })

  # --- Patient table ---
  output$patient_table <- renderDT({
    genotype_lbl <- c("1"="LQT1 (KCNQ1/IKs)", "2"="LQT2 (KCNH2-hERG/IKr)",
                       "3"="LQT3 (SCN5A/late-INa)")[[as.character(input$genotype)]]
    trigger_lbl <- c("0"="Rest/baseline", "1"="Exercise", "2"="Auditory/emotion",
                      "3"="Sleep/bradycardia")[[as.character(input$trigger)]]
    tibble(
      Field = c("Age", "Weight", "Genotype", "Dominant trigger", "Female",
                "Prior syncope/arrest", "Baseline QTc", "Baseline K+",
                "Scenario", "Horizon (d)"),
      Value = c(input$age, input$wt, genotype_lbl, trigger_lbl, input$female,
                input$prior_syncope, paste0(input$base_qtc, " ms"),
                paste0(input$serum_k, " mEq/L"), input$scenario, input$horizon_d)
    ) %>% datatable(rownames = FALSE, options = list(dom = "t"))
  })

  # --- PK ---
  output$pk_plot <- renderPlot({
    df <- results(); req(df)
    df %>% mutate(day = time/24) %>%
      select(day, conc_propranolol, conc_nadolol, conc_mexiletine, conc_qtdrugX) %>%
      pivot_longer(-day) %>%
      ggplot(aes(day, value, colour = name)) +
        geom_line(linewidth = 0.7) +
        labs(x = "Time (days)", y = "Plasma conc (mg/L)", colour = "Drug") +
        theme_minimal(base_size = 13)
  })

  # --- Ion Channel / PD ---
  output$channel_plot <- renderPlot({
    df <- results(); req(df)
    df %>% mutate(day = time/24) %>%
      select(day, GKs_idx, GKr_idx, GNaLate_idx) %>%
      pivot_longer(-day) %>%
      ggplot(aes(day, value, colour = name)) +
        geom_line(linewidth = 0.8) +
        geom_hline(yintercept = 1.0, lty = 2, colour = "grey50") +
        labs(x = "Day", y = "Relative conductance index", colour = "") +
        theme_minimal(base_size = 12)
  })
  output$symp_k_plot <- renderPlot({
    df <- results(); req(df)
    df %>% mutate(day = time/24) %>%
      select(day, SympDrive_idx, SerumK_mEqL) %>%
      pivot_longer(-day) %>%
      ggplot(aes(day, value, colour = name)) +
        geom_line(linewidth = 0.8) +
        facet_wrap(~name, scales = "free_y") +
        theme_minimal(base_size = 12) + theme(legend.position = "none")
  })
  output$qnet_plot <- renderPlot({
    df <- results(); req(df)
    ggplot(df, aes(time/24, qNet_like_metric)) +
      geom_line(colour = "#6a4c93", linewidth = 1) +
      labs(x = "Day", y = "qNet-like metric (higher = lower TdP risk)") +
      theme_minimal(base_size = 13)
  })

  # --- QTc ---
  output$qtc_plot <- renderPlot({
    df <- results(); req(df)
    ggplot(df, aes(time/24, QTc_ms)) +
      geom_line(colour = "#1f6feb", linewidth = 1) +
      geom_hline(yintercept = 500, lty = 2, colour = "#b02a37") +
      geom_hline(yintercept = 450, lty = 3, colour = "grey50") +
      labs(x = "Day", y = "QTc (ms)",
           caption = "Dashed red = high-risk threshold (500 ms); dotted = upper-normal (~450 ms)") +
      theme_minimal(base_size = 13)
  })

  # --- TdP / Arrhythmia risk ---
  output$ead_plot <- renderPlot({
    df <- results(); req(df)
    ggplot(df, aes(time/24, EAD_substrate_idx)) +
      geom_line(colour = "#9b1c1c", linewidth = 1) +
      scale_y_continuous(limits = c(0, 1)) +
      labs(x = "Day", y = "EAD-substrate index (0-1)") +
      theme_minimal(base_size = 13)
  })
  output$hazard_plot <- renderPlot({
    df <- results(); req(df)
    ggplot(df, aes(time/24, TdP_cumulative_hazard)) +
      geom_line(colour = "#e07a00", linewidth = 1) +
      labs(x = "Day", y = "Cumulative TdP hazard (integral)") +
      theme_minimal(base_size = 13)
  })
  output$tdp_prob_plot <- renderPlot({
    df <- results(); req(df)
    ggplot(df, aes(time/24, TdP_event_probability)) +
      geom_line(colour = "#cc3344", linewidth = 1) +
      scale_y_continuous(limits = c(0, 1), labels = scales::percent) +
      labs(x = "Day", y = "Cumulative P(>=1 TdP event)") +
      theme_minimal(base_size = 13)
  })

  # --- Clinical endpoints ---
  output$clin_table <- renderDT({
    df <- results(); req(df)
    horizon <- max(df$time)
    tibble(
      Endpoint = c("Final QTc (ms)", "Final EAD-substrate index",
                   "Cumulative TdP hazard", "Cumulative P(>=1 TdP event)",
                   "qNet-like metric (final)"),
      Value = c(round(tail(df$QTc_ms,1),1), round(tail(df$EAD_substrate_idx,1),3),
                round(tail(df$TdP_cumulative_hazard,1),4),
                paste0(round(100*tail(df$TdP_event_probability,1),2), "%"),
                round(tail(df$qNet_like_metric,1),3))
    ) %>% datatable(rownames = FALSE, options = list(dom = "t"))
  })

  # --- Compare ---
  output$compare_qtc_plot <- renderPlot({
    df <- all_results(); req(df)
    df %>% mutate(day = time/24) %>%
      ggplot(aes(day, QTc_ms, colour = scenario)) +
        geom_line(linewidth = 0.8) +
        geom_hline(yintercept = 500, lty = 2, colour = "grey40") +
        labs(x = "Day", y = "QTc (ms)", colour = "Scenario") +
        theme_minimal(base_size = 12) +
        theme(legend.position = "bottom")
  })
  output$compare_risk_plot <- renderPlot({
    df <- all_results(); req(df)
    df %>% mutate(day = time/24) %>%
      ggplot(aes(day, TdP_event_probability, colour = scenario)) +
        geom_line(linewidth = 0.8) +
        labs(x = "Day", y = "Cumulative P(>=1 TdP event)", colour = "Scenario") +
        theme_minimal(base_size = 12) +
        theme(legend.position = "bottom")
  })
  output$endpoint_table <- renderDT({
    df <- all_results(); req(df)
    horizon <- max(df$time)
    day30 <- min(30*24, horizon)
    df %>% filter(time %in% c(0, day30, horizon)) %>%
      mutate(Day = round(time/24)) %>%
      select(scenario, Day, QTc_ms, EAD_substrate_idx, TdP_event_probability) %>%
      distinct(scenario, Day, .keep_all = TRUE) %>%
      arrange(scenario, Day) %>%
      datatable(rownames = FALSE, options = list(pageLength = 30))
  })

  # --- Biomarkers / Risk stratification ---
  output$schwartz_result <- renderPrint({
    score <- schwartz_score(
      qtc = input$base_qtc, twave_notch = input$sw_notch,
      torsades_hist = input$sw_tdp_hist, syncope = input$prior_syncope,
      congenital_deaf = input$sw_deaf, resting_hr_low = input$sw_low_hr,
      female = input$female, family_lqts = input$sw_fam_lqts,
      family_scd = input$sw_fam_scd
    )
    risk_cat <- if (score >= 3.5) "High probability of LQTS" else
                if (score >= 1.5) "Intermediate probability" else "Low probability"
    cat("Schwartz score (illustrative approximation):", score, "\n")
    cat("Risk category:", risk_cat, "\n")
  })
  output$risk_biomarker_table <- renderDT({
    tibble(
      Biomarker = c("QTc (Bazett)", "QTc (Fridericia)", "T-wave morphology",
                    "Exercise-stress QTc response", "Genotype", "Sex",
                    "Prior syncope/arrest", "Schwartz score"),
      `Typical high-risk pattern` = c(">500 ms", ">480 ms",
                                       "Broad-based (LQT1) / notched low-amplitude (LQT2) / late peaked (LQT3)",
                                       "Failure to shorten / paradoxical prolongation (LQT1)",
                                       "LQT2 or LQT3", "Female (esp. LQT2 post-puberty)",
                                       "Present", ">=3.5")
    ) %>% datatable(rownames = FALSE, options = list(dom = "t"))
  })
}

shinyApp(ui, server)
