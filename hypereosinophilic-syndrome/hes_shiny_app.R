## ============================================================================
## Hypereosinophilic Syndrome (HES) QSP — Shiny Dashboard
## ----------------------------------------------------------------------------
## 8 tabs: Patient · Drug PK · Kinase/IL-5/Marrow PD · Eosinophil & Tissue Burden ·
##         Cardiac & Clinical Endpoints · Biomarkers · Scenario comparison · References
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

`%||%` <- function(a, b) if (is.null(a)) b else a

## ---------- Lazy model load ----------
MODEL_PATH <- "hes_mrgsolve_model.R"
get_model <- function() {
  if (!exists(".HES_MOD", envir = .GlobalEnv) || is.null(.GlobalEnv$.HES_MOD)) {
    assign(".HES_MOD", mread_cache("hes", project = "."), envir = .GlobalEnv)
  }
  .GlobalEnv$.HES_MOD
}

## ---------- Scenario builders ----------
SCENARIO_LIST <- c(
  "1. Untreated FIP1L1-PDGFRA+ M-HES (natural history)",
  "2. FIP1L1-PDGFRA+ M-HES + Imatinib 100 mg QD",
  "3. Imatinib-resistant (T674I) -> switch to Hydroxyurea",
  "4. Idiopathic HES + Prednisone (steroid-responsive)",
  "5. Idiopathic HES steroid-refractory + Mepolizumab 300 mg Q4W",
  "6. Lymphocytic-variant HES + Benralizumab 30 mg Q4W",
  "7. PDGFRB-rearranged HES + low-dose Imatinib 100 mg QD",
  "8. Cardiac Loeffler endocarditis + high-dose steroid pulse",
  "9. JAK2-rearranged HES + Ruxolitinib 20 mg BID",
  "10. Steroid taper + Mepolizumab steroid-sparing combination"
)

scenario_params <- function(scenario) {
  switch(scenario,
    "1. Untreated FIP1L1-PDGFRA+ M-HES (natural history)" =
      list(GENOTYPE = 1, T674I_MUT = 0, CARDIAC_INVOLVED = 0),
    "2. FIP1L1-PDGFRA+ M-HES + Imatinib 100 mg QD" =
      list(GENOTYPE = 1, T674I_MUT = 0, CARDIAC_INVOLVED = 0),
    "3. Imatinib-resistant (T674I) -> switch to Hydroxyurea" =
      list(GENOTYPE = 1, T674I_MUT = 1, CARDIAC_INVOLVED = 0),
    "4. Idiopathic HES + Prednisone (steroid-responsive)" =
      list(GENOTYPE = 0, T674I_MUT = 0, CARDIAC_INVOLVED = 0),
    "5. Idiopathic HES steroid-refractory + Mepolizumab 300 mg Q4W" =
      list(GENOTYPE = 0, T674I_MUT = 0, CARDIAC_INVOLVED = 0),
    "6. Lymphocytic-variant HES + Benralizumab 30 mg Q4W" =
      list(GENOTYPE = 0, T674I_MUT = 0, CARDIAC_INVOLVED = 0),
    "7. PDGFRB-rearranged HES + low-dose Imatinib 100 mg QD" =
      list(GENOTYPE = 2, T674I_MUT = 0, CARDIAC_INVOLVED = 0),
    "8. Cardiac Loeffler endocarditis + high-dose steroid pulse" =
      list(GENOTYPE = 0, T674I_MUT = 0, CARDIAC_INVOLVED = 1, BASE_CARDIAC = 0.35),
    "9. JAK2-rearranged HES + Ruxolitinib 20 mg BID" =
      list(GENOTYPE = 3, T674I_MUT = 0, CARDIAC_INVOLVED = 0),
    "10. Steroid taper + Mepolizumab steroid-sparing combination" =
      list(GENOTYPE = 0, T674I_MUT = 0, CARDIAC_INVOLVED = 0)
  )
}

build_events <- function(scenario, horizon_h) {
  nd  <- ceiling(horizon_h / 24)
  n4w <- ceiling(horizon_h / (24 * 28))
  if (scenario == "1. Untreated FIP1L1-PDGFRA+ M-HES (natural history)") {
    ev(amt = 0, cmt = "GUT_IMA", time = 0)
  } else if (scenario == "2. FIP1L1-PDGFRA+ M-HES + Imatinib 100 mg QD") {
    ev(amt = 100, cmt = "GUT_IMA", ii = 24, addl = nd)
  } else if (scenario == "3. Imatinib-resistant (T674I) -> switch to Hydroxyurea") {
    seq(ev(amt = 400, cmt = "GUT_IMA", ii = 24, addl = 29),
        ev(amt = 1000, cmt = "GUT_HU", ii = 24, time = 30*24, addl = nd))
  } else if (scenario == "4. Idiopathic HES + Prednisone (steroid-responsive)") {
    ev(amt = 60, cmt = "GUT_PRED", ii = 24, addl = nd)
  } else if (scenario == "5. Idiopathic HES steroid-refractory + Mepolizumab 300 mg Q4W") {
    ev(amt = 300, cmt = "DEPOT_MEPO", ii = 24*28, addl = n4w)
  } else if (scenario == "6. Lymphocytic-variant HES + Benralizumab 30 mg Q4W") {
    ev(amt = 30, cmt = "DEPOT_BENRA", ii = 24*28, addl = n4w)
  } else if (scenario == "7. PDGFRB-rearranged HES + low-dose Imatinib 100 mg QD") {
    ev(amt = 100, cmt = "GUT_IMA", ii = 24, addl = nd)
  } else if (scenario == "8. Cardiac Loeffler endocarditis + high-dose steroid pulse") {
    ev(amt = 100, cmt = "GUT_PRED", ii = 24, addl = nd)
  } else if (scenario == "9. JAK2-rearranged HES + Ruxolitinib 20 mg BID") {
    ev(amt = 20, cmt = "GUT_RUX", ii = 12, addl = 2*nd)
  } else if (scenario == "10. Steroid taper + Mepolizumab steroid-sparing combination") {
    seq(ev(amt = 40, cmt = "GUT_PRED", ii = 24, addl = 13),
        ev(amt = 300, cmt = "DEPOT_MEPO", ii = 24*28, time = 14*24, addl = n4w))
  } else {
    ev(amt = 0, cmt = "GUT_IMA")
  }
}

run_sim <- function(scenario, horizon_h, params) {
  mod <- get_model()
  ev_set <- build_events(scenario, horizon_h)
  sc_par <- scenario_params(scenario)
  par <- modifyList(list(
    AGE = params$age, WT = params$wt, BASE_AEC = params$base_aec
  ), sc_par)
  mod %>% param(par) %>% mrgsim(events = ev_set, end = horizon_h, delta = 1) %>%
    as_tibble() %>% mutate(scenario = scenario)
}

## ---------- UI ----------
ui <- dashboardPage(
  dashboardHeader(title = "HES QSP"),
  dashboardSidebar(
    sidebarMenu(
      menuItem("1. Patient profile",           tabName = "patient", icon = icon("user")),
      menuItem("2. Drug PK",                   tabName = "pk",      icon = icon("syringe")),
      menuItem("3. Kinase/IL-5/Marrow PD",      tabName = "pd",      icon = icon("dna")),
      menuItem("4. Eosinophil & tissue burden", tabName = "eos",     icon = icon("circle-notch")),
      menuItem("5. Cardiac & clinical endpoints", tabName = "clin",  icon = icon("heart-pulse")),
      menuItem("6. Biomarkers",                tabName = "bio",     icon = icon("flask")),
      menuItem("7. Scenario comparison",       tabName = "compare", icon = icon("layer-group")),
      menuItem("8. References",                tabName = "refs",    icon = icon("book"))
    ),
    hr(),
    selectInput("scenario", "Scenario:", SCENARIO_LIST,
                selected = "2. FIP1L1-PDGFRA+ M-HES + Imatinib 100 mg QD"),
    sliderInput("horizon_d", "Simulation horizon (days):", 7, 365, 180, step = 1),
    sliderInput("age",       "Age (years):", 18, 85, 42),
    sliderInput("wt",        "Weight (kg):",  40, 130, 75),
    sliderInput("base_aec",  "Baseline untreated AEC (cells/uL):", 1500, 30000, 8000, step = 500),
    actionButton("run", "Run simulation", icon = icon("play"),
                 style = "color:#fff;background:#1f6feb")
  ),
  dashboardBody(
    tabItems(
      tabItem("patient",
        fluidRow(
          box(width = 12, title = "Patient/genotype profile summary", status = "primary",
              solidHeader = TRUE,
              DTOutput("patient_table"),
              br(),
              p(strong("Model note:"),
                "Genotype (FIP1L1-PDGFRA+, PDGFRB-rearranged, JAK2-rearranged, or ",
                "idiopathic/lymphocytic-variant) sets the baseline constitutive ",
                "kinase-signaling or IL-5 drive; treatment sensitivity (imatinib, ",
                "ruxolitinib, corticosteroids, anti-IL-5-axis biologics) is genotype-",
                "matched as in clinical practice.")
          )
        ),
        fluidRow(
          box(width = 12, title = "About this dashboard", status = "info",
              p("This dashboard runs the mrgsolve QSP model for hypereosinophilic ",
                "syndrome (HES). Pick a genotype/treatment scenario in the left ",
                "panel, adjust the patient profile, then press ",
                strong("Run simulation"), " to update plots."),
              p("Each tab focuses on a layer of the model: PK, kinase/IL-5/marrow ",
                "pharmacodynamics, circulating/tissue eosinophil burden, cardiac ",
                "(Loeffler endocarditis) and composite clinical endpoints, ",
                "biomarkers (ECP, troponin), scenario comparison, and references.")
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
          box(width = 6, title = "Kinase-signaling index (FIP1L1-PDGFRA/PDGFRB/JAK2)", plotOutput("kin_plot", height = 360)),
          box(width = 6, title = "Endogenous IL-5 index", plotOutput("il5_plot", height = 360))
        ),
        fluidRow(
          box(width = 12, title = "Marrow eosinophil-production drive", plotOutput("marrow_plot", height = 360))
        )
      ),

      tabItem("eos",
        fluidRow(
          box(width = 6, title = "Absolute eosinophil count (AEC, log scale)", plotOutput("aec_plot", height = 380)),
          box(width = 6, title = "Tissue eosinophil burden (relative)", plotOutput("tissue_plot", height = 380))
        )
      ),

      tabItem("clin",
        fluidRow(
          box(width = 6, title = "Cardiac (Loeffler) damage index (0-1)", plotOutput("cardiac_plot", height = 360)),
          box(width = 6, title = "Composite HES symptom/flare score (0-10)", plotOutput("symptom_plot", height = 360))
        ),
        fluidRow(
          box(width = 12, title = "Flare events (symptom score >= 4)", DTOutput("flare_table"))
        )
      ),

      tabItem("bio",
        fluidRow(
          box(width = 6, title = "Serum ECP (ug/L)", plotOutput("ecp_plot", height = 360)),
          box(width = 6, title = "Cardiac troponin (ng/L)", plotOutput("trop_plot", height = 360))
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
              p("See ", code("hes_references.md"), " in this directory for the full",
                "curated list (30+ PubMed-linked references)."),
              tags$ul(
                tags$li("Klion AD. 2022 Blood — HES diagnosis/classification review."),
                tags$li("Cools J, et al. 2003 N Engl J Med — FIP1L1-PDGFRA discovery, imatinib sensitivity."),
                tags$li("Cools J, et al. 2004 Cancer Cell — T674I imatinib-resistance mutation."),
                tags$li("Roufosse F, et al. 2020 J Allergy Clin Immunol — mepolizumab HES pivotal trial."),
                tags$li("Kuang FL, et al. 2022 N Engl J Med — mepolizumab HES relapse-prevention trial."),
                tags$li("Kuang FL, et al. 2019 J Allergy Clin Immunol — benralizumab open-label HES trial."),
                tags$li("Parrillo JE, et al. 1979 Ann Intern Med — Loeffler endocarditis natural history.")
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
    p <- list(age = input$age, wt = input$wt, base_aec = input$base_aec)
    results(run_sim(input$scenario, input$horizon_d * 24, p))
  }, ignoreNULL = FALSE)

  observeEvent(input$run_all, {
    showNotification("Running 10 scenarios…", type = "message", duration = 1)
    p <- list(age = input$age, wt = input$wt, base_aec = input$base_aec)
    out <- lapply(SCENARIO_LIST, function(sc) run_sim(sc, input$horizon_d * 24, p))
    all_results(bind_rows(out))
  })

  # --- Patient table ---
  output$patient_table <- renderDT({
    sc_par <- scenario_params(input$scenario)
    tibble(
      Field = c("Age", "Weight", "Baseline AEC", "Scenario", "Genotype code",
                "T674I mutation", "Cardiac involved at baseline", "Horizon (d)"),
      Value = c(input$age, input$wt, input$base_aec, input$scenario,
                sc_par$GENOTYPE %||% NA, sc_par$T674I_MUT %||% 0,
                sc_par$CARDIAC_INVOLVED %||% 0, input$horizon_d)
    ) %>% datatable(rownames = FALSE, options = list(dom = "t"))
  })

  # --- PK ---
  output$pk_plot <- renderPlot({
    df <- results(); req(df)
    df %>% select(time, Cp_imatinib, Cp_prednisolone, Cp_mepolizumab,
                  Cp_benralizumab, Cp_hydroxyurea, Cp_ruxolitinib) %>%
      pivot_longer(-time) %>%
      ggplot(aes(time/24, value, colour = name)) +
        geom_line(linewidth = 0.7) +
        scale_y_continuous(trans = "log1p") +
        labs(x = "Time (days)", y = "Plasma conc (mg/L)", colour = "Drug") +
        theme_minimal(base_size = 13)
  })

  # --- Kinase/IL-5/Marrow PD ---
  output$kin_plot <- renderPlot({
    df <- results(); req(df)
    ggplot(df, aes(time/24, KIN)) +
      geom_line(colour = "#264653", linewidth = 1) +
      labs(x = "Day", y = "Kinase-signaling index (AU)") +
      theme_minimal(base_size = 13)
  })
  output$il5_plot <- renderPlot({
    df <- results(); req(df)
    ggplot(df, aes(time/24, IL5)) +
      geom_line(colour = "#e76f51", linewidth = 1) +
      labs(x = "Day", y = "Endogenous IL-5 index (AU)") +
      theme_minimal(base_size = 13)
  })
  output$marrow_plot <- renderPlot({
    df <- results(); req(df)
    ggplot(df, aes(time/24, MARROW)) +
      geom_line(colour = "#8d6a9f", linewidth = 1) +
      labs(x = "Day", y = "Marrow eosinophil-production drive (relative)") +
      theme_minimal(base_size = 13)
  })

  # --- Eosinophil & tissue burden ---
  output$aec_plot <- renderPlot({
    df <- results(); req(df)
    ggplot(df, aes(time/24, AEC)) +
      geom_line(colour = "#2a9d8f", linewidth = 1) +
      geom_hline(yintercept = 1500, lty = 2, colour = "grey40") +
      scale_y_log10() +
      labs(x = "Day", y = "AEC (cells/uL, log scale)",
           caption = "Dashed line = HES diagnostic threshold (1500/uL)") +
      theme_minimal(base_size = 13)
  })
  output$tissue_plot <- renderPlot({
    df <- results(); req(df)
    ggplot(df, aes(time/24, TISSUE)) +
      geom_line(colour = "#9b1c1c", linewidth = 1) +
      labs(x = "Day", y = "Tissue eosinophil burden (relative to baseline)") +
      theme_minimal(base_size = 13)
  })

  # --- Cardiac & clinical endpoints ---
  output$cardiac_plot <- renderPlot({
    df <- results(); req(df)
    ggplot(df, aes(time/24, CARDIAC)) +
      geom_line(colour = "#cc3344", linewidth = 1) +
      scale_y_continuous(limits = c(0, 1)) +
      labs(x = "Day", y = "Cardiac (Loeffler) damage index") +
      theme_minimal(base_size = 13)
  })
  output$symptom_plot <- renderPlot({
    df <- results(); req(df)
    ggplot(df, aes(time/24, SYMPTOM)) +
      geom_line(colour = "#6d4c41", linewidth = 1) +
      geom_hline(yintercept = 4, lty = 2, colour = "grey40") +
      scale_y_continuous(limits = c(0, 10)) +
      labs(x = "Day", y = "Composite HES symptom/flare score (0-10)",
           caption = "Dashed line = flare threshold (>=4)") +
      theme_minimal(base_size = 13)
  })
  output$flare_table <- renderDT({
    df <- results(); req(df)
    df %>% filter(Flare_flag == 1) %>%
      select(time, SYMPTOM, AEC, CARDIAC) %>%
      mutate(Day = round(time/24, 1)) %>%
      select(Day, SYMPTOM, AEC, CARDIAC) %>%
      datatable(rownames = FALSE, options = list(pageLength = 10))
  })

  # --- Biomarkers ---
  output$ecp_plot <- renderPlot({
    df <- results(); req(df)
    ggplot(df, aes(time/24, ECP)) +
      geom_line(colour = "#e9c46a", linewidth = 1) +
      labs(x = "Day", y = "Serum ECP (ug/L)") +
      theme_minimal(base_size = 13)
  })
  output$trop_plot <- renderPlot({
    df <- results(); req(df)
    ggplot(df, aes(time/24, TROP)) +
      geom_line(colour = "#457b9d", linewidth = 1) +
      labs(x = "Day", y = "Cardiac troponin (ng/L)") +
      theme_minimal(base_size = 13)
  })

  # --- Compare ---
  output$compare_plot <- renderPlot({
    df <- all_results(); req(df)
    ggplot(df, aes(time/24, AEC, colour = scenario)) +
      geom_line(linewidth = 0.8) +
      geom_hline(yintercept = 1500, lty = 2, colour = "grey50") +
      scale_y_log10() +
      labs(x = "Day", y = "AEC (cells/uL, log scale)", colour = "Scenario") +
      theme_minimal(base_size = 13) +
      theme(legend.position = "bottom")
  })
  output$endpoint_table <- renderDT({
    df <- all_results(); req(df)
    horizon <- max(df$time)
    mid <- horizon/2
    df %>% filter(time %in% c(0, mid, horizon) | abs(time - mid) < 0.5) %>%
      mutate(Day = round(time/24, 1)) %>%
      select(scenario, Day, AEC, TISSUE, CARDIAC, SYMPTOM, ECP, TROP) %>%
      distinct(scenario, Day, .keep_all = TRUE) %>%
      arrange(scenario, Day) %>%
      datatable(rownames = FALSE, options = list(pageLength = 30))
  })
}

shinyApp(ui, server)
