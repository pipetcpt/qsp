## ============================================================================
## Vasomotor Symptoms of Menopause (VMS) QSP — Shiny Dashboard
## ----------------------------------------------------------------------------
## 8 tabs: Patient · Drug PK · HPG/KNDy axis · Thermoregulation · Hot-flush
##         endpoint · Sleep/Mood · Scenario comparison · Biomarkers/Safety
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
MODEL_PATH <- "vms_mrgsolve_model.R"
get_model <- function() {
  if (!exists(".VMS_MOD", envir = .GlobalEnv) || is.null(.GlobalEnv$.VMS_MOD)) {
    assign(".VMS_MOD", mread_cache("vms", project = "."), envir = .GlobalEnv)
  }
  .GlobalEnv$.VMS_MOD
}

## ---------- Scenario builders ----------
build_events <- function(scenario, horizon) {
  if (scenario == "Untreated") {
    ev(amt = 0, cmt = "CEN_FEZ")
  } else if (scenario == "Fezolinetant 45 mg QD") {
    ev(amt = 45, cmt = "GUT_FEZ", ii = 1, addl = horizon)
  } else if (scenario == "Fezolinetant 30 mg QD") {
    ev(amt = 30, cmt = "GUT_FEZ", ii = 1, addl = horizon)
  } else if (scenario == "Elinzanetant 120 mg QD") {
    ev(amt = 120, cmt = "GUT_ELZ", ii = 1, addl = horizon)
  } else if (scenario == "Oral estradiol 1 mg QD") {
    ev(amt = 1, cmt = "GUT_E2O", ii = 1, addl = horizon)
  } else if (scenario == "Transdermal E2 patch 0.05 mg/d") {
    ev(amt = 3.5, cmt = "RES_E2P", ii = 7, addl = ceiling(horizon / 7))
  } else if (scenario == "Paroxetine 7.5 mg QD") {
    ev(amt = 7.5, cmt = "GUT_PRX", ii = 1, addl = horizon)
  } else if (scenario == "Venlafaxine ER 75 mg QD") {
    ev(amt = 75, cmt = "GUT_VEN", ii = 1, addl = horizon)
  } else if (scenario == "Gabapentin 900 mg/d (TID)") {
    ev(amt = 300, cmt = "GUT_GBP", ii = 0.333, addl = 3 * horizon)
  } else if (scenario == "Clonidine 0.1 mg QD") {
    ev(amt = 0.1, cmt = "GUT_CLN", ii = 1, addl = horizon)
  } else {
    ev(amt = 0, cmt = "CEN_FEZ")
  }
}

run_sim <- function(scenario, horizon, params) {
  mod <- get_model()
  ev_set <- build_events(scenario, horizon)
  par <- list(
    AGE        = params$age,
    WT         = params$wt,
    BMI        = params$bmi,
    MENO_STAGE = params$meno_stage,
    SURG_MENO  = ifelse(params$surg_meno, 1, 0),
    BASE_HF    = params$base_hf,
    BASE_SEV   = params$base_sev,
    BASE_BMD   = params$base_bmd,
    BASE_LDL   = params$base_ldl,
    BASE_PSQI  = params$base_psqi
  )
  mod %>% param(par) %>% mrgsim(events = ev_set, end = horizon, delta = 1) %>%
    as_tibble() %>% mutate(scenario = scenario)
}

SCENARIO_LIST <- c(
  "Untreated", "Fezolinetant 45 mg QD", "Fezolinetant 30 mg QD",
  "Elinzanetant 120 mg QD", "Oral estradiol 1 mg QD",
  "Transdermal E2 patch 0.05 mg/d", "Paroxetine 7.5 mg QD",
  "Venlafaxine ER 75 mg QD", "Gabapentin 900 mg/d (TID)",
  "Clonidine 0.1 mg QD"
)

## ---------- UI ----------
ui <- dashboardPage(
  dashboardHeader(title = "VMS (Menopause) QSP"),
  dashboardSidebar(
    sidebarMenu(
      menuItem("1. Patient profile",     tabName = "patient", icon = icon("user")),
      menuItem("2. Drug PK",             tabName = "pk",      icon = icon("syringe")),
      menuItem("3. HPG / KNDy axis",     tabName = "kndy",    icon = icon("brain")),
      menuItem("4. Thermoregulation",    tabName = "thermo",  icon = icon("temperature-high")),
      menuItem("5. Hot-flush endpoint",  tabName = "hf",      icon = icon("fire")),
      menuItem("6. Sleep / Mood",        tabName = "sleep",   icon = icon("bed")),
      menuItem("7. Scenario comparison", tabName = "compare", icon = icon("chart-line")),
      menuItem("8. Biomarkers / Safety", tabName = "safety",  icon = icon("shield-alt"))
    ),
    hr(),
    selectInput("scenario", "Scenario:", SCENARIO_LIST, selected = "Fezolinetant 45 mg QD"),
    sliderInput("horizon",   "Simulation horizon (days):", 28, 730, 168, step = 7),
    sliderInput("age",       "Age (years):",                40,  60,  51),
    sliderInput("wt",        "Weight (kg):",                40, 120,  65),
    sliderInput("bmi",       "BMI:",                        16,  45,  26),
    selectInput("meno_stage","Menopausal stage:",
                choices = c("Late perimenopause"=1, "Early postmenopause"=2,
                            "Late postmenopause"=3), selected = 2),
    checkboxInput("surg_meno","Surgical/abrupt menopause (BSO)", FALSE),
    sliderInput("base_hf",   "Baseline HF/day (moderate-severe):", 2, 15, 8),
    sliderInput("base_sev",  "Baseline HF severity (0-3):",       0,  3, 2.0, step=0.1),
    sliderInput("base_bmd",  "Baseline lumbar BMD (g/cm^2):",   0.7,1.3, 1.00, step = 0.01),
    sliderInput("base_ldl",  "Baseline LDL-C (mg/dL):",          70, 200, 120),
    sliderInput("base_psqi", "Baseline PSQI sleep score:",        0,  21,  10),
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
              p(strong("Menopausal-stage legend:"),
                "Late perimenopause = irregular cycles, rising FSH; ",
                "Early postmenopause = <5-6 y since FMP, highest VMS burden; ",
                "Late postmenopause = >6 y since FMP, VMS may persist (SWAN median 7.4 y).")
          )
        ),
        fluidRow(
          box(width = 12, title = "About this dashboard", status = "info",
              p("This dashboard runs the mrgsolve QSP model for menopausal",
                "vasomotor symptoms (VMS). Pick a treatment scenario in the left",
                "panel, adjust the patient profile, then press ",
                strong("Run simulation"), " to update plots."),
              p("Each tab focuses on a layer of the model: PK, HPG-axis/KNDy",
                "neuron pharmacodynamics, thermoneutral-zone narrowing,",
                "hot-flush frequency/severity, sleep & mood, scenario",
                "comparison, and safety (BMD, LDL-C).")
          )
        )
      ),

      tabItem("pk",
        fluidRow(
          box(width = 12, title = "Plasma drug exposure (selected scenario)",
              status = "primary", solidHeader = TRUE, plotOutput("pk_plot", height = 480))
        )
      ),

      tabItem("kndy",
        fluidRow(
          box(width = 6, title = "Follicle pool & Estradiol", plotOutput("foll_e2", 360)),
          box(width = 6, title = "FSH & KNDy/NKB tone", plotOutput("fsh_kndy", 360))
        )
      ),

      tabItem("thermo",
        fluidRow(
          box(width = 12, title = "Thermoneutral-zone (TNZ) half-width (°C)",
              status = "primary", solidHeader = TRUE,
              plotOutput("tnz_plot", 380))
        )
      ),

      tabItem("hf",
        fluidRow(
          box(width = 6, title = "Hot-flush frequency (events/day)",
              plotOutput("hf_freq_plot", 360)),
          box(width = 6, title = "Hot-flush severity (0-3)",
              plotOutput("hf_sev_plot", 360))
        )
      ),

      tabItem("sleep",
        fluidRow(
          box(width = 6, title = "PSQI sleep-quality index",
              plotOutput("psqi_plot", 360)),
          box(width = 6, title = "Mood / depressive-symptom composite",
              plotOutput("mood_plot", 360))
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
        )
      ),

      tabItem("safety",
        fluidRow(
          box(width = 6, title = "Bone (CTX resorption marker & lumbar BMD)",
              plotOutput("bone_plot", 360)),
          box(width = 6, title = "LDL cholesterol (mg/dL)",
              plotOutput("ldl_plot",  360))
        ),
        fluidRow(
          box(width = 12, title = "Endpoint summary table (Day-0, Day-84, Day-end)",
              status = "info",
              DTOutput("endpoint_table"))
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
      age = input$age, wt = input$wt, bmi = input$bmi,
      meno_stage = as.numeric(input$meno_stage), surg_meno = input$surg_meno,
      base_hf = input$base_hf, base_sev = input$base_sev,
      base_bmd = input$base_bmd, base_ldl = input$base_ldl,
      base_psqi = input$base_psqi
    )
    results(run_sim(input$scenario, input$horizon, p))
  }, ignoreNULL = FALSE)

  observeEvent(input$run_all, {
    showNotification("Running 10 scenarios…", type = "message", duration = 1)
    p <- list(
      age = input$age, wt = input$wt, bmi = input$bmi,
      meno_stage = as.numeric(input$meno_stage), surg_meno = input$surg_meno,
      base_hf = input$base_hf, base_sev = input$base_sev,
      base_bmd = input$base_bmd, base_ldl = input$base_ldl,
      base_psqi = input$base_psqi
    )
    out <- lapply(SCENARIO_LIST, function(sc) run_sim(sc, input$horizon, p))
    all_results(bind_rows(out))
  })

  # --- Patient table ---
  output$patient_table <- renderDT({
    tibble(
      Field = c("Age", "Weight", "BMI", "Menopausal stage", "Surgical menopause",
                "Baseline HF/day", "Baseline HF severity",
                "Baseline BMD (g/cm²)", "Baseline LDL-C (mg/dL)",
                "Baseline PSQI", "Scenario", "Horizon (d)"),
      Value = c(input$age, input$wt, input$bmi,
                names(which(c("1"="Late perimenopause","2"="Early postmenopause",
                               "3"="Late postmenopause") == as.character(input$meno_stage))),
                input$surg_meno, input$base_hf, input$base_sev,
                input$base_bmd, input$base_ldl, input$base_psqi,
                input$scenario, input$horizon)
    ) %>% datatable(rownames = FALSE, options = list(dom = "t"))
  })

  # --- PK ---
  output$pk_plot <- renderPlot({
    df <- results(); req(df)
    df %>% select(time, conc_fezolinetant, conc_elinzanetant, conc_E2_oral,
                  conc_E2_patch, conc_paroxetine, conc_venlafaxine, conc_ODV,
                  conc_gabapentin, conc_clonidine) %>%
      pivot_longer(-time) %>%
      filter(value > 0 | name == "conc_fezolinetant") %>%
      ggplot(aes(time, value, colour = name)) +
        geom_line(linewidth = 0.7) +
        scale_y_continuous(trans = "log1p") +
        labs(x = "Time (days)", y = "Plasma conc (mg/L)", colour = "Drug") +
        theme_minimal(base_size = 13)
  })

  # --- KNDy / HPG ---
  output$foll_e2 <- renderPlot({
    df <- results(); req(df)
    df %>% select(time, Follicle_pool, E2_pgmL) %>%
      pivot_longer(-time) %>%
      ggplot(aes(time, value, colour = name)) +
        geom_line(linewidth = 0.8) +
        facet_wrap(~name, scales = "free_y") +
        theme_minimal(base_size = 12) + theme(legend.position = "none")
  })
  output$fsh_kndy <- renderPlot({
    df <- results(); req(df)
    df %>% select(time, FSH_IUL, KNDy_tone) %>%
      pivot_longer(-time) %>%
      ggplot(aes(time, value, colour = name)) +
        geom_line(linewidth = 0.8) +
        facet_wrap(~name, scales = "free_y") +
        theme_minimal(base_size = 12) + theme(legend.position = "none")
  })

  # --- Thermoregulation ---
  output$tnz_plot <- renderPlot({
    df <- results(); req(df)
    ggplot(df, aes(time, TNZ_degC)) +
      geom_line(colour = "#e07a5f", linewidth = 1) +
      geom_hline(yintercept = 2.8, lty = 2, colour = "grey50") +
      labs(x = "Day", y = "Thermoneutral-zone half-width (°C)",
           caption = "Dashed line = healthy premenopausal reference (~2.8°C)") +
      theme_minimal(base_size = 13)
  })

  # --- Hot flush ---
  output$hf_freq_plot <- renderPlot({
    df <- results(); req(df)
    ggplot(df, aes(time, HF_freq_perday)) +
      geom_line(colour = "#cc3344", linewidth = 1) +
      labs(x = "Day", y = "Hot flushes / day") +
      theme_minimal(base_size = 13)
  })
  output$hf_sev_plot <- renderPlot({
    df <- results(); req(df)
    ggplot(df, aes(time, HF_severity)) +
      geom_line(colour = "#9b1c1c", linewidth = 1) +
      scale_y_continuous(limits = c(0, 3)) +
      labs(x = "Day", y = "HF severity (0-3)") +
      theme_minimal(base_size = 13)
  })

  # --- Sleep / Mood ---
  output$psqi_plot <- renderPlot({
    df <- results(); req(df)
    ggplot(df, aes(time, PSQI_sleep)) +
      geom_line(colour = "#264653", linewidth = 1) +
      scale_y_continuous(limits = c(0, 21)) +
      labs(x = "Day", y = "PSQI score (higher = worse sleep)") +
      theme_minimal(base_size = 13)
  })
  output$mood_plot <- renderPlot({
    df <- results(); req(df)
    ggplot(df, aes(time, Mood_score)) +
      geom_line(colour = "#6d4c41", linewidth = 1) +
      labs(x = "Day", y = "Mood composite (0-27, higher = worse)") +
      theme_minimal(base_size = 13)
  })

  # --- Compare ---
  output$compare_plot <- renderPlot({
    df <- all_results(); req(df)
    df %>% select(time, scenario, HF_freq_perday, TNZ_degC, PSQI_sleep, E2_pgmL) %>%
      pivot_longer(-c(time, scenario)) %>%
      ggplot(aes(time, value, colour = scenario)) +
        geom_line(linewidth = 0.6, alpha = 0.85) +
        facet_wrap(~name, scales = "free_y") +
        labs(x = "Day", y = NULL, colour = NULL) +
        theme_minimal(base_size = 12) +
        theme(legend.position = "bottom")
  })

  # --- Safety ---
  output$bone_plot <- renderPlot({
    df <- results(); req(df)
    df %>% select(time, CTX_resorption, BMD_gcm2) %>%
      pivot_longer(-time) %>%
      ggplot(aes(time, value, colour = name)) +
        geom_line(linewidth = 0.8) +
        facet_wrap(~name, scales = "free_y") +
        theme_minimal(base_size = 12) + theme(legend.position = "none")
  })
  output$ldl_plot <- renderPlot({
    df <- results(); req(df)
    ggplot(df, aes(time, LDL_mgdL)) +
      geom_line(colour = "#8a6d3b", linewidth = 1) +
      labs(x = "Day", y = "LDL cholesterol (mg/dL)") +
      theme_minimal(base_size = 13)
  })

  output$endpoint_table <- renderDT({
    df <- results(); req(df)
    tp <- c(0, 84, max(df$time))
    df %>% filter(time %in% tp) %>%
      select(Day = time, HF_freq_perday, HF_severity, TNZ_degC, PSQI_sleep,
             E2_pgmL, BMD_gcm2, LDL_mgdL) %>%
      mutate(across(-Day, ~ round(.x, 2))) %>%
      datatable(rownames = FALSE,
                options = list(dom = "t",
                               columnDefs = list(list(className='dt-center', targets="_all"))))
  })
}

shinyApp(ui, server)
