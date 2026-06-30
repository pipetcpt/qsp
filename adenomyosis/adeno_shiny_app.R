## ============================================================================
## Adenomyosis QSP — Shiny Dashboard
## ----------------------------------------------------------------------------
## 8 tabs: Patient · Drug PK · HPO axis · Lesion biology · Pain endpoint ·
##         Bleeding endpoint · Scenario comparison · Biomarkers/Safety
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
MODEL_PATH <- "adeno_mrgsolve_model.R"
get_model <- function() {
  if (!exists(".ADENO_MOD", envir = .GlobalEnv) || is.null(.GlobalEnv$.ADENO_MOD)) {
    assign(".ADENO_MOD", mread_cache("adeno", project = "."), envir = .GlobalEnv)
  }
  .GlobalEnv$.ADENO_MOD
}

## ---------- Scenario builders ----------
build_events <- function(scenario, horizon, addback_flag, sprm_flag) {
  ev_combined <- NULL
  add <- function(x) {
    ev_combined <<- if (is.null(ev_combined)) x else c(ev_combined, x)
  }
  if (scenario == "Untreated") {
    add(ev(amt = 0, cmt = "CEN_REL"))
  } else if (scenario == "NSAID + TXA") {
    add(ev(amt = 400, cmt = "GUT_IBU", ii = 0.25, addl = floor(horizon / 0.25)))
    menses_starts <- seq(0, horizon, 28)
    for (m in menses_starts)
      add(ev(amt = 1000, cmt = "GUT_TXA", time = m, ii = 0.25, addl = 19))
  } else if (scenario == "Dienogest 2 mg QD") {
    add(ev(amt = 2, cmt = "GUT_DNG", ii = 1, addl = horizon))
  } else if (scenario == "LNG-IUS 52 mg") {
    add(ev(amt = 52, cmt = "RES_LNG", time = 0))
  } else if (scenario == "Leuprolide 3.75 mg q4w") {
    add(ev(amt = 3.75, cmt = "DEP_LEU", ii = 28, addl = ceiling(horizon / 28)))
  } else if (scenario == "Relugolix-CT 40 mg QD") {
    add(ev(amt = 40, cmt = "GUT_REL", ii = 1, addl = horizon))
  } else if (scenario == "Elagolix 200 mg BID") {
    add(ev(amt = 200, cmt = "GUT_ELA", ii = 0.5, addl = floor(horizon / 0.5)))
  } else if (scenario == "Letrozole 2.5 mg QD") {
    add(ev(amt = 2.5, cmt = "GUT_LET", ii = 1, addl = horizon))
  } else if (scenario == "Hysterectomy") {
    add(ev(amt = 0, cmt = "CEN_REL"))
  } else if (scenario == "UAE") {
    add(ev(amt = 0, cmt = "CEN_REL"))
  } else {
    add(ev(amt = 0, cmt = "CEN_REL"))
  }
  ev_combined
}

run_sim <- function(scenario, horizon, params) {
  mod <- get_model()
  ev_set <- build_events(scenario, horizon, params$addback, params$sprm)
  par <- list(
    AGE        = params$age,
    WT         = params$wt,
    BMI        = params$bmi,
    PARITY     = params$parity,
    BASE_JZ    = params$base_jz,
    BASE_PBAC  = params$base_pbac,
    BASE_VAS   = params$base_vas,
    BASE_HB    = params$base_hb,
    BASE_BMD   = params$base_bmd,
    ADENO_SEV  = params$severity,
    ADDBACK_ON = ifelse(params$addback, 1, 0),
    SPRM_ON    = ifelse(params$sprm,   1, 0),
    HYS_ON     = ifelse(scenario == "Hysterectomy", 1, 0),
    UAE_ON     = ifelse(scenario == "UAE",          1, 0),
    HIFU_ON    = 0,
    ABLAT_ON   = 0
  )
  mod %>% param(par) %>% mrgsim(events = ev_set, end = horizon, delta = 1) %>%
    as_tibble() %>% mutate(scenario = scenario)
}

SCENARIO_LIST <- c(
  "Untreated", "NSAID + TXA", "Dienogest 2 mg QD",
  "LNG-IUS 52 mg", "Leuprolide 3.75 mg q4w",
  "Relugolix-CT 40 mg QD", "Elagolix 200 mg BID",
  "Letrozole 2.5 mg QD", "UAE", "Hysterectomy"
)

## ---------- UI ----------
ui <- dashboardPage(
  dashboardHeader(title = "Adenomyosis QSP"),
  dashboardSidebar(
    sidebarMenu(
      menuItem("1. Patient profile",      tabName = "patient", icon = icon("user")),
      menuItem("2. Drug PK",              tabName = "pk",      icon = icon("syringe")),
      menuItem("3. HPO axis PD",          tabName = "hpo",     icon = icon("brain")),
      menuItem("4. Lesion biology",       tabName = "lesion",  icon = icon("dna")),
      menuItem("5. Pain endpoint",        tabName = "pain",    icon = icon("bolt")),
      menuItem("6. HMB / bleeding",       tabName = "bleed",   icon = icon("tint")),
      menuItem("7. Scenario comparison",  tabName = "compare", icon = icon("chart-line")),
      menuItem("8. Biomarkers / Safety",  tabName = "safety",  icon = icon("shield-alt"))
    ),
    hr(),
    selectInput("scenario", "Scenario:", SCENARIO_LIST, selected = "Dienogest 2 mg QD"),
    sliderInput("horizon",  "Simulation horizon (days):", 90, 720, 365, step = 30),
    sliderInput("age",      "Age (years):",                18,  55,  35),
    sliderInput("wt",       "Weight (kg):",                40,  120, 65),
    sliderInput("bmi",      "BMI:",                        16,  45,  24),
    sliderInput("parity",   "Parity:",                      0,   5,  1),
    sliderInput("base_jz",  "Baseline JZ thickness (mm):",  6,  22, 12),
    sliderInput("base_pbac","Baseline PBAC:",             100, 600, 250),
    sliderInput("base_vas", "Baseline dysmenorrhea VAS:",   0,  10,   7),
    sliderInput("base_hb",  "Baseline Hb (g/dL):",          7,  15,  11),
    sliderInput("base_bmd", "Baseline BMD (g/cm²):",      0.7, 1.3,  1.05, step = 0.01),
    sliderInput("severity", "Disease severity (1-3):",      1,   3,   2),
    checkboxInput("addback", "Hormonal add-back (E2 1 mg / NETA 0.5 mg)", FALSE),
    checkboxInput("sprm",    "SPRM modifier (mifepristone/ulipristal)",   FALSE),
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
              p(strong("Severity legend:"),
                "1 = mild (JZ 8–10 mm, PBAC ~150, VAS 3–4); ",
                "2 = moderate (JZ 11–14 mm, PBAC 200–300, VAS 5–7); ",
                "3 = severe (JZ ≥15 mm, PBAC ≥300, VAS 8–10).")
          )
        ),
        fluidRow(
          box(width = 12, title = "About this dashboard", status = "info",
              p("This dashboard runs the mrgsolve QSP model for adenomyosis.",
                "Pick a treatment scenario in the left panel, adjust the patient",
                "profile, then press ", strong("Run simulation"), " to update plots."),
              p("Each tab focuses on a layer of the model: PK, HPO-axis pharmacodynamics,",
                "local lesion biology (aromatase/PGE2/TGF-β), pain endpoint (VAS),",
                "bleeding endpoint (PBAC + Hb), scenario comparison, and safety",
                "(BMD, E2, EHP-30 QoL).")
          )
        )
      ),

      tabItem("pk",
        fluidRow(
          box(width = 12, title = "Plasma drug exposure (selected scenario)",
              status = "primary", solidHeader = TRUE, plotOutput("pk_plot", height = 480))
        )
      ),

      tabItem("hpo",
        fluidRow(
          box(width = 6, title = "Gonadotropins (FSH, LH)", plotOutput("fsh_lh", 360)),
          box(width = 6, title = "Estradiol & Progesterone", plotOutput("e2_p4", 360))
        )
      ),

      tabItem("lesion",
        fluidRow(
          box(width = 6, title = "Aromatase · COX-2/PGE2 · TGF-β",
              plotOutput("lesion_PD", 360)),
          box(width = 6, title = "Local E2, fibrosis, VEGF, NGF",
              plotOutput("local_E2", 360))
        )
      ),

      tabItem("pain",
        fluidRow(
          box(width = 12, title = "Dysmenorrhea VAS (0-10)",
              status = "primary", solidHeader = TRUE,
              plotOutput("vas_plot", 380))
        ),
        fluidRow(
          box(width = 12, title = "EHP-30 quality-of-life score (lower = better)",
              plotOutput("qol_plot", 320))
        )
      ),

      tabItem("bleed",
        fluidRow(
          box(width = 6, title = "PBAC heavy menstrual bleeding",
              plotOutput("pbac_plot", 360)),
          box(width = 6, title = "Hemoglobin trajectory",
              plotOutput("hb_plot", 360))
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
          box(width = 6, title = "BMD trajectory (lumbar, g/cm²)",
              plotOutput("bmd_plot", 360)),
          box(width = 6, title = "Junctional-zone thickness (mm)",
              plotOutput("jz_plot",  360))
        ),
        fluidRow(
          box(width = 12, title = "Endpoint summary table (Day-0, Day-180, Day-end)",
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
      parity = input$parity,
      base_jz = input$base_jz, base_pbac = input$base_pbac,
      base_vas = input$base_vas, base_hb = input$base_hb,
      base_bmd = input$base_bmd, severity = input$severity,
      addback = input$addback, sprm = input$sprm
    )
    results(run_sim(input$scenario, input$horizon, p))
  })

  observeEvent(input$run_all, {
    showNotification("Running 10 scenarios…", type = "message", duration = 1)
    p <- list(
      age = input$age, wt = input$wt, bmi = input$bmi,
      parity = input$parity,
      base_jz = input$base_jz, base_pbac = input$base_pbac,
      base_vas = input$base_vas, base_hb = input$base_hb,
      base_bmd = input$base_bmd, severity = input$severity,
      addback = input$addback, sprm = input$sprm
    )
    out <- lapply(SCENARIO_LIST, function(sc) run_sim(sc, input$horizon, p))
    all_results(bind_rows(out))
  })

  # --- Patient table ---
  output$patient_table <- renderDT({
    tibble(
      Field = c("Age", "Weight", "BMI", "Parity",
                "Baseline JZ (mm)", "Baseline PBAC", "Baseline VAS",
                "Baseline Hb (g/dL)", "Baseline BMD (g/cm²)", "Severity",
                "Add-back", "SPRM", "Scenario", "Horizon (d)"),
      Value = c(input$age, input$wt, input$bmi, input$parity,
                input$base_jz, input$base_pbac, input$base_vas,
                input$base_hb, input$base_bmd, input$severity,
                input$addback, input$sprm, input$scenario, input$horizon)
    ) %>% datatable(rownames = FALSE, options = list(dom = "t"))
  })

  # --- PK ---
  output$pk_plot <- renderPlot({
    df <- results(); req(df)
    df %>% select(time, conc_relugolix, conc_elagolix, conc_leuprolide,
                  conc_dienogest, conc_LNG, conc_letrozole,
                  conc_ibuprofen, conc_TXA) %>%
      pivot_longer(-time) %>%
      filter(value > 0 | name == "conc_relugolix") %>%
      ggplot(aes(time, value, colour = name)) +
        geom_line(linewidth = 0.7) +
        scale_y_continuous(trans = "log1p") +
        labs(x = "Time (days)", y = "Plasma conc (mg/L)", colour = "Drug") +
        theme_minimal(base_size = 13)
  })

  # --- HPO ---
  output$fsh_lh <- renderPlot({
    df <- results(); req(df)
    df %>% select(time, FSH_IUL, LH_IUL) %>%
      pivot_longer(-time) %>%
      ggplot(aes(time, value, colour = name)) +
        geom_line(linewidth = 0.8) +
        labs(x = "Day", y = "IU/L", colour = NULL) +
        theme_minimal(base_size = 13)
  })
  output$e2_p4 <- renderPlot({
    df <- results(); req(df)
    df %>% select(time, E2_pgmL, P4_ngmL) %>%
      pivot_longer(-time) %>%
      ggplot(aes(time, value, colour = name)) +
        geom_line(linewidth = 0.8) +
        labs(x = "Day", y = "Conc.", colour = NULL) +
        theme_minimal(base_size = 13)
  })

  # --- Lesion ---
  output$lesion_PD <- renderPlot({
    df <- results(); req(df)
    df %>% select(time, Aromatase, COX2_PGE2, TGFb) %>%
      pivot_longer(-time) %>%
      ggplot(aes(time, value, colour = name)) +
        geom_line(linewidth = 0.8) +
        labs(x = "Day", y = "Relative activity", colour = NULL) +
        theme_minimal(base_size = 13)
  })
  output$local_E2 <- renderPlot({
    df <- results(); req(df)
    df %>% select(time, LocalE2_pgmL, Fibrosis_score, VEGF_score, NGF_score) %>%
      pivot_longer(-time) %>%
      ggplot(aes(time, value, colour = name)) +
        geom_line(linewidth = 0.8) +
        facet_wrap(~name, scales = "free_y") +
        theme_minimal(base_size = 12) + theme(legend.position = "none")
  })

  # --- Pain ---
  output$vas_plot <- renderPlot({
    df <- results(); req(df)
    ggplot(df, aes(time, VAS_pain)) +
      geom_line(colour = "#cc3344", linewidth = 1) +
      scale_y_continuous(limits = c(0,10)) +
      labs(x = "Day", y = "Dysmenorrhea VAS") +
      theme_minimal(base_size = 13)
  })
  output$qol_plot <- renderPlot({
    df <- results(); req(df)
    ggplot(df, aes(time, EHP30_QoL)) +
      geom_line(colour = "#0b7285", linewidth = 1) +
      labs(x = "Day", y = "EHP-30 score (lower=better)") +
      theme_minimal(base_size = 13)
  })

  # --- Bleeding ---
  output$pbac_plot <- renderPlot({
    df <- results(); req(df)
    ggplot(df, aes(time, PBAC_score)) +
      geom_line(colour = "#9b1c1c", linewidth = 1) +
      geom_hline(yintercept = 100, lty = 2, colour = "grey50") +
      labs(x = "Day", y = "PBAC score (≥100 = HMB)") +
      theme_minimal(base_size = 13)
  })
  output$hb_plot <- renderPlot({
    df <- results(); req(df)
    ggplot(df, aes(time, Hb_gdL)) +
      geom_line(colour = "#0f5132", linewidth = 1) +
      geom_hline(yintercept = c(12, 8), lty = 2, colour = c("grey50","red")) +
      labs(x = "Day", y = "Hemoglobin (g/dL)") +
      theme_minimal(base_size = 13)
  })

  # --- Compare ---
  output$compare_plot <- renderPlot({
    df <- all_results(); req(df)
    df %>% select(time, scenario, VAS_pain, PBAC_score, JZ_mm, E2_pgmL) %>%
      pivot_longer(-c(time, scenario)) %>%
      ggplot(aes(time, value, colour = scenario)) +
        geom_line(linewidth = 0.6, alpha = 0.85) +
        facet_wrap(~name, scales = "free_y") +
        labs(x = "Day", y = NULL, colour = NULL) +
        theme_minimal(base_size = 12) +
        theme(legend.position = "bottom")
  })

  # --- Safety ---
  output$bmd_plot <- renderPlot({
    df <- results(); req(df)
    ggplot(df, aes(time, BMD_gcm2)) +
      geom_line(colour = "#8a6d3b", linewidth = 1) +
      labs(x = "Day", y = "Lumbar BMD (g/cm²)") +
      theme_minimal(base_size = 13)
  })
  output$jz_plot <- renderPlot({
    df <- results(); req(df)
    ggplot(df, aes(time, JZ_mm)) +
      geom_line(colour = "#264653", linewidth = 1) +
      geom_hline(yintercept = c(8, 12), lty = 2, colour = "grey60") +
      labs(x = "Day", y = "Junctional-zone thickness (mm)") +
      theme_minimal(base_size = 13)
  })

  output$endpoint_table <- renderDT({
    df <- results(); req(df)
    tp <- c(0, 180, max(df$time))
    df %>% filter(time %in% tp) %>%
      select(Day = time, VAS_pain, PBAC_score, JZ_mm, Hb_gdL,
             E2_pgmL, BMD_gcm2, EHP30_QoL) %>%
      mutate(across(-Day, ~ round(.x, 2))) %>%
      datatable(rownames = FALSE,
                options = list(dom = "t",
                               columnDefs = list(list(className='dt-center', targets="_all"))))
  })
}

shinyApp(ui, server)
