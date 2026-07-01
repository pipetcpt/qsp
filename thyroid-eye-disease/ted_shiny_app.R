## ============================================================================
## Thyroid Eye Disease (TED / Graves' Orbitopathy) QSP — Shiny Dashboard
## ----------------------------------------------------------------------------
## 8 tabs: Patient profile · Drug PK · Immune/Fibroblast PD · Orbital tissue
##         remodeling · Clinical endpoints (CAS/Proptosis/GO-QoL) · Scenario
##         comparison · Biomarkers & safety · References
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
get_model <- function() {
  if (!exists(".TED_MOD", envir = .GlobalEnv) || is.null(.GlobalEnv$.TED_MOD)) {
    assign(".TED_MOD", mread_cache("ted", project = "."), envir = .GlobalEnv)
  }
  .GlobalEnv$.TED_MOD
}

SCENARIO_LIST <- c(
  "Untreated natural history",
  "IV methylprednisolone EUGOGO pulse",
  "Oral prednisone taper",
  "Teprotumumab OPTIC regimen",
  "Rituximab 1000mg d1,d15",
  "Tocilizumab 8mg/kg q4w x4",
  "Selenium 200mcg/day (mild TED)",
  "Teprotumumab + smoking cessation",
  "IVMP salvage -> Teprotumumab (sequential)"
)

## ---------- Scenario event builder (mAb doses are weight-based) ----------
build_events <- function(scenario, wt) {
  if (scenario == "Untreated natural history") {
    ev(amt = 0, cmt = "CEN_TEP", time = 0)

  } else if (scenario == "IV methylprednisolone EUGOGO pulse") {
    seq(ev(amt = 500, cmt = "GC_C",  ii = 168, addl = 5, time = 0),
        ev(amt = 500, cmt = "CUMGC", ii = 168, addl = 5, time = 0),
        ev(amt = 250, cmt = "GC_C",  ii = 168, addl = 5, time = 6*168),
        ev(amt = 250, cmt = "CUMGC", ii = 168, addl = 5, time = 6*168))

  } else if (scenario == "Oral prednisone taper") {
    doses <- c(60, 40, 30, 20, 15, 10)
    starts <- seq(0, by = 14*24, length.out = length(doses))
    ev_list <- lapply(seq_along(doses), function(i) {
      seq(ev(amt = doses[i], cmt = "GUT_GC", ii = 24, addl = 13, time = starts[i]),
          ev(amt = doses[i], cmt = "CUMGC",  ii = 24, addl = 13, time = starts[i]))
    })
    Reduce(seq, ev_list)

  } else if (scenario == "Teprotumumab OPTIC regimen") {
    seq(ev(amt = 10*wt, cmt = "CEN_TEP", time = 0),
        ev(amt = 20*wt, cmt = "CEN_TEP", ii = 21*24, addl = 6, time = 21*24))

  } else if (scenario == "Rituximab 1000mg d1,d15") {
    ev(amt = 1000, cmt = "CEN_RTX", ii = 14*24, addl = 1)

  } else if (scenario == "Tocilizumab 8mg/kg q4w x4") {
    ev(amt = 8*wt, cmt = "CEN_TCZ", ii = 28*24, addl = 3)

  } else if (scenario == "Selenium 200mcg/day (mild TED)") {
    ev(amt = 0, cmt = "CEN_TEP", time = 0)

  } else if (scenario == "Teprotumumab + smoking cessation") {
    seq(ev(amt = 10*wt, cmt = "CEN_TEP", time = 0),
        ev(amt = 20*wt, cmt = "CEN_TEP", ii = 21*24, addl = 6, time = 21*24))

  } else if (scenario == "IVMP salvage -> Teprotumumab (sequential)") {
    seq(ev(amt = 500,   cmt = "GC_C",   ii = 168, addl = 5, time = 0),
        ev(amt = 500,   cmt = "CUMGC",  ii = 168, addl = 5, time = 0),
        ev(amt = 10*wt, cmt = "CEN_TEP", time = 6*168),
        ev(amt = 20*wt, cmt = "CEN_TEP", ii = 21*24, addl = 6, time = 6*168 + 21*24))

  } else {
    ev(amt = 0, cmt = "CEN_TEP", time = 0)
  }
}

flags_for_scenario <- function(scenario) {
  list(
    selenium = as.integer(scenario == "Selenium 200mcg/day (mild TED)"),
    smoke_cessation = as.integer(scenario == "Teprotumumab + smoking cessation"),
    severity0 = if (scenario == "Selenium 200mcg/day (mild TED)") 0 else NULL
  )
}

run_sim <- function(scenario, horizon_h, params) {
  mod <- get_model()
  ev_set <- build_events(scenario, params$wt)
  fl <- flags_for_scenario(scenario)
  par <- list(
    AGE         = params$age,
    WT          = params$wt,
    SMOKING     = params$smoking,
    CAS0        = params$cas0,
    SEVERITY0   = if (!is.null(fl$severity0)) fl$severity0 else params$severity0,
    DURATION_MO = params$duration_mo,
    SELENIUM        = fl$selenium,
    SMOKE_CESSATION = fl$smoke_cessation
  )
  mod %>% param(par) %>% mrgsim(events = ev_set, end = horizon_h, delta = 24) %>%
    as_tibble() %>% mutate(scenario = scenario)
}

## ---------- UI ----------
ui <- dashboardPage(
  dashboardHeader(title = "TED QSP"),
  dashboardSidebar(
    sidebarMenu(
      menuItem("1. Patient profile",          tabName = "patient", icon = icon("user")),
      menuItem("2. Drug PK",                  tabName = "pk",      icon = icon("syringe")),
      menuItem("3. Immune / Fibroblast PD",   tabName = "pd",      icon = icon("dna")),
      menuItem("4. Orbital tissue remodeling",tabName = "remodel", icon = icon("eye")),
      menuItem("5. Clinical endpoints",       tabName = "clin",    icon = icon("clipboard-list")),
      menuItem("6. Scenario comparison",      tabName = "compare", icon = icon("layer-group")),
      menuItem("7. Biomarkers & safety",      tabName = "safety",  icon = icon("shield-alt")),
      menuItem("8. References",               tabName = "refs",    icon = icon("book"))
    ),
    hr(),
    selectInput("scenario", "Scenario:", SCENARIO_LIST, selected = "Teprotumumab OPTIC regimen"),
    sliderInput("horizon_mo", "Simulation horizon (months):", 3, 24, 12, step = 1),
    sliderInput("age",        "Age (years):",  18, 80, 45),
    sliderInput("wt",         "Weight (kg):",  40, 130, 70),
    selectInput("smoking",    "Smoking status:", choices = c("Non-smoker"=0, "Current smoker"=1), selected = 0),
    sliderInput("cas0",       "Baseline CAS (0-7):", 0, 7, 4),
    selectInput("severity0",  "Baseline EUGOGO severity:",
                choices = c("Mild"=0, "Moderate-severe"=1, "Sight-threatening"=2), selected = 1),
    sliderInput("duration_mo","Disease duration at baseline (months):", 0, 36, 6),
    actionButton("run", "Run simulation", icon = icon("play"),
                 style = "color:#fff;background:#1f6feb")
  ),
  dashboardBody(
    tabItems(
      tabItem("patient",
        fluidRow(
          box(width = 12, title = "Patient profile summary", status = "primary",
              solidHeader = TRUE, DTOutput("patient_table"), br(),
              p(strong("EUGOGO severity legend:"),
                "Mild = minor CAS/proptosis, tolerable diplopia; Moderate-severe = ",
                "sufficient impact on daily life to justify immunosuppression/",
                "teprotumumab; Sight-threatening = dysthyroid optic neuropathy (DON) ",
                "or corneal breakdown requiring urgent high-dose IVMP +/- decompression.")
          )
        ),
        fluidRow(
          box(width = 12, title = "About this dashboard", status = "info",
              p("Runs the mrgsolve QSP model for thyroid eye disease (TED/Graves' ",
                "orbitopathy): TSHR-TRAb autoimmunity cross-reacting with orbital ",
                "fibroblast IGF-1R drives T/B-cell infiltration, cytokine-mediated ",
                "fibroblast activation, hyaluronan/GAG edema, adipogenesis and late ",
                "fibrosis, producing proptosis, diplopia and (rarely) compressive ",
                "optic neuropathy."),
              p("Pick a treatment scenario and patient profile in the left panel, ",
                "then press ", strong("Run simulation"), ". Each tab focuses on a ",
                "model layer: PK, immune/fibroblast pharmacodynamics, orbital tissue ",
                "remodeling, clinical activity/proptosis/QoL endpoints, scenario ",
                "comparison, safety, and references.")
          )
        )
      ),

      tabItem("pk",
        fluidRow(
          box(width = 12, title = "Plasma mAb concentration & glucocorticoid signal",
              status = "primary", solidHeader = TRUE, plotOutput("pk_plot", height = 480))
        )
      ),

      tabItem("pd",
        fluidRow(
          box(width = 6, title = "Receptor occupancy (IGF-1R / IL-6R)", plotOutput("occ_plot", 360)),
          box(width = 6, title = "B-cell pool & TRAb titer", plotOutput("bcell_plot", 360))
        ),
        fluidRow(
          box(width = 6, title = "T-cell activation & cytokine drive", plotOutput("th_cyt_plot", 360)),
          box(width = 6, title = "Orbital fibroblast activation & hyaluronan", plotOutput("fib_ha_plot", 360))
        )
      ),

      tabItem("remodel",
        fluidRow(
          box(width = 6, title = "Orbital fat volume (mL-equivalent)", plotOutput("fat_plot", 360)),
          box(width = 6, title = "EOM volume/edema (mL-equivalent)", plotOutput("mus_plot", 360))
        ),
        fluidRow(
          box(width = 12, title = "Orbital fibrosis index (largely irreversible)", plotOutput("fib_idx_plot", 360))
        )
      ),

      tabItem("clin",
        fluidRow(
          box(width = 6, title = "Clinical Activity Score (CAS, 0-7)", plotOutput("cas_plot", 360)),
          box(width = 6, title = "Proptosis (Hertel, mm)", plotOutput("proptosis_plot", 360))
        ),
        fluidRow(
          box(width = 6, title = "GO-QoL score (0-100, higher=better)", plotOutput("qol_plot", 360)),
          box(width = 6, title = "EUGOGO 'active disease' flag (CAS>=3)", plotOutput("active_plot", 360))
        )
      ),

      tabItem("compare",
        fluidRow(
          box(width = 12, title = "Scenario comparison panel",
              status = "warning", solidHeader = TRUE,
              p("Runs all nine built-in scenarios with the current patient profile; ",
                "press the button below."),
              actionButton("run_all", "Run all scenarios", icon = icon("rocket"),
                           style = "color:#fff;background:#0f5132"),
              br(), br(),
              plotOutput("compare_plot", height = 600)
          )
        ),
        fluidRow(
          box(width = 12, title = "Endpoint summary table (baseline / month-3 / month-end)",
              status = "info", DTOutput("endpoint_table"))
        )
      ),

      tabItem("safety",
        fluidRow(
          box(width = 6, title = "Cumulative teprotumumab AUC & hearing-AE risk", plotOutput("hear_plot", 360)),
          box(width = 6, title = "Cumulative glucocorticoid dose & Cushingoid-AE risk", plotOutput("cush_plot", 360))
        ),
        fluidRow(
          box(width = 12, title = "Safety flags", status = "danger",
              DTOutput("safety_table"))
        )
      ),

      tabItem("refs",
        fluidRow(
          box(width = 12, title = "Key references", status = "primary", solidHeader = TRUE,
              p("See ", code("ted_references.md"), " in this directory for the full ",
                "curated list (30+ PubMed-linked references)."),
              tags$ul(
                tags$li("Douglas RS, et al. 2020 N Engl J Med — OPTIC trial (teprotumumab)."),
                tags$li("Kahaly GJ, et al. 2021 Lancet Diabetes Endocrinol — OPTIC-X extension."),
                tags$li("Bartalena L, et al. 2021 Eur Thyroid J — EUGOGO consensus statement."),
                tags$li("Salvi M, et al. 2015 J Clin Endocrinol Metab — rituximab, active disease."),
                tags$li("Stan MN, et al. 2015 J Clin Endocrinol Metab — rituximab RCT (no benefit)."),
                tags$li("Perez-Moreiras JV, et al. 2018 Am J Ophthalmol — tocilizumab RCT."),
                tags$li("Marcocci C, et al. 2011 N Engl J Med — selenium, mild TED."),
                tags$li("Bahn RS. 2010 N Engl J Med — TSHR-IGF-1R crosstalk mechanistic review.")
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

  cur_params <- function() list(
    age = input$age, wt = input$wt, smoking = as.numeric(input$smoking),
    cas0 = input$cas0, severity0 = as.numeric(input$severity0),
    duration_mo = input$duration_mo
  )

  observeEvent(input$run, {
    showNotification("Running mrgsolve simulation…", type = "message", duration = 1)
    results(run_sim(input$scenario, input$horizon_mo * 730, cur_params()))
  }, ignoreNULL = FALSE)

  observeEvent(input$run_all, {
    showNotification("Running 9 scenarios…", type = "message", duration = 1)
    out <- lapply(SCENARIO_LIST, function(sc) run_sim(sc, input$horizon_mo * 730, cur_params()))
    all_results(bind_rows(out))
  })

  output$patient_table <- renderDT({
    sev_lbl <- c("0"="Mild","1"="Moderate-severe","2"="Sight-threatening")[[as.character(input$severity0)]]
    smoke_lbl <- c("0"="Non-smoker","1"="Current smoker")[[as.character(input$smoking)]]
    tibble(
      Field = c("Age","Weight","Smoking","Baseline CAS","Baseline severity",
                "Disease duration (mo)","Scenario","Horizon (mo)"),
      Value = c(input$age, input$wt, smoke_lbl, input$cas0, sev_lbl,
                input$duration_mo, input$scenario, input$horizon_mo)
    ) %>% datatable(rownames = FALSE, options = list(dom = "t"))
  })

  output$pk_plot <- renderPlot({
    df <- results(); req(df)
    df %>% mutate(month = time/730) %>%
      select(month, conc_teprotumumab, conc_rituximab, conc_tocilizumab, GC_signal) %>%
      pivot_longer(-month) %>%
      ggplot(aes(month, value, colour = name)) +
        geom_line(linewidth = 0.8) +
        facet_wrap(~name, scales = "free_y") +
        labs(x = "Time (months)", y = "Concentration (mg/L) / GC signal", colour = "") +
        theme_minimal(base_size = 12) + theme(legend.position = "none")
  })

  output$occ_plot <- renderPlot({
    df <- results(); req(df)
    df %>% mutate(month = time/730) %>%
      select(month, IGF1R_occupancy, IL6R_occupancy) %>%
      pivot_longer(-month) %>%
      ggplot(aes(month, value, colour = name)) +
        geom_line(linewidth = 0.9) +
        scale_y_continuous(limits = c(0,1)) +
        labs(x = "Month", y = "Receptor occupancy (fraction)", colour = "") +
        theme_minimal(base_size = 13)
  })
  output$bcell_plot <- renderPlot({
    df <- results(); req(df)
    df %>% mutate(month = time/730) %>%
      select(month, BCell_idx, TRAb_idx) %>%
      pivot_longer(-month) %>%
      ggplot(aes(month, value, colour = name)) +
        geom_line(linewidth = 0.9) +
        geom_hline(yintercept = 1.0, lty = 2, colour = "grey50") +
        labs(x = "Month", y = "Relative index (1.0=baseline)", colour = "") +
        theme_minimal(base_size = 13)
  })
  output$th_cyt_plot <- renderPlot({
    df <- results(); req(df)
    df %>% mutate(month = time/730) %>%
      select(month, ThCell_idx, Cytokine_idx) %>%
      pivot_longer(-month) %>%
      ggplot(aes(month, value, colour = name)) +
        geom_line(linewidth = 0.9) +
        labs(x = "Month", y = "Relative index (1.0=baseline)", colour = "") +
        theme_minimal(base_size = 13)
  })
  output$fib_ha_plot <- renderPlot({
    df <- results(); req(df)
    df %>% mutate(month = time/730) %>%
      select(month, Fibroblast_idx, Hyaluronan_idx) %>%
      pivot_longer(-month) %>%
      ggplot(aes(month, value, colour = name)) +
        geom_line(linewidth = 0.9) +
        labs(x = "Month", y = "Relative index (1.0=baseline)", colour = "") +
        theme_minimal(base_size = 13)
  })

  output$fat_plot <- renderPlot({
    df <- results(); req(df)
    ggplot(df, aes(time/730, OrbitalFat_mL)) +
      geom_line(colour = "#e07a5f", linewidth = 1) +
      labs(x = "Month", y = "Orbital fat index (mL-equivalent)") +
      theme_minimal(base_size = 13)
  })
  output$mus_plot <- renderPlot({
    df <- results(); req(df)
    ggplot(df, aes(time/730, EOMVolume_mL)) +
      geom_line(colour = "#3d5a80", linewidth = 1) +
      labs(x = "Month", y = "EOM volume/edema index (mL-equivalent)") +
      theme_minimal(base_size = 13)
  })
  output$fib_idx_plot <- renderPlot({
    df <- results(); req(df)
    ggplot(df, aes(time/730, Fibrosis_idx)) +
      geom_line(colour = "#6d4c41", linewidth = 1) +
      scale_y_continuous(limits = c(0,1)) +
      labs(x = "Month", y = "Orbital fibrosis index (0-1)") +
      theme_minimal(base_size = 13)
  })

  output$cas_plot <- renderPlot({
    df <- results(); req(df)
    ggplot(df, aes(time/730, CAS_score)) +
      geom_line(colour = "#9b1c1c", linewidth = 1) +
      geom_hline(yintercept = 3, lty = 2, colour = "grey50") +
      scale_y_continuous(limits = c(0,7)) +
      labs(x = "Month", y = "CAS (0-7)", caption = "Dashed line = EUGOGO 'active disease' threshold (>=3)") +
      theme_minimal(base_size = 13)
  })
  output$proptosis_plot <- renderPlot({
    df <- results(); req(df)
    ggplot(df, aes(time/730, Proptosis_mm)) +
      geom_line(colour = "#264653", linewidth = 1) +
      labs(x = "Month", y = "Proptosis (Hertel, mm)") +
      theme_minimal(base_size = 13)
  })
  output$qol_plot <- renderPlot({
    df <- results(); req(df)
    ggplot(df, aes(time/730, GOQOL_score)) +
      geom_line(colour = "#2a9d8f", linewidth = 1) +
      scale_y_continuous(limits = c(0,100)) +
      labs(x = "Month", y = "GO-QoL score (0-100)") +
      theme_minimal(base_size = 13)
  })
  output$active_plot <- renderPlot({
    df <- results(); req(df)
    ggplot(df, aes(time/730, EUGOGO_active_flag)) +
      geom_step(colour = "#c0392b", linewidth = 1) +
      scale_y_continuous(limits = c(0,1), breaks = c(0,1)) +
      labs(x = "Month", y = "Active disease flag (CAS>=3)") +
      theme_minimal(base_size = 13)
  })

  output$compare_plot <- renderPlot({
    df <- all_results(); req(df)
    df %>% mutate(month = time/730) %>%
      ggplot(aes(month, Proptosis_mm, colour = scenario)) +
        geom_line(linewidth = 0.8) +
        labs(x = "Month", y = "Proptosis (mm)", colour = "Scenario") +
        theme_minimal(base_size = 13) +
        theme(legend.position = "bottom")
  })
  output$endpoint_table <- renderDT({
    df <- all_results(); req(df)
    horizon <- max(df$time)
    mo3 <- min(3*730, horizon)
    df %>% filter(time %in% c(0, mo3, horizon)) %>%
      mutate(Month = round(time/730, 1)) %>%
      select(scenario, Month, CAS_score, Proptosis_mm, GOQOL_score) %>%
      distinct(scenario, Month, .keep_all = TRUE) %>%
      arrange(scenario, Month) %>%
      datatable(rownames = FALSE, options = list(pageLength = 30))
  })

  output$hear_plot <- renderPlot({
    df <- results(); req(df)
    df %>% mutate(month = time/730) %>%
      select(month, AUC_TEP_mgh_L, HearingAE_risk) %>%
      pivot_longer(-month) %>%
      ggplot(aes(month, value, colour = name)) +
        geom_line(linewidth = 0.9) +
        facet_wrap(~name, scales = "free_y") +
        theme_minimal(base_size = 12) + theme(legend.position = "none")
  })
  output$cush_plot <- renderPlot({
    df <- results(); req(df)
    df %>% mutate(month = time/730) %>%
      select(month, CumGC_mg, CushingoidAE_risk) %>%
      pivot_longer(-month) %>%
      ggplot(aes(month, value, colour = name)) +
        geom_line(linewidth = 0.9) +
        facet_wrap(~name, scales = "free_y") +
        theme_minimal(base_size = 12) + theme(legend.position = "none")
  })
  output$safety_table <- renderDT({
    df <- results(); req(df)
    df %>% filter(time == max(time)) %>%
      select(AUC_TEP_mgh_L, HearingAE_risk, CumGC_mg, CushingoidAE_risk, HepatotoxFlag) %>%
      datatable(rownames = FALSE, options = list(dom = "t"))
  })
}

shinyApp(ui, server)
