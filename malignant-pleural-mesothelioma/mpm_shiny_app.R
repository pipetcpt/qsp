## ============================================================================
## Malignant Pleural Mesothelioma (MPM) QSP — Shiny Dashboard
## ----------------------------------------------------------------------------
## 8 tabs: Patient · Drug PK · Molecular pathway · Immune/TME & Effusion ·
##         Clinical endpoints · Scenario comparison · Safety · Biomarkers
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
MODEL_PATH <- "mpm_mrgsolve_model.R"
get_model <- function() {
  if (!exists(".MPM_MOD", envir = .GlobalEnv) || is.null(.GlobalEnv$.MPM_MOD)) {
    assign(".MPM_MOD", mread_cache("mpm", project = "."), envir = .GlobalEnv)
  }
  .GlobalEnv$.MPM_MOD
}

## ---------- Scenario builders ----------
build_events <- function(scenario, horizon) {
  if (scenario == "Untreated / BSC") {
    ev(amt = 0, cmt = "CIS_C1")
  } else if (scenario == "Cisplatin+Pemetrexed (EMPHACIS)") {
    seq(ev(amt = 135, cmt = "CIS_C1", ii = 21, addl = 5),
        ev(amt = 900, cmt = "PEM_C1", ii = 21, addl = 5))
  } else if (scenario == "Cisplatin+Pemetrexed+Bevacizumab (MAPS)") {
    seq(ev(amt = 135, cmt = "CIS_C1", ii = 21, addl = 5),
        ev(amt = 900, cmt = "PEM_C1", ii = 21, addl = 5),
        ev(amt = 1050, cmt = "BEV_C1", ii = 21, addl = 11))
  } else if (scenario == "Nivolumab+Ipilimumab 1L (CheckMate 743)") {
    seq(ev(amt = 240, cmt = "NIVO_C1", ii = 14, addl = 25),
        ev(amt = 63,  cmt = "IPI_C1",  ii = 42, addl = 8))
  } else if (scenario == "Nivolumab monotherapy 2L (CONFIRM)") {
    ev(amt = 240, cmt = "NIVO_C1", ii = 14, addl = 25)
  } else if (scenario == "Nivolumab+Ipilimumab 2L (MAPS2)") {
    seq(ev(amt = 180, cmt = "NIVO_C1", ii = 14, addl = 25),
        ev(amt = 63,  cmt = "IPI_C1",  ii = 42, addl = 8))
  } else if (scenario == "TTFields + Cisplatin/Pemetrexed (STELLAR)") {
    seq(ev(amt = 135, cmt = "CIS_C1", ii = 21, addl = 5),
        ev(amt = 900, cmt = "PEM_C1", ii = 21, addl = 5))
  } else if (scenario == "Anetumab ravtansine (mesothelin-high)") {
    ev(amt = 455, cmt = "ADC_C1", ii = 21, addl = 11)
  } else if (scenario == "Rucaparib (BAP1-deficient, MiST1)") {
    ev(amt = 600, cmt = "GUT_RUCA", ii = 0.5, addl = 2 * horizon)
  } else if (scenario == "Surgery (P/D) + chemo + talc pleurodesis") {
    seq(ev(amt = 135, cmt = "CIS_C1", ii = 21, addl = 5),
        ev(amt = 900, cmt = "PEM_C1", ii = 21, addl = 5))
  } else {
    ev(amt = 0, cmt = "CIS_C1")
  }
}

run_sim <- function(scenario, horizon, params) {
  mod <- get_model()
  ev_set <- build_events(scenario, horizon)
  par <- list(
    HIST_SARC_FRAC = params$sarc_frac,
    BAP1_DEFICIENT = ifelse(params$bap1_def, 1, 0),
    MSLN_HIGH      = ifelse(params$msln_high, 1, 0),
    PDL1_HIGH      = ifelse(params$pdl1_high, 1, 0),
    ECOG0          = params$ecog0,
    TTF_ACTIVE     = ifelse(scenario == "TTFields + Cisplatin/Pemetrexed (STELLAR)", 1, 0)
  )
  out <- mod %>% param(par) %>% mrgsim(events = ev_set, end = horizon, delta = 1) %>%
    as_tibble() %>% mutate(scenario = scenario)
  if (scenario == "Surgery (P/D) + chemo + talc pleurodesis") {
    out <- out %>%
      mutate(Tumor_epithelioid = ifelse(time >= 1, Tumor_epithelioid * 0.15, Tumor_epithelioid),
             Tumor_sarcomatoid = ifelse(time >= 1, Tumor_sarcomatoid * 0.15, Tumor_sarcomatoid),
             Tumor_total = Tumor_epithelioid + Tumor_sarcomatoid,
             Effusion_vol = ifelse(time >= 1, Effusion_vol * 0.45, Effusion_vol))
  }
  out
}

SCENARIO_LIST <- c(
  "Untreated / BSC", "Cisplatin+Pemetrexed (EMPHACIS)",
  "Cisplatin+Pemetrexed+Bevacizumab (MAPS)",
  "Nivolumab+Ipilimumab 1L (CheckMate 743)",
  "Nivolumab monotherapy 2L (CONFIRM)",
  "Nivolumab+Ipilimumab 2L (MAPS2)",
  "TTFields + Cisplatin/Pemetrexed (STELLAR)",
  "Anetumab ravtansine (mesothelin-high)",
  "Rucaparib (BAP1-deficient, MiST1)",
  "Surgery (P/D) + chemo + talc pleurodesis"
)

## ---------- UI ----------
ui <- dashboardPage(
  dashboardHeader(title = "MPM QSP"),
  dashboardSidebar(
    sidebarMenu(
      menuItem("1. Patient profile",     tabName = "patient", icon = icon("user")),
      menuItem("2. Drug PK",             tabName = "pk",      icon = icon("syringe")),
      menuItem("3. Molecular pathway",   tabName = "pathway", icon = icon("dna")),
      menuItem("4. Immune/TME & effusion", tabName = "tme",   icon = icon("shield-virus")),
      menuItem("5. Clinical endpoints",  tabName = "clin",    icon = icon("chart-line")),
      menuItem("6. Scenario comparison", tabName = "compare", icon = icon("balance-scale")),
      menuItem("7. Safety",              tabName = "safety",  icon = icon("triangle-exclamation")),
      menuItem("8. Biomarkers",          tabName = "bio",     icon = icon("flask"))
    ),
    hr(),
    selectInput("scenario", "Scenario:", SCENARIO_LIST,
                selected = "Nivolumab+Ipilimumab 1L (CheckMate 743)"),
    sliderInput("horizon",   "Simulation horizon (days):", 28, 730, 252, step = 7),
    sliderInput("sarc_frac", "Sarcomatoid fraction (0=pure epithelioid):", 0, 1, 0.20, step = 0.05),
    checkboxInput("bap1_def", "BAP1/BRCA1-deficient tumor", FALSE),
    checkboxInput("msln_high","Mesothelin-high tumor", TRUE),
    checkboxInput("pdl1_high","PD-L1-high tumor", FALSE),
    sliderInput("ecog0",     "Baseline ECOG PS:", 0, 3, 1),
    actionButton("run", "Run simulation", icon = icon("play"),
                 style = "color:#fff;background:#1f6feb"),
    br(), br(),
    actionButton("run_all", "Run all 10 scenarios", icon = icon("rocket"),
                 style = "color:#fff;background:#0f5132")
  ),
  dashboardBody(
    tabItems(
      tabItem("patient",
        fluidRow(
          box(width = 12, title = "Patient / tumor profile summary", status = "primary",
              solidHeader = TRUE, DTOutput("patient_table"))
        ),
        fluidRow(
          box(width = 12, title = "About this dashboard", status = "info",
              p("This dashboard runs the mrgsolve QSP model for malignant",
                "pleural mesothelioma (MPM). Pick a treatment scenario and",
                "tumor-biology profile in the left panel, then press ",
                strong("Run simulation"), " to update plots."),
              p("Each tab focuses on a layer of the model: drug PK, the",
                "BAP1/CDKN2A/NF2-Hippo-YAP molecular pathway driving the",
                "two-clone (epithelioid/sarcomatoid) tumor-burden network,",
                "immune microenvironment & pleural effusion, clinical",
                "endpoints (tumor volume/OS), scenario comparison, safety",
                "(myelosuppression/renal/irAE), and mesothelin biomarker.")
          )
        )
      ),

      tabItem("pk",
        fluidRow(
          box(width = 12, title = "Plasma drug exposure (selected scenario)",
              status = "primary", solidHeader = TRUE, plotOutput("pk_plot", height = 480))
        )
      ),

      tabItem("pathway",
        fluidRow(
          box(width = 6, title = "Tumor burden by histologic clone",
              plotOutput("clone_plot", 360)),
          box(width = 6, title = "Total tumor burden vs. carrying capacity",
              plotOutput("tumor_total_plot", 360))
        ),
        fluidRow(
          box(width = 12, title = "Pathway note", status = "info",
              p("Epithelioid clones grow more slowly but retain higher",
                "chemo-, IO-, and ADC-sensitivity; sarcomatoid clones grow",
                "faster and are relatively drug-resistant (EMT/Hippo-YAP",
                "activation), reflecting the histologic-subtype prognosis",
                "gradient seen clinically."))
        )
      ),

      tabItem("tme",
        fluidRow(
          box(width = 6, title = "Pleural effusion volume", plotOutput("eff_plot", 360)),
          box(width = 6, title = "ECOG performance status", plotOutput("ecog_plot", 360))
        )
      ),

      tabItem("clin",
        fluidRow(
          box(width = 6, title = "Tumor volume (relative units)",
              plotOutput("tumvol_plot", 360)),
          box(width = 6, title = "Survival probability (cumulative-hazard model)",
              plotOutput("surv_plot", 360))
        ),
        fluidRow(
          box(width = 12, title = "Endpoint summary table (Day-0, mid, end)",
              status = "info", DTOutput("endpoint_table"))
        )
      ),

      tabItem("compare",
        fluidRow(
          box(width = 12, title = "Scenario comparison panel",
              status = "warning", solidHeader = TRUE,
              p("Runs all ten built-in scenarios with the current tumor-",
                "biology profile; press the button in the left panel."),
              plotOutput("compare_plot", height = 600)
          )
        )
      ),

      tabItem("safety",
        fluidRow(
          box(width = 6, title = "Circulating ANC (myelosuppression)",
              plotOutput("anc_plot", 360)),
          box(width = 6, title = "Creatinine clearance (nephrotoxicity)",
              plotOutput("crcl_plot", 360))
        )
      ),

      tabItem("bio",
        fluidRow(
          box(width = 12, title = "Serum mesothelin (SMRP) biomarker trajectory",
              status = "primary", solidHeader = TRUE, plotOutput("smrp_plot", 380))
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
    p <- list(sarc_frac = input$sarc_frac, bap1_def = input$bap1_def,
              msln_high = input$msln_high, pdl1_high = input$pdl1_high,
              ecog0 = input$ecog0)
    results(run_sim(input$scenario, input$horizon, p))
  }, ignoreNULL = FALSE)

  observeEvent(input$run_all, {
    showNotification("Running 10 scenarios…", type = "message", duration = 1)
    p <- list(sarc_frac = input$sarc_frac, bap1_def = input$bap1_def,
              msln_high = input$msln_high, pdl1_high = input$pdl1_high,
              ecog0 = input$ecog0)
    out <- lapply(SCENARIO_LIST, function(sc) run_sim(sc, input$horizon, p))
    all_results(bind_rows(out))
  })

  # --- Patient table ---
  output$patient_table <- renderDT({
    tibble(
      Field = c("Sarcomatoid fraction", "BAP1/BRCA1-deficient", "Mesothelin-high",
                "PD-L1-high", "Baseline ECOG", "Scenario", "Horizon (d)"),
      Value = c(input$sarc_frac, input$bap1_def, input$msln_high,
                input$pdl1_high, input$ecog0, input$scenario, input$horizon)
    ) %>% datatable(rownames = FALSE, options = list(dom = "t"))
  })

  # --- PK ---
  output$pk_plot <- renderPlot({
    df <- results(); req(df)
    df %>% select(time, conc_cisplatin, conc_pemetrexed, conc_bevacizumab,
                  conc_nivolumab, conc_ipilimumab, conc_anetumab, conc_rucaparib) %>%
      pivot_longer(-time) %>%
      filter(value > 0) %>%
      ggplot(aes(time, value, colour = name)) +
        geom_line(linewidth = 0.7) +
        scale_y_continuous(trans = "log1p") +
        labs(x = "Time (days)", y = "Plasma conc (mg/L)", colour = "Drug") +
        theme_minimal(base_size = 13)
  })

  # --- Pathway / tumor clones ---
  output$clone_plot <- renderPlot({
    df <- results(); req(df)
    df %>% select(time, Tumor_epithelioid, Tumor_sarcomatoid) %>%
      pivot_longer(-time) %>%
      ggplot(aes(time, value, colour = name)) +
        geom_line(linewidth = 0.9) +
        labs(x = "Day", y = "Tumor burden (a.u.)", colour = NULL) +
        theme_minimal(base_size = 12)
  })
  output$tumor_total_plot <- renderPlot({
    df <- results(); req(df)
    ggplot(df, aes(time, Tumor_total)) +
      geom_line(colour = "#9b1c1c", linewidth = 1) +
      labs(x = "Day", y = "Total tumor burden (a.u.)") +
      theme_minimal(base_size = 13)
  })

  # --- TME / effusion ---
  output$eff_plot <- renderPlot({
    df <- results(); req(df)
    ggplot(df, aes(time, Effusion_vol)) +
      geom_line(colour = "#1f6feb", linewidth = 1) +
      labs(x = "Day", y = "Pleural effusion volume (a.u.)") +
      theme_minimal(base_size = 13)
  })
  output$ecog_plot <- renderPlot({
    df <- results(); req(df)
    ggplot(df, aes(time, ECOG_score)) +
      geom_line(colour = "#6d4c41", linewidth = 1) +
      scale_y_continuous(limits = c(0, 4)) +
      labs(x = "Day", y = "ECOG performance status") +
      theme_minimal(base_size = 13)
  })

  # --- Clinical endpoints ---
  output$tumvol_plot <- renderPlot({
    df <- results(); req(df)
    ggplot(df, aes(time, Tumor_total)) +
      geom_line(colour = "#264653", linewidth = 1) +
      labs(x = "Day", y = "Tumor volume (a.u.)") +
      theme_minimal(base_size = 13)
  })
  output$surv_plot <- renderPlot({
    df <- results(); req(df)
    ggplot(df, aes(time, Survival_prob)) +
      geom_line(colour = "#0f5132", linewidth = 1) +
      scale_y_continuous(limits = c(0, 1)) +
      labs(x = "Day", y = "Modeled survival probability") +
      theme_minimal(base_size = 13)
  })
  output$endpoint_table <- renderDT({
    df <- results(); req(df)
    tp <- c(0, round(max(df$time) / 2), max(df$time))
    df %>% filter(time %in% tp) %>%
      select(Day = time, Tumor_total, Effusion_vol, SMRP_biomarker,
             ANC_circulating, CrCl_mLmin, ECOG_score, Survival_prob) %>%
      mutate(across(-Day, ~ round(.x, 3))) %>%
      datatable(rownames = FALSE,
                options = list(dom = "t",
                               columnDefs = list(list(className = "dt-center", targets = "_all"))))
  })

  # --- Compare ---
  output$compare_plot <- renderPlot({
    df <- all_results(); req(df)
    df %>% select(time, scenario, Tumor_total, Survival_prob, Effusion_vol, SMRP_biomarker) %>%
      pivot_longer(-c(time, scenario)) %>%
      ggplot(aes(time, value, colour = scenario)) +
        geom_line(linewidth = 0.6, alpha = 0.85) +
        facet_wrap(~name, scales = "free_y") +
        labs(x = "Day", y = NULL, colour = NULL) +
        theme_minimal(base_size = 12) +
        theme(legend.position = "bottom")
  })

  # --- Safety ---
  output$anc_plot <- renderPlot({
    df <- results(); req(df)
    ggplot(df, aes(time, ANC_circulating)) +
      geom_line(colour = "#c0392b", linewidth = 1) +
      geom_hline(yintercept = 1.5, lty = 2, colour = "grey50") +
      labs(x = "Day", y = "ANC (10^9/L)",
           caption = "Dashed line = grade-3 neutropenia threshold (1.5x10^9/L)") +
      theme_minimal(base_size = 13)
  })
  output$crcl_plot <- renderPlot({
    df <- results(); req(df)
    ggplot(df, aes(time, CrCl_mLmin)) +
      geom_line(colour = "#8a6d3b", linewidth = 1) +
      labs(x = "Day", y = "Creatinine clearance (mL/min)") +
      theme_minimal(base_size = 13)
  })

  # --- Biomarkers ---
  output$smrp_plot <- renderPlot({
    df <- results(); req(df)
    ggplot(df, aes(time, SMRP_biomarker)) +
      geom_line(colour = "#5c6370", linewidth = 1) +
      labs(x = "Day", y = "SMRP biomarker (a.u.)") +
      theme_minimal(base_size = 13)
  })
}

shinyApp(ui, server)
