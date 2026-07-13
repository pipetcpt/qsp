## ============================================================================
## Dry Eye Disease (DED) QSP — Shiny Dashboard
## ----------------------------------------------------------------------------
## 8 tabs: Patient · Drug exposure · Tear-film/PD · Inflammatory cascade ·
##         Clinical endpoints · Neurosensory · Scenario comparison ·
##         Biomarkers/Safety
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
MODEL_PATH <- "ded_mrgsolve_model.R"
get_model <- function() {
  if (!exists(".DED_MOD", envir = .GlobalEnv) || is.null(.GlobalEnv$.DED_MOD)) {
    assign(".DED_MOD", mread_cache("ded", project = "."), envir = .GlobalEnv)
  }
  .GlobalEnv$.DED_MOD
}

## ---------- Scenario builders ----------
build_events <- function(scenario, horizon) {
  if (scenario == "Untreated") {
    ev(amt = 0, cmt = "CSA_TF")
  } else if (scenario == "Cyclosporine 0.05% BID (Restasis)") {
    ev(amt = 1.0, cmt = "CSA_TF", ii = 0.5, addl = 2 * horizon)
  } else if (scenario == "Cyclosporine 0.09% BID (Cequa)") {
    ev(amt = 1.6, cmt = "CSA_TF", ii = 0.5, addl = 2 * horizon)
  } else if (scenario == "Lifitegrast 5% BID (Xiidra)") {
    ev(amt = 1.2, cmt = "LIF_TF", ii = 0.5, addl = 2 * horizon)
  } else if (scenario == "Loteprednol 0.25% QID x 2wk induction") {
    ev(amt = 1.0, cmt = "LOT_TF", ii = 0.25, addl = 4 * 14 - 1)
  } else if (scenario == "Perfluorohexyloctane QID (Miebo)") {
    ev(amt = 1.3, cmt = "PFHO_TF", ii = 0.25, addl = 4 * horizon)
  } else if (scenario == "Varenicline nasal spray BID (Tyrvaya)") {
    ev(amt = 1.0, cmt = "VAR_NAS", ii = 0.5, addl = 2 * horizon)
  } else if (scenario == "Diquafosol 3% 6x/day") {
    ev(amt = 0.8, cmt = "DQS_TF", ii = 1 / 6, addl = 6 * horizon)
  } else if (scenario == "Combo: Loteprednol induction + Cyclosporine maintenance") {
    seq(ev(amt = 1.0, cmt = "LOT_TF", ii = 0.25, addl = 4 * 14 - 1),
        ev(amt = 1.0, cmt = "CSA_TF", ii = 0.5, addl = 2 * horizon, time = 0))
  } else {
    ev(amt = 0, cmt = "CSA_TF")
  }
}

run_sim <- function(scenario, horizon, params) {
  mod <- get_model()
  ev_set <- build_events(scenario, horizon)
  par <- list(
    MGD_SEV         = params$mgd_sev,
    ADDE_SEV        = params$adde_sev,
    AUTOIMMUNE      = ifelse(params$autoimmune, 1, 0),
    REFRACTIVE_SURG = ifelse(params$refractive, 1, 0),
    BASE_OSDI       = params$base_osdi
  )
  mod %>% param(par) %>% mrgsim(events = ev_set, end = horizon, delta = 1) %>%
    as_tibble() %>% mutate(scenario = scenario)
}

SCENARIO_LIST <- c(
  "Untreated", "Cyclosporine 0.05% BID (Restasis)", "Cyclosporine 0.09% BID (Cequa)",
  "Lifitegrast 5% BID (Xiidra)", "Loteprednol 0.25% QID x 2wk induction",
  "Perfluorohexyloctane QID (Miebo)", "Varenicline nasal spray BID (Tyrvaya)",
  "Diquafosol 3% 6x/day", "Combo: Loteprednol induction + Cyclosporine maintenance"
)

## ---------- UI ----------
ui <- dashboardPage(
  dashboardHeader(title = "Dry Eye Disease QSP"),
  dashboardSidebar(
    sidebarMenu(
      menuItem("1. Patient profile",       tabName = "patient", icon = icon("user")),
      menuItem("2. Drug exposure",         tabName = "pk",      icon = icon("eye-dropper")),
      menuItem("3. Tear-film / PD",        tabName = "tearfilm",icon = icon("water")),
      menuItem("4. Inflammatory cascade",  tabName = "inflam",  icon = icon("fire-flame-curved")),
      menuItem("5. Clinical endpoints",    tabName = "clin",    icon = icon("notes-medical")),
      menuItem("6. Neurosensory",          tabName = "neuro",   icon = icon("bolt")),
      menuItem("7. Scenario comparison",   tabName = "compare", icon = icon("chart-line")),
      menuItem("8. Biomarkers / Safety",   tabName = "safety",  icon = icon("shield-alt"))
    ),
    hr(),
    selectInput("scenario", "Scenario:", SCENARIO_LIST,
                selected = "Cyclosporine 0.05% BID (Restasis)"),
    sliderInput("horizon", "Simulation horizon (days):", 14, 180, 84, step = 7),
    sliderInput("mgd_sev",  "Evaporative MGD severity (0-1):", 0, 1, 0.6, step = 0.05),
    sliderInput("adde_sev", "Aqueous-deficient (lacrimal) severity (0-1):", 0, 1, 0.3, step = 0.05),
    checkboxInput("autoimmune", "Sjögren/autoimmune-associated DED", FALSE),
    checkboxInput("refractive", "Prior refractive surgery (LASIK/PRK)", FALSE),
    sliderInput("base_osdi", "Baseline OSDI (0-100):", 10, 90, 45),
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
              p(strong("Phenotype legend:"),
                "Evaporative DED (MGD-driven lipid-layer deficiency) and aqueous-",
                "deficient DED (lacrimal-gland/Sjögren-driven) commonly coexist; ",
                "severity sliders combine additively up to a ceiling.")
          )
        ),
        fluidRow(
          box(width = 12, title = "About this dashboard", status = "info",
              p("This dashboard runs the mrgsolve QSP model for dry eye disease",
                "(DED), linking tear-film hyperosmolarity, meibomian/lacrimal-gland",
                "dysfunction, innate/adaptive ocular-surface inflammation, epithelial",
                "barrier breakdown, and corneal neurosensory dysfunction to clinical",
                "endpoints (OSDI, TBUT, corneal staining, Schirmer)."),
              p("Pick a treatment scenario in the left panel, adjust patient severity,",
                "then press ", strong("Run simulation"), " to update plots.")
          )
        )
      ),

      tabItem("pk",
        fluidRow(
          box(width = 12, title = "Local ocular-surface drug exposure (selected scenario)",
              status = "primary", solidHeader = TRUE, plotOutput("pk_plot", height = 480))
        )
      ),

      tabItem("tearfilm",
        fluidRow(
          box(width = 6, title = "Tear osmolarity & lipid-layer quality", plotOutput("osm_lipid_plot", 360)),
          box(width = 6, title = "Aqueous production & mucin/goblet density", plotOutput("aqp_mucin_plot", 360))
        )
      ),

      tabItem("inflam",
        fluidRow(
          box(width = 6, title = "Dendritic-cell activation & IL-17A tone", plotOutput("dc_il17_plot", 360)),
          box(width = 6, title = "MMP-9 biomarker & epithelial barrier integrity", plotOutput("mmp9_barrier_plot", 360))
        )
      ),

      tabItem("clin",
        fluidRow(
          box(width = 6, title = "TBUT & Schirmer test", plotOutput("tbut_schirmer_plot", 360)),
          box(width = 6, title = "Corneal staining & OSDI symptom score", plotOutput("stain_osdi_plot", 360))
        )
      ),

      tabItem("neuro",
        fluidRow(
          box(width = 6, title = "Corneal sub-basal nerve density", plotOutput("nerve_plot", 360)),
          box(width = 6, title = "Neuropathic ocular pain score", plotOutput("pain_plot", 360))
        )
      ),

      tabItem("compare",
        fluidRow(
          box(width = 12, title = "Scenario comparison panel",
              status = "warning", solidHeader = TRUE,
              p("Runs all nine built-in scenarios with the current patient profile;",
                "press the button below."),
              actionButton("run_all", "Run all scenarios", icon = icon("rocket"),
                           style = "color:#fff;background:#0f5132"),
              br(), br(),
              plotOutput("compare_plot", height = 640)
          )
        )
      ),

      tabItem("safety",
        fluidRow(
          box(width = 6, title = "Loteprednol cumulative exposure -> IOP-risk score", plotOutput("iop_plot", 360)),
          box(width = 6, title = "Evaporation rate (barrier restoration)", plotOutput("evap_plot", 360))
        ),
        fluidRow(
          box(width = 12, title = "Endpoint summary table (Day-0, Day-28, Day-end)",
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
      mgd_sev = input$mgd_sev, adde_sev = input$adde_sev,
      autoimmune = input$autoimmune, refractive = input$refractive,
      base_osdi = input$base_osdi
    )
    results(run_sim(input$scenario, input$horizon, p))
  }, ignoreNULL = FALSE)

  observeEvent(input$run_all, {
    showNotification("Running 9 scenarios…", type = "message", duration = 1)
    p <- list(
      mgd_sev = input$mgd_sev, adde_sev = input$adde_sev,
      autoimmune = input$autoimmune, refractive = input$refractive,
      base_osdi = input$base_osdi
    )
    out <- lapply(SCENARIO_LIST, function(sc) run_sim(sc, input$horizon, p))
    all_results(bind_rows(out))
  })

  # --- Patient table ---
  output$patient_table <- renderDT({
    tibble(
      Field = c("MGD severity", "ADDE severity", "Autoimmune/Sjögren",
                "Prior refractive surgery", "Baseline OSDI", "Scenario", "Horizon (d)"),
      Value = c(input$mgd_sev, input$adde_sev, input$autoimmune,
                input$refractive, input$base_osdi, input$scenario, input$horizon)
    ) %>% datatable(rownames = FALSE, options = list(dom = "t"))
  })

  # --- Drug exposure ---
  output$pk_plot <- renderPlot({
    df <- results(); req(df)
    df %>% select(time, CsA_effect, Lifitegrast_effect, Loteprednol_effect,
                  PFHO_effect, Varenicline_effect, Diquafosol_effect) %>%
      pivot_longer(-time) %>%
      ggplot(aes(time, value, colour = name)) +
        geom_line(linewidth = 0.8) +
        labs(x = "Time (days)", y = "Fractional drug effect (0-1)", colour = "Drug effect") +
        theme_minimal(base_size = 13)
  })

  # --- Tear film ---
  output$osm_lipid_plot <- renderPlot({
    df <- results(); req(df)
    df %>% select(time, Tear_osmolarity, Lipid_layer) %>%
      pivot_longer(-time) %>%
      ggplot(aes(time, value, colour = name)) +
        geom_line(linewidth = 0.8) +
        facet_wrap(~name, scales = "free_y") +
        theme_minimal(base_size = 12) + theme(legend.position = "none")
  })
  output$aqp_mucin_plot <- renderPlot({
    df <- results(); req(df)
    df %>% select(time, Aqueous_drive, Mucin_density) %>%
      pivot_longer(-time) %>%
      ggplot(aes(time, value, colour = name)) +
        geom_line(linewidth = 0.8) +
        facet_wrap(~name, scales = "free_y") +
        theme_minimal(base_size = 12) + theme(legend.position = "none")
  })

  # --- Inflammatory cascade ---
  output$dc_il17_plot <- renderPlot({
    df <- results(); req(df)
    df %>% select(time, DC_activation, IL17_tone) %>%
      pivot_longer(-time) %>%
      ggplot(aes(time, value, colour = name)) +
        geom_line(linewidth = 0.8) +
        facet_wrap(~name, scales = "free_y") +
        theme_minimal(base_size = 12) + theme(legend.position = "none")
  })
  output$mmp9_barrier_plot <- renderPlot({
    df <- results(); req(df)
    df %>% select(time, MMP9_biomarker, Barrier_integrity) %>%
      pivot_longer(-time) %>%
      ggplot(aes(time, value, colour = name)) +
        geom_line(linewidth = 0.8) +
        facet_wrap(~name, scales = "free_y") +
        theme_minimal(base_size = 12) + theme(legend.position = "none")
  })

  # --- Clinical endpoints ---
  output$tbut_schirmer_plot <- renderPlot({
    df <- results(); req(df)
    df %>% select(time, TBUT_s, Schirmer_mm) %>%
      pivot_longer(-time) %>%
      ggplot(aes(time, value, colour = name)) +
        geom_line(linewidth = 0.8) +
        facet_wrap(~name, scales = "free_y") +
        theme_minimal(base_size = 12) + theme(legend.position = "none")
  })
  output$stain_osdi_plot <- renderPlot({
    df <- results(); req(df)
    df %>% select(time, Corneal_staining, OSDI_score) %>%
      pivot_longer(-time) %>%
      ggplot(aes(time, value, colour = name)) +
        geom_line(linewidth = 0.8) +
        facet_wrap(~name, scales = "free_y") +
        theme_minimal(base_size = 12) + theme(legend.position = "none")
  })

  # --- Neurosensory ---
  output$nerve_plot <- renderPlot({
    df <- results(); req(df)
    ggplot(df, aes(time, Nerve_density)) +
      geom_line(colour = "#2a9d8f", linewidth = 1) +
      labs(x = "Day", y = "Corneal nerve density (fibers/mm^2-equiv.)") +
      theme_minimal(base_size = 13)
  })
  output$pain_plot <- renderPlot({
    df <- results(); req(df)
    ggplot(df, aes(time, Neuropathic_pain)) +
      geom_line(colour = "#9b1c1c", linewidth = 1) +
      scale_y_continuous(limits = c(0, 10)) +
      labs(x = "Day", y = "Neuropathic ocular pain (0-10)") +
      theme_minimal(base_size = 13)
  })

  # --- Compare ---
  output$compare_plot <- renderPlot({
    df <- all_results(); req(df)
    df %>% select(time, scenario, OSDI_score, TBUT_s, Corneal_staining, Schirmer_mm) %>%
      pivot_longer(-c(time, scenario)) %>%
      ggplot(aes(time, value, colour = scenario)) +
        geom_line(linewidth = 0.6, alpha = 0.85) +
        facet_wrap(~name, scales = "free_y") +
        labs(x = "Day", y = NULL, colour = NULL) +
        theme_minimal(base_size = 12) +
        theme(legend.position = "bottom")
  })

  # --- Safety ---
  output$iop_plot <- renderPlot({
    df <- results(); req(df)
    ggplot(df, aes(time, IOP_risk_score)) +
      geom_line(colour = "#8a6d3b", linewidth = 1) +
      labs(x = "Day", y = "Cumulative IOP-risk score (loteprednol exposure)") +
      theme_minimal(base_size = 13)
  })
  output$evap_plot <- renderPlot({
    df <- results(); req(df)
    ggplot(df, aes(time, Evap_rate)) +
      geom_line(colour = "#264653", linewidth = 1) +
      labs(x = "Day", y = "Relative evaporation rate (a.u.)") +
      theme_minimal(base_size = 13)
  })

  output$endpoint_table <- renderDT({
    df <- results(); req(df)
    tp <- c(0, 28, max(df$time))
    df %>% filter(time %in% tp) %>%
      select(Day = time, OSDI_score, TBUT_s, Schirmer_mm, Corneal_staining,
             Tear_osmolarity, MMP9_biomarker) %>%
      mutate(across(-Day, ~ round(.x, 2))) %>%
      datatable(rownames = FALSE,
                options = list(dom = "t",
                               columnDefs = list(list(className='dt-center', targets="_all"))))
  })
}

shinyApp(ui, server)
