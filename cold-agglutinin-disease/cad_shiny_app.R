## ============================================================================
## Cold Agglutinin Disease (CAD) QSP — Shiny Dashboard
## ----------------------------------------------------------------------------
## 8 tabs: Patient profile · Drug PK · Complement pathway PD · Hemolysis
##         biomarkers · Clinical endpoints · Scenario comparison ·
##         Cold-exposure challenge · References
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
  if (!exists(".CAD_MOD", envir = .GlobalEnv) || is.null(.GlobalEnv$.CAD_MOD)) {
    assign(".CAD_MOD", mread_cache("cad", project = "."), envir = .GlobalEnv)
  }
  .GlobalEnv$.CAD_MOD
}

## ---------- Scenario builders ----------
build_events <- function(scenario, horizon) {
  if (scenario == "Untreated natural history") {
    ev(amt = 0, cmt = "RTX_CENT")
  } else if (scenario == "Rituximab monotherapy (weekly x4)") {
    ev(amt = 1.0, cmt = "RTX_CENT", ii = 7, addl = 3)
  } else if (scenario == "Rituximab + Bendamustine (RB, q28d x4)") {
    seq(ev(amt = 1.0, cmt = "RTX_CENT", ii = 7, addl = 3),
        ev(amt = 1.0, cmt = "BEN_CENT", ii = 28, addl = 3, time = 0))
  } else if (scenario == "Sutimlimab (load d0/d7 + biweekly)") {
    ev(amt = 1.3, cmt = "SUT_CENT", time = 0, ii = 7, addl = 1) +
      ev(amt = 1.1, cmt = "SUT_CENT", time = 14, ii = 14, addl = max(floor((horizon - 14) / 14), 0))
  } else if (scenario == "Pegcetacoplan SC (biweekly, investigational)") {
    ev(amt = 1.2, cmt = "PEG_DEPOT", ii = 14, addl = max(floor(horizon / 14), 0))
  } else if (scenario == "Sutimlimab discontinued at day 90 (relapse)") {
    ev(amt = 1.3, cmt = "SUT_CENT", time = 0, ii = 7, addl = 1) +
      ev(amt = 1.1, cmt = "SUT_CENT", time = 14, ii = 14, addl = max(floor((min(90, horizon) - 14) / 14), 0))
  } else if (scenario == "Acute cold-exposure crisis (on sutimlimab)") {
    ev(amt = 1.3, cmt = "SUT_CENT", time = 0, ii = 7, addl = 1) +
      ev(amt = 1.1, cmt = "SUT_CENT", time = 14, ii = 14, addl = max(floor((horizon - 14) / 14), 0)) +
      ev(amt = 3.0, cmt = "COLD_PULSE", time = min(60, horizon / 2))
  } else {
    ev(amt = 0, cmt = "RTX_CENT")
  }
}

run_sim <- function(scenario, horizon, params) {
  mod <- get_model()
  ev_set <- build_events(scenario, horizon)
  par <- list(
    SEV            = params$sev,
    SECONDARY      = ifelse(params$secondary, 1, 0),
    AMB_COLD       = params$amb_cold,
    COLD_AVOIDANCE = ifelse(scenario == "Untreated natural history" && params$cold_avoidance, 1,
                             ifelse(params$cold_avoidance, 1, 0)),
    BASE_HB        = params$base_hb
  )
  mod %>% param(par) %>% mrgsim(events = ev_set, end = horizon, delta = 0.5) %>%
    as_tibble() %>% mutate(scenario = scenario)
}

SCENARIO_LIST <- c(
  "Untreated natural history", "Rituximab monotherapy (weekly x4)",
  "Rituximab + Bendamustine (RB, q28d x4)", "Sutimlimab (load d0/d7 + biweekly)",
  "Pegcetacoplan SC (biweekly, investigational)",
  "Sutimlimab discontinued at day 90 (relapse)",
  "Acute cold-exposure crisis (on sutimlimab)"
)

## ---------- UI ----------
ui <- dashboardPage(
  dashboardHeader(title = "Cold Agglutinin Disease QSP"),
  dashboardSidebar(
    sidebarMenu(
      menuItem("1. Patient profile",       tabName = "patient", icon = icon("user")),
      menuItem("2. Drug PK",               tabName = "pk",      icon = icon("pills")),
      menuItem("3. Complement pathway PD", tabName = "pd",      icon = icon("dna")),
      menuItem("4. Hemolysis biomarkers",  tabName = "bio",     icon = icon("vial")),
      menuItem("5. Clinical endpoints",    tabName = "clin",    icon = icon("notes-medical")),
      menuItem("6. Scenario comparison",   tabName = "compare", icon = icon("chart-line")),
      menuItem("7. Cold-exposure challenge", tabName = "cold",  icon = icon("temperature-low")),
      menuItem("8. References",            tabName = "refs",   icon = icon("book"))
    ),
    hr(),
    selectInput("scenario", "Scenario:", SCENARIO_LIST,
                selected = "Sutimlimab (load d0/d7 + biweekly)"),
    sliderInput("horizon", "Simulation horizon (days):", 14, 365, 180, step = 1),
    sliderInput("sev", "Clonal disease severity (0-1):", 0, 1, 0.65, step = 0.05),
    checkboxInput("secondary", "Secondary CAD (infection-triggered clone)", FALSE),
    sliderInput("base_hb", "Presenting hemoglobin (g/dL):", 5, 13, 9, step = 0.5),
    sliderInput("amb_cold", "Ambient cold-exposure level (0-1):", 0, 1, 0.5, step = 0.05),
    checkboxInput("cold_avoidance", "Behavioral cold avoidance counseling", FALSE),
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
              p(strong("Note:"), "Cold Agglutinin Disease (CAD) is characteristically",
                strong(" steroid-refractory"), "— unlike warm autoimmune hemolytic anemia,",
                "corticosteroids are not an effective treatment arm and are intentionally",
                "omitted from this model's scenarios.")
          )
        ),
        fluidRow(
          box(width = 12, title = "About this dashboard", status = "info",
              p("This dashboard runs the mrgsolve QSP model for Cold Agglutinin Disease,",
                "linking clonal marrow B-cell/lymphoplasmacytic IgM cold-agglutinin",
                "production, cold-dependent RBC binding, classical complement pathway",
                "activation (C1s-C4-C3), C3b-mediated extravascular (Kupffer-cell)",
                "hemolysis, C5b-9 (MAC)-mediated intravascular hemolysis, and downstream",
                "hemolysis biomarkers."),
              p("Pick a treatment scenario in the left panel, adjust patient severity,",
                "then press ", strong("Run simulation"), " to update plots.")
          )
        )
      ),

      tabItem("pk",
        fluidRow(
          box(width = 12, title = "Drug PD-effect trajectories (selected scenario)",
              status = "primary", solidHeader = TRUE, plotOutput("pk_plot", height = 480))
        )
      ),

      tabItem("pd",
        fluidRow(
          box(width = 6, title = "Classical pathway: C1s activity, C4, C3", plotOutput("cp_plot", height = 360)),
          box(width = 6, title = "C3b opsonization & C5b-9 MAC", plotOutput("mac_plot", height = 360))
        ),
        fluidRow(
          box(width = 12, title = "Clonal burden & IgM cold-agglutinin titer", plotOutput("clone_igm_plot", height = 360))
        )
      ),

      tabItem("bio",
        fluidRow(
          box(width = 6, title = "Hemoglobin & reticulocytes", plotOutput("hb_retic_plot", height = 360)),
          box(width = 6, title = "Indirect bilirubin & LDH", plotOutput("bili_ldh_plot", height = 360))
        ),
        fluidRow(
          box(width = 12, title = "Haptoglobin (consumed by intravascular hemolysis)",
              plotOutput("hapto_plot", height = 320))
        )
      ),

      tabItem("clin",
        fluidRow(
          box(width = 12, title = "Endpoint summary table (Day-0, Day-30, Day-90, Day-end)",
              status = "info", DTOutput("endpoint_table"))
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

      tabItem("cold",
        fluidRow(
          box(width = 12, title = "Acute cold-exposure hemolytic crisis challenge",
              status = "danger", solidHeader = TRUE,
              p("Select the ", strong("Acute cold-exposure crisis (on sutimlimab)"),
                "scenario in the left panel to superimpose a decaying cold-exposure",
                "pulse at mid-horizon on a sutimlimab background, and observe the",
                "transient MAC/hemolysis spike and Hb dip."),
              plotOutput("cold_plot", height = 420)
          )
        )
      ),

      tabItem("refs",
        fluidRow(
          box(width = 12, title = "Key references", status = "info",
              p("See ", code("cad_references.md"), " in this directory for the full",
                "30+ item bibliography (mechanism, clonal biology, sutimlimab/",
                "pegcetacoplan trials, rituximab-bendamustine, cold-exposure biology)."),
              tags$a(href = "cad_references.md", "Open cad_references.md", target = "_blank")
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
      sev = input$sev, secondary = input$secondary, base_hb = input$base_hb,
      amb_cold = input$amb_cold, cold_avoidance = input$cold_avoidance
    )
    results(run_sim(input$scenario, input$horizon, p))
  }, ignoreNULL = FALSE)

  observeEvent(input$run_all, {
    showNotification("Running 7 scenarios…", type = "message", duration = 1)
    p <- list(
      sev = input$sev, secondary = input$secondary, base_hb = input$base_hb,
      amb_cold = input$amb_cold, cold_avoidance = input$cold_avoidance
    )
    out <- lapply(SCENARIO_LIST, function(sc) run_sim(sc, input$horizon, p))
    all_results(bind_rows(out))
  })

  # --- Patient table ---
  output$patient_table <- renderDT({
    tibble(
      Field = c("Clonal disease severity", "Secondary CAD", "Presenting Hb (g/dL)",
                "Ambient cold-exposure level", "Cold-avoidance counseling",
                "Scenario", "Horizon (d)"),
      Value = c(input$sev, input$secondary, input$base_hb, input$amb_cold,
                input$cold_avoidance, input$scenario, input$horizon)
    ) %>% datatable(rownames = FALSE, options = list(dom = "t"))
  })

  # --- Drug PK/PD effect ---
  output$pk_plot <- renderPlot({
    df <- results(); req(df)
    df %>% select(time, Rituximab_effect, Bendamustine_effect,
                  Sutimlimab_C1s_block, Pegcetacoplan_C3_block) %>%
      pivot_longer(-time) %>%
      ggplot(aes(time, value, colour = name)) +
        geom_line(linewidth = 0.8) +
        labs(x = "Time (days)", y = "Fractional drug effect (0-1)", colour = "Drug effect") +
        theme_minimal(base_size = 13)
  })

  # --- Complement pathway PD ---
  output$cp_plot <- renderPlot({
    df <- results(); req(df)
    df %>% select(time, Classical_pathway_flux, C4_level, C3_level) %>%
      pivot_longer(-time) %>%
      ggplot(aes(time, value, colour = name)) +
        geom_line(linewidth = 0.8) +
        facet_wrap(~name, scales = "free_y") +
        theme_minimal(base_size = 12) + theme(legend.position = "none")
  })
  output$mac_plot <- renderPlot({
    df <- results(); req(df)
    df %>% select(time, C3b_opsonization, MAC_level) %>%
      pivot_longer(-time) %>%
      ggplot(aes(time, value, colour = name)) +
        geom_line(linewidth = 0.9) +
        labs(x = "Day", y = "a.u.", colour = NULL) +
        theme_minimal(base_size = 13)
  })
  output$clone_igm_plot <- renderPlot({
    df <- results(); req(df)
    df %>% select(time, Clonal_burden, IgM_titer) %>%
      pivot_longer(-time) %>%
      ggplot(aes(time, value, colour = name)) +
        geom_line(linewidth = 0.9) +
        facet_wrap(~name, scales = "free_y") +
        theme_minimal(base_size = 12) + theme(legend.position = "none")
  })

  # --- Hemolysis biomarkers ---
  output$hb_retic_plot <- renderPlot({
    df <- results(); req(df)
    df %>% select(time, Hemoglobin, Reticulocyte_pct) %>%
      pivot_longer(-time) %>%
      ggplot(aes(time, value, colour = name)) +
        geom_line(linewidth = 0.9) +
        facet_wrap(~name, scales = "free_y") +
        theme_minimal(base_size = 12) + theme(legend.position = "none")
  })
  output$bili_ldh_plot <- renderPlot({
    df <- results(); req(df)
    df %>% select(time, Indirect_bilirubin, LDH_level) %>%
      pivot_longer(-time) %>%
      ggplot(aes(time, value, colour = name)) +
        geom_line(linewidth = 0.9) +
        facet_wrap(~name, scales = "free_y") +
        theme_minimal(base_size = 12) + theme(legend.position = "none")
  })
  output$hapto_plot <- renderPlot({
    df <- results(); req(df)
    ggplot(df, aes(time, Haptoglobin_level)) +
      geom_line(colour = "#8a6d3b", linewidth = 1) +
      labs(x = "Day", y = "Haptoglobin (mg/dL)") +
      theme_minimal(base_size = 13)
  })

  # --- Clinical endpoints table ---
  output$endpoint_table <- renderDT({
    df <- results(); req(df)
    tp <- unique(c(0, min(30, max(df$time)), min(90, max(df$time)), max(df$time)))
    df %>% filter(time %in% tp) %>%
      select(Day = time, Hemoglobin, Reticulocyte_pct, Indirect_bilirubin,
             LDH_level, Haptoglobin_level, Clonal_burden, IgM_titer) %>%
      mutate(across(-Day, ~ round(.x, 2))) %>%
      datatable(rownames = FALSE,
                options = list(dom = "t",
                               columnDefs = list(list(className = 'dt-center', targets = "_all"))))
  })

  # --- Compare ---
  output$compare_plot <- renderPlot({
    df <- all_results(); req(df)
    df %>% select(time, scenario, Hemoglobin, LDH_level, MAC_level, Clonal_burden) %>%
      pivot_longer(-c(time, scenario)) %>%
      ggplot(aes(time, value, colour = scenario)) +
        geom_line(linewidth = 0.6, alpha = 0.85) +
        facet_wrap(~name, scales = "free_y") +
        labs(x = "Day", y = NULL, colour = NULL) +
        theme_minimal(base_size = 12) +
        theme(legend.position = "bottom")
  })

  # --- Cold-exposure challenge ---
  output$cold_plot <- renderPlot({
    df <- results(); req(df)
    df %>% select(time, Hemoglobin, MAC_level, Agglutination_signal) %>%
      pivot_longer(-time) %>%
      ggplot(aes(time, value, colour = name)) +
        geom_line(linewidth = 0.9) +
        facet_wrap(~name, scales = "free_y") +
        theme_minimal(base_size = 12) + theme(legend.position = "none")
  })
}

shinyApp(ui, server)
