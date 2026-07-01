## ============================================================================
## Ménière's Disease QSP — Shiny Dashboard
## ----------------------------------------------------------------------------
## 8 tabs: Patient profile · Drug PK · Inner-ear PD/Hydrops · Clinical
##         endpoints (vertigo/hearing/tinnitus/DHI) · Scenario comparison ·
##         Biomarkers (ECoG) · Safety · References
## Dependencies: shiny, shinydashboard, mrgsolve, dplyr, tidyr, ggplot2, DT
## ----------------------------------------------------------------------------
## For research / education only. Not a substitute for clinical judgment.
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
  if (!exists(".MED_MOD", envir = .GlobalEnv) || is.null(.GlobalEnv$.MED_MOD)) {
    assign(".MED_MOD", mread_cache("med", project = "."), envir = .GlobalEnv)
  }
  .GlobalEnv$.MED_MOD
}

SCENARIO_LIST <- c(
  "Untreated natural history",
  "Betahistine standard-dose (16mg TID)",
  "Betahistine high-dose (48mg TID, BEMED)",
  "Diuretic (HCTZ+triamterene, once daily)",
  "Intratympanic dexamethasone series (wkly x4)",
  "Low-dose intratympanic gentamicin (titration)",
  "High-dose intratympanic gentamicin (fixed, destructive)",
  "Combination: betahistine high-dose + diuretic"
)

## ---------- Scenario event builder ----------
build_events <- function(scenario) {
  if (scenario == "Untreated natural history") {
    ev(amt = 0, cmt = "GUT_BETA", time = 0)

  } else if (scenario == "Betahistine standard-dose (16mg TID)") {
    ev(amt = 16, cmt = "GUT_BETA", ii = 8, addl = 999)

  } else if (scenario == "Betahistine high-dose (48mg TID, BEMED)") {
    ev(amt = 48, cmt = "GUT_BETA", ii = 8, addl = 999)

  } else if (scenario == "Diuretic (HCTZ+triamterene, once daily)") {
    ev(amt = 25, cmt = "GUT_DIUR", ii = 24, addl = 999)

  } else if (scenario == "Intratympanic dexamethasone series (wkly x4)") {
    ev(amt = 4, cmt = "ME_DEX", ii = 168, addl = 3)

  } else if (scenario == "Low-dose intratympanic gentamicin (titration)") {
    ev(amt = 10, cmt = "ME_GENT", ii = 720, addl = 2)

  } else if (scenario == "High-dose intratympanic gentamicin (fixed, destructive)") {
    ev(amt = 26.7, cmt = "ME_GENT", ii = 84, addl = 3)

  } else if (scenario == "Combination: betahistine high-dose + diuretic") {
    seq(ev(amt = 48, cmt = "GUT_BETA", ii = 8, addl = 999),
        ev(amt = 25, cmt = "GUT_DIUR", ii = 24, addl = 999))

  } else {
    ev(amt = 0, cmt = "GUT_BETA", time = 0)
  }
}

run_sim <- function(scenario, horizon_h, params) {
  mod <- get_model()
  ev_set <- build_events(scenario)
  par <- list(
    AGE         = params$age,
    WT          = params$wt,
    STAGE0      = params$stage0,
    DURATION_YR = params$duration_yr,
    BILATERAL   = params$bilateral
  )
  mod %>% param(par) %>% mrgsim(events = ev_set, end = horizon_h, delta = 6) %>%
    as_tibble() %>% mutate(scenario = scenario)
}

## ---------- UI ----------
ui <- dashboardPage(
  dashboardHeader(title = "Ménière's QSP"),
  dashboardSidebar(
    sidebarMenu(
      menuItem("1. Patient profile",         tabName = "patient", icon = icon("user")),
      menuItem("2. Drug PK",                 tabName = "pk",      icon = icon("syringe")),
      menuItem("3. Inner-ear PD / Hydrops",  tabName = "pd",      icon = icon("wave-square")),
      menuItem("4. Clinical endpoints",      tabName = "clin",    icon = icon("clipboard-list")),
      menuItem("5. Scenario comparison",     tabName = "compare", icon = icon("layer-group")),
      menuItem("6. Biomarkers (ECoG)",       tabName = "bio",     icon = icon("chart-line")),
      menuItem("7. Safety",                  tabName = "safety",  icon = icon("shield-alt")),
      menuItem("8. References",              tabName = "refs",    icon = icon("book"))
    ),
    hr(),
    selectInput("scenario", "Scenario:", SCENARIO_LIST, selected = "Betahistine standard-dose (16mg TID)"),
    sliderInput("horizon_mo", "Simulation horizon (months):", 1, 36, 12, step = 1),
    sliderInput("age",        "Age (years):", 18, 85, 55),
    sliderInput("wt",         "Weight (kg):", 40, 130, 70),
    sliderInput("stage0",     "Baseline AAO-HNS stage (I-IV):", 1, 4, 2, step = 1),
    sliderInput("duration_yr","Disease duration at baseline (years):", 0, 20, 3),
    selectInput("bilateral",  "Bilateral disease:", choices = c("Unilateral"=0, "Bilateral"=1), selected = 0),
    actionButton("run", "Run simulation", icon = icon("play"),
                 style = "color:#fff;background:#1f6feb")
  ),
  dashboardBody(
    tabItems(
      tabItem("patient",
        fluidRow(
          box(width = 12, title = "Patient profile summary", status = "primary",
              solidHeader = TRUE, DTOutput("patient_table"), br(),
              p(strong("AAO-HNS stage legend:"),
                "Stage I: 4-tone PTA <=25 dB; Stage II: 26-40 dB; Stage III: 41-70 dB; ",
                "Stage IV: >70 dB (Committee on Hearing and Equilibrium, 1995).")
          )
        ),
        fluidRow(
          box(width = 12, title = "About this dashboard", status = "info",
              p("Runs the mrgsolve QSP model for Ménière's disease: endolymphatic ",
                "sac dysfunction and stria vascularis ion-transport imbalance drive ",
                "endolymphatic hydrops, episodic Reissner's-membrane rupture with ",
                "perilymph K+ intoxication, and progressive cochlear/vestibular ",
                "hair-cell injury, producing episodic vertigo, fluctuating hearing ",
                "loss, tinnitus, aural fullness and handicap (DHI)."),
              p("Pick a treatment scenario and patient profile in the left panel, ",
                "then press ", strong("Run simulation"), ". Each tab focuses on a ",
                "model layer: PK, inner-ear pharmacodynamics/hydrops, clinical ",
                "endpoints, scenario comparison, biomarkers, and safety.")
          )
        )
      ),

      tabItem("pk",
        fluidRow(
          box(width = 12, title = "Plasma / perilymph drug concentrations",
              status = "primary", solidHeader = TRUE, plotOutput("pk_plot", height = 480))
        )
      ),

      tabItem("pd",
        fluidRow(
          box(width = 6, title = "Endolymphatic hydrops index", plotOutput("hydrops_plot", 360)),
          box(width = 6, title = "Perilymph K+ intoxication index (episodic)", plotOutput("kintox_plot", 360))
        ),
        fluidRow(
          box(width = 6, title = "Cochlear (OHC) & vestibular hair-cell viability", plotOutput("hc_plot", 360)),
          box(width = 6, title = "Central vestibular compensation index", plotOutput("comp_plot", 360))
        )
      ),

      tabItem("clin",
        fluidRow(
          box(width = 6, title = "Vertigo attack frequency (spells/month)", plotOutput("vert_plot", 360)),
          box(width = 6, title = "Pure-tone average (PTA, dB HL)", plotOutput("pta_plot", 360))
        ),
        fluidRow(
          box(width = 6, title = "Tinnitus severity (0-10)", plotOutput("tinn_plot", 360)),
          box(width = 6, title = "Dizziness Handicap Inventory (DHI, 0-100) & QoL", plotOutput("dhi_plot", 360))
        )
      ),

      tabItem("compare",
        fluidRow(
          box(width = 12, title = "Scenario comparison panel",
              status = "warning", solidHeader = TRUE,
              p("Runs all eight built-in scenarios with the current patient profile; ",
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

      tabItem("bio",
        fluidRow(
          box(width = 6, title = "ECoG SP/AP ratio (hydrops biomarker)", plotOutput("ecog_plot", 360)),
          box(width = 6, title = "AAO-HNS stage over time", plotOutput("stage_plot", 360))
        ),
        fluidRow(
          box(width = 12, title = "Tumarkin otolithic drop-attack risk flag", plotOutput("tumarkin_plot", 360))
        )
      ),

      tabItem("safety",
        fluidRow(
          box(width = 6, title = "Cumulative perilymph gentamicin exposure", plotOutput("gent_plot", 360)),
          box(width = 6, title = "Cochlear vs vestibular injury selectivity", plotOutput("selectivity_plot", 360))
        ),
        fluidRow(
          box(width = 12, title = "Safety flags at end of horizon", status = "danger",
              DTOutput("safety_table"))
        )
      ),

      tabItem("refs",
        fluidRow(
          box(width = 12, title = "Key references", status = "primary", solidHeader = TRUE,
              p("See ", code("med_references.md"), " in this directory for the full ",
                "curated list (30+ PubMed-linked references)."),
              tags$ul(
                tags$li("Committee on Hearing and Equilibrium, AAO-HNS. 1995 Otolaryngol Head Neck Surg — diagnostic/staging criteria."),
                tags$li("Basura GJ, et al. 2020 Otolaryngol Head Neck Surg — AAO-HNS Clinical Practice Guideline."),
                tags$li("Nauta JJP. 2014 Eur Arch Otorhinolaryngol — betahistine meta-analysis."),
                tags$li("Adrion C, et al. 2016 BMJ — BEMED trial (high-dose betahistine)."),
                tags$li("James AL, Burton MJ. 2001 Cochrane Database Syst Rev — diuretics."),
                tags$li("Phillips JS, Westerberg B. 2011 Cochrane Database Syst Rev — intratympanic steroids."),
                tags$li("Patel M, et al. 2016 Lancet — low- vs high-dose intratympanic gentamicin RCT."),
                tags$li("Sajjadi H, Paparella MM. 2008 Lancet — Ménière's disease review.")
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
    age = input$age, wt = input$wt, stage0 = input$stage0,
    duration_yr = input$duration_yr, bilateral = as.numeric(input$bilateral)
  )

  observeEvent(input$run, {
    showNotification("Running mrgsolve simulation…", type = "message", duration = 1)
    results(run_sim(input$scenario, input$horizon_mo * 730, cur_params()))
  }, ignoreNULL = FALSE)

  observeEvent(input$run_all, {
    showNotification("Running 8 scenarios…", type = "message", duration = 1)
    out <- lapply(SCENARIO_LIST, function(sc) run_sim(sc, input$horizon_mo * 730, cur_params()))
    all_results(bind_rows(out))
  })

  output$patient_table <- renderDT({
    bilat_lbl <- c("0"="Unilateral","1"="Bilateral")[[as.character(input$bilateral)]]
    tibble(
      Field = c("Age","Weight","Baseline stage","Disease duration (yr)",
                "Laterality","Scenario","Horizon (mo)"),
      Value = c(input$age, input$wt, paste0("Stage ", as.roman(input$stage0)),
                input$duration_yr, bilat_lbl, input$scenario, input$horizon_mo)
    ) %>% datatable(rownames = FALSE, options = list(dom = "t"))
  })

  output$pk_plot <- renderPlot({
    df <- results(); req(df)
    df %>% mutate(month = time/730) %>%
      select(month, conc_betahistine_parent, conc_betahistine_metabolite,
             conc_diuretic, peri_dexamethasone, peri_gentamicin) %>%
      pivot_longer(-month) %>%
      ggplot(aes(month, value, colour = name)) +
        geom_line(linewidth = 0.8) +
        facet_wrap(~name, scales = "free_y") +
        labs(x = "Time (months)", y = "Concentration (mg/L, arbitrary units for metabolite)", colour = "") +
        theme_minimal(base_size = 12) + theme(legend.position = "none")
  })

  output$hydrops_plot <- renderPlot({
    df <- results(); req(df)
    ggplot(df, aes(time/730, Hydrops_idx)) +
      geom_line(colour = "#264653", linewidth = 1) +
      geom_hline(yintercept = 1.0, lty = 2, colour = "grey50") +
      labs(x = "Month", y = "Hydrops index (1.0=baseline)") +
      theme_minimal(base_size = 13)
  })
  output$kintox_plot <- renderPlot({
    df <- results(); req(df)
    ggplot(df, aes(time/730, Kintox_idx)) +
      geom_line(colour = "#c0392b", linewidth = 0.8) +
      labs(x = "Month", y = "K+ intoxication index (episodic)") +
      theme_minimal(base_size = 13)
  })
  output$hc_plot <- renderPlot({
    df <- results(); req(df)
    df %>% mutate(month = time/730) %>%
      select(month, OHC_viability, Vestibular_viability) %>%
      pivot_longer(-month) %>%
      ggplot(aes(month, value, colour = name)) +
        geom_line(linewidth = 0.9) +
        scale_y_continuous(limits = c(0,1)) +
        labs(x = "Month", y = "Viability (1.0=intact)", colour = "") +
        theme_minimal(base_size = 13)
  })
  output$comp_plot <- renderPlot({
    df <- results(); req(df)
    ggplot(df, aes(time/730, Compensation_idx)) +
      geom_line(colour = "#2a9d8f", linewidth = 1) +
      labs(x = "Month", y = "Central vestibular compensation index") +
      theme_minimal(base_size = 13)
  })

  output$vert_plot <- renderPlot({
    df <- results(); req(df)
    ggplot(df, aes(time/730, VertigoFreq_permo)) +
      geom_line(colour = "#9b1c1c", linewidth = 1) +
      labs(x = "Month", y = "Vertigo spells/month") +
      theme_minimal(base_size = 13)
  })
  output$pta_plot <- renderPlot({
    df <- results(); req(df)
    ggplot(df, aes(time/730, PTA_dB)) +
      geom_line(colour = "#3d5a80", linewidth = 1) +
      labs(x = "Month", y = "PTA (dB HL, higher=worse)") +
      theme_minimal(base_size = 13)
  })
  output$tinn_plot <- renderPlot({
    df <- results(); req(df)
    ggplot(df, aes(time/730, Tinnitus_severity)) +
      geom_line(colour = "#e07a5f", linewidth = 1) +
      scale_y_continuous(limits = c(0,10)) +
      labs(x = "Month", y = "Tinnitus severity (0-10)") +
      theme_minimal(base_size = 13)
  })
  output$dhi_plot <- renderPlot({
    df <- results(); req(df)
    df %>% mutate(month = time/730) %>%
      select(month, DHI_score, QoL_score) %>%
      pivot_longer(-month) %>%
      ggplot(aes(month, value, colour = name)) +
        geom_line(linewidth = 0.9) +
        scale_y_continuous(limits = c(0,100)) +
        labs(x = "Month", y = "Score (0-100)", colour = "") +
        theme_minimal(base_size = 13)
  })

  output$compare_plot <- renderPlot({
    df <- all_results(); req(df)
    df %>% mutate(month = time/730) %>%
      ggplot(aes(month, VertigoFreq_permo, colour = scenario)) +
        geom_line(linewidth = 0.8) +
        labs(x = "Month", y = "Vertigo spells/month", colour = "Scenario") +
        theme_minimal(base_size = 13) +
        theme(legend.position = "bottom")
  })
  output$endpoint_table <- renderDT({
    df <- all_results(); req(df)
    horizon <- max(df$time)
    mo3 <- min(3*730, horizon)
    df %>% filter(time %in% c(0, mo3, horizon)) %>%
      mutate(Month = round(time/730, 1)) %>%
      select(scenario, Month, VertigoFreq_permo, PTA_dB, DHI_score, QoL_score) %>%
      distinct(scenario, Month, .keep_all = TRUE) %>%
      arrange(scenario, Month) %>%
      datatable(rownames = FALSE, options = list(pageLength = 30))
  })

  output$ecog_plot <- renderPlot({
    df <- results(); req(df)
    ggplot(df, aes(time/730, ECoG_SPAP_ratio)) +
      geom_line(colour = "#6d4c41", linewidth = 1) +
      geom_hline(yintercept = 0.40, lty = 2, colour = "grey50") +
      labs(x = "Month", y = "ECoG SP/AP ratio", caption = "Dashed line = common hydrops-suggestive threshold (~0.40)") +
      theme_minimal(base_size = 13)
  })
  output$stage_plot <- renderPlot({
    df <- results(); req(df)
    ggplot(df, aes(time/730, AAOHNS_stage_current)) +
      geom_step(colour = "#264653", linewidth = 1) +
      scale_y_continuous(limits = c(1,4), breaks = 1:4) +
      labs(x = "Month", y = "AAO-HNS stage") +
      theme_minimal(base_size = 13)
  })
  output$tumarkin_plot <- renderPlot({
    df <- results(); req(df)
    ggplot(df, aes(time/730, Tumarkin_risk_flag)) +
      geom_step(colour = "#c0392b", linewidth = 1) +
      scale_y_continuous(limits = c(0,1), breaks = c(0,1)) +
      labs(x = "Month", y = "Tumarkin drop-attack risk flag") +
      theme_minimal(base_size = 13)
  })

  output$gent_plot <- renderPlot({
    df <- results(); req(df)
    ggplot(df, aes(time/730, peri_gentamicin)) +
      geom_line(colour = "#8e44ad", linewidth = 1) +
      labs(x = "Month", y = "Perilymph gentamicin (mg/L)") +
      theme_minimal(base_size = 13)
  })
  output$selectivity_plot <- renderPlot({
    df <- results(); req(df)
    df %>% mutate(month = time/730) %>%
      select(month, OHC_viability, Vestibular_viability) %>%
      pivot_longer(-month) %>%
      ggplot(aes(month, 1 - value, colour = name)) +
        geom_line(linewidth = 0.9) +
        labs(x = "Month", y = "Fractional injury (1 - viability)", colour = "") +
        theme_minimal(base_size = 13)
  })
  output$safety_table <- renderDT({
    df <- results(); req(df)
    df %>% filter(time == max(time)) %>%
      select(peri_gentamicin, OHC_viability, Vestibular_viability,
             Tumarkin_risk_flag, PTA_dB, DHI_score) %>%
      datatable(rownames = FALSE, options = list(dom = "t"))
  })
}

shinyApp(ui, server)
