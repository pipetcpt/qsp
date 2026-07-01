## ============================================================================
## Acute Lymphoblastic Leukemia (ALL) QSP — Shiny Dashboard
## ----------------------------------------------------------------------------
## 8 tabs: Patient profile · Drug PK · Disease biology (blast/MRD) ·
##         Hematologic toxicity · Immunotherapy & CRS · Clinical endpoints ·
##         Scenario comparison · Biomarkers/Safety
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
## Builds the mrgsolve model by sourcing all_mrgsolve_model.R (which defines
## `mod` via mcode()) and caches the compiled object for the session.
get_model <- function() {
  if (!exists(".ALL_MOD", envir = .GlobalEnv) || is.null(.GlobalEnv$.ALL_MOD)) {
    model_env <- new.env()
    source("all_mrgsolve_model.R", local = model_env)
    assign(".ALL_MOD", model_env$mod, envir = .GlobalEnv)
  }
  .GlobalEnv$.ALL_MOD
}

## ---------- Scenario event builders ----------
build_induction_events <- function(bsa) {
  e_vcr <- ev(time = c(0, 7, 14, 21), amt = min(1.5 * bsa, 2.0), cmt = "VCR_cent")
  e_dex <- ev(time = seq(0, 27, by = 1), amt = 6 * bsa, cmt = "DEX_gut")
  e_asp <- ev(time = c(3, 17), amt = 2500 * bsa, cmt = "ASP_cent")
  e_it  <- ev(time = c(1, 8, 29), amt = 12, cmt = "MTX_CSF")
  e_vcr + e_dex + e_asp + e_it
}

build_events <- function(scenario, bsa) {
  if (scenario == "Untreated (natural history)") {
    ev(time = 0, amt = 0, cmt = "VCR_cent")
  } else if (scenario == "Pediatric SR B-ALL induction") {
    build_induction_events(bsa)
  } else if (scenario == "Pediatric HR B-ALL + HD-MTX") {
    build_induction_events(bsa) + ev(time = c(28, 42, 56, 70), amt = 5000 * bsa, cmt = "MTX_cent")
  } else if (scenario == "Adult Ph-neg B-ALL (hyper-CVAD-like)") {
    ev(time = c(0, 3), amt = 2.0, cmt = "VCR_cent") +
      ev(time = c(0,1,2,3,10,11,12,13), amt = 40, cmt = "DEX_gut") +
      ev(time = 15, amt = 1000 * bsa, cmt = "MTX_cent")
  } else if (scenario == "Ph+ ALL + Dasatinib") {
    ev(time = seq(0, 41, by = 1), amt = 40, cmt = "DEX_gut") +
      ev(time = seq(0, 83, by = 1), amt = 140, cmt = "TKI_gut")
  } else if (scenario == "Ph+ ALL (T315I) + Ponatinib") {
    ev(time = seq(0, 41, by = 1), amt = 40, cmt = "DEX_gut") +
      ev(time = seq(0, 83, by = 1), amt = 45, cmt = "TKI_gut")
  } else if (scenario == "R/R B-ALL + Blinatumomab (TOWER)") {
    ev(time = seq(0, 6.75, by = 0.25), amt = 9/4, cmt = "BLIN_cent") +
      ev(time = seq(7, 27.75, by = 0.25), amt = 28/4, cmt = "BLIN_cent")
  } else if (scenario == "R/R B-ALL + Inotuzumab (INO-VATE)") {
    ev(time = c(0, 7, 14), amt = c(0.8, 0.5, 0.5) * bsa, cmt = "INO_cent")
  } else if (scenario == "R/R B-ALL + CD19 CAR-T (ELIANA)") {
    ev(time = 0, amt = 5.0, cmt = "CART_blood")
  } else if (scenario == "Maintenance 6-MP (TPMT poor-metabolizer)") {
    ev(time = seq(0, 83, by = 1), amt = 50 * bsa, cmt = "MP_gut")
  } else {
    ev(time = 0, amt = 0, cmt = "VCR_cent")
  }
}

SWITCHES <- list(
  "Untreated (natural history)"              = list(),
  "Pediatric SR B-ALL induction"              = list(use_VCR=1, use_DEX=1, use_ASP=1, use_MTX=1),
  "Pediatric HR B-ALL + HD-MTX"               = list(use_VCR=1, use_DEX=1, use_ASP=1, use_MTX=1, k_grow=0.11),
  "Adult Ph-neg B-ALL (hyper-CVAD-like)"       = list(use_VCR=1, use_DEX=1, use_MTX=1, k_grow=0.10),
  "Ph+ ALL + Dasatinib"                        = list(use_DEX=1, use_TKI=1, T315I=0, k_grow=0.12),
  "Ph+ ALL (T315I) + Ponatinib"                = list(use_DEX=1, use_TKI=1, T315I=1, Ponatinib_active=1, k_grow=0.12),
  "R/R B-ALL + Blinatumomab (TOWER)"           = list(use_BLIN=1, k_grow=0.09),
  "R/R B-ALL + Inotuzumab (INO-VATE)"          = list(use_INO=1, k_grow=0.09),
  "R/R B-ALL + CD19 CAR-T (ELIANA)"            = list(use_CART=1, k_grow=0.09),
  "Maintenance 6-MP (TPMT poor-metabolizer)"   = list(use_MP=1, TPMT_mult=0.05, k_grow=0.03, BM_blast_init=2.0)
)

SCENARIO_LIST <- names(SWITCHES)

run_sim <- function(scenario, horizon, params) {
  mod <- get_model()
  ev_set <- build_events(scenario, params$bsa)
  sw <- SWITCHES[[scenario]]
  par <- c(
    list(
      BM_blast_init = params$bm_blast_init,
      k_grow        = params$k_grow,
      GR_resist     = params$gr_resist,
      TPMT_mult     = params$tpmt_mult,
      NUDT15_mult   = params$nudt15_mult,
      CD19_intact   = ifelse(params$cd19_intact, 1, 0),
      CD22_intact   = ifelse(params$cd22_intact, 1, 0),
      ASP_immunog   = params$asp_immunog
    ),
    sw
  )
  par <- par[!duplicated(names(par), fromLast = TRUE)]
  mod %>% param(par) %>% mrgsim(events = ev_set, end = horizon, delta = 0.25, obsonly = TRUE) %>%
    as_tibble() %>% mutate(scenario = scenario)
}

## ---------- UI ----------
ui <- dashboardPage(
  dashboardHeader(title = "ALL (Acute Lymphoblastic Leukemia) QSP"),
  dashboardSidebar(
    sidebarMenu(
      menuItem("1. Patient profile",       tabName = "patient", icon = icon("user")),
      menuItem("2. Drug PK",                tabName = "pk",      icon = icon("syringe")),
      menuItem("3. Disease biology",        tabName = "bio",     icon = icon("dna")),
      menuItem("4. Hematologic toxicity",   tabName = "heme",    icon = icon("tint")),
      menuItem("5. Immunotherapy & CRS",    tabName = "immuno",  icon = icon("shield-virus")),
      menuItem("6. Clinical endpoints",     tabName = "clin",    icon = icon("clipboard-check")),
      menuItem("7. Scenario comparison",    tabName = "compare", icon = icon("chart-line")),
      menuItem("8. Biomarkers / Safety",    tabName = "safety",  icon = icon("flask"))
    ),
    hr(),
    selectInput("scenario", "Treatment scenario:", SCENARIO_LIST,
                selected = "Pediatric SR B-ALL induction"),
    sliderInput("horizon", "Simulation horizon (days):", 14, 180, 84, step = 7),
    sliderInput("bsa", "Body surface area (m^2):", 0.5, 2.2, 0.9, step = 0.05),
    sliderInput("bm_blast_init", "Baseline BM blast (%):", 0, 95, 80),
    sliderInput("k_grow", "Net blast growth rate (1/day):", 0.02, 0.20, 0.08, step = 0.01),
    sliderInput("gr_resist", "Steroid-resistant fraction (NR3C1):", 0, 1, 0, step = 0.05),
    selectInput("tpmt_geno", "TPMT/NUDT15 genotype:",
                choices = c("Normal metabolizer" = "normal",
                            "Intermediate metabolizer" = "intermediate",
                            "Poor metabolizer" = "poor"), selected = "normal"),
    checkboxInput("cd19_intact", "CD19 antigen intact", TRUE),
    checkboxInput("cd22_intact", "CD22 antigen intact", TRUE),
    sliderInput("asp_immunog", "Anti-asparaginase Ab (silent inactivation):", 0, 1, 0, step = 0.05),
    actionButton("run", "Run simulation", icon = icon("play"),
                 style = "color:#fff;background:#1f6feb")
  ),
  dashboardBody(
    tabItems(
      tabItem("patient",
        fluidRow(
          box(width = 12, title = "Patient / disease profile summary", status = "primary",
              solidHeader = TRUE, DTOutput("patient_table"), br(),
              p(strong("Model scope:"), "Precursor B-/T-ALL leukemogenesis (BCR-ABL1,",
                "ETV6-RUNX1, KMT2A-r, hyperdiploidy, Ph-like, NOTCH1), bone-marrow niche,",
                "survival signaling, glucocorticoid apoptosis pathway, and 9 therapeutic",
                "agents spanning conventional chemotherapy, TKI, BiTE, ADC, and CD19 CAR-T."))
        ),
        fluidRow(
          box(width = 12, title = "About this dashboard", status = "info",
              p("Choose a treatment scenario and adjust patient/disease parameters in the",
                "left panel, then press", strong("Run simulation"), "to update all plots.",
                "Each tab isolates one layer of the QSP model: drug PK, leukemic blast/MRD",
                "dynamics, myelosuppression, immunotherapy toxicity (CRS), clinical",
                "response endpoints, cross-regimen comparison, and pharmacogenomic/safety",
                "biomarkers."))
        )
      ),

      tabItem("pk",
        fluidRow(
          box(width = 6, title = "Cytotoxic/targeted agent plasma concentration",
              status = "primary", solidHeader = TRUE, plotOutput("pk_plot1", height = 380)),
          box(width = 6, title = "Immunotherapy plasma / cellular exposure",
              status = "primary", solidHeader = TRUE, plotOutput("pk_plot2", height = 380))
        )
      ),

      tabItem("bio",
        fluidRow(
          box(width = 6, title = "BM blast % (logistic growth vs multi-agent kill)",
              plotOutput("blast_plot", 360)),
          box(width = 6, title = "MRD (log10 scale) & CNS sanctuary burden",
              plotOutput("mrd_plot", 360))
        )
      ),

      tabItem("heme",
        fluidRow(
          box(width = 6, title = "ANC — Friberg transit model",
              plotOutput("anc_plot", 360)),
          box(width = 6, title = "Platelets",
              plotOutput("plt_plot", 360))
        )
      ),

      tabItem("immuno",
        fluidRow(
          box(width = 6, title = "IL-6 (CRS marker)",
              plotOutput("il6_plot", 360)),
          box(width = 6, title = "CD19 CAR-T cell kinetics (blood vs marrow)",
              plotOutput("cart_plot", 360))
        )
      ),

      tabItem("clin",
        fluidRow(
          box(width = 12, title = "CR / MRD-negativity / relapse-risk summary",
              status = "success", solidHeader = TRUE, DTOutput("endpoint_table"))
        )
      ),

      tabItem("compare",
        fluidRow(
          box(width = 12, title = "Scenario comparison panel", status = "warning",
              solidHeader = TRUE,
              p("Runs all 10 built-in regimens with the current patient/disease profile."),
              actionButton("run_all", "Run all scenarios", icon = icon("rocket"),
                           style = "color:#fff;background:#0f5132"),
              br(), br(), plotOutput("compare_plot", height = 640))
        )
      ),

      tabItem("safety",
        fluidRow(
          box(width = 6, title = "Plasma asparagine (asparaginase PD)",
              plotOutput("asn_plot", 360)),
          box(width = 6, title = "TGN (thiopurine active metabolite, TPMT/NUDT15)",
              plotOutput("tgn_plot", 360))
        )
      )
    )
  )
)

## ---------- Server ----------
server <- function(input, output, session) {

  results <- reactiveVal(NULL)
  all_results <- reactiveVal(NULL)

  geno_mult <- function(g) switch(g, normal = 1.0, intermediate = 0.30, poor = 0.05)

  gather_params <- function() {
    list(
      bsa           = input$bsa,
      bm_blast_init = input$bm_blast_init,
      k_grow        = input$k_grow,
      gr_resist     = input$gr_resist,
      tpmt_mult     = geno_mult(input$tpmt_geno),
      nudt15_mult   = geno_mult(input$tpmt_geno),
      cd19_intact   = input$cd19_intact,
      cd22_intact   = input$cd22_intact,
      asp_immunog   = input$asp_immunog
    )
  }

  observeEvent(input$run, {
    showNotification("Running mrgsolve simulation...", type = "message", duration = 1)
    results(run_sim(input$scenario, input$horizon, gather_params()))
  }, ignoreNULL = FALSE)

  observeEvent(input$run_all, {
    showNotification("Running 10 scenarios...", type = "message", duration = 1)
    p <- gather_params()
    out <- lapply(SCENARIO_LIST, function(sc) run_sim(sc, input$horizon, p))
    all_results(bind_rows(out))
  })

  output$patient_table <- renderDT({
    tibble(
      Field = c("BSA (m^2)", "Baseline BM blast (%)", "Growth rate (1/day)",
                "Steroid-resistant fraction", "TPMT/NUDT15 genotype",
                "CD19 intact", "CD22 intact", "Anti-ASP antibody",
                "Scenario", "Horizon (d)"),
      Value = c(input$bsa, input$bm_blast_init, input$k_grow, input$gr_resist,
                input$tpmt_geno, input$cd19_intact, input$cd22_intact,
                input$asp_immunog, input$scenario, input$horizon)
    ) %>% datatable(rownames = FALSE, options = list(dom = "t"))
  })

  output$pk_plot1 <- renderPlot({
    df <- results(); if (is.null(df)) return(NULL)
    long <- pivot_longer(df, c(VCR_C1, DEX_C1, MTX_C1, TKI_C1, MP_C1),
                          names_to = "analyte", values_to = "conc")
    ggplot(long, aes(time, conc, color = analyte)) + geom_line(linewidth = 1) +
      facet_wrap(~analyte, scales = "free_y") +
      labs(x = "Time (days)", y = "Concentration") + theme_bw()
  })

  output$pk_plot2 <- renderPlot({
    df <- results(); if (is.null(df)) return(NULL)
    long <- pivot_longer(df, c(INO_C1, BLIN_C1, ASP_C1), names_to = "analyte", values_to = "conc")
    ggplot(long, aes(time, conc, color = analyte)) + geom_line(linewidth = 1) +
      facet_wrap(~analyte, scales = "free_y") +
      labs(x = "Time (days)", y = "Concentration") + theme_bw()
  })

  output$blast_plot <- renderPlot({
    df <- results(); if (is.null(df)) return(NULL)
    ggplot(df, aes(time, BM_blast_pct)) + geom_line(color = "#E41A1C", linewidth = 1.1) +
      geom_hline(yintercept = 5, linetype = "dashed", color = "darkgreen") +
      labs(x = "Time (days)", y = "BM Blast (%)") + theme_bw()
  })

  output$mrd_plot <- renderPlot({
    df <- results(); if (is.null(df)) return(NULL)
    long <- pivot_longer(df, c(MRD_log, CNS_blast_idx), names_to = "metric", values_to = "value")
    ggplot(long, aes(time, value, color = metric)) + geom_line(linewidth = 1) +
      facet_wrap(~metric, scales = "free_y") +
      labs(x = "Time (days)", y = "Value") + theme_bw()
  })

  output$anc_plot <- renderPlot({
    df <- results(); if (is.null(df)) return(NULL)
    ggplot(df, aes(time, ANC_out)) + geom_line(color = "#377EB8", linewidth = 1) +
      geom_hline(yintercept = 0.5, linetype = "dashed", color = "orange") +
      labs(x = "Time (days)", y = "ANC (x10^9/L)") + theme_bw()
  })

  output$plt_plot <- renderPlot({
    df <- results(); if (is.null(df)) return(NULL)
    ggplot(df, aes(time, PLT_out)) + geom_line(color = "#4DAF4A", linewidth = 1) +
      geom_hline(yintercept = 50, linetype = "dashed", color = "red") +
      labs(x = "Time (days)", y = "PLT (x10^9/L)") + theme_bw()
  })

  output$il6_plot <- renderPlot({
    df <- results(); if (is.null(df)) return(NULL)
    ggplot(df, aes(time, IL6_out)) + geom_line(color = "#FF7F00", linewidth = 1) +
      geom_hline(yintercept = c(20,100,500,1000), linetype = "dotted", color = "grey40") +
      labs(x = "Time (days)", y = "IL-6 (pg/mL)") + theme_bw()
  })

  output$cart_plot <- renderPlot({
    df <- results(); if (is.null(df)) return(NULL)
    long <- pivot_longer(df, c(CART_blood_out, CART_tissue_out), names_to = "compartment", values_to = "cells")
    ggplot(long, aes(time, cells, color = compartment)) + geom_line(linewidth = 1) +
      labs(x = "Time (days)", y = "CAR-T cells/uL") + theme_bw()
  })

  output$endpoint_table <- renderDT({
    df <- results(); if (is.null(df)) return(NULL)
    last_row <- df %>% slice_tail(n = 1)
    tibble(
      Endpoint = c("BM blast % (end)", "MRD log10 (end)", "CR achieved",
                   "MRD-negative", "ANC nadir", "PLT nadir", "Peak CRS grade"),
      Value = c(round(last_row$BM_blast_pct, 2), round(last_row$MRD_log, 2),
                ifelse(last_row$CR_flag == 1, "Yes", "No"),
                ifelse(last_row$MRDneg_flag == 1, "Yes", "No"),
                round(min(df$ANC_out), 2), round(min(df$PLT_out), 1),
                max(df$CRS_grade))
    ) %>% datatable(rownames = FALSE, options = list(dom = "t"))
  })

  output$compare_plot <- renderPlot({
    df <- all_results(); if (is.null(df)) return(NULL)
    p1 <- ggplot(df, aes(time, BM_blast_pct, color = scenario)) + geom_line(linewidth = 0.9) +
      labs(title = "BM Blast %", x = "Days", y = "%") + theme_bw() + theme(legend.position = "none")
    p2 <- ggplot(df, aes(time, MRD_log, color = scenario)) + geom_line(linewidth = 0.9) +
      labs(title = "MRD (log10)", x = "Days", y = "log10") + theme_bw() + theme(legend.position = "none")
    p3 <- ggplot(df, aes(time, ANC_out, color = scenario)) + geom_line(linewidth = 0.9) +
      labs(title = "ANC", x = "Days", y = "x10^9/L") + theme_bw() + theme(legend.position = "bottom", legend.title = element_blank())
    gridExtra_ok <- requireNamespace("patchwork", quietly = TRUE)
    if (gridExtra_ok) {
      library(patchwork)
      (p1 / p2 / p3)
    } else {
      p1
    }
  })

  output$asn_plot <- renderPlot({
    df <- results(); if (is.null(df)) return(NULL)
    ggplot(df, aes(time, ASP_C1)) + geom_line(color = "#984EA3", linewidth = 1) +
      labs(x = "Time (days)", y = "Asparaginase activity (IU/mL)") + theme_bw()
  })

  output$tgn_plot <- renderPlot({
    df <- results(); if (is.null(df)) return(NULL)
    ggplot(df, aes(time, TGN)) + geom_line(color = "#A65628", linewidth = 1) +
      labs(x = "Time (days)", y = "TGN (arbitrary units)") + theme_bw()
  })
}

shinyApp(ui, server)
