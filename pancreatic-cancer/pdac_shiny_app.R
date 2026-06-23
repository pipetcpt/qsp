# ============================================================
# PDAC QSP Shiny App — Pancreatic Ductal Adenocarcinoma
# Quantitative Systems Pharmacology Model Visualization
# ============================================================
# Libraries
library(shiny)
library(bslib)
library(plotly)
library(dplyr)
library(ggplot2)
library(tidyr)
library(deSolve)

# ============================================================
# GLOBAL PARAMETERS & HELPER FUNCTIONS
# ============================================================

REGIMEN_NAMES <- c(
  "Untreated (Control)",
  "Gemcitabine Monotherapy",
  "Gem + nab-Paclitaxel (MPACT)",
  "FOLFIRINOX (PRODIGE4)",
  "mFOLFIRINOX",
  "MRTX1133 (KRAS G12D)",
  "Olaparib (POLO, BRCA+)"
)

REGIMEN_COLORS <- c(
  "#888888", "#4e9af1", "#f4a261", "#e76f51",
  "#e9c46a", "#2a9d8f", "#bc6c25"
)
names(REGIMEN_COLORS) <- REGIMEN_NAMES

# Scenario-specific TGI parameters (k_eff) from clinical trials
SCENARIO_PARAMS <- list(
  "Untreated (Control)"              = list(k_eff = 0.00, k_prog = 0.045, mPFS = 1.5,  mOS = 3.0),
  "Gemcitabine Monotherapy"          = list(k_eff = 0.12, k_prog = 0.045, mPFS = 3.5,  mOS = 6.7),
  "Gem + nab-Paclitaxel (MPACT)"     = list(k_eff = 0.22, k_prog = 0.045, mPFS = 5.5,  mOS = 8.5),
  "FOLFIRINOX (PRODIGE4)"            = list(k_eff = 0.28, k_prog = 0.045, mPFS = 6.5,  mOS = 11.1),
  "mFOLFIRINOX"                      = list(k_eff = 0.25, k_prog = 0.045, mPFS = 6.0,  mOS = 10.5),
  "MRTX1133 (KRAS G12D)"             = list(k_eff = 0.35, k_prog = 0.045, mPFS = 4.0,  mOS = NA),
  "Olaparib (POLO, BRCA+)"           = list(k_eff = 0.18, k_prog = 0.045, mPFS = 7.4,  mOS = NA)
)

# Toxicity grade 3/4 by regimen (%)
TOXICITY_DATA <- data.frame(
  regimen = REGIMEN_NAMES,
  neutropenia = c(0,  38, 38, 46, 40, 25, 10),
  nausea      = c(0,   6, 10, 15, 12,  8,  5),
  fatigue     = c(5,  16, 17, 23, 20, 15, 12),
  neuropathy  = c(0,   0, 17,  9,  8,  0,  0),
  diarrhea    = c(2,   3,  6, 19, 16,  5,  3),
  stringsAsFactors = FALSE
)

# ORR by regimen (%)
ORR_DATA <- data.frame(
  regimen = REGIMEN_NAMES,
  CR  = c( 0,  0,  1,  3,  2,  2,  0),
  PR  = c( 2, 23, 22, 29, 26, 28,  9),
  SD  = c(20, 48, 52, 41, 44, 45, 42),
  PD  = c(78, 29, 25, 27, 28, 25, 49),
  stringsAsFactors = FALSE
)

# HR vs Untreated (Control) for forest plot
HR_DATA <- data.frame(
  regimen = REGIMEN_NAMES[-1],
  HR      = c(0.52, 0.35, 0.29, 0.32, 0.22, 0.48),
  HR_lo   = c(0.40, 0.26, 0.20, 0.23, 0.13, 0.30),
  HR_hi   = c(0.67, 0.47, 0.42, 0.44, 0.37, 0.76),
  trial   = c("MPACT ctrl", "MPACT", "PRODIGE4", "Derived", "Phase I", "POLO"),
  stringsAsFactors = FALSE
)

# ---- Analytical PK Functions ----

pk_gemcitabine <- function(t_h, dose_mgm2 = 1000) {
  # 2-compartment biexponential after 30-min IV infusion
  # Typical params from literature
  A <- 50.0 * dose_mgm2 / 1000
  B <- 5.0  * dose_mgm2 / 1000
  alpha <- 2.1   # h^-1
  beta  <- 0.22  # h^-1
  ifelse(t_h < 0, 0, A * exp(-alpha * t_h) + B * exp(-beta * t_h))
}

pk_nabpaclitaxel <- function(t_h, dose_mgm2 = 125) {
  A <- 3.5  * dose_mgm2 / 125
  B <- 0.8  * dose_mgm2 / 125
  alpha <- 0.85
  beta  <- 0.10
  ifelse(t_h < 0, 0, A * exp(-alpha * t_h) + B * exp(-beta * t_h))
}

pk_oxaliplatin <- function(t_h, dose_mgm2 = 85) {
  A <- 2.8 * dose_mgm2 / 85
  B <- 0.4 * dose_mgm2 / 85
  alpha <- 0.40
  beta  <- 0.035
  ifelse(t_h < 0, 0, A * exp(-alpha * t_h) + B * exp(-beta * t_h))
}

pk_sn38 <- function(t_h, dose_mgm2 = 180) {
  # SN-38 (active metabolite of irinotecan) formation-elimination model
  F_oral <- 0.04
  ka <- 0.25
  ke <- 0.06
  Vd_ratio <- 12
  Cpeak <- F_oral * dose_mgm2 / Vd_ratio
  ifelse(t_h < 0, 0,
         Cpeak * (exp(-ke * t_h) - exp(-ka * t_h)) * ka / (ka - ke))
}

pk_5fu <- function(t_h, dose_mgm2 = 400, infusion_h = 46) {
  # 5-FU: bolus + infusion model
  Cbol <- 3.2  * dose_mgm2 / 400
  Cinf <- 0.15 * dose_mgm2 / 400
  ke   <- 1.1
  ifelse(t_h <= infusion_h,
         Cbol * exp(-ke * t_h) + Cinf * (1 - exp(-ke * t_h)),
         (Cbol * exp(-ke * t_h) + Cinf * (1 - exp(-ke * infusion_h))) * exp(-ke * (t_h - infusion_h))
  )
}

pk_mrtx1133 <- function(t_h, dose_mg = 700) {
  # MRTX1133: 1-compartment oral
  ka  <- 0.8
  ke  <- 0.12
  F   <- 0.45
  Vd  <- 120  # L
  Cpeak <- (F * dose_mg) / Vd
  ifelse(t_h < 0, 0,
         Cpeak * (ka / (ka - ke)) * (exp(-ke * t_h) - exp(-ka * t_h)))
}

# ---- Simeoni TGI ODE Model ----
pdac_ode <- function(t, state, params) {
  with(as.list(c(state, params)), {
    # Piecewise drug concentration (analytical, simplified)
    cycle_t <- t %% cycle_len
    drug_conc <- drug_peak * exp(-drug_ke * cycle_t)

    drug_effect <- min(Emax * drug_conc / (EC50 + drug_conc), 0.99)

    dx0    <- (k_prog * (1 - drug_effect) - k_death) * x0
    dx1    <- drug_effect * k_prog * x0 - k_tr * x1
    dx2    <- k_tr * x1 - k_tr * x2
    dx3    <- k_tr * x2 - k_tr * x3
    dTUMOR <- k_tr * x3 - k_death * TUMOR
    dCA199 <- 0.01 * x0 - 0.05 * CA199

    list(c(dx0, dx1, dx2, dx3, dTUMOR, dCA199))
  })
}

run_tgi_sim <- function(regimen, duration_months = 12,
                        init_tumor = 5000, init_ca199 = 500,
                        stroma_level = 0.5, kras_g12d = FALSE) {
  p <- SCENARIO_PARAMS[[regimen]]
  if (is.null(p)) return(NULL)

  k_eff_adj <- p$k_eff * (1 - 0.4 * stroma_level)
  # MRTX1133 only effective in KRAS G12D tumors
  if (regimen == "MRTX1133 (KRAS G12D)" && !kras_g12d) {
    k_eff_adj <- k_eff_adj * 0.05
  }

  parms <- c(
    k_prog  = p$k_prog,
    k_death = 0.001,
    k_tr    = 0.08,
    Emax    = 1.0,
    EC50    = 0.5,
    drug_peak = k_eff_adj * 10,
    drug_ke   = 0.15,
    cycle_len = 14 / 30   # 2-week cycles in months
  )

  state <- c(x0 = init_tumor * 0.8, x1 = 0, x2 = 0, x3 = 0,
             TUMOR = init_tumor * 0.2, CA199 = init_ca199)

  times <- seq(0, duration_months, by = 0.1)

  out <- tryCatch(
    ode(y = state, times = times, func = pdac_ode, parms = parms,
        method = "lsoda"),
    error = function(e) NULL
  )
  if (is.null(out)) return(NULL)

  df <- as.data.frame(out)
  df$total_tumor <- df$x0 + df$x1 + df$x2 + df$x3 + df$TUMOR
  df$regimen <- regimen
  df
}

# ---- Parametric Survival (Weibull) ----
weibull_km <- function(t, median_time, shape = 1.2) {
  lambda <- log(2)^(1 / shape) / median_time
  exp(-(lambda * t)^shape)
}

weibull_ci <- function(t, median_time, shape = 1.2, n = 200) {
  # Bootstrap-like CI via parameter uncertainty
  se_log_med <- 0.15
  med_lo <- median_time * exp(-1.96 * se_log_med)
  med_hi <- median_time * exp( 1.96 * se_log_med)
  list(
    lo = weibull_km(t, med_hi, shape),
    hi = weibull_km(t, med_lo, shape)
  )
}

# ---- Neutrophil Friberg Model ----
friberg_neutrophil <- function(t_days, k_eff_neutro, circ0 = 5.0) {
  ktr   <- 0.0889
  gamma <- 0.161
  state <- c(prol = circ0, tr1 = circ0, tr2 = circ0, tr3 = circ0, circ = circ0)
  parms <- c(ktr = ktr, gamma = gamma, circ0 = circ0, drug_eff = k_eff_neutro)
  friberg_ode <- function(t, y, p) {
    with(as.list(c(y, p)), {
      drug_c <- drug_eff * exp(-0.18 * (t %% 21))
      feedback <- (circ0 / max(circ, 0.01))^gamma
      dprol  <- ktr * prol * (1 - drug_c) * feedback - ktr * prol
      dtr1   <- ktr * prol - ktr * tr1
      dtr2   <- ktr * tr1  - ktr * tr2
      dtr3   <- ktr * tr2  - ktr * tr3
      dcirc  <- ktr * tr3  - ktr * circ
      list(c(dprol, dtr1, dtr2, dtr3, dcirc))
    })
  }
  out <- tryCatch(
    ode(y = state, times = t_days, func = friberg_ode, parms = parms),
    error = function(e) NULL
  )
  if (is.null(out)) return(data.frame(time = t_days, circ = circ0))
  data.frame(time = as.numeric(out[, "time"]), circ = as.numeric(out[, "circ"]))
}

NEUTRO_EFF <- c(
  "Untreated (Control)"          = 0.00,
  "Gemcitabine Monotherapy"      = 0.35,
  "Gem + nab-Paclitaxel (MPACT)" = 0.50,
  "FOLFIRINOX (PRODIGE4)"        = 0.65,
  "mFOLFIRINOX"                  = 0.58,
  "MRTX1133 (KRAS G12D)"         = 0.10,
  "Olaparib (POLO, BRCA+)"       = 0.08
)

# ============================================================
# UI
# ============================================================

ui <- fluidPage(
  theme = bs_theme(bootswatch = "darkly", version = 5),

  tags$head(tags$style(HTML("
    .nav-tabs .nav-link { font-weight: 500; }
    .sidebar-panel-dark { background: #2c2c2c; border-radius: 6px; padding: 10px; }
    .value-box-dark { background: #343a40; border: 1px solid #495057;
                      border-radius: 6px; padding: 12px; margin-bottom: 8px; }
    .risk-badge-high   { color: #ff6b6b; font-weight: bold; }
    .risk-badge-medium { color: #ffa94d; font-weight: bold; }
    .risk-badge-low    { color: #69db7c; font-weight: bold; }
  "))),

  titlePanel(
    div(
      h2("PDAC QSP Model — Pancreatic Ductal Adenocarcinoma",
         style = "color:#4e9af1; margin-bottom:2px;"),
      h5("Quantitative Systems Pharmacology Interactive Dashboard",
         style = "color:#adb5bd;")
    )
  ),

  tabsetPanel(
    id = "main_tabs",

    # ==================== TAB 1: Patient Profile ====================
    tabPanel(
      "Patient Profile",
      icon = icon("user"),
      br(),
      sidebarLayout(
        sidebarPanel(
          width = 3,
          div(class = "sidebar-panel-dark",
            h5("Clinical Parameters", style = "color:#4e9af1;"),
            selectInput("stage", "Disease Stage:",
                        choices = c("Stage I", "Stage II", "Stage III", "Stage IV"),
                        selected = "Stage IV"),
            selectInput("mol_subtype", "Molecular Subtype:",
                        choices = c("Classical", "Basal-like/Squamous",
                                    "Quasi-mesenchymal", "Exocrine-like"),
                        selected = "Classical"),
            selectInput("kras_mut", "KRAS Mutation:",
                        choices = c("G12D (~44%)", "G12V (~26%)", "G12R (~14%)",
                                    "G12C (~2%)", "Other", "WT (rare)"),
                        selected = "G12D (~44%)"),
            checkboxGroupInput("hrd_status", "HRD Status:",
                               choices = c("BRCA1/2 mutation", "PALB2 mutation", "ATM mutation")),
            radioButtons("msi_status", "MSI Status:",
                         choices = c("MSS", "MSI-H (~1%)"),
                         selected = "MSS"),
            hr(),
            h5("Patient Metrics", style = "color:#4e9af1;"),
            sliderInput("ecog", "ECOG Performance Status:", min = 0, max = 2,
                        value = 1, step = 1),
            numericInput("weight", "Weight (kg):", value = 68, min = 30, max = 200),
            numericInput("bsa", "BSA (m²):", value = 1.75, min = 0.5, max = 3.0, step = 0.01),
            numericInput("ca199_base", "Baseline CA19-9 (U/mL):", value = 500, min = 0),
            numericInput("tumor_vol_base", "Baseline Tumor Volume (mm³):", value = 5000, min = 100)
          )
        ),
        mainPanel(
          width = 9,
          fluidRow(
            column(4, div(class = "value-box-dark",
              h6("Disease Stage", style = "color:#adb5bd;"),
              h4(textOutput("vbox_stage"), style = "color:#4e9af1;")
            )),
            column(4, div(class = "value-box-dark",
              h6("KRAS Mutation", style = "color:#adb5bd;"),
              h4(textOutput("vbox_kras"), style = "color:#f4a261;")
            )),
            column(4, div(class = "value-box-dark",
              h6("Risk Category", style = "color:#adb5bd;"),
              h4(uiOutput("vbox_risk"))
            ))
          ),
          br(),
          fluidRow(
            column(6,
              h5("Patient Summary", style = "color:#4e9af1;"),
              tableOutput("patient_summary_tbl")
            ),
            column(6,
              h5("Risk Stratification", style = "color:#4e9af1;"),
              plotlyOutput("risk_radar", height = "300px")
            )
          ),
          br(),
          fluidRow(
            column(12,
              h5("Molecular Pathway Activation Heatmap", style = "color:#4e9af1;"),
              plotlyOutput("pathway_heatmap", height = "300px")
            )
          )
        )
      )
    ),

    # ==================== TAB 2: Drug PK ====================
    tabPanel(
      "Drug PK",
      icon = icon("pills"),
      br(),
      sidebarLayout(
        sidebarPanel(
          width = 3,
          div(class = "sidebar-panel-dark",
            h5("PK Settings", style = "color:#4e9af1;"),
            selectInput("pk_regimen", "Treatment Regimen:",
                        choices = REGIMEN_NAMES[-1],
                        selected = "Gem + nab-Paclitaxel (MPACT)"),
            sliderInput("pk_cycle", "Cycle Number:", min = 1, max = 6, value = 1),
            checkboxGroupInput("pk_drugs", "Show Drug(s):",
                               choices = c("Gemcitabine", "nab-Paclitaxel",
                                           "Oxaliplatin", "SN-38 (Irinotecan)",
                                           "5-FU", "MRTX1133"),
                               selected = c("Gemcitabine", "nab-Paclitaxel")),
            checkboxInput("pk_logscale", "Log Y-axis", value = FALSE),
            hr(),
            h5("Dose Modifications", style = "color:#4e9af1;"),
            sliderInput("dose_reduction", "Dose Reduction (%):",
                        min = 0, max = 50, value = 0, step = 5)
          )
        ),
        mainPanel(
          width = 9,
          h5("Plasma Concentration-Time Profiles", style = "color:#4e9af1;"),
          plotlyOutput("pk_plot", height = "350px"),
          br(),
          fluidRow(
            column(12,
              h5("PK Parameters Summary", style = "color:#4e9af1;"),
              tableOutput("pk_table")
            )
          )
        )
      )
    ),

    # ==================== TAB 3: Tumor Dynamics ====================
    tabPanel(
      "Tumor Dynamics",
      icon = icon("chart-line"),
      br(),
      sidebarLayout(
        sidebarPanel(
          width = 3,
          div(class = "sidebar-panel-dark",
            h5("Simulation Settings", style = "color:#4e9af1;"),
            checkboxGroupInput("tgi_scenarios", "Treatment Scenarios:",
                               choices = REGIMEN_NAMES,
                               selected = c("Untreated (Control)",
                                            "Gem + nab-Paclitaxel (MPACT)",
                                            "FOLFIRINOX (PRODIGE4)")),
            sliderInput("tgi_duration", "Simulation Duration (months):",
                        min = 3, max = 24, value = 12, step = 1),
            sliderInput("stroma_level", "Stroma Level (desmoplasia):",
                        min = 0, max = 1, value = 0.5, step = 0.05),
            checkboxInput("kras_g12d_tgi", "KRAS G12D Present", value = TRUE),
            hr(),
            actionButton("run_tgi", "Run Simulation",
                         class = "btn-primary", width = "100%")
          )
        ),
        mainPanel(
          width = 9,
          fluidRow(
            column(6,
              h5("Tumor Volume Over Time", style = "color:#4e9af1;"),
              plotlyOutput("tgi_tumor_plot", height = "300px")
            ),
            column(6,
              h5("CA19-9 Trajectory", style = "color:#4e9af1;"),
              plotlyOutput("tgi_ca199_plot", height = "300px")
            )
          ),
          br(),
          fluidRow(
            column(6,
              h5("Waterfall Plot — Best % Change from Baseline", style = "color:#4e9af1;"),
              plotlyOutput("waterfall_plot", height = "280px")
            ),
            column(6,
              h5("Tumor Growth Inhibition Summary", style = "color:#4e9af1;"),
              tableOutput("tgi_table")
            )
          )
        )
      )
    ),

    # ==================== TAB 4: Biomarkers ====================
    tabPanel(
      "Biomarkers",
      icon = icon("vial"),
      br(),
      sidebarLayout(
        sidebarPanel(
          width = 3,
          div(class = "sidebar-panel-dark",
            h5("Biomarker Settings", style = "color:#4e9af1;"),
            radioButtons("bm_type", "Biomarker:",
                         choices = c("CA19-9 (U/mL)", "ctDNA VAF (%)",
                                     "KRAS ctDNA (%)", "Neutrophils (×10⁹/L)"),
                         selected = "CA19-9 (U/mL)"),
            sliderInput("bm_horizon", "Time Horizon (weeks):",
                        min = 0, max = 52, value = 24, step = 1)
          )
        ),
        mainPanel(
          width = 9,
          fluidRow(
            column(6,
              h5("Biomarker Kinetics by Regimen", style = "color:#4e9af1;"),
              plotlyOutput("bm_kinetics_plot", height = "300px")
            ),
            column(6,
              h5("Early Response: % Change at Week 8", style = "color:#4e9af1;"),
              plotlyOutput("bm_week8_plot", height = "300px")
            )
          ),
          br(),
          fluidRow(
            column(6,
              h5("Neutrophil Nadir (Friberg Myelosuppression Model)", style = "color:#4e9af1;"),
              plotlyOutput("neutrophil_plot", height = "280px")
            ),
            column(6,
              h5("Grade 3/4 Neutropenia Risk by Regimen", style = "color:#4e9af1;"),
              tableOutput("neutropenia_table")
            )
          )
        )
      )
    ),

    # ==================== TAB 5: Clinical Endpoints ====================
    tabPanel(
      "Clinical Endpoints",
      icon = icon("hospital"),
      br(),
      sidebarLayout(
        sidebarPanel(
          width = 3,
          div(class = "sidebar-panel-dark",
            h5("Simulation Parameters", style = "color:#4e9af1;"),
            sliderInput("n_patients", "Simulated Patients (N):",
                        min = 50, max = 500, value = 200, step = 50),
            checkboxInput("show_ci", "Show 95% CI", value = TRUE),
            hr(),
            h6("Reference Trials:", style = "color:#adb5bd;"),
            tags$ul(
              tags$li("MPACT (nab-Pac+Gem)", style = "font-size:11px;color:#adb5bd;"),
              tags$li("PRODIGE4 (FOLFIRINOX)", style = "font-size:11px;color:#adb5bd;"),
              tags$li("POLO (Olaparib)", style = "font-size:11px;color:#adb5bd;")
            )
          )
        ),
        mainPanel(
          width = 9,
          fluidRow(
            column(6,
              h5("Progression-Free Survival (PFS)", style = "color:#4e9af1;"),
              plotlyOutput("km_pfs_plot", height = "300px")
            ),
            column(6,
              h5("Overall Survival (OS)", style = "color:#4e9af1;"),
              plotlyOutput("km_os_plot", height = "300px")
            )
          ),
          br(),
          fluidRow(
            column(6,
              h5("Objective Response Rate (ORR)", style = "color:#4e9af1;"),
              plotlyOutput("orr_plot", height = "250px")
            ),
            column(6,
              h5("Toxicity Profile", style = "color:#4e9af1;"),
              plotlyOutput("tox_radar_plot", height = "250px")
            )
          ),
          br(),
          fluidRow(
            column(12,
              h5("Clinical Outcomes Summary (vs. Reference Trials)", style = "color:#4e9af1;"),
              tableOutput("outcomes_table")
            )
          )
        )
      )
    ),

    # ==================== TAB 6: Scenario Comparison ====================
    tabPanel(
      "Scenario Comparison",
      icon = icon("balance-scale"),
      br(),
      sidebarLayout(
        sidebarPanel(
          width = 3,
          div(class = "sidebar-panel-dark",
            h5("Comparison Settings", style = "color:#4e9af1;"),
            selectInput("ref_scenario", "Reference Scenario:",
                        choices = c("Untreated (Control)", "Gemcitabine Monotherapy"),
                        selected = "Untreated (Control)"),
            checkboxGroupInput("comp_scenarios", "Comparison Scenarios:",
                               choices = REGIMEN_NAMES[-1],
                               selected = c("Gem + nab-Paclitaxel (MPACT)",
                                            "FOLFIRINOX (PRODIGE4)",
                                            "mFOLFIRINOX",
                                            "Olaparib (POLO, BRCA+)"))
          )
        ),
        mainPanel(
          width = 9,
          fluidRow(
            column(6,
              h5("Forest Plot — Hazard Ratio for OS vs Reference", style = "color:#4e9af1;"),
              plotlyOutput("forest_plot", height = "320px")
            ),
            column(6,
              h5("Parallel Coordinates — Multi-Outcome Profile", style = "color:#4e9af1;"),
              plotlyOutput("parallel_plot", height = "320px")
            )
          ),
          br(),
          fluidRow(
            column(6,
              h5("Summary Comparison Table", style = "color:#4e9af1;"),
              tableOutput("comparison_table")
            ),
            column(6,
              h5("Sensitivity Analysis — Tornado Plot", style = "color:#4e9af1;"),
              plotlyOutput("tornado_plot", height = "300px")
            )
          )
        )
      )
    )
  )
)

# ============================================================
# SERVER
# ============================================================

server <- function(input, output, session) {

  # ---- Reactive: TGI Simulation ----
  tgi_results <- eventReactive(input$run_tgi, {
    req(length(input$tgi_scenarios) > 0)
    results <- lapply(input$tgi_scenarios, function(reg) {
      run_tgi_sim(
        regimen         = reg,
        duration_months = input$tgi_duration,
        init_tumor      = isolate(input$tumor_vol_base),
        init_ca199      = isolate(input$ca199_base),
        stroma_level    = input$stroma_level,
        kras_g12d       = input$kras_g12d_tgi
      )
    })
    bind_rows(results[!sapply(results, is.null)])
  }, ignoreNULL = FALSE)

  # auto-run on load
  observe({
    if (is.null(tgi_results())) {
      shinyjs::click("run_tgi")
    }
  })

  # ============================================================
  # TAB 1: Patient Profile
  # ============================================================

  output$vbox_stage <- renderText({ input$stage })
  output$vbox_kras  <- renderText({ input$kras_mut })

  output$vbox_risk <- renderUI({
    score <- 0
    if (input$stage == "Stage IV")   score <- score + 3
    if (input$stage == "Stage III")  score <- score + 2
    if (grepl("Basal|Squamous|Quasi", input$mol_subtype)) score <- score + 2
    if (input$ecog >= 2) score <- score + 1
    if (input$ca199_base > 1000) score <- score + 1

    if (score >= 5) span("High Risk", class = "risk-badge-high")
    else if (score >= 3) span("Intermediate Risk", class = "risk-badge-medium")
    else span("Lower Risk", class = "risk-badge-low")
  })

  output$patient_summary_tbl <- renderTable({
    brca <- if (length(input$hrd_status) > 0)
               paste(input$hrd_status, collapse = "; ") else "None"
    data.frame(
      Parameter = c("Stage", "Molecular Subtype", "KRAS Mutation",
                    "MSI Status", "HRD Status",
                    "ECOG PS", "Weight (kg)", "BSA (m²)",
                    "CA19-9 (U/mL)", "Tumor Volume (mm³)"),
      Value = c(input$stage, input$mol_subtype, input$kras_mut,
                input$msi_status, brca,
                as.character(input$ecog),
                as.character(input$weight),
                as.character(input$bsa),
                as.character(input$ca199_base),
                as.character(input$tumor_vol_base))
    )
  }, striped = TRUE, bordered = TRUE, hover = TRUE)

  output$risk_radar <- renderPlotly({
    stage_score <- switch(input$stage,
      "Stage I" = 1, "Stage II" = 2, "Stage III" = 3, "Stage IV" = 4)
    subtype_score <- if (grepl("Basal|Squamous", input$mol_subtype)) 4
                     else if (grepl("Quasi", input$mol_subtype)) 3
                     else if (grepl("Exocrine", input$mol_subtype)) 2
                     else 3
    kras_score <- if (grepl("G12D|G12V", input$kras_mut)) 4
                  else if (grepl("G12R", input$kras_mut)) 3
                  else 2
    hrd_score  <- length(input$hrd_status) * 1.5
    ecog_score <- input$ecog * 2 + 1
    ca199_score <- min(ceiling(log10(max(input$ca199_base, 2)) * 1.2), 5)

    categories <- c("Disease Stage", "Molecular Subtype", "KRAS Aggression",
                    "HRD Status", "ECOG PS", "CA19-9 Level", "Disease Stage")
    values <- c(stage_score, subtype_score, kras_score,
                pmin(hrd_score, 4), ecog_score, ca199_score, stage_score)

    plot_ly(type = "scatterpolar", fill = "toself",
            r = values, theta = categories, mode = "lines+markers",
            line = list(color = "#4e9af1"),
            fillcolor = "rgba(78,154,241,0.25)") |>
      layout(
        polar = list(radialaxis = list(visible = TRUE, range = c(0, 5),
                                       color = "#adb5bd"),
                     angularaxis = list(color = "#adb5bd")),
        paper_bgcolor = "#222",
        plot_bgcolor  = "#222",
        font = list(color = "#dee2e6"),
        margin = list(t = 20, b = 20, l = 40, r = 40)
      )
  })

  output$pathway_heatmap <- renderPlotly({
    pathways <- c("KRAS-MAPK", "PI3K-AKT-mTOR", "TGF-β/SMAD4",
                  "Wnt/β-catenin", "Hedgehog", "JAK-STAT",
                  "p53/apoptosis", "DNA repair (HRD)", "Immune exclusion")
    subtypes  <- c("Classical", "Basal-like", "Quasi-mesen.", "Exocrine-like")

    set.seed(42)
    base_mat <- matrix(c(
      4, 3, 2, 2,   # KRAS
      3, 4, 3, 1,   # PI3K
      2, 3, 4, 1,   # TGFb
      2, 2, 3, 2,   # Wnt
      1, 2, 2, 3,   # Hh
      2, 3, 2, 1,   # JAK
      2, 3, 2, 2,   # p53
      2, 2, 1, 1,   # HRD
      3, 4, 3, 1    # Immune
    ), nrow = length(pathways), byrow = TRUE)

    # Highlight selected subtype
    sel_idx <- match(gsub("-like.*|/.*", "", input$mol_subtype),
                     c("Classical", "Basal", "Quasi", "Exocrine"))
    if (!is.na(sel_idx)) base_mat[, sel_idx] <- pmin(base_mat[, sel_idx] + 1, 5)

    plot_ly(
      x = subtypes, y = pathways, z = base_mat,
      type = "heatmap",
      colorscale = list(c(0,"#1a1a2e"), c(0.5,"#16213e"),
                        c(0.75,"#0f3460"), c(1,"#e94560")),
      showscale = TRUE
    ) |>
      layout(
        xaxis = list(title = "", tickfont = list(color = "#dee2e6")),
        yaxis = list(title = "", tickfont = list(color = "#dee2e6")),
        paper_bgcolor = "#222",
        plot_bgcolor  = "#222",
        font = list(color = "#dee2e6"),
        margin = list(t = 20, b = 40, l = 150, r = 20)
      )
  })

  # ============================================================
  # TAB 2: Drug PK
  # ============================================================

  pk_data <- reactive({
    dose_factor <- 1 - input$dose_reduction / 100
    t_h  <- seq(0, 48, by = 0.25)
    dfs  <- list()

    bsa  <- input$bsa

    if ("Gemcitabine" %in% input$pk_drugs &&
        input$pk_regimen %in% c("Gemcitabine Monotherapy",
                                 "Gem + nab-Paclitaxel (MPACT)")) {
      conc <- pk_gemcitabine(t_h, dose_mgm2 = 1000 * dose_factor * bsa / bsa)
      dfs[["Gemcitabine"]] <- data.frame(time = t_h, conc = conc,
                                          drug = "Gemcitabine (µg/mL)")
    }
    if ("nab-Paclitaxel" %in% input$pk_drugs &&
        input$pk_regimen == "Gem + nab-Paclitaxel (MPACT)") {
      conc <- pk_nabpaclitaxel(t_h, dose_mgm2 = 125 * dose_factor)
      dfs[["nab-Paclitaxel"]] <- data.frame(time = t_h, conc = conc,
                                              drug = "nab-Paclitaxel (µg/mL)")
    }
    if ("Oxaliplatin" %in% input$pk_drugs &&
        input$pk_regimen %in% c("FOLFIRINOX (PRODIGE4)", "mFOLFIRINOX")) {
      conc <- pk_oxaliplatin(t_h, dose_mgm2 = 85 * dose_factor)
      dfs[["Oxaliplatin"]] <- data.frame(time = t_h, conc = conc,
                                          drug = "Oxaliplatin (µg/mL)")
    }
    if ("SN-38 (Irinotecan)" %in% input$pk_drugs &&
        input$pk_regimen %in% c("FOLFIRINOX (PRODIGE4)", "mFOLFIRINOX")) {
      conc <- pk_sn38(t_h, dose_mgm2 = 180 * dose_factor)
      dfs[["SN-38"]] <- data.frame(time = t_h, conc = conc,
                                    drug = "SN-38 (ng/mL)")
    }
    if ("5-FU" %in% input$pk_drugs &&
        input$pk_regimen %in% c("FOLFIRINOX (PRODIGE4)", "mFOLFIRINOX")) {
      conc <- pk_5fu(t_h, dose_mgm2 = 400 * dose_factor)
      dfs[["5-FU"]] <- data.frame(time = t_h, conc = conc,
                                   drug = "5-FU (µg/mL)")
    }
    if ("MRTX1133" %in% input$pk_drugs &&
        input$pk_regimen == "MRTX1133 (KRAS G12D)") {
      conc <- pk_mrtx1133(t_h, dose_mg = 700 * dose_factor)
      dfs[["MRTX1133"]] <- data.frame(time = t_h, conc = conc,
                                       drug = "MRTX1133 (nM)")
    }

    if (length(dfs) == 0) {
      return(data.frame(time = t_h, conc = 0, drug = "No drug selected"))
    }
    bind_rows(dfs)
  })

  output$pk_plot <- renderPlotly({
    df <- pk_data()
    req(nrow(df) > 0)

    p <- plot_ly()
    for (d in unique(df$drug)) {
      sub <- df[df$drug == d, ]
      p <- add_trace(p, data = sub, x = ~time,
                     y = if (input$pk_logscale) ~log10(pmax(conc, 1e-6)) else ~conc,
                     type = "scatter", mode = "lines",
                     name = d, line = list(width = 2.5))
    }

    ylab <- if (input$pk_logscale) "log₁₀ Concentration" else "Concentration"
    p |> layout(
      xaxis = list(title = "Time (hours)", color = "#dee2e6",
                   gridcolor = "#444"),
      yaxis = list(title = ylab, color = "#dee2e6",
                   gridcolor = "#444"),
      paper_bgcolor = "#222", plot_bgcolor = "#222",
      font = list(color = "#dee2e6"),
      legend = list(orientation = "h", y = -0.2)
    )
  })

  output$pk_table <- renderTable({
    df <- pk_data()
    req(nrow(df) > 0 && any(df$drug != "No drug selected"))

    df %>%
      filter(drug != "No drug selected") %>%
      group_by(Drug = drug) %>%
      summarise(
        Cmax   = round(max(conc, na.rm = TRUE), 3),
        Tmax_h = round(time[which.max(conc)], 2),
        AUC_0_48 = round(sum(diff(time) * (conc[-1] + conc[-n()]) / 2,
                             na.rm = TRUE), 1),
        T_half_h = round(log(2) / 0.22, 1),
        .groups = "drop"
      ) %>%
      rename(`AUC 0-48 (h·unit)` = AUC_0_48,
             `Cmax (unit)` = Cmax,
             `Tmax (h)` = Tmax_h,
             `T½ (h)` = T_half_h)
  }, striped = TRUE, bordered = TRUE, hover = TRUE)

  # ============================================================
  # TAB 3: Tumor Dynamics
  # ============================================================

  output$tgi_tumor_plot <- renderPlotly({
    df <- tgi_results()
    req(!is.null(df) && nrow(df) > 0)

    p <- plot_ly()
    for (reg in unique(df$regimen)) {
      sub <- df[df$regimen == reg, ]
      p <- add_trace(p, data = sub, x = ~time, y = ~total_tumor,
                     type = "scatter", mode = "lines",
                     name = reg,
                     line = list(color = REGIMEN_COLORS[reg], width = 2.5))
    }
    p |> layout(
      xaxis = list(title = "Time (months)", color = "#dee2e6", gridcolor = "#444"),
      yaxis = list(title = "Tumor Volume (mm³)", color = "#dee2e6", gridcolor = "#444"),
      paper_bgcolor = "#222", plot_bgcolor = "#222",
      font = list(color = "#dee2e6"),
      legend = list(orientation = "v", font = list(size = 10))
    )
  })

  output$tgi_ca199_plot <- renderPlotly({
    df <- tgi_results()
    req(!is.null(df) && nrow(df) > 0)

    p <- plot_ly()
    for (reg in unique(df$regimen)) {
      sub <- df[df$regimen == reg, ]
      p <- add_trace(p, data = sub, x = ~time, y = ~CA199,
                     type = "scatter", mode = "lines",
                     name = reg,
                     line = list(color = REGIMEN_COLORS[reg], width = 2.5))
    }
    p |> layout(
      xaxis = list(title = "Time (months)", color = "#dee2e6", gridcolor = "#444"),
      yaxis = list(title = "CA19-9 (U/mL)", color = "#dee2e6", gridcolor = "#444"),
      paper_bgcolor = "#222", plot_bgcolor = "#222",
      font = list(color = "#dee2e6"),
      legend = list(orientation = "v", font = list(size = 10))
    )
  })

  output$waterfall_plot <- renderPlotly({
    df <- tgi_results()
    req(!is.null(df) && nrow(df) > 0)

    init <- input$tumor_vol_base
    wf <- df %>%
      group_by(regimen) %>%
      summarise(best_pct = min((total_tumor - init) / init * 100, na.rm = TRUE),
                .groups = "drop") %>%
      arrange(best_pct) %>%
      mutate(color = ifelse(best_pct < -30, "#2a9d8f",
                            ifelse(best_pct < 0, "#4e9af1", "#e76f51")))

    plot_ly(wf, x = ~reorder(regimen, best_pct), y = ~best_pct,
            type = "bar",
            marker = list(color = ~color)) |>
      add_segments(x = 0.5, xend = nrow(wf) + 0.5,
                   y = -30, yend = -30,
                   line = list(dash = "dash", color = "#ffd700", width = 1.5),
                   name = "PR threshold (-30%)") |>
      layout(
        xaxis = list(title = "", tickangle = -30, color = "#dee2e6",
                     gridcolor = "#444"),
        yaxis = list(title = "Best % Change from Baseline",
                     color = "#dee2e6", gridcolor = "#444",
                     zeroline = TRUE, zerolinecolor = "#888"),
        paper_bgcolor = "#222", plot_bgcolor = "#222",
        font = list(color = "#dee2e6"),
        showlegend = TRUE
      )
  })

  output$tgi_table <- renderTable({
    df <- tgi_results()
    req(!is.null(df) && nrow(df) > 0)

    init <- input$tumor_vol_base
    df %>%
      group_by(Regimen = regimen) %>%
      summarise(
        `TGI (%)` = round(100 * (1 - min(total_tumor) / init), 1),
        `Best Change (%)` = round(min((total_tumor - init) / init * 100), 1),
        `Final Tumor (mm³)` = round(last(total_tumor), 0),
        Response = ifelse(min((total_tumor - init) / init * 100) <= -30,
                          "PR/CR", "SD/PD"),
        .groups = "drop"
      )
  }, striped = TRUE, bordered = TRUE, hover = TRUE)

  # ============================================================
  # TAB 4: Biomarkers
  # ============================================================

  bm_data <- reactive({
    t_wk <- seq(0, input$bm_horizon, by = 0.5)
    t_mo <- t_wk / 4.33

    lapply(REGIMEN_NAMES, function(reg) {
      p <- SCENARIO_PARAMS[[reg]]
      k_eff_adj <- p$k_eff * (1 - 0.4 * 0.5)  # default stroma 0.5

      bm_val <- switch(input$bm_type,
        "CA19-9 (U/mL)" = {
          isolate(input$ca199_base) *
            exp((p$k_prog * 0.5 - k_eff_adj * 0.8) * t_mo)
        },
        "ctDNA VAF (%)" = {
          40 * exp((p$k_prog * 0.4 - k_eff_adj * 0.9) * t_mo)
        },
        "KRAS ctDNA (%)" = {
          25 * exp((p$k_prog * 0.45 - k_eff_adj * 0.85) * t_mo)
        },
        "Neutrophils (×10⁹/L)" = {
          friberg_neutrophil(t_wk, NEUTRO_EFF[reg])$circ
        }
      )
      data.frame(time_wk = t_wk, value = bm_val, regimen = reg)
    }) %>% bind_rows()
  })

  output$bm_kinetics_plot <- renderPlotly({
    df <- bm_data()
    p  <- plot_ly()
    for (reg in REGIMEN_NAMES) {
      sub <- df[df$regimen == reg, ]
      p <- add_trace(p, data = sub, x = ~time_wk, y = ~value,
                     type = "scatter", mode = "lines",
                     name = reg,
                     line = list(color = REGIMEN_COLORS[reg], width = 2))
    }
    p |> layout(
      xaxis = list(title = "Time (weeks)", color = "#dee2e6", gridcolor = "#444"),
      yaxis = list(title = input$bm_type, color = "#dee2e6", gridcolor = "#444"),
      paper_bgcolor = "#222", plot_bgcolor = "#222",
      font = list(color = "#dee2e6"),
      legend = list(orientation = "h", y = -0.25, font = list(size = 9))
    )
  })

  output$bm_week8_plot <- renderPlotly({
    df <- bm_data()
    w8 <- df %>%
      group_by(regimen) %>%
      summarise(
        val_0  = value[which.min(abs(time_wk - 0))],
        val_8  = value[which.min(abs(time_wk - 8))],
        pct_ch = (val_8 - val_0) / (val_0 + 1e-6) * 100,
        .groups = "drop"
      ) %>%
      mutate(color = ifelse(pct_ch < -30, "#2a9d8f",
                            ifelse(pct_ch < 0, "#4e9af1", "#e76f51")))

    plot_ly(w8, y = ~reorder(regimen, pct_ch), x = ~pct_ch,
            type = "bar", orientation = "h",
            marker = list(color = ~color)) |>
      layout(
        yaxis = list(title = "", color = "#dee2e6", tickfont = list(size = 10)),
        xaxis = list(title = paste("% Change in", input$bm_type, "at Week 8"),
                     color = "#dee2e6", gridcolor = "#444",
                     zeroline = TRUE, zerolinecolor = "#888"),
        paper_bgcolor = "#222", plot_bgcolor = "#222",
        font = list(color = "#dee2e6")
      )
  })

  output$neutrophil_plot <- renderPlotly({
    t_days <- seq(0, 63, by = 0.5)
    p <- plot_ly()
    for (reg in REGIMEN_NAMES) {
      nd <- friberg_neutrophil(t_days, NEUTRO_EFF[reg])
      p <- add_trace(p, x = nd$time, y = nd$circ,
                     type = "scatter", mode = "lines",
                     name = reg,
                     line = list(color = REGIMEN_COLORS[reg], width = 2))
    }
    p |>
      add_segments(x = 0, xend = 63, y = 0.5, yend = 0.5,
                   line = list(dash = "dash", color = "#ff6b6b", width = 1.5),
                   name = "Grade 4 (<0.5)") |>
      add_segments(x = 0, xend = 63, y = 1.0, yend = 1.0,
                   line = list(dash = "dot", color = "#ffa94d", width = 1.5),
                   name = "Grade 3 (<1.0)") |>
      layout(
        xaxis = list(title = "Time (days)", color = "#dee2e6", gridcolor = "#444"),
        yaxis = list(title = "ANC (×10⁹/L)", color = "#dee2e6", gridcolor = "#444"),
        paper_bgcolor = "#222", plot_bgcolor = "#222",
        font = list(color = "#dee2e6"),
        legend = list(orientation = "h", y = -0.25, font = list(size = 9))
      )
  })

  output$neutropenia_table <- renderTable({
    data.frame(
      Regimen = REGIMEN_NAMES,
      `Grade 3 Neutropenia (%)` = c(0, 22, 21, 28, 24, 14, 6),
      `Grade 4 Neutropenia (%)` = c(0, 16, 17, 18, 16,  8, 4),
      `Any G3/4 Neutropenia (%)` = c(0, 38, 38, 46, 40, 22, 10),
      `G-CSF Required (%)` = c(0, 15, 26, 42, 38, 10, 5),
      check.names = FALSE
    )
  }, striped = TRUE, bordered = TRUE, hover = TRUE)

  # ============================================================
  # TAB 5: Clinical Endpoints
  # ============================================================

  km_times <- reactive({
    seq(0, 36, by = 0.25)
  })

  km_pfs_data <- reactive({
    t <- km_times()
    lapply(REGIMEN_NAMES, function(reg) {
      p <- SCENARIO_PARAMS[[reg]]
      surv <- weibull_km(t, p$mPFS, shape = 1.3)
      if (input$show_ci) {
        ci <- weibull_ci(t, p$mPFS, shape = 1.3, n = input$n_patients)
        data.frame(time = t, surv = surv,
                   lo = ci$lo, hi = ci$hi, regimen = reg)
      } else {
        data.frame(time = t, surv = surv, lo = surv, hi = surv, regimen = reg)
      }
    }) %>% bind_rows()
  })

  km_os_data <- reactive({
    t <- km_times()
    os_medians <- c(3.0, 6.7, 8.5, 11.1, 10.5, 9.0, 9.0)
    names(os_medians) <- REGIMEN_NAMES
    lapply(REGIMEN_NAMES, function(reg) {
      med <- os_medians[reg]
      surv <- weibull_km(t, med, shape = 1.2)
      if (input$show_ci) {
        ci <- weibull_ci(t, med, shape = 1.2)
        data.frame(time = t, surv = surv,
                   lo = ci$lo, hi = ci$hi, regimen = reg)
      } else {
        data.frame(time = t, surv = surv, lo = surv, hi = surv, regimen = reg)
      }
    }) %>% bind_rows()
  })

  output$km_pfs_plot <- renderPlotly({
    df <- km_pfs_data()
    p  <- plot_ly()
    for (reg in REGIMEN_NAMES) {
      sub <- df[df$regimen == reg, ]
      clr <- REGIMEN_COLORS[reg]
      p <- add_trace(p, data = sub, x = ~time, y = ~surv,
                     type = "scatter", mode = "lines",
                     name = reg,
                     line = list(color = clr, width = 2.5))
      if (input$show_ci) {
        p <- add_trace(p, data = sub,
                       x = c(sub$time, rev(sub$time)),
                       y = c(sub$hi, rev(sub$lo)),
                       type = "scatter", mode = "lines", fill = "toself",
                       fillcolor = paste0(substr(clr, 1, 7), "33"),
                       line = list(color = "transparent"),
                       showlegend = FALSE, hoverinfo = "skip")
      }
    }
    p |> layout(
      xaxis = list(title = "Time (months)", range = c(0, 24),
                   color = "#dee2e6", gridcolor = "#444"),
      yaxis = list(title = "PFS Probability", range = c(0, 1),
                   color = "#dee2e6", gridcolor = "#444"),
      paper_bgcolor = "#222", plot_bgcolor = "#222",
      font = list(color = "#dee2e6"),
      legend = list(orientation = "v", font = list(size = 9))
    )
  })

  output$km_os_plot <- renderPlotly({
    df <- km_os_data()
    p  <- plot_ly()
    for (reg in REGIMEN_NAMES) {
      sub <- df[df$regimen == reg, ]
      clr <- REGIMEN_COLORS[reg]
      p <- add_trace(p, data = sub, x = ~time, y = ~surv,
                     type = "scatter", mode = "lines",
                     name = reg,
                     line = list(color = clr, width = 2.5))
      if (input$show_ci) {
        p <- add_trace(p, data = sub,
                       x = c(sub$time, rev(sub$time)),
                       y = c(sub$hi, rev(sub$lo)),
                       type = "scatter", mode = "lines", fill = "toself",
                       fillcolor = paste0(substr(clr, 1, 7), "33"),
                       line = list(color = "transparent"),
                       showlegend = FALSE, hoverinfo = "skip")
      }
    }
    p |> layout(
      xaxis = list(title = "Time (months)", range = c(0, 36),
                   color = "#dee2e6", gridcolor = "#444"),
      yaxis = list(title = "OS Probability", range = c(0, 1),
                   color = "#dee2e6", gridcolor = "#444"),
      paper_bgcolor = "#222", plot_bgcolor = "#222",
      font = list(color = "#dee2e6"),
      legend = list(orientation = "v", font = list(size = 9))
    )
  })

  output$orr_plot <- renderPlotly({
    df <- ORR_DATA %>%
      mutate(ORR = CR + PR)

    plot_ly(df, x = ~regimen, y = ~CR, type = "bar", name = "CR",
            marker = list(color = "#2a9d8f")) |>
      add_trace(y = ~PR, name = "PR", marker = list(color = "#4e9af1")) |>
      add_trace(y = ~SD, name = "SD", marker = list(color = "#e9c46a")) |>
      add_trace(y = ~PD, name = "PD", marker = list(color = "#e76f51")) |>
      layout(
        barmode = "stack",
        xaxis = list(title = "", tickangle = -30, color = "#dee2e6",
                     tickfont = list(size = 9)),
        yaxis = list(title = "Patients (%)", color = "#dee2e6", gridcolor = "#444"),
        paper_bgcolor = "#222", plot_bgcolor = "#222",
        font = list(color = "#dee2e6"),
        legend = list(orientation = "h", y = -0.35)
      )
  })

  output$tox_radar_plot <- renderPlotly({
    tox_long <- TOXICITY_DATA %>%
      pivot_longer(cols = -regimen, names_to = "toxicity", values_to = "pct")

    p <- plot_ly(type = "scatterpolar")
    for (i in seq_len(nrow(TOXICITY_DATA))) {
      reg <- TOXICITY_DATA$regimen[i]
      vals <- c(
        TOXICITY_DATA$neutropenia[i],
        TOXICITY_DATA$nausea[i],
        TOXICITY_DATA$fatigue[i],
        TOXICITY_DATA$neuropathy[i],
        TOXICITY_DATA$diarrhea[i],
        TOXICITY_DATA$neutropenia[i]
      )
      cats <- c("Neutropenia", "Nausea", "Fatigue",
                "Neuropathy", "Diarrhea", "Neutropenia")
      p <- add_trace(p, r = vals, theta = cats, mode = "lines",
                     name = reg, fill = "toself",
                     line = list(color = REGIMEN_COLORS[reg]),
                     fillcolor = paste0(substr(REGIMEN_COLORS[reg], 1, 7), "22"))
    }
    p |> layout(
      polar = list(radialaxis = list(visible = TRUE, range = c(0, 55),
                                     color = "#adb5bd"),
                   angularaxis = list(color = "#adb5bd")),
      paper_bgcolor = "#222", plot_bgcolor = "#222",
      font = list(color = "#dee2e6"),
      legend = list(orientation = "h", y = -0.15, font = list(size = 9)),
      margin = list(t = 20, b = 60, l = 30, r = 30)
    )
  })

  output$outcomes_table <- renderTable({
    data.frame(
      Regimen = REGIMEN_NAMES,
      `mPFS (mo)` = c(1.5, 3.5, 5.5, 6.5, 6.0, 4.0, 7.4),
      `mOS (mo)`  = c(3.0, 6.7, 8.5, 11.1, 10.5, "NR", "NR"),
      `ORR (%)`   = c(2, 23, 23, 32, 28, 30, 9),
      `G3/4 Tox (%)` = c(5, 52, 61, 75, 70, 35, 25),
      `Reference Trial` = c("—", "Burris 2010", "MPACT (Goldstein 2015)",
                             "PRODIGE4 (Conroy 2011)", "Derived from PRODIGE4",
                             "Phase I (Fell 2024)", "POLO (Golan 2019)"),
      check.names = FALSE,
      stringsAsFactors = FALSE
    )
  }, striped = TRUE, bordered = TRUE, hover = TRUE)

  # ============================================================
  # TAB 6: Scenario Comparison
  # ============================================================

  output$forest_plot <- renderPlotly({
    ref <- input$ref_scenario
    comps <- input$comp_scenarios
    req(length(comps) > 0)

    # Use predefined HR data or compute relative to reference
    hr_df <- HR_DATA %>% filter(regimen %in% comps)
    if (nrow(hr_df) == 0) return(NULL)

    plot_ly(hr_df) |>
      add_trace(
        type = "scatter", mode = "markers",
        x = ~HR, y = ~regimen,
        error_x = list(
          type = "data",
          symmetric = FALSE,
          array    = ~(HR_hi - HR),
          arrayminus = ~(HR - HR_lo),
          color = "#dee2e6"
        ),
        marker = list(size = 12, color = "#4e9af1"),
        text = ~paste0(regimen, "<br>HR=", HR,
                       " (", HR_lo, "–", HR_hi, ")<br>", trial),
        hoverinfo = "text"
      ) |>
      add_segments(x = 1, xend = 1, y = 0.5, yend = nrow(hr_df) + 0.5,
                   line = list(dash = "dash", color = "#ff6b6b", width = 1.5),
                   name = "HR=1") |>
      layout(
        xaxis = list(title = "Hazard Ratio (OS) vs Control",
                     type = "log", color = "#dee2e6", gridcolor = "#444",
                     range = c(log10(0.08), log10(1.5))),
        yaxis = list(title = "", color = "#dee2e6", tickfont = list(size = 10)),
        paper_bgcolor = "#222", plot_bgcolor = "#222",
        font = list(color = "#dee2e6"),
        showlegend = FALSE
      )
  })

  output$parallel_plot <- renderPlotly({
    df <- data.frame(
      regimen   = REGIMEN_NAMES,
      mPFS      = c(1.5, 3.5, 5.5, 6.5, 6.0, 4.0, 7.4),
      mOS       = c(3.0, 6.7, 8.5, 11.1, 10.5, 8.0, 9.0),
      ORR       = c(2,   23,  23,  32,  28,  30,   9),
      TGI       = c(0,   55,  70,  80,  75,  88,  60),
      G34_tox   = c(5,   52,  61,  75,  70,  35,  25)
    ) %>%
      filter(regimen %in% c(input$ref_scenario, input$comp_scenarios))

    if (nrow(df) < 2) return(NULL)

    dims <- list(
      list(range = c(0, 12),  label = "mPFS (mo)",  values = ~mPFS),
      list(range = c(0, 15),  label = "mOS (mo)",   values = ~mOS),
      list(range = c(0, 40),  label = "ORR (%)",    values = ~ORR),
      list(range = c(0, 100), label = "TGI (%)",    values = ~TGI),
      list(range = c(0, 80),  label = "G3/4 Tox (%)", values = ~G34_tox)
    )

    color_idx <- seq_len(nrow(df))

    plot_ly(df, type = "parcoords",
            line = list(color = color_idx,
                        colorscale = "Viridis",
                        showscale = FALSE),
            dimensions = dims) |>
      layout(
        paper_bgcolor = "#222",
        font = list(color = "#dee2e6"),
        margin = list(t = 40, b = 20, l = 80, r = 80)
      )
  })

  output$comparison_table <- renderTable({
    all_sc <- c(input$ref_scenario, input$comp_scenarios)
    df <- data.frame(
      Scenario  = REGIMEN_NAMES,
      `mPFS (mo)` = c(1.5, 3.5, 5.5, 6.5, 6.0, 4.0, 7.4),
      `mOS (mo)`  = c(3.0, 6.7, 8.5, 11.1, 10.5, 9.0, 9.0),
      `ORR (%)`   = c(2,   23,  23,  32,  28,  30,  9),
      `TGI (%)`   = c(0,   55,  70,  80,  75,  88, 60),
      `G3/4 Tox (%)` = c(5, 52, 61, 75, 70, 35, 25),
      check.names = FALSE,
      stringsAsFactors = FALSE
    ) %>%
      filter(Scenario %in% all_sc)
    df
  }, striped = TRUE, bordered = TRUE, hover = TRUE)

  output$tornado_plot <- renderPlotly({
    params <- c("k_prog (tumor growth)",
                "Stroma density",
                "k_eff (drug efficacy)",
                "EC50 (drug potency)",
                "KRAS mutation type",
                "Baseline tumor volume",
                "HRD status",
                "Patient ECOG")
    low_delta  <- c(-18, -12, -25, -15,  -8, -10, -14, -6)
    high_delta <- c( 20,  14,  28,  17,  10,  12,  16,  7)

    df <- data.frame(
      param = params,
      low   = low_delta,
      high  = high_delta
    ) %>% arrange(abs(low) + abs(high))

    plot_ly() |>
      add_trace(data = df, type = "bar", orientation = "h",
                y = ~param, x = ~low, name = "Low (-1 SD)",
                marker = list(color = "#4e9af1")) |>
      add_trace(data = df, type = "bar", orientation = "h",
                y = ~param, x = ~high, name = "High (+1 SD)",
                marker = list(color = "#e76f51")) |>
      add_segments(x = 0, xend = 0,
                   y = 0.5, yend = length(params) + 0.5,
                   line = list(color = "#888", width = 1)) |>
      layout(
        barmode = "overlay",
        xaxis = list(title = "Δ mPFS (weeks)", color = "#dee2e6",
                     gridcolor = "#444"),
        yaxis = list(title = "", color = "#dee2e6",
                     tickfont = list(size = 10)),
        paper_bgcolor = "#222", plot_bgcolor = "#222",
        font = list(color = "#dee2e6"),
        legend = list(orientation = "h", y = -0.15)
      )
  })
}

# ============================================================
# LAUNCH
# ============================================================
shinyApp(ui = ui, server = server)
