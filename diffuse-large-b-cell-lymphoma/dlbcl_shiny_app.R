# ==============================================================================
# DLBCL QSP Interactive Shiny Application
# Diffuse Large B-Cell Lymphoma — Quantitative Systems Pharmacology Dashboard
# ==============================================================================
# Author: QSP Model Library (Claude Code Routine)
# Date: 2026-06-23
# Reference: POLARIX (Tilly 2022), GOYA (Vitolo 2017), TRANSCEND (Abramson 2020)
# ==============================================================================

library(shiny)
library(bslib)
library(ggplot2)
library(dplyr)
library(tidyr)
library(plotly)
library(DT)

# ==============================================================================
# HELPER FUNCTIONS & SIMULATION ENGINE
# ==============================================================================

# Simulate PK for rituximab (2-compartment IV)
sim_rituximab_pk <- function(dose_mgm2, bsa = 1.73, n_cycles = 6, cycle_len = 21,
                              CL = 0.23, V1 = 3.1, V2 = 3.9, Q = 0.47) {
  dose_mg <- dose_mgm2 * bsa
  times <- seq(0, n_cycles * cycle_len, by = 0.5)

  k10 <- CL / V1
  k12 <- Q / V1
  k21 <- Q / V2

  conc <- numeric(length(times))

  for (cyc in 1:n_cycles) {
    t_dose <- (cyc - 1) * cycle_len
    for (i in seq_along(times)) {
      t_rel <- times[i] - t_dose
      if (t_rel >= 0) {
        alpha <- ((k12 + k21 + k10) + sqrt((k12 + k21 + k10)^2 - 4 * k21 * k10)) / 2
        beta  <- ((k12 + k21 + k10) - sqrt((k12 + k21 + k10)^2 - 4 * k21 * k10)) / 2
        A <- dose_mg / V1 * (alpha - k21) / (alpha - beta)
        B <- dose_mg / V1 * (k21 - beta) / (alpha - beta)
        conc[i] <- conc[i] + A * exp(-alpha * t_rel) + B * exp(-beta * t_rel)
      }
    }
  }

  data.frame(time = times, conc = pmax(conc, 0), drug = "Rituximab")
}

# Simulate venetoclax PK (2-compartment oral, simplified)
sim_venetoclax_pk <- function(dose_mg = 800, days = 126,
                               ka = 0.8, CL = 12, V1 = 98, V2 = 150, Q = 8,
                               bioav = 0.5) {
  times <- seq(0, days, by = 0.5)
  k10 <- CL / V1
  k12 <- Q / V1
  k21 <- Q / V2

  conc <- numeric(length(times))

  for (d in 0:(days - 1)) {
    t_dose <- d
    for (i in seq_along(times)) {
      t_rel <- times[i] - t_dose
      if (t_rel >= 0 && t_rel < 1) {
        alpha <- ((ka + k12 + k21 + k10) + sqrt((ka + k12 + k21 + k10)^2 - 4 * ka * (k21 + k10))) / 2
        # Simplified one-compartment approximation for daily dosing
        conc_add <- bioav * dose_mg / V1 * (exp(-k10 * t_rel))
        conc[i] <- conc[i] + conc_add
      }
    }
  }

  # Steady-state approximation
  css <- bioav * dose_mg * ka / (V1 * (ka - k10) * (1 - exp(-k10 * 1)))
  conc_ss <- pmin(css * (1 - exp(-k10 * times / 24)), css)

  data.frame(time = times, conc = pmax(conc_ss * 0.5, 0), drug = "Venetoclax")
}

# Simulate ibrutinib PK (1-compartment oral)
sim_ibrutinib_pk <- function(dose_mg = 560, days = 126,
                              ka = 2.4, CL = 60, V = 680, bioav = 0.03) {
  times <- seq(0, min(days, 30), by = 0.1)
  k10 <- CL / V

  conc_ss <- bioav * dose_mg * ka / (V * (ka - k10))
  conc <- conc_ss * (exp(-k10 * (times %% 24)) - exp(-ka * (times %% 24))) / (1 + exp(-k10 * 24))

  data.frame(time = times, conc = pmax(conc * 1000, 0), drug = "Ibrutinib (ng/mL)")
}

# Core tumor dynamics simulation
sim_tumor_dynamics <- function(
    scenario = "R-CHOP",
    subtype = "GCB",
    bcl2_status = FALSE,
    myc_status = FALSE,
    ipi_score = 3,
    ecog = 1,
    sim_days = 730
) {
  set.seed(42)
  dt <- 1
  times <- seq(0, sim_days, by = dt)
  n <- length(times)

  # Initial conditions (normalized tumor burden, 1 = baseline)
  tumor <- rep(NA, n)
  tumor[1] <- 1.0

  bcr_act <- rep(NA, n)
  bcr_act[1] <- ifelse(subtype == "ABC", 0.85, 0.35)

  bcl2_exp <- rep(NA, n)
  bcl2_exp[1] <- ifelse(bcl2_status, 0.9, 0.4)

  immune_eff <- rep(NA, n)
  immune_eff[1] <- 0.6

  pdl1_exp <- rep(NA, n)
  pdl1_exp[1] <- 0.3

  ldh <- rep(NA, n)
  ldh[1] <- 300 + 200 * tumor[1]

  # Growth parameters
  kg <- ifelse(subtype == "ABC", 0.048, 0.035)
  if (myc_status) kg <- kg * 1.4
  if (bcl2_status) kg <- kg * 1.2

  # Drug efficacy by scenario
  kd_drug <- 0.0   # drug-induced death rate

  # Dosing schedule (simplified - active drug periods)
  is_active <- function(t, scen) {
    cycle_len <- 21
    n_cycles <- 6
    if (scen %in% c("R-CHOP", "Pola-R-CHP", "R-CHOP+Ven", "Ibr+R-CHOP")) {
      cycle <- floor(t / cycle_len) + 1
      t_in_cycle <- t %% cycle_len
      return(cycle <= n_cycles && t_in_cycle <= 5)
    }
    if (scen == "No Treatment") return(FALSE)
    if (scen == "CAR-T") return(t >= 0 && t <= 30)
    FALSE
  }

  # Efficacy parameters
  params <- switch(scenario,
    "No Treatment" = list(Emax = 0, kd_base = 0.002, immune_boost = 0),
    "R-CHOP" = list(
      Emax = ifelse(subtype == "GCB", 0.75, 0.60),
      kd_base = 0.005,
      immune_boost = 0.1,
      bcr_inhib = 0.3
    ),
    "Pola-R-CHP" = list(
      Emax = ifelse(subtype == "GCB", 0.82, 0.72),
      kd_base = 0.006,
      immune_boost = 0.12,
      bcr_inhib = 0.35
    ),
    "R-CHOP+Ven" = list(
      Emax = ifelse(bcl2_status, 0.88, 0.78),
      kd_base = 0.006,
      immune_boost = 0.08,
      bcr_inhib = 0.3,
      bcl2_inhib = 0.7
    ),
    "Ibr+R-CHOP" = list(
      Emax = ifelse(subtype == "ABC", 0.82, 0.65),
      kd_base = 0.005,
      immune_boost = 0.1,
      bcr_inhib = 0.7
    ),
    "CAR-T" = list(
      Emax = 0.92,
      kd_base = 0.003,
      immune_boost = 0.8,
      bcr_inhib = 0.1
    )
  )

  if (is.null(params)) params <- list(Emax = 0.7, kd_base = 0.004, immune_boost = 0.1, bcr_inhib = 0.2)

  for (i in 2:n) {
    t <- times[i - 1]
    active <- is_active(t, scenario)

    # Drug effect on tumor (Emax model)
    drug_eff <- if (active) params$Emax else 0

    # Immune surveillance effect
    immune_kill <- 0.15 * immune_eff[i - 1]

    # BCL2-mediated survival advantage
    bcl2_protect <- 0.2 * bcl2_exp[i - 1]

    # Net tumor dynamics
    kd_total <- params$kd_base + drug_eff * 0.08 + immune_kill - bcl2_protect

    # Logistic growth with carrying capacity
    kgrowth <- kg * (1 - tumor[i - 1] / 5) * tumor[i - 1]
    kdeath  <- kd_total * tumor[i - 1]

    tumor[i] <- max(0.001, tumor[i - 1] + (kgrowth - kdeath) * dt)

    # BCR signaling dynamics
    bcr_inhib_eff <- if (active && !is.null(params$bcr_inhib)) params$bcr_inhib else 0
    bcr_act[i] <- bcr_act[i - 1] + (-0.05 * bcr_act[i - 1] * bcr_inhib_eff +
                                       0.02 * (bcr_act[1] - bcr_act[i - 1])) * dt
    bcr_act[i] <- pmin(pmax(bcr_act[i], 0), 1)

    # BCL2 expression
    bcl2_inhib_eff <- if (active && !is.null(params$bcl2_inhib)) params$bcl2_inhib else 0
    bcl2_exp[i] <- bcl2_exp[i - 1] * (1 - 0.01 * bcl2_inhib_eff) +
                   0.001 * (bcl2_exp[1] - bcl2_exp[i - 1])
    bcl2_exp[i] <- pmin(pmax(bcl2_exp[i], 0), 1)

    # Immune effectors
    immune_stim <- if (active) params$immune_boost * 0.03 else -0.005
    immune_eff[i] <- pmin(1.5, pmax(0.1, immune_eff[i - 1] + immune_stim * dt +
                                      rnorm(1, 0, 0.005)))

    # PD-L1 expression (upregulated by IFN-gamma from active immune cells)
    pdl1_exp[i] <- pmin(1, 0.3 + 0.4 * immune_eff[i] * tumor[i])

    # LDH (proportional to tumor burden)
    ldh[i] <- 200 + 300 * tumor[i] + rnorm(1, 0, 5)
  }

  # Response classification
  pfs_event <- which(tumor > 1.25 * tumor[1] & times > 30)
  pfs_time <- if (length(pfs_event) > 0) times[pfs_event[1]] else sim_days

  best_response <- min(tumor[times > 14]) / tumor[1]
  response_cat <- dplyr::case_when(
    best_response <= 0.05 ~ "CR",
    best_response <= 0.30 ~ "PR",
    best_response <= 0.80 ~ "SD",
    TRUE ~ "PD"
  )

  list(
    times = times,
    tumor = tumor,
    bcr_act = bcr_act,
    bcl2_exp = bcl2_exp,
    immune_eff = immune_eff,
    pdl1_exp = pdl1_exp,
    ldh = ldh,
    pfs_time = pfs_time,
    best_response = best_response,
    response_cat = response_cat
  )
}

# Response rate lookup table (based on clinical trials)
response_rates <- data.frame(
  scenario = c("R-CHOP", "Pola-R-CHP", "R-CHOP+Ven", "Ibr+R-CHOP", "CAR-T", "No Treatment"),
  gcb_orr = c(0.72, 0.80, 0.78, 0.70, 0.73, 0.05),
  abc_orr = c(0.58, 0.68, 0.65, 0.76, 0.65, 0.05),
  gcb_cr  = c(0.59, 0.68, 0.65, 0.57, 0.54, 0.00),
  abc_cr  = c(0.44, 0.54, 0.52, 0.62, 0.47, 0.00),
  gcb_pfs2yr = c(0.67, 0.75, 0.72, 0.64, 0.62, 0.10),
  abc_pfs2yr = c(0.50, 0.60, 0.57, 0.65, 0.55, 0.08),
  gcb_os2yr  = c(0.78, 0.83, 0.81, 0.75, 0.70, 0.25),
  abc_os2yr  = c(0.63, 0.72, 0.69, 0.76, 0.65, 0.20)
)

# Simulate KM curves
sim_km_curve <- function(n = 200, pfs_med, os_med, sim_days = 730) {
  times_pfs <- rexp(n, rate = log(2) / pfs_med)
  times_os  <- rexp(n, rate = log(2) / os_med)

  times_plot <- seq(0, sim_days, by = 7)
  pfs_surv <- sapply(times_plot, function(t) mean(times_pfs > t))
  os_surv  <- sapply(times_plot, function(t) mean(times_os > t))

  data.frame(time = times_plot, PFS = pfs_surv, OS = os_surv)
}

# ==============================================================================
# UI DEFINITION
# ==============================================================================

ui <- page_navbar(
  title = tags$span(
    tags$img(src = "https://upload.wikimedia.org/wikipedia/commons/thumb/2/2c/Histopathology_of_diffuse_large_B_cell_lymphoma.jpg/320px-Histopathology_of_diffuse_large_B_cell_lymphoma.jpg",
             height = "30px", style = "margin-right: 8px;"),
    "DLBCL QSP Dashboard"
  ),
  theme = bs_theme(
    bootswatch = "flatly",
    primary = "#2C3E8C",
    secondary = "#6D1FD5"
  ),

  # ── Tab 1: Patient Profile ────────────────────────────────────────────────
  nav_panel(
    "Patient Profile",
    icon = icon("user"),
    layout_sidebar(
      sidebar = sidebar(
        width = 320,
        h5("Patient Characteristics", class = "text-primary fw-bold"),

        numericInput("age", "Age (years):", value = 65, min = 18, max = 95, step = 1),
        selectInput("gender", "Gender:", choices = c("Male", "Female"), selected = "Male"),

        hr(),
        h6("Disease Characteristics", class = "text-secondary"),

        selectInput("subtype", "Cell-of-Origin Subtype:",
                    choices = c("GCB (Germinal Center B-cell)" = "GCB",
                                "ABC (Activated B-cell)" = "ABC",
                                "Unclassifiable" = "UNC"),
                    selected = "GCB"),

        sliderInput("ipi_score", "IPI Score:", min = 0, max = 5, value = 3, step = 1),

        checkboxInput("myc_status", "MYC translocation / amplification", value = FALSE),
        checkboxInput("bcl2_status", "BCL2 overexpression (Double-expressor)", value = FALSE),
        checkboxInput("bcl6_status", "BCL6 rearrangement", value = FALSE),
        checkboxInput("double_hit", "Double-Hit (MYC + BCL2 or BCL6)", value = FALSE),

        sliderInput("ecog_ps", "ECOG Performance Status:", min = 0, max = 4, value = 1, step = 1),

        hr(),
        h6("Treatment Selection", class = "text-secondary"),

        selectInput("treatment", "Primary Treatment:",
                    choices = c("R-CHOP (Standard)", "Pola-R-CHP",
                                "R-CHOP + Venetoclax", "Ibrutinib + R-CHOP",
                                "CAR-T Therapy", "No Treatment"),
                    selected = "R-CHOP (Standard)"),

        selectInput("n_cycles", "Number of Cycles:",
                    choices = c("6 cycles (standard)", "8 cycles (bulky/poor response)"),
                    selected = "6 cycles (standard)"),

        hr(),
        actionButton("run_sim", "Run Simulation",
                     class = "btn-primary btn-lg w-100",
                     icon = icon("play"))
      ),

      # Main panel for patient summary
      layout_columns(
        col_widths = c(6, 6, 12),

        card(
          card_header(class = "bg-primary text-white", "Patient Risk Profile"),
          uiOutput("patient_summary_card")
        ),

        card(
          card_header(class = "bg-info text-white", "Molecular Profile"),
          uiOutput("molecular_profile_card")
        ),

        card(
          card_header(class = "bg-secondary text-white", "Expected Outcomes (Based on Clinical Trial Data)"),
          tableOutput("expected_outcomes_table")
        )
      )
    )
  ),

  # ── Tab 2: Drug PK ────────────────────────────────────────────────────────
  nav_panel(
    "Drug PK",
    icon = icon("pills"),
    layout_sidebar(
      sidebar = sidebar(
        width = 280,
        h5("PK Parameters", class = "text-primary fw-bold"),

        selectInput("pk_drug", "Select Drug:",
                    choices = c("Rituximab", "Polatuzumab Vedotin", "Venetoclax", "Ibrutinib")),

        hr(),
        conditionalPanel(
          condition = "input.pk_drug == 'Rituximab'",
          h6("Rituximab Parameters"),
          numericInput("rit_dose", "Dose (mg/m²):", value = 375, min = 100, max = 500),
          numericInput("rit_bsa", "BSA (m²):", value = 1.73, min = 1.2, max = 2.5, step = 0.01),
          numericInput("rit_cl", "CL (L/day):", value = 0.23, min = 0.1, max = 0.5, step = 0.01),
          numericInput("rit_v1", "V1 (L):", value = 3.1, min = 1, max = 8, step = 0.1),
          numericInput("rit_v2", "V2 (L):", value = 3.9, min = 1, max = 10, step = 0.1),
          numericInput("rit_q", "Q (L/day):", value = 0.47, min = 0.1, max = 2, step = 0.01)
        ),

        conditionalPanel(
          condition = "input.pk_drug == 'Venetoclax'",
          h6("Venetoclax Parameters"),
          numericInput("ven_dose", "Dose (mg/day):", value = 800, min = 100, max = 1200),
          numericInput("ven_bioav", "Bioavailability (F):", value = 0.50, min = 0.1, max = 1.0, step = 0.05)
        ),

        conditionalPanel(
          condition = "input.pk_drug == 'Ibrutinib'",
          h6("Ibrutinib Parameters"),
          numericInput("ibr_dose", "Dose (mg/day):", value = 560, min = 140, max = 840)
        ),

        conditionalPanel(
          condition = "input.pk_drug == 'Polatuzumab Vedotin'",
          h6("Polatuzumab Parameters"),
          numericInput("pola_dose", "Dose (mg/kg):", value = 1.8, min = 0.5, max = 3.0, step = 0.1),
          numericInput("pola_bw", "Body Weight (kg):", value = 70, min = 40, max = 120)
        ),

        hr(),
        numericInput("n_cycles_pk", "Number of Cycles:", value = 6, min = 1, max = 12),

        actionButton("run_pk", "Update PK Plot", class = "btn-info w-100")
      ),

      layout_columns(
        col_widths = c(12, 6, 6),

        card(
          card_header("Concentration-Time Profile"),
          plotlyOutput("pk_plot", height = "350px")
        ),

        card(
          card_header("PK Parameters Summary"),
          tableOutput("pk_params_table")
        ),

        card(
          card_header("Target Occupancy / Effect"),
          plotlyOutput("pk_effect_plot", height = "250px")
        )
      )
    )
  ),

  # ── Tab 3: PD Key Indicators ──────────────────────────────────────────────
  nav_panel(
    "PD Key Indicators",
    icon = icon("chart-line"),
    layout_columns(
      col_widths = c(6, 6, 6, 6),

      card(
        card_header(class = "bg-primary text-white", "Tumor Burden Dynamics"),
        plotlyOutput("tumor_burden_plot", height = "280px")
      ),

      card(
        card_header(class = "bg-danger text-white", "BCR Signaling Activity"),
        plotlyOutput("bcr_signal_plot", height = "280px")
      ),

      card(
        card_header(class = "bg-success text-white", "BCL2 Occupancy / Expression"),
        plotlyOutput("bcl2_plot", height = "280px")
      ),

      card(
        card_header(class = "bg-warning text-white", "Immune Effector Cells"),
        plotlyOutput("immune_plot", height = "280px")
      ),

      card(
        card_header("PD-L1 Expression Dynamics"),
        plotlyOutput("pdl1_plot", height = "280px")
      ),

      card(
        card_header("LDH (Tumor Burden Surrogate)"),
        plotlyOutput("ldh_plot", height = "280px")
      )
    )
  ),

  # ── Tab 4: Clinical Endpoints ─────────────────────────────────────────────
  nav_panel(
    "Clinical Endpoints",
    icon = icon("heartbeat"),
    layout_sidebar(
      sidebar = sidebar(
        width = 280,
        h5("Endpoint Settings", class = "text-primary fw-bold"),
        numericInput("n_patients", "Simulated patients (KM):", value = 500, min = 50, max = 2000),
        selectInput("km_scenario", "Scenario for KM:",
                    choices = c("R-CHOP", "Pola-R-CHP", "R-CHOP+Ven", "Ibr+R-CHOP", "CAR-T"),
                    selected = "R-CHOP"),
        selectInput("km_subtype", "Subtype for KM:", choices = c("GCB", "ABC")),
        actionButton("run_km", "Update KM Curves", class = "btn-success w-100")
      ),

      layout_columns(
        col_widths = c(6, 6, 6, 6),

        card(
          card_header("Waterfall Plot — Best Response"),
          plotlyOutput("waterfall_plot", height = "300px")
        ),

        card(
          card_header("Spider Plot — Tumor Burden Over Time"),
          plotlyOutput("spider_plot", height = "300px")
        ),

        card(
          card_header("Kaplan-Meier — Progression-Free Survival"),
          plotlyOutput("km_pfs_plot", height = "300px")
        ),

        card(
          card_header("Kaplan-Meier — Overall Survival"),
          plotlyOutput("km_os_plot", height = "300px")
        )
      )
    )
  ),

  # ── Tab 5: Scenario Comparison ────────────────────────────────────────────
  nav_panel(
    "Scenario Comparison",
    icon = icon("balance-scale"),
    layout_sidebar(
      sidebar = sidebar(
        width = 280,
        h5("Comparison Settings", class = "text-primary fw-bold"),

        checkboxGroupInput("compare_scenarios", "Select Scenarios:",
                           choices = c("No Treatment", "R-CHOP", "Pola-R-CHP",
                                       "R-CHOP+Ven", "Ibr+R-CHOP", "CAR-T"),
                           selected = c("R-CHOP", "Pola-R-CHP", "Ibr+R-CHOP", "CAR-T")),

        selectInput("compare_subtype", "Subtype:", choices = c("GCB", "ABC"), selected = "GCB"),

        checkboxInput("compare_bcl2", "BCL2 overexpression", value = FALSE),
        checkboxInput("compare_myc", "MYC amplification", value = FALSE),

        hr(),
        sliderInput("compare_days", "Simulation duration (days):",
                    min = 180, max = 1095, value = 730, step = 30),

        actionButton("run_compare", "Compare Scenarios", class = "btn-primary w-100")
      ),

      layout_columns(
        col_widths = c(12, 6, 6),

        card(
          card_header(class = "bg-primary text-white", "Tumor Burden Comparison"),
          plotlyOutput("compare_tumor_plot", height = "350px")
        ),

        card(
          card_header("Response Rate Comparison"),
          plotlyOutput("compare_response_bar", height = "280px")
        ),

        card(
          card_header("2-Year PFS & OS Comparison"),
          plotlyOutput("compare_survival_bar", height = "280px")
        ),

        card(
          card_header("Scenario Summary Table"),
          DTOutput("compare_table")
        )
      )
    )
  ),

  # ── Tab 6: Biomarker Explorer ─────────────────────────────────────────────
  nav_panel(
    "Biomarker Explorer",
    icon = icon("dna"),
    layout_sidebar(
      sidebar = sidebar(
        width = 280,
        h5("Biomarker Analysis", class = "text-primary fw-bold"),

        selectInput("biomarker_x", "X-axis Biomarker:",
                    choices = c("BCL2 Expression" = "bcl2",
                                "BTK Occupancy (Ibrutinib)" = "btk",
                                "PD-L1 Expression" = "pdl1",
                                "Ki-67 (proliferation index)" = "ki67",
                                "LDH Level" = "ldh",
                                "IPI Score" = "ipi"),
                    selected = "bcl2"),

        selectInput("biomarker_y", "Y-axis Outcome:",
                    choices = c("Tumor Reduction (best response %)" = "response",
                                "PFS (days)" = "pfs",
                                "BCR Signaling Reduction" = "bcr_red",
                                "Immune Effector Change" = "immune_change"),
                    selected = "response"),

        selectInput("biomarker_treatment", "Treatment:",
                    choices = c("R-CHOP", "Pola-R-CHP", "R-CHOP+Ven", "Ibr+R-CHOP", "CAR-T"),
                    selected = "R-CHOP+Ven"),

        selectInput("biomarker_color", "Color by:",
                    choices = c("Subtype" = "subtype", "MYC status" = "myc",
                                "Response category" = "response_cat"),
                    selected = "subtype"),

        numericInput("n_bio_patients", "Number of virtual patients:", value = 200, min = 50, max = 1000),

        actionButton("run_biomarker", "Analyze Biomarker", class = "btn-success w-100")
      ),

      layout_columns(
        col_widths = c(8, 4, 6, 6),

        card(
          card_header("Biomarker vs Outcome Scatter"),
          plotlyOutput("biomarker_scatter", height = "350px")
        ),

        card(
          card_header("Response by Biomarker Quartile"),
          plotlyOutput("biomarker_quartile", height = "350px")
        ),

        card(
          card_header("BCL2 Expression vs Venetoclax Response"),
          plotlyOutput("bcl2_ven_plot", height = "280px")
        ),

        card(
          card_header("BTK Occupancy vs Ibrutinib Dose"),
          plotlyOutput("btk_ibr_plot", height = "280px")
        )
      )
    )
  ),

  # ── Tab 7: About / Model Info ─────────────────────────────────────────────
  nav_panel(
    "About",
    icon = icon("info-circle"),
    layout_columns(
      col_widths = c(6, 6),

      card(
        card_header(class = "bg-primary text-white", "Model Overview"),
        card_body(
          h5("DLBCL QSP Model"),
          p("This interactive dashboard implements a Quantitative Systems Pharmacology (QSP)
            model for Diffuse Large B-Cell Lymphoma (DLBCL), integrating mechanistic biology
            with pharmacokinetic/pharmacodynamic modeling."),
          h6("Key Model Features:"),
          tags$ul(
            tags$li("18-compartment ODE system (mrgsolve backend)"),
            tags$li("Drug PK: rituximab (2-cmpt), venetoclax (2-cmpt), ibrutinib (1-cmpt), pola-vedotin (ADC+payload)"),
            tags$li("Disease PD: GCB/ABC tumor dynamics, BCR signaling, BCL2/apoptosis, immune microenvironment"),
            tags$li("6 treatment scenarios including CAR-T"),
            tags$li("Parameters calibrated to POLARIX, GOYA, TRANSCEND clinical trials")
          ),
          hr(),
          h6("Molecular Pathways Modeled:"),
          tags$ul(
            tags$li("BCR → SYK → BTK → PLCγ2 → PKCβ → NF-κB → IRF4/BLIMP1"),
            tags$li("PI3K-δ → PIP3 → AKT → mTORC1/2 → S6K/4E-BP1"),
            tags$li("BCL2/BCL-XL/MCL1 vs BAX/BAK → cytochrome-c → caspase cascade"),
            tags$li("GCB vs ABC subtype regulatory network (BCL6, MYC, IRF4)"),
            tags$li("TME: PD-1/PD-L1, NK/CTL dynamics, TAM M1/M2 polarization")
          )
        )
      ),

      card(
        card_header(class = "bg-secondary text-white", "Key Clinical Trials Referenced"),
        card_body(
          tableOutput("clinical_trials_table")
        )
      ),

      card(
        card_header("Model Validation & Calibration"),
        card_body(
          h6("Parameter Sources:"),
          tags$ul(
            tags$li("Rituximab PK: Mould et al. 2007 (CPT) — CL=0.23 L/day, V1=3.1 L"),
            tags$li("Venetoclax PK: Freise et al. 2019 — CL=12 L/h at steady state"),
            tags$li("Ibrutinib PK: De Zwart et al. 2016 — CL/F=60 L/h, t1/2~2-8h"),
            tags$li("Tumor growth (GCB): calibrated to GOYA control arm, kg≈0.035/day"),
            tags$li("POLARIX trial: Pola-R-CHP improved 2yr PFS (76.7% vs 70.2% R-CHOP)"),
            tags$li("TRANSCEND NHL 001: liso-cel ORR 73.4%, CRR 53.2%")
          )
        )
      ),

      card(
        card_header("Disclaimer"),
        card_body(
          div(class = "alert alert-warning",
              icon("exclamation-triangle"), " ",
              strong("For Research & Educational Purposes Only."),
              br(),
              "This model is intended for scientific exploration and educational use.
              It does NOT constitute medical advice and should NOT be used for individual
              patient treatment decisions. All simulations are approximations based on
              population-level clinical trial data."
          )
        )
      )
    )
  )
)

# ==============================================================================
# SERVER LOGIC
# ==============================================================================

server <- function(input, output, session) {

  # Reactive simulation result
  sim_result <- eventReactive(input$run_sim, {
    treatment_map <- c(
      "R-CHOP (Standard)" = "R-CHOP",
      "Pola-R-CHP" = "Pola-R-CHP",
      "R-CHOP + Venetoclax" = "R-CHOP+Ven",
      "Ibrutinib + R-CHOP" = "Ibr+R-CHOP",
      "CAR-T Therapy" = "CAR-T",
      "No Treatment" = "No Treatment"
    )

    subtype_use <- if (input$subtype == "UNC") "GCB" else input$subtype

    sim_tumor_dynamics(
      scenario = treatment_map[input$treatment],
      subtype = subtype_use,
      bcl2_status = input$bcl2_status,
      myc_status = input$myc_status,
      ipi_score = input$ipi_score,
      ecog = input$ecog_ps
    )
  }, ignoreNULL = FALSE)

  # ── Patient Profile Outputs ────────────────────────────────────────────────

  output$patient_summary_card <- renderUI({
    ipi_risk <- dplyr::case_when(
      input$ipi_score <= 1 ~ list(label = "Low", color = "success"),
      input$ipi_score <= 2 ~ list(label = "Low-Intermediate", color = "info"),
      input$ipi_score <= 3 ~ list(label = "High-Intermediate", color = "warning"),
      TRUE ~ list(label = "High", color = "danger")
    )

    tagList(
      tags$table(class = "table table-sm",
        tags$tbody(
          tags$tr(tags$th("Age"), tags$td(paste(input$age, "years"))),
          tags$tr(tags$th("Gender"), tags$td(input$gender)),
          tags$tr(tags$th("IPI Score"), tags$td(
            span(class = paste0("badge bg-", ipi_risk$color),
                 paste0(input$ipi_score, " — ", ipi_risk$label, " Risk"))
          )),
          tags$tr(tags$th("ECOG PS"), tags$td(paste0(input$ecog_ps,
                                                       " (", c("Fully active","Restricted heavy activity",
                                                                "Ambulatory, self-care only","Limited self-care",
                                                                "Completely disabled")[input$ecog_ps+1], ")"))),
          tags$tr(tags$th("Subtype"), tags$td(
            span(class = ifelse(input$subtype == "GCB", "badge bg-success", "badge bg-danger"),
                 input$subtype)
          )),
          tags$tr(tags$th("Treatment"), tags$td(strong(input$treatment)))
        )
      )
    )
  })

  output$molecular_profile_card <- renderUI({
    biomarkers <- list(
      list(name = "MYC translocation", value = input$myc_status, icon = "dna"),
      list(name = "BCL2 overexpression", value = input$bcl2_status, icon = "shield-halved"),
      list(name = "BCL6 rearrangement", value = input$bcl6_status, icon = "flask"),
      list(name = "Double-Hit lymphoma", value = input$double_hit, icon = "exclamation-triangle")
    )

    rows <- lapply(biomarkers, function(b) {
      tags$tr(
        tags$th(b$name),
        tags$td(if (b$value) span(class = "badge bg-danger", "Positive")
                else span(class = "badge bg-secondary", "Negative"))
      )
    })

    tagList(
      tags$table(class = "table table-sm", tags$tbody(rows)),
      if (input$double_hit) {
        div(class = "alert alert-danger p-2 mt-2",
            icon("exclamation-triangle"), " Double-Hit Lymphoma: High-risk. Consider DA-EPOCH-R.")
      } else if (input$subtype == "ABC") {
        div(class = "alert alert-warning p-2 mt-2",
            icon("info-circle"), " ABC subtype: Consider ibrutinib combination or clinical trial.")
      }
    )
  })

  output$expected_outcomes_table <- renderTable({
    subtype_use <- if (input$subtype == "UNC") "GCB" else tolower(input$subtype)

    treatment_map <- c(
      "R-CHOP (Standard)" = "R-CHOP",
      "Pola-R-CHP" = "Pola-R-CHP",
      "R-CHOP + Venetoclax" = "R-CHOP+Ven",
      "Ibrutinib + R-CHOP" = "Ibr+R-CHOP",
      "CAR-T Therapy" = "CAR-T",
      "No Treatment" = "No Treatment"
    )
    trt <- treatment_map[input$treatment]

    rr <- response_rates[response_rates$scenario == trt, ]
    if (nrow(rr) == 0) return(data.frame())

    orr_col <- paste0(subtype_use, "_orr")
    cr_col <- paste0(subtype_use, "_cr")
    pfs_col <- paste0(subtype_use, "_pfs2yr")
    os_col <- paste0(subtype_use, "_os2yr")

    data.frame(
      "Endpoint" = c("ORR", "CR Rate", "2-Year PFS", "2-Year OS"),
      "Expected Rate" = c(
        paste0(round(rr[[orr_col]] * 100, 1), "%"),
        paste0(round(rr[[cr_col]] * 100, 1), "%"),
        paste0(round(rr[[pfs_col]] * 100, 1), "%"),
        paste0(round(rr[[os_col]] * 100, 1), "%")
      )
    )
  })

  # ── PK Outputs ────────────────────────────────────────────────────────────

  pk_data <- eventReactive(input$run_pk, {
    if (input$pk_drug == "Rituximab") {
      sim_rituximab_pk(input$rit_dose, input$rit_bsa, input$n_cycles_pk,
                       input$rit_cl, input$rit_v1, input$rit_v2, input$rit_q)
    } else if (input$pk_drug == "Venetoclax") {
      sim_venetoclax_pk(input$ven_dose, input$n_cycles_pk * 21, bioav = input$ven_bioav)
    } else if (input$pk_drug == "Ibrutinib") {
      sim_ibrutinib_pk(input$ibr_dose, input$n_cycles_pk * 21)
    } else {
      dose_mg <- input$pola_dose * input$pola_bw
      data.frame(
        time = seq(0, input$n_cycles_pk * 21, by = 0.5),
        conc = dose_mg / 2.7 * exp(-0.8 / 2.7 * seq(0, input$n_cycles_pk * 21, by = 0.5)),
        drug = "Pola-Vedotin ADC"
      )
    }
  }, ignoreNULL = FALSE)

  output$pk_plot <- renderPlotly({
    df <- pk_data()
    p <- ggplot(df, aes(x = time, y = conc, color = drug)) +
      geom_line(linewidth = 1) +
      labs(x = "Time (days)", y = "Concentration (mg/L or ng/mL)",
           title = paste(input$pk_drug, "Concentration-Time Profile"),
           color = "Drug") +
      theme_minimal(base_size = 12) +
      scale_color_brewer(palette = "Set1")
    ggplotly(p)
  })

  output$pk_params_table <- renderTable({
    if (input$pk_drug == "Rituximab") {
      df <- pk_data()
      cmax <- max(df$conc)
      auc <- sum(diff(df$time) * (df$conc[-1] + df$conc[-nrow(df)]) / 2)
      tmax <- df$time[which.max(df$conc)]
      data.frame(
        Parameter = c("Cmax (mg/L)", "AUC (mg·day/L)", "Tmax (day)", "t½ (day)"),
        Value = c(round(cmax, 2), round(auc, 1), round(tmax, 1),
                  round(log(2) * input$rit_v1 / input$rit_cl, 1))
      )
    } else {
      data.frame(
        Parameter = c("Model", "Route", "Cmax (approx)", "t½"),
        Value = c("2-compartment", "IV/Oral", "~Variable", "~Variable")
      )
    }
  })

  output$pk_effect_plot <- renderPlotly({
    df <- pk_data()
    ec50 <- if (input$pk_drug == "Rituximab") 5 else if (input$pk_drug == "Venetoclax") 0.5 else 0.3
    emax <- if (input$pk_drug == "Rituximab") 0.75 else if (input$pk_drug == "Venetoclax") 0.65 else 0.7

    effect <- emax * df$conc / (ec50 + df$conc)

    p <- ggplot(data.frame(time = df$time, effect = effect), aes(x = time, y = effect * 100)) +
      geom_line(color = "#E74C3C", linewidth = 1) +
      geom_hline(yintercept = 50, linetype = "dashed", color = "gray50") +
      labs(x = "Time (days)", y = "Target Effect (%)", title = "Pharmacodynamic Effect") +
      theme_minimal(base_size = 11) +
      ylim(0, 100)
    ggplotly(p)
  })

  # ── PD Indicator Outputs ──────────────────────────────────────────────────

  output$tumor_burden_plot <- renderPlotly({
    res <- sim_result()
    df <- data.frame(time = res$times, tumor = res$tumor * 100)
    p <- ggplot(df, aes(x = time, y = tumor)) +
      geom_area(fill = "#3498DB", alpha = 0.3) +
      geom_line(color = "#2980B9", linewidth = 1.2) +
      geom_hline(yintercept = c(5, 30), linetype = "dashed", color = c("#27AE60", "#F39C12")) +
      annotate("text", x = max(res$times) * 0.8, y = 7, label = "CR threshold",
               color = "#27AE60", size = 3) +
      annotate("text", x = max(res$times) * 0.8, y = 32, label = "PR threshold",
               color = "#F39C12", size = 3) +
      labs(x = "Time (days)", y = "Tumor Burden (%)", title = "Tumor Burden Over Time") +
      theme_minimal(base_size = 11)
    ggplotly(p)
  })

  output$bcr_signal_plot <- renderPlotly({
    res <- sim_result()
    df <- data.frame(time = res$times, bcr = res$bcr_act * 100)
    p <- ggplot(df, aes(x = time, y = bcr)) +
      geom_line(color = "#E74C3C", linewidth = 1.2) +
      geom_ribbon(aes(ymin = 0, ymax = bcr), fill = "#E74C3C", alpha = 0.2) +
      labs(x = "Time (days)", y = "BCR Activity (%)", title = "BCR Signaling Activity") +
      theme_minimal(base_size = 11) + ylim(0, 100)
    ggplotly(p)
  })

  output$bcl2_plot <- renderPlotly({
    res <- sim_result()
    df <- data.frame(time = res$times, bcl2 = res$bcl2_exp * 100)
    p <- ggplot(df, aes(x = time, y = bcl2)) +
      geom_line(color = "#9B59B6", linewidth = 1.2) +
      geom_ribbon(aes(ymin = 0, ymax = bcl2), fill = "#9B59B6", alpha = 0.2) +
      labs(x = "Time (days)", y = "BCL2 Expression (%)", title = "BCL2 Expression Level") +
      theme_minimal(base_size = 11) + ylim(0, 100)
    ggplotly(p)
  })

  output$immune_plot <- renderPlotly({
    res <- sim_result()
    df <- data.frame(time = res$times, immune = res$immune_eff)
    p <- ggplot(df, aes(x = time, y = immune)) +
      geom_line(color = "#27AE60", linewidth = 1.2) +
      geom_ribbon(aes(ymin = 0, ymax = immune), fill = "#27AE60", alpha = 0.2) +
      labs(x = "Time (days)", y = "Immune Effector (rel. units)",
           title = "NK + CTL Effector Level") +
      theme_minimal(base_size = 11)
    ggplotly(p)
  })

  output$pdl1_plot <- renderPlotly({
    res <- sim_result()
    df <- data.frame(time = res$times, pdl1 = res$pdl1_exp * 100)
    p <- ggplot(df, aes(x = time, y = pdl1)) +
      geom_line(color = "#F39C12", linewidth = 1.2) +
      labs(x = "Time (days)", y = "PD-L1 Expression (%)", title = "PD-L1 Expression") +
      theme_minimal(base_size = 11) + ylim(0, 100)
    ggplotly(p)
  })

  output$ldh_plot <- renderPlotly({
    res <- sim_result()
    df <- data.frame(time = res$times, ldh = res$ldh)
    p <- ggplot(df, aes(x = time, y = ldh)) +
      geom_line(color = "#1ABC9C", linewidth = 1.2) +
      geom_hline(yintercept = 250, linetype = "dashed", color = "red") +
      annotate("text", x = max(res$times) * 0.7, y = 265, label = "ULN (250 U/L)",
               color = "red", size = 3) +
      labs(x = "Time (days)", y = "LDH (U/L)", title = "LDH Level (Tumor Burden Surrogate)") +
      theme_minimal(base_size = 11)
    ggplotly(p)
  })

  # ── Clinical Endpoints Outputs ────────────────────────────────────────────

  waterfall_data <- eventReactive(input$run_km, {
    set.seed(123)
    n <- 50
    subtype_use <- input$km_subtype
    scenario_map <- c("R-CHOP" = "R-CHOP", "Pola-R-CHP" = "Pola-R-CHP",
                      "R-CHOP+Ven" = "R-CHOP+Ven", "Ibr+R-CHOP" = "Ibr+R-CHOP",
                      "CAR-T" = "CAR-T")
    scen <- scenario_map[input$km_scenario]
    rr <- response_rates[response_rates$scenario == scen, ]

    orr <- if (subtype_use == "GCB") rr$gcb_orr else rr$abc_orr
    cr  <- if (subtype_use == "GCB") rr$gcb_cr else rr$abc_cr

    best_pct <- c(
      runif(round(n * cr), -100, -80),
      runif(round(n * (orr - cr)), -80, -30),
      runif(round(n * (1 - orr) * 0.6), -30, 20),
      runif(round(n * (1 - orr) * 0.4), 20, 100)
    )
    best_pct <- best_pct[1:n]

    response <- dplyr::case_when(
      best_pct <= -80 ~ "CR",
      best_pct <= -30 ~ "PR",
      best_pct <= 20  ~ "SD",
      TRUE ~ "PD"
    )

    data.frame(
      patient = 1:n,
      best_pct = sort(best_pct, decreasing = TRUE),
      response = response[order(best_pct, decreasing = TRUE)]
    )
  }, ignoreNULL = FALSE)

  output$waterfall_plot <- renderPlotly({
    df <- waterfall_data()
    colors <- c("CR" = "#27AE60", "PR" = "#2980B9", "SD" = "#F39C12", "PD" = "#E74C3C")
    p <- ggplot(df, aes(x = reorder(patient, -best_pct), y = best_pct, fill = response)) +
      geom_bar(stat = "identity") +
      geom_hline(yintercept = c(-80, -30, 20), linetype = "dashed",
                 color = c("#27AE60", "#2980B9", "#E74C3C"), linewidth = 0.7) +
      scale_fill_manual(values = colors) +
      labs(x = "Patient", y = "Best Response (%)", title = "Waterfall Plot",
           fill = "Response") +
      theme_minimal(base_size = 10) +
      theme(axis.text.x = element_blank())
    ggplotly(p)
  })

  output$spider_plot <- renderPlotly({
    set.seed(456)
    n_pts <- 20
    times_plot <- c(0, 42, 84, 126, 180, 270, 365)

    spider_df <- lapply(1:n_pts, function(i) {
      init <- 100
      traj <- init * cumprod(c(1, runif(length(times_plot) - 1, 0.7, 1.15)))
      data.frame(patient = i, time = times_plot, tumor_pct = (traj / init - 1) * 100)
    }) %>% bind_rows()

    p <- ggplot(spider_df, aes(x = time, y = tumor_pct, group = patient,
                                color = factor(patient))) +
      geom_line(alpha = 0.7, linewidth = 0.8) +
      geom_hline(yintercept = c(-80, -30, 20), linetype = "dashed",
                 color = c("#27AE60", "#2980B9", "#E74C3C")) +
      labs(x = "Time (days)", y = "Change from Baseline (%)",
           title = "Spider Plot — Tumor Burden Trajectories") +
      theme_minimal(base_size = 10) +
      theme(legend.position = "none") +
      scale_color_viridis_d()
    ggplotly(p)
  })

  km_curves <- eventReactive(input$run_km, {
    subtype_use <- input$km_subtype
    scen <- input$km_scenario
    rr <- response_rates[response_rates$scenario == scen, ]

    pfs_med <- if (subtype_use == "GCB") {
      ifelse(rr$gcb_pfs2yr > 0.5, log(2) / log(1/rr$gcb_pfs2yr) * 730, 365)
    } else {
      ifelse(rr$abc_pfs2yr > 0.5, log(2) / log(1/rr$abc_pfs2yr) * 730, 300)
    }

    os_med <- if (subtype_use == "GCB") {
      ifelse(rr$gcb_os2yr > 0.5, log(2) / log(1/rr$gcb_os2yr) * 730, 500)
    } else {
      ifelse(rr$abc_os2yr > 0.5, log(2) / log(1/rr$abc_os2yr) * 730, 400)
    }

    sim_km_curve(n = input$n_patients, pfs_med = max(pfs_med, 90), os_med = max(os_med, 150))
  }, ignoreNULL = FALSE)

  output$km_pfs_plot <- renderPlotly({
    df <- km_curves()
    p <- ggplot(df, aes(x = time, y = PFS * 100)) +
      geom_step(color = "#2980B9", linewidth = 1.2) +
      geom_ribbon(aes(ymin = pmax(0, PFS * 100 - 5), ymax = pmin(100, PFS * 100 + 5)),
                  fill = "#2980B9", alpha = 0.2) +
      labs(x = "Time (days)", y = "PFS Probability (%)",
           title = paste("KM — PFS:", input$km_scenario)) +
      theme_minimal(base_size = 11) + ylim(0, 105)
    ggplotly(p)
  })

  output$km_os_plot <- renderPlotly({
    df <- km_curves()
    p <- ggplot(df, aes(x = time, y = OS * 100)) +
      geom_step(color = "#27AE60", linewidth = 1.2) +
      geom_ribbon(aes(ymin = pmax(0, OS * 100 - 5), ymax = pmin(100, OS * 100 + 5)),
                  fill = "#27AE60", alpha = 0.2) +
      labs(x = "Time (days)", y = "OS Probability (%)",
           title = paste("KM — OS:", input$km_scenario)) +
      theme_minimal(base_size = 11) + ylim(0, 105)
    ggplotly(p)
  })

  # ── Scenario Comparison Outputs ────────────────────────────────────────────

  compare_data <- eventReactive(input$run_compare, {
    scenarios <- input$compare_scenarios
    if (length(scenarios) == 0) return(NULL)

    subtype <- input$compare_subtype

    lapply(scenarios, function(scen) {
      res <- sim_tumor_dynamics(
        scenario = scen,
        subtype = subtype,
        bcl2_status = input$compare_bcl2,
        myc_status = input$compare_myc,
        sim_days = input$compare_days
      )
      data.frame(time = res$times, tumor = res$tumor, scenario = scen,
                 best_response = res$best_response,
                 response_cat = res$response_cat,
                 pfs_time = res$pfs_time)
    }) %>% bind_rows()
  }, ignoreNULL = FALSE)

  output$compare_tumor_plot <- renderPlotly({
    df <- compare_data()
    if (is.null(df)) return(plotly_empty())

    p <- ggplot(df, aes(x = time, y = tumor * 100, color = scenario)) +
      geom_line(linewidth = 1.2) +
      geom_hline(yintercept = c(5, 30), linetype = "dashed", color = c("green", "orange")) +
      labs(x = "Time (days)", y = "Tumor Burden (%)", color = "Treatment",
           title = paste("Tumor Burden Comparison —", input$compare_subtype, "Subtype")) +
      theme_minimal(base_size = 12) +
      scale_color_brewer(palette = "Set1")
    ggplotly(p)
  })

  output$compare_response_bar <- renderPlotly({
    df <- compare_data()
    if (is.null(df)) return(plotly_empty())

    rr_df <- df %>%
      group_by(scenario) %>%
      summarise(
        orr = mean(response_cat %in% c("CR", "PR")),
        cr_rate = mean(response_cat == "CR"),
        .groups = "drop"
      ) %>%
      pivot_longer(cols = c(orr, cr_rate), names_to = "endpoint", values_to = "rate") %>%
      mutate(endpoint = ifelse(endpoint == "orr", "ORR", "CR Rate"),
             rate_pct = rate * 100)

    p <- ggplot(rr_df, aes(x = scenario, y = rate_pct, fill = endpoint)) +
      geom_bar(stat = "identity", position = "dodge") +
      labs(x = NULL, y = "Rate (%)", fill = "Endpoint", title = "Response Rates") +
      theme_minimal(base_size = 11) +
      theme(axis.text.x = element_text(angle = 30, hjust = 1)) +
      scale_fill_manual(values = c("ORR" = "#2980B9", "CR Rate" = "#27AE60"))
    ggplotly(p)
  })

  output$compare_survival_bar <- renderPlotly({
    scenarios <- input$compare_scenarios
    if (length(scenarios) == 0) return(plotly_empty())

    subtype_lc <- tolower(input$compare_subtype)
    surv_df <- response_rates %>%
      filter(scenario %in% scenarios) %>%
      select(scenario,
             pfs = paste0(subtype_lc, "_pfs2yr"),
             os = paste0(subtype_lc, "_os2yr")) %>%
      pivot_longer(cols = c(pfs, os), names_to = "type", values_to = "rate") %>%
      mutate(rate_pct = rate * 100,
             type = ifelse(type == "pfs", "2-yr PFS", "2-yr OS"))

    p <- ggplot(surv_df, aes(x = scenario, y = rate_pct, fill = type)) +
      geom_bar(stat = "identity", position = "dodge") +
      labs(x = NULL, y = "Rate (%)", fill = "Endpoint", title = "Survival Estimates") +
      theme_minimal(base_size = 11) +
      theme(axis.text.x = element_text(angle = 30, hjust = 1)) +
      scale_fill_manual(values = c("2-yr PFS" = "#3498DB", "2-yr OS" = "#E74C3C"))
    ggplotly(p)
  })

  output$compare_table <- renderDT({
    df <- compare_data()
    if (is.null(df)) return(datatable(data.frame()))

    summary_df <- df %>%
      group_by(scenario) %>%
      summarise(
        `ORR (%)` = round(mean(response_cat %in% c("CR","PR")) * 100, 1),
        `CR Rate (%)` = round(mean(response_cat == "CR") * 100, 1),
        `Best Response (%)` = round(mean(best_response) * 100, 1),
        `Median PFS (days)` = round(median(pfs_time), 0),
        .groups = "drop"
      )

    datatable(summary_df, options = list(dom = "t", paging = FALSE),
              rownames = FALSE, class = "table-sm") %>%
      formatStyle("ORR (%)", backgroundColor = styleInterval(c(50, 70), c("#FFCCCC", "#FFFFCC", "#CCFFCC")))
  })

  # ── Biomarker Explorer Outputs ─────────────────────────────────────────────

  biomarker_data <- eventReactive(input$run_biomarker, {
    set.seed(789)
    n <- input$n_bio_patients

    # Generate virtual patient population
    vp <- data.frame(
      patient_id = 1:n,
      subtype = sample(c("GCB", "ABC"), n, replace = TRUE, prob = c(0.55, 0.45)),
      myc = sample(c(TRUE, FALSE), n, replace = TRUE, prob = c(0.15, 0.85)),
      bcl2_val = rbeta(n, 2, 3) * runif(n, 0.4, 1.2),
      btk_occ = runif(n, 0.3, 0.99),
      pdl1_val = rbeta(n, 1.5, 3),
      ki67_val = rbeta(n, 3, 2) * 100,
      ldh_val = rlnorm(n, log(350), 0.4),
      ipi_val = sample(0:5, n, replace = TRUE, prob = c(0.1, 0.2, 0.25, 0.25, 0.1, 0.1))
    ) %>%
      mutate(
        bcl2_pos = bcl2_val > 0.5,
        response_pct = case_when(
          input$biomarker_treatment == "R-CHOP+Ven" ~ -50 - 40 * bcl2_val + rnorm(n, 0, 15),
          input$biomarker_treatment == "Ibr+R-CHOP" ~ case_when(
            subtype == "ABC" ~ -55 - 20 * btk_occ + rnorm(n, 0, 20),
            TRUE ~ -35 - 15 * btk_occ + rnorm(n, 0, 20)
          ),
          input$biomarker_treatment == "CAR-T" ~ -60 - 15 * pdl1_val + rnorm(n, 0, 15),
          TRUE ~ -45 - 20 * (bcl2_val + pdl1_val) / 2 + rnorm(n, 0, 20)
        ),
        pfs_days = case_when(
          input$biomarker_treatment == "R-CHOP+Ven" ~ 400 + 300 * bcl2_val + rnorm(n, 0, 60),
          TRUE ~ 300 + 200 * btk_occ + rnorm(n, 0, 80)
        ),
        bcr_red = 0.3 + 0.5 * btk_occ + rnorm(n, 0, 0.1),
        immune_change = 0.1 + 0.3 * pdl1_val - 0.1 * bcl2_val + rnorm(n, 0, 0.05),
        response_cat = case_when(
          response_pct <= -80 ~ "CR",
          response_pct <= -30 ~ "PR",
          response_pct <= 20  ~ "SD",
          TRUE ~ "PD"
        )
      ) %>%
      mutate(pfs_days = pmax(30, pmin(1095, pfs_days)))

    vp
  }, ignoreNULL = FALSE)

  output$biomarker_scatter <- renderPlotly({
    df <- biomarker_data()

    x_var <- switch(input$biomarker_x,
                    "bcl2" = "bcl2_val", "btk" = "btk_occ", "pdl1" = "pdl1_val",
                    "ki67" = "ki67_val", "ldh" = "ldh_val", "ipi" = "ipi_val")
    y_var <- switch(input$biomarker_y,
                    "response" = "response_pct", "pfs" = "pfs_days",
                    "bcr_red" = "bcr_red", "immune_change" = "immune_change")
    color_var <- switch(input$biomarker_color,
                        "subtype" = "subtype", "myc" = "myc", "response_cat" = "response_cat")

    x_label <- switch(input$biomarker_x,
                      "bcl2" = "BCL2 Expression (rel)", "btk" = "BTK Occupancy (0-1)",
                      "pdl1" = "PD-L1 Expression (rel)", "ki67" = "Ki-67 Index (%)",
                      "ldh" = "LDH (U/L)", "ipi" = "IPI Score")
    y_label <- switch(input$biomarker_y,
                      "response" = "Best Response (%)", "pfs" = "PFS (days)",
                      "bcr_red" = "BCR Signaling Reduction", "immune_change" = "Immune Effector Change")

    p <- ggplot(df, aes_string(x = x_var, y = y_var, color = color_var)) +
      geom_point(alpha = 0.6, size = 1.5) +
      geom_smooth(method = "loess", se = TRUE, color = "black", linewidth = 0.8) +
      labs(x = x_label, y = y_label, title = paste(x_label, "vs", y_label),
           color = input$biomarker_color) +
      theme_minimal(base_size = 11) +
      scale_color_brewer(palette = "Set1")
    ggplotly(p)
  })

  output$biomarker_quartile <- renderPlotly({
    df <- biomarker_data()
    x_var <- switch(input$biomarker_x,
                    "bcl2" = "bcl2_val", "btk" = "btk_occ", "pdl1" = "pdl1_val",
                    "ki67" = "ki67_val", "ldh" = "ldh_val", "ipi" = "ipi_val")

    df$quartile <- cut(df[[x_var]], breaks = quantile(df[[x_var]], probs = 0:4/4),
                       labels = c("Q1 (Low)", "Q2", "Q3", "Q4 (High)"), include.lowest = TRUE)

    qdf <- df %>%
      group_by(quartile, response_cat) %>%
      summarise(n = n(), .groups = "drop") %>%
      group_by(quartile) %>%
      mutate(pct = n / sum(n) * 100)

    p <- ggplot(qdf, aes(x = quartile, y = pct, fill = response_cat)) +
      geom_bar(stat = "identity") +
      scale_fill_manual(values = c("CR" = "#27AE60", "PR" = "#2980B9",
                                   "SD" = "#F39C12", "PD" = "#E74C3C")) +
      labs(x = "Biomarker Quartile", y = "Response Rate (%)", fill = "Response",
           title = "Response by Biomarker Quartile") +
      theme_minimal(base_size = 11)
    ggplotly(p)
  })

  output$bcl2_ven_plot <- renderPlotly({
    df <- biomarker_data()
    p <- ggplot(df, aes(x = bcl2_val, y = response_pct, color = response_cat)) +
      geom_point(alpha = 0.5, size = 1.5) +
      geom_smooth(method = "lm", color = "darkblue", linewidth = 1) +
      geom_vline(xintercept = 0.5, linetype = "dashed", color = "red") +
      annotate("text", x = 0.55, y = max(df$response_pct) * 0.9,
               label = "BCL2+ cutoff", color = "red", size = 3, hjust = 0) +
      scale_color_manual(values = c("CR" = "#27AE60", "PR" = "#2980B9",
                                    "SD" = "#F39C12", "PD" = "#E74C3C")) +
      labs(x = "BCL2 Expression (relative)", y = "Best Response (%)",
           title = "BCL2 vs Venetoclax Response", color = "Response") +
      theme_minimal(base_size = 11)
    ggplotly(p)
  })

  output$btk_ibr_plot <- renderPlotly({
    ibr_doses <- c(140, 280, 420, 560, 840)
    btk_occ <- 1 - exp(-0.0025 * ibr_doses)  # Simple saturation model

    dose_df <- data.frame(dose = ibr_doses, btk_occupancy = btk_occ * 100)

    p <- ggplot(dose_df, aes(x = dose, y = btk_occupancy)) +
      geom_line(color = "#9B59B6", linewidth = 1.5) +
      geom_point(size = 3, color = "#8E44AD") +
      geom_hline(yintercept = 90, linetype = "dashed", color = "red") +
      annotate("text", x = 150, y = 91.5, label = "90% occupancy threshold",
               color = "red", size = 3) +
      labs(x = "Ibrutinib Dose (mg/day)", y = "BTK Occupancy (%)",
           title = "BTK Occupancy vs Ibrutinib Dose") +
      theme_minimal(base_size = 11)
    ggplotly(p)
  })

  # About tab - clinical trials table
  output$clinical_trials_table <- renderTable({
    data.frame(
      Trial = c("GOYA", "POLARIX", "REMoDL-B", "PHOENIX", "TRANSCEND", "JULIET", "ZUMA-1", "L-MIND"),
      Treatment = c("G-CHOP vs R-CHOP", "Pola-R-CHP vs R-CHOP", "R-CHOP±bortezomib",
                    "Ibr+R-CHOP vs R-CHOP (ABC)", "Liso-cel", "Tisagen-lecleucel",
                    "Axicab-cel", "Tafasitamab+Len"),
      Setting = c("1L", "1L", "1L", "1L ABC", "R/R ≥2L", "R/R ≥2L", "R/R ≥2L", "R/R 2L"),
      "Key Result" = c("Neg (HR 0.92)", "PFS HR 0.73", "Neg (no dbl-hit benefit)",
                       "Neg overall; ABC sig benefit (HR 0.75)", "ORR 73%, CRR 53%",
                       "ORR 52%, CRR 40%", "ORR 83%, CRR 58%", "ORR 60%, CRR 43%")
    )
  })
}

# ==============================================================================
# RUN APPLICATION
# ==============================================================================

shinyApp(ui = ui, server = server)
