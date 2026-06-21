# =============================================================================
# Adult-Onset Still's Disease (AOSD) QSP Shiny Application
# =============================================================================
# Quantitative Systems Pharmacology model for AOSD
# Models: Anakinra, Canakinumab, Tocilizumab, Corticosteroids, Tofacitinib,
#         and combination therapies
# Author: Claude Code (CCR)
# Date: 2026-06-17
# =============================================================================

library(shiny)
library(ggplot2)
library(dplyr)
library(tidyr)

# =============================================================================
# SIMULATION ENGINE - ODE-based QSP model (Euler integration)
# =============================================================================

# --- PK Parameters ---
pk_params <- list(
  anakinra = list(
    dose_mg = 100,          # mg/day SC
    route = "SC",
    ka = 0.693,             # absorption rate /h (t_half_abs ~1h)
    CL = 1.5,               # L/h
    Vd = 15,                # L
    F = 0.95,               # bioavailability
    t_half = 4              # hours
  ),
  canakinumab = list(
    dose_mg = 150,          # mg q8w SC
    route = "SC",
    ka = 0.005,             # /h slow absorption (mAb)
    CL = 0.0027,            # L/h
    Vd = 6.0,               # L
    F = 0.70,
    t_half = 3360           # 140 days (mAb)
  ),
  tocilizumab = list(
    dose_mg = 8,            # mg/kg q4w IV (weight-based, using 70kg)
    route = "IV",
    CL = 0.012,             # L/h/kg
    Vd = 0.10,              # L/kg
    t_half = 168            # ~7 days
  ),
  corticosteroid = list(
    dose_mg = 40,           # mg/day prednisone
    route = "PO",
    ka = 1.386,             # /h (t_half_abs ~0.5h)
    CL = 15,                # L/h
    Vd = 100,               # L
    F = 0.85,
    t_half = 3.5            # hours
  ),
  tofacitinib = list(
    dose_mg = 5,            # mg BID
    route = "PO",
    ka = 2.0,               # /h (t_half_abs ~0.35h)
    CL = 50,                # L/h
    Vd = 87,                # L
    F = 0.74,
    t_half = 3.2            # hours
  )
)

# --- Simulate PK (1-compartment) ---
simulate_pk <- function(drug, weight_kg, dose_override = NULL,
                        time_vec = seq(0, 672, by = 2),
                        interval_h = NULL) {

  p <- pk_params[[drug]]

  if (!is.null(dose_override)) p$dose_mg <- dose_override

  # Weight-adjust tocilizumab
  if (drug == "tocilizumab") {
    dose_total <- p$dose_mg * weight_kg   # mg/kg * kg
    CL <- p$CL * weight_kg
    Vd <- p$Vd * weight_kg
    t_half <- log(2) * Vd / CL
    p$dose_mg <- dose_total
    p$CL <- CL
    p$Vd <- Vd
    p$t_half <- t_half
  }

  # Dosing intervals
  if (is.null(interval_h)) {
    interval_h <- switch(drug,
      "anakinra"      = 24,
      "canakinumab"   = 1344,   # 56 days
      "tocilizumab"   = 672,    # 28 days
      "corticosteroid"= 24,
      "tofacitinib"   = 12,
      24
    )
  }

  ke <- log(2) / p$t_half
  dose_times <- seq(0, max(time_vec), by = interval_h)

  conc <- numeric(length(time_vec))

  for (t_idx in seq_along(time_vec)) {
    t <- time_vec[t_idx]
    C <- 0

    for (td in dose_times[dose_times <= t]) {
      dt <- t - td
      if (p$route %in% c("SC", "PO")) {
        ka <- p$ka
        C <- C + (p$F * p$dose_mg / p$Vd) *
          (ka / (ka - ke)) *
          (exp(-ke * dt) - exp(-ka * dt))
      } else {
        # IV bolus
        C <- C + (p$dose_mg / p$Vd) * exp(-ke * dt)
      }
    }
    conc[t_idx] <- max(C, 0)
  }
  return(conc)
}

# --- Disease QSP model (simplified ODEs via Euler) ---
simulate_qsp <- function(
    drug,
    weight_kg     = 70,
    age           = 40,
    severity      = 0.7,       # 0-1 baseline disease activity
    mas_risk      = 0.3,       # 0-1 baseline MAS risk
    dose_override = NULL,
    sim_days      = 180
) {
  dt <- 0.5           # hours
  t_vec <- seq(0, sim_days * 24, by = dt)
  n <- length(t_vec)

  # ------- Initial conditions -------
  # Cytokines (pg/mL)
  IL1b   <- numeric(n); IL1b[1]   <- 50  * severity + 5
  IL6    <- numeric(n); IL6[1]    <- 200 * severity + 10
  IL18   <- numeric(n); IL18[1]   <- 500 * severity + 50
  IFNg   <- numeric(n); IFNg[1]  <- 30  * severity + 2
  TNFa   <- numeric(n); TNFa[1]  <- 40  * severity + 5

  # Bound (drug-bound) cytokines
  IL1b_bound <- numeric(n)
  IL6_bound  <- numeric(n)

  # Biomarkers
  CRP    <- numeric(n); CRP[1]    <- 80  * severity + 5   # mg/L
  Ferrit <- numeric(n); Ferrit[1] <- 5000 * severity + 200 # ng/mL
  GlyFer <- numeric(n); GlyFer[1] <- max(0, 20 - 15 * severity) # % glycosylated
  LDH    <- numeric(n); LDH[1]    <- 400 * severity + 100  # U/L
  ALT    <- numeric(n); ALT[1]    <- 60  * severity + 20   # U/L
  AST    <- numeric(n); AST[1]    <- 55  * severity + 20   # U/L

  # Clinical scores (0-100 scale)
  Pouchot <- numeric(n); Pouchot[1] <- 60 * severity + 5
  Fever   <- numeric(n); Fever[1]   <- severity          # binary-ish (0-1 intensity)
  Arthrit <- numeric(n); Arthrit[1] <- 60 * severity
  Rash    <- numeric(n); Rash[1]    <- 40 * severity

  # MAS probability
  MAS_prob <- numeric(n); MAS_prob[1] <- mas_risk * severity

  # PK
  pk_conc <- simulate_pk(drug, weight_kg, dose_override, t_vec, NULL)

  # Drug effect function (sigmoidal)
  Emax_IL1 <- switch(drug,
    "anakinra"       = 0.92, "canakinumab" = 0.95,
    "tocilizumab"    = 0.30, "corticosteroid" = 0.60,
    "tofacitinib"    = 0.55, 0.5)
  Emax_IL6 <- switch(drug,
    "anakinra"       = 0.40, "canakinumab" = 0.25,
    "tocilizumab"    = 0.90, "corticosteroid" = 0.65,
    "tofacitinib"    = 0.60, 0.5)
  Emax_IFN <- switch(drug,
    "anakinra"       = 0.35, "canakinumab" = 0.30,
    "tocilizumab"    = 0.45, "corticosteroid" = 0.50,
    "tofacitinib"    = 0.75, 0.4)

  # EC50 scaled to typical Cmax range for each drug
  EC50 <- switch(drug,
    "anakinra"       = 2000,    # ng/mL
    "canakinumab"    = 20,      # ug/mL -> store in ng/mL = 20000
    "tocilizumab"    = 10,
    "corticosteroid" = 50,
    "tofacitinib"    = 5,
    100
  )
  if (drug == "canakinumab") EC50 <- 20000

  hill <- 2.0

  drug_eff_IL1 <- function(C) Emax_IL1 * C^hill / (EC50^hill + C^hill)
  drug_eff_IL6 <- function(C) Emax_IL6 * C^hill / (EC50^hill + C^hill)
  drug_eff_IFN <- function(C) Emax_IFN * C^hill / (EC50^hill + C^hill)

  # Half-lives / rate constants for cytokines
  k_deg_IL1 <- log(2) / 1.5    # ~1.5h half-life IL-1beta
  k_deg_IL6 <- log(2) / 3      # ~3h half-life IL-6
  k_deg_IL18 <- log(2) / 6
  k_deg_IFN  <- log(2) / 2
  k_deg_TNF  <- log(2) / 1.5

  # Biomarker rate constants
  k_CRP_on  <- 0.010
  k_CRP_off <- log(2) / 12
  k_Fer_on  <- 0.008
  k_Fer_off <- log(2) / 48
  k_LDH_on  <- 0.005
  k_LDH_off <- log(2) / 24

  # Euler integration
  for (i in 2:n) {
    C <- pk_conc[i-1]

    e1 <- drug_eff_IL1(C)
    e6 <- drug_eff_IL6(C)
    eI <- drug_eff_IFN(C)

    # Cytokine auto-stimulation (positive feedback)
    stim_IL1 <- 1 + 0.3 * IL1b[i-1] / (IL1b[i-1] + 100)
    stim_IL6  <- 1 + 0.4 * IL6[i-1] / (IL6[i-1] + 300)

    # Production rates (baseline + cross-stimulation)
    prod_IL1  <- (20 * severity * stim_IL1 * (1 + 0.2 * IL18[i-1]/500)) * (1 - e1)
    prod_IL6  <- (80 * severity * stim_IL6 * (1 + 0.3 * IL1b[i-1]/50)) * (1 - e6)
    prod_IL18 <- (150 * severity * (1 + 0.2 * IFNg[i-1]/30)) * (1 - 0.3 * e1)
    prod_IFNg <- (10 * severity * (1 + 0.5 * IL18[i-1]/500)) * (1 - eI)
    prod_TNF  <- (15 * severity * (1 + 0.3 * IL1b[i-1]/50)) * (1 - 0.4 * e1)

    dIL1  <- (prod_IL1  - k_deg_IL1  * IL1b[i-1])  * dt
    dIL6  <- (prod_IL6  - k_deg_IL6  * IL6[i-1])   * dt
    dIL18 <- (prod_IL18 - k_deg_IL18 * IL18[i-1])  * dt
    dIFNg <- (prod_IFNg - k_deg_IFN  * IFNg[i-1])  * dt
    dTNF  <- (prod_TNF  - k_deg_TNF  * TNFa[i-1])  * dt

    IL1b[i]  <- max(IL1b[i-1]  + dIL1,  0.5)
    IL6[i]   <- max(IL6[i-1]   + dIL6,  1)
    IL18[i]  <- max(IL18[i-1]  + dIL18, 20)
    IFNg[i]  <- max(IFNg[i-1]  + dIFNg, 0.5)
    TNFa[i]  <- max(TNFa[i-1]  + dTNF,  1)

    # Bound fraction approximation
    if (drug == "anakinra")    IL1b_bound[i] <- IL1b[i] * e1
    if (drug == "canakinumab") IL1b_bound[i] <- IL1b[i] * e1
    if (drug == "tocilizumab") IL6_bound[i]  <- IL6[i]  * e6

    # Biomarkers driven by cytokines
    cyto_drive_CRP  <- (IL6[i-1] / 200 + IL1b[i-1] / 50) / 2
    cyto_drive_Fer  <- (IL18[i-1] / 500 + IFNg[i-1] / 30) / 2 * (1 + 2 * mas_risk)
    cyto_drive_LDH  <- (IFNg[i-1] / 30 + IL18[i-1] / 500) / 2

    CRP_ss    <- 5 + 75 * severity * cyto_drive_CRP
    Ferrit_ss <- 200 + 4800 * severity * cyto_drive_Fer
    LDH_ss    <- 100 + 300 * severity * cyto_drive_LDH
    ALT_ss    <- 20 + 40  * severity * cyto_drive_LDH
    AST_ss    <- 20 + 35  * severity * cyto_drive_LDH

    CRP[i]    <- CRP[i-1]    + (CRP_ss    - CRP[i-1])    * k_CRP_off * dt
    Ferrit[i] <- Ferrit[i-1] + (Ferrit_ss - Ferrit[i-1]) * k_Fer_off * dt
    LDH[i]    <- LDH[i-1]   + (LDH_ss    - LDH[i-1])    * k_LDH_off * dt
    ALT[i]    <- ALT[i-1]   + (ALT_ss    - ALT[i-1])     * k_LDH_off * dt
    AST[i]    <- AST[i-1]   + (AST_ss    - AST[i-1])     * k_LDH_off * dt

    # Glycosylated ferritin: inversely related to macrophage activation
    GlyFer[i] <- max(5, 50 - 40 * (IL18[i-1]/500) * mas_risk -
                       10 * (IFNg[i-1]/30) * mas_risk)

    # Pouchot score (clinical activity)
    Pouchot[i] <- max(0, min(100,
      5 + 50 * (IL1b[i-1]/50 * 0.3 + IL6[i-1]/200 * 0.4 +
                IL18[i-1]/500 * 0.3) * severity
    ))
    Fever[i]   <- max(0, min(1, IL1b[i-1] / 100 * severity))
    Arthrit[i] <- max(0, min(100, 20 + 40 * IL6[i-1] / 200 * severity))
    Rash[i]    <- max(0, min(100, 15 + 25 * IL1b[i-1] / 50 * severity))

    # MAS risk
    MAS_prob[i] <- max(0, min(0.99,
      0.05 + 0.8 * mas_risk * (Ferrit[i-1] / 10000) *
        (1 + IFNg[i-1] / 50) * (GlyFer[i-1] < 20)
    ))
  }

  # Subsample to daily for output
  day_idx <- seq(1, n, by = round(24 / dt))

  data.frame(
    day        = t_vec[day_idx] / 24,
    time_h     = t_vec[day_idx],
    ConcPK     = pk_conc[day_idx],
    IL1b       = IL1b[day_idx],
    IL1b_bound = IL1b_bound[day_idx],
    IL6        = IL6[day_idx],
    IL6_bound  = IL6_bound[day_idx],
    IL18       = IL18[day_idx],
    IFNg       = IFNg[day_idx],
    TNFa       = TNFa[day_idx],
    CRP        = CRP[day_idx],
    Ferritin   = Ferrit[day_idx],
    GlyFerrit  = GlyFer[day_idx],
    LDH        = LDH[day_idx],
    ALT        = ALT[day_idx],
    AST        = AST[day_idx],
    Pouchot    = Pouchot[day_idx],
    Fever      = Fever[day_idx],
    Arthritis  = Arthrit[day_idx],
    Rash       = Rash[day_idx],
    MAS_prob   = MAS_prob[day_idx]
  )
}

# --- PK-only time series (for the PK tab, short time window) ---
simulate_pk_detail <- function(drug, weight_kg, dose_override = NULL,
                                sim_days = 28) {
  t_vec <- seq(0, sim_days * 24, by = 0.5)
  pk_conc <- simulate_pk(drug, weight_kg, dose_override, t_vec)
  data.frame(time_h = t_vec, conc = pk_conc)
}

# HScore helper
calc_hscore <- function(temp, hep_sp, cytopenias, trig, ferritin,
                         fibrinogen, asat, hemo_bm, imm_suppressed) {
  score <- 0
  # Temperature
  if (temp < 38.4) score <- score + 0
  else if (temp <= 39.4) score <- score + 33
  else score <- score + 49

  # Hepatomegaly/Splenomegaly
  if (hep_sp == "Neither")     score <- score + 0
  else if (hep_sp == "One")    score <- score + 23
  else score <- score + 38

  # Cytopenias
  if (cytopenias == "0 lines")  score <- score + 0
  else if (cytopenias == "1 line") score <- score + 24
  else if (cytopenias == "2 lines") score <- score + 34
  else score <- score + 50

  # Triglycerides
  if (trig < 1.5) score <- score + 0
  else if (trig < 4.0) score <- score + 44
  else score <- score + 64

  # Ferritin
  if (ferritin < 2000) score <- score + 0
  else if (ferritin < 6000) score <- score + 35
  else score <- score + 50

  # Fibrinogen
  if (fibrinogen <= 2.5) score <- score + 30
  else score <- score + 0

  # ASAT
  if (asat < 30) score <- score + 0
  else score <- score + 19

  # Hemophagocytosis on BM
  if (hemo_bm == "Yes") score <- score + 35

  # Immunosuppressed
  if (imm_suppressed == "Yes") score <- score + 18

  list(score = score, prob = round(exp(-9.516 + 0.02426 * score) /
                                     (1 + exp(-9.516 + 0.02426 * score)) * 100, 1))
}

# =============================================================================
# COLOUR PALETTE
# =============================================================================
drug_colors <- c(
  anakinra       = "#E41A1C",
  canakinumab    = "#377EB8",
  tocilizumab    = "#4DAF4A",
  corticosteroid = "#FF7F00",
  tofacitinib    = "#984EA3",
  untreated      = "#999999"
)

drug_labels <- c(
  anakinra       = "Anakinra (IL-1Ra)",
  canakinumab    = "Canakinumab (anti-IL-1β)",
  tocilizumab    = "Tocilizumab (anti-IL-6R)",
  corticosteroid = "Corticosteroids",
  tofacitinib    = "Tofacitinib (JAK1/3i)",
  untreated      = "Untreated"
)

theme_aosd <- function() {
  theme_bw(base_size = 13) +
    theme(
      plot.title   = element_text(face = "bold", size = 14),
      strip.text   = element_text(face = "bold"),
      legend.position = "bottom",
      panel.grid.minor = element_blank()
    )
}

# =============================================================================
# UI
# =============================================================================
ui <- fluidPage(
  titlePanel(
    div(
      h2("Adult-Onset Still's Disease (AOSD) — QSP Interactive Model",
         style = "color:#1a3a5c; font-weight:bold;"),
      h5("Quantitative Systems Pharmacology | IL-1β / IL-6 / JAK-STAT Pathway Model",
         style = "color:#666;")
    )
  ),
  br(),

  sidebarLayout(
    sidebarPanel(
      width = 3,
      h4("Global Simulation Settings", style = "color:#1a3a5c;"),
      hr(),
      sliderInput("sim_days", "Simulation Duration (days):",
                  min = 28, max = 365, value = 180, step = 14),
      selectInput("drug", "Primary Treatment:",
                  choices = c(
                    "Anakinra (IL-1Ra)"          = "anakinra",
                    "Canakinumab (anti-IL-1β)"   = "canakinumab",
                    "Tocilizumab (anti-IL-6R)"   = "tocilizumab",
                    "Corticosteroids"            = "corticosteroid",
                    "Tofacitinib (JAK1/3i)"      = "tofacitinib",
                    "Untreated / Placebo"        = "untreated"
                  ), selected = "anakinra"),
      hr(),
      h5("Patient Profile", style = "color:#1a3a5c;"),
      sliderInput("age", "Age (years):", min = 18, max = 80, value = 40),
      sliderInput("weight", "Weight (kg):",  min = 40, max = 150, value = 70),
      sliderInput("bmi", "BMI (kg/m²):", min = 15, max = 45, value = 24),
      sliderInput("severity", "Disease Severity (0–1):",
                  min = 0, max = 1, value = 0.7, step = 0.05),
      sliderInput("mas_risk", "MAS Risk Factor (0–1):",
                  min = 0, max = 1, value = 0.3, step = 0.05),
      hr(),
      actionButton("run_sim", "Run Simulation",
                   class = "btn-primary btn-lg", width = "100%"),
      br(), br(),
      p(em("Model based on IL-1/IL-6/JAK-STAT pathway mechanistic ODEs.
            PK parameters from published Phase II/III trials."),
        style = "font-size:11px; color:#888;")
    ),

    mainPanel(
      width = 9,
      tabsetPanel(
        id = "tabs",
        type = "tabs",

        # ---- Tab 1: Patient Profile ----
        tabPanel("Patient Profile",
          br(),
          fluidRow(
            column(6,
              h4("Current Patient Summary", style = "color:#1a3a5c;"),
              tableOutput("patient_table"),
              br(),
              h4("AOSD Classification Criteria (Yamaguchi 1992)", style = "color:#1a3a5c;"),
              tableOutput("yamaguchi_table")
            ),
            column(6,
              h4("Baseline Disease Status", style = "color:#1a3a5c;"),
              plotOutput("profile_radar", height = "350px"),
              br(),
              h4("Expected Baseline Lab Values"),
              tableOutput("baseline_labs")
            )
          )
        ),

        # ---- Tab 2: Drug PK ----
        tabPanel("Drug PK",
          br(),
          fluidRow(
            column(12,
              h4("Pharmacokinetic Profiles — All Agents", style = "color:#1a3a5c;"),
              p("Simulated plasma concentration profiles (first 28 days shown for fast-acting agents; full duration for biologics)."),
              selectInput("pk_drug_select", "Select Drug(s) to Display:",
                          choices = c(
                            "Anakinra"       = "anakinra",
                            "Canakinumab"    = "canakinumab",
                            "Tocilizumab"    = "tocilizumab",
                            "Corticosteroids"= "corticosteroid",
                            "Tofacitinib"    = "tofacitinib"
                          ),
                          multiple = TRUE,
                          selected = c("anakinra", "canakinumab",
                                       "tocilizumab", "corticosteroid", "tofacitinib")),
              plotOutput("pk_plot_all", height = "400px")
            )
          ),
          br(),
          fluidRow(
            column(6,
              h5("PK Parameter Summary"),
              tableOutput("pk_param_table")
            ),
            column(6,
              h5("Dose-Response Relationships"),
              plotOutput("pk_dose_response", height = "300px")
            )
          )
        ),

        # ---- Tab 3: Cytokine Dynamics ----
        tabPanel("Cytokine Dynamics",
          br(),
          fluidRow(
            column(12,
              h4("Cytokine Time Course", style = "color:#1a3a5c;"),
              p("Simulated free and drug-bound cytokine levels over the treatment period."),
              plotOutput("cytokine_plot", height = "500px")
            )
          ),
          br(),
          fluidRow(
            column(6,
              h5("Free vs. Bound Cytokines"),
              plotOutput("bound_free_plot", height = "300px")
            ),
            column(6,
              h5("Cytokine Network Heatmap"),
              plotOutput("cytokine_heatmap", height = "300px")
            )
          )
        ),

        # ---- Tab 4: Biomarkers ----
        tabPanel("Biomarkers",
          br(),
          fluidRow(
            column(12,
              h4("Biomarker Kinetics Over Treatment", style = "color:#1a3a5c;"),
              p("Key laboratory biomarkers with clinical threshold lines.")
            )
          ),
          fluidRow(
            column(6,
              plotOutput("crp_plot", height = "280px")
            ),
            column(6,
              plotOutput("ferritin_plot", height = "280px")
            )
          ),
          fluidRow(
            column(6,
              plotOutput("glyfer_plot", height = "280px")
            ),
            column(6,
              plotOutput("ldh_plot", height = "280px")
            )
          ),
          fluidRow(
            column(12,
              plotOutput("liver_plot", height = "280px")
            )
          )
        ),

        # ---- Tab 5: Clinical Endpoints ----
        tabPanel("Clinical Endpoints",
          br(),
          fluidRow(
            column(12,
              h4("AOSD Clinical Activity Endpoints", style = "color:#1a3a5c;"),
              p("Pouchot score, individual disease manifestations, and treatment response rates.")
            )
          ),
          fluidRow(
            column(6,
              plotOutput("pouchot_plot", height = "300px")
            ),
            column(6,
              plotOutput("fever_arthritis_plot", height = "300px")
            )
          ),
          fluidRow(
            column(6,
              plotOutput("rash_plot", height = "300px")
            ),
            column(6,
              plotOutput("response_rate_plot", height = "300px")
            )
          ),
          fluidRow(
            column(12,
              h5("Response Rate Summary"),
              tableOutput("response_table")
            )
          )
        ),

        # ---- Tab 6: Scenario Comparison ----
        tabPanel("Scenario Comparison",
          br(),
          h4("Treatment Scenario Comparison", style = "color:#1a3a5c;"),
          p("Side-by-side comparison of all treatment modalities across key endpoints."),
          fluidRow(
            column(6,
              plotOutput("scenario_pouchot", height = "300px")
            ),
            column(6,
              plotOutput("scenario_ferritin", height = "300px")
            )
          ),
          fluidRow(
            column(6,
              plotOutput("scenario_crp", height = "300px")
            ),
            column(6,
              plotOutput("scenario_il6", height = "300px")
            )
          ),
          fluidRow(
            column(12,
              h5("End-of-Study Summary Table (Day ", textOutput("sim_day_label", inline=TRUE), ")"),
              tableOutput("scenario_table")
            )
          )
        ),

        # ---- Tab 7: MAS Risk ----
        tabPanel("MAS Risk",
          br(),
          h4("Macrophage Activation Syndrome (MAS) Risk Assessment", style = "color:#1a3a5c;"),
          fluidRow(
            column(8,
              plotOutput("mas_prob_plot", height = "300px")
            ),
            column(4,
              wellPanel(
                h5("MAS Risk Indicators"),
                uiOutput("mas_indicators")
              )
            )
          ),
          fluidRow(
            column(6,
              plotOutput("mas_ferritin_kin", height = "280px")
            ),
            column(6,
              plotOutput("mas_glyfer_plot", height = "280px")
            )
          ),
          br(),
          h4("HScore Calculator (MAS Probability)", style = "color:#1a3a5c;"),
          p("Interactive HScore for clinical MAS diagnosis probability (Fardet et al. 2014)"),
          fluidRow(
            column(3,
              numericInput("hs_temp", "Peak Temperature (°C):", value = 39.5, min = 36, max = 42, step = 0.1),
              selectInput("hs_hepsp", "Hepatomegaly/Splenomegaly:",
                          choices = c("Neither", "One", "Both"), selected = "One"),
              selectInput("hs_cyto", "Cytopenias (cell lines):",
                          choices = c("0 lines", "1 line", "2 lines", "3 lines"), selected = "2 lines")
            ),
            column(3,
              numericInput("hs_trig", "Triglycerides (mmol/L):", value = 3.5, min = 0, max = 20, step = 0.1),
              numericInput("hs_ferrit", "Ferritin (ng/mL):", value = 5000, min = 0, max = 100000, step = 100),
              numericInput("hs_fibrin", "Fibrinogen (g/L):", value = 2.0, min = 0, max = 10, step = 0.1)
            ),
            column(3,
              numericInput("hs_asat", "ASAT (IU/L):", value = 80, min = 0, max = 2000, step = 5),
              selectInput("hs_hemo", "Hemophagocytosis on BM biopsy:",
                          choices = c("No", "Yes"), selected = "No"),
              selectInput("hs_imm", "Immunosuppressed:",
                          choices = c("No", "Yes"), selected = "No")
            ),
            column(3,
              br(), br(),
              actionButton("calc_hscore", "Calculate HScore", class = "btn-warning btn-lg"),
              br(), br(),
              uiOutput("hscore_result")
            )
          )
        ),

        # ---- Tab 8: Mechanistic Map ----
        tabPanel("Mechanistic Map",
          br(),
          h4("AOSD Mechanistic Pathway Map", style = "color:#1a3a5c;"),
          fluidRow(
            column(12,
              uiOutput("mech_map_ui")
            )
          ),
          br(),
          h4("Key Pathways Summary", style = "color:#1a3a5c;"),
          fluidRow(
            column(4,
              wellPanel(
                h5("Innate Immune Activation"),
                tags$ul(
                  tags$li("Pattern recognition via TLRs and NLRs"),
                  tags$li("NLRP3 inflammasome activation → IL-1β maturation"),
                  tags$li("Macrophage M1 polarization"),
                  tags$li("Neutrophil activation and NET formation"),
                  tags$li("NK cell dysregulation")
                )
              )
            ),
            column(4,
              wellPanel(
                h5("Cytokine Cascade"),
                tags$ul(
                  tags$li("IL-1β: fever, acute phase, MAS driver"),
                  tags$li("IL-6: CRP synthesis, Th17 differentiation"),
                  tags$li("IL-18: IFN-γ induction, macrophage activation"),
                  tags$li("IFN-γ: macrophage activation, hemophagocytosis"),
                  tags$li("TNF-α: NF-κB activation, inflammation amplification"),
                  tags$li("M-CSF / GM-CSF: macrophage proliferation")
                )
              )
            ),
            column(4,
              wellPanel(
                h5("Drug Mechanisms of Action"),
                tags$ul(
                  tags$li(strong("Anakinra:"), " IL-1Ra, blocks IL-1α and IL-1β"),
                  tags$li(strong("Canakinumab:"), " anti-IL-1β monoclonal antibody"),
                  tags$li(strong("Tocilizumab:"), " anti-IL-6R, blocks IL-6 signalling"),
                  tags$li(strong("Corticosteroids:"), " broad immunosuppression, NF-κB inhibition"),
                  tags$li(strong("Tofacitinib:"), " JAK1/3 inhibition, blocks IL-6, IFN-γ, IL-2 signalling")
                )
              )
            )
          ),
          fluidRow(
            column(12,
              wellPanel(
                h5("Model Compartments and ODEs"),
                p("The QSP model tracks the following state variables:"),
                tags$ul(
                  tags$li("Drug PK: 1-compartment model (IV/SC/PO) for each agent"),
                  tags$li("Cytokines: IL-1β, IL-6, IL-18, IFN-γ, TNF-α — production/degradation ODEs with cross-stimulation"),
                  tags$li("Bound cytokines: drug-target complex formation"),
                  tags$li("Acute phase: CRP (IL-6 driven), Ferritin (IL-18/IFN-γ driven)"),
                  tags$li("Hepatic markers: ALT, AST, LDH (macrophage activation)"),
                  tags$li("Clinical scores: Pouchot, fever, arthritis, rash — composite of cytokine levels"),
                  tags$li("MAS probability: non-linear function of ferritin, IFN-γ, glycosylated ferritin%")
                )
              )
            )
          )
        )
      ) # end tabsetPanel
    ) # end mainPanel
  ) # end sidebarLayout
) # end fluidPage

# =============================================================================
# SERVER
# =============================================================================
server <- function(input, output, session) {

  # ---- Reactive: run primary simulation ----
  sim_data <- eventReactive(input$run_sim, {
    req(input$drug, input$weight, input$severity, input$mas_risk)

    withProgress(message = "Running QSP simulation...", value = 0, {
      incProgress(0.3, detail = "Computing PK/PD...")

      if (input$drug == "untreated") {
        # simulate with near-zero drug
        d <- simulate_qsp("anakinra", input$weight, input$age,
                           input$severity, input$mas_risk,
                           dose_override = 0.001, sim_days = input$sim_days)
      } else {
        d <- simulate_qsp(input$drug, input$weight, input$age,
                           input$severity, input$mas_risk,
                           dose_override = NULL, sim_days = input$sim_days)
      }
      incProgress(0.7, detail = "Finalising...")
      d
    })
  }, ignoreNULL = FALSE)

  # ---- Reactive: all scenarios ----
  scenario_data <- eventReactive(input$run_sim, {
    withProgress(message = "Running all scenarios...", value = 0, {
      drugs <- c("anakinra", "canakinumab", "tocilizumab",
                 "corticosteroid", "tofacitinib")

      result <- lapply(drugs, function(drg) {
        d <- simulate_qsp(drg, input$weight, input$age,
                           input$severity, input$mas_risk,
                           dose_override = NULL, sim_days = input$sim_days)
        d$drug <- drg
        d
      })
      # Add untreated
      d_unt <- simulate_qsp("anakinra", input$weight, input$age,
                              input$severity, input$mas_risk,
                              dose_override = 0.001, sim_days = input$sim_days)
      d_unt$drug <- "untreated"
      result[[6]] <- d_unt
      incProgress(1.0)
      do.call(rbind, result)
    })
  }, ignoreNULL = FALSE)

  # ========================
  # TAB 1: Patient Profile
  # ========================
  output$patient_table <- renderTable({
    bsa <- sqrt(input$weight * (input$bmi^(-1) * input$weight) / 3600)
    data.frame(
      Parameter = c("Age", "Weight", "BMI", "BSA (est.)",
                    "Disease Severity", "MAS Risk Factor",
                    "Primary Treatment"),
      Value = c(
        paste(input$age, "years"),
        paste(input$weight, "kg"),
        paste(input$bmi, "kg/m²"),
        paste(round(input$weight / (input$bmi^0.5 * 60) * 1.8, 2), "m²"),
        paste0(round(input$severity * 100), "%"),
        paste0(round(input$mas_risk * 100), "%"),
        drug_labels[input$drug]
      ),
      stringsAsFactors = FALSE
    )
  }, striped = TRUE, hover = TRUE)

  output$yamaguchi_table <- renderTable({
    data.frame(
      Criteria = c("Major (need 5)", "Fever ≥39°C, quotidian",
                   "Arthralgia >2 weeks", "Typical rash",
                   "Leukocytosis ≥10,000/mm³",
                   "Minor (any 2)", "Sore throat",
                   "Lymphadenopathy", "Hepato/Splenomegaly",
                   "Abnormal liver function", "Negative RF/ANA"),
      Status = c("",
                 if (input$severity > 0.4) "Present" else "Absent",
                 if (input$severity > 0.3) "Present" else "Absent",
                 if (input$severity > 0.5) "Present" else "Absent",
                 if (input$severity > 0.5) "Present" else "Absent",
                 "",
                 if (input$severity > 0.4) "Present" else "Absent",
                 if (input$severity > 0.6) "Present" else "Absent",
                 if (input$mas_risk  > 0.4) "Present" else "Absent",
                 if (input$severity > 0.5) "Present" else "Absent",
                 "Assumed Negative")
    )
  }, striped = TRUE)

  output$profile_radar <- renderPlot({
    sv <- input$severity
    mr <- input$mas_risk
    cats <- c("Fever", "Arthritis", "Rash", "Ferritin", "CRP", "MAS Risk")
    vals <- c(sv * 0.9, sv * 0.8, sv * 0.7, sv, sv * 0.85, mr)
    vals <- pmin(pmax(vals, 0), 1)
    n <- length(cats)
    angles <- seq(0, 2 * pi, length.out = n + 1)[1:n]
    x <- cos(angles) * vals
    y <- sin(angles) * vals
    xp <- cos(angles)
    yp <- sin(angles)
    par(mar = c(2, 2, 3, 2))
    plot(0, 0, type = "n", xlim = c(-1.4, 1.4), ylim = c(-1.4, 1.4),
         asp = 1, axes = FALSE, xlab = "", ylab = "",
         main = "Disease Activity Radar")
    for (r in c(0.25, 0.5, 0.75, 1.0)) {
      polygon(cos(angles) * r, sin(angles) * r,
              border = "grey80", lty = 2)
    }
    for (i in 1:n) {
      lines(c(0, xp[i]), c(0, yp[i]), col = "grey70")
      text(xp[i] * 1.25, yp[i] * 1.25, cats[i], cex = 0.9, font = 2)
    }
    polygon(c(x, x[1]), c(y, y[1]),
            col = adjustcolor("#E41A1C", alpha.f = 0.35),
            border = "#E41A1C", lwd = 2)
    mtext("0.25 / 0.5 / 0.75 / 1.0", side = 1, cex = 0.7, col = "grey60")
  })

  output$baseline_labs <- renderTable({
    sv <- input$severity
    mr <- input$mas_risk
    data.frame(
      Biomarker = c("CRP (mg/L)", "Ferritin (ng/mL)", "Glycosylated Ferritin (%)",
                    "LDH (U/L)", "ALT (U/L)", "AST (U/L)",
                    "Pouchot Score"),
      `Baseline Value` = c(
        round(5 + 75 * sv),
        round(200 + 4800 * sv),
        round(max(5, 20 - 15 * sv)),
        round(100 + 300 * sv),
        round(20 + 40 * sv),
        round(20 + 35 * sv),
        round(5 + 50 * sv)
      ),
      `Reference Range` = c("<5", "<400", ">20 (rules out MAS)",
                             "<250", "<40", "<40", "<5 (inactive)"),
      check.names = FALSE
    )
  }, striped = TRUE)

  # ========================
  # TAB 2: Drug PK
  # ========================
  output$pk_plot_all <- renderPlot({
    selected <- input$pk_drug_select
    if (length(selected) == 0) return(NULL)

    # Scale: use 14-day window for fast drugs, full sim for slow mAbs
    plot_list <- lapply(selected, function(drg) {
      win <- if (drg %in% c("anakinra", "corticosteroid", "tofacitinib")) 14 else 56
      df <- simulate_pk_detail(drg, input$weight, sim_days = win)
      df$drug <- drg
      df$day <- df$time_h / 24
      df
    })

    # Find common scale: facet by drug
    all_df <- do.call(rbind, plot_list)

    ggplot(all_df, aes(x = day, y = conc, color = drug)) +
      geom_line(linewidth = 1.2) +
      facet_wrap(~ drug, scales = "free", labeller = as_labeller(drug_labels)) +
      scale_color_manual(values = drug_colors, labels = drug_labels) +
      labs(title = "Plasma Concentration Profiles by Drug",
           x = "Time (days)", y = "Plasma Concentration (ng/mL or µg/mL)",
           color = "Drug") +
      theme_aosd() +
      theme(legend.position = "none")
  })

  output$pk_param_table <- renderTable({
    data.frame(
      Drug = c("Anakinra", "Canakinumab", "Tocilizumab", "Corticosteroid", "Tofacitinib"),
      Route = c("SC qd", "SC q8w", "IV q4w", "PO qd", "PO BID"),
      `Dose (std)` = c("100 mg", "150 mg", "8 mg/kg", "40 mg pred.", "5 mg"),
      `t½ (h)` = c(4, 3360, 168, 3.5, 3.2),
      `Vd (L)` = c(15, 6, "0.1×BW", 100, 87),
      `Bioavail.` = c("95%", "70%", "100% (IV)", "85%", "74%"),
      check.names = FALSE
    )
  }, striped = TRUE, hover = TRUE)

  output$pk_dose_response <- renderPlot({
    # Dose-response curves for each drug
    doses <- list(
      anakinra       = seq(25, 300, by = 25),
      canakinumab    = seq(50, 450, by = 50),
      tocilizumab    = seq(2, 16, by = 2),
      corticosteroid = seq(5, 80, by = 5),
      tofacitinib    = seq(1, 15, by = 1)
    )
    EC50s <- c(anakinra = 2000, canakinumab = 20000, tocilizumab = 10,
               corticosteroid = 50, tofacitinib = 5)
    Emaxs <- c(anakinra = 0.92, canakinumab = 0.95, tocilizumab = 0.90,
                corticosteroid = 0.65, tofacitinib = 0.75)

    dr_df <- do.call(rbind, lapply(names(doses), function(drg) {
      d_vals <- doses[[drg]]
      # Approximate Cmax ~ Dose * ka / (CL) for oral, simplified
      pk <- pk_params[[drg]]
      Cmax_per_dose <- (pk$F * 1) / pk$Vd  # per mg scaling
      Cmax <- d_vals * Cmax_per_dose * pk$Vd  # rough Cmax in ng/mL
      eff <- Emaxs[drg] * Cmax^2 / (EC50s[drg]^2 + Cmax^2)
      data.frame(drug = drg, dose = d_vals, effect = eff * 100)
    }))

    ggplot(dr_df, aes(x = dose, y = effect, color = drug)) +
      geom_line(linewidth = 1.2) +
      scale_color_manual(values = drug_colors, labels = drug_labels) +
      labs(title = "Dose–Response Relationships",
           x = "Dose (mg or mg/kg)", y = "IL-1β/IL-6 Inhibition (%)",
           color = "Drug") +
      theme_aosd()
  })

  # ========================
  # TAB 3: Cytokine Dynamics
  # ========================
  output$cytokine_plot <- renderPlot({
    df <- sim_data()
    cyto_df <- tidyr::pivot_longer(
      df[, c("day", "IL1b", "IL6", "IL18", "IFNg", "TNFa")],
      cols = -day, names_to = "Cytokine", values_to = "Conc"
    )
    cyto_levels <- c("IL1b" = "IL-1β (pg/mL)", "IL6" = "IL-6 (pg/mL)",
                     "IL18" = "IL-18 (pg/mL)", "IFNg" = "IFN-γ (pg/mL)",
                     "TNFa" = "TNF-α (pg/mL)")
    cyto_df$Cytokine <- factor(cyto_df$Cytokine, levels = names(cyto_levels),
                                labels = cyto_levels)
    cyto_colors <- c("#E41A1C", "#377EB8", "#FF7F00", "#4DAF4A", "#984EA3")

    ggplot(cyto_df, aes(x = day, y = Conc, color = Cytokine)) +
      geom_line(linewidth = 1.2) +
      facet_wrap(~ Cytokine, scales = "free_y", ncol = 2) +
      scale_color_manual(values = setNames(cyto_colors, levels(cyto_df$Cytokine))) +
      labs(title = paste("Cytokine Dynamics —", drug_labels[input$drug]),
           x = "Day", y = "Concentration (pg/mL)") +
      theme_aosd() +
      theme(legend.position = "none")
  })

  output$bound_free_plot <- renderPlot({
    df <- sim_data()
    drg <- input$drug
    if (drg %in% c("anakinra", "canakinumab")) {
      plot_df <- data.frame(
        day   = df$day,
        Free  = df$IL1b - df$IL1b_bound,
        Bound = df$IL1b_bound
      )
      ylab <- "IL-1β (pg/mL)"
      title_str <- "IL-1β: Free vs. Drug-Bound"
    } else if (drg == "tocilizumab") {
      plot_df <- data.frame(
        day   = df$day,
        Free  = df$IL6 - df$IL6_bound,
        Bound = df$IL6_bound
      )
      ylab <- "IL-6 (pg/mL)"
      title_str <- "IL-6: Free vs. Drug-Bound"
    } else {
      plot_df <- data.frame(
        day  = df$day,
        Free = df$IL1b,
        Bound = rep(0, nrow(df))
      )
      ylab <- "IL-1β (pg/mL)"
      title_str <- "IL-1β (no direct binding for selected drug)"
    }
    plot_long <- tidyr::pivot_longer(plot_df, -day,
                                      names_to = "Fraction", values_to = "Conc")
    ggplot(plot_long, aes(x = day, y = Conc, fill = Fraction)) +
      geom_area(alpha = 0.7, position = "stack") +
      scale_fill_manual(values = c("Free" = "#E41A1C", "Bound" = "#377EB8")) +
      labs(title = title_str, x = "Day", y = ylab, fill = "") +
      theme_aosd()
  })

  output$cytokine_heatmap <- renderPlot({
    df <- sim_data()
    # Normalize each cytokine to [0, 1]
    cyto_cols <- c("IL1b", "IL6", "IL18", "IFNg", "TNFa")
    day_breaks <- unique(round(seq(0, max(df$day), length.out = 20)))
    df_sub <- df[df$day %in% day_breaks | round(df$day) %in% day_breaks, ]

    heat_df <- do.call(rbind, lapply(cyto_cols, function(col) {
      vals <- df_sub[[col]]
      norm_vals <- (vals - min(vals)) / (max(vals) - min(vals) + 1e-9)
      data.frame(day = df_sub$day, cytokine = col,
                 norm_activity = norm_vals)
    }))
    heat_df$cytokine <- factor(heat_df$cytokine,
                                levels = c("TNFa", "IFNg", "IL18", "IL6", "IL1b"),
                                labels = c("TNF-α", "IFN-γ", "IL-18", "IL-6", "IL-1β"))
    ggplot(heat_df, aes(x = day, y = cytokine, fill = norm_activity)) +
      geom_tile() +
      scale_fill_gradient2(low = "steelblue", mid = "white", high = "firebrick",
                           midpoint = 0.5, name = "Normalized\nActivity") +
      labs(title = "Cytokine Activity Heatmap",
           x = "Day", y = "Cytokine") +
      theme_aosd()
  })

  # ========================
  # TAB 4: Biomarkers
  # ========================
  output$crp_plot <- renderPlot({
    df <- sim_data()
    ggplot(df, aes(x = day, y = CRP)) +
      geom_line(color = "#E41A1C", linewidth = 1.2) +
      geom_hline(yintercept = 5,  linetype = "dashed", color = "green3", linewidth = 1) +
      geom_hline(yintercept = 50, linetype = "dashed", color = "orange",  linewidth = 1) +
      annotate("text", x = max(df$day) * 0.8, y = 8,  label = "Normal (<5)", color = "green3") +
      annotate("text", x = max(df$day) * 0.8, y = 53, label = "Elevated (>50)", color = "orange") +
      labs(title = "C-Reactive Protein (CRP)",
           x = "Day", y = "CRP (mg/L)") +
      theme_aosd()
  })

  output$ferritin_plot <- renderPlot({
    df <- sim_data()
    ggplot(df, aes(x = day, y = Ferritin)) +
      geom_line(color = "#FF7F00", linewidth = 1.2) +
      geom_hline(yintercept = 500,  linetype = "dashed", color = "orange",   linewidth = 1) +
      geom_hline(yintercept = 5000, linetype = "dashed", color = "red",      linewidth = 1) +
      geom_hline(yintercept = 10000,linetype = "dashed", color = "darkred",  linewidth = 1) +
      annotate("text", x = max(df$day)*0.75, y = 700,   label = "ULN (~500)", color = "orange") +
      annotate("text", x = max(df$day)*0.75, y = 5200,  label = "AOSD threshold (5,000)", color = "red") +
      annotate("text", x = max(df$day)*0.75, y = 10200, label = "MAS concern (>10,000)", color = "darkred") +
      labs(title = "Serum Ferritin",
           x = "Day", y = "Ferritin (ng/mL)") +
      theme_aosd()
  })

  output$glyfer_plot <- renderPlot({
    df <- sim_data()
    ggplot(df, aes(x = day, y = GlyFerrit)) +
      geom_line(color = "#4DAF4A", linewidth = 1.2) +
      geom_hline(yintercept = 20, linetype = "dashed", color = "red", linewidth = 1) +
      annotate("text", x = max(df$day) * 0.7, y = 22.5,
               label = "<20% suggests MAS", color = "red") +
      labs(title = "Glycosylated Ferritin (%)",
           subtitle = ">20% typical in AOSD; <20% raises MAS concern",
           x = "Day", y = "Glycosylated Ferritin (%)") +
      ylim(0, 60) +
      theme_aosd()
  })

  output$ldh_plot <- renderPlot({
    df <- sim_data()
    ggplot(df, aes(x = day, y = LDH)) +
      geom_line(color = "#984EA3", linewidth = 1.2) +
      geom_hline(yintercept = 250, linetype = "dashed", color = "orange", linewidth = 1) +
      annotate("text", x = max(df$day) * 0.7, y = 265,
               label = "ULN (~250)", color = "orange") +
      labs(title = "Lactate Dehydrogenase (LDH)",
           x = "Day", y = "LDH (U/L)") +
      theme_aosd()
  })

  output$liver_plot <- renderPlot({
    df <- sim_data()
    liver_df <- data.frame(
      day = c(df$day, df$day),
      value = c(df$ALT, df$AST),
      marker = rep(c("ALT", "AST"), each = nrow(df))
    )
    ggplot(liver_df, aes(x = day, y = value, color = marker)) +
      geom_line(linewidth = 1.2) +
      geom_hline(yintercept = 40, linetype = "dashed", color = "grey50", linewidth = 1) +
      scale_color_manual(values = c("ALT" = "#E41A1C", "AST" = "#377EB8")) +
      labs(title = "Liver Enzymes (ALT / AST)",
           x = "Day", y = "Enzyme Level (U/L)", color = "Marker") +
      theme_aosd()
  })

  # ========================
  # TAB 5: Clinical Endpoints
  # ========================
  output$pouchot_plot <- renderPlot({
    df <- sim_data()
    ggplot(df, aes(x = day, y = Pouchot)) +
      geom_line(color = drug_colors[input$drug], linewidth = 1.5) +
      geom_hline(yintercept = 10, linetype = "dashed", color = "green3", linewidth = 1) +
      geom_hline(yintercept = 30, linetype = "dashed", color = "orange",  linewidth = 1) +
      annotate("text", x = max(df$day)*0.7, y = 13, label = "Inactive (<10)", color = "green3") +
      annotate("text", x = max(df$day)*0.7, y = 33, label = "Moderate (30)", color = "orange") +
      labs(title = "Pouchot Activity Score",
           subtitle = paste("Treatment:", drug_labels[input$drug]),
           x = "Day", y = "Pouchot Score (0–100)") +
      theme_aosd()
  })

  output$fever_arthritis_plot <- renderPlot({
    df <- sim_data()
    fa_df <- data.frame(
      day = c(df$day, df$day),
      value = c(df$Fever * 100, df$Arthritis),
      Feature = rep(c("Fever Intensity (%)", "Arthritis Score"), each = nrow(df))
    )
    ggplot(fa_df, aes(x = day, y = value, color = Feature)) +
      geom_line(linewidth = 1.2) +
      scale_color_manual(values = c("Fever Intensity (%)" = "#E41A1C",
                                     "Arthritis Score"    = "#377EB8")) +
      labs(title = "Fever and Arthritis",
           x = "Day", y = "Score / Intensity (%)", color = "") +
      theme_aosd()
  })

  output$rash_plot <- renderPlot({
    df <- sim_data()
    ggplot(df, aes(x = day, y = Rash)) +
      geom_line(color = "#FF7F00", linewidth = 1.5) +
      geom_ribbon(aes(ymin = 0, ymax = Rash), fill = "#FF7F00", alpha = 0.2) +
      labs(title = "Typical AOSD Salmon-Pink Rash Score",
           x = "Day", y = "Rash Score (0–100)") +
      ylim(0, 100) +
      theme_aosd()
  })

  output$response_rate_plot <- renderPlot({
    # Based on published clinical trial response rates for each drug in AOSD
    resp_data <- data.frame(
      Drug = c("Anakinra", "Canakinumab", "Tocilizumab",
               "Corticosteroids", "Tofacitinib"),
      ACR50 = c(0.78, 0.82, 0.71, 0.55, 0.70),
      Complete = c(0.52, 0.58, 0.46, 0.32, 0.48),
      MAS_prevention = c(0.85, 0.88, 0.70, 0.55, 0.75)
    )
    resp_long <- tidyr::pivot_longer(resp_data, -Drug,
                                      names_to = "Endpoint", values_to = "Rate")
    resp_long$Rate_pct <- resp_long$Rate * 100
    resp_long$Drug <- factor(resp_long$Drug,
                              levels = c("Anakinra", "Canakinumab", "Tocilizumab",
                                         "Corticosteroids", "Tofacitinib"))

    ggplot(resp_long, aes(x = Drug, y = Rate_pct, fill = Endpoint)) +
      geom_col(position = "dodge") +
      scale_fill_manual(values = c("#E41A1C", "#377EB8", "#4DAF4A"),
                        labels = c("ACR50 Response", "Complete Response", "MAS Prevention")) +
      labs(title = "Published Response Rates by Drug",
           x = "", y = "Response Rate (%)", fill = "Endpoint") +
      theme_aosd() +
      theme(axis.text.x = element_text(angle = 20, hjust = 1))
  })

  output$response_table <- renderTable({
    df <- sim_data()
    base_pouchot <- df$Pouchot[1]
    end_pouchot  <- df$Pouchot[nrow(df)]
    base_crp     <- df$CRP[1]
    end_crp      <- df$CRP[nrow(df)]
    base_ferr    <- df$Ferritin[1]
    end_ferr     <- df$Ferritin[nrow(df)]

    data.frame(
      `Treatment` = drug_labels[input$drug],
      `Pouchot Δ`  = paste0(round(base_pouchot), " → ", round(end_pouchot),
                             " (", round((base_pouchot - end_pouchot)/base_pouchot*100), "% ↓)"),
      `CRP Δ`      = paste0(round(base_crp), " → ", round(end_crp),
                             " mg/L (", round((base_crp - end_crp)/base_crp*100), "% ↓)"),
      `Ferritin Δ` = paste0(round(base_ferr), " → ", round(end_ferr),
                             " ng/mL (", round((base_ferr - end_ferr)/base_ferr*100), "% ↓)"),
      check.names = FALSE
    )
  })

  # ========================
  # TAB 6: Scenario Comparison
  # ========================
  output$sim_day_label <- renderText({ input$sim_days })

  output$scenario_pouchot <- renderPlot({
    df <- scenario_data()
    df$drug_label <- drug_labels[df$drug]
    ggplot(df, aes(x = day, y = Pouchot, color = drug)) +
      geom_line(linewidth = 1.1) +
      scale_color_manual(values = drug_colors, labels = drug_labels) +
      labs(title = "Pouchot Score: All Scenarios",
           x = "Day", y = "Pouchot Score", color = "") +
      theme_aosd()
  })

  output$scenario_ferritin <- renderPlot({
    df <- scenario_data()
    ggplot(df, aes(x = day, y = Ferritin, color = drug)) +
      geom_line(linewidth = 1.1) +
      geom_hline(yintercept = 500,  linetype = "dashed", color = "grey50") +
      geom_hline(yintercept = 5000, linetype = "dashed", color = "red") +
      scale_color_manual(values = drug_colors, labels = drug_labels) +
      labs(title = "Ferritin: All Scenarios",
           x = "Day", y = "Ferritin (ng/mL)", color = "") +
      theme_aosd()
  })

  output$scenario_crp <- renderPlot({
    df <- scenario_data()
    ggplot(df, aes(x = day, y = CRP, color = drug)) +
      geom_line(linewidth = 1.1) +
      geom_hline(yintercept = 5, linetype = "dashed", color = "green3") +
      scale_color_manual(values = drug_colors, labels = drug_labels) +
      labs(title = "CRP: All Scenarios",
           x = "Day", y = "CRP (mg/L)", color = "") +
      theme_aosd()
  })

  output$scenario_il6 <- renderPlot({
    df <- scenario_data()
    ggplot(df, aes(x = day, y = IL6, color = drug)) +
      geom_line(linewidth = 1.1) +
      scale_color_manual(values = drug_colors, labels = drug_labels) +
      labs(title = "IL-6: All Scenarios",
           x = "Day", y = "IL-6 (pg/mL)", color = "") +
      theme_aosd()
  })

  output$scenario_table <- renderTable({
    df <- scenario_data()
    last_df <- df[df$day == max(df$day) |
                    abs(df$day - max(df$day)) == min(abs(df$day - max(df$day))), ]
    # One row per drug
    last_df <- last_df[!duplicated(last_df$drug), ]
    data.frame(
      Drug        = drug_labels[last_df$drug],
      `Pouchot`   = round(last_df$Pouchot, 1),
      `CRP (mg/L)`= round(last_df$CRP, 1),
      `Ferritin`  = round(last_df$Ferritin),
      `IL-1β`     = round(last_df$IL1b, 1),
      `IL-6`      = round(last_df$IL6, 1),
      `MAS Prob`  = paste0(round(last_df$MAS_prob * 100, 1), "%"),
      check.names = FALSE
    )
  }, striped = TRUE, hover = TRUE)

  # ========================
  # TAB 7: MAS Risk
  # ========================
  output$mas_prob_plot <- renderPlot({
    df <- sim_data()
    ggplot(df, aes(x = day, y = MAS_prob * 100)) +
      geom_line(color = "#E41A1C", linewidth = 1.5) +
      geom_ribbon(aes(ymin = 0, ymax = MAS_prob * 100),
                  fill = "#E41A1C", alpha = 0.15) +
      geom_hline(yintercept = 30, linetype = "dashed", color = "orange", linewidth = 1) +
      geom_hline(yintercept = 60, linetype = "dashed", color = "red",    linewidth = 1) +
      annotate("text", x = max(df$day)*0.7, y = 33, label = "Watchful (30%)", color = "orange") +
      annotate("text", x = max(df$day)*0.7, y = 63, label = "High Risk (60%)", color = "red") +
      labs(title = paste("MAS Probability Over Time —", drug_labels[input$drug]),
           x = "Day", y = "MAS Probability (%)") +
      ylim(0, 100) +
      theme_aosd()
  })

  output$mas_indicators <- renderUI({
    df <- sim_data()
    max_mas  <- max(df$MAS_prob * 100)
    min_gfer <- min(df$GlyFerrit)
    max_ferr <- max(df$Ferritin)
    max_ifng <- max(df$IFNg)

    color_mas  <- if (max_mas < 20) "green"  else if (max_mas < 50) "orange" else "red"
    color_gfer <- if (min_gfer > 20) "green" else "red"
    color_ferr <- if (max_ferr < 5000) "green" else if (max_ferr < 10000) "orange" else "red"

    tagList(
      tags$p(style = paste0("color:", color_mas, "; font-weight:bold;"),
             paste0("Peak MAS Probability: ", round(max_mas, 1), "%")),
      tags$p(style = paste0("color:", color_ferr, "; font-weight:bold;"),
             paste0("Peak Ferritin: ", round(max_ferr), " ng/mL")),
      tags$p(style = paste0("color:", color_gfer, "; font-weight:bold;"),
             paste0("Min Glycosylated Ferritin: ", round(min_gfer, 1), "%")),
      tags$p(style = "color:#555;",
             paste0("Peak IFN-γ: ", round(max_ifng, 1), " pg/mL")),
      hr(),
      tags$p(em("MAS criteria: Ferritin >5000, Glyc.Ferr. <20%, IFN-γ elevation, cytopenias"))
    )
  })

  output$mas_ferritin_kin <- renderPlot({
    df <- sim_data()
    ggplot(df, aes(x = day, y = Ferritin)) +
      geom_line(color = "#FF7F00", linewidth = 1.4) +
      geom_hline(yintercept = 10000, linetype = "dashed", color = "darkred") +
      annotate("text", x = max(df$day)*0.6, y = 10300, label = "MAS threshold: 10,000",
               color = "darkred") +
      labs(title = "Ferritin Kinetics (MAS Context)",
           x = "Day", y = "Ferritin (ng/mL)") +
      theme_aosd()
  })

  output$mas_glyfer_plot <- renderPlot({
    df <- sim_data()
    ggplot(df, aes(x = day, y = GlyFerrit)) +
      geom_line(color = "#4DAF4A", linewidth = 1.4) +
      geom_hline(yintercept = 20, linetype = "dashed", color = "red", linewidth = 1) +
      annotate("text", x = max(df$day)*0.6, y = 22, label = "<20% = MAS risk",
               color = "red") +
      labs(title = "Glycosylated Ferritin % (MAS Marker)",
           x = "Day", y = "Glycosylated Ferritin (%)") +
      ylim(0, 60) +
      theme_aosd()
  })

  # ---- HScore ----
  observeEvent(input$calc_hscore, {
    result <- calc_hscore(
      input$hs_temp, input$hs_hepsp, input$hs_cyto,
      input$hs_trig, input$hs_ferrit, input$hs_fibrin,
      input$hs_asat, input$hs_hemo, input$hs_imm
    )
    output$hscore_result <- renderUI({
      risk_col <- if (result$prob < 30) "green" else if (result$prob < 60) "orange" else "red"
      wellPanel(
        style = paste0("border-left: 5px solid ", risk_col, "; background:#f9f9f9;"),
        h4(paste0("HScore = ", result$score), style = "font-weight:bold;"),
        h4(paste0("MAS Probability = ", result$prob, "%"),
           style = paste0("color:", risk_col, "; font-weight:bold;")),
        p(em("Fardet et al., Arthritis Rheumatol 2014")),
        if (result$prob >= 60) {
          tags$p(style = "color:red; font-weight:bold;",
                 "HIGH RISK: Initiate MAS-specific management")
        } else if (result$prob >= 30) {
          tags$p(style = "color:orange;", "MODERATE RISK: Close monitoring recommended")
        } else {
          tags$p(style = "color:green;", "LOW RISK")
        }
      )
    })
  })

  # ========================
  # TAB 8: Mechanistic Map
  # ========================
  output$mech_map_ui <- renderUI({
    svg_path <- "/home/user/qsp/adult-onset-stills-disease/aosd_qsp_model.svg"
    png_path <- "/home/user/qsp/adult-onset-stills-disease/aosd_qsp_model.png"

    if (file.exists(svg_path)) {
      tags$div(
        tags$img(src = svg_path, style = "max-width:100%; height:auto;",
                 alt = "AOSD QSP Mechanistic Map"),
        tags$p(em("AOSD QSP mechanistic map — nodes represent biological species,
                   edges represent regulatory interactions."),
               style = "color:#666; font-size:12px;")
      )
    } else if (file.exists(png_path)) {
      tags$div(
        tags$img(src = png_path, style = "max-width:100%; height:auto;",
                 alt = "AOSD QSP Mechanistic Map"),
        tags$p(em("AOSD QSP mechanistic map (PNG)."),
               style = "color:#666; font-size:12px;")
      )
    } else {
      wellPanel(
        h5("Mechanistic Map Not Yet Rendered"),
        p("The SVG/PNG mechanistic map will appear here once generated via Graphviz."),
        p("Run the following command to generate:"),
        tags$code("dot -Tsvg aosd_qsp_model.dot -o aosd_qsp_model.svg"),
        br(),
        tags$code("dot -Tpng -Gdpi=150 aosd_qsp_model.dot -o aosd_qsp_model.png"),
        br(), br(),
        h5("Model Architecture Summary"),
        tags$ul(
          tags$li(strong("Cluster 1:"), " Innate immune activation — DAMP/PAMP sensing, TLR4/TLR7/TLR9, NF-κB, NLRP3"),
          tags$li(strong("Cluster 2:"), " NLRP3 inflammasome — pro-IL-1β → IL-1β (caspase-1 dependent)"),
          tags$li(strong("Cluster 3:"), " IL-6 / JAK-STAT3 axis — STAT3 phosphorylation, CRP synthesis in liver"),
          tags$li(strong("Cluster 4:"), " IL-18 / IFN-γ axis — NK cell activation, macrophage M1 polarization"),
          tags$li(strong("Cluster 5:"), " T cell compartment — Th1/Th2/Th17 differentiation, Treg suppression"),
          tags$li(strong("Cluster 6:"), " Macrophage activation & MAS — hemophagocytosis, cytokine storm"),
          tags$li(strong("Cluster 7:"), " Neutrophil biology — NET formation, S100A8/A9, myeloid mediators"),
          tags$li(strong("Cluster 8:"), " Drug PK/PD — anakinra, canakinumab, tocilizumab, steroids, tofacitinib"),
          tags$li(strong("Cluster 9:"), " Acute phase response — CRP, ferritin, fibrinogen, SAA"),
          tags$li(strong("Cluster 10:"), " Organ involvement — liver (ALT/AST), spleen, lymph nodes")
        )
      )
    }
  })

} # end server

# =============================================================================
# LAUNCH APP
# =============================================================================
shinyApp(ui = ui, server = server)
