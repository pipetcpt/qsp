# =============================================================================
# NSCLC QSP Model Dashboard — Shiny Application
# Non-Small Cell Lung Cancer Quantitative Systems Pharmacology
# =============================================================================
# Run with: shiny::runApp("nsclc_shiny_app.R")
# All data embedded; no external files required.
# =============================================================================

library(shiny)
library(shinydashboard)
library(ggplot2)
library(dplyr)
library(tidyr)
library(DT)
library(plotly)
library(scales)

# =============================================================================
# SECTION 1: STANDALONE SIMULATION FUNCTIONS (defined before server)
# =============================================================================

# -----------------------------------------------------------------------------
# PK Simulation: 2-compartment oral / 1-compartment IV
# Returns data.frame with columns: time_h, conc, cycle
# -----------------------------------------------------------------------------
simulate_pk <- function(drug, dose_mg, freq_h, n_cycles, weight_kg,
                        egfr_renal, child_pugh) {

  # Drug-specific PK parameter library
  pk_params <- list(
    "Osimertinib" = list(
      model    = "oral_2cmt",
      ka       = 0.50,   # h^-1 absorption rate constant
      CL       = 14.3,   # L/h apparent clearance
      V1       = 986,    # L central volume
      Q        = 20.0,   # L/h inter-compartmental clearance
      V2       = 1000,   # L peripheral volume
      F        = 0.70,   # oral bioavailability
      typical_dose = 80
    ),
    "Alectinib" = list(
      model    = "oral_2cmt",
      ka       = 0.30,
      CL       = 12.0,
      V1       = 475,
      Q        = 15.0,
      V2       = 500,
      F        = 0.37,
      typical_dose = 600
    ),
    "Sotorasib" = list(
      model    = "oral_2cmt",
      ka       = 1.20,
      CL       = 18.0,
      V1       = 90,
      Q        = 10.0,
      V2       = 100,
      F        = 0.90,
      typical_dose = 960
    ),
    "Pembrolizumab" = list(
      model    = "iv_1cmt",
      CL       = 0.22,
      V1       = 3.5,
      Q        = 0,
      V2       = 0,
      F        = 1.0,
      typical_dose = 200
    ),
    "Cisplatin" = list(
      model    = "iv_1cmt",
      CL       = 20.0,
      V1       = 16.0,
      Q        = 0,
      V2       = 0,
      F        = 1.0,
      typical_dose = 75
    ),
    "Pemetrexed" = list(
      model    = "iv_1cmt",
      CL       = 4.5,
      V1       = 9.0,
      Q        = 0,
      V2       = 0,
      F        = 1.0,
      typical_dose = 500
    ),
    "Carboplatin" = list(
      model    = "iv_1cmt",
      CL       = 5.0,
      V1       = 23.0,
      Q        = 0,
      V2       = 0,
      F        = 1.0,
      typical_dose = 450
    )
  )

  p <- pk_params[[drug]]
  if (is.null(p)) return(NULL)

  # Apply organ-function-based dose adjustments
  CL_adj <- p$CL
  if (egfr_renal < 60) CL_adj <- CL_adj * 0.75   # moderate/severe renal impairment
  if (child_pugh == "B") CL_adj <- CL_adj * 0.70  # moderate hepatic impairment
  if (child_pugh == "C") CL_adj <- CL_adj * 0.50  # severe hepatic impairment

  dt    <- 0.25  # hours — Euler step size
  t_end <- freq_h  # simulate one full dosing interval per cycle
  times <- seq(0, t_end, by = dt)

  all_data <- list()

  for (cyc in seq_len(n_cycles)) {

    if (p$model == "oral_2cmt") {
      # Three-state Euler: [Depot (A1), Central (A2), Peripheral (A3)]
      Dose_adj <- dose_mg * p$F
      state    <- c(Dose_adj, 0, 0)
      conc_vec <- numeric(length(times))

      for (i in seq_along(times)) {
        conc_vec[i] <- state[2] / p$V1  # concentration in central compartment

        dA1 <- -p$ka * state[1]
        dA2 <-  p$ka * state[1] - (CL_adj / p$V1) * state[2] -
                (p$Q / p$V1) * state[2] + (p$Q / p$V2) * state[3]
        dA3 <- (p$Q / p$V1) * state[2] - (p$Q / p$V2) * state[3]

        state <- state + c(dA1, dA2, dA3) * dt
        state <- pmax(state, 0)
      }

    } else {
      # IV 1-compartment: instantaneous bolus → mono-exponential decay
      Dose_adj <- dose_mg
      conc_vec <- numeric(length(times))
      C0       <- Dose_adj / p$V1

      for (i in seq_along(times)) {
        conc_vec[i] <- C0 * exp(-(CL_adj / p$V1) * times[i])
      }
    }

    # Approximate steady-state accumulation from prior cycles
    if (cyc > 1) {
      C_trough  <- tail(conc_vec, 1)
      conc_vec  <- conc_vec + C_trough * 0.5
    }

    time_global <- times + (cyc - 1) * freq_h

    all_data[[cyc]] <- data.frame(
      time_h = time_global,
      conc   = conc_vec,
      cycle  = cyc
    )
  }

  do.call(rbind, all_data)
}

# Compute PK metrics (Cmax, Cmin, AUC, t1/2) from a single-cycle time series
compute_pk_metrics <- function(times, concs) {
  Cmax <- max(concs,  na.rm = TRUE)
  Cmin <- min(concs,  na.rm = TRUE)

  # Linear trapezoidal AUC
  n   <- length(times)
  auc <- sum((concs[-1] + concs[-n]) / 2 * diff(times))

  # Terminal half-life from log-linear regression of latter 50% of interval
  idx_term <- which(times > (max(times) * 0.5))
  if (length(idx_term) > 3 && all(concs[idx_term] > 0)) {
    fit   <- tryCatch(
      lm(log(concs[idx_term]) ~ times[idx_term]),
      error = function(e) NULL
    )
    slope  <- if (!is.null(fit)) coef(fit)[2] else -0.05
    t_half <- if (!is.na(slope) && slope < 0) log(2) / (-slope) else NA_real_
  } else {
    t_half <- NA_real_
  }

  list(Cmax = Cmax, Cmin = Cmin, AUC = auc, t_half = t_half)
}

# -----------------------------------------------------------------------------
# Tumor Dynamics ODE (Euler): S+R two-population model
#
#   dS/dt = lambda * S - delta * D * S   (drug-sensitive cells)
#   dR/dt = kappa * S  + lambda * R      (drug-resistant cells)
#   D = 1 (treatment on) or 0 (untreated)
#
# Returns data.frame: time_days, time_months, S_sld, R_sld, total_sld
# -----------------------------------------------------------------------------
simulate_tumor <- function(baseline_sld_mm, lambda, delta, kappa,
                           sim_months, treatment_active = TRUE) {
  dt    <- 0.1   # days
  t_end <- sim_months * 30
  times <- seq(0, t_end, by = dt)
  n     <- length(times)

  S <- baseline_sld_mm * 0.9  # 90% sensitive at baseline
  R <- baseline_sld_mm * 0.1  # 10% pre-existing resistance

  S_vec <- numeric(n)
  R_vec <- numeric(n)

  for (i in seq_len(n)) {
    S_vec[i] <- S
    R_vec[i] <- R

    D  <- if (treatment_active) 1 else 0
    dS <- lambda * S - delta * D * S
    dR <- kappa * S  + lambda * R

    S <- S + dS * dt
    R <- R + dR * dt
    S <- max(S, 0)
    R <- max(R, 0)
  }

  data.frame(
    time_days   = times,
    time_months = times / 30,
    S_sld       = S_vec,
    R_sld       = R_vec,
    total_sld   = S_vec + R_vec
  )
}

# Exponential PFS survival curve: S(t) = exp(-lambda_pfs * t)
simulate_pfs <- function(median_pfs_months, max_months = 36) {
  lambda_pfs <- log(2) / median_pfs_months
  t_months   <- seq(0, max_months, by = 0.5)
  pfs_prob   <- exp(-lambda_pfs * t_months)
  data.frame(time_months = t_months, pfs_prob = pfs_prob)
}

# CEA biomarker simulation with treatment-induced decay and resistance rebound
simulate_cea <- function(drug, sim_months) {
  times   <- seq(0, sim_months, by = 0.1)
  CEA0    <- 85   # ng/mL baseline
  decay_r <- switch(drug,
    "No Treatment"       = -0.005,  # slight rise without treatment
    "Osimertinib"        =  0.200,
    "Pembrolizumab mono" =  0.150,
    "Platinum Doublet"   =  0.120
  )

  set.seed(42)
  cea_vals <- sapply(times, function(t) {
    base_val <- CEA0 * exp(-decay_r * t)
    # Resistance-mediated rebound after month 6 for treated patients
    rebound  <- if (drug != "No Treatment" && t > 6) 0.8 * exp(0.05 * (t - 6)) else 0
    noise    <- rnorm(1, 0, 2)
    max(base_val + rebound + noise, 0.5)
  })
  data.frame(time_months = times, CEA = cea_vals)
}

# ctDNA allele frequency simulation with rebound after ~month 5
simulate_ctdna <- function(drug, sim_months) {
  times   <- seq(0, sim_months, by = 0.1)
  AF0     <- 5.0   # percent
  decay_r <- switch(drug,
    "No Treatment"       = 0.00,
    "Osimertinib"        = 0.25,
    "Pembrolizumab mono" = 0.18,
    "Platinum Doublet"   = 0.15
  )

  set.seed(99)
  af_vals <- sapply(times, function(t) {
    base_val <- AF0 * exp(-decay_r * t) + 0.2
    rebound  <- if (drug != "No Treatment" && t > 5) 0.3 * exp(0.04 * (t - 5)) else 0
    noise    <- rnorm(1, 0, 0.15)
    max(base_val + rebound + noise, 0.1)
  })
  data.frame(time_months = times, AF = af_vals)
}

# ANC simulation: 21-day chemotherapy cycle with nadir-recovery pattern
simulate_anc <- function(drug, n_cycles_bm) {
  chemo_flag <- grepl("Chemo|Platinum|Cisplatin|Pemetrexed|Carboplatin", drug)
  ANC_base   <- 4.5   # 10^9/L normal baseline
  cycle_days <- 21
  total_days <- n_cycles_bm * cycle_days

  times    <- seq(0, total_days, by = 1)
  set.seed(7)
  anc_vals <- sapply(times, function(t) {
    day_in_cycle <- t %% cycle_days
    if (chemo_flag) {
      if (day_in_cycle <= 14) {
        # Decline to nadir around day 10-14
        nadir_fraction <- 0.30
        phase          <- day_in_cycle / 14
        anc_drop       <- ANC_base * (1 - nadir_fraction) * sin(pi * phase)
        anc            <- ANC_base - anc_drop
      } else {
        # Recovery phase (day 14-21)
        recovery_phase <- (day_in_cycle - 14) / 7
        anc <- ANC_base * (0.30 + 0.70 * min(recovery_phase, 1))
      }
    } else {
      anc <- ANC_base + rnorm(1, 0, 0.2)
    }
    max(anc + rnorm(1, 0, 0.1), 0.1)
  })
  data.frame(time_days = times, ANC = anc_vals)
}

# -----------------------------------------------------------------------------
# Drug-specific toxicity profile tables
# -----------------------------------------------------------------------------
get_toxicity_profile <- function(drug) {
  profiles <- list(

    "Osimertinib" = data.frame(
      AE           = c("Rash/acneiform dermatitis","Diarrhea","Paronychia",
                       "Stomatitis","Interstitial lung disease",
                       "QTc prolongation","Decreased appetite","Fatigue"),
      Grade_1_2    = c("66%","42%","31%","20%","2%","1%","18%","21%"),
      Grade_3_plus = c("2%","1%","1%","0%","3%","<1%","1%","<1%"),
      Management   = c("Topical steroids/tetracycline","Loperamide",
                       "Antibiotics/emollients","Supportive care",
                       "Drug hold/systemic steroids","ECG monitoring",
                       "Nutritional support","Supportive care"),
      stringsAsFactors = FALSE
    ),

    "Alectinib" = data.frame(
      AE           = c("Constipation","Fatigue","Myalgia","Peripheral edema",
                       "Elevated CPK","Bradycardia","Photosensitivity","Nausea"),
      Grade_1_2    = c("33%","26%","24%","23%","11%","8%","9%","14%"),
      Grade_3_plus = c("0%","2%","0%","1%","5%","<1%","0%","0%"),
      Management   = c("Stool softeners","Rest","Analgesics",
                       "Elevation/diuretics","Monitor CPK; hold if G3",
                       "ECG; adjust concurrent drugs","Sunscreen/protective clothing",
                       "Antiemetics"),
      stringsAsFactors = FALSE
    ),

    "Sotorasib" = data.frame(
      AE           = c("Diarrhea","Nausea","Elevated AST/ALT","Fatigue",
                       "Musculoskeletal pain","Cough","Constipation","Rash"),
      Grade_1_2    = c("46%","26%","17%","23%","15%","14%","12%","9%"),
      Grade_3_plus = c("11%","4%","15%","4%","1%","1%","1%","1%"),
      Management   = c("Loperamide; hold for G3","Antiemetics","LFT monitoring; hold for G3",
                       "Supportive","Analgesics","Cough suppressants",
                       "Stool softeners","Topical steroids"),
      stringsAsFactors = FALSE
    ),

    "Pembrolizumab mono" = data.frame(
      AE           = c("Fatigue","Rash","Diarrhea","Pruritus",
                       "Hypothyroidism","Pneumonitis","Immune hepatitis","Immune colitis"),
      Grade_1_2    = c("18%","13%","11%","10%","8%","4%","2%","2%"),
      Grade_3_plus = c("2%","1%","1%","0%","0%","3%","1%","1%"),
      Management   = c("Supportive","Topical/systemic steroids","Supportive/antidiarrheals",
                       "Antihistamines","Levothyroxine","High-dose corticosteroids",
                       "Hold + prednisone 0.5-1 mg/kg","Hold + prednisone 1-2 mg/kg"),
      stringsAsFactors = FALSE
    ),

    "Pembrolizumab+Chemo" = data.frame(
      AE           = c("Nausea","Anemia","Neutropenia","Fatigue",
                       "Alopecia","Peripheral neuropathy","Pneumonitis","Thrombocytopenia"),
      Grade_1_2    = c("43%","66%","44%","34%","38%","10%","5%","35%"),
      Grade_3_plus = c("5%","19%","28%","6%","3%","2%","4%","8%"),
      Management   = c("Antiemetics","EPO/transfusion","G-CSF/dose reduction",
                       "Supportive","Supportive","Gabapentin/B6 supplementation",
                       "Corticosteroids; hold ICI","Platelet support/dose delay"),
      stringsAsFactors = FALSE
    ),

    "Platinum Doublet" = data.frame(
      AE           = c("Nausea/vomiting","Neutropenia","Anemia","Thrombocytopenia",
                       "Renal toxicity","Peripheral neuropathy","Ototoxicity","Fatigue"),
      Grade_1_2    = c("55%","48%","65%","38%","10%","12%","5%","30%"),
      Grade_3_plus = c("12%","25%","15%","8%","5%","3%","2%","8%"),
      Management   = c("Antiemetics + aggressive hydration","G-CSF; delay cycle",
                       "EPO/transfusion","Delay/reduce dose",
                       "Pre/post hydration; amifostine","Gabapentin; consider carboplatin switch",
                       "Audiometry; dose reduce","Supportive care"),
      stringsAsFactors = FALSE
    )
  )
  profiles[[drug]]
}

# Logistic approximation of irAE probabilities based on PD-L1 TPS
compute_irae_probs <- function(pdl1_tps) {
  base_risk <- 0.10 + (pdl1_tps / 100) * 0.15
  data.frame(
    irAE = c("Pneumonitis","Colitis","Hepatitis","Thyroiditis",
             "Hypophysitis","Adrenal insufficiency","Nephritis","Dermatitis"),
    Grade_1_2_prob = round(c(0.80, 0.60, 0.50, 1.20, 0.30, 0.40, 0.30, 0.90) *
                             base_risk * 100, 1),
    Grade_3_4_prob = round(c(0.25, 0.15, 0.10, 0.05, 0.05, 0.05, 0.05, 0.10) *
                             base_risk * 100, 1),
    stringsAsFactors = FALSE
  )
}

# =============================================================================
# SECTION 2: UI
# =============================================================================

ui <- dashboardPage(
  skin = "blue",

  # ---------------------------------------------------------------------------
  # Header
  # ---------------------------------------------------------------------------
  dashboardHeader(
    title = "NSCLC QSP Dashboard",
    titleWidth = 260
  ),

  # ---------------------------------------------------------------------------
  # Sidebar
  # ---------------------------------------------------------------------------
  dashboardSidebar(
    width = 260,
    sidebarMenu(
      id = "sidebar",
      menuItem("Patient Profile",     tabName = "patient",    icon = icon("user")),
      menuItem("Pharmacokinetics",    tabName = "pk",         icon = icon("chart-line")),
      menuItem("Tumor Dynamics",      tabName = "tumor",      icon = icon("circle-dot")),
      menuItem("Biomarkers",          tabName = "biomarkers", icon = icon("vial")),
      menuItem("Scenario Comparison", tabName = "scenarios",  icon = icon("sliders")),
      menuItem("Toxicity",            tabName = "toxicity",   icon = icon("triangle-exclamation"))
    ),
    hr(),
    div(style = "padding: 10px; color: #ccc; font-size: 11px;",
        "NSCLC QSP Model v1.0",
        br(), "For research use only",
        br(), "Parameters from published trials"
    )
  ),

  # ---------------------------------------------------------------------------
  # Body
  # ---------------------------------------------------------------------------
  dashboardBody(

    # Custom CSS for traffic lights, callout boxes, etc.
    tags$head(tags$style(HTML("
      .content-wrapper { background-color: #f4f6f9; }
      .box-header { border-bottom: 2px solid #3c8dbc; }
      .traffic-light-green  {
        background: #27ae60; color: white; padding: 7px 12px;
        border-radius: 5px; font-weight: bold; text-align: center;
        margin: 3px 0; display: block;
      }
      .traffic-light-yellow {
        background: #f39c12; color: white; padding: 7px 12px;
        border-radius: 5px; font-weight: bold; text-align: center;
        margin: 3px 0; display: block;
      }
      .traffic-light-red {
        background: #c0392b; color: white; padding: 7px 12px;
        border-radius: 5px; font-weight: bold; text-align: center;
        margin: 3px 0; display: block;
      }
      .rec-box {
        background: #eaf4fb; border-left: 4px solid #3c8dbc;
        padding: 12px 16px; border-radius: 4px;
        margin-top: 10px; font-size: 13px;
      }
      .warn-box {
        background: #fef9e7; border-left: 4px solid #f39c12;
        padding: 12px 16px; border-radius: 4px;
        margin-top: 10px; font-size: 13px;
      }
    "))),

    tabItems(

      # -----------------------------------------------------------------------
      # TAB 1: Patient Profile
      # -----------------------------------------------------------------------
      tabItem(tabName = "patient",
        fluidRow(
          box(title = "Patient Demographics & Molecular Profile",
              width = 4, status = "primary", solidHeader = TRUE,
              numericInput("age", "Age (years)", value = 62, min = 18, max = 90),
              selectInput("ecog", "ECOG Performance Status",
                          choices = c("0","1","2","3"), selected = "1"),
              selectInput("stage", "Disease Stage",
                          choices = c("IIIA","IIIB","IVA","IVB"), selected = "IVA"),
              selectInput("histology", "Histology",
                          choices = c("Adenocarcinoma","Squamous Cell",
                                      "Large Cell NEC")),
              hr(),
              h5("Molecular Biomarkers"),
              selectInput("egfr_mut", "EGFR Mutation",
                          choices = c("None","Exon 19 del","L858R","T790M")),
              checkboxInput("alk_fusion", "ALK Fusion Positive",  value = FALSE),
              checkboxInput("kras_g12c",  "KRAS G12C Positive",  value = FALSE)
          ),
          box(title = "Tumor & Immune Parameters",
              width = 4, status = "primary", solidHeader = TRUE,
              sliderInput("pdl1_tps", "PD-L1 TPS (%)",
                          min = 0, max = 100, value = 55, step = 1),
              sliderInput("tmb", "TMB (mut/Mb)",
                          min = 0, max = 50, value = 10, step = 1),
              selectInput("smoking", "Smoking Status",
                          choices = c("Never","Former","Current")),
              numericInput("pack_years", "Pack-Years (if applicable)",
                           value = 30, min = 0, max = 200),
              numericInput("baseline_tumor",
                           "Baseline Tumor Size — SLD (mm)",
                           value = 45, min = 5, max = 200),
              helpText("SLD = sum of longest diameters per RECIST 1.1")
          ),
          box(title = "Clinical Summary Indicators",
              width = 4, status = "info", solidHeader = TRUE,
              valueBoxOutput("vbox_ecog",  width = 12),
              valueBoxOutput("vbox_pdl1",  width = 12),
              valueBoxOutput("vbox_tmb",   width = 12)
          )
        ),
        fluidRow(
          box(title = "Patient Characteristics Summary",
              width = 6, status = "primary", solidHeader = TRUE,
              DTOutput("patient_table")
          ),
          box(title = "Treatment Eligibility Matrix",
              width = 6, status = "warning", solidHeader = TRUE,
              DTOutput("eligibility_table")
          )
        ),
        fluidRow(
          box(title = "Baseline Tumor Assessment (RECIST Waterfall — Baseline)",
              width = 12, status = "primary", solidHeader = TRUE,
              plotlyOutput("baseline_waterfall", height = "250px")
          )
        )
      ),

      # -----------------------------------------------------------------------
      # TAB 2: Pharmacokinetics
      # -----------------------------------------------------------------------
      tabItem(tabName = "pk",
        fluidRow(
          box(title = "Drug & Dosing Parameters",
              width = 3, status = "primary", solidHeader = TRUE,
              selectInput("drug", "Drug",
                          choices = c("Osimertinib","Alectinib","Sotorasib",
                                      "Pembrolizumab","Cisplatin",
                                      "Pemetrexed","Carboplatin"),
                          selected = "Osimertinib"),
              numericInput("dose_mg", "Dose (mg)",
                           value = 80, min = 10, max = 1200),
              selectInput("dose_freq", "Dosing Frequency",
                          choices = c("Once daily","Twice daily","Q3W","Q21D"),
                          selected = "Once daily"),
              hr(),
              h5("Patient Organ Function"),
              numericInput("body_weight", "Body Weight (kg)",
                           value = 70, min = 30, max = 150),
              numericInput("egfr_renal",
                           "eGFR (mL/min/1.73m²)",
                           value = 80, min = 10, max = 150),
              selectInput("child_pugh", "Child-Pugh Class",
                          choices = c("A","B","C"), selected = "A"),
              numericInput("n_cycles_pk", "Number of Cycles",
                           value = 4, min = 1, max = 12),
              actionButton("sim_pk", "Simulate PK",
                           class = "btn-primary btn-block",
                           icon  = icon("play"))
          ),
          box(title = "Concentration-Time Profile (Multi-Cycle)",
              width = 9, status = "primary", solidHeader = TRUE,
              plotlyOutput("pk_plot", height = "380px"),
              div(class = "rec-box", uiOutput("pk_dose_recommendation"))
          )
        ),
        fluidRow(
          box(title = "PK Parameter Summary by Cycle",
              width = 12, status = "primary", solidHeader = TRUE,
              DTOutput("pk_table")
          )
        )
      ),

      # -----------------------------------------------------------------------
      # TAB 3: Tumor Dynamics
      # -----------------------------------------------------------------------
      tabItem(tabName = "tumor",
        fluidRow(
          box(title = "Tumor Growth-Kill Model Parameters",
              width = 3, status = "primary", solidHeader = TRUE,
              selectInput("tumor_drug", "Treatment Regimen",
                          choices = c("No Treatment","Osimertinib","Alectinib",
                                      "Sotorasib","Pembrolizumab mono",
                                      "Pembrolizumab+Chemo","Platinum Doublet"),
                          selected = "Osimertinib"),
              sliderInput("sim_months", "Simulation Duration (months)",
                          min = 3, max = 24, value = 12, step = 1),
              numericInput("tumor_baseline_sz", "Baseline SLD (mm)",
                           value = 45, min = 5, max = 200),
              hr(),
              h5("ODE Model Parameters"),
              sliderInput("growth_rate", "Tumor Growth Rate λ",
                          min = 0.001, max = 0.05, value = 0.01, step = 0.001),
              sliderInput("kill_rate", "Drug Kill Rate δ",
                          min = 0.000, max = 0.10, value = 0.04, step = 0.001),
              sliderInput("resistance_rate", "Resistance Emergence κ",
                          min = 0.001, max = 0.02, value = 0.005, step = 0.001),
              actionButton("sim_tumor", "Simulate Tumor",
                           class = "btn-primary btn-block",
                           icon  = icon("play"))
          ),
          box(title = "Tumor SLD Dynamics Over Time",
              width = 9, status = "primary", solidHeader = TRUE,
              plotlyOutput("tumor_timecourse", height = "360px")
          )
        ),
        fluidRow(
          box(title = "RECIST Waterfall — Best Response by Timepoint",
              width = 6, status = "success", solidHeader = TRUE,
              plotlyOutput("recist_waterfall", height = "300px")
          ),
          box(title = "Spider Plot — Parameter Sensitivity (±20%)",
              width = 6, status = "info", solidHeader = TRUE,
              plotlyOutput("spider_plot", height = "300px")
          )
        ),
        fluidRow(
          box(title = "Response Metrics Summary",
              width = 12, status = "primary", solidHeader = TRUE,
              DTOutput("tumor_summary_table")
          )
        )
      ),

      # -----------------------------------------------------------------------
      # TAB 4: Biomarkers
      # -----------------------------------------------------------------------
      tabItem(tabName = "biomarkers",
        fluidRow(
          box(title = "Biomarker Simulation Settings",
              width = 3, status = "primary", solidHeader = TRUE,
              selectInput("bm_drug", "Treatment",
                          choices = c("No Treatment","Osimertinib",
                                      "Pembrolizumab mono","Platinum Doublet"),
                          selected = "Osimertinib"),
              sliderInput("bm_months", "Simulation Duration (months)",
                          min = 3, max = 18, value = 12, step = 1),
              actionButton("sim_bm", "Simulate Biomarkers",
                           class = "btn-primary btn-block",
                           icon  = icon("play")),
              hr(),
              div(class = "rec-box",
                  h5("Clinical Reference PFS"),
                  p(strong("Osimertinib:"), " 18.9 mo (FLAURA)"),
                  p(strong("Pembrolizumab:"), " 10.3 mo (KEYNOTE-024)"),
                  p(strong("Platinum doublet:"), " 5.7 mo (historical)"),
                  p(strong("No treatment:"), " 4.0 mo")
              )
          ),
          box(title = "Four-Panel Biomarker Trajectories",
              width = 9, status = "primary", solidHeader = TRUE,
              plotlyOutput("bm_panel_plot", height = "520px")
          )
        ),
        fluidRow(
          box(title = "PFS Probability at Key Timepoints",
              width = 6, status = "info", solidHeader = TRUE,
              DTOutput("pfs_table")
          ),
          box(title = "Resistance Emergence Interpretation",
              width = 6, status = "warning", solidHeader = TRUE,
              uiOutput("resistance_text")
          )
        )
      ),

      # -----------------------------------------------------------------------
      # TAB 5: Scenario Comparison
      # -----------------------------------------------------------------------
      tabItem(tabName = "scenarios",
        fluidRow(
          box(title = "Regimen Selection",
              width = 3, status = "primary", solidHeader = TRUE,
              checkboxGroupInput("scenarios", "Select Regimens to Compare",
                choices  = c("No Treatment","Osimertinib","Alectinib",
                             "Sotorasib","Pembrolizumab mono",
                             "Pembrolizumab+Chemo","Platinum Doublet"),
                selected = c("No Treatment","Osimertinib","Pembrolizumab mono")
              ),
              hr(),
              sliderInput("comp_months", "Comparison Duration (months)",
                          min = 6, max = 24, value = 18, step = 1),
              numericInput("comp_baseline", "Baseline SLD (mm)",
                           value = 45, min = 5, max = 200),
              actionButton("sim_scenarios", "Compare Scenarios",
                           class = "btn-primary btn-block",
                           icon  = icon("play"))
          ),
          box(title = "Tumor Volume Comparison — All Regimens",
              width = 9, status = "primary", solidHeader = TRUE,
              plotlyOutput("scenario_timecourse", height = "380px")
          )
        ),
        fluidRow(
          box(title = "Efficacy Summary Table",
              width = 7, status = "primary", solidHeader = TRUE,
              DTOutput("scenario_table")
          ),
          box(title = "Forest Plot — Hazard Ratio vs No Treatment",
              width = 5, status = "info", solidHeader = TRUE,
              plotlyOutput("forest_plot", height = "360px")
          )
        )
      ),

      # -----------------------------------------------------------------------
      # TAB 6: Toxicity
      # -----------------------------------------------------------------------
      tabItem(tabName = "toxicity",
        fluidRow(
          box(title = "Toxicity Assessment Parameters",
              width = 3, status = "primary", solidHeader = TRUE,
              selectInput("tox_drug", "Treatment",
                          choices = c("Osimertinib","Alectinib","Sotorasib",
                                      "Pembrolizumab mono","Pembrolizumab+Chemo",
                                      "Platinum Doublet"),
                          selected = "Pembrolizumab+Chemo"),
              numericInput("tox_cycles", "Number of Cycles",
                           value = 6, min = 1, max = 12),
              hr(),
              h5("Patient Risk Factors"),
              sliderInput("rash_grade", "Observed Rash Grade (EGFR TKI)",
                          min = 0, max = 4, value = 1, step = 1),
              sliderInput("diarrhea_grade", "Observed Diarrhea Grade",
                          min = 0, max = 4, value = 1, step = 1),
              sliderInput("pdl1_for_tox", "PD-L1 TPS for irAE Risk (%)",
                          min = 0, max = 100, value = 55, step = 1),
              actionButton("sim_tox", "Assess Toxicity",
                           class = "btn-danger btn-block",
                           icon  = icon("play"))
          ),
          box(title = "Hematologic Toxicity Over Cycles",
              width = 9, status = "danger", solidHeader = TRUE,
              plotlyOutput("hematox_plot", height = "360px")
          )
        ),
        fluidRow(
          box(title = "Drug-Specific Toxicity Profile",
              width = 7, status = "warning", solidHeader = TRUE,
              DTOutput("tox_table")
          ),
          box(title = "Organ System Traffic Light",
              width = 5, status = "primary", solidHeader = TRUE,
              uiOutput("traffic_light")
          )
        ),
        fluidRow(
          box(title = "irAE Probability Estimates (by PD-L1 TPS)",
              width = 6, status = "warning", solidHeader = TRUE,
              DTOutput("irae_table")
          ),
          box(title = "Dose Modification Recommendations",
              width = 6, status = "info", solidHeader = TRUE,
              uiOutput("dose_mod_text")
          )
        )
      )  # end tabItem toxicity

    )  # end tabItems
  )  # end dashboardBody
)  # end dashboardPage

# =============================================================================
# SECTION 3: SERVER
# =============================================================================

server <- function(input, output, session) {

  # ---------------------------------------------------------------------------
  # Auto-update dose defaults when drug changes
  # ---------------------------------------------------------------------------
  observeEvent(input$drug, {
    defaults <- list(
      "Osimertinib"   = list(dose = 80,  freq = "Once daily"),
      "Alectinib"     = list(dose = 600, freq = "Twice daily"),
      "Sotorasib"     = list(dose = 960, freq = "Once daily"),
      "Pembrolizumab" = list(dose = 200, freq = "Q3W"),
      "Cisplatin"     = list(dose = 75,  freq = "Q3W"),
      "Pemetrexed"    = list(dose = 500, freq = "Q3W"),
      "Carboplatin"   = list(dose = 450, freq = "Q3W")
    )
    d <- defaults[[input$drug]]
    if (!is.null(d)) {
      updateNumericInput(session, "dose_mg",   value    = d$dose)
      updateSelectInput(session,  "dose_freq", selected = d$freq)
    }
  })

  # Auto-suggest kill rate when tumor_drug changes
  drug_kill_rates <- c(
    "No Treatment"        = 0.000,
    "Osimertinib"         = 0.055,
    "Alectinib"           = 0.060,
    "Sotorasib"           = 0.045,
    "Pembrolizumab mono"  = 0.040,
    "Pembrolizumab+Chemo" = 0.048,
    "Platinum Doublet"    = 0.032
  )
  observeEvent(input$tumor_drug, {
    kr <- drug_kill_rates[input$tumor_drug]
    if (!is.na(kr)) updateSliderInput(session, "kill_rate", value = kr)
  })

  # ---------------------------------------------------------------------------
  # TAB 1 outputs
  # ---------------------------------------------------------------------------

  output$vbox_ecog <- renderValueBox({
    ecog_val   <- as.integer(input$ecog)
    subtitles  <- c("Fully active","Ambulatory, restricted strenuous",
                    "Self-care only","Limited self-care")
    clr        <- c("green","yellow","orange","red")[ecog_val + 1]
    clr        <- if (clr == "orange") "red" else clr
    valueBox(
      value    = paste0("ECOG PS: ", input$ecog),
      subtitle = subtitles[ecog_val + 1],
      icon     = icon("person-walking"),
      color    = clr
    )
  })

  output$vbox_pdl1 <- renderValueBox({
    tps    <- input$pdl1_tps
    cat_lbl <- if (tps >= 50) "High (≥50%)" else
               if (tps >= 1)  "Low (1–49%)" else "Negative (<1%)"
    clr    <- if (tps >= 50) "green" else if (tps >= 1) "yellow" else "red"
    valueBox(
      value    = paste0("PD-L1: ", tps, "%"),
      subtitle = cat_lbl,
      icon     = icon("microscope"),
      color    = clr
    )
  })

  output$vbox_tmb <- renderValueBox({
    tmb_val <- input$tmb
    cat_lbl <- if (tmb_val >= 10) "TMB-High (≥10)" else "TMB-Low (<10)"
    clr     <- if (tmb_val >= 10) "green" else "yellow"
    valueBox(
      value    = paste0("TMB: ", tmb_val, " mut/Mb"),
      subtitle = cat_lbl,
      icon     = icon("dna"),
      color    = clr
    )
  })

  output$patient_table <- renderDT({
    df <- data.frame(
      Parameter = c("Age","ECOG PS","Stage","Histology","EGFR Mutation",
                    "ALK Fusion","KRAS G12C","PD-L1 TPS","TMB",
                    "Smoking Status","Pack-Years","Baseline Tumor SLD"),
      Value = c(
        paste0(input$age, " yr"),
        paste0("PS ", input$ecog),
        input$stage,
        input$histology,
        input$egfr_mut,
        ifelse(input$alk_fusion, "Positive", "Negative"),
        ifelse(input$kras_g12c,  "Positive", "Negative"),
        paste0(input$pdl1_tps, "%"),
        paste0(input$tmb, " mut/Mb"),
        input$smoking,
        input$pack_years,
        paste0(input$baseline_tumor, " mm")
      ),
      stringsAsFactors = FALSE
    )
    datatable(df,
              options  = list(dom = "t", pageLength = 15),
              rownames = FALSE,
              class    = "compact stripe")
  })

  output$eligibility_table <- renderDT({
    egfr_pos  <- input$egfr_mut %in% c("Exon 19 del","L858R","T790M")
    alk_pos   <- input$alk_fusion
    kras_pos  <- input$kras_g12c
    pdl1_high <- input$pdl1_tps >= 50
    sq_cell   <- input$histology == "Squamous Cell"
    no_driver <- !egfr_pos && !alk_pos && !kras_pos

    df <- data.frame(
      Drug = c("Osimertinib","Alectinib","Sotorasib",
               "Pembrolizumab mono","Pembrolizumab+Chemo",
               "Platinum Doublet","Bevacizumab add-on"),
      Biomarker_Requirement = c(
        "EGFR mut+ (Ex19del/L858R/T790M)",
        "ALK fusion+",
        "KRAS G12C+",
        "PD-L1 TPS ≥50%",
        "Any (no EGFR/ALK driver preferred)",
        "No restriction (standard 2nd line)",
        "Non-squamous histology"
      ),
      Eligible = c(
        ifelse(egfr_pos,    "YES", "NO"),
        ifelse(alk_pos,     "YES", "NO"),
        ifelse(kras_pos,    "YES", "NO"),
        ifelse(pdl1_high,   "YES", "Conditional"),
        ifelse(no_driver,   "YES", "Consider"),
        "YES",
        ifelse(!sq_cell, "YES", "NO")
      ),
      Line = c("1L / 2L (T790M)","1L","2L+","1L","1L","1L / 2L","1L add-on"),
      stringsAsFactors = FALSE
    )

    datatable(df,
              options  = list(dom = "t", pageLength = 10),
              rownames = FALSE,
              class    = "compact stripe") |>
      formatStyle("Eligible",
        backgroundColor = styleEqual(
          c("YES","NO","Conditional","Consider"),
          c("#d4efdf","#fadbd8","#fef9e7","#d6eaf8")
        )
      )
  })

  output$baseline_waterfall <- renderPlotly({
    tryCatch({
      df <- data.frame(
        patient = "Current Patient",
        pct_chg = 0,
        status  = "Baseline"
      )
      p <- ggplot(df, aes(x = patient, y = pct_chg, fill = status)) +
        geom_col(width = 0.25) +
        geom_hline(yintercept = -30, linetype = "dashed", color = "#27ae60", linewidth = 0.7) +
        geom_hline(yintercept =  20, linetype = "dashed", color = "#e74c3c", linewidth = 0.7) +
        annotate("text", x = 1.35, y = -32, label = "PR threshold (-30%)",
                 color = "#27ae60", size = 3) +
        annotate("text", x = 1.35, y =  22, label = "PD threshold (+20%)",
                 color = "#e74c3c", size = 3) +
        scale_fill_manual(values = c("Baseline" = "#3498db")) +
        labs(x = NULL, y = "% Change from Baseline",
             title = paste0("Baseline SLD: ", input$baseline_tumor, " mm — Awaiting treatment response")) +
        theme_minimal(base_size = 12) +
        theme(legend.position = "none") +
        ylim(-100, 60)
      ggplotly(p, tooltip = c("x","y"))
    }, error = function(e) {
      plotly::plot_ly() |> layout(title = paste("Render error:", conditionMessage(e)))
    })
  })

  # ---------------------------------------------------------------------------
  # TAB 2: Pharmacokinetics
  # ---------------------------------------------------------------------------

  freq_to_hours <- function(freq_str) {
    switch(freq_str,
      "Once daily"  = 24,
      "Twice daily" = 12,
      "Q3W"         = 504,
      "Q21D"        = 504,
      24  # default
    )
  }

  pk_data <- eventReactive(input$sim_pk, {
    validate(
      need(input$dose_mg > 0,     "Dose must be > 0 mg"),
      need(input$n_cycles_pk > 0, "Cycles must be ≥ 1"),
      need(input$egfr_renal > 0,  "eGFR must be > 0")
    )
    freq_h <- freq_to_hours(input$dose_freq)
    simulate_pk(
      drug       = input$drug,
      dose_mg    = input$dose_mg,
      freq_h     = freq_h,
      n_cycles   = input$n_cycles_pk,
      weight_kg  = input$body_weight,
      egfr_renal = input$egfr_renal,
      child_pugh = input$child_pugh
    )
  })

  output$pk_plot <- renderPlotly({
    tryCatch({
      df <- pk_data()
      req(df)
      p <- ggplot(df, aes(x = time_h, y = conc, color = factor(cycle))) +
        geom_line(linewidth = 0.9) +
        scale_color_viridis_d(name = "Cycle", option = "D") +
        labs(x = "Time (hours)", y = "Plasma Concentration (mg/L)",
             title = paste(input$drug, "— Multi-Cycle PK Simulation")) +
        theme_minimal(base_size = 12) +
        theme(legend.position = "right")
      ggplotly(p, tooltip = c("x","y","colour"))
    }, error = function(e) {
      plotly::plot_ly() |> layout(title = paste("Error:", conditionMessage(e)))
    })
  })

  output$pk_table <- renderDT({
    tryCatch({
      df     <- pk_data()
      req(df)
      freq_h <- freq_to_hours(input$dose_freq)

      metrics_list <- lapply(unique(df$cycle), function(cyc) {
        sub_df  <- df[df$cycle == cyc, ]
        t_local <- sub_df$time_h - min(sub_df$time_h)
        m       <- compute_pk_metrics(t_local, sub_df$conc)
        data.frame(
          Cycle  = cyc,
          Cmax   = round(m$Cmax,   3),
          Cmin   = round(m$Cmin,   4),
          AUC    = round(m$AUC,    1),
          t_half = round(m$t_half, 1)
        )
      })
      out_df        <- do.call(rbind, metrics_list)
      names(out_df) <- c("Cycle", "Cmax (mg/L)", "Cmin (mg/L)",
                         "AUC (mg·h/L)", "t½ (h)")
      datatable(out_df,
                options  = list(dom = "t", pageLength = 12),
                rownames = FALSE,
                class    = "compact stripe") |>
        formatRound(columns = 2:5, digits = 3)
    }, error = function(e) {
      datatable(data.frame(Error = conditionMessage(e)))
    })
  })

  output$pk_dose_recommendation <- renderUI({
    msgs <- character(0)
    if (input$egfr_renal < 30)
      msgs <- c(msgs, "<b>Severe renal impairment (eGFR &lt;30 mL/min):</b> CL adjusted by ~50%; significant drug accumulation expected — consider dose reduction and close monitoring.")
    else if (input$egfr_renal < 60)
      msgs <- c(msgs, "<b>Moderate renal impairment (eGFR 30-59 mL/min):</b> CL reduced ~25%; dose adjustment may be warranted depending on drug; monitor for toxicity.")
    if (input$child_pugh == "B")
      msgs <- c(msgs, "<b>Child-Pugh B:</b> Hepatic CL reduced ~30%; standard dose generally used but close LFT monitoring required.")
    if (input$child_pugh == "C")
      msgs <- c(msgs, "<b>Child-Pugh C:</b> Hepatic CL reduced ~50%; dose reduction strongly recommended; avoid hepatotoxic agents if possible.")
    if (length(msgs) == 0)
      msgs <- "No dose adjustments required based on current organ function inputs. Normal renal/hepatic function assumed."
    HTML(paste0("<b>Dose Adjustment Guidance:</b><br>", paste(msgs, collapse = "<br><br>")))
  })

  # ---------------------------------------------------------------------------
  # TAB 3: Tumor Dynamics
  # ---------------------------------------------------------------------------

  tumor_data <- eventReactive(input$sim_tumor, {
    validate(
      need(input$tumor_baseline_sz > 0, "Baseline SLD must be > 0"),
      need(input$sim_months > 0,        "Simulation duration must be > 0")
    )
    treatment_on <- input$tumor_drug != "No Treatment"
    simulate_tumor(
      baseline_sld_mm  = input$tumor_baseline_sz,
      lambda           = input$growth_rate,
      delta            = input$kill_rate,
      kappa            = input$resistance_rate,
      sim_months       = input$sim_months,
      treatment_active = treatment_on
    )
  })

  output$tumor_timecourse <- renderPlotly({
    tryCatch({
      df <- tumor_data()
      req(df)

      df_long <- tidyr::pivot_longer(
        df[, c("time_months","S_sld","R_sld","total_sld")],
        cols      = c("S_sld","R_sld","total_sld"),
        names_to  = "Component",
        values_to = "SLD_mm"
      )
      df_long$Component <- factor(df_long$Component,
        levels = c("S_sld","R_sld","total_sld"),
        labels = c("Sensitive Cells","Resistant Cells","Total SLD")
      )

      baseline_val <- df$total_sld[1]

      p <- ggplot(df_long, aes(x = time_months, y = SLD_mm,
                               color = Component, linetype = Component)) +
        geom_line(linewidth = 1.1) +
        scale_color_manual(values = c("#2980b9","#e74c3c","#2c3e50")) +
        scale_linetype_manual(values = c("solid","dashed","solid")) +
        geom_hline(yintercept = baseline_val * 0.70,
                   linetype = "dotted", color = "#27ae60", alpha = 0.8) +
        geom_hline(yintercept = baseline_val * 1.20,
                   linetype = "dotted", color = "#e74c3c", alpha = 0.8) +
        annotate("text", x = max(df$time_months) * 0.85,
                 y = baseline_val * 0.70 - 2,
                 label = "PR threshold", color = "#27ae60", size = 3) +
        annotate("text", x = max(df$time_months) * 0.85,
                 y = baseline_val * 1.20 + 2,
                 label = "PD threshold", color = "#e74c3c", size = 3) +
        labs(x = "Time (months)", y = "SLD (mm)",
             title = paste("Tumor Dynamics —", input$tumor_drug),
             color = NULL, linetype = NULL) +
        theme_minimal(base_size = 12) +
        theme(legend.position = "bottom")
      ggplotly(p, tooltip = c("x","y","colour"))
    }, error = function(e) {
      plotly::plot_ly() |> layout(title = paste("Error:", conditionMessage(e)))
    })
  })

  output$recist_waterfall <- renderPlotly({
    tryCatch({
      df           <- tumor_data()
      req(df)
      baseline_sld <- df$total_sld[1]

      months_check <- c(2, 4, 6, 8, 10, 12)[c(2, 4, 6, 8, 10, 12) <= input$sim_months]
      if (length(months_check) == 0) months_check <- input$sim_months

      wf_data <- lapply(months_check, function(m) {
        t_days  <- m * 30
        idx     <- which.min(abs(df$time_days - t_days))
        pct_chg <- (df$total_sld[idx] - baseline_sld) / baseline_sld * 100
        data.frame(label = paste0("Mo ", m), pct_chg = round(pct_chg, 1))
      })
      wf_df <- do.call(rbind, wf_data)
      wf_df$color_cat <- ifelse(wf_df$pct_chg <= -30, "PR/CR",
                          ifelse(wf_df$pct_chg >= 20, "PD", "SD"))

      p <- ggplot(wf_df, aes(x = reorder(label, pct_chg),
                             y = pct_chg, fill = color_cat)) +
        geom_col() +
        geom_hline(yintercept = -30, linetype = "dashed", color = "#27ae60") +
        geom_hline(yintercept =  20, linetype = "dashed", color = "#e74c3c") +
        scale_fill_manual(values = c("PR/CR" = "#27ae60","SD" = "#f39c12","PD" = "#e74c3c")) +
        coord_flip() +
        labs(x = NULL, y = "% Change from Baseline", fill = "RECIST",
             title = "Waterfall: Best Response by Timepoint") +
        theme_minimal(base_size = 11)
      ggplotly(p, tooltip = c("x","y","fill"))
    }, error = function(e) {
      plotly::plot_ly() |> layout(title = paste("Error:", conditionMessage(e)))
    })
  })

  output$spider_plot <- renderPlotly({
    tryCatch({
      df_base      <- tumor_data()
      req(df_base)
      baseline_sld <- df_base$total_sld[1]
      treatment_on <- input$tumor_drug != "No Treatment"

      # Five perturbation scenarios: ±20% in growth and kill
      perturb_sets <- list(
        list(lam = input$growth_rate,        del = input$kill_rate,        label = "Base case"),
        list(lam = input$growth_rate * 1.2,  del = input$kill_rate,        label = "+20% Growth"),
        list(lam = input$growth_rate * 0.8,  del = input$kill_rate,        label = "-20% Growth"),
        list(lam = input$growth_rate,        del = input$kill_rate * 1.2,  label = "+20% Kill"),
        list(lam = input$growth_rate,        del = input$kill_rate * 0.8,  label = "-20% Kill")
      )

      all_spider <- lapply(perturb_sets, function(ps) {
        sim_df     <- simulate_tumor(
          baseline_sld_mm  = input$tumor_baseline_sz,
          lambda           = ps$lam,
          delta            = ps$del,
          kappa            = input$resistance_rate,
          sim_months       = input$sim_months,
          treatment_active = treatment_on
        )
        months_seq <- seq(0, input$sim_months, by = 1)
        sampled    <- lapply(months_seq, function(m) {
          idx  <- which.min(abs(sim_df$time_months - m))
          pct  <- (sim_df$total_sld[idx] - baseline_sld) / baseline_sld * 100
          data.frame(time_months = m, pct_chg = round(pct, 1), scenario = ps$label)
        })
        do.call(rbind, sampled)
      })

      spider_df <- do.call(rbind, all_spider)

      p <- ggplot(spider_df, aes(x = time_months, y = pct_chg,
                                 color = scenario, group = scenario)) +
        geom_line(linewidth = 0.9) +
        geom_hline(yintercept = -30, linetype = "dashed",
                   color = "#27ae60", alpha = 0.7) +
        geom_hline(yintercept =  20, linetype = "dashed",
                   color = "#e74c3c", alpha = 0.7) +
        scale_color_viridis_d(option = "C") +
        labs(x = "Time (months)", y = "% Change from Baseline",
             title = "Spider Plot — Parameter Sensitivity",
             color = NULL) +
        theme_minimal(base_size = 11) +
        theme(legend.position = "bottom",
              legend.text     = element_text(size = 9))
      ggplotly(p, tooltip = c("x","y","colour"))
    }, error = function(e) {
      plotly::plot_ly() |> layout(title = paste("Error:", conditionMessage(e)))
    })
  })

  output$tumor_summary_table <- renderDT({
    tryCatch({
      df           <- tumor_data()
      req(df)
      baseline_sld <- df$total_sld[1]

      get_pct_at <- function(months) {
        if (months > input$sim_months) return(NA_real_)
        idx <- which.min(abs(df$time_months - months))
        round((df$total_sld[idx] - baseline_sld) / baseline_sld * 100, 1)
      }

      min_idx   <- which.min(df$total_sld)
      bor_pct   <- round((df$total_sld[min_idx] - baseline_sld) / baseline_sld * 100, 1)
      bor_label <- if (bor_pct <= -100) "CR" else if (bor_pct <= -30) "PR" else
                   if (bor_pct >= 20) "PD" else "SD"
      time_best <- round(df$time_months[min_idx], 1)

      pd_idx  <- which(df$total_sld > baseline_sld * 1.20)
      est_pfs <- if (length(pd_idx) > 0) round(df$time_months[pd_idx[1]], 1) else
                 paste0(">", input$sim_months)

      summary_df <- data.frame(
        Metric = c(
          "Best Overall Response (BOR)",
          "% Change at Best Response",
          "Time to Best Response (months)",
          "Estimated PFS (months)",
          "% Change at Month 3",
          "% Change at Month 6",
          "% Change at Month 12"
        ),
        Value = c(
          bor_label,
          paste0(bor_pct, "%"),
          time_best,
          est_pfs,
          paste0(get_pct_at(3),  "%"),
          paste0(get_pct_at(6),  "%"),
          paste0(get_pct_at(12), "%")
        ),
        stringsAsFactors = FALSE
      )
      datatable(summary_df,
                options  = list(dom = "t", pageLength = 10),
                rownames = FALSE,
                class    = "compact stripe")
    }, error = function(e) {
      datatable(data.frame(Error = conditionMessage(e)))
    })
  })

  # ---------------------------------------------------------------------------
  # TAB 4: Biomarkers
  # ---------------------------------------------------------------------------

  bm_data <- eventReactive(input$sim_bm, {
    validate(need(input$bm_months > 0, "Duration must be > 0"))
    drug <- input$bm_drug
    dur  <- input$bm_months

    pfs_medians <- c(
      "No Treatment"       = 4.0,
      "Osimertinib"        = 18.9,
      "Pembrolizumab mono" = 10.3,
      "Platinum Doublet"   = 5.7
    )
    med_pfs <- pfs_medians[drug]

    list(
      cea     = simulate_cea(drug, dur),
      ctdna   = simulate_ctdna(drug, dur),
      anc     = simulate_anc(drug, ceiling(dur / 0.7)),
      pfs     = simulate_pfs(med_pfs, max_months = max(36, dur)),
      med_pfs = med_pfs,
      drug    = drug
    )
  })

  output$bm_panel_plot <- renderPlotly({
    tryCatch({
      bm <- bm_data()
      req(bm)

      # Panel 1: CEA
      p1 <- plot_ly(bm$cea, x = ~time_months, y = ~CEA,
                    type = "scatter", mode = "lines",
                    line = list(color = "#2980b9", width = 2),
                    name = "CEA (ng/mL)") |>
        layout(
          xaxis = list(title = "Time (months)"),
          yaxis = list(title = "CEA (ng/mL)"),
          title = list(text = "CEA Trajectory", font = list(size = 12))
        )

      # Panel 2: ctDNA
      p2 <- plot_ly(bm$ctdna, x = ~time_months, y = ~AF,
                    type = "scatter", mode = "lines",
                    line = list(color = "#e74c3c", width = 2),
                    name = "ctDNA AF (%)") |>
        layout(
          xaxis = list(title = "Time (months)"),
          yaxis = list(title = "ctDNA AF (%)"),
          title = list(text = "ctDNA Allele Frequency", font = list(size = 12))
        )

      # Panel 3: ANC
      bm$anc$time_months_anc <- bm$anc$time_days / 30
      anc_sub <- bm$anc[bm$anc$time_months_anc <= input$bm_months, ]
      p3 <- plot_ly(anc_sub, x = ~time_months_anc, y = ~ANC,
                    type = "scatter", mode = "lines",
                    line = list(color = "#27ae60", width = 2),
                    name = "ANC (10⁹/L)") |>
        add_segments(x = 0, xend = input$bm_months, y = 1.5, yend = 1.5,
                     line = list(color = "red", dash = "dash", width = 1),
                     showlegend = FALSE) |>
        layout(
          xaxis = list(title = "Time (months)"),
          yaxis = list(title = "ANC (10⁹/L)"),
          title = list(text = "ANC Trajectory", font = list(size = 12))
        )

      # Panel 4: PFS Kaplan-Meier style
      pfs_sub <- bm$pfs[bm$pfs$time_months <= max(36, input$bm_months), ]
      p4 <- plot_ly(pfs_sub, x = ~time_months, y = ~pfs_prob,
                    type = "scatter", mode = "lines",
                    line  = list(color = "#8e44ad", width = 2),
                    name  = "PFS Probability") |>
        add_segments(x = 0, xend = max(pfs_sub$time_months),
                     y = 0.5, yend = 0.5,
                     line = list(color = "grey50", dash = "dash", width = 1),
                     showlegend = FALSE) |>
        layout(
          xaxis = list(title = "Time (months)"),
          yaxis = list(title = "PFS Probability", range = c(0, 1)),
          title = list(text = paste("PFS Curve (median:", bm$med_pfs, "mo)"),
                       font = list(size = 12))
        )

      subplot(p1, p2, p3, p4,
              nrows   = 2,
              titleX  = TRUE,
              titleY  = TRUE,
              shareX  = FALSE,
              shareY  = FALSE,
              margin  = 0.08) |>
        layout(
          title      = paste("Biomarker Panel —", bm$drug),
          showlegend = FALSE
        )
    }, error = function(e) {
      plotly::plot_ly() |> layout(title = paste("Error:", conditionMessage(e)))
    })
  })

  output$pfs_table <- renderDT({
    tryCatch({
      bm     <- bm_data()
      req(bm)
      pfs_df <- bm$pfs

      rows <- lapply(c(6, 12, 18, 24), function(t) {
        idx  <- which.min(abs(pfs_df$time_months - t))
        prob <- round(pfs_df$pfs_prob[idx] * 100, 1)
        data.frame(
          Timepoint   = paste0(t, " months"),
          PFS_Prob    = paste0(prob, "%"),
          Median_PFS  = paste0(bm$med_pfs, " months"),
          stringsAsFactors = FALSE
        )
      })
      out_df <- do.call(rbind, rows)
      names(out_df) <- c("Timepoint", "PFS Probability", "Median PFS")
      datatable(out_df,
                options  = list(dom = "t"),
                rownames = FALSE,
                class    = "compact stripe")
    }, error = function(e) {
      datatable(data.frame(Error = conditionMessage(e)))
    })
  })

  output$resistance_text <- renderUI({
    tryCatch({
      bm   <- bm_data()
      req(bm)
      drug <- bm$drug

      content <- switch(drug,
        "Osimertinib" = div(
          h5("Resistance Mechanisms — Osimertinib"),
          p("After a median ~18 months, on-target and off-target resistance mechanisms emerge:"),
          tags$ul(
            tags$li(strong("On-target:"), " EGFR C797S (cis or trans with T790M), EGFR amplification"),
            tags$li(strong("Off-target:"), " MET amplification, HER2 amplification, KRAS mutations"),
            tags$li(strong("Histologic:"), " Small-cell transformation (~5-10% of cases)")
          ),
          p("ctDNA rebound (as modeled) typically precedes radiologic progression by 4-8 weeks."),
          div(class = "rec-box",
              strong("Recommendation:"), " Repeat liquid biopsy at progression. If C797S detected:
              consider 4th-generation EGFR TKI (BLU-945) in trials or platinum doublet.")
        ),
        "Pembrolizumab mono" = div(
          h5("Resistance Mechanisms — Pembrolizumab"),
          p("Primary resistance (~55% of PD-L1 high NSCLC):"),
          tags$ul(
            tags$li("JAK1/JAK2 loss-of-function mutations (impairs IFN-γ signaling)"),
            tags$li("Beta-2-microglobulin mutations (MHC I antigen presentation loss)"),
            tags$li("STK11/KEAP1 mutations (associated with non-response)")
          ),
          p("Acquired resistance via T-cell exhaustion; upregulation of TIM-3, LAG-3, TIGIT."),
          div(class = "rec-box",
              strong("Recommendation:"), " Consider dual ICI (pembrolizumab + ipilimumab) or switch
              to platinum doublet. Anti-LAG-3 (relatlimab) trials open for this setting.")
        ),
        "Platinum Doublet" = div(
          h5("Resistance Mechanisms — Platinum Doublet"),
          tags$ul(
            tags$li("DNA repair upregulation: ERCC1 overexpression, XPC polymorphisms"),
            tags$li("Drug efflux: ATP7B/ATP7A copper transporters"),
            tags$li("Epithelial-mesenchymal transition (EMT)"),
            tags$li("TP53 loss-of-function mutations")
          ),
          p("Median duration of platinum response 4-6 months; CEA rebound is an early indicator."),
          div(class = "rec-box",
              strong("2nd line options:"), " Docetaxel ± ramucirumab (REVEL trial);
              atezolizumab (OAK trial); nintedanib + docetaxel (adenocarcinoma).")
        ),
        div(
          h5("Natural Disease History"),
          p("Without treatment, NSCLC follows exponential growth kinetics (Gompertzian at large volumes)."),
          p("Median OS untreated stage IV NSCLC: ~4 months historically."),
          p("Spontaneous regressions are exceedingly rare (<0.1%).")
        )
      )
      content
    }, error = function(e) {
      div(class = "warn-box", paste("Error:", conditionMessage(e)))
    })
  })

  # ---------------------------------------------------------------------------
  # TAB 5: Scenario Comparison
  # ---------------------------------------------------------------------------

  # Literature-derived scenario parameters
  scenario_params <- list(
    "No Treatment"        = list(lambda = 0.012, delta = 0.000, ORR = 0,  mPFS = 4.0,  HR = 1.00),
    "Osimertinib"         = list(lambda = 0.010, delta = 0.055, ORR = 80, mPFS = 18.9, HR = 0.46),
    "Alectinib"           = list(lambda = 0.009, delta = 0.060, ORR = 83, mPFS = 25.7, HR = 0.47),
    "Sotorasib"           = list(lambda = 0.011, delta = 0.045, ORR = 37, mPFS = 6.8,  HR = 0.63),
    "Pembrolizumab mono"  = list(lambda = 0.010, delta = 0.040, ORR = 45, mPFS = 10.3, HR = 0.50),
    "Pembrolizumab+Chemo" = list(lambda = 0.010, delta = 0.048, ORR = 48, mPFS = 8.8,  HR = 0.49),
    "Platinum Doublet"    = list(lambda = 0.011, delta = 0.032, ORR = 35, mPFS = 5.7,  HR = 0.71)
  )

  scenario_colors <- c(
    "No Treatment"        = "#7f8c8d",
    "Osimertinib"         = "#2980b9",
    "Alectinib"           = "#27ae60",
    "Sotorasib"           = "#e67e22",
    "Pembrolizumab mono"  = "#8e44ad",
    "Pembrolizumab+Chemo" = "#c0392b",
    "Platinum Doublet"    = "#16a085"
  )

  scenario_data <- eventReactive(input$sim_scenarios, {
    validate(
      need(length(input$scenarios) >= 1, "Select at least one scenario"),
      need(input$comp_baseline > 0,      "Baseline SLD must be > 0")
    )

    results <- lapply(input$scenarios, function(sc) {
      p   <- scenario_params[[sc]]
      sim <- simulate_tumor(
        baseline_sld_mm  = input$comp_baseline,
        lambda           = p$lambda,
        delta            = p$delta,
        kappa            = 0.005,
        sim_months       = input$comp_months,
        treatment_active = (sc != "No Treatment")
      )
      sim$scenario <- sc
      sim
    })
    list(
      all_sims = do.call(rbind, results),
      selected = input$scenarios
    )
  })

  output$scenario_timecourse <- renderPlotly({
    tryCatch({
      dat <- scenario_data()
      req(dat)
      df  <- dat$all_sims

      p <- ggplot(df, aes(x = time_months, y = total_sld,
                          color = scenario, group = scenario)) +
        geom_line(linewidth = 1.1) +
        scale_color_manual(values = scenario_colors, name = "Regimen") +
        labs(x = "Time (months)", y = "Total SLD (mm)",
             title = "Tumor Volume Comparison — All Selected Regimens") +
        theme_minimal(base_size = 12) +
        theme(legend.position = "right")
      ggplotly(p, tooltip = c("x","y","colour"))
    }, error = function(e) {
      plotly::plot_ly() |> layout(title = paste("Error:", conditionMessage(e)))
    })
  })

  output$scenario_table <- renderDT({
    tryCatch({
      dat <- scenario_data()
      req(dat)

      ref_trials <- c(
        "No Treatment"        = "—",
        "Osimertinib"         = "FLAURA (HR=0.46)",
        "Alectinib"           = "ALEX (HR=0.47)",
        "Sotorasib"           = "CodeBreaK 100",
        "Pembrolizumab mono"  = "KEYNOTE-024 (HR=0.50)",
        "Pembrolizumab+Chemo" = "KEYNOTE-189 (HR=0.49)",
        "Platinum Doublet"    = "Historical control"
      )

      rows <- lapply(dat$selected, function(sc) {
        p   <- scenario_params[[sc]]
        dcr <- min(p$ORR + 20, 95)
        data.frame(
          Regimen   = sc,
          ORR       = paste0(p$ORR, "%"),
          DCR       = paste0(dcr, "%"),
          mPFS      = paste0(p$mPFS, " mo"),
          HR        = p$HR,
          Reference = ref_trials[sc],
          stringsAsFactors = FALSE
        )
      })
      out_df        <- do.call(rbind, rows)
      names(out_df) <- c("Regimen","ORR","DCR","Median PFS","HR vs Control","Reference Trial")

      datatable(out_df,
                options  = list(dom = "t", pageLength = 10),
                rownames = FALSE,
                class    = "compact stripe") |>
        formatStyle("HR vs Control",
          backgroundColor = styleInterval(
            c(0.50, 0.70),
            c("#d4efdf","#d6eaf8","#fadbd8")
          )
        )
    }, error = function(e) {
      datatable(data.frame(Error = conditionMessage(e)))
    })
  })

  output$forest_plot <- renderPlotly({
    tryCatch({
      dat <- scenario_data()
      req(dat)

      non_ctrl <- dat$selected[dat$selected != "No Treatment"]
      if (length(non_ctrl) == 0)
        return(plotly::plot_ly() |>
               layout(title = "Select at least one non-'No Treatment' regimen"))

      forest_df <- do.call(rbind, lapply(non_ctrl, function(sc) {
        p  <- scenario_params[[sc]]
        hr <- p$HR
        data.frame(
          drug   = sc,
          HR     = hr,
          HR_lo  = max(hr - 0.15, 0.10),
          HR_hi  = min(hr + 0.15, 0.99),
          stringsAsFactors = FALSE
        )
      }))
      forest_df$y_pos <- seq_len(nrow(forest_df))

      # Benchmark annotations from pivotal trials
      benchmarks <- data.frame(
        drug = c("FLAURA (Osimertinib)","ALEX (Alectinib)","KEYNOTE-024 (Pembro)"),
        HR   = c(0.46, 0.47, 0.50),
        y    = c(-0.5, -1.2, -1.9),
        stringsAsFactors = FALSE
      )

      clrs <- sapply(forest_df$drug, function(d) {
        v <- scenario_colors[[d]]
        if (is.null(v)) "#333" else v
      })

      p <- plot_ly() |>
        add_segments(
          data  = forest_df,
          x     = ~HR_lo, xend = ~HR_hi,
          y     = ~y_pos, yend = ~y_pos,
          line  = list(width = 3, color = "black"),
          showlegend = FALSE
        ) |>
        add_markers(
          data      = forest_df,
          x         = ~HR,
          y         = ~y_pos,
          marker    = list(size = 14, color = clrs, symbol = "diamond"),
          text      = ~paste0(drug, "<br>HR = ", HR,
                              " [", round(HR_lo, 2), "–", round(HR_hi, 2), "]"),
          hoverinfo = "text",
          showlegend = FALSE
        ) |>
        add_vline(x = 1.0,
                  line = list(color = "black", dash = "dash", width = 1)) |>
        # Benchmark reference markers
        add_markers(
          data       = benchmarks,
          x          = ~HR, y = ~y,
          marker     = list(size = 10, color = "#aaa", symbol = "star"),
          text       = ~paste0(drug, " HR=", HR),
          hoverinfo  = "text",
          showlegend = FALSE
        ) |>
        layout(
          title  = "Forest Plot — HR vs No Treatment (diamonds = simulation; stars = trials)",
          xaxis  = list(title = "Hazard Ratio (95% CI)", range = c(0.1, 1.3)),
          yaxis  = list(
            title    = "",
            tickvals = forest_df$y_pos,
            ticktext = forest_df$drug,
            showgrid = FALSE
          ),
          margin = list(l = 180)
        )
      p
    }, error = function(e) {
      plotly::plot_ly() |> layout(title = paste("Error:", conditionMessage(e)))
    })
  })

  # ---------------------------------------------------------------------------
  # TAB 6: Toxicity
  # ---------------------------------------------------------------------------

  tox_sim_data <- eventReactive(input$sim_tox, {
    validate(need(input$tox_cycles > 0, "Cycles must be ≥ 1"))
    drug       <- input$tox_drug
    n_cyc      <- input$tox_cycles
    chemo_flag <- grepl("Chemo|Platinum", drug)
    cycle_days <- 21
    total_days <- n_cyc * cycle_days
    times      <- seq(0, total_days, by = 1)

    # ANC
    anc_full <- simulate_anc(drug, n_cyc)
    anc_df   <- data.frame(
      time_days = times,
      ANC       = anc_full$ANC[seq_along(times)]
    )

    # Platelet count simulation
    set.seed(17)
    plt_vals <- sapply(times, function(t) {
      day_in_cycle <- t %% cycle_days
      if (chemo_flag) {
        plt_base <- 220
        if (day_in_cycle <= 16) {
          phase <- day_in_cycle / 16
          drop  <- plt_base * 0.4 * sin(pi * phase)
          val   <- plt_base - drop
        } else {
          rp  <- min((day_in_cycle - 16) / 5, 1)
          val <- 220 * (0.60 + 0.40 * rp)
        }
      } else {
        val <- 220
      }
      max(val + rnorm(1, 0, 5), 10)
    })

    # Hemoglobin simulation — cumulative decline with chemo
    set.seed(23)
    hgb_vals <- sapply(times, function(t) {
      cycle_num <- floor(t / cycle_days) + 1
      if (chemo_flag) {
        base_hgb <- 13.5 - (cycle_num - 1) * 0.25
        day_drop <- sin(pi * (t %% cycle_days) / cycle_days) * 0.5
        val      <- base_hgb - day_drop
      } else {
        val <- 13.5
      }
      max(val + rnorm(1, 0, 0.1), 5)
    })

    list(
      anc_df    = anc_df,
      plt_df    = data.frame(time_days = times, PLT = plt_vals),
      hgb_df    = data.frame(time_days = times, HGB = hgb_vals),
      tox_table = get_toxicity_profile(drug),
      irae_tbl  = compute_irae_probs(input$pdl1_for_tox),
      drug      = drug,
      n_cyc     = n_cyc
    )
  })

  output$hematox_plot <- renderPlotly({
    tryCatch({
      td <- tox_sim_data()
      req(td)

      # ANC panel
      p1 <- plot_ly(td$anc_df, x = ~time_days, y = ~ANC,
                    type = "scatter", mode = "lines",
                    line = list(color = "#2980b9", width = 1.5),
                    name = "ANC") |>
        add_segments(x = 0, xend = max(td$anc_df$time_days),
                     y = 1.5, yend = 1.5,
                     line = list(color = "#e74c3c", dash = "dash", width = 1),
                     showlegend = FALSE) |>
        add_segments(x = 0, xend = max(td$anc_df$time_days),
                     y = 0.5, yend = 0.5,
                     line = list(color = "darkred", dash = "dot", width = 1),
                     showlegend = FALSE) |>
        layout(xaxis = list(title = "Day"),
               yaxis = list(title = "ANC (10⁹/L)"),
               title = list(text = "ANC", font = list(size = 12)))

      # Platelet panel
      p2 <- plot_ly(td$plt_df, x = ~time_days, y = ~PLT,
                    type = "scatter", mode = "lines",
                    line = list(color = "#e74c3c", width = 1.5),
                    name = "Platelets") |>
        add_segments(x = 0, xend = max(td$plt_df$time_days),
                     y = 100, yend = 100,
                     line = list(color = "#e74c3c", dash = "dash", width = 1),
                     showlegend = FALSE) |>
        add_segments(x = 0, xend = max(td$plt_df$time_days),
                     y = 50, yend = 50,
                     line = list(color = "darkred", dash = "dot", width = 1),
                     showlegend = FALSE) |>
        layout(xaxis = list(title = "Day"),
               yaxis = list(title = "PLT (10⁹/L)"),
               title = list(text = "Platelet Count", font = list(size = 12)))

      # Hemoglobin panel
      p3 <- plot_ly(td$hgb_df, x = ~time_days, y = ~HGB,
                    type = "scatter", mode = "lines",
                    line = list(color = "#27ae60", width = 1.5),
                    name = "Hgb") |>
        add_segments(x = 0, xend = max(td$hgb_df$time_days),
                     y = 10, yend = 10,
                     line = list(color = "#f39c12", dash = "dash", width = 1),
                     showlegend = FALSE) |>
        add_segments(x = 0, xend = max(td$hgb_df$time_days),
                     y = 8, yend = 8,
                     line = list(color = "#e74c3c", dash = "dot", width = 1),
                     showlegend = FALSE) |>
        layout(xaxis = list(title = "Day"),
               yaxis = list(title = "Hgb (g/dL)"),
               title = list(text = "Hemoglobin", font = list(size = 12)))

      subplot(p1, p2, p3,
              nrows  = 1,
              titleX = TRUE,
              titleY = TRUE,
              shareX = FALSE,
              shareY = FALSE) |>
        layout(
          title      = paste("Hematologic Toxicity —", td$drug,
                             "(", td$n_cyc, "cycles)"),
          showlegend = FALSE,
          annotations = list(list(
            x = 0.5, y = 1.07, xref = "paper", yref = "paper",
            text = "Dashed = G3 threshold | Dotted = G4 threshold",
            showarrow = FALSE, font = list(size = 10, color = "#666")
          ))
        )
    }, error = function(e) {
      plotly::plot_ly() |> layout(title = paste("Error:", conditionMessage(e)))
    })
  })

  output$tox_table <- renderDT({
    tryCatch({
      td <- tox_sim_data()
      req(td)
      df <- td$tox_table
      names(df) <- c("Adverse Event", "G1-2 Frequency", "G3+ Frequency", "Management")
      datatable(df,
                options  = list(dom = "t", pageLength = 15, scrollX = TRUE),
                rownames = FALSE,
                class    = "compact stripe")
    }, error = function(e) {
      datatable(data.frame(Error = conditionMessage(e)))
    })
  })

  output$traffic_light <- renderUI({
    tryCatch({
      td   <- tox_sim_data()
      req(td)
      drug       <- td$drug
      chemo_flag <- grepl("Chemo|Platinum", drug)
      egfr_tki   <- drug %in% c("Osimertinib","Alectinib","Sotorasib")
      ici_flag   <- grepl("Pembrolizumab", drug)

      # Assign expected toxicity grade per organ system
      organs <- data.frame(
        System = c("Hematologic","Gastrointestinal","Skin",
                   "Pulmonary","Hepatic","Renal","Neurologic"),
        Grade  = c(
          if (chemo_flag) "G3-4" else if (egfr_tki || ici_flag) "G0-1" else "G0-1",
          if (chemo_flag || egfr_tki) "G2" else if (ici_flag) "G2" else "G0-1",
          if (egfr_tki) "G2" else if (ici_flag) "G2" else "G0-1",
          if (ici_flag) "G2" else "G0-1",
          if (ici_flag || drug == "Sotorasib") "G2" else if (chemo_flag) "G2" else "G0-1",
          if (grepl("Cisplatin|Platinum", drug)) "G2" else if (chemo_flag) "G2" else "G0-1",
          if (chemo_flag) "G2" else "G0-1"
        ),
        stringsAsFactors = FALSE
      )

      grade_class <- function(g) {
        switch(g,
          "G0-1" = "traffic-light-green",
          "G2"   = "traffic-light-yellow",
          "G3-4" = "traffic-light-red",
          "traffic-light-green"
        )
      }

      boxes <- lapply(seq_len(nrow(organs)), function(i) {
        span(class = grade_class(organs$Grade[i]),
             paste0(organs$System[i], ": ", organs$Grade[i]))
      })

      div(
        h5("Organ System Toxicity Risk"),
        p(style = "font-size: 11px; color: #555; margin-bottom: 8px;",
          "Based on drug class and typical toxicity profiles from registration trials."),
        do.call(tagList, boxes),
        br(),
        div(style = "font-size: 10px; color: #777;",
            span(class = "traffic-light-green",  style = "display:inline-block;padding:2px 6px;font-size:10px;", "G0-1"),
            " = Low risk  ",
            span(class = "traffic-light-yellow", style = "display:inline-block;padding:2px 6px;font-size:10px;", "G2"),
            " = Moderate risk  ",
            span(class = "traffic-light-red",    style = "display:inline-block;padding:2px 6px;font-size:10px;", "G3-4"),
            " = High risk"
        )
      )
    }, error = function(e) {
      div(class = "warn-box", paste("Error:", conditionMessage(e)))
    })
  })

  output$irae_table <- renderDT({
    tryCatch({
      td <- tox_sim_data()
      req(td)
      df <- td$irae_tbl
      names(df) <- c("irAE", "G1-2 Probability (%)", "G3-4 Probability (%)")
      datatable(df,
                options  = list(dom = "t", pageLength = 10),
                rownames = FALSE,
                class    = "compact stripe",
                caption  = paste0("Estimated irAE risk at PD-L1 TPS ",
                                  input$pdl1_for_tox, "%")) |>
        formatStyle("G3-4 Probability (%)",
          backgroundColor = styleInterval(c(3, 8),
                                          c("#d4efdf","#fef9e7","#fadbd8")))
    }, error = function(e) {
      datatable(data.frame(Error = conditionMessage(e)))
    })
  })

  output$dose_mod_text <- renderUI({
    tryCatch({
      td   <- tox_sim_data()
      req(td)
      drug <- td$drug

      content <- switch(drug,
        "Osimertinib" = div(class = "rec-box",
          h5("Osimertinib Dose Modifications (TAGRISSO PI)"),
          tags$ul(
            tags$li(strong("Rash G3:"), " Withhold; restart 40 mg QD if ≤G2 within 3 wks; if not, discontinue"),
            tags$li(strong("ILD G1:"),  " Consider interruption; close monitoring required"),
            tags$li(strong("ILD ≥G2:"), " Permanently discontinue"),
            tags$li(strong("QTcF >500 ms:"), " Withhold; restart 40 mg if QTcF ≤480 ms"),
            tags$li(strong("Diarrhea G3+:"), " Withhold; restart 40 mg if ≤G1 within 3 wks"),
            tags$li(strong("Keratitis:"), " Hold; ophthalmology referral; restart 40 mg if resolved")
          )
        ),
        "Alectinib" = div(class = "rec-box",
          h5("Alectinib Dose Modifications (ALECENSA PI)"),
          tags$ul(
            tags$li(strong("CPK ≥5×ULN or myopathy G3:"), " Hold; restart 450 mg BID when ≤G1"),
            tags$li(strong("Bradycardia (symptom.):"), " Hold; adjust concurrent rate-lowering drugs; restart 450 mg BID"),
            tags$li(strong("Hepatotoxicity G3 (≥5×ULN):"), " Hold; restart 450 mg BID if ≤G1 within 3 wks"),
            tags$li(strong("Hepatotoxicity G4:"), " Permanently discontinue"),
            tags$li(strong("Photosensitivity G3+:"), " Hold; restart 450 mg BID; counsel sun avoidance")
          )
        ),
        "Sotorasib" = div(class = "rec-box",
          h5("Sotorasib Dose Modifications (LUMAKRAS PI)"),
          tags$ul(
            tags$li(strong("Hepatotoxicity G3 (>5×ULN):"), " Hold; reduce to 480 mg QD if resolved ≤G1"),
            tags$li(strong("Hepatotoxicity G4:"), " Permanently discontinue"),
            tags$li(strong("Diarrhea G3+:"), " Hold; restart 480 mg if ≤G1; avoid with proton pump inhibitors (reduce absorption)"),
            tags$li(strong("ILD/Pneumonitis any grade:"), " Permanently discontinue"),
            tags$li(strong("2nd dose reduction:"), " 240 mg QD; 3rd = discontinue")
          )
        ),
        "Pembrolizumab mono" = div(class = "rec-box",
          h5("Pembrolizumab irAE Dose Modifications (KEYTRUDA PI)"),
          tags$ul(
            tags$li(strong("Pneumonitis G2:"), " Hold; prednisone 1-2 mg/kg/day → taper over ≥4 weeks"),
            tags$li(strong("Pneumonitis G3-4:"), " Permanently discontinue; IV methylprednisolone 1-2 mg/kg/day"),
            tags$li(strong("Colitis G2:"), " Hold; consider budesonide or prednisone 1 mg/kg"),
            tags$li(strong("Colitis G3-4:"), " Permanently discontinue; IV methylprednisolone + infliximab if refractory"),
            tags$li(strong("Hepatitis G3 (>5×ULN):"), " Hold; prednisone 1-2 mg/kg + mycophenolate"),
            tags$li(strong("Hepatitis G4:"), " Permanently discontinue")
          )
        ),
        "Pembrolizumab+Chemo" = div(class = "rec-box",
          h5("Pembrolizumab+Chemo Modifications (KEYNOTE-189)"),
          tags$ul(
            tags$li(strong("Febrile neutropenia:"), " Hold chemotherapy; G-CSF; restart at 75% pemetrexed/carboplatin"),
            tags$li(strong("Thrombocytopenia G4:"), " Hold pemetrexed; restart at 75% dose if recovered"),
            tags$li(strong("Renal impairment (CrCl <45):"), " Pemetrexed contraindicated; switch regimen"),
            tags$li(strong("irAE:"), " Manage pembrolizumab as per single-agent ICI guidance above"),
            tags$li(strong("Neuropathy G3+:"), " Reduce/discontinue cisplatin; switch to carboplatin")
          )
        ),
        "Platinum Doublet" = div(class = "rec-box",
          h5("Platinum Doublet Dose Modifications"),
          tags$ul(
            tags$li(strong("Creatinine ≥1.5×ULN:"), " Hold cisplatin; consider carboplatin substitution; maintain hydration 3L/day"),
            tags$li(strong("Neuropathy G2:"), " Hold cisplatin; restart at 75% if resolved to ≤G1"),
            tags$li(strong("Neuropathy G3:"), " Permanently reduce or discontinue cisplatin; switch to carboplatin"),
            tags$li(strong("Neutropenia G4 or FN:"), " Reduce next dose 25%; prophylactic G-CSF for subsequent cycles"),
            tags$li(strong("Ototoxicity:"), " Serial audiometry; dose reduce or switch to carboplatin/oxaliplatin"),
            tags$li(strong("Pemetrexed folic acid:"), " Ensure folic acid 400-1000 mcg daily + B12 1000 mcg q9wks (mandatory)")
          )
        ),
        div(class = "warn-box", p("No specific modifications available for selected drug."))
      )
      content
    }, error = function(e) {
      div(class = "warn-box", paste("Error:", conditionMessage(e)))
    })
  })

}  # end server

# =============================================================================
# SECTION 4: LAUNCH
# =============================================================================

shiny::shinyApp(ui = ui, server = server)
