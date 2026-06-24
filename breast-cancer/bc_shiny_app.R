# =============================================================================
# Breast Cancer QSP Interactive Dashboard
# bc_shiny_app.R
#
# A fully self-contained Shiny application for quantitative systems
# pharmacology (QSP) simulation of breast cancer drug therapies.
#
# All PK/PD simulations use base-R Euler integration — mrgsolve is NOT required.
#
# Run with: shiny::runApp("bc_shiny_app.R")
# =============================================================================

library(shiny)
library(shinydashboard)
library(ggplot2)
library(dplyr)
library(plotly)

# =============================================================================
# HELPER FUNCTIONS — PK / PD / Survival simulation
# =============================================================================

# -----------------------------------------------------------------------------
# 1-compartment oral PK via Euler integration
# Returns data.frame(time, Cp, dose_event)
# -----------------------------------------------------------------------------
sim_pk_1cpt <- function(dose, ka, Vd, CL, F_oral = 1, ii = 24,
                        n_doses = 1, t_end = 336, dt = 0.5) {
  times <- seq(0, t_end, by = dt)
  n     <- length(times)

  Adepot <- 0
  Cp_vec <- numeric(n)

  for (i in seq_len(n)) {
    t <- times[i]

    # Dose events
    dose_event <- 0
    if (n_doses == 1) {
      if (i == 1) { Adepot <- Adepot + F_oral * dose; dose_event <- 1 }
    } else {
      for (d in seq_len(n_doses)) {
        dose_time <- (d - 1) * ii
        if (abs(t - dose_time) < dt / 2) { Adepot <- Adepot + F_oral * dose; dose_event <- 1 }
      }
    }

    # Euler step
    dAdepot <- -ka * Adepot
    Cp      <- if (Vd > 0) Adepot / Vd else 0          # approximate: depot → central directly
    # Better 2-state: depot + central
    # We do it properly below:
    Adepot  <- max(0, Adepot + dAdepot * dt)
    Cp_vec[i] <- max(0, Adepot / Vd)                    # simplified; reset below with proper ODE
  }

  # ---- proper 2-compartment Euler (depot + central) ----
  Adepot  <- 0
  Acentral <- 0
  Cp_out  <- numeric(n)
  de_out  <- integer(n)

  for (i in seq_len(n)) {
    t <- times[i]

    for (d in seq_len(n_doses)) {
      dose_time <- (d - 1) * ii
      if (abs(t - dose_time) < dt / 2) {
        Adepot <- Adepot + F_oral * dose
        de_out[i] <- 1L
      }
    }

    dAdepot   <- -ka * Adepot
    dAcentral <-  ka * Adepot - (CL / Vd) * Acentral

    Adepot   <- max(0, Adepot   + dAdepot   * dt)
    Acentral <- max(0, Acentral + dAcentral * dt)

    Cp_out[i] <- Acentral / Vd
  }

  data.frame(time = times, Cp = Cp_out, dose_event = de_out)
}

# -----------------------------------------------------------------------------
# Drug PK parameter library
# Returns list(ka, Vd_L, CL_L_h, F, ii_h, label_unit)
# -----------------------------------------------------------------------------
drug_pk_params <- function(drug_name, dose_frac = 1.0,
                           renal = "Normal", hepatic = "Normal") {
  params <- switch(drug_name,
    "Palbociclib 125mg" = list(
      dose = 125 * dose_frac, ka = 0.3, Vd = 2583, CL = 101.4 / 24,
      F = 0.46, ii = 24, n_doses = 1, label = "Palbociclib (ng/mL)"
    ),
    "Letrozole 2.5mg" = list(
      dose = 2500 * dose_frac, ka = 0.7, Vd = 187, CL = 2.1 / 24,
      F = 0.99, ii = 24, n_doses = 1, label = "Letrozole (ng/mL)"
    ),
    "Ribociclib 600mg" = list(
      dose = 600 * dose_frac, ka = 0.5, Vd = 1090, CL = 25.8 / 24,
      F = 0.56, ii = 24, n_doses = 1, label = "Ribociclib (ng/mL)"
    ),
    "Trastuzumab 6mg/kg" = list(
      dose = 6000 * dose_frac, ka = 0.02, Vd = 3.0, CL = 0.225 / 24,
      F = 1.0, ii = 504, n_doses = 1, label = "Trastuzumab (µg/mL)"
    ),
    "Olaparib 300mg BID" = list(
      dose = 300 * dose_frac, ka = 1.2, Vd = 158, CL = 7.3 / 24,
      F = 0.77, ii = 12, n_doses = 2, label = "Olaparib (ng/mL)"
    ),
    "Abemaciclib 150mg BID" = list(
      dose = 150 * dose_frac, ka = 0.4, Vd = 690, CL = 18.9 / 24,
      F = 0.45, ii = 12, n_doses = 2, label = "Abemaciclib (ng/mL)"
    ),
    "Pembrolizumab 200mg q3w" = list(
      dose = 200000 * dose_frac, ka = 0.008, Vd = 3.8, CL = 0.21 / 24,
      F = 1.0, ii = 504, n_doses = 1, label = "Pembrolizumab (µg/mL)"
    )
  )

  # Renal / hepatic adjustments (simplified scalar)
  renal_adj <- switch(renal,
    "Normal" = 1.0, "Mild (eGFR 60-89)" = 0.9,
    "Moderate (eGFR 30-59)" = 0.75, "Severe (eGFR <30)" = 0.6, 1.0
  )
  hepatic_adj <- switch(hepatic,
    "Normal" = 1.0, "Child-Pugh A" = 0.85,
    "Child-Pugh B" = 0.65, "Child-Pugh C" = 0.45, 1.0
  )

  params$CL <- params$CL * renal_adj * hepatic_adj
  params
}

# -----------------------------------------------------------------------------
# Tumor volume ODE (Simeoni / TGI-like) via Euler
# Returns data.frame(time_wk, vol_cm3, ki67, cdk_inhib)
# -----------------------------------------------------------------------------
sim_tumor <- function(regimen, weeks = 52, init_diam_cm = 3,
                      esr1_mut = FALSE, pik3ca_mut = FALSE,
                      baseline_ki67 = 30) {

  # Initial volume (sphere): V = (4/3)*pi*(d/2)^3
  V0 <- (4 / 3) * pi * (init_diam_cm / 2)^3

  # Regimen parameters: (growth_rate/day, kill_rate/day, resistance_half_life_days)
  reg_params <- switch(regimen,
    "Letrozole alone"                   = list(kg = 0.012, kk = 0.018, t_res = 365),
    "Palbociclib + Letrozole"           = list(kg = 0.012, kk = 0.045, t_res = 500),
    "Ribociclib + Letrozole"            = list(kg = 0.012, kk = 0.042, t_res = 480),
    "Abemaciclib + Letrozole"           = list(kg = 0.012, kk = 0.046, t_res = 510),
    "Trastuzumab + Pertuzumab + Chemo"  = list(kg = 0.015, kk = 0.060, t_res = 400),
    "Pembrolizumab + Chemo (TNBC)"      = list(kg = 0.018, kk = 0.055, t_res = 300),
    "Olaparib (BRCAm)"                  = list(kg = 0.014, kk = 0.050, t_res = 350),
    list(kg = 0.012, kk = 0.030, t_res = 400)
  )

  kg  <- reg_params$kg
  kk  <- reg_params$kk
  t_r <- reg_params$t_res

  # Mutation penalties
  if (esr1_mut)  kk <- kk * 0.60
  if (pik3ca_mut) kk <- kk * 0.75

  # Euler integration (daily steps)
  days <- seq(0, weeks * 7, by = 1)
  V    <- numeric(length(days))
  V[1] <- V0
  Vmax <- V0 * 50   # cap

  for (i in seq(2, length(days))) {
    t   <- days[i]
    # Resistance grows over time (exponential decay of kill rate)
    kk_t <- kk * exp(-t / t_r)
    dV   <- (kg - kk_t) * V[i - 1]
    V[i] <- max(1e-6, min(V[i - 1] + dV, Vmax))
  }

  vol_cm3 <- V
  diam    <- (vol_cm3 * 3 / (4 * pi))^(1 / 3) * 2  # diameter cm

  # Ki-67 proxy: drops with treatment, rebounds with resistance
  ki67_t <- baseline_ki67 * (V / V0)^0.4
  ki67_t <- pmax(2, pmin(ki67_t, 100))

  # CDK4/6 inhibition proxy (for CDK4/6i regimens)
  cdk_inhib <- if (grepl("Palbociclib|Ribociclib|Abemaciclib", regimen)) {
    pmax(0, 80 * exp(-days / t_r) * (1 - days / (weeks * 7 + 1)))
  } else {
    rep(0, length(days))
  }

  data.frame(
    time_day = days,
    time_wk  = days / 7,
    vol_cm3  = vol_cm3,
    diam_cm  = diam,
    ki67     = ki67_t,
    cdk_inhib = pmin(100, pmax(0, cdk_inhib))
  )
}

# -----------------------------------------------------------------------------
# Immune microenvironment ODE
# Compartments: TIL (CD8), Treg, tumor antigen load, PD-L1 expression
# -----------------------------------------------------------------------------
sim_immune <- function(immuno_drug, baseline_til, pdl1_cps, weeks,
                       brca_mut = FALSE, subtype = "TNBC") {

  til0   <- switch(baseline_til,
    "Low (<10%)"   = 5,
    "Medium (10-30%)" = 20,
    "High (>30%)"  = 45, 20
  )
  treg0  <- til0 * 0.3
  ag0    <- 100   # arbitrary antigen units

  # Drug effects
  pd1_block  <- immuno_drug %in% c("Pembrolizumab", "Atezolizumab")
  ctla4_block <- immuno_drug == "Ipilimumab"

  days <- seq(0, weeks * 7, by = 1)
  TIL  <- numeric(length(days))
  Treg <- numeric(length(days))
  Ag   <- numeric(length(days))

  TIL[1]  <- til0
  Treg[1] <- treg0
  Ag[1]   <- ag0

  for (i in seq(2, length(days))) {
    t <- days[i]

    # PD-1/PD-L1 blockade: boosts TIL proliferation, reduces exhaustion
    pd1_boost   <- if (pd1_block)  1 + 0.8 * (pdl1_cps / 100) else 1.0
    ctla4_boost <- if (ctla4_block) 1.3 else 1.0
    treg_supp   <- if (ctla4_block) 0.6 else 1.0

    # BRCA mutation → higher mutational burden → more antigen
    ag_load <- if (brca_mut) 1.3 else 1.0
    # TNBC has higher immunogenicity
    subtype_boost <- if (subtype == "TNBC") 1.2 else 1.0

    dTIL  <-  0.05 * pd1_boost * ctla4_boost * subtype_boost * TIL[i-1] -
               0.03 * TIL[i-1] - 0.002 * TIL[i-1] * Treg[i-1] -
               0.001 * TIL[i-1] * Ag[i-1] / (Ag[i-1] + 50)

    dTreg <- -0.02 * treg_supp * Treg[i-1] + 0.01 * Treg[i-1]

    dAg   <- -0.04 * ag_load * TIL[i-1] * Ag[i-1] / (Ag[i-1] + 30) +
               0.005 * Ag[i-1]

    TIL[i]  <- max(0, TIL[i-1]  + dTIL)
    Treg[i] <- max(0, Treg[i-1] + dTreg)
    Ag[i]   <- max(0, Ag[i-1]   + dAg)
  }

  # pCR probability (logistic based on TIL at week 12 relative to baseline)
  til_week12 <- TIL[min(85, length(TIL))]
  log_odds <- -1.5 + 0.04 * til_week12 + 0.02 * pdl1_cps +
               if (brca_mut) 0.5 else 0 +
               if (subtype == "TNBC") 0.3 else 0
  pcr_prob <- 1 / (1 + exp(-log_odds))

  list(
    traj = data.frame(
      day  = days,
      week = days / 7,
      TIL  = TIL,
      Treg = Treg,
      Ag   = Ag
    ),
    pcr_prob = round(pcr_prob * 100, 1)
  )
}

# -----------------------------------------------------------------------------
# Weibull PFS survival curve
# S(t) = exp(-(t/scale)^shape), scale = median / (log(2)^(1/shape))
# -----------------------------------------------------------------------------
weibull_surv <- function(t, median_pfs, shape) {
  scale <- median_pfs / (log(2)^(1 / shape))
  exp(-(t / scale)^shape)
}

pfs_params <- list(
  "Letrozole alone"                     = list(median = 14.5, shape = 1.2),
  "Palbociclib + Letrozole"             = list(median = 27.6, shape = 1.3),
  "Ribociclib + Letrozole"              = list(median = 25.3, shape = 1.2),
  "Abemaciclib + Letrozole"             = list(median = 28.2, shape = 1.3),
  "Trastuzumab + Pertuzumab + Docetaxel"= list(median = 18.7, shape = 1.2),
  "Pembrolizumab + Chemo (TNBC)"        = list(median =  9.7, shape = 1.0),
  "Olaparib (BRCAm)"                    = list(median =  7.0, shape = 1.1)
)

# ORR lookup (approximate published values, %)
orr_lookup <- c(
  "Letrozole alone"                      = 21,
  "Palbociclib + Letrozole"              = 42,
  "Ribociclib + Letrozole"               = 41,
  "Abemaciclib + Letrozole"              = 48,
  "Trastuzumab + Pertuzumab + Docetaxel" = 56,
  "Pembrolizumab + Chemo (TNBC)"         = 53,
  "Olaparib (BRCAm)"                     = 60
)

# -----------------------------------------------------------------------------
# Resistance & biomarker simulation
# Returns list of ctDNA, ANC, LVEF trajectories
# -----------------------------------------------------------------------------
sim_resistance <- function(resist_mut, scenario, baseline_ctdna_vaf,
                           baseline_lvef, anthracycline, weeks = 52) {

  days <- seq(0, weeks * 7, by = 1)
  n    <- length(days)

  # ctDNA dynamics — VAF rises with resistance
  resist_boost <- switch(resist_mut,
    "None"           = 0.0,
    "ESR1 Y537S"     = 0.8,
    "ESR1 D538G"     = 0.7,
    "PIK3CA H1047R"  = 0.5,
    "PTEN loss"      = 0.6,
    0.0
  )
  ctdna <- numeric(n)
  ctdna[1] <- baseline_ctdna_vaf
  for (i in seq(2, n)) {
    decay  <- 0.008
    growth <- 0.003 * resist_boost * (days[i] / 365)
    ctdna[i] <- max(0, ctdna[i - 1] * (1 - decay) + growth * ctdna[i - 1] + 0.01)
    ctdna[i] <- min(ctdna[i], 80)
  }

  # ANC nadir (CDK4/6i causes neutropenia in ~80%)
  anc_base <- 4.5   # normal 1.8–7.7 × 10^9/L
  anc <- numeric(n)
  anc[1] <- anc_base
  for (i in seq(2, n)) {
    # Cyclic pattern: nadir at day 15, recovery by day 28 per cycle
    cycle_day <- (days[i] - 1) %% 28 + 1
    nadir_factor <- if (cycle_day <= 15) {
      1 - 0.55 * sin(pi * cycle_day / 15)
    } else {
      1 - 0.55 * sin(pi * (28 - cycle_day) / 13)
    }
    anc[i] <- max(0.2, anc_base * nadir_factor + rnorm(1, 0, 0.05))
  }

  # LVEF (cardiotoxicity — relevant for anthracyclines & trastuzumab)
  lvef <- numeric(n)
  lvef[1] <- baseline_lvef
  lvef_drop_rate <- if (anthracycline) 0.004 else 0.001
  for (i in seq(2, n)) {
    noise  <- rnorm(1, 0, 0.1)
    lvef[i] <- max(20, lvef[i - 1] - lvef_drop_rate + noise)
  }

  list(
    days  = days,
    ctdna = ctdna,
    anc   = anc,
    lvef  = lvef
  )
}

# =============================================================================
# UI
# =============================================================================

ui <- dashboardPage(
  skin = "red",

  dashboardHeader(
    title = "Breast Cancer QSP Dashboard",
    titleWidth = 300
  ),

  dashboardSidebar(
    width = 230,
    sidebarMenu(
      id = "tabs",
      menuItem("Patient Profile",      tabName = "tab_patient",  icon = icon("user-md")),
      menuItem("Drug PK Profiles",     tabName = "tab_pk",       icon = icon("pills")),
      menuItem("Tumor Dynamics & PD",  tabName = "tab_tumor",    icon = icon("dna")),
      menuItem("Immune Microenvironment", tabName = "tab_immune", icon = icon("shield-halved")),
      menuItem("Clinical Endpoints",   tabName = "tab_clinical", icon = icon("chart-line")),
      menuItem("Scenario & Biomarkers",tabName = "tab_scenario", icon = icon("flask"))
    ),
    br(),
    div(style = "padding:10px; color:#ccc; font-size:11px;",
        "QSP Model v1.0 | 2025",
        br(), "Euler ODE simulation",
        br(), "Parameters from published trials")
  ),

  dashboardBody(
    tags$head(tags$style(HTML("
      .content-wrapper, .right-side { background-color: #f4f4f4; }
      .box.box-solid.box-danger > .box-header { color:#fff; background:#c0392b; }
      .box.box-solid.box-primary > .box-header { color:#fff; background:#2980b9; }
      .summary-card { background:#fff; border-left:5px solid #c0392b;
                      padding:15px; border-radius:4px; margin-bottom:10px; }
      .summary-card h4 { color:#c0392b; margin-top:0; }
      .metric-badge { display:inline-block; background:#c0392b; color:#fff;
                      border-radius:3px; padding:2px 8px; font-size:13px; }
      .metric-ok   { background:#27ae60; }
      .metric-warn { background:#e67e22; }
    "))),

    tabItems(

      # -----------------------------------------------------------------------
      # TAB 1: Patient Profile & Disease Classification
      # -----------------------------------------------------------------------
      tabItem(tabName = "tab_patient",
        fluidRow(
          box(title = "Patient Characteristics", status = "danger", solidHeader = TRUE,
              width = 4, collapsible = FALSE,
              selectInput("pt_subtype", "Molecular Subtype",
                choices = c("ER+/HER2- (Luminal A/B)", "HER2+ ER+ (Luminal B HER2+)",
                            "HER2+ ER- (HER2-enriched)", "TNBC"),
                selected = "ER+/HER2- (Luminal A/B)"),
              sliderInput("pt_age", "Age (years)", 30, 80, 55, step = 1),
              radioButtons("pt_meno", "Menopausal Status",
                choices = c("Pre-menopausal", "Peri-menopausal", "Post-menopausal"),
                selected = "Post-menopausal"),
              sliderInput("pt_er", "ER Expression (%)", 0, 100, 75, step = 1),
              sliderInput("pt_pr", "PR Expression (%)", 0, 100, 60, step = 1)
          ),
          box(title = "Biomarker Profile", status = "danger", solidHeader = TRUE,
              width = 4, collapsible = FALSE,
              selectInput("pt_her2", "HER2 Status",
                choices = c("0 (Negative)", "1+ (Negative)", "2+ (Equivocal/FISH-)",
                            "2+ (Equivocal/FISH+)", "3+ (Positive)"),
                selected = "0 (Negative)"),
              sliderInput("pt_ki67", "Ki-67 Proliferation Index (%)", 0, 80, 25, step = 1),
              sliderInput("pt_pdl1", "PD-L1 CPS Score", 0, 100, 10, step = 1),
              sliderInput("pt_oncotype", "Oncotype DX Recurrence Score", 0, 100, 22, step = 1),
              selectInput("pt_stage", "Clinical Stage",
                choices = c("Stage I", "Stage II", "Stage III", "Stage IV (metastatic)"),
                selected = "Stage II"),
              checkboxInput("pt_brca", "Germline BRCA1/2 Mutation", value = FALSE)
          ),
          box(title = "Patient Summary", status = "primary", solidHeader = TRUE,
              width = 4, collapsible = FALSE,
              uiOutput("patient_summary_card")
          )
        ),
        fluidRow(
          box(title = "Biomarker Radar Chart", status = "primary", solidHeader = TRUE,
              width = 6,
              plotlyOutput("biomarker_radar", height = 320)),
          box(title = "Risk Stratification", status = "primary", solidHeader = TRUE,
              width = 6,
              plotlyOutput("risk_bar", height = 320))
        )
      ),

      # -----------------------------------------------------------------------
      # TAB 2: Drug PK Profiles
      # -----------------------------------------------------------------------
      tabItem(tabName = "tab_pk",
        fluidRow(
          box(title = "PK Simulation Parameters", status = "danger", solidHeader = TRUE,
              width = 3,
              selectInput("pk_drug", "Drug",
                choices = c("Palbociclib 125mg", "Letrozole 2.5mg",
                            "Ribociclib 600mg", "Trastuzumab 6mg/kg",
                            "Olaparib 300mg BID", "Abemaciclib 150mg BID",
                            "Pembrolizumab 200mg q3w"),
                selected = "Palbociclib 125mg"),
              sliderInput("pk_dose_frac", "Relative Dose (fraction)", 0.5, 1.0, 1.0, step = 0.05),
              sliderInput("pk_days", "Simulation Duration (days)", 7, 90, 28, step = 1),
              selectInput("pk_renal", "Renal Function",
                choices = c("Normal", "Mild (eGFR 60-89)",
                            "Moderate (eGFR 30-59)", "Severe (eGFR <30)"),
                selected = "Normal"),
              selectInput("pk_hepatic", "Hepatic Function",
                choices = c("Normal", "Child-Pugh A", "Child-Pugh B", "Child-Pugh C"),
                selected = "Normal"),
              numericInput("pk_cv_pct", "Inter-individual CV% (variability)", 30, 5, 60, step = 5),
              actionButton("pk_run", "Run Simulation", icon = icon("play"),
                           class = "btn-danger btn-block")
          ),
          box(title = "PK Concentration-Time Profile", status = "primary", solidHeader = TRUE,
              width = 9,
              plotlyOutput("pk_plot", height = 420),
              hr(),
              fluidRow(
                valueBoxOutput("pk_cmax", width = 3),
                valueBoxOutput("pk_tmax", width = 3),
                valueBoxOutput("pk_auc",  width = 3),
                valueBoxOutput("pk_cmin", width = 3)
              )
          )
        ),
        fluidRow(
          box(title = "Dose Proportionality", status = "primary", solidHeader = TRUE,
              width = 6,
              plotlyOutput("pk_dose_prop", height = 280)),
          box(title = "Organ Function Impact on AUC", status = "primary", solidHeader = TRUE,
              width = 6,
              plotlyOutput("pk_organ_impact", height = 280))
        )
      ),

      # -----------------------------------------------------------------------
      # TAB 3: Tumor Dynamics & PD
      # -----------------------------------------------------------------------
      tabItem(tabName = "tab_tumor",
        fluidRow(
          box(title = "Tumor Model Parameters", status = "danger", solidHeader = TRUE,
              width = 3,
              selectInput("td_regimen", "Treatment Regimen",
                choices = c("Letrozole alone",
                            "Palbociclib + Letrozole",
                            "Ribociclib + Letrozole",
                            "Abemaciclib + Letrozole",
                            "Trastuzumab + Pertuzumab + Chemo",
                            "Pembrolizumab + Chemo (TNBC)",
                            "Olaparib (BRCAm)"),
                selected = "Palbociclib + Letrozole"),
              sliderInput("td_weeks", "Simulation Duration (weeks)", 12, 104, 52, step = 4),
              sliderInput("td_init_diam", "Initial Tumor Diameter (cm)", 1, 10, 3, step = 0.5),
              sliderInput("td_ki67", "Baseline Ki-67 (%)", 5, 80, 30, step = 1),
              checkboxInput("td_esr1", "ESR1 Mutation (resistance)", value = FALSE),
              checkboxInput("td_pik3ca", "PIK3CA Mutation (resistance)", value = FALSE),
              br(),
              div(style = "font-size:12px; color:#777;",
                  "Model: modified Simeoni TGI with",
                  "exponential resistance decay")
          ),
          box(title = "Tumor Volume Dynamics", status = "primary", solidHeader = TRUE,
              width = 9,
              plotlyOutput("tumor_vol_plot", height = 280)),
        ),
        fluidRow(
          box(title = "Ki-67 Proliferation Index Over Time", status = "primary", solidHeader = TRUE,
              width = 4,
              plotlyOutput("ki67_plot", height = 260)),
          box(title = "CDK4/6 Inhibition (%)", status = "primary", solidHeader = TRUE,
              width = 4,
              plotlyOutput("cdk_plot", height = 260)),
          box(title = "PD Metrics Summary", status = "primary", solidHeader = TRUE,
              width = 4,
              tableOutput("pd_table"),
              br(),
              uiOutput("response_badge"))
        )
      ),

      # -----------------------------------------------------------------------
      # TAB 4: Immune Microenvironment
      # -----------------------------------------------------------------------
      tabItem(tabName = "tab_immune",
        fluidRow(
          box(title = "Immunotherapy Parameters", status = "danger", solidHeader = TRUE,
              width = 3,
              selectInput("imm_drug", "Immunotherapy Agent",
                choices = c("None", "Pembrolizumab", "Atezolizumab", "Ipilimumab"),
                selected = "Pembrolizumab"),
              selectInput("imm_til", "Baseline TIL Level",
                choices = c("Low (<10%)", "Medium (10-30%)", "High (>30%)"),
                selected = "Medium (10-30%)"),
              sliderInput("imm_pdl1", "PD-L1 CPS Score", 0, 100, 15, step = 1),
              sliderInput("imm_weeks", "Simulation Duration (weeks)", 8, 52, 24, step = 4),
              checkboxInput("imm_brca", "BRCA Mutation (↑TMB)", value = FALSE),
              selectInput("imm_subtype", "Tumor Subtype",
                choices = c("TNBC", "HER2+", "ER+/HER2-"),
                selected = "TNBC"),
              actionButton("imm_run", "Simulate Immune Response",
                           icon = icon("play"), class = "btn-danger btn-block")
          ),
          box(title = "Immune Cell Dynamics", status = "primary", solidHeader = TRUE,
              width = 9,
              plotlyOutput("immune_plot", height = 340),
              hr(),
              fluidRow(
                valueBoxOutput("pcr_box",   width = 4),
                valueBoxOutput("til_box",   width = 4),
                valueBoxOutput("treg_box",  width = 4)
              )
          )
        ),
        fluidRow(
          box(title = "pCR Probability vs PD-L1 CPS", status = "primary", solidHeader = TRUE,
              width = 6,
              plotlyOutput("pcr_pdl1_plot", height = 280)),
          box(title = "Trial Benchmark Comparison", status = "primary", solidHeader = TRUE,
              width = 6,
              plotlyOutput("trial_bench_plot", height = 280))
        )
      ),

      # -----------------------------------------------------------------------
      # TAB 5: Clinical Endpoint Simulation
      # -----------------------------------------------------------------------
      tabItem(tabName = "tab_clinical",
        fluidRow(
          box(title = "Endpoint Simulation Parameters", status = "danger", solidHeader = TRUE,
              width = 3,
              checkboxGroupInput("clin_regimens", "Regimens to Compare",
                choices = names(pfs_params),
                selected = c("Letrozole alone",
                             "Palbociclib + Letrozole",
                             "Ribociclib + Letrozole")),
              sliderInput("clin_months", "Follow-up Duration (months)", 12, 60, 36, step = 1),
              selectInput("clin_subtype", "Patient Subtype",
                choices = c("ER+/HER2-", "HER2+", "TNBC", "BRCAm"),
                selected = "ER+/HER2-"),
              numericInput("clin_n_patients", "Hypothetical cohort size", 200, 50, 1000, step = 50),
              actionButton("clin_run", "Simulate Endpoints",
                           icon = icon("play"), class = "btn-danger btn-block")
          ),
          box(title = "Progression-Free Survival Curves (Weibull)", status = "primary",
              solidHeader = TRUE, width = 9,
              plotlyOutput("pfs_plot", height = 380))
        ),
        fluidRow(
          box(title = "Overall Response Rate (ORR)", status = "primary", solidHeader = TRUE,
              width = 5,
              plotlyOutput("orr_plot", height = 280)),
          box(title = "Endpoint Comparison Table", status = "primary", solidHeader = TRUE,
              width = 7,
              tableOutput("endpoint_table"))
        )
      ),

      # -----------------------------------------------------------------------
      # TAB 6: Scenario Comparison & Biomarkers
      # -----------------------------------------------------------------------
      tabItem(tabName = "tab_scenario",
        fluidRow(
          box(title = "Scenario Parameters", status = "danger", solidHeader = TRUE,
              width = 3,
              selectInput("sc_resist", "Resistance Mutation",
                choices = c("None", "ESR1 Y537S", "ESR1 D538G",
                            "PIK3CA H1047R", "PTEN loss"),
                selected = "None"),
              selectInput("sc_compare", "Comparison Scenario",
                choices = c("Letrozole alone",
                            "Palbociclib + Letrozole",
                            "Ribociclib + Letrozole",
                            "Abemaciclib + Letrozole",
                            "Olaparib (BRCAm)"),
                selected = "Palbociclib + Letrozole"),
              sliderInput("sc_ctdna_vaf", "Baseline ctDNA VAF (%)", 0, 50, 5, step = 0.5),
              sliderInput("sc_lvef", "Baseline LVEF (%)", 40, 70, 62, step = 1),
              checkboxInput("sc_anthracycline", "Anthracycline-containing Regimen", value = FALSE),
              sliderInput("sc_weeks", "Monitoring Duration (weeks)", 12, 52, 24, step = 4),
              actionButton("sc_run", "Run Scenarios", icon = icon("play"),
                           class = "btn-danger btn-block")
          ),
          box(title = "Resistance Impact on Tumor Dynamics", status = "primary",
              solidHeader = TRUE, width = 9,
              plotlyOutput("resist_plot", height = 300))
        ),
        fluidRow(
          box(title = "ctDNA VAF Dynamics", status = "primary", solidHeader = TRUE,
              width = 4,
              plotlyOutput("ctdna_plot", height = 250)),
          box(title = "Toxicity: ANC & Neutropenia Risk", status = "primary",
              solidHeader = TRUE, width = 4,
              plotlyOutput("anc_plot", height = 250)),
          box(title = "Cardiac Safety: LVEF Monitoring", status = "primary",
              solidHeader = TRUE, width = 4,
              plotlyOutput("lvef_plot", height = 250))
        ),
        fluidRow(
          box(title = "Treatment Response Waterfall Chart", status = "primary",
              solidHeader = TRUE, width = 12,
              plotlyOutput("waterfall_plot", height = 280))
        )
      )
    )   # end tabItems
  )   # end dashboardBody
)   # end dashboardPage

# =============================================================================
# SERVER
# =============================================================================

server <- function(input, output, session) {

  set.seed(42)   # reproducible random variability

  # ===========================================================================
  # TAB 1: Patient Profile
  # ===========================================================================

  # Helper: derive clinical metrics from patient inputs
  patient_metrics <- reactive({
    subtype   <- input$pt_subtype
    age       <- input$pt_age
    er        <- input$pt_er
    pr        <- input$pt_pr
    ki67      <- input$pt_ki67
    pdl1      <- input$pt_pdl1
    oncotype  <- input$pt_oncotype
    her2_str  <- input$pt_her2
    stage     <- input$pt_stage
    brca      <- input$pt_brca

    her2_pos  <- grepl("FISH\\+|3\\+", her2_str)

    # Risk category
    risk <- if (grepl("TNBC", subtype)) "High"
            else if (her2_pos || ki67 > 30 || oncotype > 25) "Intermediate-High"
            else if (ki67 < 15 && oncotype < 18) "Low"
            else "Intermediate"

    # Recommended treatment (simplified NCCN-like)
    tx <- if (grepl("TNBC", subtype)) {
      if (brca) "Pembrolizumab + Chemo (neoadjuvant); Olaparib (maintenance)"
      else if (pdl1 >= 10) "Pembrolizumab + nab-Paclitaxel / Gem-Carbo"
      else "Anthracycline + Taxane-based chemotherapy"
    } else if (her2_pos && grepl("ER\\+", subtype)) {
      "Pertuzumab + Trastuzumab + Docetaxel → Palbociclib + Letrozole"
    } else if (her2_pos) {
      "Pertuzumab + Trastuzumab + Docetaxel"
    } else {
      # ER+/HER2-
      if (oncotype > 25 || ki67 > 30) "Palbociclib + Letrozole (or Ribociclib + Letrozole)"
      else if (oncotype <= 18) "Endocrine therapy alone (Letrozole/Anastrozole)"
      else "Abemaciclib + Letrozole (Monarchies E criteria)"
    }

    # pCR probability estimate (logistic model, simplified)
    pcr_logit <- -2.5 +
      (if (grepl("TNBC", subtype)) 1.2 else 0) +
      (if (her2_pos) 0.9 else 0) +
      (if (pdl1 > 10) 0.4 else 0) +
      (if (brca) 0.6 else 0) +
      (if (ki67 > 30) 0.3 else 0) +
      (if (stage == "Stage II") -0.2 else if (stage == "Stage III") 0.1 else -0.5)
    pcr <- round(100 / (1 + exp(-pcr_logit)), 1)

    # 5-year survival estimate
    surv5 <- if (stage == "Stage I") 99
              else if (stage == "Stage II") {
                if (grepl("TNBC", subtype)) 77 else 86
              } else if (stage == "Stage III") {
                if (grepl("TNBC", subtype)) 65 else 72
              } else 28    # Stage IV

    list(risk = risk, tx = tx, pcr = pcr, surv5 = surv5)
  })

  output$patient_summary_card <- renderUI({
    m <- patient_metrics()
    risk_col <- switch(m$risk,
      "Low" = "metric-ok", "Intermediate" = "metric-warn",
      "Intermediate-High" = "metric-warn", "High" = "metric-badge", "metric-badge")

    div(
      div(class = "summary-card",
        h4("Risk Category"),
        span(class = paste("metric-badge", risk_col), m$risk)
      ),
      div(class = "summary-card",
        h4("Recommended Treatment"),
        p(style = "font-size:13px;", m$tx)
      ),
      div(class = "summary-card",
        h4("Estimated pCR Probability"),
        span(class = "metric-badge", paste0(m$pcr, "%"))
      ),
      div(class = "summary-card",
        h4("5-Year Survival Estimate"),
        span(class = paste("metric-badge", if (m$surv5 > 80) "metric-ok" else "metric-warn"),
             paste0(m$surv5, "%"))
      )
    )
  })

  output$biomarker_radar <- renderPlotly({
    vals <- c(input$pt_er, input$pt_pr, input$pt_ki67,
              input$pt_pdl1, input$pt_oncotype,
              min(100, input$pt_age * 1.25))
    cats <- c("ER%", "PR%", "Ki-67%", "PD-L1 CPS",
              "Oncotype DX", "Age-adj Risk")
    df <- data.frame(category = cats, value = vals)

    plot_ly(
      type = "scatterpolar",
      r    = c(vals, vals[1]),
      theta = c(cats, cats[1]),
      fill = "toself",
      fillcolor = "rgba(192,57,43,0.25)",
      line = list(color = "#c0392b")
    ) |>
      layout(
        polar = list(radialaxis = list(visible = TRUE, range = c(0, 100))),
        showlegend = FALSE,
        margin = list(l = 40, r = 40, t = 20, b = 20)
      )
  })

  output$risk_bar <- renderPlotly({
    m <- patient_metrics()
    # Relative contribution scores
    scores <- c(
      "Subtype"        = if (grepl("TNBC", input$pt_subtype)) 40 else if (grepl("HER2\\+", input$pt_subtype)) 30 else 10,
      "Ki-67"          = input$pt_ki67 * 0.5,
      "Oncotype DX"    = input$pt_oncotype * 0.4,
      "Stage"          = c("Stage I"=5,"Stage II"=15,"Stage III"=35,"Stage IV (metastatic)"=60)[input$pt_stage],
      "PD-L1 CPS"      = input$pt_pdl1 * 0.2,
      "BRCA Mutation"  = if (input$pt_brca) 20 else 0
    )
    df <- data.frame(
      Factor = names(scores),
      Score  = as.numeric(scores)
    )
    plot_ly(df, x = ~Score, y = ~Factor, type = "bar", orientation = "h",
            marker = list(color = "#c0392b")) |>
      layout(
        xaxis = list(title = "Risk Contribution Score"),
        yaxis = list(title = "", categoryorder = "total ascending"),
        margin = list(l = 120, t = 20, b = 40, r = 20)
      )
  })

  # ===========================================================================
  # TAB 2: Drug PK Profiles
  # ===========================================================================

  pk_data <- eventReactive(input$pk_run, {
    params <- drug_pk_params(input$pk_drug, input$pk_dose_frac,
                             input$pk_renal, input$pk_hepatic)
    t_end  <- input$pk_days * 24
    n_d    <- max(1, round(input$pk_days * 24 / params$ii))

    df <- sim_pk_1cpt(
      dose    = params$dose,
      ka      = params$ka,
      Vd      = params$Vd,
      CL      = params$CL,
      F_oral  = params$F,
      ii      = params$ii,
      n_doses = n_d,
      t_end   = t_end,
      dt      = 0.25
    )
    df$label <- params$label

    # Monte Carlo variability band (log-normal IIV)
    cv  <- input$pk_cv_pct / 100
    set.seed(123)
    sims <- lapply(1:50, function(i) {
      cl_i  <- params$CL  * exp(rnorm(1, 0, cv))
      vd_i  <- params$Vd  * exp(rnorm(1, 0, cv * 0.5))
      s     <- sim_pk_1cpt(params$dose, params$ka, vd_i, cl_i,
                           params$F, params$ii, n_d, t_end, dt = 0.5)
      s$Cp
    })
    sim_mat <- do.call(cbind, sims)
    df$Cp_lo <- apply(sim_mat, 1, quantile, 0.05)
    df$Cp_hi <- apply(sim_mat, 1, quantile, 0.95)

    df
  }, ignoreNULL = FALSE)

  observeEvent(input$pk_drug, {
    # Trigger re-run when drug changes
    shinyjs::click("pk_run")
  }, ignoreInit = TRUE)

  output$pk_plot <- renderPlotly({
    df <- pk_data()
    if (is.null(df)) return(NULL)

    p <- plot_ly(df, x = ~time / 24) |>
      add_ribbons(ymin = ~Cp_lo, ymax = ~Cp_hi,
                  fillcolor = "rgba(192,57,43,0.15)",
                  line = list(color = "transparent"),
                  name = "90% PI") |>
      add_lines(y = ~Cp, line = list(color = "#c0392b", width = 2),
                name = df$label[1]) |>
      layout(
        xaxis = list(title = "Time (days)"),
        yaxis = list(title = df$label[1]),
        hovermode = "x unified",
        legend  = list(orientation = "h"),
        margin  = list(t = 20, b = 50, l = 60, r = 20)
      )
    p
  })

  pk_stats <- reactive({
    df <- pk_data()
    if (is.null(df)) return(list(cmax = 0, tmax = 0, auc = 0, cmin = 0))
    # Steady-state window: last dosing interval
    # Use last 10% of simulation
    n   <- nrow(df)
    ss  <- df[round(n * 0.85):n, ]
    cmax <- max(ss$Cp)
    tmax <- ss$time[which.max(ss$Cp)] / 24
    auc  <- sum(diff(ss$time) * (head(ss$Cp, -1) + tail(ss$Cp, -1)) / 2)
    cmin <- min(ss$Cp)
    list(cmax = round(cmax, 2), tmax = round(tmax, 2),
         auc  = round(auc, 1),  cmin = round(cmin, 2))
  })

  output$pk_cmax <- renderValueBox({
    s <- pk_stats()
    valueBox(s$cmax, "Cmax (SS)", icon = icon("arrow-up"), color = "red")
  })
  output$pk_tmax <- renderValueBox({
    s <- pk_stats()
    valueBox(paste0(s$tmax, "d"), "Tmax (SS)", icon = icon("clock"), color = "orange")
  })
  output$pk_auc <- renderValueBox({
    s <- pk_stats()
    valueBox(s$auc, "AUC (SS)", icon = icon("chart-area"), color = "blue")
  })
  output$pk_cmin <- renderValueBox({
    s <- pk_stats()
    valueBox(s$cmin, "Ctrough (SS)", icon = icon("arrow-down"), color = "green")
  })

  output$pk_dose_prop <- renderPlotly({
    params <- drug_pk_params(input$pk_drug, 1.0, "Normal", "Normal")
    fracs  <- seq(0.25, 1.5, by = 0.25)
    aucs   <- vapply(fracs, function(f) {
      d2 <- sim_pk_1cpt(params$dose * f, params$ka, params$Vd, params$CL,
                        params$F, params$ii, 7, 168, dt = 0.5)
      ss <- d2[round(nrow(d2) * 0.8):nrow(d2), ]
      sum(diff(ss$time) * (head(ss$Cp, -1) + tail(ss$Cp, -1)) / 2)
    }, numeric(1))

    plot_ly(x = fracs, y = aucs, type = "scatter", mode = "lines+markers",
            line = list(color = "#2980b9"), marker = list(color = "#c0392b", size = 8)) |>
      layout(xaxis = list(title = "Relative Dose Fraction"),
             yaxis = list(title = "AUC (SS)"),
             title = list(text = "Dose Proportionality", font = list(size = 13)),
             margin = list(t = 40, b = 50, l = 60, r = 20))
  })

  output$pk_organ_impact <- renderPlotly({
    params  <- drug_pk_params(input$pk_drug, 1.0, "Normal", "Normal")
    combos  <- list(
      c("Normal",             "Normal"),
      c("Mild (eGFR 60-89)",  "Normal"),
      c("Moderate (eGFR 30-59)", "Normal"),
      c("Severe (eGFR <30)",  "Normal"),
      c("Normal",             "Child-Pugh A"),
      c("Normal",             "Child-Pugh B"),
      c("Normal",             "Child-Pugh C")
    )
    labels <- vapply(combos, function(x) paste(x[1], "/", x[2]), character(1))

    aucs <- vapply(combos, function(x) {
      p2 <- drug_pk_params(input$pk_drug, 1.0, x[1], x[2])
      d2 <- sim_pk_1cpt(p2$dose, p2$ka, p2$Vd, p2$CL, p2$F, p2$ii, 7, 168, 0.5)
      ss <- d2[round(nrow(d2) * 0.8):nrow(d2), ]
      sum(diff(ss$time) * (head(ss$Cp, -1) + tail(ss$Cp, -1)) / 2)
    }, numeric(1))

    auc_ref <- aucs[1]
    pct_change <- round((aucs / auc_ref - 1) * 100, 1)
    bar_col    <- ifelse(pct_change > 0, "#c0392b", "#27ae60")

    plot_ly(x = labels, y = pct_change, type = "bar",
            marker = list(color = bar_col)) |>
      layout(
        xaxis = list(title = "", tickangle = -30),
        yaxis = list(title = "AUC change vs Normal (%)"),
        title = list(text = "Organ Function Impact", font = list(size = 13)),
        margin = list(t = 40, b = 110, l = 60, r = 20)
      )
  })

  # ===========================================================================
  # TAB 3: Tumor Dynamics & PD
  # ===========================================================================

  tumor_data <- reactive({
    sim_tumor(
      regimen     = input$td_regimen,
      weeks       = input$td_weeks,
      init_diam_cm = input$td_init_diam,
      esr1_mut    = input$td_esr1,
      pik3ca_mut  = input$td_pik3ca,
      baseline_ki67 = input$td_ki67
    )
  })

  output$tumor_vol_plot <- renderPlotly({
    df <- tumor_data()
    plot_ly(df, x = ~time_wk) |>
      add_lines(y = ~diam_cm, line = list(color = "#c0392b", width = 2),
                name = "Tumor Diameter (cm)") |>
      add_lines(y = ~sqrt(vol_cm3 / pi) * 2, line = list(color = "#2980b9",
                width = 1, dash = "dash"), name = "RECIST diameter proxy") |>
      layout(
        xaxis  = list(title = "Time (weeks)"),
        yaxis  = list(title = "Diameter (cm)"),
        legend = list(orientation = "h"),
        shapes = list(
          list(type = "line", x0 = 0, x1 = max(df$time_wk),
               y0 = input$td_init_diam * 0.7, y1 = input$td_init_diam * 0.7,
               line = list(color = "gray", dash = "dot", width = 1)),
          list(type = "line", x0 = 0, x1 = max(df$time_wk),
               y0 = input$td_init_diam * 1.2, y1 = input$td_init_diam * 1.2,
               line = list(color = "#e74c3c", dash = "dot", width = 1))
        ),
        annotations = list(
          list(x = max(df$time_wk) * 0.95, y = input$td_init_diam * 0.7,
               text = "PR threshold (−30%)", showarrow = FALSE,
               font = list(size = 10, color = "gray")),
          list(x = max(df$time_wk) * 0.95, y = input$td_init_diam * 1.2,
               text = "PD threshold (+20%)", showarrow = FALSE,
               font = list(size = 10, color = "#e74c3c"))
        ),
        margin = list(t = 20, b = 50, l = 60, r = 20)
      )
  })

  output$ki67_plot <- renderPlotly({
    df <- tumor_data()
    plot_ly(df, x = ~time_wk, y = ~ki67, type = "scatter", mode = "lines",
            line = list(color = "#8e44ad", width = 2)) |>
      layout(
        xaxis = list(title = "Time (weeks)"),
        yaxis = list(title = "Ki-67 (%)", range = c(0, 85)),
        margin = list(t = 20, b = 50, l = 50, r = 20)
      )
  })

  output$cdk_plot <- renderPlotly({
    df <- tumor_data()
    plot_ly(df, x = ~time_wk, y = ~cdk_inhib, type = "scatter", mode = "lines",
            fill = "tozeroy", fillcolor = "rgba(41,128,185,0.2)",
            line = list(color = "#2980b9", width = 2)) |>
      layout(
        xaxis = list(title = "Time (weeks)"),
        yaxis = list(title = "CDK4/6 Inhibition (%)", range = c(0, 100)),
        margin = list(t = 20, b = 50, l = 50, r = 20)
      )
  })

  output$pd_table <- renderTable({
    df <- tumor_data()
    last <- tail(df, 1)
    init_vol <- df$vol_cm3[1]
    min_vol  <- min(df$vol_cm3)
    nadir_wk <- df$time_wk[which.min(df$vol_cm3)]

    pct_change <- round((last$diam_cm - input$td_init_diam) / input$td_init_diam * 100, 1)
    pct_nadir  <- round((sqrt(min_vol / pi) * 2 - input$td_init_diam) / input$td_init_diam * 100, 1)

    data.frame(
      Metric = c("Final diameter (cm)", "Diameter change (%)",
                 "Min diameter (cm)", "Min change (%)",
                 "Nadir at week", "Final Ki-67 (%)"),
      Value  = c(round(last$diam_cm, 2), pct_change,
                 round(sqrt(min_vol / pi) * 2, 2), pct_nadir,
                 round(nadir_wk, 1), round(last$ki67, 1))
    )
  }, striped = TRUE, hover = TRUE, bordered = TRUE, align = "lr")

  output$response_badge <- renderUI({
    df <- tumor_data()
    pct <- (tail(df$diam_cm, 1) - input$td_init_diam) / input$td_init_diam * 100
    resp <- if (pct <= -30) list("Partial/Complete Response", "green")
            else if (pct >= 20)  list("Progressive Disease",     "red")
            else                  list("Stable Disease",           "orange")
    div(style = paste0("background:", resp[[2]], "; color:#fff; padding:8px 14px;",
                       " border-radius:4px; text-align:center; font-weight:bold;"),
        resp[[1]])
  })

  # ===========================================================================
  # TAB 4: Immune Microenvironment
  # ===========================================================================

  immune_result <- eventReactive(input$imm_run, {
    sim_immune(input$imm_drug, input$imm_til, input$imm_pdl1,
               input$imm_weeks, input$imm_brca, input$imm_subtype)
  }, ignoreNULL = FALSE)

  output$immune_plot <- renderPlotly({
    res <- immune_result()
    df  <- res$traj

    plot_ly(df, x = ~week) |>
      add_lines(y = ~TIL,  name = "CD8+ TIL (%)",
                line = list(color = "#27ae60", width = 2)) |>
      add_lines(y = ~Treg, name = "Treg (%)",
                line = list(color = "#e74c3c", width = 2, dash = "dash")) |>
      add_lines(y = ~Ag / 2, name = "Antigen load (scaled)",
                line = list(color = "#8e44ad", width = 2, dash = "dot")) |>
      layout(
        xaxis  = list(title = "Time (weeks)"),
        yaxis  = list(title = "Relative abundance (%)"),
        legend = list(orientation = "h"),
        margin = list(t = 20, b = 50, l = 60, r = 20)
      )
  })

  output$pcr_box <- renderValueBox({
    res <- immune_result()
    col <- if (res$pcr_prob >= 40) "green" else if (res$pcr_prob >= 20) "yellow" else "red"
    valueBox(paste0(res$pcr_prob, "%"), "Estimated pCR Probability",
             icon = icon("bullseye"), color = col)
  })
  output$til_box <- renderValueBox({
    res <- immune_result()
    df  <- res$traj
    valueBox(round(tail(df$TIL, 1), 1), "Final CD8+ TIL (%)",
             icon = icon("shield-halved"), color = "green")
  })
  output$treg_box <- renderValueBox({
    res <- immune_result()
    df  <- res$traj
    valueBox(round(tail(df$Treg, 1), 1), "Final Treg (%)",
             icon = icon("circle-minus"), color = "orange")
  })

  output$pcr_pdl1_plot <- renderPlotly({
    pdl1_vals <- seq(0, 100, by = 5)
    pcr_probs <- vapply(pdl1_vals, function(p) {
      r <- sim_immune(input$imm_drug, input$imm_til, p,
                      input$imm_weeks, input$imm_brca, input$imm_subtype)
      r$pcr_prob
    }, numeric(1))

    plot_ly(x = pdl1_vals, y = pcr_probs, type = "scatter", mode = "lines+markers",
            line = list(color = "#c0392b", width = 2),
            marker = list(color = "#c0392b", size = 6)) |>
      layout(
        xaxis = list(title = "PD-L1 CPS Score"),
        yaxis = list(title = "Estimated pCR (%)", range = c(0, 100)),
        title = list(text = "pCR vs PD-L1 CPS", font = list(size = 13)),
        margin = list(t = 40, b = 50, l = 60, r = 20)
      )
  })

  output$trial_bench_plot <- renderPlotly({
    # Published pCR benchmarks
    trials <- data.frame(
      trial = c("KEYNOTE-522\n(pembro+chemo,TNBC)",
                "IMpassion031\n(atezo+chemo,TNBC)",
                "GeparNuevo\n(durva+chemo,TNBC)",
                "NeoSphere\n(pertuz+trast,HER2+)",
                "I-SPY2\n(pembro,HER2-/HR+)"),
      pcr   = c(64.8, 57.6, 53.4, 45.8, 34.0),
      lower = c(59.0, 49.5, 44.0, 36.0, 26.0),
      upper = c(70.6, 65.7, 62.8, 55.7, 42.0)
    )
    model_pcr <- immune_result()$pcr_prob

    plot_ly(trials, x = ~trial, y = ~pcr, type = "bar",
            error_y = list(type = "data",
                           array    = trials$upper - trials$pcr,
                           arrayminus = trials$pcr - trials$lower,
                           color = "gray"),
            marker = list(color = "#2980b9"), name = "Published pCR") |>
      add_trace(x = "Model\nPrediction", y = model_pcr, type = "bar",
                marker = list(color = "#c0392b"), name = "Model") |>
      layout(
        xaxis  = list(title = ""),
        yaxis  = list(title = "pCR Rate (%)", range = c(0, 85)),
        legend = list(orientation = "h"),
        margin = list(t = 20, b = 80, l = 60, r = 20)
      )
  })

  # ===========================================================================
  # TAB 5: Clinical Endpoints
  # ===========================================================================

  endpoint_data <- eventReactive(input$clin_run, {
    regs    <- input$clin_regimens
    months  <- input$clin_months
    t_vec   <- seq(0, months, by = 0.5)

    pfs_list <- lapply(regs, function(r) {
      p <- pfs_params[[r]]
      s <- weibull_surv(t_vec, p$median, p$shape)
      data.frame(time = t_vec, surv = s, regimen = r)
    })
    pfs_df <- do.call(rbind, pfs_list)
    list(pfs = pfs_df, regs = regs, months = months)
  }, ignoreNULL = FALSE)

  output$pfs_plot <- renderPlotly({
    d    <- endpoint_data()
    pfs  <- d$pfs
    regs <- d$regs

    palette <- c("#c0392b","#2980b9","#27ae60","#8e44ad",
                 "#e67e22","#16a085","#2c3e50")
    p <- plot_ly()
    for (i in seq_along(regs)) {
      sub <- pfs[pfs$regimen == regs[i], ]
      p   <- add_lines(p, data = sub, x = ~time, y = ~surv * 100,
                       name = regs[i],
                       line = list(color = palette[i %% length(palette) + 1], width = 2))
    }
    p |>
      layout(
        xaxis  = list(title = "Time (months)"),
        yaxis  = list(title = "PFS Probability (%)", range = c(0, 102)),
        shapes = list(
          list(type = "line", x0 = 0, x1 = d$months,
               y0 = 50, y1 = 50,
               line = list(color = "gray", dash = "dot", width = 1))
        ),
        legend = list(orientation = "v", x = 1.01, y = 1),
        margin = list(t = 20, b = 50, l = 60, r = 160)
      )
  })

  output$orr_plot <- renderPlotly({
    regs <- input$clin_regimens
    orrs <- vapply(regs, function(r) orr_lookup[r], numeric(1))
    df   <- data.frame(regimen = regs, orr = orrs)
    df   <- df[order(df$orr, decreasing = TRUE), ]

    plot_ly(df, x = ~orr, y = ~reorder(regimen, orr), type = "bar",
            orientation = "h",
            marker = list(color = "#c0392b")) |>
      layout(
        xaxis  = list(title = "ORR (%)", range = c(0, 80)),
        yaxis  = list(title = ""),
        margin = list(t = 20, b = 50, l = 200, r = 40)
      )
  })

  output$endpoint_table <- renderTable({
    regs <- input$clin_regimens
    rows <- lapply(regs, function(r) {
      p        <- pfs_params[[r]]
      med_pfs  <- p$median
      # OS estimate: approx 1.6x PFS for HR+ / 1.3x for TNBC
      mult_os  <- if (grepl("TNBC|Pembro", r)) 1.3 else 1.6
      med_os   <- round(med_pfs * mult_os, 1)
      orr      <- orr_lookup[r]
      # 2-year PFS rate
      pfs_24   <- round(weibull_surv(24, p$median, p$shape) * 100, 1)

      data.frame(
        Regimen      = r,
        `Median PFS (mo)` = med_pfs,
        `2-yr PFS (%)`    = pfs_24,
        `Est. Median OS (mo)` = med_os,
        `ORR (%)` = orr,
        stringsAsFactors = FALSE, check.names = FALSE
      )
    })
    do.call(rbind, rows)
  }, striped = TRUE, hover = TRUE, bordered = TRUE)

  # ===========================================================================
  # TAB 6: Scenario Comparison & Biomarkers
  # ===========================================================================

  scenario_data <- eventReactive(input$sc_run, {
    # Tumor dynamics for baseline vs resistance scenario
    base_reg <- input$sc_compare
    esr1   <- grepl("ESR1", input$sc_resist)
    pik3ca <- grepl("PIK3CA|PTEN", input$sc_resist)

    td_base <- sim_tumor(base_reg, input$sc_weeks, 3, FALSE, FALSE, 30)
    td_res  <- sim_tumor(base_reg, input$sc_weeks, 3, esr1, pik3ca, 30)

    biomarkers <- sim_resistance(
      resist_mut      = input$sc_resist,
      scenario        = input$sc_compare,
      baseline_ctdna_vaf = input$sc_ctdna_vaf,
      baseline_lvef   = input$sc_lvef,
      anthracycline   = input$sc_anthracycline,
      weeks           = input$sc_weeks
    )

    list(base = td_base, res = td_res, bio = biomarkers)
  }, ignoreNULL = FALSE)

  output$resist_plot <- renderPlotly({
    d  <- scenario_data()
    plot_ly() |>
      add_lines(data = d$base, x = ~time_wk, y = ~diam_cm,
                name = paste("No resistance:", input$sc_compare),
                line = list(color = "#27ae60", width = 2)) |>
      add_lines(data = d$res, x = ~time_wk, y = ~diam_cm,
                name = paste(input$sc_resist, "mutation"),
                line = list(color = "#c0392b", width = 2, dash = "dash")) |>
      layout(
        xaxis  = list(title = "Time (weeks)"),
        yaxis  = list(title = "Tumor Diameter (cm)"),
        legend = list(orientation = "h"),
        margin = list(t = 20, b = 50, l = 60, r = 20)
      )
  })

  output$ctdna_plot <- renderPlotly({
    d   <- scenario_data()
    bio <- d$bio
    days <- bio$days

    plot_ly(x = days / 7, y = bio$ctdna, type = "scatter", mode = "lines",
            fill = "tozeroy", fillcolor = "rgba(192,57,43,0.15)",
            line = list(color = "#c0392b", width = 2)) |>
      layout(
        xaxis = list(title = "Time (weeks)"),
        yaxis = list(title = "ctDNA VAF (%)"),
        margin = list(t = 20, b = 50, l = 50, r = 20)
      )
  })

  output$anc_plot <- renderPlotly({
    d   <- scenario_data()
    bio <- d$bio

    anc_df <- data.frame(
      week = bio$days / 7,
      ANC  = bio$anc
    )
    # Color by grade
    col_vec <- ifelse(bio$anc < 0.5, "#c0392b",
               ifelse(bio$anc < 1.0, "#e67e22",
               ifelse(bio$anc < 1.5, "#f1c40f", "#27ae60")))

    plot_ly(anc_df, x = ~week, y = ~ANC, type = "scatter", mode = "lines",
            line = list(color = "#2980b9", width = 1.5)) |>
      add_lines(y = rep(1.0, nrow(anc_df)),
                line = list(color = "#c0392b", dash = "dot", width = 1),
                name = "Grade 3 threshold") |>
      add_lines(y = rep(0.5, nrow(anc_df)),
                line = list(color = "#8e44ad", dash = "dot", width = 1),
                name = "Grade 4 threshold") |>
      layout(
        xaxis  = list(title = "Time (weeks)"),
        yaxis  = list(title = "ANC (×10⁹/L)"),
        legend = list(orientation = "h", y = -0.3),
        margin = list(t = 20, b = 70, l = 50, r = 20)
      )
  })

  output$lvef_plot <- renderPlotly({
    d   <- scenario_data()
    bio <- d$bio

    plot_ly(x = bio$days / 7, y = bio$lvef, type = "scatter", mode = "lines",
            fill = "tozeroy", fillcolor = "rgba(39,174,96,0.1)",
            line = list(color = "#27ae60", width = 2)) |>
      add_lines(x = bio$days / 7, y = rep(53, length(bio$days)),
                line = list(color = "#c0392b", dash = "dot", width = 1),
                name = "LVEF 53% (LLN)") |>
      layout(
        xaxis  = list(title = "Time (weeks)"),
        yaxis  = list(title = "LVEF (%)", range = c(15, 75)),
        legend = list(orientation = "h", y = -0.3),
        margin = list(t = 20, b = 70, l = 50, r = 20)
      )
  })

  output$waterfall_plot <- renderPlotly({
    set.seed(77)
    n   <- 40
    reg <- input$sc_compare
    p   <- pfs_params[[reg]]

    # Simulate patient-level % change from baseline
    # Based on Weibull model and variability
    responses <- rnorm(n,
                       mean = if (p$median > 20) -35 else if (p$median > 12) -20 else -10,
                       sd   = 25)
    responses  <- pmax(-100, pmin(100, responses))

    # Apply resistance penalty
    resist_pen <- switch(input$sc_resist,
      "None"          = 0,
      "ESR1 Y537S"    = 15,
      "ESR1 D538G"    = 12,
      "PIK3CA H1047R" = 10,
      "PTEN loss"     = 18,
      0
    )
    responses <- responses + resist_pen

    df <- data.frame(
      patient  = seq_len(n),
      change   = sort(responses),
      response = cut(sort(responses),
                     breaks = c(-Inf, -30, 20, Inf),
                     labels = c("Response (≤−30%)", "Stable (>−30% to +20%)", "PD (>+20%)"))
    )

    pal_resp <- c("Response (≤−30%)" = "#27ae60",
                  "Stable (>−30% to +20%)" = "#f39c12",
                  "PD (>+20%)"             = "#c0392b")

    plot_ly(df, x = ~patient, y = ~change, type = "bar",
            color = ~response, colors = pal_resp) |>
      layout(
        xaxis  = list(title = "Patient (sorted by response)"),
        yaxis  = list(title = "Tumor Size Change from Baseline (%)",
                      range = c(-105, 105),
                      zeroline = TRUE),
        shapes = list(
          list(type="line", x0=0, x1=n+1, y0=-30, y1=-30,
               line=list(color="#27ae60", dash="dot", width=1)),
          list(type="line", x0=0, x1=n+1, y0= 20, y1= 20,
               line=list(color="#c0392b", dash="dot", width=1))
        ),
        barmode = "overlay",
        legend  = list(orientation = "h", y = -0.25),
        margin  = list(t = 20, b = 80, l = 60, r = 20)
      )
  })

}   # end server

# =============================================================================
# RUN APP
# =============================================================================
shinyApp(ui = ui, server = server)
