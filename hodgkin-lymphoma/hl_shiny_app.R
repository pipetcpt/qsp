# ============================================================
# Hodgkin Lymphoma QSP Shiny Application
# Quantitative Systems Pharmacology Dashboard
# ============================================================
# Tabs:
#   1. Patient Profile & Disease Staging
#   2. Pharmacokinetics (ABVD / BV-AVD / Pembrolizumab)
#   3. Tumor Dynamics & Response
#   4. Clinical Endpoints (PFS/OS/PET)
#   5. Treatment Scenario Comparison
#   6. Biomarkers & Toxicity
# ============================================================

library(shiny)
library(ggplot2)
library(dplyr)
library(tidyr)

# ── optional packages with graceful fallback ──────────────────
has_DT        <- requireNamespace("DT",        quietly = TRUE)
has_scales    <- requireNamespace("scales",    quietly = TRUE)
has_gridExtra <- requireNamespace("gridExtra", quietly = TRUE)

if (has_DT)        library(DT)
if (has_scales)    library(scales)
if (has_gridExtra) library(gridExtra)

# ============================================================
# Helper colour palette
# ============================================================
reg_colors <- c(
  "ABVD"                 = "#2196F3",
  "BV-AVD"               = "#4CAF50",
  "Escalated BEACOPP"    = "#F44336",
  "Pembrolizumab (R/R)"  = "#FF9800"
)

# ============================================================
# PK simulation functions (one-compartment IV bolus or Bateman)
# ============================================================

# Generic two-compartment approximation: biexponential decline
pk_biexp <- function(time, dose, Vd, alpha, beta, f_alpha = 0.7) {
  A <- f_alpha      * dose / Vd
  B <- (1 - f_alpha)* dose / Vd
  conc <- A * exp(-alpha * time) + B * exp(-beta * time)
  pmax(conc, 0)
}

# Antibody-drug conjugate (BV): slow distribution, Fc recycling
pk_adc <- function(time, dose_mg_kg, weight, CL = 0.50, Vc = 3.5, k12 = 0.18, k21 = 0.09) {
  dose_ug <- dose_mg_kg * weight * 1000
  Vd <- Vc * weight  # L
  lambda1 <- ((CL/Vd + k12 + k21) + sqrt((CL/Vd + k12 + k21)^2 - 4*(CL/Vd)*k21)) / 2
  lambda2 <- ((CL/Vd + k12 + k21) - sqrt((CL/Vd + k12 + k21)^2 - 4*(CL/Vd)*k21)) / 2
  A <- (lambda1 - k21) / (lambda1 - lambda2) * dose_ug / Vd
  B <- (k21 - lambda2) / (lambda1 - lambda2) * dose_ug / Vd
  A * exp(-lambda1 * time) + B * exp(-lambda2 * time)
}

# Pembrolizumab: ~3-week dosing, long t1/2 ~25 days
pk_pembro <- function(time, dose_mg, weight, CL_L_day = 0.22, Vd_L = 7.5) {
  dose_ug <- dose_mg * 1000
  Vd_body <- Vd_L * (weight / 70)^0.75
  ke <- CL_L_day / Vd_body
  dose_ug / Vd_body * exp(-ke * time)
}

# ============================================================
# Tumor volume ODE (Gompertz + drug kill)
# ============================================================
sim_tumor <- function(cycles, regimen = "ABVD",
                      stage = "III", ips = 3,
                      dose_scale = 1.0, weight = 70,
                      bulky = FALSE) {
  n_days  <- cycles * 28
  t       <- seq(0, n_days, by = 1)
  V0      <- if (stage %in% c("I","II")) 80 else 180   # cm³ baseline
  if (bulky) V0 <- V0 * 1.6
  K       <- 600     # carrying capacity cm³
  lambda  <- 0.012   # Gompertz growth

  # Drug kill parameters per regimen
  kill_params <- list(
    "ABVD"              = list(delta_max = 0.045, IC50 = 0.5, gamma = 1.2),
    "BV-AVD"            = list(delta_max = 0.060, IC50 = 0.4, gamma = 1.3),
    "Escalated BEACOPP" = list(delta_max = 0.075, IC50 = 0.3, gamma = 1.4),
    "Pembrolizumab (R/R)"= list(delta_max = 0.035, IC50 = 0.6, gamma = 1.1)
  )
  p <- kill_params[[regimen]]

  # IPS penalty (higher IPS → slower response)
  ips_factor <- 1 - 0.05 * pmin(ips, 7)

  V <- numeric(length(t))
  V[1] <- V0

  for (i in seq_along(t)[-1]) {
    day <- t[i]
    # Cycle day (0-based within cycle)
    cyc_day <- day %% 28

    # Drug effect: peaks at day 1 and 15 of each cycle (ABVD/BV-AVD schedule)
    drug_effect_raw <- if (cyc_day <= 2) {
      exp(-0.8 * cyc_day)
    } else if (cyc_day >= 14 && cyc_day <= 16) {
      exp(-0.8 * (cyc_day - 14))
    } else {
      0.05 * exp(-0.12 * pmin(cyc_day, 28 - cyc_day))
    }

    drug_conc_norm <- drug_effect_raw * dose_scale * ips_factor
    delta <- p$delta_max * (drug_conc_norm^p$gamma) / (p$IC50^p$gamma + drug_conc_norm^p$gamma)

    # Gompertz growth
    growth <- lambda * log(K / pmax(V[i-1], 0.001)) * V[i-1]
    death  <- delta * V[i-1]

    V[i] <- pmax(V[i-1] + growth - death, 0.001)
  }

  data.frame(day = t, volume_cm3 = V, regimen = regimen)
}

# ============================================================
# ANC Friberg model (myelosuppression)
# ============================================================
sim_anc <- function(cycles, regimen = "ABVD", dose_scale = 1.0, weight = 70) {
  n_days <- cycles * 28
  t      <- seq(0, n_days, by = 0.5)

  ANC0    <- 5.0   # 10^9/L baseline
  ktr     <- 0.15  # transit rate (days^-1)
  gamma   <- 0.17  # feedback exponent

  # Myelosuppression slope (Emax model)
  slope_map <- c("ABVD"=0.8, "BV-AVD"=0.9, "Escalated BEACOPP"=1.4, "Pembrolizumab (R/R)"=0.2)
  slope <- slope_map[regimen] * dose_scale

  n_transit <- 4
  # States: prol, transit1, transit2, transit3, circ
  state <- c(ANC0, ANC0, ANC0, ANC0, ANC0)

  anc_out <- numeric(length(t))
  anc_out[1] <- ANC0

  dt <- 0.5
  for (i in seq_along(t)[-1]) {
    day     <- t[i]
    cyc_day <- day %% 28

    # Drug exposure pulse
    drug_exp <- if (cyc_day <= 1 || (cyc_day >= 14 && cyc_day <= 15)) {
      slope * exp(-0.5 * pmin(cyc_day, abs(cyc_day - 14)))
    } else 0

    E_drug <- drug_exp

    feedback <- (ANC0 / pmax(state[5], 0.01))^gamma

    dprol <- ktr * state[1] * (feedback * (1 - E_drug) - 1)
    dt1   <- ktr * (state[1] - state[2])
    dt2   <- ktr * (state[2] - state[3])
    dt3   <- ktr * (state[3] - state[4])
    dcirc <- ktr * (state[4] - state[5])

    state <- pmax(state + dt * c(dprol, dt1, dt2, dt3, dcirc), 0.01)
    anc_out[i] <- state[5]
  }

  data.frame(day = t, ANC = anc_out, regimen = regimen)
}

# ============================================================
# PFS survival curve (parametric Weibull)
# ============================================================
sim_pfs <- function(t_years, regimen = "ABVD", stage = "III", ips = 3) {
  # Weibull scale (lambda) and shape (k) calibrated to trial data
  params <- list(
    "ABVD"              = list(lambda = 1/7, k = 1.3),
    "BV-AVD"            = list(lambda = 1/10, k = 1.4),
    "Escalated BEACOPP" = list(lambda = 1/9, k = 1.5),
    "Pembrolizumab (R/R)"= list(lambda = 1/4, k = 1.1)
  )
  p <- params[[regimen]]

  # Adjust for IPS
  lambda_adj <- p$lambda * (1 + 0.04 * pmax(ips - 2, 0))

  # Weibull survival: S(t) = exp(-(lambda*t)^k)
  S <- exp(-(lambda_adj * t_years)^p$k)
  pmin(pmax(S, 0), 1)
}

# ============================================================
# LDH normalization (exponential decay)
# ============================================================
sim_ldh <- function(t_days, ldh_baseline, regimen = "ABVD", dose_scale = 1.0) {
  rate_map <- c("ABVD"=0.04, "BV-AVD"=0.05, "Escalated BEACOPP"=0.06, "Pembrolizumab (R/R)"=0.03)
  rate <- rate_map[regimen] * dose_scale
  ldh_normal <- 250  # ULN U/L

  ldh <- ldh_normal + (ldh_baseline - ldh_normal) * exp(-rate * t_days)
  pmax(ldh, ldh_normal * 0.6)
}

# ============================================================
# TARC (CCL17) biomarker dynamics
# ============================================================
sim_tarc <- function(t_days, baseline_tarc = 5000, regimen = "ABVD", dose_scale = 1.0) {
  rate_map <- c("ABVD"=0.05, "BV-AVD"=0.07, "Escalated BEACOPP"=0.06, "Pembrolizumab (R/R)"=0.04)
  rate <- rate_map[regimen] * dose_scale
  normal_tarc <- 450  # pg/mL upper normal

  tarc <- normal_tarc + (baseline_tarc - normal_tarc) * exp(-rate * t_days)
  pmax(tarc, normal_tarc * 0.3)
}

# ============================================================
# Waterfall plot data generator
# ============================================================
gen_waterfall <- function(regimen, n_patients = 30, seed = 42) {
  set.seed(seed)

  # Response distribution based on clinical trial data
  resp_params <- list(
    "ABVD"              = list(mean_change = -55, sd = 28),
    "BV-AVD"            = list(mean_change = -65, sd = 25),
    "Escalated BEACOPP" = list(mean_change = -70, sd = 22),
    "Pembrolizumab (R/R)"= list(mean_change = -50, sd = 35)
  )
  p <- resp_params[[regimen]]

  changes <- rnorm(n_patients, mean = p$mean_change, sd = p$sd)
  changes <- pmax(pmin(changes, 150), -100)  # clamp to realistic range

  data.frame(
    patient_id   = seq_len(n_patients),
    pct_change   = changes,
    response     = case_when(
      changes <= -50 ~ "CR/PR",
      changes <= 0   ~ "SD",
      TRUE           ~ "PD"
    ),
    regimen      = regimen
  ) |> arrange(pct_change)
}

# ============================================================
# Treatment comparison table
# ============================================================
treatment_table <- data.frame(
  Regimen            = c("ABVD", "BV-AVD", "Escalated BEACOPP", "Pembrolizumab (R/R)"),
  Setting            = c("1L Advanced", "1L Advanced", "1L Advanced (unfav.)", "R/R"),
  ORR_pct            = c(76, 86, 90, 69),
  CR_rate_pct        = c(73, 83, 88, 22),
  PFS_5yr_pct        = c(60, 82, 70, 27),
  OS_5yr_pct         = c(82, 87, 85, 59),
  FN_risk_pct        = c(18, 19, 45, 5),
  Grade3_4_neuro_pct = c(6, 25, 8, 3),
  Pulm_toxicity_pct  = c(12, 3, 5, 2),
  Key_Trial          = c("HD10/HD14/ECHELON-1", "ECHELON-1", "HD18/GHSG", "KEYNOTE-087"),
  stringsAsFactors   = FALSE
)

# ============================================================
# UI
# ============================================================
ui <- fluidPage(
  title = "Hodgkin Lymphoma QSP Dashboard",

  tags$head(
    tags$style(HTML("
      body { font-family: 'Segoe UI', Arial, sans-serif; background: #f5f7fa; }
      .navbar { background-color: #1a237e; }
      .well { background-color: #fff; border: 1px solid #e0e0e0; box-shadow: 0 1px 3px rgba(0,0,0,.1); }
      h3 { color: #1a237e; border-bottom: 2px solid #e8eaf6; padding-bottom: 6px; }
      h4 { color: #283593; }
      .metric-box { background:#e8eaf6; border-radius:8px; padding:12px; text-align:center; margin:4px; }
      .metric-box .value { font-size:1.8em; font-weight:bold; color:#1a237e; }
      .metric-box .label { font-size:.85em; color:#555; }
      .tab-content { padding-top: 10px; }
    "))
  ),

  titlePanel(
    div(
      style = "background: linear-gradient(135deg,#1a237e,#283593); color:white; padding:16px 24px; border-radius:8px; margin-bottom:16px;",
      h2("Hodgkin Lymphoma — QSP Interactive Dashboard", style="margin:0; font-size:1.5em;"),
      p("Mechanistic PK/PD · Tumor Dynamics · Clinical Endpoints · Biomarker Simulation", style="margin:4px 0 0; opacity:.85; font-size:.9em;")
    )
  ),

  tabsetPanel(id = "main_tabs", type = "tabs",

    # ──────────────────────────────────────────────────────────
    # TAB 1: Patient Profile & Disease Staging
    # ──────────────────────────────────────────────────────────
    tabPanel("Patient Profile & Staging",
      sidebarLayout(
        sidebarPanel(width = 3,
          h4("Patient Demographics"),
          sliderInput("pt_age",   "Age (years)",    min=18, max=80, value=32),
          sliderInput("pt_wt",    "Weight (kg)",    min=40, max=130, value=70),
          selectInput("pt_sex",   "Sex",            choices=c("Male","Female"), selected="Male"),
          hr(),
          h4("Disease Characteristics"),
          selectInput("hl_stage", "Ann Arbor Stage",
                      choices=c("I","II","III","IV"), selected="III"),
          selectInput("hl_subtype","Histologic Subtype",
                      choices=c("NSHL (Nodular Sclerosis)"="NSHL",
                                "MCHL (Mixed Cellularity)"="MCHL",
                                "LRHL (Lymphocyte-Rich)"="LRHL",
                                "LDHL (Lymphocyte-Depleted)"="LDHL"), selected="NSHL"),
          checkboxInput("bulky",    "Bulky Disease (≥10 cm)",    value=FALSE),
          checkboxInput("ebv_pos",  "EBV Positive",              value=FALSE),
          checkboxInput("b_syms",   "B-Symptoms Present",        value=TRUE),
          hr(),
          h4("IPS Factors (check if present)"),
          checkboxInput("ips_alb",  "Albumin < 4 g/dL",   value=FALSE),
          checkboxInput("ips_hgb",  "Hgb < 10.5 g/dL",    value=FALSE),
          checkboxInput("ips_male", "Male sex",            value=TRUE),
          checkboxInput("ips_age",  "Age ≥ 45",            value=FALSE),
          checkboxInput("ips_stg4", "Stage IV",            value=FALSE),
          checkboxInput("ips_wbc",  "WBC ≥ 15,000/μL",    value=FALSE),
          checkboxInput("ips_lymp", "Lymphocyte < 8%",     value=FALSE),
          hr(),
          sliderInput("ldh_base", "Baseline LDH (U/L)", min=100, max=1500, value=450, step=10)
        ),
        mainPanel(width = 9,
          fluidRow(
            column(4, div(class="metric-box", div(class="value", textOutput("ips_score")), div(class="label","IPS Score"))),
            column(4, div(class="metric-box", div(class="value", textOutput("risk_group")), div(class="label","Risk Group"))),
            column(4, div(class="metric-box", div(class="value", textOutput("ldh_uln")), div(class="label","LDH ×ULN")))
          ),
          br(),
          fluidRow(
            column(6, plotOutput("spider_ips", height="350px")),
            column(6, plotOutput("stage_diagram", height="350px"))
          ),
          br(),
          fluidRow(
            column(12, plotOutput("subtype_pie", height="280px"))
          )
        )
      )
    ),

    # ──────────────────────────────────────────────────────────
    # TAB 2: Pharmacokinetics
    # ──────────────────────────────────────────────────────────
    tabPanel("Pharmacokinetics",
      sidebarLayout(
        sidebarPanel(width = 3,
          h4("Regimen Selection"),
          selectInput("pk_regimen", "Treatment Regimen",
                      choices = c("ABVD","BV-AVD","Pembrolizumab (R/R)"), selected="ABVD"),
          hr(),
          h4("ABVD Doses (% standard)"),
          sliderInput("dose_dox",  "Doxorubicin (25 mg/m²)",   min=50, max=120, value=100),
          sliderInput("dose_bleo", "Bleomycin (10 U/m²)",      min=50, max=120, value=100),
          sliderInput("dose_vinb", "Vinblastine (6 mg/m²)",    min=50, max=120, value=100),
          sliderInput("dose_dtic", "Dacarbazine (375 mg/m²)",  min=50, max=120, value=100),
          hr(),
          h4("BV + MMAE"),
          sliderInput("dose_bv",   "BV dose (mg/kg)",          min=0.5, max=2.4, value=1.8, step=0.1),
          hr(),
          h4("Pembrolizumab"),
          sliderInput("dose_pembro","Pembro dose (mg)",         min=100, max=400, value=200, step=50),
          hr(),
          sliderInput("pk_wt",     "Patient weight (kg)",      min=40, max=130, value=70)
        ),
        mainPanel(width = 9,
          fluidRow(
            column(4, div(class="metric-box", div(class="value", textOutput("pk_cmax")),  div(class="label","Cₐₐₐ (primary drug)"))),
            column(4, div(class="metric-box", div(class="value", textOutput("pk_auc")),   div(class="label","AUC₀₋₁₄ days·μg/mL"))),
            column(4, div(class="metric-box", div(class="value", textOutput("pk_thalf")), div(class="label","t₁₂ (hours)")))
          ),
          br(),
          plotOutput("pk_conc_plot", height = "380px"),
          br(),
          plotOutput("pk_multiDrug", height = "300px")
        )
      )
    ),

    # ──────────────────────────────────────────────────────────
    # TAB 3: Tumor Dynamics & Response
    # ──────────────────────────────────────────────────────────
    tabPanel("Tumor Dynamics & Response",
      sidebarLayout(
        sidebarPanel(width = 3,
          h4("Simulation Settings"),
          selectInput("td_regimen", "Regimen",
                      choices=c("ABVD","BV-AVD","Escalated BEACOPP","Pembrolizumab (R/R)"), selected="ABVD"),
          selectInput("td_stage",  "Stage",  choices=c("I","II","III","IV"), selected="III"),
          sliderInput("td_ips",    "IPS Score",  min=0, max=7, value=3, step=1),
          sliderInput("td_cycles", "No. Cycles", min=2, max=8, value=6, step=1),
          sliderInput("td_dose",   "Dose Intensity (%)", min=50, max=120, value=100),
          checkboxInput("td_bulky","Bulky Disease", value=FALSE),
          hr(),
          h4("Waterfall Plot"),
          sliderInput("wf_n", "No. Patients (simulation)", min=10, max=60, value=30),
          numericInput("wf_seed", "Random Seed", value=42, min=1, max=999)
        ),
        mainPanel(width = 9,
          plotOutput("tumor_vol_plot", height = "320px"),
          br(),
          fluidRow(
            column(6, plotOutput("waterfall_plot", height="280px")),
            column(6, plotOutput("ldh_curve",       height="280px"))
          ),
          br(),
          fluidRow(
            column(6, div(class="metric-box",
              div(class="value", textOutput("deauville_score")),
              div(class="label","Predicted Deauville Score at Cycle 2 PET")
            )),
            column(6, div(class="metric-box",
              div(class="value", textOutput("interim_pet_resp")),
              div(class="label","Interim PET-2 Response")
            ))
          )
        )
      )
    ),

    # ──────────────────────────────────────────────────────────
    # TAB 4: Clinical Endpoints
    # ──────────────────────────────────────────────────────────
    tabPanel("Clinical Endpoints",
      sidebarLayout(
        sidebarPanel(width = 3,
          h4("Survival Simulation"),
          checkboxGroupInput("ep_regimens", "Show Regimens",
                             choices  = c("ABVD","BV-AVD","Escalated BEACOPP","Pembrolizumab (R/R)"),
                             selected = c("ABVD","BV-AVD")),
          selectInput("ep_stage", "Stage",  choices=c("I","II","III","IV"), selected="III"),
          sliderInput("ep_ips",   "IPS Score", min=0, max=7, value=3, step=1),
          hr(),
          h4("PET-CT Scenarios"),
          selectInput("ep_pet2_result", "Interim PET-2 Result",
                      choices=c("Negative (Deauville 1-2)"="neg",
                                "Positive (Deauville 3-5)"="pos"), selected="neg"),
          hr(),
          h4("Response Distribution"),
          selectInput("ep_resp_reg", "Regimen for Response Pie",
                      choices=c("ABVD","BV-AVD","Escalated BEACOPP","Pembrolizumab (R/R)"), selected="ABVD")
        ),
        mainPanel(width = 9,
          fluidRow(
            column(6, plotOutput("pfs_curve", height="320px")),
            column(6, plotOutput("os_curve",  height="320px"))
          ),
          br(),
          fluidRow(
            column(6, plotOutput("response_pie", height="260px")),
            column(6, plotOutput("pet2_plot",    height="260px"))
          )
        )
      )
    ),

    # ──────────────────────────────────────────────────────────
    # TAB 5: Treatment Scenario Comparison
    # ──────────────────────────────────────────────────────────
    tabPanel("Treatment Scenario Comparison",
      fluidPage(
        h3("Side-by-Side Regimen Comparison"),
        fluidRow(
          column(3,
            h4("Simulation Parameters"),
            sliderInput("sc_ips",    "IPS Score",         min=0, max=7,   value=3, step=1),
            selectInput("sc_stage",  "Stage",             choices=c("I","II","III","IV"), selected="III"),
            sliderInput("sc_cycles", "Cycles to Simulate",min=2, max=8,   value=6, step=1),
            sliderInput("sc_dose",   "Dose Intensity (%)",min=60, max=120, value=100),
            checkboxInput("sc_bulky","Bulky Disease",     value=FALSE)
          ),
          column(9,
            plotOutput("sc_tumor_overlay", height="320px"),
            br(),
            plotOutput("sc_anc_overlay",   height="260px")
          )
        ),
        br(),
        h3("Efficacy & Safety Summary"),
        if (has_DT) DT::dataTableOutput("regimen_table") else tableOutput("regimen_table_base"),
        br(),
        plotOutput("efficacy_radar", height="320px")
      )
    ),

    # ──────────────────────────────────────────────────────────
    # TAB 6: Biomarkers & Toxicity
    # ──────────────────────────────────────────────────────────
    tabPanel("Biomarkers & Toxicity",
      sidebarLayout(
        sidebarPanel(width = 3,
          h4("Biomarker Settings"),
          selectInput("bm_regimen","Regimen",
                      choices=c("ABVD","BV-AVD","Escalated BEACOPP","Pembrolizumab (R/R)"), selected="ABVD"),
          sliderInput("bm_cycles","Cycles",    min=2, max=8, value=6, step=1),
          sliderInput("bm_dose",  "Dose (%)",  min=50, max=120, value=100),
          sliderInput("bm_wt",    "Weight (kg)",min=40, max=130, value=70),
          hr(),
          h4("Baseline Biomarkers"),
          sliderInput("bm_ldh",   "LDH (U/L)",  min=100, max=1500, value=450, step=10),
          sliderInput("bm_tarc",  "TARC (pg/mL)", min=200, max=15000, value=5000, step=100),
          sliderInput("bm_anc0",  "Baseline ANC (×10⁹/L)", min=2, max=15, value=5, step=0.5),
          sliderInput("bm_plt0",  "Baseline PLT (×10⁹/L)", min=100, max=600, value=220, step=10),
          hr(),
          h4("CD30 Expression"),
          sliderInput("bm_cd30",  "CD30 Expression (%)", min=10, max=100, value=85, step=5)
        ),
        mainPanel(width = 9,
          fluidRow(
            column(6, plotOutput("bm_ldh_plot",  height="260px")),
            column(6, plotOutput("bm_tarc_plot", height="260px"))
          ),
          br(),
          fluidRow(
            column(6, plotOutput("bm_anc_plot",  height="260px")),
            column(6, plotOutput("bm_plt_plot",  height="260px"))
          ),
          br(),
          fluidRow(
            column(6, plotOutput("bm_cd30_plot", height="240px")),
            column(6,
              div(class="metric-box", style="margin-top:30px;",
                div(class="value", textOutput("bsym_resolution")),
                div(class="label","Predicted B-Symptom Resolution (day)")
              ),
              br(),
              div(class="metric-box",
                div(class="value", textOutput("anc_nadir")),
                div(class="label","Predicted ANC Nadir (×10⁹/L)")
              ),
              br(),
              div(class="metric-box",
                div(class="value", textOutput("bm_fn_risk")),
                div(class="label","Febrile Neutropenia Risk (%)")
              )
            )
          )
        )
      )
    )
  ) # end tabsetPanel
)

# ============================================================
# SERVER
# ============================================================
server <- function(input, output, session) {

  # ── IPS calculation ──────────────────────────────────────────
  ips_val <- reactive({
    sum(c(input$ips_alb, input$ips_hgb, input$ips_male,
          input$ips_age,  input$ips_stg4, input$ips_wbc, input$ips_lymp))
  })

  output$ips_score  <- renderText(ips_val())
  output$risk_group <- renderText({
    ips <- ips_val()
    if (ips <= 2) "Favorable" else if (ips <= 4) "Intermediate" else "Unfavorable"
  })
  output$ldh_uln <- renderText(sprintf("%.1f×", input$ldh_base / 250))

  # IPS spider chart
  output$spider_ips <- renderPlot({
    factors <- c("Albumin<4","Hgb<10.5","Male","Age≥45","Stage IV","WBC≥15k","Lymph<8%")
    vals    <- c(input$ips_alb, input$ips_hgb, input$ips_male,
                 input$ips_age, input$ips_stg4, input$ips_wbc, input$ips_lymp) * 1
    df <- data.frame(factor=factors, value=vals, stringsAsFactors=FALSE)
    df$factor <- factor(df$factor, levels=factors)

    ggplot(df, aes(x=factor, y=value, fill=factor(value))) +
      geom_col(width=0.6, show.legend=FALSE) +
      scale_fill_manual(values=c("0"="#90CAF9","1"="#F44336")) +
      scale_y_continuous(breaks=c(0,1), limits=c(0,1.3)) +
      labs(title=sprintf("IPS Factors (Score = %d/7)", ips_val()),
           x=NULL, y="Present (1) / Absent (0)") +
      theme_minimal(base_size=12) +
      theme(axis.text.x=element_text(angle=30, hjust=1, size=9),
            plot.title=element_text(colour="#1a237e", face="bold"))
  })

  # Stage diagram (Ann Arbor schematic)
  output$stage_diagram <- renderPlot({
    stage_num <- match(input$hl_stage, c("I","II","III","IV"))
    body_regions <- data.frame(
      region=c("Cervical\n(above diaphragm)","Mediastinal\n(above diaphragm)",
               "Axillary\n(above diaphragm)","Splenic\n(below diaphragm)",
               "Para-aortic\n(below diaphragm)","Inguinal\n(below diaphragm)"),
      involved=c(1, stage_num>=2, stage_num>=2, stage_num>=3, stage_num>=3, stage_num>=4),
      y=c(3,3,3,1,1,1), x=c(1,2,3,1,2,3),
      stringsAsFactors=FALSE
    )
    ggplot(body_regions, aes(x=x, y=y, fill=factor(involved), label=region)) +
      geom_tile(colour="white", linewidth=2, width=0.85, height=0.85) +
      geom_text(size=3, lineheight=0.9) +
      scale_fill_manual(values=c("0"="#E3F2FD","1"="#F44336"),
                        labels=c("Not Involved","Involved"), name="Status") +
      scale_y_continuous(breaks=c(1,3), labels=c("Below Diaphragm","Above Diaphragm")) +
      labs(title=sprintf("Ann Arbor Stage %s — %s", input$hl_stage, input$hl_subtype),
           x=NULL, y=NULL) +
      theme_minimal(base_size=11) +
      theme(axis.text.x=element_blank(), axis.ticks.x=element_blank(),
            plot.title=element_text(colour="#1a237e", face="bold"),
            legend.position="bottom")
  })

  # Subtype distribution pie
  output$subtype_pie <- renderPlot({
    subtype_prevalence <- data.frame(
      Subtype=c("NSHL","MCHL","LRHL","LDHL"),
      Prevalence=c(65,25,5,5),
      stringsAsFactors=FALSE
    )
    subtype_prevalence$highlight <- subtype_prevalence$Subtype == input$hl_subtype

    ggplot(subtype_prevalence, aes(x="", y=Prevalence, fill=Subtype,
                                   alpha=highlight)) +
      geom_col(width=1, colour="white", linewidth=1.5) +
      coord_polar("y") +
      scale_alpha_manual(values=c("FALSE"=0.5,"TRUE"=1), guide="none") +
      scale_fill_manual(values=c(NSHL="#2196F3",MCHL="#FF9800",LRHL="#4CAF50",LDHL="#9C27B0")) +
      geom_text(aes(label=paste0(Subtype,"\n",Prevalence,"%")),
                position=position_stack(vjust=0.5), size=3.5) +
      labs(title="Classical HL Histologic Subtypes (population prevalence)",
           fill="Subtype") +
      theme_void(base_size=12) +
      theme(plot.title=element_text(colour="#1a237e", face="bold", hjust=0.5))
  })

  # ── TAB 2: PK ────────────────────────────────────────────────
  pk_data <- reactive({
    reg <- input$pk_regimen
    wt  <- input$pk_wt
    t   <- seq(0, 14, by=0.1)  # one cycle half (14 days)

    if (reg == "ABVD") {
      dox_dose  <- 25 * (1.73) * (input$dose_dox/100)   # mg (using BSA≈1.73)
      bleo_dose <- 10 * (1.73) * (input$dose_bleo/100)
      vinb_dose <- 6  * (1.73) * (input$dose_vinb/100)
      dtic_dose <- 375* (1.73) * (input$dose_dtic/100)

      bind_rows(
        data.frame(time=t, conc=pk_biexp(t, dox_dose,  Vd=30*wt/70, alpha=0.35, beta=0.045), drug="Doxorubicin"),
        data.frame(time=t, conc=pk_biexp(t, bleo_dose, Vd=18*wt/70, alpha=0.28, beta=0.035), drug="Bleomycin"),
        data.frame(time=t, conc=pk_biexp(t, vinb_dose, Vd=40*wt/70, alpha=0.45, beta=0.025), drug="Vinblastine"),
        data.frame(time=t, conc=pk_biexp(t, dtic_dose, Vd=22*wt/70, alpha=0.55, beta=0.055), drug="Dacarbazine")
      )
    } else if (reg == "BV-AVD") {
      bv_conc <- pk_adc(t, input$dose_bv, wt) / 1000   # ng/mL → μg/mL
      mmae_conc <- bv_conc * 0.08 * exp(-0.25*t)       # MMAE: ~8% payload release, faster clearance
      bind_rows(
        data.frame(time=t, conc=bv_conc,   drug="BV (antibody)"),
        data.frame(time=t, conc=mmae_conc, drug="MMAE (payload)")
      )
    } else {
      p_conc <- pk_pembro(t, input$dose_pembro, wt)
      data.frame(time=t, conc=p_conc/1000, drug="Pembrolizumab")
    }
  })

  output$pk_cmax <- renderText({
    d <- pk_data()
    sprintf("%.2f μg/mL", max(d$conc[d$drug == d$drug[1]], na.rm=TRUE))
  })
  output$pk_auc <- renderText({
    d <- pk_data() |> filter(drug == drug[1])
    auc <- sum(diff(d$time) * (head(d$conc,-1) + tail(d$conc,-1))/2)
    sprintf("%.1f", auc)
  })
  output$pk_thalf <- renderText({
    reg <- input$pk_regimen
    t12 <- switch(reg,
      "ABVD"               = 30,
      "BV-AVD"             = 108,
      "Pembrolizumab (R/R)"= 600
    )
    sprintf("~%d h", t12)
  })

  output$pk_conc_plot <- renderPlot({
    d <- pk_data()
    ggplot(d, aes(x=time, y=conc, colour=drug)) +
      geom_line(linewidth=1.2) +
      scale_colour_brewer(palette="Set1", name="Drug") +
      scale_y_log10(labels=function(x) sprintf("%.3g", x)) +
      labs(title=paste("PK Concentration Profile —", input$pk_regimen),
           x="Time (days)", y="Concentration (μg/mL, log scale)") +
      theme_minimal(base_size=12) +
      theme(legend.position="right",
            plot.title=element_text(colour="#1a237e", face="bold"))
  })

  output$pk_multiDrug <- renderPlot({
    # PK summary bar chart (Cmax normalised)
    df_summary <- pk_data() |>
      group_by(drug) |>
      summarise(Cmax=max(conc, na.rm=TRUE),
                AUC=sum(diff(time)*(head(conc,-1)+tail(conc,-1))/2), .groups="drop")

    df_long <- df_summary |>
      pivot_longer(c(Cmax,AUC), names_to="metric", values_to="value") |>
      group_by(metric) |>
      mutate(rel_value = value / max(value))

    ggplot(df_long, aes(x=drug, y=rel_value, fill=drug)) +
      geom_col(show.legend=FALSE) +
      facet_wrap(~metric, scales="free_y") +
      scale_fill_brewer(palette="Set1") +
      labs(title="Relative PK Metrics by Drug",
           x=NULL, y="Relative Value (normalised to max)") +
      theme_minimal(base_size=11) +
      theme(axis.text.x=element_text(angle=30, hjust=1),
            plot.title=element_text(colour="#1a237e", face="bold"))
  })

  # ── TAB 3: Tumor Dynamics ────────────────────────────────────
  td_sim <- reactive({
    sim_tumor(input$td_cycles, input$td_regimen,
              input$td_stage, input$td_ips,
              input$td_dose/100, 70, input$td_bulky)
  })

  output$tumor_vol_plot <- renderPlot({
    d    <- td_sim()
    v0   <- d$volume_cm3[1]
    d$pct_change <- (d$volume_cm3 - v0) / v0 * 100

    # Mark cycle boundaries
    cycle_days <- seq(0, max(d$day), by=28)

    ggplot(d, aes(x=day, y=volume_cm3)) +
      geom_vline(xintercept=cycle_days, linetype="dashed", colour="#bdbdbd", linewidth=0.5) +
      geom_line(colour=reg_colors[input$td_regimen], linewidth=1.5) +
      geom_hline(yintercept=v0*0.5, linetype="dotted", colour="#F44336", linewidth=0.8) +
      annotate("text", x=5, y=v0*0.52, label="50% reduction threshold", size=3.2, colour="#F44336") +
      scale_x_continuous(breaks=cycle_days, labels=paste0("C",seq_along(cycle_days)-1)) +
      labs(title=sprintf("Tumor Volume Dynamics — %s (Stage %s, IPS %d)",
                         input$td_regimen, input$td_stage, input$td_ips),
           x="Treatment Cycle", y="Tumor Volume (cm³)") +
      theme_minimal(base_size=12) +
      theme(plot.title=element_text(colour="#1a237e", face="bold"))
  })

  output$waterfall_plot <- renderPlot({
    df <- gen_waterfall(input$td_regimen, input$wf_n, input$wf_seed)
    df$patient_id <- factor(df$patient_id, levels=df$patient_id)

    ggplot(df, aes(x=patient_id, y=pct_change, fill=response)) +
      geom_col(width=0.85) +
      geom_hline(yintercept=c(-50,20), linetype="dashed", colour=c("#2196F3","#F44336")) +
      scale_fill_manual(values=c("CR/PR"="#4CAF50","SD"="#FF9800","PD"="#F44336"), name="Response") +
      labs(title=paste("Waterfall Plot — Interim PET\n", input$td_regimen),
           x="Patient (ranked)", y="Change in Tumor Volume (%)") +
      theme_minimal(base_size=11) +
      theme(axis.text.x=element_blank(), axis.ticks.x=element_blank(),
            plot.title=element_text(colour="#1a237e", face="bold"),
            legend.position="top")
  })

  output$ldh_curve <- renderPlot({
    t    <- seq(0, input$td_cycles*28, by=1)
    ldh  <- sim_ldh(t, input$ldh_base, input$td_regimen, input$td_dose/100)
    df   <- data.frame(day=t, LDH=ldh)

    ggplot(df, aes(x=day, y=LDH)) +
      geom_line(colour="#9C27B0", linewidth=1.3) +
      geom_hline(yintercept=250, linetype="dashed", colour="#F44336") +
      annotate("text", x=5, y=260, label="ULN 250 U/L", size=3, colour="#F44336") +
      labs(title="LDH Normalization Curve",
           x="Day", y="LDH (U/L)") +
      theme_minimal(base_size=11) +
      theme(plot.title=element_text(colour="#1a237e", face="bold"))
  })

  deauville_reactive <- reactive({
    d  <- td_sim()
    v0 <- d$volume_cm3[1]
    # Volume at day 56 (end of cycle 2)
    v56_idx <- which.min(abs(d$day - 56))
    v56 <- d$volume_cm3[v56_idx]
    pct_red <- (v0 - v56) / v0 * 100

    if      (pct_red > 90) 1
    else if (pct_red > 70) 2
    else if (pct_red > 50) 3
    else if (pct_red > 0)  4
    else                   5
  })

  output$deauville_score  <- renderText(deauville_reactive())
  output$interim_pet_resp <- renderText({
    ds <- deauville_reactive()
    if (ds <= 2) "Negative (Excellent Prognosis)" else if (ds==3) "Equivocal" else "Positive (Consider Escalation)"
  })

  # ── TAB 4: Clinical Endpoints ─────────────────────────────────
  output$pfs_curve <- renderPlot({
    t_years <- seq(0, 10, by=0.1)
    reg_list <- input$ep_regimens
    if (length(reg_list) == 0) reg_list <- "ABVD"

    df <- bind_rows(lapply(reg_list, function(r) {
      data.frame(time_years=t_years,
                 PFS=sim_pfs(t_years, r, input$ep_stage, input$ep_ips),
                 regimen=r)
    }))

    ggplot(df, aes(x=time_years, y=PFS*100, colour=regimen)) +
      geom_line(linewidth=1.3) +
      geom_hline(yintercept=50, linetype="dotted", colour="#9E9E9E") +
      scale_colour_manual(values=reg_colors, name="Regimen") +
      scale_y_continuous(limits=c(0,100), breaks=seq(0,100,20)) +
      labs(title=sprintf("Progression-Free Survival (Stage %s, IPS %d)", input$ep_stage, input$ep_ips),
           x="Time (years)", y="PFS (%)") +
      theme_minimal(base_size=12) +
      theme(plot.title=element_text(colour="#1a237e", face="bold"),
            legend.position="bottom")
  })

  output$os_curve <- renderPlot({
    t_years <- seq(0, 10, by=0.1)
    reg_list <- input$ep_regimens
    if (length(reg_list) == 0) reg_list <- "ABVD"

    # OS ~ slightly higher and slower declining than PFS
    os_params <- list(
      "ABVD"              = list(lambda=1/10, k=1.2),
      "BV-AVD"            = list(lambda=1/12, k=1.3),
      "Escalated BEACOPP" = list(lambda=1/11, k=1.3),
      "Pembrolizumab (R/R)"= list(lambda=1/6,  k=1.1)
    )

    df <- bind_rows(lapply(reg_list, function(r) {
      p <- os_params[[r]]
      S_os <- exp(-(p$lambda * t_years)^p$k)
      data.frame(time_years=t_years, OS=pmin(S_os,1)*100, regimen=r)
    }))

    ggplot(df, aes(x=time_years, y=OS, colour=regimen)) +
      geom_line(linewidth=1.3) +
      scale_colour_manual(values=reg_colors, name="Regimen") +
      scale_y_continuous(limits=c(0,100)) +
      labs(title="Overall Survival (simulated)", x="Time (years)", y="OS (%)") +
      theme_minimal(base_size=12) +
      theme(plot.title=element_text(colour="#1a237e", face="bold"),
            legend.position="bottom")
  })

  output$response_pie <- renderPlot({
    reg <- input$ep_resp_reg
    resp_data <- list(
      "ABVD"              = c(CR=73, PR=5, SD=10, PD=12),
      "BV-AVD"            = c(CR=83, PR=3, SD=7,  PD=7),
      "Escalated BEACOPP" = c(CR=88, PR=2, SD=5,  PD=5),
      "Pembrolizumab (R/R)"= c(CR=22, PR=47, SD=16, PD=15)
    )
    rv  <- resp_data[[reg]]
    df  <- data.frame(response=names(rv), pct=as.numeric(rv), stringsAsFactors=FALSE)
    df$response <- factor(df$response, levels=c("CR","PR","SD","PD"))

    ggplot(df, aes(x="", y=pct, fill=response)) +
      geom_col(width=1, colour="white") +
      coord_polar("y") +
      geom_text(aes(label=paste0(response,"\n",pct,"%")),
                position=position_stack(vjust=0.5), size=3.5, fontface="bold") +
      scale_fill_manual(values=c(CR="#4CAF50",PR="#8BC34A",SD="#FF9800",PD="#F44336")) +
      labs(title=sprintf("Best Response — %s\n(ORR = %d%%)", reg, rv["CR"]+rv["PR"])) +
      theme_void(base_size=11) +
      theme(plot.title=element_text(colour="#1a237e", face="bold", hjust=0.5),
            legend.position="none")
  })

  output$pet2_plot <- renderPlot({
    t_years <- seq(0, 8, by=0.1)
    pet2    <- input$ep_pet2_result
    # PET-2 negative → ~90% 5yr PFS; positive → ~60% depending on escalation
    if (pet2 == "neg") {
      S_pfs  <- exp(-(t_years/9)^1.5)
      label  <- "PET-2 Negative (Deauville 1-2)\n~88% 5-yr PFS"
      col    <- "#4CAF50"
    } else {
      S_pfs  <- exp(-(t_years/5)^1.4)
      label  <- "PET-2 Positive (Deauville 3-5)\n~60% 5-yr PFS"
      col    <- "#F44336"
    }
    df <- data.frame(t=t_years, pfs=pmin(S_pfs*100,100))

    ggplot(df, aes(x=t, y=pfs)) +
      geom_line(colour=col, linewidth=1.5) +
      geom_hline(yintercept=50, linetype="dotted", colour="#9E9E9E") +
      annotate("text", x=1, y=95, label=label, hjust=0, colour=col, size=3.5, fontface="bold") +
      labs(title="PET-2 Response Impact on PFS",
           x="Time (years)", y="PFS (%)") +
      scale_y_continuous(limits=c(0,100)) +
      theme_minimal(base_size=11) +
      theme(plot.title=element_text(colour="#1a237e", face="bold"))
  })

  # ── TAB 5: Scenario Comparison ────────────────────────────────
  output$sc_tumor_overlay <- renderPlot({
    regs <- c("ABVD","BV-AVD","Escalated BEACOPP","Pembrolizumab (R/R)")
    df   <- bind_rows(lapply(regs, function(r) {
      sim_tumor(input$sc_cycles, r, input$sc_stage, input$sc_ips,
                input$sc_dose/100, 70, input$sc_bulky)
    }))
    v0_map <- df |> group_by(regimen) |> summarise(v0=first(volume_cm3))
    df     <- df |> left_join(v0_map, by="regimen") |>
                    mutate(pct_change = (volume_cm3 - v0) / v0 * 100)

    cycle_days <- seq(0, input$sc_cycles*28, by=28)

    ggplot(df, aes(x=day, y=pct_change, colour=regimen)) +
      geom_hline(yintercept=0, colour="#9E9E9E") +
      geom_vline(xintercept=cycle_days, linetype="dashed", colour="#e0e0e0") +
      geom_line(linewidth=1.2) +
      scale_colour_manual(values=reg_colors, name="Regimen") +
      scale_x_continuous(breaks=cycle_days, labels=paste0("C",seq_along(cycle_days)-1)) +
      labs(title=sprintf("Tumor Volume Change — All Regimens (Stage %s, IPS %d)", input$sc_stage, input$sc_ips),
           x="Treatment Cycle", y="Tumor Volume Change (%)") +
      theme_minimal(base_size=12) +
      theme(legend.position="bottom",
            plot.title=element_text(colour="#1a237e", face="bold"))
  })

  output$sc_anc_overlay <- renderPlot({
    regs <- c("ABVD","BV-AVD","Escalated BEACOPP","Pembrolizumab (R/R)")
    df   <- bind_rows(lapply(regs, function(r) {
      sim_anc(input$sc_cycles, r, input$sc_dose/100, 70)
    }))

    ggplot(df, aes(x=day, y=ANC, colour=regimen)) +
      geom_line(linewidth=1.0, alpha=0.85) +
      geom_hline(yintercept=0.5, linetype="dashed", colour="#F44336") +
      annotate("text", x=2, y=0.55, label="ANC 0.5 (severe neutropenia)", size=3, colour="#F44336") +
      scale_colour_manual(values=reg_colors, name="Regimen") +
      labs(title="Absolute Neutrophil Count (ANC) — Friberg Model",
           x="Day", y="ANC (×10⁹/L)") +
      theme_minimal(base_size=11) +
      theme(legend.position="bottom",
            plot.title=element_text(colour="#1a237e", face="bold"))
  })

  output$regimen_table <- if (has_DT) {
    DT::renderDataTable({
      colnames(treatment_table) <- c("Regimen","Setting","ORR (%)","CR Rate (%)","5-yr PFS (%)","5-yr OS (%)","FN Risk (%)","G3-4 Neuro (%)","Pulm Tox (%)","Key Trial")
      DT::datatable(treatment_table, options=list(pageLength=5, scrollX=TRUE),
                    rownames=FALSE) |>
        DT::formatStyle("ORR (%)", background=DT::styleColorBar(c(0,100),"#90CAF9")) |>
        DT::formatStyle("5-yr PFS (%)", background=DT::styleColorBar(c(0,100),"#A5D6A7"))
    })
  } else {
    renderTable({
      colnames(treatment_table) <- c("Regimen","Setting","ORR (%)","CR Rate (%)","5-yr PFS (%)","5-yr OS (%)","FN Risk (%)","G3-4 Neuro (%)","Pulm Tox (%)","Key Trial")
      treatment_table
    })
  }

  output$regimen_table_base <- renderTable(treatment_table)

  output$efficacy_radar <- renderPlot({
    metrics <- c("ORR","CR Rate","5-yr PFS","5-yr OS","Safety\n(100-FN%)")
    df_radar <- data.frame(
      regimen = rep(treatment_table$Regimen, each=5),
      metric  = rep(metrics, times=nrow(treatment_table)),
      value   = c(
        rbind(treatment_table$ORR_pct,
              treatment_table$CR_rate_pct,
              treatment_table$PFS_5yr_pct,
              treatment_table$OS_5yr_pct,
              100 - treatment_table$FN_risk_pct)
      ),
      stringsAsFactors=FALSE
    )

    ggplot(df_radar, aes(x=metric, y=value, fill=regimen, group=regimen)) +
      geom_col(position="dodge", width=0.7) +
      scale_fill_manual(values=reg_colors, name="Regimen") +
      scale_y_continuous(limits=c(0,105)) +
      geom_text(aes(label=paste0(round(value),"%")), position=position_dodge(0.7),
                vjust=-0.3, size=3) +
      labs(title="Efficacy & Safety Comparison (Clinical Trial Data)",
           x=NULL, y="Percent (%)") +
      theme_minimal(base_size=11) +
      theme(legend.position="bottom",
            plot.title=element_text(colour="#1a237e", face="bold"))
  })

  # ── TAB 6: Biomarkers & Toxicity ─────────────────────────────
  t_days_bm <- reactive(seq(0, input$bm_cycles*28, by=1))

  output$bm_ldh_plot <- renderPlot({
    t   <- t_days_bm()
    ldh <- sim_ldh(t, input$bm_ldh, input$bm_regimen, input$bm_dose/100)
    df  <- data.frame(day=t, LDH=ldh)

    ggplot(df, aes(x=day, y=LDH)) +
      geom_line(colour="#7B1FA2", linewidth=1.3) +
      geom_hline(yintercept=250, linetype="dashed", colour="#F44336") +
      geom_ribbon(aes(ymin=250, ymax=pmax(LDH,250)), fill="#E1BEE7", alpha=0.4) +
      annotate("text", x=max(t)*0.7, y=265, label="ULN", colour="#F44336", size=3.5) +
      labs(title="LDH Dynamics", x="Day", y="LDH (U/L)") +
      theme_minimal(base_size=11) +
      theme(plot.title=element_text(colour="#1a237e", face="bold"))
  })

  output$bm_tarc_plot <- renderPlot({
    t    <- t_days_bm()
    tarc <- sim_tarc(t, input$bm_tarc, input$bm_regimen, input$bm_dose/100)
    df   <- data.frame(day=t, TARC=tarc)

    ggplot(df, aes(x=day, y=TARC)) +
      geom_line(colour="#D32F2F", linewidth=1.3) +
      geom_hline(yintercept=450, linetype="dashed", colour="#F44336") +
      geom_ribbon(aes(ymin=450, ymax=pmax(TARC,450)), fill="#FFCDD2", alpha=0.4) +
      annotate("text", x=max(t)*0.7, y=480, label="Normal <450 pg/mL", colour="#F44336", size=3.2) +
      labs(title="TARC (CCL17) — RS Cell Biomarker", x="Day", y="TARC (pg/mL)") +
      theme_minimal(base_size=11) +
      theme(plot.title=element_text(colour="#1a237e", face="bold"))
  })

  output$bm_anc_plot <- renderPlot({
    anc_df <- sim_anc(input$bm_cycles, input$bm_regimen, input$bm_dose/100, input$bm_wt)
    # Scale by patient baseline
    anc_df$ANC <- anc_df$ANC * (input$bm_anc0 / 5.0)

    ggplot(anc_df, aes(x=day, y=ANC)) +
      geom_line(colour="#1976D2", linewidth=1.2) +
      geom_hline(yintercept=c(0.5,1.0), linetype="dashed",
                 colour=c("#F44336","#FF9800")) +
      annotate("text", x=5, y=0.55, label="Severe (<0.5)", colour="#F44336", size=3) +
      annotate("text", x=5, y=1.1,  label="Moderate (<1.0)",colour="#FF9800",size=3) +
      labs(title="ANC Dynamics (Friberg Model)", x="Day", y="ANC (×10⁹/L)") +
      scale_y_continuous(limits=c(0, NA)) +
      theme_minimal(base_size=11) +
      theme(plot.title=element_text(colour="#1a237e", face="bold"))
  })

  output$bm_plt_plot <- renderPlot({
    t   <- t_days_bm()
    # Platelet: similar Friberg structure but slower recovery
    plt0 <- input$bm_plt0
    slope_map <- c("ABVD"=0.4, "BV-AVD"=0.55, "Escalated BEACOPP"=0.85, "Pembrolizumab (R/R)"=0.1)
    slope <- slope_map[input$bm_regimen] * (input$bm_dose/100)

    PLT <- numeric(length(t))
    PLT[1] <- plt0
    for (i in seq_along(t)[-1]) {
      day <- t[i]
      cd  <- day %% 28
      drug_pulse <- if (cd <= 1 || (cd >= 14 && cd <= 15)) slope * exp(-0.5*pmin(cd,abs(cd-14))) else 0
      E <- drug_pulse
      feedback <- (plt0 / pmax(PLT[i-1], 1))^0.2
      growth <- 0.05 * log(plt0 / pmax(PLT[i-1],1)) * PLT[i-1]
      death  <- E * PLT[i-1]
      PLT[i] <- pmax(PLT[i-1] + growth - death, 5)
    }
    df <- data.frame(day=t, PLT=PLT)

    ggplot(df, aes(x=day, y=PLT)) +
      geom_line(colour="#388E3C", linewidth=1.2) +
      geom_hline(yintercept=c(50,100), linetype="dashed",
                 colour=c("#F44336","#FF9800")) +
      annotate("text", x=5, y=58,  label="Grade 3 (<50)",  colour="#F44336", size=3) +
      annotate("text", x=5, y=108, label="Grade 2 (<100)", colour="#FF9800", size=3) +
      labs(title="Platelet Dynamics", x="Day", y="PLT (×10⁹/L)") +
      theme_minimal(base_size=11) +
      theme(plot.title=element_text(colour="#1a237e", face="bold"))
  })

  output$bm_cd30_plot <- renderPlot({
    t      <- t_days_bm()
    cd30_0 <- input$bm_cd30
    # CD30 surface expression decreases with BV (direct target), stable with ABVD
    rate_map <- c("ABVD"=0.005, "BV-AVD"=0.06, "Escalated BEACOPP"=0.01, "Pembrolizumab (R/R)"=0.02)
    rate <- rate_map[input$bm_regimen]
    cd30 <- cd30_0 * exp(-rate * t)
    df   <- data.frame(day=t, CD30=cd30)

    ggplot(df, aes(x=day, y=CD30)) +
      geom_line(colour="#F57C00", linewidth=1.3) +
      geom_hline(yintercept=10, linetype="dashed", colour="#9E9E9E") +
      annotate("text", x=max(t)*0.7, y=12, label="BV sensitivity threshold ~10%", size=3, colour="#9E9E9E") +
      labs(title="CD30 Surface Expression on RS Cells", x="Day", y="CD30 Expression (%)") +
      scale_y_continuous(limits=c(0,100)) +
      theme_minimal(base_size=11) +
      theme(plot.title=element_text(colour="#1a237e", face="bold"))
  })

  output$bsym_resolution <- renderText({
    # B-symptom resolution: exponential, depends on regimen
    rate_map <- c("ABVD"=0.05, "BV-AVD"=0.07, "Escalated BEACOPP"=0.09, "Pembrolizumab (R/R)"=0.03)
    rate <- rate_map[input$bm_regimen] * (input$bm_dose/100)
    day  <- ceiling(-log(0.1) / rate)  # 90% resolution
    sprintf("~Day %d", day)
  })

  output$anc_nadir <- renderText({
    anc_df <- sim_anc(input$bm_cycles, input$bm_regimen, input$bm_dose/100, input$bm_wt)
    nadir  <- min(anc_df$ANC) * (input$bm_anc0 / 5.0)
    sprintf("%.2f ×10⁹/L", nadir)
  })

  output$bm_fn_risk <- renderText({
    fn_base <- c("ABVD"=18, "BV-AVD"=19, "Escalated BEACOPP"=45, "Pembrolizumab (R/R)"=5)
    risk    <- fn_base[input$bm_regimen] * (input$bm_dose/100)^1.2
    sprintf("%.1f%%", pmin(risk, 60))
  })
}

# ============================================================
# Launch
# ============================================================
shinyApp(ui = ui, server = server)
