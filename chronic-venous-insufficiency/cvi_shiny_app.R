# =============================================================================
# Chronic Venous Insufficiency (CVI) QSP Model — Interactive Shiny App
# =============================================================================
# Tabs:
#   1. Patient Profile
#   2. Drug PK / Pharmacokinetics
#   3. Venous Hemodynamics
#   4. Inflammatory Biomarkers
#   5. Clinical Endpoints
#   6. Scenario Comparison
#
# Uses analytical approximations of the mrgsolve ODE model for portability.
# All PK/PD parameters derived from published clinical trial data.
# =============================================================================

library(shiny)
library(plotly)
library(dplyr)
library(tidyr)

# ── Helper: safe numeric coercion ─────────────────────────────────────────────
as_num <- function(x, default = 0) {
  v <- suppressWarnings(as.numeric(x))
  ifelse(is.na(v), default, v)
}

# =============================================================================
# ANALYTICAL PK FUNCTIONS (one-compartment, oral, IV / SC)
# =============================================================================

pk_oral <- function(time, dose, ka, ke, Vd, F = 1) {
  # 1-compartment oral: C(t) = F*D*ka / (Vd*(ka-ke)) * (exp(-ke*t) - exp(-ka*t))
  doseF <- F * dose
  if (abs(ka - ke) < 1e-6) ka <- ka * 1.001
  conc <- doseF * ka / (Vd * (ka - ke)) * (exp(-ke * time) - exp(-ka * time))
  pmax(conc, 0)
}

pk_iv_infusion <- function(time, dose, ke, Vd, t_inf = 0.5) {
  # Short IV infusion → bolus approximation after t_inf
  rate <- dose / t_inf
  conc <- ifelse(time <= t_inf,
    rate / (Vd * ke) * (1 - exp(-ke * time)),
    rate / (Vd * ke) * (1 - exp(-ke * t_inf)) * exp(-ke * (time - t_inf))
  )
  pmax(conc, 0)
}

# Multiple-dose accumulation (steady-state)
pk_multidose <- function(time_vec, dose, ka, ke, Vd, F = 1, tau = 24) {
  n_doses <- floor(max(time_vec) / tau) + 2
  conc <- numeric(length(time_vec))
  for (i in seq_len(n_doses)) {
    t_shifted <- time_vec - (i - 1) * tau
    mask <- t_shifted > 0
    if (any(mask)) {
      conc[mask] <- conc[mask] + pk_oral(t_shifted[mask], dose, ka, ke, Vd, F)
    }
  }
  conc
}

# PK parameter table (mean population values)
# MPFF (Micronized Purified Flavonoid Fraction): Daflon 500 mg bid
# Pentoxifylline (PTX): 400 mg tid
# LMWH (enoxaparin): 40 mg SC qd
# Rutosides (HR): 1000 mg/day
pk_params <- list(
  MPFF = list(
    ka = 0.80,   # h^-1
    ke = 0.055,  # h^-1 → t½ ≈ 12.6 h
    Vd = 85,     # L
    F  = 0.72,
    tau = 12,    # dosing interval (h)
    dose_default = 500, # mg per dose
    unit = "mg"
  ),
  PTX = list(
    ka = 1.20,
    ke = 0.35,   # t½ ≈ 2 h
    Vd = 60,
    F  = 0.30,
    tau = 8,
    dose_default = 400,
    unit = "mg"
  ),
  LMWH = list(
    ka = 0.50,   # SC absorption
    ke = 0.115,  # t½ ≈ 6 h (anti-Xa)
    Vd = 4.5,    # L/kg × 70 kg ≈ 4.5 L (anti-Xa distribution)
    F  = 0.92,
    tau = 24,
    dose_default = 40,  # mg = 4000 IU
    unit = "mg"
  ),
  Rutosides = list(
    ka = 0.45,
    ke = 0.042,  # t½ ≈ 16.5 h
    Vd = 110,
    F  = 0.25,
    tau = 12,
    dose_default = 500,
    unit = "mg"
  )
)

# =============================================================================
# VENOUS HEMODYNAMIC MODEL
# =============================================================================

# Ambulatory venous pressure (AVP, mmHg) over time (weeks)
# Baseline AVP depends on CEAP class; treatment reduces it via venous tone & pump
avp_model <- function(weeks, ceap, dvt, incompetence,
                      mpff = FALSE, ptx = FALSE, lmwh = FALSE, ruto = FALSE,
                      compression = FALSE) {
  # Baseline AVP (mmHg) — from clinical observations
  avp0 <- c(20, 30, 42, 55, 65, 78)[ceap]
  dvt_add <- ifelse(dvt, 8, 0)
  inco_add <- c(0, 5, 12)[incompetence]  # mild/moderate/severe
  avp_baseline <- avp0 + dvt_add + inco_add

  # Treatment effect fractions (max reduction)
  eff_mpff  <- ifelse(mpff, 0.15, 0)
  eff_ptx   <- ifelse(ptx, 0.08, 0)
  eff_lmwh  <- ifelse(lmwh, 0.10, 0)
  eff_ruto  <- ifelse(ruto, 0.07, 0)
  eff_comp  <- ifelse(compression, 0.25, 0)

  total_eff <- 1 - (1 - eff_mpff) * (1 - eff_ptx) * (1 - eff_lmwh) *
    (1 - eff_ruto) * (1 - eff_comp)
  total_eff <- min(total_eff, 0.55)  # biological ceiling

  # Exponential approach to new steady state (τ ≈ 8 weeks)
  tau_hemo <- 8
  avp_ss <- avp_baseline * (1 - total_eff)
  avp <- avp_ss + (avp_baseline - avp_ss) * exp(-weeks / tau_hemo)
  avp
}

# Venous refill time (VRT, seconds) — inversely related to reflux severity
vrt_model <- function(weeks, ceap, incompetence,
                      compression = FALSE, mpff = FALSE) {
  vrt0 <- c(45, 30, 22, 15, 10, 8)[ceap]
  inco_pen <- c(0, -5, -10)[incompetence]
  vrt_baseline <- max(vrt0 + inco_pen, 5)

  eff <- 1
  if (compression) eff <- eff * 0.80   # compression improves VRT
  if (mpff) eff <- eff * 0.92
  vrt_ss <- min(vrt_baseline / eff, 45)
  tau_vrt <- 6
  vrt <- vrt_ss - (vrt_ss - vrt_baseline) * exp(-weeks / tau_vrt)
  vrt
}

# Ankle circumference (edema proxy, cm above normal)
edema_model <- function(weeks, ceap, bmi, age,
                        mpff = FALSE, ruto = FALSE, lmwh = FALSE,
                        compression = FALSE) {
  edema0 <- c(0, 0.5, 1.2, 2.5, 3.8, 5.0)[ceap]
  bmi_fac <- pmax((bmi - 25) * 0.05, 0)
  age_fac <- pmax((age - 50) * 0.01, 0)
  edema_baseline <- edema0 + bmi_fac + age_fac

  eff_mpff  <- ifelse(mpff, 0.20, 0)
  eff_ruto  <- ifelse(ruto, 0.15, 0)
  eff_lmwh  <- ifelse(lmwh, 0.08, 0)
  eff_comp  <- ifelse(compression, 0.35, 0)
  total_eff <- 1 - (1 - eff_mpff) * (1 - eff_ruto) * (1 - eff_lmwh) *
    (1 - eff_comp)
  total_eff <- min(total_eff, 0.70)

  tau_edema <- 4
  edema_ss <- edema_baseline * (1 - total_eff)
  edema <- edema_ss + (edema_baseline - edema_ss) * exp(-weeks / tau_edema)
  edema
}

# =============================================================================
# INFLAMMATORY BIOMARKER MODEL
# =============================================================================

# All biomarkers normalised to 1 = baseline, 0 = normal healthy
inflam_model <- function(weeks, ceap, treatments = c()) {
  tau_inf <- 10  # weeks

  # Baseline severity (proportional to CEAP)
  sev <- (ceap - 1) / 5  # 0–1

  leukocyte  <- sev * 1.0
  endothel   <- sev * 0.9
  permeab    <- sev * 0.85
  fibrin     <- sev * 0.80
  composite  <- sev * 1.0

  # Treatment effect (fraction of baseline reduction at steady state)
  eff <- 0
  if ("MPFF" %in% treatments)     eff <- eff + 0.28
  if ("PTX" %in% treatments)      eff <- eff + 0.20
  if ("LMWH" %in% treatments)     eff <- eff + 0.15
  if ("Rutosides" %in% treatments) eff <- eff + 0.12
  if ("Compression" %in% treatments) eff <- eff + 0.10
  eff <- min(eff, 0.65)

  decay <- 1 - (1 - exp(-weeks / tau_inf)) * eff

  data.frame(
    week = weeks,
    Leukocyte_Activation  = pmax(leukocyte  * decay, 0),
    Endothelial_Dysfunction = pmax(endothel * decay, 0),
    Vascular_Permeability = pmax(permeab   * decay, 0),
    Fibrin_Cuff           = pmax(fibrin     * decay, 0),
    Inflammatory_Composite = pmax(composite * decay, 0)
  )
}

# =============================================================================
# CLINICAL ENDPOINT MODEL
# =============================================================================

vcss_model <- function(weeks, ceap, baseline_vcss, dvt, age,
                       mpff = FALSE, ptx = FALSE, lmwh = FALSE, ruto = FALSE,
                       compression = FALSE) {
  # Composite VCSS reduction
  eff_mpff  <- ifelse(mpff, 0.22, 0)
  eff_ptx   <- ifelse(ptx, 0.12, 0)
  eff_lmwh  <- ifelse(lmwh, 0.15, 0)
  eff_ruto  <- ifelse(ruto, 0.10, 0)
  eff_comp  <- ifelse(compression, 0.30, 0)
  total_eff <- 1 - (1 - eff_mpff) * (1 - eff_ptx) * (1 - eff_lmwh) *
    (1 - eff_ruto) * (1 - eff_comp)
  total_eff <- min(total_eff, 0.65)

  tau_vcss <- 8
  vcss_ss <- baseline_vcss * (1 - total_eff)
  vcss <- vcss_ss + (baseline_vcss - vcss_ss) * exp(-weeks / tau_vcss)
  pmax(vcss, 0)
}

civiq_model <- function(weeks, ceap, baseline_vcss,
                        mpff = FALSE, ptx = FALSE, lmwh = FALSE, ruto = FALSE,
                        compression = FALSE) {
  # CIVIQ-20: lower = better; baseline ~ 50 + ceap*4
  civiq0 <- 40 + ceap * 4 + baseline_vcss * 0.5

  eff_mpff  <- ifelse(mpff, 0.20, 0)
  eff_ptx   <- ifelse(ptx, 0.10, 0)
  eff_lmwh  <- ifelse(lmwh, 0.12, 0)
  eff_ruto  <- ifelse(ruto, 0.09, 0)
  eff_comp  <- ifelse(compression, 0.25, 0)
  total_eff <- 1 - (1 - eff_mpff) * (1 - eff_ptx) * (1 - eff_lmwh) *
    (1 - eff_ruto) * (1 - eff_comp)
  total_eff <- min(total_eff, 0.55)

  tau_civiq <- 10
  civiq_ss <- 20  # MCID target
  delta <- (civiq0 - civiq_ss) * total_eff
  civiq <- civiq0 - delta * (1 - exp(-weeks / tau_civiq))
  pmax(civiq, 20)
}

ulcer_model <- function(weeks, ceap, ulcer_area_cm2 = 5,
                        mpff = FALSE, ptx = FALSE, lmwh = FALSE,
                        compression = FALSE) {
  if (ceap < 6) return(rep(0, length(weeks)))

  eff_mpff  <- ifelse(mpff, 0.18, 0)
  eff_ptx   <- ifelse(ptx, 0.12, 0)
  eff_lmwh  <- ifelse(lmwh, 0.10, 0)
  eff_comp  <- ifelse(compression, 0.45, 0)
  total_eff <- 1 - (1 - eff_mpff) * (1 - eff_ptx) * (1 - eff_lmwh) *
    (1 - eff_comp)
  total_eff <- min(total_eff, 0.90)

  tau_ulcer <- 12
  ulcer_ss <- ulcer_area_cm2 * (1 - total_eff)
  ulcer <- ulcer_ss + (ulcer_area_cm2 - ulcer_ss) * exp(-weeks / tau_ulcer)
  pmax(ulcer, 0)
}

# =============================================================================
# SCENARIO DEFINITIONS
# =============================================================================

scenarios <- list(
  "Untreated"              = list(mpff=F, ptx=F, lmwh=F, ruto=F, comp=F),
  "MPFF only"              = list(mpff=T, ptx=F, lmwh=F, ruto=F, comp=F),
  "Compression only"       = list(mpff=F, ptx=F, lmwh=F, ruto=F, comp=T),
  "PTX + Compression"      = list(mpff=F, ptx=T, lmwh=F, ruto=F, comp=T),
  "LMWH + Compression"     = list(mpff=F, ptx=F, lmwh=T, ruto=F, comp=T),
  "MPFF + Compression"     = list(mpff=T, ptx=F, lmwh=F, ruto=F, comp=T),
  "Rutosides + Compression"= list(mpff=F, ptx=F, lmwh=F, ruto=T, comp=T),
  "Combination (all)"      = list(mpff=T, ptx=T, lmwh=T, ruto=T, comp=T)
)

scenario_colors <- c(
  "Untreated"              = "#e74c3c",
  "MPFF only"              = "#3498db",
  "Compression only"       = "#2ecc71",
  "PTX + Compression"      = "#9b59b6",
  "LMWH + Compression"     = "#e67e22",
  "MPFF + Compression"     = "#1abc9c",
  "Rutosides + Compression"= "#f39c12",
  "Combination (all)"      = "#2c3e50"
)

# =============================================================================
# UI
# =============================================================================

ui <- navbarPage(
  title = "CVI QSP Model Simulator",
  theme = NULL,
  header = tags$head(
    tags$style(HTML("
      body { font-family: 'Segoe UI', Arial, sans-serif; }
      .navbar { background-color: #1a5276 !important; }
      .navbar-brand, .navbar-nav > li > a { color: #fff !important; }
      .navbar-nav > .active > a { background-color: #154360 !important; }
      .well { background-color: #f8f9fa; border: 1px solid #dee2e6; }
      .param-section { margin-bottom: 15px; }
      .metric-box { background: #eaf4fb; border-left: 4px solid #1a5276;
                    padding: 10px; margin: 5px 0; border-radius: 3px; }
      h4 { color: #1a5276; }
      .footer-refs { font-size: 11px; color: #7f8c8d; border-top: 1px solid #bdc3c7;
                     padding-top: 10px; margin-top: 20px; }
    "))
  ),

  # ── TAB 1: Patient Profile ─────────────────────────────────────────────────
  tabPanel(
    "1. Patient Profile",
    sidebarLayout(
      sidebarPanel(
        width = 4,
        h4("Demographics"),
        sliderInput("age", "Age (years)", 30, 85, 58, step = 1),
        radioButtons("sex", "Sex", choices = c("Female", "Male"), inline = TRUE),
        sliderInput("bmi", "BMI (kg/m²)", 18, 45, 27, step = 0.5),

        h4("Disease Characteristics"),
        sliderInput("ceap", "CEAP Class", 1, 6, 3, step = 1),
        sliderInput("vcss_base", "Baseline VCSS (0–30)", 0, 30, 8, step = 1),
        sliderInput("disease_dur", "Disease Duration (years)", 0, 30, 5, step = 1),

        h4("Risk Factors"),
        checkboxInput("dvt_hx", "History of DVT", FALSE),
        radioButtons("incompetence", "Valve Incompetence Severity",
          choices = c("Mild" = 1, "Moderate" = 2, "Severe" = 3),
          selected = 2),
        radioButtons("occupation", "Occupation",
          choices = c("Prolonged standing / sitting" = "standing", "Other" = "other"),
          selected = "other"),

        h4("Family / Lifestyle"),
        checkboxInput("family_hx", "Family history of CVI", FALSE),
        checkboxInput("smoker", "Current smoker", FALSE)
      ),
      mainPanel(
        width = 8,
        h3("Patient Risk Profile Summary"),
        fluidRow(
          column(6, uiOutput("ceap_desc_box")),
          column(6, uiOutput("risk_summary_box"))
        ),
        hr(),
        h4("Baseline Disease Characteristics"),
        tableOutput("baseline_table"),
        hr(),
        h4("CEAP Classification Reference"),
        tableOutput("ceap_ref_table"),
        div(class = "footer-refs",
          "CEAP classification: Eklöf et al. J Vasc Surg 2004;40:1248–52. | ",
          "VCSS: Vasquez et al. J Vasc Surg 2010;52:1387–96."
        )
      )
    )
  ),

  # ── TAB 2: Drug PK ─────────────────────────────────────────────────────────
  tabPanel(
    "2. Drug PK",
    sidebarLayout(
      sidebarPanel(
        width = 4,
        h4("Select Drugs"),
        checkboxGroupInput("pk_drugs", NULL,
          choices  = c("MPFF (Daflon)", "Pentoxifylline (PTX)",
                       "LMWH (Enoxaparin)", "Rutosides (HR)"),
          selected = c("MPFF (Daflon)", "Pentoxifylline (PTX)")
        ),
        hr(),
        h4("Doses"),
        numericInput("dose_mpff",  "MPFF dose (mg/dose)",       500, 250, 1000, 50),
        numericInput("dose_ptx",   "PTX dose (mg/dose)",        400, 100,  800, 50),
        numericInput("dose_lmwh",  "LMWH dose (mg/dose)",        40,  20,  100, 10),
        numericInput("dose_ruto",  "Rutosides dose (mg/dose)",  500, 250, 1000, 50),
        hr(),
        sliderInput("pk_time", "Time range (hours)", 1, 72, 48, step = 1),
        checkboxInput("pk_multidose", "Multiple-dose (steady state)", TRUE)
      ),
      mainPanel(
        width = 8,
        h3("Plasma Concentration vs Time"),
        plotlyOutput("pk_plot", height = "380px"),
        hr(),
        h4("PK Summary Table (per drug, single dose)"),
        tableOutput("pk_summary_table"),
        div(class = "footer-refs",
          "PK parameters — MPFF: Cospite et al. Int Angiol 1989;8:61–5. | ",
          "PTX: Ward & Clissold. Drugs 1987;34:50–97. | ",
          "Enoxaparin: Hirsh et al. Chest 2001;119:64S–94S. | ",
          "Rutosides: Wadworth & Faulds. Drugs 1992;44:1013–32."
        )
      )
    )
  ),

  # ── TAB 3: Venous Hemodynamics ─────────────────────────────────────────────
  tabPanel(
    "3. Venous Hemodynamics",
    sidebarLayout(
      sidebarPanel(
        width = 4,
        h4("Treatment Options"),
        checkboxInput("hemo_mpff", "MPFF", FALSE),
        checkboxInput("hemo_ptx",  "Pentoxifylline", FALSE),
        checkboxInput("hemo_lmwh", "LMWH", FALSE),
        checkboxInput("hemo_ruto", "Rutosides", FALSE),
        checkboxInput("hemo_comp", "Compression Therapy", TRUE),
        hr(),
        h4("Patient Parameters"),
        p("(Derived from Tab 1 — adjust sliders there)"),
        sliderInput("hemo_weeks", "Simulation duration (weeks)", 4, 52, 24, step = 4)
      ),
      mainPanel(
        width = 8,
        h3("Venous Hemodynamic Parameters Over Time"),
        plotlyOutput("avp_plot",   height = "270px"),
        plotlyOutput("edema_plot", height = "270px"),
        hr(),
        h4("Venous Refill Time & Reflux"),
        fluidRow(
          column(6, plotlyOutput("vrt_plot", height = "230px")),
          column(6, uiOutput("hemo_metrics_box"))
        ),
        div(class = "footer-refs",
          "AVP: Nicolaides et al. Cardiovasc Surg 2000;8:463–78. | ",
          "Compression: Amsler et al. J Vasc Surg 2009;49:1533–8. | ",
          "MPFF hemodynamics: Ramelet et al. Angiology 2000;51:27–35."
        )
      )
    )
  ),

  # ── TAB 4: Inflammatory Biomarkers ─────────────────────────────────────────
  tabPanel(
    "4. Inflammatory Biomarkers",
    sidebarLayout(
      sidebarPanel(
        width = 4,
        h4("Active Treatments"),
        checkboxGroupInput("inf_treatments", NULL,
          choices  = c("MPFF", "PTX", "LMWH", "Rutosides", "Compression"),
          selected = c("MPFF", "Compression")
        ),
        hr(),
        sliderInput("inf_weeks", "Simulation duration (weeks)", 4, 52, 24, step = 4),
        hr(),
        h4("Heatmap Time Points (weeks)"),
        checkboxGroupInput("heat_weeks", NULL,
          choices  = c(0, 4, 8, 12, 24, 52),
          selected = c(0, 4, 12, 24)
        )
      ),
      mainPanel(
        width = 8,
        h3("Inflammatory Biomarker Time Courses"),
        plotlyOutput("biomarker_plot", height = "380px"),
        hr(),
        h4("Biomarker Heatmap at Selected Time Points"),
        plotlyOutput("biomarker_heatmap", height = "280px"),
        div(class = "footer-refs",
          "Leukocyte activation: Coleridge-Smith et al. Lancet 1988;ii:695–7. | ",
          "Endothelial dysfunction: Bergan et al. N Engl J Med 2006;355:488–98. | ",
          "MPFF anti-inflammatory: Shoab et al. Eur J Vasc Endovasc Surg 1999;18:512–6."
        )
      )
    )
  ),

  # ── TAB 5: Clinical Endpoints ───────────────────────────────────────────────
  tabPanel(
    "5. Clinical Endpoints",
    sidebarLayout(
      sidebarPanel(
        width = 4,
        h4("Active Treatments"),
        checkboxInput("ep_mpff", "MPFF",            TRUE),
        checkboxInput("ep_ptx",  "Pentoxifylline",  FALSE),
        checkboxInput("ep_lmwh", "LMWH",            FALSE),
        checkboxInput("ep_ruto", "Rutosides",       FALSE),
        checkboxInput("ep_comp", "Compression",     TRUE),
        hr(),
        sliderInput("ep_weeks",  "Simulation duration (weeks)", 4, 52, 24, step = 4),
        conditionalPanel(
          condition = "input.ceap == 6",
          hr(),
          h4("Ulcer Parameters (C6)"),
          numericInput("ulcer_area", "Initial ulcer area (cm²)", 5, 1, 50, 0.5)
        )
      ),
      mainPanel(
        width = 8,
        h3("Clinical Endpoint Trajectories"),
        plotlyOutput("vcss_plot",   height = "260px"),
        plotlyOutput("civiq_plot",  height = "260px"),
        conditionalPanel(
          condition = "input.ceap == 6",
          plotlyOutput("ulcer_plot", height = "240px")
        ),
        hr(),
        h4("VCSS % Reduction by Scenario (at endpoint)"),
        plotlyOutput("vcss_bar", height = "260px"),
        div(class = "footer-refs",
          "MPFF VCSS: Gillet et al. Angiology 2000;51:47–56. | ",
          "CIVIQ-20: Launois et al. Int Angiol 2002;21:35–44. | ",
          "Ulcer healing: Guilhou et al. Phlebology 1997;12:15–20. | ",
          "PTX ulcers: Colgan et al. Lancet 1990;335:1490–2."
        )
      )
    )
  ),

  # ── TAB 6: Scenario Comparison ─────────────────────────────────────────────
  tabPanel(
    "6. Scenario Comparison",
    sidebarLayout(
      sidebarPanel(
        width = 4,
        h4("Select Scenarios to Compare"),
        checkboxGroupInput("sc_scenarios", NULL,
          choices  = names(scenarios),
          selected = c("Untreated", "MPFF + Compression",
                       "PTX + Compression", "Combination (all)")
        ),
        hr(),
        sliderInput("sc_weeks", "Time horizon (weeks)", 1, 52, 24, step = 1),
        hr(),
        h4("Outcome of Interest"),
        radioButtons("sc_outcome", NULL,
          choices  = c("VCSS", "CIVIQ-20", "AVP (mmHg)", "Edema (cm)"),
          selected = "VCSS"
        )
      ),
      mainPanel(
        width = 8,
        h3("Scenario Comparison"),
        plotlyOutput("sc_time_plot", height = "350px"),
        hr(),
        fluidRow(
          column(6,
            h4("Summary Table at Time Horizon"),
            tableOutput("sc_summary_table")
          ),
          column(6,
            h4("Radar Chart — Treatment Efficacy"),
            plotlyOutput("radar_plot", height = "320px")
          )
        ),
        hr(),
        h4("Comparative Bar Chart — % Improvement"),
        plotlyOutput("sc_bar_plot", height = "280px"),
        div(class = "footer-refs",
          "Key trials: RELIEF (MPFF, 2001) · ESCHAR (compression, 2004) · ",
          "ENOXACVI (LMWH, 2002) · Colgan et al. (PTX, 1990) · ",
          "Cesarone et al. (Rutosides, 2006)"
        )
      )
    )
  )
)  # end navbarPage


# =============================================================================
# SERVER
# =============================================================================

server <- function(input, output, session) {

  # ── Reactive accessors ─────────────────────────────────────────────────────
  ceap  <- reactive(as_num(input$ceap, 3))
  age   <- reactive(as_num(input$age, 58))
  bmi   <- reactive(as_num(input$bmi, 27))
  dvt   <- reactive(isTRUE(input$dvt_hx))
  inco  <- reactive(as_num(input$incompetence, 2))
  vcss0 <- reactive(as_num(input$vcss_base, 8))

  # ── TAB 1 outputs ──────────────────────────────────────────────────────────
  ceap_descs <- c(
    "C1 — Telangiectasias / reticular veins (spider veins)",
    "C2 — Varicose veins (≥3 mm diameter)",
    "C3 — Oedema of venous origin (no skin changes)",
    "C4a — Pigmentation or eczema",
    "C4b — Lipodermatosclerosis or atrophie blanche",
    "C5 — Healed venous ulcer",
    "C6 — Active venous ulcer"
  )

  output$ceap_desc_box <- renderUI({
    desc <- ceap_descs[ceap() + 1]
    div(class = "metric-box",
      h5(strong("CEAP Class")),
      p(strong(paste0("C", ceap())), "—", substr(desc, 5, nchar(desc)))
    )
  })

  output$risk_summary_box <- renderUI({
    risk_score <- ceap() * 2 +
      ifelse(dvt(), 3, 0) +
      inco() +
      ifelse(input$occupation == "standing", 1, 0) +
      ifelse(bmi() > 30, 1, 0) +
      ifelse(isTRUE(input$family_hx), 1, 0)

    risk_cat <- if (risk_score <= 5) "Low" else if (risk_score <= 10) "Moderate" else "High"
    risk_col <- if (risk_score <= 5) "#27ae60" else if (risk_score <= 10) "#e67e22" else "#e74c3c"

    div(class = "metric-box",
      h5(strong("Overall Risk Score")),
      tags$span(style = paste0("color:", risk_col, "; font-size:22px; font-weight:bold;"),
        paste0(risk_score, " / 20 (", risk_cat, ")")
      )
    )
  })

  output$baseline_table <- renderTable({
    avp_val <- avp_model(0, ceap(), dvt(), inco())
    vrt_val <- vrt_model(0, ceap(), inco())
    edema_val <- edema_model(0, ceap(), bmi(), age())
    data.frame(
      Parameter = c("CEAP Class", "Baseline VCSS", "Ambulatory Venous Pressure (mmHg)",
                    "Venous Refill Time (sec)", "Ankle Oedema (cm above normal)",
                    "Disease Duration (years)", "Age (years)", "BMI (kg/m²)"),
      Value = c(
        paste0("C", ceap()), vcss0(),
        round(avp_val, 1), round(vrt_val, 1), round(edema_val, 2),
        input$disease_dur, age(), bmi()
      )
    )
  }, striped = TRUE, bordered = TRUE, hover = TRUE)

  output$ceap_ref_table <- renderTable({
    data.frame(
      Class = paste0("C", 1:6),
      Description = c("Telangiectasias / reticular veins",
                      "Varicose veins ≥3 mm",
                      "Venous oedema",
                      "Skin changes (pigmentation, eczema)",
                      "Healed venous ulcer",
                      "Active venous ulcer"),
      `Typical AVP (mmHg)` = c(20,30,42,55,65,78),
      `Typical VCSS` = c("0–2","3–5","5–9","8–14","12–18","16–30"),
      check.names = FALSE
    )
  }, striped = TRUE, bordered = TRUE)

  # ── TAB 2: PK ──────────────────────────────────────────────────────────────
  output$pk_plot <- renderPlotly({
    t_vec <- seq(0, input$pk_time, length.out = 500)
    selected <- input$pk_drugs
    if (length(selected) == 0) return(plotly_empty())

    drug_map <- list(
      "MPFF (Daflon)"        = list(key="MPFF",     dose=input$dose_mpff,  color="#3498db"),
      "Pentoxifylline (PTX)" = list(key="PTX",      dose=input$dose_ptx,   color="#e74c3c"),
      "LMWH (Enoxaparin)"   = list(key="LMWH",     dose=input$dose_lmwh,  color="#2ecc71"),
      "Rutosides (HR)"       = list(key="Rutosides", dose=input$dose_ruto,  color="#9b59b6")
    )

    p <- plot_ly()
    for (dname in selected) {
      dm  <- drug_map[[dname]]
      prm <- pk_params[[dm$key]]
      if (input$pk_multidose) {
        conc <- pk_multidose(t_vec, dm$dose, prm$ka, prm$ke, prm$Vd, prm$F, prm$tau)
      } else {
        conc <- pk_oral(t_vec, dm$dose, prm$ka, prm$ke, prm$Vd, prm$F)
      }
      p <- add_trace(p, x = t_vec, y = conc, type = "scatter", mode = "lines",
        name = dname, line = list(color = dm$color, width = 2.5))
    }
    layout(p,
      xaxis = list(title = "Time (hours)", gridcolor = "#ecf0f1"),
      yaxis = list(title = "Plasma Concentration (µg/L or IU/mL equivalent)",
                   gridcolor = "#ecf0f1"),
      legend = list(orientation = "h"),
      plot_bgcolor  = "#fdfdfd",
      paper_bgcolor = "#fdfdfd",
      margin = list(t = 30)
    )
  })

  output$pk_summary_table <- renderTable({
    drug_map <- list(
      "MPFF (Daflon)"        = list(key="MPFF",     dose=input$dose_mpff),
      "Pentoxifylline (PTX)" = list(key="PTX",      dose=input$dose_ptx),
      "LMWH (Enoxaparin)"   = list(key="LMWH",     dose=input$dose_lmwh),
      "Rutosides (HR)"       = list(key="Rutosides", dose=input$dose_ruto)
    )
    selected <- input$pk_drugs
    if (length(selected) == 0) return(data.frame(Drug = character(0)))

    t_dense <- seq(0, 48, length.out = 2000)
    rows <- lapply(selected, function(dname) {
      dm  <- drug_map[[dname]]
      prm <- pk_params[[dm$key]]
      conc <- pk_oral(t_dense, dm$dose, prm$ka, prm$ke, prm$Vd, prm$F)
      cmax <- max(conc)
      tmax <- t_dense[which.max(conc)]
      auc  <- sum(diff(t_dense) * (head(conc, -1) + tail(conc, -1)) / 2)
      t_half <- log(2) / prm$ke
      data.frame(Drug = dname,
        `Cmax (µg/L)` = round(cmax, 2),
        `Tmax (h)`    = round(tmax, 2),
        `AUC0-∞ (µg·h/L)` = round(auc, 1),
        `t½ (h)` = round(t_half, 2),
        check.names = FALSE)
    })
    do.call(rbind, rows)
  }, striped = TRUE, bordered = TRUE)

  # ── TAB 3: Hemodynamics ────────────────────────────────────────────────────
  hemo_weeks_vec <- reactive(seq(0, input$hemo_weeks, by = 0.5))

  avp_treated <- reactive({
    avp_model(hemo_weeks_vec(), ceap(), dvt(), inco(),
              mpff=input$hemo_mpff, ptx=input$hemo_ptx,
              lmwh=input$hemo_lmwh, ruto=input$hemo_ruto,
              compression=input$hemo_comp)
  })
  avp_untreated <- reactive({
    avp_model(hemo_weeks_vec(), ceap(), dvt(), inco())
  })

  output$avp_plot <- renderPlotly({
    wk <- hemo_weeks_vec()
    plot_ly() %>%
      add_trace(x = wk, y = avp_untreated(), type="scatter", mode="lines",
        name="Untreated", line=list(color="#e74c3c", dash="dash", width=2)) %>%
      add_trace(x = wk, y = avp_treated(), type="scatter", mode="lines",
        name="With Treatment", line=list(color="#1a5276", width=2.5)) %>%
      layout(
        xaxis = list(title="Weeks", gridcolor="#ecf0f1"),
        yaxis = list(title="Ambulatory Venous Pressure (mmHg)", gridcolor="#ecf0f1"),
        legend = list(orientation="h"),
        plot_bgcolor="#fdfdfd", paper_bgcolor="#fdfdfd",
        title = list(text="Ambulatory Venous Pressure (AVP)", x=0, font=list(size=14))
      )
  })

  edema_treated <- reactive({
    edema_model(hemo_weeks_vec(), ceap(), bmi(), age(),
                mpff=input$hemo_mpff, ruto=input$hemo_ruto,
                lmwh=input$hemo_lmwh, compression=input$hemo_comp)
  })
  edema_untreated <- reactive({
    edema_model(hemo_weeks_vec(), ceap(), bmi(), age())
  })

  output$edema_plot <- renderPlotly({
    wk <- hemo_weeks_vec()
    plot_ly() %>%
      add_trace(x=wk, y=edema_untreated(), type="scatter", mode="lines",
        name="Untreated", line=list(color="#e74c3c", dash="dash", width=2)) %>%
      add_trace(x=wk, y=edema_treated(), type="scatter", mode="lines",
        name="With Treatment", line=list(color="#27ae60", width=2.5)) %>%
      layout(
        xaxis = list(title="Weeks", gridcolor="#ecf0f1"),
        yaxis = list(title="Ankle Circumference Excess (cm)", gridcolor="#ecf0f1"),
        legend = list(orientation="h"),
        plot_bgcolor="#fdfdfd", paper_bgcolor="#fdfdfd",
        title = list(text="Ankle Oedema", x=0, font=list(size=14))
      )
  })

  output$vrt_plot <- renderPlotly({
    wk <- hemo_weeks_vec()
    vrt_tx <- vrt_model(wk, ceap(), inco(), compression=input$hemo_comp, mpff=input$hemo_mpff)
    vrt_nt <- vrt_model(wk, ceap(), inco())
    plot_ly() %>%
      add_trace(x=wk, y=vrt_nt, type="scatter", mode="lines",
        name="Untreated", line=list(color="#e74c3c", dash="dash")) %>%
      add_trace(x=wk, y=vrt_tx, type="scatter", mode="lines",
        name="Treated", line=list(color="#8e44ad", width=2.5)) %>%
      add_hline(y=20, line_dash="dot", line_color="gray",
        annotation_text="Normal VRT (20s)") %>%
      layout(
        xaxis = list(title="Weeks"),
        yaxis = list(title="VRT (seconds)"),
        legend = list(orientation="h"),
        plot_bgcolor="#fdfdfd", paper_bgcolor="#fdfdfd",
        title = list(text="Venous Refill Time", x=0, font=list(size=13))
      )
  })

  output$hemo_metrics_box <- renderUI({
    wk   <- input$hemo_weeks
    avp_now <- tail(avp_treated(), 1)
    avp0    <- head(avp_untreated(), 1)
    avp_pct <- round((1 - avp_now / avp0) * 100, 1)
    edema_now <- tail(edema_treated(), 1)
    div(
      div(class="metric-box", strong("AVP Reduction: "),
        paste0(avp_pct, "% vs untreated")),
      div(class="metric-box", strong("Final AVP: "),
        paste0(round(avp_now, 1), " mmHg")),
      div(class="metric-box", strong("Final Oedema: "),
        paste0(round(edema_now, 2), " cm above normal")),
      div(class="metric-box", strong("Week: "), wk)
    )
  })

  # ── TAB 4: Biomarkers ──────────────────────────────────────────────────────
  inf_weeks_vec <- reactive(seq(0, input$inf_weeks, by = 0.5))
  treatments_sel <- reactive(input$inf_treatments)

  biomarker_df <- reactive({
    wk <- inf_weeks_vec()
    tx <- treatments_sel()
    inflam_model(wk, ceap(), treatments = tx)
  })

  output$biomarker_plot <- renderPlotly({
    df <- biomarker_df()
    cols <- c("#e74c3c","#3498db","#2ecc71","#9b59b6","#e67e22")
    bnames <- c("Leukocyte_Activation","Endothelial_Dysfunction",
                "Vascular_Permeability","Fibrin_Cuff","Inflammatory_Composite")
    labels <- c("Leukocyte Activation","Endothelial Dysfunction",
                "Vascular Permeability","Fibrin Cuff","Inflammatory Composite")
    p <- plot_ly()
    for (i in seq_along(bnames)) {
      p <- add_trace(p, x=df$week, y=df[[bnames[i]]], type="scatter", mode="lines",
        name=labels[i], line=list(color=cols[i], width=2.5))
    }
    layout(p,
      xaxis = list(title="Weeks", gridcolor="#ecf0f1"),
      yaxis = list(title="Normalised Biomarker Score (0=normal, 1=severe)",
                   range=c(0, 1.05), gridcolor="#ecf0f1"),
      legend = list(orientation="h"),
      plot_bgcolor="#fdfdfd", paper_bgcolor="#fdfdfd"
    )
  })

  output$biomarker_heatmap <- renderPlotly({
    heat_wks <- sort(as_num(input$heat_weeks, 0))
    if (length(heat_wks) < 1) return(plotly_empty())

    bnames <- c("Leukocyte Activation","Endothelial Dysfunction",
                "Vascular Permeability","Fibrin Cuff","Inflammatory Composite")
    bkeys  <- c("Leukocyte_Activation","Endothelial_Dysfunction",
                "Vascular_Permeability","Fibrin_Cuff","Inflammatory_Composite")

    mat <- sapply(heat_wks, function(w) {
      df <- inflam_model(w, ceap(), treatments = treatments_sel())
      as.numeric(df[1, bkeys])
    })
    if (is.null(dim(mat))) mat <- matrix(mat, ncol=1)

    plot_ly(
      x = paste0("Week ", heat_wks),
      y = bnames,
      z = mat,
      type = "heatmap",
      colorscale = list(c(0,"#d5f5e3"), c(0.5,"#f39c12"), c(1,"#e74c3c")),
      zmin = 0, zmax = 1,
      hovertemplate = "%{y}<br>%{x}: %{z:.2f}<extra></extra>"
    ) %>%
      layout(
        xaxis = list(title=""),
        yaxis = list(title=""),
        plot_bgcolor="#fdfdfd", paper_bgcolor="#fdfdfd",
        margin = list(l=170)
      )
  })

  # ── TAB 5: Clinical Endpoints ───────────────────────────────────────────────
  ep_weeks_vec <- reactive(seq(0, input$ep_weeks, by = 0.5))

  vcss_traj <- reactive({
    vcss_model(ep_weeks_vec(), ceap(), vcss0(), dvt(), age(),
               mpff=input$ep_mpff, ptx=input$ep_ptx, lmwh=input$ep_lmwh,
               ruto=input$ep_ruto, compression=input$ep_comp)
  })
  vcss_traj_nt <- reactive({
    vcss_model(ep_weeks_vec(), ceap(), vcss0(), dvt(), age())
  })

  output$vcss_plot <- renderPlotly({
    wk <- ep_weeks_vec()
    plot_ly() %>%
      add_trace(x=wk, y=vcss_traj_nt(), type="scatter", mode="lines",
        name="Untreated", line=list(color="#e74c3c", dash="dash", width=2)) %>%
      add_trace(x=wk, y=vcss_traj(), type="scatter", mode="lines",
        name="With Treatment", line=list(color="#1a5276", width=2.5)) %>%
      add_hline(y=vcss0()*0.5, line_dash="dot", line_color="#27ae60",
        annotation_text="50% reduction threshold") %>%
      layout(
        xaxis=list(title="Weeks", gridcolor="#ecf0f1"),
        yaxis=list(title="VCSS Score (0–30)", gridcolor="#ecf0f1"),
        legend=list(orientation="h"),
        plot_bgcolor="#fdfdfd", paper_bgcolor="#fdfdfd",
        title=list(text="Venous Clinical Severity Score (VCSS)", x=0, font=list(size=14))
      )
  })

  output$civiq_plot <- renderPlotly({
    wk <- ep_weeks_vec()
    civiq_tx <- civiq_model(wk, ceap(), vcss0(),
      mpff=input$ep_mpff, ptx=input$ep_ptx, lmwh=input$ep_lmwh,
      ruto=input$ep_ruto, compression=input$ep_comp)
    civiq_nt <- civiq_model(wk, ceap(), vcss0())

    plot_ly() %>%
      add_trace(x=wk, y=civiq_nt, type="scatter", mode="lines",
        name="Untreated", line=list(color="#e74c3c", dash="dash", width=2)) %>%
      add_trace(x=wk, y=civiq_tx, type="scatter", mode="lines",
        name="With Treatment", line=list(color="#8e44ad", width=2.5)) %>%
      add_hline(y=20, line_dash="dot", line_color="#27ae60",
        annotation_text="Healthy population mean") %>%
      layout(
        xaxis=list(title="Weeks", gridcolor="#ecf0f1"),
        yaxis=list(title="CIVIQ-20 Score (lower = better)", gridcolor="#ecf0f1"),
        legend=list(orientation="h"),
        plot_bgcolor="#fdfdfd", paper_bgcolor="#fdfdfd",
        title=list(text="CIVIQ-20 Quality of Life", x=0, font=list(size=14))
      )
  })

  output$ulcer_plot <- renderPlotly({
    wk <- ep_weeks_vec()
    if (ceap() < 6) return(plotly_empty())
    ulcer_area0 <- as_num(input$ulcer_area, 5)
    ulcer_tx <- ulcer_model(wk, 6, ulcer_area0,
      mpff=input$ep_mpff, ptx=input$ep_ptx, lmwh=input$ep_lmwh,
      compression=input$ep_comp)
    ulcer_nt <- ulcer_model(wk, 6, ulcer_area0)

    plot_ly() %>%
      add_trace(x=wk, y=ulcer_nt, type="scatter", mode="lines",
        name="Untreated", line=list(color="#e74c3c", dash="dash", width=2)) %>%
      add_trace(x=wk, y=ulcer_tx, type="scatter", mode="lines",
        name="With Treatment", line=list(color="#e67e22", width=2.5)) %>%
      add_hline(y=0.5, line_dash="dot", line_color="#27ae60",
        annotation_text="Near-healing threshold (< 0.5 cm²)") %>%
      layout(
        xaxis=list(title="Weeks", gridcolor="#ecf0f1"),
        yaxis=list(title="Ulcer Area (cm²)", gridcolor="#ecf0f1"),
        legend=list(orientation="h"),
        plot_bgcolor="#fdfdfd", paper_bgcolor="#fdfdfd",
        title=list(text="Venous Ulcer Area (C6)", x=0, font=list(size=14))
      )
  })

  output$vcss_bar <- renderPlotly({
    end_wk <- input$ep_weeks
    sc_names <- names(scenarios)
    pct_red <- sapply(sc_names, function(sn) {
      sc <- scenarios[[sn]]
      vcss_end <- tail(vcss_model(c(0, end_wk), ceap(), vcss0(), dvt(), age(),
        mpff=sc$mpff, ptx=sc$ptx, lmwh=sc$lmwh,
        ruto=sc$ruto, compression=sc$comp), 1)
      round((1 - vcss_end / max(vcss0(), 0.01)) * 100, 1)
    })
    cols <- unname(scenario_colors[sc_names])
    plot_ly(x=sc_names, y=pct_red, type="bar",
      marker=list(color=cols),
      text=paste0(pct_red, "%"), textposition="outside") %>%
      layout(
        xaxis=list(title="", tickangle=-30),
        yaxis=list(title="VCSS % Reduction at Week", range=c(0, 80)),
        plot_bgcolor="#fdfdfd", paper_bgcolor="#fdfdfd",
        title=list(text=paste("VCSS % Reduction at Week", end_wk), x=0, font=list(size=13))
      )
  })

  # ── TAB 6: Scenario Comparison ─────────────────────────────────────────────
  sc_selected <- reactive({
    sel <- input$sc_scenarios
    if (length(sel) == 0) "Untreated" else sel
  })

  sc_weeks_vec <- reactive(seq(0, input$sc_weeks, by = 0.5))

  outcome_fn <- reactive({
    switch(input$sc_outcome,
      "VCSS"      = function(wk, sc) vcss_model(wk, ceap(), vcss0(), dvt(), age(),
                      mpff=sc$mpff, ptx=sc$ptx, lmwh=sc$lmwh,
                      ruto=sc$ruto, compression=sc$comp),
      "CIVIQ-20"  = function(wk, sc) civiq_model(wk, ceap(), vcss0(),
                      mpff=sc$mpff, ptx=sc$ptx, lmwh=sc$lmwh,
                      ruto=sc$ruto, compression=sc$comp),
      "AVP (mmHg)"= function(wk, sc) avp_model(wk, ceap(), dvt(), inco(),
                      mpff=sc$mpff, ptx=sc$ptx, lmwh=sc$lmwh,
                      ruto=sc$ruto, compression=sc$comp),
      "Edema (cm)"= function(wk, sc) edema_model(wk, ceap(), bmi(), age(),
                      mpff=sc$mpff, ruto=sc$ruto, lmwh=sc$lmwh, compression=sc$comp)
    )
  })

  output$sc_time_plot <- renderPlotly({
    wk <- sc_weeks_vec()
    fn <- outcome_fn()
    p  <- plot_ly()
    for (sn in sc_selected()) {
      sc   <- scenarios[[sn]]
      yval <- fn(wk, sc)
      p <- add_trace(p, x=wk, y=yval, type="scatter", mode="lines",
        name=sn, line=list(color=scenario_colors[sn], width=2.5))
    }
    layout(p,
      xaxis=list(title="Weeks", gridcolor="#ecf0f1"),
      yaxis=list(title=input$sc_outcome, gridcolor="#ecf0f1"),
      legend=list(orientation="h", y=-0.25),
      plot_bgcolor="#fdfdfd", paper_bgcolor="#fdfdfd",
      margin=list(b=80)
    )
  })

  output$sc_summary_table <- renderTable({
    end_wk <- input$sc_weeks
    fn     <- outcome_fn()
    rows <- lapply(sc_selected(), function(sn) {
      sc  <- scenarios[[sn]]
      y0  <- fn(0,      sc)[1]
      y_end <- fn(end_wk, sc)
      y_end <- if (length(y_end) == 1) y_end else tail(y_end, 1)
      pct <- round((y0 - y_end) / max(y0, 0.001) * 100, 1)
      data.frame(Scenario=sn,
        `Baseline` = round(y0, 2),
        `At Endpoint` = round(y_end, 2),
        `% Change` = pct,
        check.names=FALSE)
    })
    do.call(rbind, rows)
  }, striped=TRUE, bordered=TRUE)

  output$radar_plot <- renderPlotly({
    end_wk <- input$sc_weeks
    dims   <- c("VCSS Reduction", "CIVIQ Improvement", "AVP Reduction",
                "Oedema Reduction", "Inflammation")
    sn_sel <- sc_selected()

    p <- plot_ly(type="scatterpolar", fill="toself")
    for (sn in sn_sel) {
      sc <- scenarios[[sn]]
      vcss_end  <- tail(vcss_model(c(0, end_wk), ceap(), vcss0(), dvt(), age(),
        mpff=sc$mpff, ptx=sc$ptx, lmwh=sc$lmwh, ruto=sc$ruto, compression=sc$comp), 1)
      civiq_end <- tail(civiq_model(c(0, end_wk), ceap(), vcss0(),
        mpff=sc$mpff, ptx=sc$ptx, lmwh=sc$lmwh, ruto=sc$ruto, compression=sc$comp), 1)
      avp_end   <- tail(avp_model(c(0, end_wk), ceap(), dvt(), inco(),
        mpff=sc$mpff, ptx=sc$ptx, lmwh=sc$lmwh, ruto=sc$ruto, compression=sc$comp), 1)
      edema_end <- tail(edema_model(c(0, end_wk), ceap(), bmi(), age(),
        mpff=sc$mpff, ruto=sc$ruto, lmwh=sc$lmwh, compression=sc$comp), 1)
      inf_end   <- tail(inflam_model(c(0, end_wk), ceap(),
        treatments = c(
          if (sc$mpff) "MPFF", if (sc$ptx)  "PTX",
          if (sc$lmwh) "LMWH", if (sc$ruto) "Rutosides",
          if (sc$comp) "Compression"
        ))$Inflammatory_Composite, 1)

      vcss0_v   <- max(vcss0(), 0.01)
      avp0_v    <- avp_model(0, ceap(), dvt(), inco())[1]
      edema0_v  <- edema_model(0, ceap(), bmi(), age())[1]
      civiq0_v  <- civiq_model(0, ceap(), vcss0())[1]

      r_vals <- c(
        pmax((1 - vcss_end  / vcss0_v)   * 100, 0),
        pmax((civiq0_v - civiq_end) / max(civiq0_v - 20, 0.01) * 100, 0),
        pmax((1 - avp_end   / max(avp0_v,   0.01)) * 100, 0),
        pmax((1 - edema_end / max(edema0_v, 0.01)) * 100, 0),
        pmax((1 - inf_end   / max((ceap()-1)/5, 0.01)) * 100, 0)
      )
      r_vals <- c(r_vals, r_vals[1])  # close polygon
      dims_closed <- c(dims, dims[1])

      p <- add_trace(p,
        r    = pmin(r_vals, 100),
        theta= dims_closed,
        name = sn,
        line = list(color=scenario_colors[sn]),
        fillcolor = paste0(gsub("#","",scenario_colors[sn]), "33")
      )
    }
    layout(p,
      polar = list(radialaxis = list(visible=TRUE, range=c(0,100), ticksuffix="%")),
      legend = list(orientation="h", y=-0.15),
      paper_bgcolor="#fdfdfd",
      margin = list(t=30)
    )
  })

  output$sc_bar_plot <- renderPlotly({
    end_wk <- input$sc_weeks
    fn     <- outcome_fn()
    sn_sel <- sc_selected()
    pct_vals <- sapply(sn_sel, function(sn) {
      sc <- scenarios[[sn]]
      y0  <- fn(0, sc)[1]
      y_end_vec <- fn(end_wk, sc)
      y_end <- if (length(y_end_vec) == 1) y_end_vec else tail(y_end_vec, 1)
      round((y0 - y_end) / max(y0, 0.001) * 100, 1)
    })
    cols <- unname(scenario_colors[sn_sel])
    plot_ly(x=sn_sel, y=pct_vals, type="bar",
      marker=list(color=cols),
      text=paste0(pct_vals, "%"), textposition="outside") %>%
      layout(
        xaxis=list(title="", tickangle=-30),
        yaxis=list(title=paste0(input$sc_outcome, " % Improvement"), range=c(0, 100)),
        plot_bgcolor="#fdfdfd", paper_bgcolor="#fdfdfd",
        title=list(text=paste("% Improvement in", input$sc_outcome, "at Week", end_wk),
                   x=0, font=list(size=13)),
        margin=list(b=120)
      )
  })

}  # end server

# =============================================================================
# RUN
# =============================================================================
shinyApp(ui = ui, server = server)
