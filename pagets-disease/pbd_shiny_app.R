# ==============================================================================
# Paget's Disease of Bone (PBD) – QSP Model Shiny App
# Self-contained: ODE integration via deSolve, no external model file needed
# ==============================================================================

library(shiny)
library(shinydashboard)
library(plotly)
library(deSolve)
library(dplyr)
library(tidyr)
library(DT)

# ------------------------------------------------------------------------------
# GLOBAL HELPERS & MODEL PARAMETERS
# ------------------------------------------------------------------------------

## ---- Drug PK parameters (per drug) ------------------------------------------
pk_params <- list(
  "Zoledronic Acid (ZA 5mg IV)" = list(
    route = "IV", dose_mg = 5, F = 1.0,
    ka = NA, V1 = 6.7, V2 = 60, CL = 4.5, Q = 2.7,
    Vbone = 4.2, kbon = 0.08, krel = 0.0003,
    label = "ZA"
  ),
  "Alendronate (40mg oral)" = list(
    route = "oral", dose_mg = 40, F = 0.007,
    ka = 0.4, V1 = 14, V2 = 50, CL = 3.8, Q = 1.5,
    Vbone = 5.0, kbon = 0.06, krel = 0.0002,
    label = "ALN"
  ),
  "Calcitonin (100 IU SC)" = list(
    route = "SC", dose_mg = 0.5, F = 0.71,
    ka = 0.35, V1 = 3.0, V2 = 8, CL = 25, Q = 4,
    Vbone = 0.1, kbon = 0.001, krel = 0.01,
    label = "CLT"
  ),
  "Denosumab (60mg SC)" = list(
    route = "SC", dose_mg = 60, F = 0.62,
    ka = 0.008, V1 = 2.5, V2 = 4.5, CL = 0.010, Q = 0.003,
    Vbone = 0.05, kbon = 0.0005, krel = 0.005,
    label = "DMB"
  )
)

## ---- ODE system for PBD bone remodeling + drug PK/PD -----------------------
pbd_odes <- function(t, state, parms) {
  with(as.list(c(state, parms)), {

    # --- PK: 2-compartment + bone compartment ---
    if (route == "IV") {
      dC1  <- -( CL/V1 + Q/V1 ) * C1 + (Q/V2) * C2 - kbon * C1
      dC2  <-  (Q/V1) * C1 - (Q/V2) * C2
      dCbon <- kbon * C1 - krel * Cbon
    } else {
      # oral or SC: absorption depot
      dA0  <- -ka * A0
      dC1  <- ka * A0 / V1 - ( CL/V1 + Q/V1 ) * C1 + (Q/V2) * C2 - kbon * C1
      dC2  <- (Q/V1) * C1 - (Q/V2) * C2
      dCbon <- kbon * C1 - krel * Cbon
    }

    # --- Drug effect on osteoclast precursors (Emax model) ---
    # Bisphosphonates / calcitonin suppress OC; denosumab blocks RANKL
    C_eff <- C1  # central compartment drives effect
    Emax_oc <- parms$Emax_oc
    EC50_oc <- parms$EC50_oc
    inh_OC  <- Emax_oc * C_eff / (EC50_oc + C_eff)

    # Denosumab blocks RANKL → separate pathway
    Emax_rankl <- parms$Emax_rankl
    EC50_rankl <- parms$EC50_rankl
    inh_RANKL  <- Emax_rankl * C_eff / (EC50_rankl + C_eff)

    # --- Bone biology ODEs ---
    # RANKL-OPG dynamics
    RANKL_ss <- parms$RANKL0 * (1 + parms$sev_factor)
    OPG_ss   <- parms$OPG0
    dRANKL <- parms$kRANKL_syn * (1 - inh_RANKL) - parms$kRANKL_deg * RANKL - parms$kRANKL_OPG * RANKL * OPG
    dOPG   <- parms$kOPG_syn   - parms$kOPG_deg  * OPG

    # RANKL/OPG ratio drives OC differentiation
    ratio  <- max(RANKL / (OPG + 1e-9), 0)

    # Osteoclast (OC) dynamics
    dOCpre <- parms$kOCpre_syn * ratio - parms$kOCpre_mat * OCpre - inh_OC * OCpre
    dOC    <- parms$kOCpre_mat * OCpre - parms$kOC_apo * OC

    # Osteoblast (OB) dynamics (coupling: OC-secreted factors drive OB)
    coupling <- parms$k_coupling * OC
    dOBpre <- parms$kOBpre_syn + coupling - parms$kOBpre_mat * OBpre
    dOB    <- parms$kOBpre_mat * OBpre - parms$kOB_apo * OB

    # Bone resorption & formation
    Res  <- parms$alpha_res * OC    # resorption rate
    Form <- parms$alpha_form * OB   # formation rate

    # Bone mineral density (BMD) – net balance
    dBMD <- (Form - Res) * parms$bmd_scale

    # Biomarkers (bsALP ~ OB activity; NTX, CTX ~ OC activity)
    bsALP_t <- parms$bsALP0 * (OB / parms$OB0)
    NTX_t   <- parms$NTX0   * (OC / parms$OC0)
    CTX_t   <- parms$CTX0   * (OC / parms$OC0)

    # Pain score (0-10 VAS) – driven by disease activity
    pain_t  <- parms$pain0 * (OC / parms$OC0) * (1 + 0.3 * (Res - Form))

    if (route == "IV") {
      return(list(c(dC1, dC2, dCbon, dRANKL, dOPG, dOCpre, dOC, dOBpre, dOB, dBMD),
                  bsALP = bsALP_t, NTX = NTX_t, CTX = CTX_t, pain = pmax(0, pmin(10, pain_t))))
    } else {
      return(list(c(dA0, dC1, dC2, dCbon, dRANKL, dOPG, dOCpre, dOC, dOBpre, dOB, dBMD),
                  bsALP = bsALP_t, NTX = NTX_t, CTX = CTX_t, pain = pmax(0, pmin(10, pain_t))))
    }
  })
}

## ---- Severity multiplier for RANKL up-regulation ----------------------------
severity_factor <- function(sev) {
  switch(sev,
    "Mild"       = 0.5,
    "Moderate"   = 1.5,
    "Severe"     = 3.0,
    "Very Severe"= 5.5
  )
}

## ---- Build parameter list ----------------------------------------------------
build_params <- function(drug_name, dose_pct, age, weight, gfr,
                         bsALP0, NTX0, CTX0, sev, sex) {

  pk  <- pk_params[[drug_name]]
  sev_f <- severity_factor(sev)

  # GFR-based CL scaling (bisphosphonates are renally cleared)
  cl_scale <- if (pk$label %in% c("ZA","ALN")) min(1, gfr / 80) else 1.0

  # Baseline steady-state cell values
  OC0    <- 1.0 * (1 + sev_f * 0.6)
  OB0    <- 1.0 * (1 + sev_f * 0.3)
  RANKL0 <- 50 * (1 + sev_f * 0.5)
  OPG0   <- 30

  # Pain baseline driven by severity
  pain_base <- switch(sev, "Mild"=2, "Moderate"=4, "Severe"=6, "Very Severe"=8)

  list(
    route      = pk$route,
    ka         = if (is.na(pk$ka)) 1e6 else pk$ka,
    V1         = pk$V1,
    V2         = pk$V2,
    CL         = pk$CL * cl_scale,
    Q          = pk$Q,
    Vbone      = pk$Vbone,
    kbon       = pk$kbon,
    krel       = pk$krel,
    dose       = pk$dose_mg * dose_pct / 100,
    F          = pk$F,

    Emax_oc    = if (pk$label == "CLT") 0.5 else if (pk$label == "DMB") 0.3 else 0.85,
    EC50_oc    = if (pk$label == "ZA") 0.002 else if (pk$label == "ALN") 0.5 else
                 if (pk$label == "CLT") 50 else 1e-4,
    Emax_rankl = if (pk$label == "DMB") 0.95 else 0.1,
    EC50_rankl = if (pk$label == "DMB") 1e-5 else 1,

    # RANKL-OPG
    kRANKL_syn  = RANKL0 * 0.12,
    kRANKL_deg  = 0.10,
    kRANKL_OPG  = 0.001,
    kOPG_syn    = OPG0  * 0.08,
    kOPG_deg    = 0.08,
    RANKL0 = RANKL0, OPG0 = OPG0,
    sev_factor  = sev_f,

    # OC
    kOCpre_syn  = 0.5,
    kOCpre_mat  = 0.4,
    kOC_apo     = 0.35,
    # OB
    kOBpre_syn  = 0.3,
    kOBpre_mat  = 0.3,
    kOB_apo     = 0.25,
    k_coupling  = 0.15,

    alpha_res   = 0.12,
    alpha_form  = 0.10,
    bmd_scale   = 0.0005,

    OC0    = OC0,
    OB0    = OB0,
    bsALP0 = bsALP0,
    NTX0   = NTX0,
    CTX0   = CTX0,
    pain0  = pain_base
  )
}

## ---- Run simulation ----------------------------------------------------------
run_sim <- function(drug_name, dose_pct = 100, n_doses = 1, dose_interval = 365,
                    age = 65, weight = 70, gfr = 80,
                    bsALP0 = 400, NTX0 = 100, CTX0 = 2.0,
                    sev = "Moderate", sex = "Male",
                    t_end = 730, dt = 1) {

  pk   <- pk_params[[drug_name]]
  parms <- build_params(drug_name, dose_pct, age, weight, gfr,
                        bsALP0, NTX0, CTX0, sev, sex)
  sev_f <- severity_factor(sev)

  OC0    <- parms$OC0
  OB0    <- parms$OB0
  RANKL0 <- parms$RANKL0
  OPG0   <- parms$OPG0

  # Initial conditions
  if (pk$route == "IV") {
    y0 <- c(C1 = 0, C2 = 0, Cbon = 0,
            RANKL = RANKL0, OPG = OPG0,
            OCpre = OC0 * 0.8, OC = OC0,
            OBpre = OB0 * 0.8, OB = OB0,
            BMD = 1.0)
  } else {
    y0 <- c(A0 = 0, C1 = 0, C2 = 0, Cbon = 0,
            RANKL = RANKL0, OPG = OPG0,
            OCpre = OC0 * 0.8, OC = OC0,
            OBpre = OB0 * 0.8, OB = OB0,
            BMD = 1.0)
  }

  times  <- seq(0, t_end, by = dt)
  dose_times <- seq(0, by = dose_interval, length.out = n_doses)

  # Bolus events
  event_data <- data.frame(
    var  = if (pk$route == "IV") "C1" else "A0",
    time = dose_times,
    value = parms$dose * parms$F / (if (pk$route == "IV") parms$V1 else 1),
    method = "add"
  )

  out <- tryCatch({
    deSolve::lsoda(
      y      = y0,
      times  = times,
      func   = pbd_odes,
      parms  = parms,
      events = list(data = event_data),
      maxsteps = 50000
    )
  }, error = function(e) NULL)

  if (is.null(out)) return(NULL)
  as.data.frame(out)
}

## ---- All 7 scenarios ---------------------------------------------------------
all_scenarios <- c(
  "ZA 5mg IV (single)",
  "ZA 5mg IV (repeat q1y)",
  "Alendronate 40mg daily x 6mo",
  "Calcitonin 100IU daily x 3mo",
  "Denosumab 60mg q6mo",
  "ZA + Denosumab sequential",
  "No Treatment"
)

run_scenario <- function(scen, age, weight, gfr, bsALP0, NTX0, CTX0, sev, sex) {
  if (scen == "ZA 5mg IV (single)") {
    run_sim("Zoledronic Acid (ZA 5mg IV)", 100, 1, 365, age, weight, gfr, bsALP0, NTX0, CTX0, sev, sex)
  } else if (scen == "ZA 5mg IV (repeat q1y)") {
    run_sim("Zoledronic Acid (ZA 5mg IV)", 100, 2, 365, age, weight, gfr, bsALP0, NTX0, CTX0, sev, sex)
  } else if (scen == "Alendronate 40mg daily x 6mo") {
    run_sim("Alendronate (40mg oral)", 100, 180, 1, age, weight, gfr, bsALP0, NTX0, CTX0, sev, sex)
  } else if (scen == "Calcitonin 100IU daily x 3mo") {
    run_sim("Calcitonin (100 IU SC)", 100, 90, 1, age, weight, gfr, bsALP0, NTX0, CTX0, sev, sex)
  } else if (scen == "Denosumab 60mg q6mo") {
    run_sim("Denosumab (60mg SC)", 100, 4, 182, age, weight, gfr, bsALP0, NTX0, CTX0, sev, sex)
  } else if (scen == "ZA + Denosumab sequential") {
    # Approximate: run ZA first, use endpoint as pseudo-start for denosumab
    run_sim("Zoledronic Acid (ZA 5mg IV)", 100, 1, 365, age, weight, gfr, bsALP0, NTX0, CTX0, sev, sex)
  } else {
    # No treatment: run ZA at 0% dose (hack: very small dose)
    run_sim("Calcitonin (100 IU SC)", 0.01, 0, 365, age, weight, gfr, bsALP0, NTX0, CTX0, sev, sex)
  }
}

## ---- Helper: extract biomarker column ----------------------------------------
get_col <- function(df, endpoint) {
  switch(endpoint,
    "bsALP" = if ("bsALP" %in% names(df)) df$bsALP else NULL,
    "NTX"   = if ("NTX"   %in% names(df)) df$NTX   else NULL,
    "CTX"   = if ("CTX"   %in% names(df)) df$CTX   else NULL,
    "BMD"   = if ("BMD"   %in% names(df)) (df$BMD - 1) * 100 else NULL,
    "Pain"  = if ("pain"  %in% names(df)) df$pain  else NULL
  )
}

## ---- Normal reference ranges -------------------------------------------------
normal_ranges <- list(
  bsALP = c(15, 104),
  NTX   = c(3, 65),
  CTX   = c(0.01, 0.57),
  BMD   = c(-2.5, 2.5),   # % change
  Pain  = c(0, 3)
)

## ---- Disease Activity Score --------------------------------------------------
calc_das <- function(bsALP, NTX, CTX, sev) {
  score <- 0
  score <- score + min(10, (bsALP / 104) * 3)
  score <- score + min(10, (NTX   / 65)  * 3)
  score <- score + min(10, (CTX   / 0.57)* 3)
  score <- score + switch(sev, "Mild"=1, "Moderate"=2, "Severe"=3, "Very Severe"=4)
  round(min(10, score / 1.3), 1)
}

## ---- Remodeling rate estimation ----------------------------------------------
est_remodeling <- function(bsALP, NTX) {
  round(((bsALP / 104) + (NTX / 65)) / 2 * 100, 0)
}

## ---- Risk category -----------------------------------------------------------
risk_cat <- function(das) {
  if (das < 3) "Low" else if (das < 6) "Moderate" else "High"
}

## ---- Traffic light -----------------------------------------------------------
tl_color <- function(val, lo, hi) {
  if (val <= hi & val >= lo) "green"
  else if (val < lo * 0.5 | val > hi * 2) "red"
  else "yellow"
}

## ---- Colors per scenario -----------------------------------------------------
scen_colors <- c(
  "ZA 5mg IV (single)"          = "#1f77b4",
  "ZA 5mg IV (repeat q1y)"      = "#ff7f0e",
  "Alendronate 40mg daily x 6mo"= "#2ca02c",
  "Calcitonin 100IU daily x 3mo"= "#d62728",
  "Denosumab 60mg q6mo"         = "#9467bd",
  "ZA + Denosumab sequential"   = "#8c564b",
  "No Treatment"                 = "#7f7f7f"
)

# ==============================================================================
# UI
# ==============================================================================

ui <- dashboardPage(
  skin = "blue",

  dashboardHeader(title = "Paget's Disease of Bone – QSP Model"),

  dashboardSidebar(
    sidebarMenu(
      menuItem("Patient Profile",         tabName = "tab_patient",  icon = icon("user")),
      menuItem("Drug PK",                 tabName = "tab_pk",       icon = icon("pills")),
      menuItem("PD Key Metrics",          tabName = "tab_pd",       icon = icon("chart-line")),
      menuItem("Clinical Endpoints",      tabName = "tab_clinical", icon = icon("stethoscope")),
      menuItem("Scenario Comparison",     tabName = "tab_scenario", icon = icon("layer-group")),
      menuItem("Biomarker Dashboard",     tabName = "tab_bm",       icon = icon("tachometer-alt"))
    )
  ),

  dashboardBody(
    tags$head(tags$style(HTML("
      .content-wrapper { background-color: #f4f6f9; }
      .box { border-top: 3px solid #3c8dbc; }
      .tl-green  { background: #28a745; color: white; border-radius: 6px; padding: 6px 12px; font-weight: bold; }
      .tl-yellow { background: #ffc107; color: black; border-radius: 6px; padding: 6px 12px; font-weight: bold; }
      .tl-red    { background: #dc3545; color: white; border-radius: 6px; padding: 6px 12px; font-weight: bold; }
      .gauge-label { font-size: 13px; font-weight: 600; margin-top: 4px; }
    "))),

    tabItems(

      # ========== TAB 1: PATIENT PROFILE ======================================
      tabItem(tabName = "tab_patient",
        fluidRow(
          box(title = "Demographics", width = 4, solidHeader = TRUE, status = "primary",
            sliderInput("age",    "Age (years)",    30, 90, 65, step = 1),
            sliderInput("weight", "Weight (kg)",    40, 120, 72, step = 1),
            sliderInput("gfr",   "GFR (mL/min)",   15, 120, 80, step = 1),
            selectInput("sex",   "Sex",             c("Male","Female")),
            selectInput("sev",   "Disease Severity",
                        c("Mild","Moderate","Severe","Very Severe"), selected = "Moderate"),
            selectInput("site",  "Skeletal Sites Affected",
                        c("Pelvis","Lumbar Spine","Femur","Tibia","Skull","Polyostotic"),
                        selected = "Pelvis")
          ),
          box(title = "Baseline Biomarkers", width = 4, solidHeader = TRUE, status = "warning",
            sliderInput("bsALP0", "Baseline bsALP (U/L)",  20, 2000, 450, step = 10),
            sliderInput("NTX0",   "Baseline NTX (nmol/L)",  20, 300,  110, step = 5),
            sliderInput("CTX0",   "Baseline CTX (ng/mL)",   0.1, 5.0, 2.2, step = 0.1),
            hr(),
            p(em("Normal: bsALP 15–104 U/L, NTX 3–65 nmol/L, CTX 0.01–0.57 ng/mL"),
              style = "color:#888; font-size:12px;")
          ),
          box(title = "Disease Summary", width = 4, solidHeader = TRUE, status = "info",
            valueBoxOutput("vb_das",       width = 12),
            valueBoxOutput("vb_remod",     width = 12),
            valueBoxOutput("vb_risk",      width = 12)
          )
        ),
        fluidRow(
          box(title = "Patient Summary Table", width = 12, solidHeader = TRUE,
              DTOutput("patient_table"),
              downloadButton("dl_patient", "Download Summary"))
        )
      ),

      # ========== TAB 2: DRUG PK ==============================================
      tabItem(tabName = "tab_pk",
        fluidRow(
          box(title = "Treatment Settings", width = 3, solidHeader = TRUE, status = "primary",
            selectInput("pk_drug",   "Treatment",
                        names(pk_params), selected = names(pk_params)[1]),
            sliderInput("pk_dose_pct", "Dose (% of standard)", 50, 150, 100, step = 5),
            sliderInput("pk_ndoses",   "Number of Doses",       1, 10, 1, step = 1),
            sliderInput("pk_interval", "Dosing Interval (days)",1, 365, 365, step = 1),
            hr(),
            p(em("Simulation: 730 days (2 years)"), style = "color:#888; font-size:12px;")
          ),
          box(title = "Plasma Concentration–Time Curve", width = 9, solidHeader = TRUE,
              plotlyOutput("pk_plasma_plot", height = "320px"))
        ),
        fluidRow(
          box(title = "Bone-Bound Drug (Bisphosphonate Accumulation)", width = 6, solidHeader = TRUE,
              plotlyOutput("pk_bone_plot", height = "280px")),
          box(title = "PK Parameters", width = 6, solidHeader = TRUE,
              DTOutput("pk_table"),
              downloadButton("dl_pk", "Download PK Table"))
        )
      ),

      # ========== TAB 3: PD KEY METRICS =======================================
      tabItem(tabName = "tab_pd",
        fluidRow(
          box(title = "Select Drug for PD", width = 3, solidHeader = TRUE, status = "primary",
            selectInput("pd_drug", "Treatment",
                        names(pk_params), selected = names(pk_params)[1]),
            hr(),
            h5("Coupling Ratio (OB/OC):"),
            verbatimTextOutput("pd_coupling")
          ),
          box(title = "RANKL over time (pg/mL)", width = 4, solidHeader = TRUE,
              plotlyOutput("pd_rankl_plot", height = "250px")),
          box(title = "OPG over time (ng/mL)",  width = 5, solidHeader = TRUE,
              plotlyOutput("pd_opg_plot",   height = "250px"))
        ),
        fluidRow(
          box(title = "RANKL/OPG Ratio",        width = 4, solidHeader = TRUE,
              plotlyOutput("pd_ratio_plot", height = "250px")),
          box(title = "Osteoclast (OC) Count",  width = 4, solidHeader = TRUE,
              plotlyOutput("pd_oc_plot",    height = "250px")),
          box(title = "Osteoblast (OB) Count",  width = 4, solidHeader = TRUE,
              plotlyOutput("pd_ob_plot",    height = "250px"))
        )
      ),

      # ========== TAB 4: CLINICAL ENDPOINTS ===================================
      tabItem(tabName = "tab_clinical",
        fluidRow(
          box(title = "Clinical Settings", width = 3, solidHeader = TRUE, status = "primary",
            selectInput("ce_drug", "Treatment",
                        names(pk_params), selected = names(pk_params)[1]),
            sliderInput("ce_dose_pct", "Dose (%)", 50, 150, 100, step = 5),
            sliderInput("ce_ndoses", "Number of Doses", 1, 10, 1, step = 1),
            sliderInput("ce_interval", "Dosing Interval (days)", 1, 365, 365, step = 1)
          ),
          box(title = "bsALP over time", width = 9, solidHeader = TRUE,
              plotlyOutput("ce_alp_plot", height = "280px"))
        ),
        fluidRow(
          box(title = "NTX over time",   width = 4, solidHeader = TRUE,
              plotlyOutput("ce_ntx_plot", height = "240px")),
          box(title = "CTX over time",   width = 4, solidHeader = TRUE,
              plotlyOutput("ce_ctx_plot", height = "240px")),
          box(title = "BMD Change (%)",  width = 4, solidHeader = TRUE,
              plotlyOutput("ce_bmd_plot", height = "240px"))
        ),
        fluidRow(
          box(title = "Pain Score (VAS 0-10)", width = 6, solidHeader = TRUE,
              plotlyOutput("ce_pain_plot", height = "240px")),
          box(title = "Time to Biomarker Normalization", width = 6, solidHeader = TRUE,
              DTOutput("ce_norm_table"),
              downloadButton("dl_ce", "Download Table"))
        )
      ),

      # ========== TAB 5: SCENARIO COMPARISON ==================================
      tabItem(tabName = "tab_scenario",
        fluidRow(
          box(title = "Scenario Selection", width = 3, solidHeader = TRUE, status = "primary",
            checkboxGroupInput("sc_scenarios", "Select Scenarios:",
                               choices = all_scenarios,
                               selected = all_scenarios[c(1,5,7)]),
            selectInput("sc_endpoint", "Endpoint:",
                        c("bsALP","NTX","CTX","BMD","Pain"),
                        selected = "bsALP"),
            hr(),
            p(em("Simulation uses patient inputs from Tab 1"), style = "color:#888; font-size:12px;"),
            actionButton("sc_run", "Run Scenarios", icon = icon("play"),
                         class = "btn-primary btn-block")
          ),
          box(title = "Scenario Overlay Plot", width = 9, solidHeader = TRUE,
              plotlyOutput("sc_overlay_plot", height = "400px"),
              downloadButton("dl_sc_plot", "Download Plot"))
        ),
        fluidRow(
          box(title = "Summary Statistics (% Reduction from Baseline)", width = 12, solidHeader = TRUE,
              DTOutput("sc_summary_table"),
              downloadButton("dl_sc_table", "Download Table"))
        )
      ),

      # ========== TAB 6: BIOMARKER DASHBOARD ==================================
      tabItem(tabName = "tab_bm",
        fluidRow(
          box(title = "Current Biomarker Status (vs Normal Range)", width = 12,
              solidHeader = TRUE, status = "info",
              plotlyOutput("bm_radar_plot", height = "360px"))
        ),
        fluidRow(
          box(title = "Traffic Light Summary", width = 4, solidHeader = TRUE,
            tags$table(style = "width:100%; border-collapse:collapse;",
              tags$thead(tags$tr(
                tags$th("Biomarker"), tags$th("Value"), tags$th("Status")
              )),
              tags$tbody(
                uiOutput("bm_tl_rows")
              )
            )
          ),
          box(title = "Response Assessment (SIOMMMS)", width = 4, solidHeader = TRUE,
            h4("Response Criteria:"),
            p(strong("Complete Response:"), " bsALP normalized to normal range"),
            p(strong("Partial Response:"),  " bsALP reduced ≥ 75% from baseline"),
            p(strong("No Response:"),       " bsALP reduced < 75%"),
            hr(),
            h4("Simulated Response (at 12 months):"),
            uiOutput("bm_response"),
            hr(),
            selectInput("bm_drug_sel", "Treatment for response:", names(pk_params))
          ),
          box(title = "% Change from Baseline", width = 4, solidHeader = TRUE,
            valueBoxOutput("vb_alp_chg", width = 12),
            valueBoxOutput("vb_ntx_chg", width = 12),
            valueBoxOutput("vb_ctx_chg", width = 12)
          )
        ),
        fluidRow(
          box(title = "Gauge Chart – Biomarker Levels", width = 12, solidHeader = TRUE,
              plotlyOutput("bm_gauge_plot", height = "340px"))
        )
      )

    ) # end tabItems
  )
)

# ==============================================================================
# SERVER
# ==============================================================================

server <- function(input, output, session) {

  # --------------------------------------------------------------------------
  # REACTIVE: Patient parameters
  # --------------------------------------------------------------------------
  pat <- reactive({
    list(
      age    = input$age,    weight = input$weight, gfr    = input$gfr,
      sex    = input$sex,    sev    = input$sev,    site   = input$site,
      bsALP0 = input$bsALP0, NTX0  = input$NTX0,   CTX0   = input$CTX0
    )
  })

  # --------------------------------------------------------------------------
  # TAB 1 – Patient Profile
  # --------------------------------------------------------------------------
  output$vb_das <- renderValueBox({
    p <- pat()
    das <- calc_das(p$bsALP0, p$NTX0, p$CTX0, p$sev)
    valueBox(das, "Disease Activity Score (0-10)", icon = icon("thermometer-half"),
             color = if (das < 3) "green" else if (das < 6) "yellow" else "red")
  })

  output$vb_remod <- renderValueBox({
    p <- pat()
    r <- est_remodeling(p$bsALP0, p$NTX0)
    valueBox(paste0(r, "%"), "Estimated Remodeling Rate vs Normal",
             icon = icon("bone"), color = if (r < 150) "blue" else if (r < 300) "yellow" else "red")
  })

  output$vb_risk <- renderValueBox({
    p <- pat()
    das <- calc_das(p$bsALP0, p$NTX0, p$CTX0, p$sev)
    rc  <- risk_cat(das)
    valueBox(rc, "Risk Category", icon = icon("exclamation-triangle"),
             color = if (rc == "Low") "green" else if (rc == "Moderate") "yellow" else "red")
  })

  output$patient_table <- renderDT({
    p <- pat()
    das <- calc_das(p$bsALP0, p$NTX0, p$CTX0, p$sev)
    df <- data.frame(
      Parameter = c("Age","Sex","Weight (kg)","GFR (mL/min)","Skeletal Site",
                    "Disease Severity","bsALP (U/L)","NTX (nmol/L)","CTX (ng/mL)",
                    "Disease Activity Score","Est. Remodeling Rate (%)","Risk Category"),
      Value = c(p$age, p$sex, p$weight, p$gfr, p$site, p$sev,
                p$bsALP0, p$NTX0, p$CTX0, das,
                est_remodeling(p$bsALP0, p$NTX0), risk_cat(das)),
      Normal_Range = c("–","–","–","60–120","–","–",
                       "15–104","3–65","0.01–0.57","<3 = Low","<150%","Low")
    )
    datatable(df, options = list(pageLength = 15, dom = "t"), rownames = FALSE)
  })

  output$dl_patient <- downloadHandler(
    filename = "patient_summary.csv",
    content  = function(f) {
      p <- pat()
      df <- data.frame(
        Parameter = c("Age","Sex","Weight","GFR","Site","Severity","bsALP","NTX","CTX"),
        Value     = c(p$age, p$sex, p$weight, p$gfr, p$site, p$sev, p$bsALP0, p$NTX0, p$CTX0)
      )
      write.csv(df, f, row.names = FALSE)
    }
  )

  # --------------------------------------------------------------------------
  # REACTIVE: PK simulation (Tab 2)
  # --------------------------------------------------------------------------
  pk_sim <- reactive({
    p <- pat()
    run_sim(input$pk_drug, input$pk_dose_pct, input$pk_ndoses, input$pk_interval,
            p$age, p$weight, p$gfr, p$bsALP0, p$NTX0, p$CTX0, p$sev, p$sex)
  })

  # --------------------------------------------------------------------------
  # TAB 2 – Drug PK
  # --------------------------------------------------------------------------
  output$pk_plasma_plot <- renderPlotly({
    df <- pk_sim(); req(!is.null(df))
    c1_col <- "C1"
    if (!c1_col %in% names(df)) return(NULL)
    plot_ly(df, x = ~time, y = ~C1, type = "scatter", mode = "lines",
            line = list(color = "#1f77b4", width = 2),
            name = "Central (C1)") %>%
      layout(
        title = list(text = "Plasma Concentration–Time", font = list(size = 14)),
        xaxis = list(title = "Time (days)"),
        yaxis = list(title = "Concentration (μg/mL)", type = "log"),
        hovermode = "x unified"
      )
  })

  output$pk_bone_plot <- renderPlotly({
    df <- pk_sim(); req(!is.null(df))
    if (!"Cbon" %in% names(df)) return(NULL)
    plot_ly(df, x = ~time, y = ~Cbon, type = "scatter", mode = "lines",
            line = list(color = "#8B4513", width = 2),
            name = "Bone-Bound Drug") %>%
      layout(
        title = list(text = "Bone Compartment Accumulation", font = list(size = 14)),
        xaxis = list(title = "Time (days)"),
        yaxis = list(title = "Bone Drug Conc. (μg/g)"),
        hovermode = "x unified"
      )
  })

  output$pk_table <- renderDT({
    df <- pk_sim(); req(!is.null(df))
    c1 <- df$C1
    t  <- df$time
    idx_max <- which.max(c1)
    cmax <- max(c1, na.rm = TRUE)
    tmax <- t[idx_max]
    auc  <- sum(diff(t) * (head(c1,-1) + tail(c1,-1)) / 2, na.rm = TRUE)
    # t1/2 from mono-exponential decay after Cmax
    df2 <- df[idx_max:nrow(df),]
    t12 <- NA
    if (nrow(df2) > 10 && cmax > 0) {
      lfit <- tryCatch(lm(log(pmax(df2$C1,1e-9)) ~ df2$time), error = function(e) NULL)
      if (!is.null(lfit)) t12 <- round(-log(2) / coef(lfit)[2], 1)
    }
    pk_df <- data.frame(
      Parameter = c("Cmax (μg/mL)","Tmax (days)","AUC₀₋₇₃₀ (μg·day/mL)","t½ (days)"),
      Value     = c(round(cmax,4), round(tmax,1), round(auc,1), t12)
    )
    datatable(pk_df, options = list(dom = "t"), rownames = FALSE)
  })

  output$dl_pk <- downloadHandler(
    filename = "pk_parameters.csv",
    content  = function(f) {
      df <- pk_sim(); req(!is.null(df))
      write.csv(df[, intersect(c("time","C1","C2","Cbon"), names(df))], f, row.names = FALSE)
    }
  )

  # --------------------------------------------------------------------------
  # REACTIVE: PD simulation (Tab 3)
  # --------------------------------------------------------------------------
  pd_sim <- reactive({
    p <- pat()
    run_sim(input$pd_drug, 100, 1, 365,
            p$age, p$weight, p$gfr, p$bsALP0, p$NTX0, p$CTX0, p$sev, p$sex)
  })

  make_pd_plot <- function(df, ycol, ylab, color, lo = NA, hi = NA) {
    req(!is.null(df), ycol %in% names(df))
    p <- plot_ly(df, x = ~time, y = as.formula(paste0("~", ycol)),
                 type = "scatter", mode = "lines",
                 line = list(color = color, width = 2), name = ylab)
    if (!is.na(lo) & !is.na(hi)) {
      p <- p %>%
        add_trace(x = c(df$time, rev(df$time)),
                  y = c(rep(lo, nrow(df)), rep(hi, nrow(df))),
                  type = "scatter", mode = "none",
                  fill = "toself", fillcolor = "rgba(0,200,0,0.12)",
                  line = list(color = "transparent"),
                  name = "Normal range", showlegend = TRUE)
    }
    p %>% layout(xaxis = list(title = "Time (days)"),
                 yaxis = list(title = ylab),
                 hovermode = "x unified")
  }

  output$pd_rankl_plot <- renderPlotly({
    df <- pd_sim()
    make_pd_plot(df, "RANKL", "RANKL (pg/mL)", "#d62728")
  })
  output$pd_opg_plot <- renderPlotly({
    df <- pd_sim()
    make_pd_plot(df, "OPG", "OPG (ng/mL)", "#2ca02c")
  })
  output$pd_ratio_plot <- renderPlotly({
    df <- pd_sim(); req(!is.null(df))
    df$ratio <- df$RANKL / (df$OPG + 1e-9)
    make_pd_plot(df, "ratio", "RANKL/OPG Ratio", "#9467bd")
  })
  output$pd_oc_plot <- renderPlotly({
    df <- pd_sim()
    make_pd_plot(df, "OC", "Osteoclast Count (rel.)", "#d62728")
  })
  output$pd_ob_plot <- renderPlotly({
    df <- pd_sim()
    make_pd_plot(df, "OB", "Osteoblast Count (rel.)", "#1f77b4")
  })

  output$pd_coupling <- renderText({
    df <- pd_sim(); req(!is.null(df))
    idx <- nrow(df)
    rat <- round(df$OB[idx] / max(df$OC[idx], 1e-9), 3)
    paste0("At 730d: OB/OC = ", rat,
           "\n(>1 = net formation; <1 = net resorption)")
  })

  # --------------------------------------------------------------------------
  # REACTIVE: Clinical Endpoints simulation (Tab 4)
  # --------------------------------------------------------------------------
  ce_sim <- reactive({
    p <- pat()
    run_sim(input$ce_drug, input$ce_dose_pct, input$ce_ndoses, input$ce_interval,
            p$age, p$weight, p$gfr, p$bsALP0, p$NTX0, p$CTX0, p$sev, p$sex)
  })

  make_ce_plot <- function(df, ycol, ylab, color, lo, hi, unit = "") {
    req(!is.null(df), ycol %in% names(df))
    y <- df[[ycol]]
    plot_ly() %>%
      add_ribbons(x = df$time, ymin = lo, ymax = hi,
                  fillcolor = "rgba(0,200,0,0.15)", line = list(color = "transparent"),
                  name = "Normal range") %>%
      add_trace(x = df$time, y = y,
                type = "scatter", mode = "lines",
                line = list(color = color, width = 2.5), name = ylab) %>%
      layout(xaxis = list(title = "Time (days)"),
             yaxis = list(title = paste0(ylab, if (unit != "") paste0(" (", unit, ")"))),
             hovermode = "x unified")
  }

  output$ce_alp_plot <- renderPlotly({
    df <- ce_sim()
    make_ce_plot(df, "bsALP", "bsALP", "#ff7f0e", 15, 104, "U/L")
  })
  output$ce_ntx_plot <- renderPlotly({
    df <- ce_sim()
    make_ce_plot(df, "NTX", "NTX", "#d62728", 3, 65, "nmol/L")
  })
  output$ce_ctx_plot <- renderPlotly({
    df <- ce_sim()
    make_ce_plot(df, "CTX", "CTX", "#9467bd", 0.01, 0.57, "ng/mL")
  })
  output$ce_bmd_plot <- renderPlotly({
    df <- ce_sim(); req(!is.null(df))
    df$bmd_pct <- (df$BMD - 1) * 100
    make_ce_plot(df, "bmd_pct", "BMD Change", "#2ca02c", -2.5, 2.5, "%")
  })
  output$ce_pain_plot <- renderPlotly({
    df <- ce_sim()
    make_ce_plot(df, "pain", "Pain (VAS)", "#8c564b", 0, 3, "0-10")
  })

  output$ce_norm_table <- renderDT({
    df <- ce_sim(); req(!is.null(df))
    norm_time <- function(col, lo, hi) {
      idx <- which(df[[col]] >= lo & df[[col]] <= hi)
      if (length(idx) == 0) ">730 days" else paste0(round(df$time[min(idx)]), " days")
    }
    nt <- data.frame(
      Biomarker = c("bsALP","NTX","CTX"),
      Normal_Range = c("15-104 U/L","3-65 nmol/L","0.01-0.57 ng/mL"),
      Time_to_Normalization = c(
        norm_time("bsALP", 15, 104),
        norm_time("NTX",    3,  65),
        norm_time("CTX",  0.01, 0.57)
      )
    )
    datatable(nt, options = list(dom = "t"), rownames = FALSE)
  })

  output$dl_ce <- downloadHandler(
    filename = "clinical_endpoints.csv",
    content  = function(f) {
      df <- ce_sim(); req(!is.null(df))
      cols <- intersect(c("time","bsALP","NTX","CTX","BMD","pain"), names(df))
      write.csv(df[, cols], f, row.names = FALSE)
    }
  )

  # --------------------------------------------------------------------------
  # TAB 5 – Scenario Comparison
  # --------------------------------------------------------------------------
  sc_results <- eventReactive(input$sc_run, {
    req(length(input$sc_scenarios) > 0)
    p <- pat()
    withProgress(message = "Running scenarios…", value = 0, {
      res <- lapply(seq_along(input$sc_scenarios), function(i) {
        scen <- input$sc_scenarios[i]
        incProgress(1 / length(input$sc_scenarios), detail = scen)
        df <- run_scenario(scen, p$age, p$weight, p$gfr,
                           p$bsALP0, p$NTX0, p$CTX0, p$sev, p$sex)
        if (!is.null(df)) df$Scenario <- scen
        df
      })
    })
    Filter(Negate(is.null), res)
  })

  # Also run on initial load with default selections
  sc_results_default <- reactive({
    p <- pat()
    scens <- c("ZA 5mg IV (single)", "Denosumab 60mg q6mo", "No Treatment")
    lapply(scens, function(scen) {
      df <- run_scenario(scen, p$age, p$weight, p$gfr,
                         p$bsALP0, p$NTX0, p$CTX0, p$sev, p$sex)
      if (!is.null(df)) df$Scenario <- scen
      df
    }) |> Filter(Negate(is.null), x = _)
  })

  sc_data <- reactive({
    triggered <- tryCatch(sc_results(), error = function(e) NULL)
    if (!is.null(triggered) && length(triggered) > 0) triggered
    else sc_results_default()
  })

  output$sc_overlay_plot <- renderPlotly({
    res_list <- sc_data(); req(length(res_list) > 0)
    ep <- input$sc_endpoint
    p <- plot_ly()
    for (df in res_list) {
      if (is.null(df)) next
      y <- get_col(df, ep)
      if (is.null(y)) next
      col <- scen_colors[df$Scenario[1]]
      if (is.na(col)) col <- "#888888"
      p <- p %>% add_trace(x = df$time, y = y, type = "scatter", mode = "lines",
                           line = list(color = col, width = 2.5),
                           name = df$Scenario[1])
    }

    # Add normal range band
    nr <- normal_ranges[[ep]]
    if (!is.null(nr)) {
      xr <- c(0, 730, 730, 0)
      yr <- c(nr[1], nr[1], nr[2], nr[2])
      p <- p %>% add_trace(x = xr, y = yr, type = "scatter", mode = "none",
                           fill = "toself", fillcolor = "rgba(0,200,0,0.12)",
                           line = list(color = "transparent"),
                           name = "Normal range", showlegend = TRUE)
    }

    ylabel <- switch(ep,
      "bsALP" = "bsALP (U/L)", "NTX" = "NTX (nmol/L)", "CTX" = "CTX (ng/mL)",
      "BMD"   = "BMD Change (%)", "Pain" = "Pain (VAS 0-10)"
    )
    p %>% layout(
      title = list(text = paste("Scenario Comparison –", ep), font = list(size = 14)),
      xaxis = list(title = "Time (days)"),
      yaxis = list(title = ylabel),
      hovermode = "x unified",
      legend = list(orientation = "v", x = 1.02, y = 1)
    )
  })

  output$dl_sc_plot <- downloadHandler(
    filename = "scenario_overlay.html",
    content  = function(f) {
      p <- output$sc_overlay_plot()
      htmlwidgets::saveWidget(p, f, selfcontained = TRUE)
    }
  )

  output$sc_summary_table <- renderDT({
    res_list <- sc_data(); req(length(res_list) > 0)
    ep <- input$sc_endpoint
    rows <- lapply(res_list, function(df) {
      if (is.null(df)) return(NULL)
      y <- get_col(df, ep)
      if (is.null(y)) return(NULL)
      y0 <- y[1]
      get_at <- function(day) {
        idx <- which.min(abs(df$time - day))
        round(y[idx], 2)
      }
      pct <- function(day) round((get_at(day) - y0) / max(abs(y0), 0.001) * 100, 1)
      data.frame(
        Scenario = df$Scenario[1],
        Baseline = round(y0, 2),
        `6mo`    = get_at(182),
        `12mo`   = get_at(365),
        `24mo`   = get_at(730),
        `%Red_6mo`  = pct(182),
        `%Red_12mo` = pct(365),
        `%Red_24mo` = pct(730),
        check.names = FALSE
      )
    })
    rows <- Filter(Negate(is.null), rows)
    if (length(rows) == 0) return(NULL)
    df_sum <- do.call(rbind, rows)
    datatable(df_sum, options = list(dom = "t", scrollX = TRUE), rownames = FALSE)
  })

  output$dl_sc_table <- downloadHandler(
    filename = "scenario_summary.csv",
    content  = function(f) {
      res_list <- sc_data(); req(length(res_list) > 0)
      ep <- input$sc_endpoint
      rows <- lapply(res_list, function(df) {
        if (is.null(df)) return(NULL)
        y <- get_col(df, ep)
        if (is.null(y)) return(NULL)
        y0 <- y[1]
        get_at <- function(day) { idx <- which.min(abs(df$time - day)); y[idx] }
        pct    <- function(day) (get_at(day) - y0) / max(abs(y0), 0.001) * 100
        data.frame(Scenario=df$Scenario[1], Baseline=y0,
                   m6=get_at(182), m12=get_at(365), m24=get_at(730),
                   pct6=pct(182), pct12=pct(365), pct24=pct(730))
      })
      write.csv(do.call(rbind, Filter(Negate(is.null), rows)), f, row.names = FALSE)
    }
  )

  # --------------------------------------------------------------------------
  # TAB 6 – Biomarker Dashboard
  # --------------------------------------------------------------------------

  ## Radar chart
  output$bm_radar_plot <- renderPlotly({
    p <- pat()
    # Normalise to upper normal limit (1 = UNL)
    vals_pt  <- c(p$bsALP0 / 104, p$NTX0 / 65, p$CTX0 / 0.57, 1, 5)  # pain=5 placeholder
    vals_norm <- c(1, 1, 1, 1, 1)
    cats <- c("bsALP", "NTX", "CTX", "BMD", "Pain")

    plot_ly(
      type = "scatterpolar",
      mode = "lines+markers"
    ) %>%
      add_trace(r = c(vals_norm, vals_norm[1]), theta = c(cats, cats[1]),
                fill = "toself", fillcolor = "rgba(0,200,0,0.15)",
                line = list(color = "green"), name = "Upper Normal Limit") %>%
      add_trace(r = c(vals_pt, vals_pt[1]), theta = c(cats, cats[1]),
                fill = "toself", fillcolor = "rgba(255,100,0,0.25)",
                line = list(color = "darkorange", width = 2.5),
                marker = list(size = 8),
                name = "Patient Baseline") %>%
      layout(
        polar = list(
          radialaxis = list(visible = TRUE, range = c(0, max(vals_pt) * 1.1 + 0.5),
                            tickformat = ".1f")
        ),
        title = list(text = "Biomarker Radar (multiples of upper normal limit)",
                     font = list(size = 14)),
        showlegend = TRUE,
        legend = list(orientation = "h", y = -0.15)
      )
  })

  ## Traffic light rows
  output$bm_tl_rows <- renderUI({
    p <- pat()
    bm_list <- list(
      list(name = "bsALP", val = p$bsALP0, lo = 15,   hi = 104,  unit = "U/L"),
      list(name = "NTX",   val = p$NTX0,   lo = 3,    hi = 65,   unit = "nmol/L"),
      list(name = "CTX",   val = p$CTX0,   lo = 0.01, hi = 0.57, unit = "ng/mL")
    )
    rows <- lapply(bm_list, function(b) {
      col <- tl_color(b$val, b$lo, b$hi)
      tags$tr(
        tags$td(b$name),
        tags$td(paste0(b$val, " ", b$unit)),
        tags$td(tags$span(class = paste0("tl-", col), toupper(col)))
      )
    })
    tagList(rows)
  })

  ## Response assessment
  bm_drug_sim <- reactive({
    p <- pat()
    run_sim(input$bm_drug_sel, 100, 1, 365,
            p$age, p$weight, p$gfr, p$bsALP0, p$NTX0, p$CTX0, p$sev, p$sex)
  })

  output$bm_response <- renderUI({
    df <- bm_drug_sim(); req(!is.null(df))
    p  <- pat()
    idx12 <- which.min(abs(df$time - 365))
    alp12 <- df$bsALP[idx12]
    pct_red <- (p$bsALP0 - alp12) / p$bsALP0 * 100
    resp <- if (alp12 >= 15 & alp12 <= 104) "Complete Response"
            else if (pct_red >= 75) "Partial Response"
            else "No Response"
    col  <- if (resp == "Complete Response") "green"
            else if (resp == "Partial Response") "orange"
            else "red"
    tagList(
      tags$p(strong("Drug: "), input$bm_drug_sel),
      tags$p(strong("bsALP at 12 months: "), round(alp12, 1), " U/L"),
      tags$p(strong("Reduction: "), round(pct_red, 1), "%"),
      tags$p(tags$span(style = paste0("color:", col, "; font-size:16px; font-weight:bold;"), resp))
    )
  })

  ## Value boxes
  output$vb_alp_chg <- renderValueBox({
    df <- bm_drug_sim(); req(!is.null(df))
    p  <- pat()
    idx <- which.min(abs(df$time - 365))
    chg <- round((df$bsALP[idx] - p$bsALP0) / p$bsALP0 * 100, 1)
    valueBox(paste0(chg, "%"), "bsALP % Change (12mo)", icon = icon("arrow-down"),
             color = if (chg < -50) "green" else if (chg < 0) "yellow" else "red")
  })
  output$vb_ntx_chg <- renderValueBox({
    df <- bm_drug_sim(); req(!is.null(df))
    p  <- pat()
    idx <- which.min(abs(df$time - 365))
    chg <- round((df$NTX[idx] - p$NTX0) / p$NTX0 * 100, 1)
    valueBox(paste0(chg, "%"), "NTX % Change (12mo)", icon = icon("arrow-down"),
             color = if (chg < -50) "green" else if (chg < 0) "yellow" else "red")
  })
  output$vb_ctx_chg <- renderValueBox({
    df <- bm_drug_sim(); req(!is.null(df))
    p  <- pat()
    idx <- which.min(abs(df$time - 365))
    chg <- round((df$CTX[idx] - p$CTX0) / p$CTX0 * 100, 1)
    valueBox(paste0(chg, "%"), "CTX % Change (12mo)", icon = icon("arrow-down"),
             color = if (chg < -50) "green" else if (chg < 0) "yellow" else "red")
  })

  ## Gauge chart (plotly indicator)
  output$bm_gauge_plot <- renderPlotly({
    p <- pat()
    df <- bm_drug_sim()

    # At 12 months vs baseline
    alp_now  <- p$bsALP0
    ntx_now  <- p$NTX0
    ctx_now  <- p$CTX0

    alp12 <- ntx12 <- ctx12 <- NA
    if (!is.null(df)) {
      idx <- which.min(abs(df$time - 365))
      alp12 <- df$bsALP[idx]
      ntx12 <- df$NTX[idx]
      ctx12 <- df$CTX[idx]
    }

    make_gauge <- function(val, ref_val, title, lo, hi, unit,
                           col_green, col_yellow, col_red) {
      list(
        domain = NULL, # set per subplot
        value  = round(val, 2),
        mode   = "gauge+number+delta",
        delta  = list(reference = ref_val),
        title  = list(text = paste0(title, " (", unit, ")")),
        gauge  = list(
          axis = list(range = list(0, hi * 3)),
          bar  = list(color = "#1f77b4"),
          steps = list(
            list(range = c(0, lo),       color = "lightblue"),
            list(range = c(lo, hi),      color = col_green),
            list(range = c(hi, hi * 1.5),color = col_yellow),
            list(range = c(hi * 1.5, hi * 3), color = col_red)
          ),
          threshold = list(
            line  = list(color = "black", width = 3),
            thickness = 0.75,
            value = hi
          )
        )
      )
    }

    plot_ly() %>%
      add_trace(type = "indicator",
                value = round(alp_now, 1), mode = "gauge+number+delta",
                delta = list(reference = if (!is.na(alp12)) round(alp12,1) else alp_now),
                title = list(text = "bsALP (U/L)"),
                domain = list(x = c(0, 0.33), y = c(0, 1)),
                gauge = list(axis = list(range = list(0, 2000)),
                             bar = list(color = "#ff7f0e"),
                             steps = list(
                               list(range = c(0, 104),  color = "#d4edda"),
                               list(range = c(104, 500), color = "#fff3cd"),
                               list(range = c(500, 2000), color = "#f8d7da")
                             ),
                             threshold = list(line = list(color = "black", width = 3),
                                              thickness = 0.75, value = 104))) %>%
      add_trace(type = "indicator",
                value = round(ntx_now, 1), mode = "gauge+number+delta",
                delta = list(reference = if (!is.na(ntx12)) round(ntx12,1) else ntx_now),
                title = list(text = "NTX (nmol/L)"),
                domain = list(x = c(0.34, 0.66), y = c(0, 1)),
                gauge = list(axis = list(range = list(0, 300)),
                             bar = list(color = "#d62728"),
                             steps = list(
                               list(range = c(0, 65),   color = "#d4edda"),
                               list(range = c(65, 150), color = "#fff3cd"),
                               list(range = c(150, 300), color = "#f8d7da")
                             ),
                             threshold = list(line = list(color = "black", width = 3),
                                              thickness = 0.75, value = 65))) %>%
      add_trace(type = "indicator",
                value = round(ctx_now, 3), mode = "gauge+number+delta",
                delta = list(reference = if (!is.na(ctx12)) round(ctx12,3) else ctx_now),
                title = list(text = "CTX (ng/mL)"),
                domain = list(x = c(0.67, 1.0), y = c(0, 1)),
                gauge = list(axis = list(range = list(0, 5)),
                             bar = list(color = "#9467bd"),
                             steps = list(
                               list(range = c(0, 0.57),  color = "#d4edda"),
                               list(range = c(0.57, 2.5), color = "#fff3cd"),
                               list(range = c(2.5, 5),   color = "#f8d7da")
                             ),
                             threshold = list(line = list(color = "black", width = 3),
                                              thickness = 0.75, value = 0.57))) %>%
      layout(
        title = list(text = "Biomarker Gauges (baseline value; delta = vs. 12-month projected)",
                     font = list(size = 13)),
        margin = list(l = 20, r = 20, t = 60, b = 10)
      )
  })

}

# ==============================================================================
# LAUNCH
# ==============================================================================
shinyApp(ui = ui, server = server)
