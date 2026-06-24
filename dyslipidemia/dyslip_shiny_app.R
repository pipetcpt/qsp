# ==============================================================================
# Dyslipidemia QSP Shiny Application
# Quantitative Systems Pharmacology Model for Dyslipidemia Management
# Author: Claude Code (CCR - Claude Code Routine)
# Date: 2026-06-17
# ==============================================================================
# Required packages:
#   shiny, bslib, ggplot2, plotly, DT, dplyr, tidyr, scales, shinycssloaders
# ==============================================================================

library(shiny)
library(bslib)
library(ggplot2)
library(plotly)
library(DT)
library(dplyr)
library(tidyr)
library(scales)

# ------------------------------------------------------------------------------
# HELPER FUNCTIONS & MODEL PARAMETERS
# ------------------------------------------------------------------------------

# --- Pooled Cohort Equations (PCE) for 10-year ASCVD risk ---
compute_pce_risk <- function(age, sex, race = "white", total_chol, hdl_c, sbp,
                              bp_treatment = FALSE, diabetes = FALSE, smoker = FALSE) {
  # Coefficients from Goff DC Jr et al. JACC 2014
  # Simplified white female / white male / AA female / AA male
  ln_age       <- log(age)
  ln_tchol     <- log(total_chol)
  ln_hdl       <- log(hdl_c)
  ln_sbp_tx    <- if (bp_treatment) log(sbp) else 0
  ln_sbp_notx  <- if (!bp_treatment) log(sbp) else 0
  dm           <- as.numeric(diabetes)
  smk          <- as.numeric(smoker)

  if (sex == "Female") {
    # White female
    sum_coef <- (
      -29.799 * ln_age +
      4.884   * ln_age^2 +
      13.540  * ln_tchol +
      -3.114  * ln_age * ln_tchol +
      -13.578 * ln_hdl +
      3.149   * ln_age * ln_hdl +
      2.019   * ln_sbp_tx +
      1.957   * ln_sbp_notx +
      7.574   * smk +
      -1.665  * ln_age * smk +
      0.661   * dm
    )
    baseline_surv <- 0.9665
    mean_coef     <- -29.799
  } else {
    # White male
    sum_coef <- (
      12.344  * ln_age +
      11.853  * ln_tchol +
      -2.664  * ln_age * ln_tchol +
      -7.990  * ln_hdl +
      1.769   * ln_age * ln_hdl +
      1.797   * ln_sbp_tx +
      1.764   * ln_sbp_notx +
      7.837   * smk +
      -1.795  * ln_age * smk +
      0.658   * dm
    )
    baseline_surv <- 0.9144
    mean_coef     <- 61.18
  }

  risk_10yr <- (1 - baseline_surv^exp(sum_coef - mean_coef)) * 100
  pmax(0.5, pmin(risk_10yr, 99))
}

# --- LDL-C target by risk category ---
get_ldl_target <- function(risk_category) {
  targets <- list(
    "Very High"    = list(esc = 55,  acc = 70,  label = "ESC: <55 | ACC/AHA: <70 mg/dL"),
    "High"         = list(esc = 70,  acc = 100, label = "ESC: <70 | ACC/AHA: <100 mg/dL"),
    "Intermediate" = list(esc = 100, acc = 130, label = "ESC: <100 | ACC/AHA: <130 mg/dL"),
    "Low"          = list(esc = 116, acc = 160, label = "ESC: <116 | ACC/AHA: <160 mg/dL")
  )
  targets[[risk_category]]
}

# --- Drug PK parameters ---
drug_pk_params <- list(
  "Atorvastatin 10mg"    = list(dose=10,   drug="Atorvastatin", Cmax=8.08, AUC=229, t_half=14, Tmax=1.0,  F=12,  MW=558.6, class="Statin"),
  "Atorvastatin 20mg"    = list(dose=20,   drug="Atorvastatin", Cmax=16.2, AUC=458, t_half=14, Tmax=1.0,  F=12,  MW=558.6, class="Statin"),
  "Atorvastatin 40mg"    = list(dose=40,   drug="Atorvastatin", Cmax=28.0, AUC=812, t_half=14, Tmax=1.2,  F=12,  MW=558.6, class="Statin"),
  "Atorvastatin 80mg"    = list(dose=80,   drug="Atorvastatin", Cmax=52.4, AUC=1540,t_half=14, Tmax=1.2,  F=12,  MW=558.6, class="Statin"),
  "Rosuvastatin 10mg"    = list(dose=10,   drug="Rosuvastatin", Cmax=6.25, AUC=155, t_half=19, Tmax=3.0,  F=20,  MW=481.5, class="Statin"),
  "Rosuvastatin 20mg"    = list(dose=20,   drug="Rosuvastatin", Cmax=12.0, AUC=311, t_half=19, Tmax=3.0,  F=20,  MW=481.5, class="Statin"),
  "Rosuvastatin 40mg"    = list(dose=40,   drug="Rosuvastatin", Cmax=22.6, AUC=587, t_half=19, Tmax=3.0,  F=20,  MW=481.5, class="Statin"),
  "Evolocumab 140mg Q2W" = list(dose=140,  drug="Evolocumab",   Cmax=13.8, AUC=1000,t_half=264,Tmax=72,   F=72,  MW=143500,class="PCSK9i"),
  "Evolocumab 420mg QM"  = list(dose=420,  drug="Evolocumab",   Cmax=59.0, AUC=3150,t_half=264,Tmax=120,  F=72,  MW=143500,class="PCSK9i"),
  "Alirocumab 75mg Q2W"  = list(dose=75,   drug="Alirocumab",   Cmax=6.04, AUC=684, t_half=336,Tmax=72,   F=85,  MW=146000,class="PCSK9i"),
  "Alirocumab 150mg Q2W" = list(dose=150,  drug="Alirocumab",   Cmax=11.3, AUC=1291,t_half=336,Tmax=72,   F=85,  MW=146000,class="PCSK9i"),
  "Inclisiran 284mg SC"  = list(dose=284,  drug="Inclisiran",   Cmax=508,  AUC=3420,t_half=9.4,Tmax=4,    F=NA,  MW=16100, class="siRNA"),
  "Ezetimibe 10mg"       = list(dose=10,   drug="Ezetimibe",    Cmax=3.4,  AUC=217, t_half=22, Tmax=4.0,  F=35,  MW=409.4, class="NPC1L1i"),
  "Bempedoic acid 180mg" = list(dose=180,  drug="Bempedoic acid",Cmax=11.8,AUC=236, t_half=21, Tmax=3.5,  F=NA,  MW=344.4, class="ACL-inhibitor")
)

# --- Drug PD (LDL-C reduction %) ---
drug_pd_ldl_reduction <- list(
  "Atorvastatin 10mg"    = list(ldl_pct=-37, hdl_pct=5,  tg_pct=-12, apob_pct=-34, pcsk9_pct=40),
  "Atorvastatin 20mg"    = list(ldl_pct=-43, hdl_pct=6,  tg_pct=-15, apob_pct=-39, pcsk9_pct=45),
  "Atorvastatin 40mg"    = list(ldl_pct=-49, hdl_pct=7,  tg_pct=-18, apob_pct=-44, pcsk9_pct=50),
  "Atorvastatin 80mg"    = list(ldl_pct=-55, hdl_pct=7,  tg_pct=-22, apob_pct=-49, pcsk9_pct=55),
  "Rosuvastatin 10mg"    = list(ldl_pct=-46, hdl_pct=8,  tg_pct=-15, apob_pct=-42, pcsk9_pct=50),
  "Rosuvastatin 20mg"    = list(ldl_pct=-52, hdl_pct=9,  tg_pct=-18, apob_pct=-47, pcsk9_pct=55),
  "Rosuvastatin 40mg"    = list(ldl_pct=-58, hdl_pct=10, tg_pct=-22, apob_pct=-52, pcsk9_pct=60),
  "Evolocumab 140mg Q2W" = list(ldl_pct=-60, hdl_pct=6,  tg_pct=-26, apob_pct=-56, pcsk9_pct=-98),
  "Evolocumab 420mg QM"  = list(ldl_pct=-59, hdl_pct=6,  tg_pct=-25, apob_pct=-55, pcsk9_pct=-98),
  "Alirocumab 75mg Q2W"  = list(ldl_pct=-48, hdl_pct=5,  tg_pct=-18, apob_pct=-43, pcsk9_pct=-97),
  "Alirocumab 150mg Q2W" = list(ldl_pct=-58, hdl_pct=6,  tg_pct=-24, apob_pct=-53, pcsk9_pct=-97),
  "Inclisiran 284mg SC"  = list(ldl_pct=-52, hdl_pct=4,  tg_pct=-15, apob_pct=-48, pcsk9_pct=-60),
  "Ezetimibe 10mg"       = list(ldl_pct=-18, hdl_pct=3,  tg_pct=-7,  apob_pct=-16, pcsk9_pct=0),
  "Bempedoic acid 180mg" = list(ldl_pct=-18, hdl_pct=3,  tg_pct=-8,  apob_pct=-19, pcsk9_pct=0)
)

# --- PK simulation function (1-compartment oral / SC) ---
simulate_pk <- function(drug_name, duration_weeks = 26, n_points = 500) {
  pk <- drug_pk_params[[drug_name]]
  t  <- seq(0, duration_weeks * 7, length.out = n_points)  # days

  if (pk$class %in% c("PCSK9i")) {
    # SC biologic: 1-compartment SC absorption
    # Determine dosing interval
    if (grepl("Q2W", drug_name)) {
      tau <- 14  # days
    } else if (grepl("QM", drug_name)) {
      tau <- 28
    } else {
      tau <- 14
    }
    ka   <- log(2) / 1.5  # absorption half-life ~1.5 days
    ke   <- log(2) / (pk$t_half / 24)  # convert hours to days
    V_F  <- pk$dose * 1000 / pk$Cmax  # rough volume/F in mL
    n_doses <- floor(duration_weeks * 7 / tau) + 1
    dose_times <- seq(0, (n_doses - 1) * tau, by = tau)

    conc <- numeric(n_points)
    for (td in dose_times) {
      t_after <- pmax(t - td, 0)
      conc <- conc + (pk$Cmax * ka / (ka - ke)) * (exp(-ke * t_after) - exp(-ka * t_after)) *
        (t >= td)
    }

  } else if (pk$class == "siRNA") {
    # Inclisiran: two-dose induction then Q6M
    # Doses at day 0, 90, 270, 450 (0, 90, then Q6M = 180 days)
    dose_times <- c(0, 90, 270, 450)
    dose_times <- dose_times[dose_times <= duration_weeks * 7]
    ke   <- log(2) / pk$t_half  # per day (t_half already in days)
    ka   <- log(2) / 0.25       # rapid absorption (Tmax ~4h ~ 0.17 day)
    conc <- numeric(n_points)
    for (td in dose_times) {
      t_after <- pmax(t - td, 0)
      conc <- conc + pk$Cmax * (exp(-ke * t_after) - exp(-ka * t_after)) /
        (1 - ke / ka) * (t >= td)
    }
    conc <- pmax(conc, 0)

  } else {
    # Oral drugs (statins, ezetimibe, bempedoic acid): QD dosing
    tau  <- 1  # daily
    ka   <- log(2) / (pk$Tmax / 2)  # rough ka from Tmax
    ke   <- log(2) / (pk$t_half / 24)  # convert to per-day
    n_doses <- floor(duration_weeks * 7) + 1
    dose_times <- seq(0, (n_doses - 1) * tau, by = tau)

    # Compute steady-state approximation after ~5 half-lives (~10 days for statins)
    # For display, show first 7 days detailed then steady-state envelope
    t_display <- seq(0, min(14, duration_weeks * 7), length.out = n_points)
    t <- t_display
    conc <- numeric(n_points)
    for (td in dose_times[1:min(15, length(dose_times))]) {
      t_after <- pmax(t - td, 0)
      conc <- conc + (pk$Cmax * ka / (ka - ke)) *
        (exp(-ke * t_after) - exp(-ka * t_after)) * (t >= td)
    }
  }

  data.frame(time_days = t, concentration = pmax(conc, 0), drug = drug_name)
}

# --- PD simulation: lipid time course ---
simulate_pd_lipids <- function(drug_name, baseline_ldl = 160, baseline_hdl = 45,
                                baseline_tg = 180, baseline_apob = 120,
                                baseline_pcsk9 = 300, duration_weeks = 52) {
  pd  <- drug_pd_ldl_reduction[[drug_name]]
  pk  <- drug_pk_params[[drug_name]]
  t_w <- seq(0, duration_weeks, by = 1)  # weekly

  # Effect onset kinetics: Emax model with delay
  if (pk$class == "Statin") {
    t50_effect <- 2   # weeks to 50% max effect
    hill       <- 2
  } else if (pk$class == "PCSK9i") {
    t50_effect <- 3
    hill       <- 2.5
  } else if (pk$class == "siRNA") {
    t50_effect <- 6
    hill       <- 1.5
  } else {
    t50_effect <- 2
    hill       <- 2
  }

  effect_frac <- (t_w^hill) / (t50_effect^hill + t_w^hill)

  ldl_c  <- baseline_ldl  * (1 + pd$ldl_pct  / 100 * effect_frac)
  hdl_c  <- baseline_hdl  * (1 + pd$hdl_pct  / 100 * effect_frac)
  tg     <- baseline_tg   * (1 + pd$tg_pct   / 100 * effect_frac)
  apob   <- baseline_apob * (1 + pd$apob_pct / 100 * effect_frac)

  # PCSK9 plasma: statins increase, PCSK9i decrease
  if (pk$class == "Statin") {
    pcsk9 <- baseline_pcsk9 * (1 + pd$pcsk9_pct / 100 * effect_frac)
  } else if (pk$class %in% c("PCSK9i", "siRNA")) {
    pcsk9 <- baseline_pcsk9 * (1 + pd$pcsk9_pct / 100 * effect_frac)
  } else {
    pcsk9 <- rep(baseline_pcsk9, length(t_w))
  }

  # LDL receptor density (relative): inversely related to PCSK9
  ldlr_density <- 100 * (baseline_pcsk9 / pmax(pcsk9, 1))^0.5

  non_hdl <- ldl_c + tg / 5  # simplified Friedewald

  data.frame(
    week        = t_w,
    LDL_C       = pmax(ldl_c, 10),
    HDL_C       = hdl_c,
    TG          = pmax(tg, 30),
    ApoB        = pmax(apob, 30),
    PCSK9       = pmax(pcsk9, 5),
    Non_HDL     = pmax(non_hdl, 15),
    LDLR_density= pmin(ldlr_density, 300),
    drug        = drug_name
  )
}

# --- CV risk reduction (Mendelian randomization / statin meta-analysis) ---
compute_cv_risk_reduction <- function(baseline_ldl_mmol, final_ldl_mmol,
                                       baseline_cv_risk_pct, years = 5) {
  delta_ldl_mmol <- baseline_ldl_mmol - final_ldl_mmol
  # 22% RRR per 1 mmol/L LDL reduction (CTT meta-analysis)
  rr_per_mmol    <- 0.78
  relative_risk  <- rr_per_mmol^delta_ldl_mmol

  # Expand to n years (linear approximation from 1-yr to 5-yr rate)
  annual_risk    <- 1 - (1 - baseline_cv_risk_pct / 100)^(1 / 10)
  baseline_5yr   <- (1 - (1 - annual_risk)^years) * 100
  treated_5yr    <- baseline_5yr * relative_risk
  arr            <- baseline_5yr - treated_5yr
  nnt            <- if (arr > 0) round(100 / arr) else Inf

  list(
    baseline_5yr  = round(baseline_5yr, 2),
    treated_5yr   = round(treated_5yr, 2),
    arr           = round(arr, 2),
    rrr           = round((1 - relative_risk) * 100, 1),
    nnt           = nnt,
    delta_ldl_mmol= round(delta_ldl_mmol, 2)
  )
}

# --- Statin myopathy risk score ---
compute_myopathy_risk <- function(age, sex, dose_intensity = "moderate",
                                   ckd = FALSE, hypothyroid = FALSE,
                                   cyp3a4_inhibitor = FALSE, prior_myopathy = FALSE) {
  score <- 0
  if (age >= 75)              score <- score + 2
  else if (age >= 65)         score <- score + 1
  if (sex == "Female")        score <- score + 1
  if (dose_intensity == "high") score <- score + 2
  else if (dose_intensity == "moderate") score <- score + 1
  if (ckd)                    score <- score + 2
  if (hypothyroid)            score <- score + 2
  if (cyp3a4_inhibitor)       score <- score + 3
  if (prior_myopathy)         score <- score + 4

  risk_label <- dplyr::case_when(
    score >= 8 ~ "Very High (>5%)",
    score >= 5 ~ "High (2-5%)",
    score >= 3 ~ "Moderate (0.5-2%)",
    TRUE       ~ "Low (<0.5%)"
  )
  list(score = score, label = risk_label)
}

# --- CYP3A4 interaction table ---
cyp3a4_interactions <- data.frame(
  Drug             = c("Cyclosporine", "Itraconazole", "Ketoconazole",
                       "Clarithromycin", "Erythromycin", "Amiodarone",
                       "Diltiazem", "Verapamil", "Nefazodone",
                       "HIV protease inhibitors", "Niacin (high dose)",
                       "Fibrates (gemfibrozil)", "Amlodipine"),
  Mechanism        = c("CYP3A4+P-gp", "CYP3A4", "CYP3A4",
                       "CYP3A4", "CYP3A4", "CYP3A4",
                       "CYP3A4", "CYP3A4", "CYP3A4",
                       "CYP3A4", "OATP1B1", "OATP1B1/CYP",
                       "CYP3A4 weak"),
  Statin_increase  = c("10-15x", "3-4x", "3-4x", "4-8x", "2-4x",
                        "2-3x", "2-4x", "2-4x", "4-6x",
                        "3-10x", "myopathy risk", "myopathy risk", "1.4x"),
  Risk_level       = c("Contraindicated", "Contraindicated", "Contraindicated",
                        "High", "High", "Moderate",
                        "Moderate", "Moderate", "High",
                        "High", "Moderate", "Moderate", "Low"),
  Action           = c("Use pravastatin/rosuvastatin only",
                        "Avoid atorvastatin, lovastatin, simvastatin",
                        "Avoid atorvastatin, lovastatin, simvastatin",
                        "Dose limit or switch",
                        "Dose limit or switch",
                        "Max simvastatin 20mg",
                        "Max simvastatin 20mg",
                        "Max simvastatin 20mg",
                        "Avoid simvastatin/lovastatin",
                        "Max simvastatin 20mg",
                        "Caution, monitor CK",
                        "Avoid rosuvastatin+gemfibrozil",
                        "Monitor"),
  stringsAsFactors = FALSE
)

# ============================================================================
# UI
# ============================================================================

ui <- navbarPage(
  title = div(
    tags$b("Dyslipidemia QSP"),
    tags$small(" | Quantitative Systems Pharmacology Simulator", style="color:#aaa;")
  ),
  theme = bs_theme(
    bootswatch  = "flatly",
    primary     = "#2c7bb6",
    secondary   = "#5cb85c",
    base_font   = font_google("Inter"),
    heading_font= font_google("Inter")
  ),
  windowTitle = "Dyslipidemia QSP Simulator",
  collapsible = TRUE,

  # ---- Header CSS ----
  header = tags$head(
    tags$style(HTML("
      .navbar-brand { font-size: 18px; }
      .well { background:#f8f9fa; border:1px solid #dee2e6; border-radius:6px; }
      .risk-box { padding:12px; border-radius:8px; text-align:center; margin:4px; font-weight:600; }
      .risk-low  { background:#d4edda; color:#155724; }
      .risk-int  { background:#fff3cd; color:#856404; }
      .risk-high { background:#f8d7da; color:#721c24; }
      .risk-vh   { background:#f5c6cb; color:#491217; border:2px solid #721c24; }
      .target-box{ padding:10px; border-radius:6px; background:#e8f4f8; border:1px solid #b8daff; }
      .metric-card { background:#fff; border:1px solid #dee2e6; border-radius:8px;
                      padding:15px; text-align:center; margin:5px; }
      .metric-card h4 { margin:0; font-size:22px; font-weight:700; }
      .metric-card p  { margin:0; color:#6c757d; font-size:12px; }
      .section-header { border-bottom:2px solid #2c7bb6; padding-bottom:6px;
                         margin-bottom:15px; color:#2c7bb6; font-weight:700; }
    "))
  ),

  # ==========================================================================
  # TAB 1: Patient Profile
  # ==========================================================================
  tabPanel(
    title = icon_text("user", "Patient Profile"),
    value = "patient",
    fluidRow(
      # --- Sidebar inputs ---
      column(3,
        wellPanel(
          h4(class="section-header", "Patient Demographics"),
          sliderInput("age", "Age (years)", 30, 85, 55, step=1),
          selectInput("sex", "Sex", c("Male","Female")),
          selectInput("race_input", "Race/Ethnicity",
                      c("White","African American","Other")),
          hr(),
          h4(class="section-header", "Lipid Panel"),
          sliderInput("ldl_base", "Baseline LDL-C (mg/dL)", 70, 300, 160, step=5),
          sliderInput("hdl_base", "Baseline HDL-C (mg/dL)", 20, 100, 45, step=1),
          sliderInput("tg_base",  "Baseline TG (mg/dL)",    50, 600, 180, step=10),
          numericInput("total_chol_in", "Total Cholesterol (mg/dL)", 220, 100, 400, step=5),
          hr(),
          h4(class="section-header", "Clinical Factors"),
          sliderInput("bmi",     "BMI (kg/m²)",  16, 50, 26, step=0.5),
          sliderInput("sbp_in",  "Systolic BP (mmHg)", 90, 200, 130, step=2),
          checkboxInput("bp_tx_in",    "On antihypertensive therapy?", FALSE),
          checkboxInput("diabetes_in", "Type 2 Diabetes?", FALSE),
          selectInput("smoking_in", "Smoking Status",
                      c("Never" = "never", "Former" = "former", "Current" = "current")),
          selectInput("risk_cat_in", "CV Risk Category (Manual Override)",
                      c("Low","Intermediate","High","Very High")),
          hr(),
          h4(class="section-header", "Comorbidities"),
          checkboxInput("chd_in",   "Known ASCVD / CHD?", FALSE),
          checkboxInput("ckd_in",   "Chronic Kidney Disease?", FALSE),
          checkboxInput("fh_in",    "Familial Hypercholesterolemia?", FALSE)
        )
      ),
      # --- Main panel ---
      column(9,
        fluidRow(
          column(6,
            h4(class="section-header", "10-Year ASCVD Risk (Pooled Cohort Equations)"),
            uiOutput("pce_risk_ui"),
            br(),
            h4(class="section-header", "LDL-C Treatment Targets"),
            uiOutput("ldl_target_ui")
          ),
          column(6,
            h4(class="section-header", "Patient Summary"),
            tableOutput("patient_summary_table"),
            br(),
            h4(class="section-header", "Lipid Classification"),
            uiOutput("lipid_class_ui")
          )
        ),
        hr(),
        fluidRow(
          column(12,
            h4(class="section-header", "Framingham / ESC Risk Factor Profile"),
            plotlyOutput("risk_factor_radar", height="320px")
          )
        )
      )
    )
  ),

  # ==========================================================================
  # TAB 2: Drug PK
  # ==========================================================================
  tabPanel(
    title = icon_text("flask", "Drug PK"),
    value = "pk",
    fluidRow(
      column(3,
        wellPanel(
          h4(class="section-header", "Drug Selection"),
          selectInput("pk_drug1", "Primary Drug",
                      choices = names(drug_pk_params),
                      selected = "Atorvastatin 40mg"),
          checkboxInput("pk_overlay", "Compare with second drug?", FALSE),
          conditionalPanel(
            "input.pk_overlay == true",
            selectInput("pk_drug2", "Second Drug",
                        choices = names(drug_pk_params),
                        selected = "Evolocumab 140mg Q2W")
          ),
          hr(),
          h4(class="section-header", "PK Parameters"),
          tableOutput("pk_param_table"),
          hr(),
          h4(class="section-header", "Simulation Settings"),
          sliderInput("pk_weeks", "Duration (weeks)", 2, 26, 14, step=1),
          checkboxInput("pk_log_scale", "Log-scale Y-axis?", FALSE),
          checkboxInput("pk_show_ss",   "Show SS concentration band?", TRUE)
        )
      ),
      column(9,
        fluidRow(
          column(12,
            h4(class="section-header", "Plasma Concentration-Time Profile"),
            plotlyOutput("pk_conc_plot", height = "380px")
          )
        ),
        fluidRow(
          column(6,
            h4(class="section-header", "Hepatic Exposure (Statin Active Metabolite)"),
            uiOutput("hepatic_exposure_ui"),
            plotlyOutput("pk_hepatic_plot", height = "260px")
          ),
          column(6,
            h4(class="section-header", "Drug Class PK Comparison"),
            plotlyOutput("pk_class_compare", height = "260px")
          )
        )
      )
    )
  ),

  # ==========================================================================
  # TAB 3: PD – Lipid Biomarkers
  # ==========================================================================
  tabPanel(
    title = icon_text("chart-line", "PD - Lipid Biomarkers"),
    value = "pd",
    fluidRow(
      column(3,
        wellPanel(
          h4(class="section-header", "Treatment"),
          selectInput("pd_drug", "Select Drug",
                      choices = names(drug_pd_ldl_reduction),
                      selected = "Atorvastatin 40mg"),
          hr(),
          h4(class="section-header", "Baseline Lipids"),
          p("(From Patient Profile tab)"),
          sliderInput("pd_duration", "Simulation Duration (weeks)", 12, 52, 52, step=4),
          hr(),
          h4(class="section-header", "Display Options"),
          checkboxGroupInput("pd_biomarkers", "Biomarkers to Show",
            choices = c("LDL-C"="LDL_C","HDL-C"="HDL_C","TG"="TG",
                        "ApoB"="ApoB","Non-HDL-C"="Non_HDL",
                        "PCSK9"="PCSK9","LDLR Density"="LDLR_density"),
            selected = c("LDL_C","HDL_C","TG","ApoB")),
          checkboxInput("pd_pct_change", "Show % Change from Baseline?", FALSE),
          checkboxInput("pd_show_targets", "Show ACC/AHA/ESC Target Lines?", TRUE)
        )
      ),
      column(9,
        fluidRow(
          column(12,
            h4(class="section-header", "Lipid Biomarker Time Course"),
            plotlyOutput("pd_lipid_plot", height = "380px")
          )
        ),
        fluidRow(
          column(6,
            h4(class="section-header", "PCSK9 Plasma Levels"),
            plotlyOutput("pd_pcsk9_plot", height = "260px")
          ),
          column(6,
            h4(class="section-header", "LDL Receptor Density (Relative %)"),
            plotlyOutput("pd_ldlr_plot", height = "260px")
          )
        )
      )
    )
  ),

  # ==========================================================================
  # TAB 4: Clinical Endpoints
  # ==========================================================================
  tabPanel(
    title = icon_text("heartbeat", "Clinical Endpoints"),
    value = "endpoints",
    fluidRow(
      column(3,
        wellPanel(
          h4(class="section-header", "Treatment Parameters"),
          selectInput("ep_drug", "Select Drug",
                      choices = names(drug_pd_ldl_reduction),
                      selected = "Atorvastatin 40mg"),
          sliderInput("ep_adherence", "Treatment Adherence (%)", 50, 100, 85, step=5),
          hr(),
          h4(class="section-header", "Projection Settings"),
          sliderInput("ep_years", "Projection Horizon (years)", 1, 10, 5, step=1),
          selectInput("ep_guideline", "Guideline Framework",
                      c("ACC/AHA 2019","ESC 2019","Both")),
          hr(),
          h4(class="section-header", "Plaque Model"),
          sliderInput("ep_plaque_base", "Baseline Plaque Volume (mm³)", 20, 200, 80, step=5),
          sliderInput("ep_progression", "Annual Plaque Progression Rate (%)", 0, 10, 3, step=0.5)
        )
      ),
      column(9,
        fluidRow(
          column(4,
            h4(class="section-header", "CV Risk Metrics"),
            uiOutput("cv_metrics_ui")
          ),
          column(8,
            h4(class="section-header", "Cumulative MACE Risk Over Time"),
            plotlyOutput("mace_plot", height = "320px")
          )
        ),
        fluidRow(
          column(6,
            h4(class="section-header", "Atherosclerotic Plaque Simulation"),
            plotlyOutput("plaque_plot", height = "280px")
          ),
          column(6,
            h4(class="section-header", "NNT / Absolute Risk Reduction"),
            plotlyOutput("nnt_plot", height = "280px"),
            uiOutput("nnt_summary_ui")
          )
        )
      )
    )
  ),

  # ==========================================================================
  # TAB 5: Scenario Comparison
  # ==========================================================================
  tabPanel(
    title = icon_text("balance-scale", "Scenario Comparison"),
    value = "scenarios",
    fluidRow(
      column(3,
        wellPanel(
          h4(class="section-header", "Select Scenarios"),
          p("Choose up to 6 treatment regimens to compare:"),
          lapply(1:6, function(i) {
            default_choices <- c(
              "No treatment",
              "Atorvastatin 10mg",
              "Atorvastatin 40mg",
              "Atorvastatin 80mg",
              "Rosuvastatin 40mg",
              "Evolocumab 140mg Q2W"
            )
            selectInput(
              paste0("sc_drug_", i),
              paste("Scenario", i),
              choices = c("None", "No treatment", names(drug_pd_ldl_reduction)),
              selected = if (i <= 6) default_choices[i] else "None"
            )
          }),
          hr(),
          actionButton("run_scenarios", "Run Comparison",
                       class="btn btn-primary btn-block",
                       icon = icon("play"))
        )
      ),
      column(9,
        fluidRow(
          column(12,
            h4(class="section-header", "LDL-C at Weeks 12, 26 & 52 by Scenario"),
            plotlyOutput("sc_bar_plot", height = "340px")
          )
        ),
        fluidRow(
          column(12,
            h4(class="section-header", "Scenario Comparison Table"),
            DTOutput("sc_table"),
            br(),
            h4(class="section-header", "Multi-Parameter Spider Plot"),
            plotlyOutput("sc_spider_plot", height = "380px")
          )
        )
      )
    )
  ),

  # ==========================================================================
  # TAB 6: Biomarker Dashboard
  # ==========================================================================
  tabPanel(
    title = icon_text("dashboard", "Biomarker Dashboard"),
    value = "dashboard",
    fluidRow(
      column(3,
        wellPanel(
          h4(class="section-header", "PCSK9 Genetics"),
          sliderInput("pcsk9_base",     "Baseline PCSK9 (ng/mL)", 100, 800, 300, step=10),
          selectInput("pcsk9_variant", "Genetic Variant",
                      c("Wild-type (normal)"     = "wt",
                        "GOF variant (elevated)" = "gof",
                        "LOF variant (reduced)"  = "lof")),
          hr(),
          h4(class="section-header", "Treatment"),
          selectInput("db_drug",      "Drug",
                      choices = names(drug_pd_ldl_reduction),
                      selected = "Atorvastatin 40mg"),
          hr(),
          h4(class="section-header", "Myopathy Risk Assessment"),
          checkboxInput("db_hypothyroid",  "Hypothyroidism?", FALSE),
          checkboxInput("db_cyp_inhibitor","On CYP3A4 inhibitor?", FALSE),
          checkboxInput("db_prior_myo",    "Prior statin myopathy?", FALSE),
          hr(),
          h4(class="section-header", "DDI Checker"),
          selectInput("ddi_drug", "Co-medication",
                      choices = c("None", cyp3a4_interactions$Drug))
        )
      ),
      column(9,
        fluidRow(
          column(6,
            h4(class="section-header", "Lipid Panel Trends with Guideline Targets"),
            plotlyOutput("db_lipid_trends", height = "320px")
          ),
          column(6,
            h4(class="section-header", "PCSK9 & LDLR Response"),
            plotlyOutput("db_pcsk9_ldlr", height = "320px")
          )
        ),
        fluidRow(
          column(6,
            h4(class="section-header", "Statin Myopathy Risk Score"),
            uiOutput("myopathy_risk_ui"),
            plotlyOutput("myopathy_gauge", height = "220px")
          ),
          column(6,
            h4(class="section-header", "Drug-Drug Interaction (CYP3A4)"),
            uiOutput("ddi_result_ui"),
            DTOutput("ddi_table")
          )
        )
      )
    )
  )
)

# Helper function to create icon + text label
icon_text <- function(icon_name, text) {
  tagList(icon(icon_name), text)
}

# ============================================================================
# SERVER
# ============================================================================

server <- function(input, output, session) {

  # --------------------------------------------------------------------------
  # REACTIVE: PCE risk calculation
  # --------------------------------------------------------------------------
  pce_risk_val <- reactive({
    smoker <- input$smoking_in == "current"
    compute_pce_risk(
      age        = input$age,
      sex        = input$sex,
      total_chol = input$total_chol_in,
      hdl_c      = input$hdl_base,
      sbp        = input$sbp_in,
      bp_treatment = input$bp_tx_in,
      diabetes   = input$diabetes_in,
      smoker     = smoker
    )
  })

  # Override risk category if CHD or FH
  effective_risk_cat <- reactive({
    if (input$chd_in || input$fh_in) return("Very High")
    risk <- pce_risk_val()
    if (!is.null(input$risk_cat_in)) return(input$risk_cat_in)
    if (risk >= 20)      "Very High"
    else if (risk >= 10) "High"
    else if (risk >= 5)  "Intermediate"
    else                  "Low"
  })

  # --------------------------------------------------------------------------
  # TAB 1 OUTPUTS
  # --------------------------------------------------------------------------
  output$pce_risk_ui <- renderUI({
    risk <- round(pce_risk_val(), 1)
    cat_label <- effective_risk_cat()
    css_class <- switch(cat_label,
      "Very High"    = "risk-vh",
      "High"         = "risk-high",
      "Intermediate" = "risk-int",
      "Low"          = "risk-low"
    )
    div(
      div(class = paste("risk-box", css_class),
          h2(paste0(risk, "%"), style="margin:0;font-size:48px;"),
          p("10-Year ASCVD Risk")
      ),
      br(),
      div(class = paste("risk-box", css_class),
          h4(paste("Category:", cat_label))
      ),
      if (risk >= 20) tags$p(tags$em("⚠ Very high risk: statin + possible PCSK9 inhibitor indicated"), style="color:#721c24;") else NULL
    )
  })

  output$ldl_target_ui <- renderUI({
    targets <- get_ldl_target(effective_risk_cat())
    div(
      class = "target-box",
      h4("LDL-C Targets"),
      tags$table(class="table table-sm",
        tags$thead(tags$tr(tags$th("Guideline"), tags$th("Target (mg/dL)"), tags$th("Current LDL-C"))),
        tags$tbody(
          tags$tr(
            tags$td("ESC 2019"),
            tags$td(paste0("<", targets$esc)),
            tags$td(
              style = if (input$ldl_base > targets$esc) "color:red;font-weight:700;" else "color:green;font-weight:700;",
              paste(input$ldl_base, if (input$ldl_base > targets$esc) "▲ Above target" else "✓ At target")
            )
          ),
          tags$tr(
            tags$td("ACC/AHA 2019"),
            tags$td(paste0("<", targets$acc)),
            tags$td(
              style = if (input$ldl_base > targets$acc) "color:red;font-weight:700;" else "color:green;font-weight:700;",
              paste(input$ldl_base, if (input$ldl_base > targets$acc) "▲ Above target" else "✓ At target")
            )
          )
        )
      ),
      p(class="text-muted", paste("ESC 2019 extreme-risk target: <40 mg/dL (recurrent events on max therapy)"))
    )
  })

  output$patient_summary_table <- renderTable({
    data.frame(
      Parameter = c("Age","Sex","BMI","LDL-C","HDL-C","Triglycerides","Total Chol",
                    "Systolic BP","Diabetes","Smoking","Known ASCVD","CKD","FH"),
      Value = c(
        paste(input$age, "years"),
        input$sex,
        paste0(input$bmi, " kg/m²"),
        paste0(input$ldl_base, " mg/dL"),
        paste0(input$hdl_base, " mg/dL"),
        paste0(input$tg_base, " mg/dL"),
        paste0(input$total_chol_in, " mg/dL"),
        paste0(input$sbp_in, " mmHg"),
        if (input$diabetes_in) "Yes" else "No",
        tools::toTitleCase(input$smoking_in),
        if (input$chd_in) "Yes" else "No",
        if (input$ckd_in) "Yes" else "No",
        if (input$fh_in) "Yes" else "No"
      ),
      stringsAsFactors = FALSE
    )
  }, striped=TRUE, hover=TRUE, bordered=TRUE)

  output$lipid_class_ui <- renderUI({
    ldl <- input$ldl_base
    hdl <- input$hdl_base
    tg  <- input$tg_base
    div(
      tags$ul(
        tags$li(
          style = if (ldl < 100) "color:green;" else if (ldl < 130) "color:orange;" else "color:red;",
          paste0("LDL-C: ", ldl, " mg/dL — ",
                 if (ldl < 100) "Optimal" else if (ldl < 130) "Near-optimal" else if (ldl < 160) "Borderline high" else "High")
        ),
        tags$li(
          style = if (hdl >= 60) "color:green;" else if (hdl >= 40) "color:orange;" else "color:red;",
          paste0("HDL-C: ", hdl, " mg/dL — ",
                 if (hdl >= 60) "High (protective)" else if (hdl >= 40) "Acceptable" else "Low (risk factor)")
        ),
        tags$li(
          style = if (tg < 150) "color:green;" else if (tg < 200) "color:orange;" else "color:red;",
          paste0("TG: ", tg, " mg/dL — ",
                 if (tg < 150) "Normal" else if (tg < 200) "Borderline" else if (tg < 500) "High" else "Very High")
        )
      )
    )
  })

  output$risk_factor_radar <- renderPlotly({
    # Normalize risk factors to 0-100 scale for radar
    age_norm    <- (input$age - 30) / (85 - 30) * 100
    ldl_norm    <- (input$ldl_base - 70) / (300 - 70) * 100
    hdl_norm    <- (1 - (input$hdl_base - 20) / (100 - 20)) * 100  # inverse
    tg_norm     <- (input$tg_base - 50) / (500 - 50) * 100
    bmi_norm    <- (input$bmi - 16) / (50 - 16) * 100
    sbp_norm    <- (input$sbp_in - 90) / (200 - 90) * 100
    dm_norm     <- if (input$diabetes_in) 80 else 10
    smk_norm    <- if (input$smoking_in == "current") 90 else if (input$smoking_in == "former") 40 else 5

    categories  <- c("Age","LDL-C","Low HDL","Triglycerides","BMI","SBP","Diabetes","Smoking", "Age")
    values      <- c(age_norm, ldl_norm, hdl_norm, tg_norm, bmi_norm, sbp_norm, dm_norm, smk_norm, age_norm)

    plot_ly(
      type = 'scatterpolar',
      r    = values,
      theta= categories,
      fill = 'toself',
      fillcolor = 'rgba(44,123,182,0.3)',
      line = list(color='#2c7bb6', width=2),
      name = "Patient"
    ) %>%
      layout(
        polar  = list(radialaxis = list(visible=TRUE, range=c(0,100))),
        margin = list(l=50,r=50,t=30,b=30)
      ) %>%
      config(displayModeBar=FALSE)
  })

  # --------------------------------------------------------------------------
  # TAB 2: PK OUTPUTS
  # --------------------------------------------------------------------------
  output$pk_param_table <- renderTable({
    pk <- drug_pk_params[[input$pk_drug1]]
    data.frame(
      Parameter = c("Drug Class","Dose","Cmax (ng/mL or μg/mL)","AUC",
                    "t½ (h)","Tmax (h)","Bioavailability (%)"),
      Value = c(
        pk$class, paste0(pk$dose, " mg"),
        pk$Cmax, pk$AUC,
        pk$t_half, pk$Tmax,
        if (!is.na(pk$F)) paste0(pk$F, "%") else "N/A (SC)"
      ),
      stringsAsFactors = FALSE
    )
  }, striped=TRUE, hover=TRUE)

  output$pk_conc_plot <- renderPlotly({
    df1 <- simulate_pk(input$pk_drug1, duration_weeks = input$pk_weeks)

    p <- plot_ly() %>%
      add_lines(data=df1, x=~time_days, y=~concentration,
                name = input$pk_drug1,
                line = list(color='#2c7bb6', width=2.5))

    if (input$pk_overlay && !is.null(input$pk_drug2) && input$pk_drug2 != input$pk_drug1) {
      df2 <- simulate_pk(input$pk_drug2, duration_weeks = input$pk_weeks)
      p <- p %>%
        add_lines(data=df2, x=~time_days, y=~concentration,
                  name = input$pk_drug2,
                  line = list(color='#d62728', width=2.5, dash='dash'))
    }

    y_type <- if (input$pk_log_scale) "log" else "linear"

    p %>% layout(
      xaxis = list(title = "Time (days)"),
      yaxis = list(title = "Plasma Concentration (ng/mL)", type = y_type),
      legend = list(orientation = "h"),
      hovermode = "x unified",
      margin = list(l=60,r=20,t=20,b=50)
    ) %>% config(displayModeBar=FALSE)
  })

  output$hepatic_exposure_ui <- renderUI({
    pk <- drug_pk_params[[input$pk_drug1]]
    if (pk$class != "Statin") {
      return(p(class="text-muted",
               "Hepatic active metabolite data shown for statins only. Selected drug is a ",
               strong(pk$class), "."))
    }
    # Hepatic extraction ratio: ~60-70% for most statins
    her     <- if (grepl("Atorvastatin", input$pk_drug1)) 0.65 else 0.63
    c_liver <- pk$Cmax * her * 8  # rough portal vein concentration factor
    div(
      class = "metric-card",
      p("Hepatic Active Metabolite Cmax"),
      h4(paste0(round(c_liver, 1), " nM")),
      p(class="text-muted", "Estimated portal concentration"),
      p(class="text-muted", paste0("HMG-CoA reductase IC50: ~1-10 nM | Ratio: ", round(c_liver/5, 0), "x IC50"))
    )
  })

  output$pk_hepatic_plot <- renderPlotly({
    pk <- drug_pk_params[[input$pk_drug1]]
    if (pk$class != "Statin") {
      return(plot_ly() %>% layout(title="Hepatic data: statins only"))
    }
    df <- simulate_pk(input$pk_drug1, duration_weeks = min(input$pk_weeks, 7))
    her <- 0.65
    df$hepatic <- df$concentration * her * 8

    plot_ly(df, x=~time_days, y=~hepatic, type='scatter', mode='lines',
            line=list(color='#2ca02c', width=2),
            name="Hepatic [Active metabolite]") %>%
      add_lines(x=range(df$time_days), y=c(5,5), line=list(color='red',dash='dot'),
                name="IC50 ~5 nM") %>%
      layout(
        xaxis = list(title="Time (days)"),
        yaxis = list(title="Concentration (nM)"),
        margin = list(l=60,r=20,t=20,b=50)
      ) %>% config(displayModeBar=FALSE)
  })

  output$pk_class_compare <- renderPlotly({
    comparison_drugs <- c("Atorvastatin 40mg","Rosuvastatin 20mg",
                          "Evolocumab 140mg Q2W","Inclisiran 284mg SC","Ezetimibe 10mg")
    params <- lapply(comparison_drugs, function(d) drug_pk_params[[d]])
    df_compare <- data.frame(
      Drug     = comparison_drugs,
      Class    = sapply(params, `[[`, "class"),
      t_half_d = sapply(params, function(x) x$t_half / 24),
      AUC      = sapply(params, `[[`, "AUC"),
      stringsAsFactors = FALSE
    )
    plot_ly(df_compare, x=~Drug, y=~t_half_d, type='bar',
            color=~Class,
            colors=c("Statin"="#2c7bb6","PCSK9i"="#d62728",
                     "siRNA"="#2ca02c","NPC1L1i"="#ff7f0e"),
            text=~paste0(round(t_half_d,1)," days"),
            textposition="outside") %>%
      layout(
        xaxis = list(title="", tickangle=-20),
        yaxis = list(title="Half-life (days)"),
        margin = list(l=60,r=20,t=20,b=80)
      ) %>% config(displayModeBar=FALSE)
  })

  # --------------------------------------------------------------------------
  # TAB 3: PD OUTPUTS
  # --------------------------------------------------------------------------
  pd_data <- reactive({
    simulate_pd_lipids(
      drug_name     = input$pd_drug,
      baseline_ldl  = input$ldl_base,
      baseline_hdl  = input$hdl_base,
      baseline_tg   = input$tg_base,
      baseline_apob = round(input$ldl_base * 0.75),  # rough ApoB estimate
      baseline_pcsk9= 300,
      duration_weeks= input$pd_duration
    )
  })

  output$pd_lipid_plot <- renderPlotly({
    df  <- pd_data()
    sel <- input$pd_biomarkers

    color_map <- c(LDL_C="#d62728", HDL_C="#2ca02c", TG="#ff7f0e",
                   ApoB="#9467bd", Non_HDL="#8c564b", PCSK9="#1f77b4",
                   LDLR_density="#17becf")
    label_map <- c(LDL_C="LDL-C", HDL_C="HDL-C", TG="Triglycerides",
                   ApoB="ApoB", Non_HDL="Non-HDL-C", PCSK9="PCSK9 (ng/mL)",
                   LDLR_density="LDLR Density (%)")

    sel_valid <- intersect(sel, c("LDL_C","HDL_C","TG","ApoB","Non_HDL"))
    if (length(sel_valid) == 0) sel_valid <- "LDL_C"

    p <- plot_ly()
    for (bm in sel_valid) {
      y_vals <- if (input$pd_pct_change) {
        (df[[bm]] - df[[bm]][1]) / df[[bm]][1] * 100
      } else {
        df[[bm]]
      }
      p <- p %>% add_lines(x=df$week, y=y_vals, name=label_map[bm],
                           line=list(color=color_map[bm], width=2.5))
    }

    # Target lines
    if (input$pd_show_targets && !input$pd_pct_change && "LDL_C" %in% sel_valid) {
      targets <- get_ldl_target(effective_risk_cat())
      p <- p %>%
        add_lines(x=c(0, input$pd_duration), y=c(targets$esc, targets$esc),
                  line=list(color='blue', dash='dot', width=1),
                  name=paste0("ESC target <", targets$esc)) %>%
        add_lines(x=c(0, input$pd_duration), y=c(targets$acc, targets$acc),
                  line=list(color='darkgreen', dash='dash', width=1),
                  name=paste0("ACC/AHA target <", targets$acc))
    }

    y_title <- if (input$pd_pct_change) "% Change from Baseline" else "Concentration (mg/dL)"
    p %>% layout(
      xaxis  = list(title="Week"),
      yaxis  = list(title=y_title),
      legend = list(orientation="h", y=-0.25),
      hovermode="x unified",
      margin = list(l=60,r=20,t=20,b=80)
    ) %>% config(displayModeBar=FALSE)
  })

  output$pd_pcsk9_plot <- renderPlotly({
    df <- pd_data()
    pk <- drug_pk_params[[input$pd_drug]]

    color <- if (pk$class == "Statin") "#d62728" else if (pk$class %in% c("PCSK9i","siRNA")) "#2c7bb6" else "#7f7f7f"
    note  <- if (pk$class == "Statin") "Statins ↑ PCSK9 (reactive upregulation)" else
             if (pk$class %in% c("PCSK9i","siRNA")) "PCSK9 inhibition/silencing" else
             "No direct PCSK9 effect"

    plot_ly(df, x=~week, y=~PCSK9, type='scatter', mode='lines',
            line=list(color=color, width=2.5),
            name="PCSK9") %>%
      add_lines(x=c(0, max(df$week)), y=c(300, 300),
                line=list(color='grey', dash='dot'), name="Normal range") %>%
      layout(
        title  = list(text=note, font=list(size=11)),
        xaxis  = list(title="Week"),
        yaxis  = list(title="PCSK9 (ng/mL)"),
        margin = list(l=60,r=20,t=40,b=50)
      ) %>% config(displayModeBar=FALSE)
  })

  output$pd_ldlr_plot <- renderPlotly({
    df <- pd_data()
    plot_ly(df, x=~week, y=~LDLR_density, type='scatter', mode='lines',
            fill='tozeroy', fillcolor='rgba(44,123,182,0.2)',
            line=list(color='#2c7bb6', width=2.5),
            name="LDLR Density") %>%
      add_lines(x=c(0, max(df$week)), y=c(100, 100),
                line=list(color='grey', dash='dot'), name="Baseline") %>%
      layout(
        xaxis  = list(title="Week"),
        yaxis  = list(title="LDLR Density (% of baseline)"),
        margin = list(l=60,r=20,t=20,b=50)
      ) %>% config(displayModeBar=FALSE)
  })

  # --------------------------------------------------------------------------
  # TAB 4: CLINICAL ENDPOINTS
  # --------------------------------------------------------------------------
  ep_ldl_final <- reactive({
    pd  <- drug_pd_ldl_reduction[[input$ep_drug]]
    adh <- input$ep_adherence / 100
    # Effective reduction adjusted for adherence
    eff_red <- pd$ldl_pct * adh
    baseline_ldl <- input$ldl_base
    baseline_ldl * (1 + eff_red / 100)
  })

  cv_risk_results <- reactive({
    final_ldl <- ep_ldl_final()
    base_mmol  <- input$ldl_base / 38.67
    final_mmol <- final_ldl / 38.67
    base_risk  <- pce_risk_val()
    compute_cv_risk_reduction(base_mmol, final_mmol, base_risk, years = input$ep_years)
  })

  output$cv_metrics_ui <- renderUI({
    res <- cv_risk_results()
    final_ldl <- ep_ldl_final()
    div(
      div(class="metric-card",
          p("Baseline LDL-C"), h4(paste0(input$ldl_base, " mg/dL"))),
      div(class="metric-card",
          p("Treated LDL-C"), h4(paste0(round(final_ldl,0), " mg/dL"),
                                 style="color:#2ca02c;")),
      div(class="metric-card",
          p("LDL-C Reduction"), h4(paste0(round(res$delta_ldl_mmol * 38.67, 0), " mg/dL"),
                                   style="color:#d62728;")),
      div(class="metric-card",
          p(paste0(input$ep_years,"-yr CV Risk (baseline)")),
          h4(paste0(res$baseline_5yr, "%"), style="color:#d62728;")),
      div(class="metric-card",
          p(paste0(input$ep_years,"-yr CV Risk (treated)")),
          h4(paste0(res$treated_5yr, "%"), style="color:#2ca02c;")),
      div(class="metric-card",
          p("Relative RR Reduction"), h4(paste0(res$rrr, "%"))),
      div(class="metric-card",
          p("Absolute RR Reduction"), h4(paste0(res$arr, "%"))),
      div(class="metric-card",
          p("NNT"), h4(res$nnt))
    )
  })

  output$mace_plot <- renderPlotly({
    res   <- cv_risk_results()
    years <- seq(0, input$ep_years, by=0.25)

    annual_risk_base <- 1 - (1 - pce_risk_val() / 100)^(1/10)
    annual_risk_tx   <- annual_risk_base * (res$treated_5yr / res$baseline_5yr)

    cum_base <- (1 - (1 - annual_risk_base)^years) * 100
    cum_tx   <- (1 - (1 - annual_risk_tx)^years)   * 100

    plot_ly() %>%
      add_lines(x=years, y=cum_base, name="No treatment",
                line=list(color='#d62728', width=2.5)) %>%
      add_lines(x=years, y=cum_tx, name=input$ep_drug,
                line=list(color='#2ca02c', width=2.5)) %>%
      add_ribbons(x=c(years, rev(years)),
                  y=c(cum_base, rev(cum_tx)),
                  fillcolor='rgba(214,39,40,0.15)',
                  line=list(color='transparent'),
                  name="ARR area",
                  showlegend=FALSE) %>%
      layout(
        xaxis  = list(title="Years"),
        yaxis  = list(title="Cumulative MACE Risk (%)"),
        legend = list(orientation="h"),
        hovermode="x unified",
        margin = list(l=60,r=20,t=20,b=50)
      ) %>% config(displayModeBar=FALSE)
  })

  output$plaque_plot <- renderPlotly({
    yrs    <- seq(0, max(input$ep_years, 5), by=0.5)
    p_base <- input$ep_plaque_base
    prog   <- input$ep_progression / 100

    pd  <- drug_pd_ldl_reduction[[input$ep_drug]]
    # Plaque regression rate proportional to LDL-C lowering
    regress_factor <- abs(pd$ldl_pct) / 100 * 0.4  # max 40% regression possible

    plaque_untreated <- p_base * (1 + prog)^yrs
    plaque_treated   <- p_base * (1 + prog * (1 - regress_factor * 2))^yrs

    plot_ly() %>%
      add_lines(x=yrs, y=plaque_untreated, name="Untreated",
                line=list(color='#d62728', width=2)) %>%
      add_lines(x=yrs, y=plaque_treated, name=input$ep_drug,
                line=list(color='#2c7bb6', width=2)) %>%
      add_lines(x=c(0, max(yrs)), y=c(p_base, p_base),
                line=list(color='grey', dash='dot'), name="Baseline volume") %>%
      layout(
        xaxis  = list(title="Years"),
        yaxis  = list(title="Plaque Volume (mm³)"),
        margin = list(l=60,r=20,t=20,b=50)
      ) %>% config(displayModeBar=FALSE)
  })

  output$nnt_plot <- renderPlotly({
    # NNT vs adherence
    adh_range <- seq(50, 100, by=5)
    pd  <- drug_pd_ldl_reduction[[input$ep_drug]]
    base_mmol <- input$ldl_base / 38.67
    base_risk <- pce_risk_val()

    nnt_vals <- sapply(adh_range, function(a) {
      eff_red <- pd$ldl_pct * a / 100
      final_ldl <- input$ldl_base * (1 + eff_red / 100)
      final_mmol <- final_ldl / 38.67
      res <- compute_cv_risk_reduction(base_mmol, final_mmol, base_risk, input$ep_years)
      min(res$nnt, 500)
    })

    plot_ly(x=adh_range, y=nnt_vals, type='scatter', mode='lines+markers',
            line=list(color='#2c7bb6', width=2),
            marker=list(color='#2c7bb6', size=6),
            text=paste0("NNT=", round(nnt_vals)),
            hoverinfo="text+x") %>%
      add_lines(x=c(input$ep_adherence, input$ep_adherence), y=c(0, max(nnt_vals)),
                line=list(color='red', dash='dot'), name="Current adherence") %>%
      layout(
        xaxis  = list(title="Adherence (%)"),
        yaxis  = list(title=paste0("NNT (", input$ep_years, "-yr MACE)")),
        showlegend=FALSE,
        margin = list(l=60,r=20,t=20,b=50)
      ) %>% config(displayModeBar=FALSE)
  })

  output$nnt_summary_ui <- renderUI({
    res <- cv_risk_results()
    div(class="target-box",
        tags$b("Summary: "),
        paste0("Treating ", res$nnt, " patients for ", input$ep_years,
               " years prevents 1 MACE event. ARR = ", res$arr,
               "%, RRR = ", res$rrr, "%.")
    )
  })

  # --------------------------------------------------------------------------
  # TAB 5: SCENARIO COMPARISON
  # --------------------------------------------------------------------------
  sc_drugs_selected <- eventReactive(input$run_scenarios, {
    raw <- sapply(1:6, function(i) input[[paste0("sc_drug_", i)]])
    raw[raw != "None"]
  }, ignoreNULL = FALSE)

  sc_results <- reactive({
    drugs <- sc_drugs_selected()
    if (length(drugs) == 0) return(NULL)

    lapply(drugs, function(d) {
      if (d == "No treatment") {
        data.frame(
          Drug      = "No treatment",
          LDL_12w   = input$ldl_base,
          LDL_26w   = input$ldl_base,
          LDL_52w   = input$ldl_base,
          LDL_pct   = 0,
          HDL_pct   = 0,
          TG_pct    = 0,
          ApoB_pct  = 0,
          stringsAsFactors = FALSE
        )
      } else {
        pd  <- drug_pd_ldl_reduction[[d]]
        sim <- simulate_pd_lipids(d, input$ldl_base, input$hdl_base,
                                   input$tg_base,
                                   round(input$ldl_base * 0.75), 300, 52)
        targets <- get_ldl_target(effective_risk_cat())
        data.frame(
          Drug      = d,
          LDL_12w   = round(sim$LDL_C[sim$week == 12], 1),
          LDL_26w   = round(sim$LDL_C[sim$week == 26], 1),
          LDL_52w   = round(sim$LDL_C[sim$week == 52], 1),
          LDL_pct   = pd$ldl_pct,
          HDL_pct   = pd$hdl_pct,
          TG_pct    = pd$tg_pct,
          ApoB_pct  = pd$apob_pct,
          stringsAsFactors = FALSE
        )
      }
    }) %>% bind_rows()
  })

  output$sc_table <- renderDT({
    df <- sc_results()
    if (is.null(df)) return(NULL)
    targets <- get_ldl_target(effective_risk_cat())
    df_display <- df %>%
      mutate(
        `LDL-C 12w (mg/dL)` = LDL_12w,
        `LDL-C 26w (mg/dL)` = LDL_26w,
        `LDL-C 52w (mg/dL)` = LDL_52w,
        `LDL % reduction`   = paste0(LDL_pct, "%"),
        `HDL % change`      = paste0("+", HDL_pct, "%"),
        `TG % change`       = paste0(TG_pct, "%"),
        `ApoB % change`     = paste0(ApoB_pct, "%"),
        `ESC Goal (<55)`    = ifelse(LDL_52w < targets$esc, "YES", "NO"),
        `ACC Goal (<70)`    = ifelse(LDL_52w < targets$acc, "YES", "NO")
      ) %>%
      select(Drug, `LDL-C 12w (mg/dL)`, `LDL-C 26w (mg/dL)`, `LDL-C 52w (mg/dL)`,
             `LDL % reduction`, `HDL % change`, `TG % change`,
             `ESC Goal (<55)`, `ACC Goal (<70)`)

    datatable(df_display, options=list(pageLength=10, dom='t'),
              rownames=FALSE) %>%
      formatStyle("ESC Goal (<55)", backgroundColor=styleEqual(c("YES","NO"), c("#d4edda","#f8d7da"))) %>%
      formatStyle("ACC Goal (<70)", backgroundColor=styleEqual(c("YES","NO"), c("#d4edda","#f8d7da")))
  })

  output$sc_bar_plot <- renderPlotly({
    df <- sc_results()
    if (is.null(df)) return(plot_ly() %>% layout(title="Click 'Run Comparison'"))
    targets <- get_ldl_target(effective_risk_cat())

    df_long <- df %>%
      select(Drug, LDL_12w, LDL_26w, LDL_52w) %>%
      pivot_longer(cols=-Drug, names_to="Timepoint", values_to="LDL_C") %>%
      mutate(Timepoint = recode(Timepoint, LDL_12w="Week 12", LDL_26w="Week 26", LDL_52w="Week 52"))

    colors <- c("#2c7bb6","#d62728","#2ca02c","#ff7f0e","#9467bd","#8c564b")
    drugs  <- unique(df$Drug)

    p <- plot_ly()
    for (i in seq_along(drugs)) {
      sub <- df_long %>% filter(Drug == drugs[i])
      p <- p %>% add_bars(data=sub, x=~Timepoint, y=~LDL_C,
                          name=drugs[i],
                          marker=list(color=colors[i %% length(colors) + 1]))
    }

    p %>%
      add_lines(x=c("Week 12","Week 26","Week 52"),
                y=rep(targets$esc, 3),
                line=list(color='blue', dash='dot', width=1.5),
                name=paste0("ESC <", targets$esc), inherit=FALSE) %>%
      add_lines(x=c("Week 12","Week 26","Week 52"),
                y=rep(targets$acc, 3),
                line=list(color='darkgreen', dash='dash', width=1.5),
                name=paste0("ACC <", targets$acc), inherit=FALSE) %>%
      layout(
        barmode = "group",
        xaxis   = list(title=""),
        yaxis   = list(title="LDL-C (mg/dL)"),
        legend  = list(orientation="h", y=-0.3),
        margin  = list(l=60,r=20,t=20,b=80)
      ) %>% config(displayModeBar=FALSE)
  })

  output$sc_spider_plot <- renderPlotly({
    df <- sc_results()
    if (is.null(df)) return(plot_ly() %>% layout(title="Click 'Run Comparison'"))

    categories <- c("LDL ↓","HDL ↑","TG ↓","ApoB ↓","LDL Goal", "LDL ↓")
    targets <- get_ldl_target(effective_risk_cat())
    colors  <- c("#2c7bb6","#d62728","#2ca02c","#ff7f0e","#9467bd","#8c564b")

    p <- plot_ly(type='scatterpolar')
    for (i in seq_len(nrow(df))) {
      row <- df[i, ]
      # Normalize each metric to 0-100 scale
      ldl_norm  <- pmin(abs(row$LDL_pct), 70) / 70 * 100
      hdl_norm  <- pmin(row$HDL_pct, 15) / 15 * 100
      tg_norm   <- pmin(abs(row$TG_pct), 30) / 30 * 100
      apob_norm <- pmin(abs(row$ApoB_pct), 60) / 60 * 100
      goal_norm <- max(0, (targets$acc - row$LDL_52w) / targets$acc * 100)
      vals <- c(ldl_norm, hdl_norm, tg_norm, apob_norm, goal_norm, ldl_norm)

      p <- p %>% add_trace(
        r     = vals,
        theta = categories,
        fill  = 'toself',
        fillcolor = paste0(gsub("\\)", ",0.15)", gsub("rgb", "rgba",
                     plotly::toRGB(colors[i %% length(colors) + 1]))),
                           .na.rm = TRUE),
        line  = list(color = colors[i %% length(colors) + 1], width=2),
        name  = row$Drug
      )
    }
    p %>% layout(
      polar  = list(radialaxis=list(visible=TRUE, range=c(0,100))),
      legend = list(orientation="h"),
      margin = list(l=50,r=50,t=30,b=60)
    ) %>% config(displayModeBar=FALSE)
  })

  # --------------------------------------------------------------------------
  # TAB 6: BIOMARKER DASHBOARD
  # --------------------------------------------------------------------------
  # Adjust PCSK9 by genetic variant
  pcsk9_adjusted <- reactive({
    base <- input$pcsk9_base
    switch(input$pcsk9_variant,
           "wt"  = base,
           "gof" = base * 1.8,   # GOF: ~80% higher
           "lof" = base * 0.45   # LOF: ~55% lower
    )
  })

  db_pd_data <- reactive({
    simulate_pd_lipids(
      drug_name     = input$db_drug,
      baseline_ldl  = input$ldl_base,
      baseline_hdl  = input$hdl_base,
      baseline_tg   = input$tg_base,
      baseline_apob = round(input$ldl_base * 0.75),
      baseline_pcsk9= pcsk9_adjusted(),
      duration_weeks= 52
    )
  })

  output$db_lipid_trends <- renderPlotly({
    df      <- db_pd_data()
    targets <- get_ldl_target(effective_risk_cat())

    plot_ly() %>%
      add_lines(data=df, x=~week, y=~LDL_C,   name="LDL-C",   line=list(color='#d62728', width=2.5)) %>%
      add_lines(data=df, x=~week, y=~HDL_C,   name="HDL-C",   line=list(color='#2ca02c', width=2.5)) %>%
      add_lines(data=df, x=~week, y=~Non_HDL, name="Non-HDL", line=list(color='#ff7f0e', width=2)) %>%
      add_lines(data=df, x=~week, y=~ApoB,    name="ApoB",    line=list(color='#9467bd', width=2)) %>%
      add_lines(x=c(0,52), y=c(targets$esc, targets$esc),
                line=list(color='blue', dash='dot'), name=paste0("ESC <",targets$esc)) %>%
      add_lines(x=c(0,52), y=c(targets$acc, targets$acc),
                line=list(color='darkgreen', dash='dash'), name=paste0("ACC <",targets$acc)) %>%
      add_lines(x=c(0,52), y=c(40,40),
                line=list(color='purple', dash='longdash'), name="ACC extreme <40") %>%
      layout(
        xaxis  = list(title="Week"),
        yaxis  = list(title="mg/dL"),
        legend = list(orientation="h", y=-0.3),
        hovermode="x unified",
        margin = list(l=60,r=20,t=20,b=80)
      ) %>% config(displayModeBar=FALSE)
  })

  output$db_pcsk9_ldlr <- renderPlotly({
    df <- db_pd_data()

    p <- plot_ly() %>%
      add_lines(data=df, x=~week, y=~PCSK9,
                name="PCSK9 (ng/mL)",
                yaxis="y",
                line=list(color='#1f77b4', width=2.5)) %>%
      add_lines(data=df, x=~week, y=~LDLR_density,
                name="LDLR Density (%)",
                yaxis="y2",
                line=list(color='#e377c2', width=2.5, dash='dash'))

    variant_text <- switch(input$pcsk9_variant,
                           "wt"  = "Wild-type PCSK9",
                           "gof" = "GOF variant: elevated baseline PCSK9",
                           "lof" = "LOF variant: reduced baseline PCSK9 (naturally lower LDL-C)")

    p %>% layout(
      title  = list(text=variant_text, font=list(size=11)),
      xaxis  = list(title="Week"),
      yaxis  = list(title="PCSK9 (ng/mL)", side="left"),
      yaxis2 = list(title="LDLR Density (%)", overlaying="y", side="right"),
      legend = list(orientation="h", y=-0.25),
      margin = list(l=60,r=60,t=40,b=70)
    ) %>% config(displayModeBar=FALSE)
  })

  myopathy_risk_val <- reactive({
    pk <- drug_pk_params[[input$db_drug]]
    if (pk$class != "Statin") return(list(score=0, label="N/A (not a statin)"))
    dose_intensity <- if (grepl("80mg|40mg", input$db_drug)) "high" else "moderate"
    compute_myopathy_risk(
      age             = input$age,
      sex             = input$sex,
      dose_intensity  = dose_intensity,
      ckd             = input$ckd_in,
      hypothyroid     = input$db_hypothyroid,
      cyp3a4_inhibitor= input$db_cyp_inhibitor,
      prior_myopathy  = input$db_prior_myo
    )
  })

  output$myopathy_risk_ui <- renderUI({
    res <- myopathy_risk_val()
    risk_css <- if (grepl("Very High", res$label)) "risk-vh" else
                if (grepl("^High",     res$label)) "risk-high" else
                if (grepl("Moderate",  res$label)) "risk-int" else "risk-low"
    div(
      div(class=paste("risk-box", risk_css),
          h4(res$label),
          p(paste("Risk score:", res$score, "/ 14"))
      ),
      br(),
      tags$ul(
        if (input$age >= 75) tags$li("Age ≥75: +2 pts") else NULL,
        if (input$sex == "Female") tags$li("Female sex: +1 pt") else NULL,
        if (input$ckd_in) tags$li("CKD: +2 pts") else NULL,
        if (input$db_hypothyroid) tags$li("Hypothyroidism: +2 pts") else NULL,
        if (input$db_cyp_inhibitor) tags$li("CYP3A4 inhibitor: +3 pts") else NULL,
        if (input$db_prior_myo) tags$li("Prior myopathy: +4 pts") else NULL
      )
    )
  })

  output$myopathy_gauge <- renderPlotly({
    res    <- myopathy_risk_val()
    score  <- min(res$score, 14)
    color  <- if (score >= 8) "#d62728" else if (score >= 5) "#ff7f0e" else
              if (score >= 3) "#ffd700" else "#2ca02c"

    plot_ly(
      type  = "indicator",
      mode  = "gauge+number+delta",
      value = score,
      gauge = list(
        axis    = list(range=list(0, 14), tickwidth=1),
        bar     = list(color=color),
        steps   = list(
          list(range=c(0,3),  color="#d4edda"),
          list(range=c(3,5),  color="#fff3cd"),
          list(range=c(5,8),  color="#f8d7da"),
          list(range=c(8,14), color="#f5c6cb")
        ),
        threshold = list(line=list(color="red", width=4), value=8)
      ),
      title  = list(text="Myopathy Risk Score"),
      number = list(suffix="/14")
    ) %>% layout(margin=list(l=20,r=20,t=30,b=20)) %>% config(displayModeBar=FALSE)
  })

  output$ddi_result_ui <- renderUI({
    if (input$ddi_drug == "None") {
      return(div(class="risk-box risk-low", h5("No co-medication selected")))
    }
    row <- cyp3a4_interactions[cyp3a4_interactions$Drug == input$ddi_drug, ]
    if (nrow(row) == 0) return(div(p("No interaction data.")))
    risk_css <- switch(row$Risk_level,
      "Contraindicated" = "risk-vh",
      "High"            = "risk-high",
      "Moderate"        = "risk-int",
      "Low"             = "risk-low"
    )
    div(
      div(class=paste("risk-box", risk_css),
          h5(paste("Interaction Level:", row$Risk_level)),
          p(paste("Mechanism:", row$Mechanism)),
          p(paste("Statin exposure increase:", row$Statin_increase))
      ),
      br(),
      div(class="target-box",
          tags$b("Recommended Action: "),
          row$Action
      )
    )
  })

  output$ddi_table <- renderDT({
    datatable(
      cyp3a4_interactions %>% select(Drug, Mechanism, Statin_increase, Risk_level, Action),
      options = list(pageLength=5, dom='tip', scrollX=TRUE),
      rownames = FALSE,
      colnames = c("Co-medication","Mechanism","Statin ↑","Risk","Action")
    ) %>%
      formatStyle("Risk_level",
                  backgroundColor = styleEqual(
                    c("Contraindicated","High","Moderate","Low"),
                    c("#f5c6cb","#f8d7da","#fff3cd","#d4edda")
                  ))
  })

  # --------------------------------------------------------------------------
  # Auto-trigger scenario comparison on startup
  # --------------------------------------------------------------------------
  observe({
    # Run scenarios immediately on start
    if (!is.null(input$run_scenarios) && input$run_scenarios == 0) {
      shinyjs::click("run_scenarios")
    }
  })

}

# ============================================================================
# RUN APP
# ============================================================================
shinyApp(ui = ui, server = server)
