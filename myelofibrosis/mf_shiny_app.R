# ============================================================
# Myelofibrosis QSP Shiny Application
# Quantitative Systems Pharmacology Model for Myelofibrosis
# JAK-STAT Pathway · Cytokine Network · Clinical Endpoints
# ============================================================
# Packages: shiny, bslib, ggplot2, plotly, DT
# Run: shiny::runApp("/home/user/qsp/myelofibrosis/mf_shiny_app.R")
# ============================================================

library(shiny)
library(bslib)
library(ggplot2)
library(plotly)
library(DT)

# ============================================================
# PK PARAMETERS (literature-based)
# ============================================================
pk_params <- list(
  Ruxolitinib = list(
    name = "Ruxolitinib", cl = 17.7, v = 72, ka = 2.0,
    f = 0.95, t_half = 3.0, protein_bind = 0.97,
    cmax_ref = 1467, auc_ref = 3692,  # nM, nM*h at 20mg BID
    mw = 306.4, ic50_jak1 = 3.3, ic50_jak2 = 2.8,  # nM
    ref_dose = 20, ref_interval = 12
  ),
  Fedratinib = list(
    name = "Fedratinib", cl = 12.0, v = 1770, ka = 0.8,
    f = 0.80, t_half = 41.0, protein_bind = 0.92,
    cmax_ref = 4200, auc_ref = 108000,  # nM at 400mg QD
    mw = 525.6, ic50_jak2 = 22, ic50_flt3 = 15,
    ref_dose = 400, ref_interval = 24
  ),
  Pacritinib = list(
    name = "Pacritinib", cl = 8.5, v = 2300, ka = 0.6,
    f = 0.70, t_half = 48.0, protein_bind = 0.97,
    cmax_ref = 290, auc_ref = 4600,  # nM at 200mg BID
    mw = 473.6, ic50_jak2 = 23, ic50_irak1 = 10,
    ref_dose = 200, ref_interval = 12
  ),
  Momelotinib = list(
    name = "Momelotinib", cl = 22.0, v = 105, ka = 1.5,
    f = 0.85, t_half = 6.0, protein_bind = 0.97,
    cmax_ref = 1250, auc_ref = 6300,  # nM at 200mg QD
    mw = 375.4, ic50_jak1 = 11, ic50_jak2 = 18, ic50_acvr1 = 0.5,
    ref_dose = 200, ref_interval = 24
  )
)

# ============================================================
# SIMULATION FUNCTIONS
# ============================================================

# One-compartment oral PK model (steady-state)
sim_pk_ss <- function(drug, dose_mg, interval_h, n_cycles = 5) {
  p  <- pk_params[[drug]]
  dose_nmol <- (dose_mg / p$mw) * 1e6  # nmol
  cl <- p$cl  # L/h
  v  <- p$v   # L
  ka <- p$ka  # h^-1
  ke <- cl / v
  F  <- p$f

  times <- seq(0, interval_h * n_cycles, by = 0.1)

  cp <- numeric(length(times))
  for (i in seq_along(times)) {
    t <- times[i]
    # Superposition of multiple doses
    cp_t <- 0
    for (n in 0:(n_cycles - 1)) {
      t_after <- t - n * interval_h
      if (t_after > 0) {
        cp_t <- cp_t + (F * dose_nmol * ka / (v * (ka - ke))) *
          (exp(-ke * t_after) - exp(-ka * t_after))
      }
    }
    cp[i] <- max(0, cp_t)
  }

  # Extract last dose cycle for steady-state display
  last_start <- (n_cycles - 1) * interval_h
  idx <- times >= last_start
  t_ss  <- times[idx] - last_start
  cp_ss <- cp[idx]

  cmax <- max(cp_ss)
  cmin <- min(cp_ss[t_ss > 0.5])
  auc  <- sum(diff(t_ss) * (head(cp_ss, -1) + tail(cp_ss, -1)) / 2)

  list(time = t_ss, cp = cp_ss,
       cmax = cmax, cmin = cmin, auc = auc, thalf = log(2) / (cl / v))
}

# JAK inhibition via sigmoid Emax
jak_inhibition <- function(cp, ic50, hill = 1.5) {
  100 * cp^hill / (ic50^hill + cp^hill)
}

# pSTAT3/pSTAT5 inhibition dynamics
sim_pd_biomarkers <- function(drug, dose_mg, interval_h, weeks = 52) {
  p    <- pk_params[[drug]]
  pk   <- sim_pk_ss(drug, dose_mg, interval_h, n_cycles = 3)
  times_wk <- seq(0, weeks, by = 0.5)

  # Use average Cp at steady state as approximation for PD
  cp_avg <- pk$auc / interval_h

  # Gradual VAF change: exponential + plateau
  vaf_baseline <- 0.55
  vaf_response <- switch(drug,
    Ruxolitinib = 0.18,
    Fedratinib  = 0.20,
    Pacritinib  = 0.15,
    Momelotinib = 0.16
  )
  vaf <- vaf_baseline - vaf_response * (1 - exp(-times_wk / 20))
  vaf <- pmax(vaf, vaf_baseline - vaf_response)

  ic50_s3 <- if (!is.null(p$ic50_jak2)) p$ic50_jak2 else 10
  ic50_s5 <- if (!is.null(p$ic50_jak1)) p$ic50_jak1 else 10

  pstat3_inh <- jak_inhibition(cp_avg, ic50_s3) * (1 - exp(-times_wk / 2))
  pstat5_inh <- jak_inhibition(cp_avg, ic50_s5) * (1 - exp(-times_wk / 1.5))

  # Cytokines: IL-6 and TNF-α baseline → suppression
  il6_base <- 35; tnfa_base <- 22  # pg/mL
  il6   <- il6_base  * exp(-0.025 * times_wk * pstat3_inh / 100) +
            rnorm(length(times_wk), 0, 0.5)
  tnfa  <- tnfa_base * exp(-0.020 * times_wk * pstat3_inh / 100) +
            rnorm(length(times_wk), 0, 0.3)

  data.frame(
    week       = times_wk,
    pstat3_inh = pmin(pstat3_inh, 95),
    pstat5_inh = pmin(pstat5_inh, 95),
    vaf        = pmax(vaf, 0.05),
    il6        = pmax(il6, 2),
    tnfa       = pmax(tnfa, 1)
  )
}

# Clinical endpoints simulation
sim_clinical <- function(drug, dose_mg, interval_h,
                         spleen_base, hgb_base, plt_base,
                         tss_base, weeks = 24) {

  p  <- pk_params[[drug]]
  pk <- sim_pk_ss(drug, dose_mg, interval_h, n_cycles = 3)

  times_wk <- seq(0, weeks, by = 0.5)
  cp_avg <- pk$auc / interval_h

  # Drug efficacy parameters
  eff <- switch(drug,
    Ruxolitinib = list(svr = 0.42, tss50 = 0.46, hgb_delta = -0.5, plt_mult = 0.90),
    Fedratinib  = list(svr = 0.47, tss50 = 0.40, hgb_delta = -0.8, plt_mult = 0.75),
    Pacritinib  = list(svr = 0.29, tss50 = 0.25, hgb_delta =  0.3, plt_mult = 0.85),
    Momelotinib = list(svr = 0.26, tss50 = 0.28, hgb_delta =  1.5, plt_mult = 0.88)
  )

  # Spleen volume: rapid initial response then plateau
  svol <- spleen_base * (1 - eff$svr * (1 - exp(-times_wk / 6)))
  svol <- svol + rnorm(length(times_wk), 0, spleen_base * 0.01)

  # Hemoglobin: initial drop then partial recovery
  hgb <- hgb_base + eff$hgb_delta * (1 - exp(-times_wk / 4)) +
         0.01 * times_wk * abs(eff$hgb_delta)
  hgb <- hgb + rnorm(length(times_wk), 0, 0.1)

  # Platelets
  plt <- plt_base * (eff$plt_mult + (1 - eff$plt_mult) * exp(-times_wk / 8))
  plt <- plt + rnorm(length(times_wk), 0, plt_base * 0.02)

  # TSS: gradual improvement
  tss <- tss_base * (1 - eff$tss50 * 0.7 * (1 - exp(-times_wk / 5)))
  tss <- pmax(tss + rnorm(length(times_wk), 0, 0.3), 0)

  # MF grade (0-3): very slow change
  mf_grades <- c("MF-0", "MF-1", "MF-2", "MF-3")
  mf_val <- 2.5 - 0.015 * times_wk  # continuous proxy
  mf_val <- pmax(pmin(mf_val + rnorm(length(times_wk), 0, 0.05), 3), 0)

  data.frame(
    week  = times_wk,
    svol  = pmax(svol, 50),
    hgb   = pmax(hgb, 4),
    plt   = pmax(plt, 10),
    tss   = pmax(tss, 0),
    mf    = mf_val
  )
}

# Multi-drug scenario comparison
scenario_compare <- function(spleen_base, hgb_base, plt_base, tss_base) {
  scenarios <- list(
    list(label = "No Treatment",         drug = NULL,         dose = 0),
    list(label = "Ruxolitinib 20mg BID", drug = "Ruxolitinib", dose = 20, int = 12),
    list(label = "Ruxolitinib 15mg BID", drug = "Ruxolitinib", dose = 15, int = 12),
    list(label = "Fedratinib 400mg QD",  drug = "Fedratinib",  dose = 400, int = 24),
    list(label = "Pacritinib 200mg BID", drug = "Pacritinib",  dose = 200, int = 12),
    list(label = "Momelotinib 200mg QD", drug = "Momelotinib", dose = 200, int = 24)
  )

  results <- lapply(scenarios, function(s) {
    if (is.null(s$drug)) {
      svr35 <- 0; tss50 <- 0; hgb_ch <- 0
      spleen_final <- spleen_base * 1.05
    } else {
      eff <- switch(s$drug,
        Ruxolitinib = if (s$dose == 20) list(svr = 0.42, tss50 = 0.46, hgb = -0.5)
                      else list(svr = 0.33, tss50 = 0.38, hgb = -0.3),
        Fedratinib  = list(svr = 0.47, tss50 = 0.40, hgb = -0.8),
        Pacritinib  = list(svr = 0.29, tss50 = 0.25, hgb =  0.3),
        Momelotinib = list(svr = 0.26, tss50 = 0.28, hgb =  1.5)
      )
      svr35 <- eff$svr * 100
      tss50 <- eff$tss50 * 100
      hgb_ch <- eff$hgb
      spleen_final <- spleen_base * (1 - eff$svr)
    }
    data.frame(label = s$label, svr35 = svr35, tss50 = tss50,
               hgb_ch = hgb_ch, spleen_pct = (spleen_final - spleen_base) / spleen_base * 100)
  })
  do.call(rbind, results)
}

# DIPSS Plus scoring
calc_dipss_plus <- function(age, hgb, wbc, blasts, symptoms, plt, transfuse, karyotype) {
  score <- 0
  score <- score + ifelse(age > 65, 1, 0)
  score <- score + ifelse(hgb < 10, 2, 0)
  score <- score + ifelse(wbc > 25, 1, 0)
  score <- score + ifelse(blasts >= 1, 1, 0)
  score <- score + ifelse(symptoms == "Yes", 1, 0)
  score <- score + ifelse(plt < 100, 2, 0)
  score <- score + ifelse(transfuse == "Yes", 1, 0)
  score <- score + ifelse(karyotype %in% c("Unfavorable", "Very unfavorable"), 1, 0)

  risk <- if (score == 0) "Low" else if (score <= 2) "Int-1" else if (score <= 4) "Int-2" else "High"
  list(score = score, risk = risk)
}

# ============================================================
# UI
# ============================================================
ui <- page_navbar(
  title = tags$span(
    tags$img(src = "https://img.icons8.com/color/48/000000/blood.png",
             height = "30px", style = "margin-right:8px;"),
    "Myelofibrosis QSP Dashboard"
  ),
  theme = bs_theme(
    version = 5,
    bg = "#0f1117", fg = "#e8eaf6",
    primary = "#7c4dff", secondary = "#536dfe",
    success = "#00e676", warning = "#ffab40", danger = "#ff5252",
    base_font = font_google("Inter"),
    heading_font = font_google("Rajdhani"),
    "card-bg" = "#1a1d2e",
    "card-border-color" = "#2d3561",
    "input-bg" = "#252840",
    "input-color" = "#e8eaf6",
    "input-border-color" = "#7c4dff"
  ),
  fillable = TRUE,

  # --------------------------------------------------------
  # TAB 1: Patient Profile
  # --------------------------------------------------------
  nav_panel(
    title = "Patient Profile",
    icon  = icon("user-circle"),
    layout_sidebar(
      sidebar = sidebar(
        width = 310, open = TRUE,
        bg = "#1a1d2e",
        tags$h5("Demographics", style = "color:#7c4dff; font-weight:700;"),
        numericInput("age",     "Age (years)", 67, 18, 100, 1),
        selectInput("sex",      "Sex", c("Male", "Female")),
        numericInput("dis_dur", "Disease Duration (months)", 18, 0, 360, 1),

        tags$hr(style = "border-color:#2d3561;"),
        tags$h5("Mutation Profile", style = "color:#7c4dff; font-weight:700;"),
        selectInput("mutation", "Driver Mutation",
          c("JAK2 V617F", "CALR Type 1", "CALR Type 2", "MPL W515", "Triple Negative")),
        numericInput("vaf_bl", "JAK2 VAF / CALR VAF (%)", 55, 0, 100, 1),

        tags$hr(style = "border-color:#2d3561;"),
        tags$h5("Baseline Labs", style = "color:#7c4dff; font-weight:700;"),
        numericInput("hgb_bl",  "Hemoglobin (g/dL)",     9.5, 4, 20, 0.1),
        numericInput("plt_bl",  "Platelets (×10⁹/L)",    150, 10, 1500, 10),
        numericInput("wbc_bl",  "WBC (×10⁹/L)",          8.0, 1, 100, 0.5),
        numericInput("blasts",  "Peripheral Blasts (%)",   1, 0, 30, 1),
        numericInput("svol_bl", "Spleen Volume (cm³)",    800, 100, 5000, 50),
        numericInput("tss_bl",  "Total Symptom Score",     18, 0, 50, 1),
        selectInput("constit",  "Constitutional Symptoms","Yes"),

        tags$hr(style = "border-color:#2d3561;"),
        tags$h5("Risk Factors", style = "color:#7c4dff; font-weight:700;"),
        selectInput("transfuse",  "Transfusion Dependent", c("No", "Yes")),
        selectInput("karyotype",  "Cytogenetics",
          c("Favorable", "Unfavorable", "Very unfavorable")),
        selectInput("mf_grade",   "BM Fibrosis Grade", c("MF-1","MF-2","MF-3")),
        selectInput("tx_hist",    "Prior Treatment",
          c("Treatment naive", "Prior Ruxolitinib", "Prior Fedratinib", "Multiple prior JAKi"))
      ),

      # Main panel
      layout_columns(
        col_widths = c(6, 6, 12),
        card(
          card_header("DIPSS Plus Risk Stratification"),
          card_body(
            uiOutput("dipss_result"),
            tags$br(),
            uiOutput("risk_badge")
          )
        ),
        card(
          card_header("Mutation Risk Profile"),
          card_body(plotlyOutput("mutation_radar", height = "250px"))
        ),
        card(
          card_header("Patient Summary"),
          card_body(DTOutput("patient_table"))
        )
      )
    )
  ),

  # --------------------------------------------------------
  # TAB 2: Pharmacokinetics
  # --------------------------------------------------------
  nav_panel(
    title = "Pharmacokinetics",
    icon  = icon("capsules"),
    layout_sidebar(
      sidebar = sidebar(
        width = 280, bg = "#1a1d2e",
        tags$h5("Drug Selection", style = "color:#7c4dff; font-weight:700;"),
        selectInput("pk_drug", "JAK Inhibitor",
          c("Ruxolitinib", "Fedratinib", "Pacritinib", "Momelotinib")),
        numericInput("pk_dose", "Dose (mg)", 20, 1, 800, 1),
        numericInput("pk_int",  "Dosing Interval (h)", 12, 6, 48, 6),
        selectInput("pk_unit", "Concentration Unit", c("nM", "ng/mL", "µM")),
        tags$hr(style = "border-color:#2d3561;"),
        actionButton("run_pk", "Run PK Simulation",
          class = "btn-primary w-100",
          icon  = icon("play-circle"))
      ),
      layout_columns(
        col_widths = c(8, 4),
        card(
          card_header("Steady-State Concentration-Time Profile"),
          card_body(plotlyOutput("pk_plot", height = "350px"))
        ),
        card(
          card_header("PK Parameters"),
          card_body(
            uiOutput("pk_params_ui"),
            tags$br(),
            DTOutput("pk_table")
          )
        )
      ),
      layout_columns(
        col_widths = c(6, 6),
        card(
          card_header("JAK1/JAK2 Inhibition vs Concentration"),
          card_body(plotlyOutput("jak_inh_plot", height = "280px"))
        ),
        card(
          card_header("Multi-Dose Accumulation"),
          card_body(plotlyOutput("pk_accum_plot", height = "280px"))
        )
      )
    )
  ),

  # --------------------------------------------------------
  # TAB 3: PD Biomarkers
  # --------------------------------------------------------
  nav_panel(
    title = "PD Biomarkers",
    icon  = icon("dna"),
    layout_sidebar(
      sidebar = sidebar(
        width = 260, bg = "#1a1d2e",
        tags$h5("Simulation Settings", style = "color:#7c4dff; font-weight:700;"),
        selectInput("pd_drug", "Drug", c("Ruxolitinib","Fedratinib","Pacritinib","Momelotinib")),
        numericInput("pd_dose", "Dose (mg)", 20, 1, 800, 1),
        numericInput("pd_int",  "Interval (h)", 12, 6, 48, 6),
        numericInput("pd_weeks","Duration (weeks)", 52, 4, 104, 4),
        actionButton("run_pd", "Simulate Biomarkers",
          class = "btn-primary w-100", icon = icon("chart-line"))
      ),
      layout_columns(
        col_widths = c(6, 6),
        card(
          card_header("pSTAT3 & pSTAT5 Inhibition (%)"),
          card_body(plotlyOutput("pstat_plot", height = "280px"))
        ),
        card(
          card_header("JAK2 V617F Allele Burden (VAF)"),
          card_body(plotlyOutput("vaf_plot", height = "280px"))
        )
      ),
      layout_columns(
        col_widths = c(6, 6),
        card(
          card_header("Cytokine Levels over Time"),
          card_body(plotlyOutput("cytokine_plot", height = "280px"))
        ),
        card(
          card_header("Dose-Response: JAK Inhibition"),
          card_body(plotlyOutput("dr_curve_plot", height = "280px"))
        )
      )
    )
  ),

  # --------------------------------------------------------
  # TAB 4: Clinical Endpoints
  # --------------------------------------------------------
  nav_panel(
    title = "Clinical Endpoints",
    icon  = icon("stethoscope"),
    layout_sidebar(
      sidebar = sidebar(
        width = 260, bg = "#1a1d2e",
        tags$h5("Simulation Settings", style = "color:#7c4dff; font-weight:700;"),
        selectInput("ce_drug",  "Drug", c("Ruxolitinib","Fedratinib","Pacritinib","Momelotinib")),
        numericInput("ce_dose",  "Dose (mg)", 20, 1, 800, 1),
        numericInput("ce_int",   "Interval (h)", 12, 6, 48, 6),
        numericInput("ce_weeks", "Weeks", 24, 4, 52, 4),
        actionButton("run_ce", "Simulate Endpoints",
          class = "btn-primary w-100", icon = icon("chart-area"))
      ),
      layout_columns(
        col_widths = c(6, 6),
        card(
          card_header("Spleen Volume Response (SVR35 threshold = 35%)"),
          card_body(plotlyOutput("svol_plot", height = "280px"))
        ),
        card(
          card_header("Hemoglobin Dynamics"),
          card_body(plotlyOutput("hgb_plot", height = "280px"))
        )
      ),
      layout_columns(
        col_widths = c(4, 4, 4),
        card(
          card_header("Platelet Count"),
          card_body(plotlyOutput("plt_plot", height = "250px"))
        ),
        card(
          card_header("Total Symptom Score (TSS50)"),
          card_body(plotlyOutput("tss_plot", height = "250px"))
        ),
        card(
          card_header("Bone Marrow Fibrosis Grade"),
          card_body(plotlyOutput("mf_plot", height = "250px"))
        )
      )
    )
  ),

  # --------------------------------------------------------
  # TAB 5: Treatment Comparison
  # --------------------------------------------------------
  nav_panel(
    title = "Treatment Comparison",
    icon  = icon("balance-scale"),
    layout_sidebar(
      sidebar = sidebar(
        width = 260, bg = "#1a1d2e",
        tags$h5("Select Treatments", style = "color:#7c4dff; font-weight:700;"),
        checkboxGroupInput("comp_drugs", NULL,
          choices  = c("No Treatment", "Ruxolitinib 20mg BID",
                       "Ruxolitinib 15mg BID", "Fedratinib 400mg QD",
                       "Pacritinib 200mg BID", "Momelotinib 200mg QD"),
          selected = c("No Treatment", "Ruxolitinib 20mg BID",
                       "Fedratinib 400mg QD", "Pacritinib 200mg BID",
                       "Momelotinib 200mg QD")),
        actionButton("run_comp", "Compare",
          class = "btn-primary w-100", icon = icon("chart-bar"))
      ),
      layout_columns(
        col_widths = c(6, 6),
        card(
          card_header("SVR35 & TSS50 Response Rates (%)"),
          card_body(plotlyOutput("bar_compare", height = "300px"))
        ),
        card(
          card_header("Hemoglobin Change (g/dL) at 24 weeks"),
          card_body(plotlyOutput("hgb_compare", height = "300px"))
        )
      ),
      layout_columns(
        col_widths = c(6, 6),
        card(
          card_header("Spleen Volume Waterfall Plot"),
          card_body(plotlyOutput("waterfall_plot", height = "300px"))
        ),
        card(
          card_header("Forest Plot: Treatment Effects"),
          card_body(plotlyOutput("forest_plot", height = "300px"))
        )
      )
    )
  ),

  # --------------------------------------------------------
  # TAB 6: Biomarker Dynamics
  # --------------------------------------------------------
  nav_panel(
    title = "Biomarker Dynamics",
    icon  = icon("microscope"),
    layout_sidebar(
      sidebar = sidebar(
        width = 260, bg = "#1a1d2e",
        tags$h5("Analysis Settings", style = "color:#7c4dff; font-weight:700;"),
        selectInput("bm_drug", "Drug", c("Ruxolitinib","Fedratinib","Pacritinib","Momelotinib")),
        numericInput("bm_dose", "Dose (mg)", 20, 1, 800, 1),
        numericInput("bm_int",  "Interval (h)", 12, 6, 48, 6),
        sliderInput("bm_weeks", "Weeks to Simulate", 4, 104, 52, 4),
        actionButton("run_bm", "Run Analysis",
          class = "btn-primary w-100", icon = icon("play"))
      ),
      layout_columns(
        col_widths = c(6, 6),
        card(
          card_header("VAF vs Spleen Volume Response Correlation"),
          card_body(plotlyOutput("vaf_svol_corr", height = "280px"))
        ),
        card(
          card_header("Cytokine Heatmap by Treatment"),
          card_body(plotlyOutput("cytokine_heat", height = "280px"))
        )
      ),
      layout_columns(
        col_widths = c(6, 6),
        card(
          card_header("pSTAT3/pSTAT5 Temporal Dynamics"),
          card_body(plotlyOutput("pstat_dyn", height = "280px"))
        ),
        card(
          card_header("Time to AML Transformation (KM-style)"),
          card_body(plotlyOutput("km_aml", height = "280px"))
        )
      )
    )
  ),

  # Footer
  nav_spacer(),
  nav_item(
    tags$small("MF QSP v1.0 | mrgsolve-compatible model", style = "color:#536dfe;")
  )
)

# Null-coalescing operator (defined globally so all server code can use it)
`%||%` <- function(a, b) if (!is.null(a)) a else b

# ============================================================
# SERVER
# ============================================================
server <- function(input, output, session) {

  # ---- Reactive: update defaults when drug changes ----------
  observeEvent(input$pk_drug, {
    p <- pk_params[[input$pk_drug]]
    updateNumericInput(session, "pk_dose", value = p$ref_dose)
    updateNumericInput(session, "pk_int",  value = p$ref_interval)
  })
  observeEvent(input$pd_drug, {
    p <- pk_params[[input$pd_drug]]
    updateNumericInput(session, "pd_dose", value = p$ref_dose)
    updateNumericInput(session, "pd_int",  value = p$ref_interval)
  })
  observeEvent(input$ce_drug, {
    p <- pk_params[[input$ce_drug]]
    updateNumericInput(session, "ce_dose", value = p$ref_dose)
    updateNumericInput(session, "ce_int",  value = p$ref_interval)
  })
  observeEvent(input$bm_drug, {
    p <- pk_params[[input$bm_drug]]
    updateNumericInput(session, "bm_dose", value = p$ref_dose)
    updateNumericInput(session, "bm_int",  value = p$ref_interval)
  })

  # ===========================================================
  # TAB 1: Patient Profile
  # ===========================================================
  dipss_rv <- reactive({
    calc_dipss_plus(input$age, input$hgb_bl, input$wbc_bl, input$blasts,
                    input$constit, input$plt_bl, input$transfuse, input$karyotype)
  })

  output$dipss_result <- renderUI({
    d <- dipss_rv()
    tags$div(
      tags$h3(paste("DIPSS Plus Score:", d$score),
              style = "color:#ffab40; font-weight:700;"),
      tags$p(paste("Components contributing to score:"),
             style = "font-size:0.85rem; color:#aaa;"),
      tags$ul(style = "font-size:0.85rem;",
        tags$li(paste("Age >65:", if(input$age>65) "1 pt" else "0 pt")),
        tags$li(paste("Hgb <10 g/dL:", if(input$hgb_bl<10) "2 pts" else "0 pt")),
        tags$li(paste("WBC >25 ×10⁹/L:", if(input$wbc_bl>25) "1 pt" else "0 pt")),
        tags$li(paste("Blasts ≥1%:", if(input$blasts>=1) "1 pt" else "0 pt")),
        tags$li(paste("Constitutional Sx:", if(input$constit=="Yes") "1 pt" else "0 pt")),
        tags$li(paste("Plt <100:", if(input$plt_bl<100) "2 pts" else "0 pt")),
        tags$li(paste("Transfusion dep.:", if(input$transfuse=="Yes") "1 pt" else "0 pt")),
        tags$li(paste("Unfav. karyotype:", if(input$karyotype!="Favorable") "1 pt" else "0 pt"))
      )
    )
  })

  output$risk_badge <- renderUI({
    d <- dipss_rv()
    col <- switch(d$risk,
      "Low"   = "#00e676",
      "Int-1" = "#ffeb3b",
      "Int-2" = "#ffab40",
      "High"  = "#ff5252"
    )
    tags$div(
      style = paste0("background:", col,
                     "; color:#000; padding:12px 20px; border-radius:8px;",
                     " font-weight:700; font-size:1.4rem; text-align:center;"),
      paste("Risk Category:", d$risk)
    )
  })

  output$mutation_radar <- renderPlotly({
    mut_scores <- switch(input$mutation,
      "JAK2 V617F"    = c(JAK2=9, CALR=1, MPL=1, Thrombosis=8, AML=5, BM_Fibrosis=7),
      "CALR Type 1"   = c(JAK2=1, CALR=9, MPL=1, Thrombosis=3, AML=3, BM_Fibrosis=5),
      "CALR Type 2"   = c(JAK2=1, CALR=7, MPL=1, Thrombosis=4, AML=7, BM_Fibrosis=6),
      "MPL W515"      = c(JAK2=1, CALR=1, MPL=9, Thrombosis=6, AML=4, BM_Fibrosis=6),
      "Triple Negative"= c(JAK2=1, CALR=1, MPL=1, Thrombosis=3, AML=8, BM_Fibrosis=8)
    )
    cats <- names(mut_scores)
    vals <- as.numeric(mut_scores)
    plot_ly(type = "scatterpolar", fill = "toself",
            r = c(vals, vals[1]), theta = c(cats, cats[1]),
            fillcolor = "rgba(124,77,255,0.3)",
            line = list(color = "#7c4dff", width = 2),
            hoverinfo = "r+theta") %>%
      layout(polar = list(radialaxis = list(visible = TRUE, range = c(0, 10),
                                            color = "#aaa"),
                          bgcolor = "#1a1d2e"),
             paper_bgcolor = "#1a1d2e", plot_bgcolor = "#1a1d2e",
             font = list(color = "#e8eaf6"), showlegend = FALSE,
             margin = list(l = 30, r = 30, t = 20, b = 20))
  })

  output$patient_table <- renderDT({
    d <- dipss_rv()
    df <- data.frame(
      Parameter = c("Age", "Sex", "Disease Duration", "Driver Mutation",
                    "JAK2/CALR VAF", "Spleen Volume", "Hemoglobin",
                    "Platelets", "WBC", "TSS", "MF Grade",
                    "Prior Treatment", "DIPSS Plus Score", "Risk Category"),
      Value = c(
        paste(input$age, "years"),
        input$sex,
        paste(input$dis_dur, "months"),
        input$mutation,
        paste0(input$vaf_bl, "%"),
        paste(input$svol_bl, "cm³"),
        paste(input$hgb_bl, "g/dL"),
        paste(input$plt_bl, "×10⁹/L"),
        paste(input$wbc_bl, "×10⁹/L"),
        as.character(input$tss_bl),
        input$mf_grade,
        input$tx_hist,
        as.character(d$score),
        d$risk
      )
    )
    datatable(df, options = list(dom = "t", pageLength = 20, scrollY = "300px"),
              rownames = FALSE,
              style = "bootstrap5") %>%
      formatStyle("Parameter", fontWeight = "bold", color = "#7c4dff") %>%
      formatStyle("Value", color = "#e8eaf6")
  })

  # ===========================================================
  # TAB 2: Pharmacokinetics
  # ===========================================================
  pk_rv <- eventReactive(input$run_pk, {
    req(input$pk_drug, input$pk_dose, input$pk_int)
    set.seed(42)
    sim_pk_ss(input$pk_drug, input$pk_dose, input$pk_int, n_cycles = 5)
  }, ignoreNULL = FALSE)

  unit_conv <- reactive({
    switch(input$pk_unit,
      "nM"    = list(factor = 1,    label = "nM"),
      "ng/mL" = list(factor = pk_params[[input$pk_drug]]$mw / 1000,
                     label = "ng/mL"),
      "µM"    = list(factor = 0.001, label = "µM")
    )
  })

  output$pk_plot <- renderPlotly({
    pk <- pk_rv()
    uc <- unit_conv()
    cp_display <- pk$cp * uc$factor
    ic50_j2 <- if (!is.null(pk_params[[input$pk_drug]]$ic50_jak2))
                 pk_params[[input$pk_drug]]$ic50_jak2 * uc$factor else NULL
    ic50_j1 <- if (!is.null(pk_params[[input$pk_drug]]$ic50_jak1))
                 pk_params[[input$pk_drug]]$ic50_jak1 * uc$factor else NULL

    p <- plot_ly() %>%
      add_trace(x = pk$time, y = cp_display, type = "scatter", mode = "lines",
                line = list(color = "#7c4dff", width = 2.5),
                name = paste(input$pk_drug, "Cp"),
                hovertemplate = paste0("Time: %{x:.1f}h<br>Cp: %{y:.1f} ", uc$label))

    if (!is.null(ic50_j2))
      p <- p %>% add_lines(x = range(pk$time), y = c(ic50_j2, ic50_j2),
                            line = list(color = "#ff5252", dash = "dash", width = 1.5),
                            name = "IC50 JAK2")
    if (!is.null(ic50_j1))
      p <- p %>% add_lines(x = range(pk$time), y = c(ic50_j1, ic50_j1),
                            line = list(color = "#ffab40", dash = "dash", width = 1.5),
                            name = "IC50 JAK1")

    p %>% layout(
      xaxis = list(title = "Time (h)", color = "#aaa", gridcolor = "#2d3561"),
      yaxis = list(title = paste0("Plasma Concentration (", uc$label, ")"),
                   color = "#aaa", gridcolor = "#2d3561"),
      paper_bgcolor = "#1a1d2e", plot_bgcolor = "#1a1d2e",
      font = list(color = "#e8eaf6"), legend = list(x = 0.7, y = 0.9),
      margin = list(l = 60, r = 20, t = 20, b = 50)
    )
  })

  output$pk_params_ui <- renderUI({
    pk <- pk_rv()
    uc <- unit_conv()
    tags$div(
      style = "font-size: 0.85rem; color: #aaa;",
      tags$p(paste0("Drug: ", input$pk_drug, " | Dose: ", input$pk_dose,
                    " mg q", input$pk_int, "h"))
    )
  })

  output$pk_table <- renderDT({
    pk <- pk_rv()
    uc <- unit_conv()
    df <- data.frame(
      Parameter = c("Cmax", "Cmin", "AUC₀₋τ", "t½", "Protein Binding"),
      Value = c(
        paste0(round(pk$cmax * uc$factor, 1), " ", uc$label),
        paste0(round(pk$cmin * uc$factor, 1), " ", uc$label),
        paste0(round(pk$auc * uc$factor, 0), " ", uc$label, "·h"),
        paste0(round(pk$thalf, 1), " h"),
        paste0(pk_params[[input$pk_drug]]$protein_bind * 100, "%")
      )
    )
    datatable(df, options = list(dom = "t", pageLength = 10),
              rownames = FALSE, style = "bootstrap5") %>%
      formatStyle("Parameter", fontWeight = "bold", color = "#7c4dff")
  })

  output$jak_inh_plot <- renderPlotly({
    p  <- pk_params[[input$pk_drug]]
    cp_range <- 10^seq(-1, 4, by = 0.05)  # nM
    uc <- unit_conv()

    traces <- list()
    if (!is.null(p$ic50_jak2)) {
      inh_j2 <- jak_inhibition(cp_range, p$ic50_jak2)
      traces[[1]] <- list(x = cp_range * uc$factor, y = inh_j2,
                          name = "JAK2 Inhibition", color = "#7c4dff")
    }
    if (!is.null(p$ic50_jak1)) {
      inh_j1 <- jak_inhibition(cp_range, p$ic50_jak1)
      traces[[length(traces)+1]] <- list(x = cp_range * uc$factor, y = inh_j1,
                                         name = "JAK1 Inhibition", color = "#00e676")
    }

    pl <- plot_ly()
    for (tr in traces)
      pl <- pl %>% add_trace(x = tr$x, y = tr$y, type = "scatter", mode = "lines",
                             name = tr$name, line = list(color = tr$color, width = 2))

    pl %>% layout(
      xaxis = list(title = paste0("Concentration (", uc$label, ")"), type = "log",
                   color = "#aaa", gridcolor = "#2d3561"),
      yaxis = list(title = "Inhibition (%)", range = c(0, 100),
                   color = "#aaa", gridcolor = "#2d3561"),
      paper_bgcolor = "#1a1d2e", plot_bgcolor = "#1a1d2e",
      font = list(color = "#e8eaf6"),
      shapes = list(list(type="line", x0=0, x1=1, y0=50, y1=50,
                         xref="paper", yref="y",
                         line=list(color="#ff5252", dash="dot", width=1))),
      margin = list(l=60, r=20, t=20, b=50)
    )
  })

  output$pk_accum_plot <- renderPlotly({
    pk_full <- sim_pk_ss(input$pk_drug, input$pk_dose, input$pk_int, n_cycles = 5)
    uc <- unit_conv()
    # Re-simulate full 5-dose profile
    p  <- pk_params[[input$pk_drug]]
    dose_nmol <- (input$pk_dose / p$mw) * 1e6
    cl <- p$cl; v <- p$v; ka <- p$ka; ke <- cl / v; F <- p$f
    n_cyc <- 5
    times <- seq(0, input$pk_int * n_cyc, by = 0.1)
    cp <- numeric(length(times))
    for (i in seq_along(times)) {
      t <- times[i]
      for (n in 0:(n_cyc - 1)) {
        t_after <- t - n * input$pk_int
        if (t_after > 0)
          cp[i] <- cp[i] + (F*dose_nmol*ka/(v*(ka-ke))) *
                   (exp(-ke*t_after)-exp(-ka*t_after))
      }
    }
    cp <- pmax(cp, 0)

    dose_markers <- data.frame(
      x = seq(0, input$pk_int*(n_cyc-1), by=input$pk_int),
      y = rep(0, n_cyc)
    )

    plot_ly() %>%
      add_trace(x = times, y = cp * uc$factor, type = "scatter", mode = "lines",
                line = list(color = "#536dfe", width = 2),
                name = "Cp") %>%
      add_trace(x = dose_markers$x, y = dose_markers$y + max(cp)*uc$factor*0.02,
                type = "scatter", mode = "markers",
                marker = list(symbol = "triangle-up", size = 12, color = "#00e676"),
                name = "Dose") %>%
      layout(
        xaxis = list(title = "Time (h)", color = "#aaa", gridcolor = "#2d3561"),
        yaxis = list(title = paste0("Cp (", uc$label, ")"),
                     color = "#aaa", gridcolor = "#2d3561"),
        paper_bgcolor = "#1a1d2e", plot_bgcolor = "#1a1d2e",
        font = list(color = "#e8eaf6"),
        margin = list(l=60, r=20, t=20, b=50)
      )
  })

  # ===========================================================
  # TAB 3: PD Biomarkers
  # ===========================================================
  pd_rv <- eventReactive(input$run_pd, {
    req(input$pd_drug, input$pd_dose, input$pd_int, input$pd_weeks)
    set.seed(123)
    sim_pd_biomarkers(input$pd_drug, input$pd_dose, input$pd_int, input$pd_weeks)
  }, ignoreNULL = FALSE)

  output$pstat_plot <- renderPlotly({
    df <- pd_rv()
    plot_ly(df) %>%
      add_trace(x = ~week, y = ~pstat3_inh, type = "scatter", mode = "lines",
                name = "pSTAT3 inhibition", line = list(color = "#7c4dff", width = 2)) %>%
      add_trace(x = ~week, y = ~pstat5_inh, type = "scatter", mode = "lines",
                name = "pSTAT5 inhibition", line = list(color = "#00e676", width = 2)) %>%
      layout(
        xaxis = list(title = "Week", color = "#aaa", gridcolor = "#2d3561"),
        yaxis = list(title = "Inhibition (%)", range = c(0, 100),
                     color = "#aaa", gridcolor = "#2d3561"),
        paper_bgcolor = "#1a1d2e", plot_bgcolor = "#1a1d2e",
        font = list(color = "#e8eaf6"), legend = list(x = 0.5, y = 0.15),
        margin = list(l=60, r=20, t=20, b=50)
      )
  })

  output$vaf_plot <- renderPlotly({
    df <- pd_rv()
    plot_ly(df, x = ~week, y = ~vaf * 100, type = "scatter", mode = "lines",
            line = list(color = "#ffab40", width = 2.5),
            fill = "tozeroy", fillcolor = "rgba(255,171,64,0.15)",
            hovertemplate = "Week %{x}<br>VAF: %{y:.1f}%") %>%
      add_lines(x = range(df$week), y = c(50, 50),
                line = list(color = "#ff5252", dash = "dash", width = 1.5),
                name = "Baseline VAF") %>%
      layout(
        xaxis = list(title = "Week", color = "#aaa", gridcolor = "#2d3561"),
        yaxis = list(title = "JAK2 V617F VAF (%)", color = "#aaa", gridcolor = "#2d3561"),
        paper_bgcolor = "#1a1d2e", plot_bgcolor = "#1a1d2e",
        font = list(color = "#e8eaf6"), showlegend = FALSE,
        margin = list(l=60, r=20, t=20, b=50)
      )
  })

  output$cytokine_plot <- renderPlotly({
    df <- pd_rv()
    plot_ly(df) %>%
      add_trace(x = ~week, y = ~il6, type = "scatter", mode = "lines",
                name = "IL-6 (pg/mL)", line = list(color = "#ff5252", width = 2)) %>%
      add_trace(x = ~week, y = ~tnfa, type = "scatter", mode = "lines",
                name = "TNF-α (pg/mL)", line = list(color = "#ffeb3b", width = 2)) %>%
      layout(
        xaxis = list(title = "Week", color = "#aaa", gridcolor = "#2d3561"),
        yaxis = list(title = "Cytokine Level (pg/mL)", color = "#aaa", gridcolor = "#2d3561"),
        paper_bgcolor = "#1a1d2e", plot_bgcolor = "#1a1d2e",
        font = list(color = "#e8eaf6"), legend = list(x = 0.6, y = 0.9),
        margin = list(l=60, r=20, t=20, b=50)
      )
  })

  output$dr_curve_plot <- renderPlotly({
    conc_range <- 10^seq(-1, 4, by = 0.1)
    drugs_sel  <- c("Ruxolitinib", "Fedratinib", "Pacritinib", "Momelotinib")
    colors     <- c("#7c4dff", "#00e676", "#ffab40", "#ff5252")

    pl <- plot_ly()
    for (i in seq_along(drugs_sel)) {
      dr <- drugs_sel[i]
      ic50 <- pk_params[[dr]]$ic50_jak2 %||% pk_params[[dr]]$ic50_jak1 %||% 10
      inh  <- jak_inhibition(conc_range, ic50)
      pl <- pl %>% add_trace(
        x = conc_range, y = inh, type = "scatter", mode = "lines",
        name = dr, line = list(color = colors[i], width = 2)
      )
    }

    pl %>% layout(
      xaxis = list(title = "Concentration (nM)", type = "log",
                   color = "#aaa", gridcolor = "#2d3561"),
      yaxis = list(title = "JAK2 Inhibition (%)", range = c(0, 100),
                   color = "#aaa", gridcolor = "#2d3561"),
      paper_bgcolor = "#1a1d2e", plot_bgcolor = "#1a1d2e",
      font = list(color = "#e8eaf6"),
      shapes = list(list(type="line", x0=0, x1=1, y0=50, y1=50,
                         xref="paper", yref="y",
                         line=list(color="grey", dash="dot", width=1))),
      margin = list(l=60, r=20, t=20, b=50)
    )
  })

  # ===========================================================
  # TAB 4: Clinical Endpoints
  # ===========================================================
  ce_rv <- eventReactive(input$run_ce, {
    req(input$ce_drug, input$ce_dose, input$ce_int, input$ce_weeks)
    set.seed(456)
    sim_clinical(input$ce_drug, input$ce_dose, input$ce_int,
                 input$svol_bl, input$hgb_bl, input$plt_bl,
                 input$tss_bl, input$ce_weeks)
  }, ignoreNULL = FALSE)

  output$svol_plot <- renderPlotly({
    df <- ce_rv()
    svr35_thresh <- input$svol_bl * 0.65  # 35% reduction
    plot_ly(df) %>%
      add_trace(x = ~week, y = ~svol, type = "scatter", mode = "lines",
                line = list(color = "#7c4dff", width = 2.5),
                fill = "tozeroy", fillcolor = "rgba(124,77,255,0.1)",
                name = "Spleen Volume") %>%
      add_lines(x = range(df$week), y = c(svr35_thresh, svr35_thresh),
                line = list(color = "#00e676", dash = "dash", width = 2),
                name = "SVR35 Threshold") %>%
      add_lines(x = range(df$week), y = c(input$svol_bl, input$svol_bl),
                line = list(color = "#ff5252", dash = "dot", width = 1.5),
                name = "Baseline") %>%
      layout(
        xaxis = list(title = "Week", color = "#aaa", gridcolor = "#2d3561"),
        yaxis = list(title = "Spleen Volume (cm³)", color = "#aaa", gridcolor = "#2d3561"),
        paper_bgcolor = "#1a1d2e", plot_bgcolor = "#1a1d2e",
        font = list(color = "#e8eaf6"),
        margin = list(l=60, r=20, t=20, b=50)
      )
  })

  output$hgb_plot <- renderPlotly({
    df <- ce_rv()
    plot_ly(df, x = ~week, y = ~hgb, type = "scatter", mode = "lines",
            line = list(color = "#ff5252", width = 2.5),
            fill = "tozeroy", fillcolor = "rgba(255,82,82,0.1)",
            name = "Hgb") %>%
      add_lines(x = range(df$week), y = c(10, 10),
                line = list(color = "#ffab40", dash = "dash", width = 1.5),
                name = "Anemia threshold (10 g/dL)") %>%
      layout(
        xaxis = list(title = "Week", color = "#aaa", gridcolor = "#2d3561"),
        yaxis = list(title = "Hemoglobin (g/dL)", color = "#aaa", gridcolor = "#2d3561"),
        paper_bgcolor = "#1a1d2e", plot_bgcolor = "#1a1d2e",
        font = list(color = "#e8eaf6"),
        margin = list(l=60, r=20, t=20, b=50)
      )
  })

  output$plt_plot <- renderPlotly({
    df <- ce_rv()
    plot_ly(df, x = ~week, y = ~plt, type = "scatter", mode = "lines",
            line = list(color = "#ffeb3b", width = 2)) %>%
      add_ribbons(x = ~week,
                  ymin = rep(400, nrow(df)), ymax = rep(1200, nrow(df)),
                  fillcolor = "rgba(255,87,34,0.1)",
                  line = list(color = "transparent"),
                  name = "Thrombosis risk zone") %>%
      add_lines(x = range(df$week), y = c(50, 50),
                line = list(color = "#ff5252", dash = "dash", width = 1.5),
                name = "Bleeding risk (<50)") %>%
      layout(
        xaxis = list(title = "Week", color = "#aaa", gridcolor = "#2d3561"),
        yaxis = list(title = "Platelets (×10⁹/L)", color = "#aaa", gridcolor = "#2d3561"),
        paper_bgcolor = "#1a1d2e", plot_bgcolor = "#1a1d2e",
        font = list(color = "#e8eaf6"), showlegend = FALSE,
        margin = list(l=60, r=20, t=20, b=50)
      )
  })

  output$tss_plot <- renderPlotly({
    df <- ce_rv()
    tss50_thresh <- input$tss_bl * 0.5
    plot_ly(df, x = ~week, y = ~tss, type = "scatter", mode = "lines",
            line = list(color = "#00e676", width = 2.5),
            fill = "tozeroy", fillcolor = "rgba(0,230,118,0.1)") %>%
      add_lines(x = range(df$week), y = c(tss50_thresh, tss50_thresh),
                line = list(color = "#ffab40", dash = "dash", width = 2),
                name = "TSS50 threshold") %>%
      layout(
        xaxis = list(title = "Week", color = "#aaa", gridcolor = "#2d3561"),
        yaxis = list(title = "Total Symptom Score", color = "#aaa", gridcolor = "#2d3561"),
        paper_bgcolor = "#1a1d2e", plot_bgcolor = "#1a1d2e",
        font = list(color = "#e8eaf6"), showlegend = FALSE,
        margin = list(l=60, r=20, t=20, b=50)
      )
  })

  output$mf_plot <- renderPlotly({
    df <- ce_rv()
    grade_labels <- c("MF-0", "MF-1", "MF-2", "MF-3")
    plot_ly(df, x = ~week, y = ~mf, type = "scatter", mode = "lines",
            line = list(color = "#536dfe", width = 2)) %>%
      layout(
        xaxis = list(title = "Week", color = "#aaa", gridcolor = "#2d3561"),
        yaxis = list(title = "BM Fibrosis (0=MF-0, 3=MF-3)",
                     range = c(-0.1, 3.1),
                     tickvals = 0:3, ticktext = grade_labels,
                     color = "#aaa", gridcolor = "#2d3561"),
        paper_bgcolor = "#1a1d2e", plot_bgcolor = "#1a1d2e",
        font = list(color = "#e8eaf6"),
        margin = list(l=70, r=20, t=20, b=50)
      )
  })

  # ===========================================================
  # TAB 5: Treatment Comparison
  # ===========================================================
  comp_rv <- eventReactive(input$run_comp, {
    req(input$comp_drugs)
    set.seed(789)
    df <- scenario_compare(input$svol_bl, input$hgb_bl, input$plt_bl, input$tss_bl)
    df[df$label %in% input$comp_drugs, ]
  }, ignoreNULL = FALSE)

  output$bar_compare <- renderPlotly({
    df <- comp_rv()
    df_long <- rbind(
      data.frame(label = df$label, value = df$svr35, metric = "SVR35 (%)"),
      data.frame(label = df$label, value = df$tss50, metric = "TSS50 (%)")
    )
    plot_ly(df_long, x = ~label, y = ~value, color = ~metric,
            type = "bar", barmode = "group",
            colors = c("#7c4dff", "#00e676")) %>%
      layout(
        xaxis = list(title = "", color = "#aaa", tickangle = -30),
        yaxis = list(title = "Response Rate (%)", range = c(0, 70),
                     color = "#aaa", gridcolor = "#2d3561"),
        paper_bgcolor = "#1a1d2e", plot_bgcolor = "#1a1d2e",
        font = list(color = "#e8eaf6"), legend = list(x = 0.7, y = 0.95),
        margin = list(l=60, r=20, t=20, b=80)
      )
  })

  output$hgb_compare <- renderPlotly({
    df <- comp_rv()
    colors <- ifelse(df$hgb_ch >= 0, "#00e676", "#ff5252")
    plot_ly(df, x = ~label, y = ~hgb_ch, type = "bar",
            marker = list(color = colors)) %>%
      add_lines(x = c(-0.5, nrow(df)-0.5), y = c(0, 0),
                line = list(color = "white", width = 1)) %>%
      layout(
        xaxis = list(title = "", color = "#aaa", tickangle = -30),
        yaxis = list(title = "Hgb Change (g/dL)", color = "#aaa", gridcolor = "#2d3561"),
        paper_bgcolor = "#1a1d2e", plot_bgcolor = "#1a1d2e",
        font = list(color = "#e8eaf6"), showlegend = FALSE,
        margin = list(l=60, r=20, t=20, b=80)
      )
  })

  output$waterfall_plot <- renderPlotly({
    df <- comp_rv()
    df <- df[order(df$spleen_pct), ]
    colors <- ifelse(df$spleen_pct <= -35, "#00e676",
                     ifelse(df$spleen_pct < 0, "#7c4dff", "#ff5252"))
    plot_ly(df, x = ~label, y = ~spleen_pct, type = "bar",
            marker = list(color = colors)) %>%
      add_lines(x = c(-0.5, nrow(df)-0.5), y = c(-35, -35),
                line = list(color = "#ffab40", dash = "dash", width = 2),
                name = "SVR35 threshold") %>%
      layout(
        xaxis = list(title = "", color = "#aaa", tickangle = -30,
                     categoryorder = "array", categoryarray = df$label),
        yaxis = list(title = "Spleen Volume Change (%)",
                     color = "#aaa", gridcolor = "#2d3561"),
        paper_bgcolor = "#1a1d2e", plot_bgcolor = "#1a1d2e",
        font = list(color = "#e8eaf6"), showlegend = FALSE,
        margin = list(l=60, r=20, t=20, b=80)
      )
  })

  output$forest_plot <- renderPlotly({
    df <- comp_rv()
    df <- df[df$label != "No Treatment", ]
    # Odds ratios vs No-treatment reference (approximate)
    df$or_svr  <- (df$svr35 / (100 - df$svr35)) / (0.05 / 0.95)
    df$or_lo   <- df$or_svr * 0.65
    df$or_hi   <- df$or_svr * 1.55
    df$y_pos   <- seq_len(nrow(df))

    plot_ly() %>%
      add_trace(
        x = df$or_svr, y = df$y_pos, type = "scatter", mode = "markers",
        marker = list(size = 12, color = "#7c4dff", symbol = "diamond"),
        error_x = list(
          type = "data",
          symmetric = FALSE,
          array    = df$or_hi - df$or_svr,
          arrayminus = df$or_svr - df$or_lo,
          color = "#7c4dff"
        ),
        hovertemplate = paste0(df$label, "<br>OR: %{x:.2f}")
      ) %>%
      add_lines(x = c(1, 1), y = c(0.5, nrow(df)+0.5),
                line = list(color = "#ff5252", dash = "dash", width = 1.5)) %>%
      layout(
        xaxis = list(title = "Odds Ratio for SVR35 (vs No Treatment)", type = "log",
                     color = "#aaa", gridcolor = "#2d3561"),
        yaxis = list(title = "", tickvals = df$y_pos, ticktext = df$label,
                     color = "#aaa", gridcolor = "#2d3561"),
        paper_bgcolor = "#1a1d2e", plot_bgcolor = "#1a1d2e",
        font = list(color = "#e8eaf6"), showlegend = FALSE,
        margin = list(l=170, r=20, t=20, b=50)
      )
  })

  # ===========================================================
  # TAB 6: Biomarker Dynamics
  # ===========================================================
  bm_rv <- eventReactive(input$run_bm, {
    req(input$bm_drug, input$bm_dose, input$bm_int, input$bm_weeks)
    set.seed(999)
    list(
      pd  = sim_pd_biomarkers(input$bm_drug, input$bm_dose, input$bm_int, input$bm_weeks),
      cln = sim_clinical(input$bm_drug, input$bm_dose, input$bm_int,
                         input$svol_bl, input$hgb_bl, input$plt_bl,
                         input$tss_bl, min(input$bm_weeks, 52))
    )
  }, ignoreNULL = FALSE)

  output$vaf_svol_corr <- renderPlotly({
    d <- bm_rv()
    # Match weeks between pd and cln
    wks <- intersect(d$pd$week, d$cln$week)
    pd2 <- d$pd[d$pd$week %in% wks, ]
    cl2 <- d$cln[d$cln$week %in% wks, ]

    svol_pct_chg <- (cl2$svol - input$svol_bl) / input$svol_bl * 100

    plot_ly(x = pd2$vaf * 100, y = svol_pct_chg, type = "scatter",
            mode = "markers",
            marker = list(
              size  = 8,
              color = pd2$week,
              colorscale = "Viridis",
              showscale = TRUE,
              colorbar = list(title = "Week")
            ),
            hovertemplate = "VAF: %{x:.1f}%<br>Spleen ΔVol: %{y:.1f}%") %>%
      add_lines(x = range(pd2$vaf * 100), y = c(-35, -35),
                line = list(color = "#00e676", dash = "dash"),
                name = "SVR35") %>%
      layout(
        xaxis = list(title = "JAK2 VAF (%)", color = "#aaa", gridcolor = "#2d3561"),
        yaxis = list(title = "Spleen Volume Change (%)",
                     color = "#aaa", gridcolor = "#2d3561"),
        paper_bgcolor = "#1a1d2e", plot_bgcolor = "#1a1d2e",
        font = list(color = "#e8eaf6"), showlegend = FALSE,
        margin = list(l=60, r=80, t=20, b=50)
      )
  })

  output$cytokine_heat <- renderPlotly({
    drugs <- c("Ruxolitinib", "Fedratinib", "Pacritinib", "Momelotinib")
    timepoints <- c(0, 4, 12, 24, 52)
    cytokines  <- c("IL-6", "TNF-α", "IL-1β", "TGF-β", "CXCL9", "EPO")

    set.seed(111)
    # Matrix of cytokine levels (scaled 0-1)
    mat <- matrix(NA, nrow = length(drugs), ncol = length(cytokines) * length(timepoints))
    heat_x <- c(); heat_y <- c(); heat_z <- c()

    for (di in seq_along(drugs)) {
      pd_d <- sim_pd_biomarkers(drugs[di], pk_params[[drugs[di]]]$ref_dose,
                                pk_params[[drugs[di]]]$ref_interval, 52)
      for (ci in seq_along(cytokines)) {
        for (ti in seq_along(timepoints)) {
          idx <- which.min(abs(pd_d$week - timepoints[ti]))
          base_val <- switch(cytokines[ci],
            "IL-6"  = pd_d$il6[idx],
            "TNF-α" = pd_d$tnfa[idx],
            "IL-1β" = 18 * exp(-0.015 * timepoints[ti] * pd_d$pstat3_inh[idx] / 100),
            "TGF-β" = 25 + 0.3 * timepoints[ti],
            "CXCL9" = 120 * exp(-0.02 * timepoints[ti]),
            "EPO"   = 30 * (1 - 0.05 * min(timepoints[ti], 12))
          )
          heat_x <- c(heat_x, paste0(cytokines[ci], "\nW", timepoints[ti]))
          heat_y <- c(heat_y, drugs[di])
          heat_z <- c(heat_z, base_val)
        }
      }
    }

    z_mat <- matrix(heat_z, nrow = length(drugs),
                    ncol = length(cytokines) * length(timepoints), byrow = FALSE)
    # Normalize each column (cytokine)
    for (j in seq_len(ncol(z_mat)))
      z_mat[, j] <- (z_mat[, j] - min(z_mat[, j])) /
                    (max(z_mat[, j]) - min(z_mat[, j]) + 1e-6)

    x_labels <- unique(heat_x)
    plot_ly(z = z_mat, x = x_labels, y = drugs, type = "heatmap",
            colorscale = "RdBu", reversescale = TRUE,
            hovertemplate = "Drug: %{y}<br>%{x}<br>Rel. Level: %{z:.2f}") %>%
      layout(
        xaxis = list(title = "", color = "#aaa", tickangle = -45, tickfont = list(size = 9)),
        yaxis = list(title = "", color = "#aaa"),
        paper_bgcolor = "#1a1d2e", plot_bgcolor = "#1a1d2e",
        font = list(color = "#e8eaf6"),
        margin = list(l=110, r=20, t=20, b=120)
      )
  })

  output$pstat_dyn <- renderPlotly({
    d <- bm_rv()
    df <- d$pd
    plot_ly(df) %>%
      add_trace(x = ~week, y = ~pstat3_inh, type = "scatter", mode = "lines",
                line = list(color = "#7c4dff", width = 2),
                name = "pSTAT3", fill = "tozeroy",
                fillcolor = "rgba(124,77,255,0.12)") %>%
      add_trace(x = ~week, y = ~pstat5_inh, type = "scatter", mode = "lines",
                line = list(color = "#00e676", width = 2),
                name = "pSTAT5", fill = "tozeroy",
                fillcolor = "rgba(0,230,118,0.12)") %>%
      layout(
        xaxis = list(title = "Week", color = "#aaa", gridcolor = "#2d3561"),
        yaxis = list(title = "Inhibition (%)", range = c(0, 100),
                     color = "#aaa", gridcolor = "#2d3561"),
        paper_bgcolor = "#1a1d2e", plot_bgcolor = "#1a1d2e",
        font = list(color = "#e8eaf6"), legend = list(x = 0.7, y = 0.15),
        margin = list(l=60, r=20, t=20, b=50)
      )
  })

  output$km_aml <- renderPlotly({
    # Kaplan-Meier-style curves for AML transformation by risk group
    times <- seq(0, 120, by = 1)  # months

    # Hazard rates (per month) for AML transformation by DIPSS Plus risk
    hz <- list(
      Low   = 0.001,
      "Int-1" = 0.004,
      "Int-2" = 0.009,
      High    = 0.018
    )
    colors <- c("#00e676", "#ffeb3b", "#ffab40", "#ff5252")
    names(colors) <- names(hz)

    current_risk <- dipss_rv()$risk

    pl <- plot_ly()
    for (i in seq_along(hz)) {
      risk_name <- names(hz)[i]
      surv <- exp(-hz[[i]] * times)
      lw <- if (risk_name == current_risk) 3.5 else 1.5
      pl <- pl %>% add_trace(
        x = times, y = surv * 100, type = "scatter", mode = "lines",
        name = paste(risk_name, "Risk"),
        line = list(color = colors[risk_name], width = lw),
        hovertemplate = paste0(risk_name, "<br>Month: %{x}<br>AML-free: %{y:.1f}%")
      )
    }

    # Highlight patient's risk
    pl %>% layout(
      xaxis = list(title = "Time (months)", color = "#aaa", gridcolor = "#2d3561"),
      yaxis = list(title = "AML-free Survival (%)", range = c(0, 102),
                   color = "#aaa", gridcolor = "#2d3561"),
      paper_bgcolor = "#1a1d2e", plot_bgcolor = "#1a1d2e",
      font = list(color = "#e8eaf6"),
      legend = list(x = 0.65, y = 0.95),
      annotations = list(list(
        x = 60, y = 50,
        text = paste("Your patient:", current_risk),
        showarrow = FALSE, font = list(color = "#ffab40", size = 12)
      )),
      margin = list(l=60, r=20, t=20, b=50)
    )
  })

}

# ============================================================
# LAUNCH
# ============================================================
shinyApp(ui = ui, server = server)
