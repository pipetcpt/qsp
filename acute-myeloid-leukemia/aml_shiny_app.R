# ============================================================
# AML QSP Shiny Application
# Acute Myeloid Leukemia — Quantitative Systems Pharmacology
# ============================================================
# Libraries: shiny, bslib, plotly, ggplot2, dplyr, scales
# 6 Tabs: Patient Profile | Drug PK | Leukemia Dynamics (PD) |
#         Clinical Endpoints | Scenario Comparison | Biomarker Dashboard
# ============================================================

library(shiny)
library(bslib)
library(plotly)
library(ggplot2)
library(dplyr)
library(scales)

# ============================================================
# Helper Functions
# ============================================================

#' simulate_aml: Core ODE-based AML simulation
#' @param params list of patient and treatment parameters
#' @param days numeric vector of time points
#' @return data.frame with longitudinal outputs
simulate_aml <- function(params, days = seq(0, 336, by = 1)) {

  # Unpack parameters
  age         <- params$age          %||% 65
  weight      <- params$weight       %||% 75
  wbc0        <- params$wbc0         %||% 30
  bm_blast0   <- params$bm_blast0    %||% 60
  pb_blast0   <- params$pb_blast0    %||% 40
  regimen     <- params$regimen      %||% "VEN_AZA"
  mutations   <- params$mutations    %||% character(0)
  cyto_risk   <- params$cyto_risk    %||% "Intermediate"
  flt3_itd    <- "FLT3-ITD" %in% mutations
  npm1        <- "NPM1" %in% mutations
  idh1        <- "IDH1" %in% mutations
  idh2        <- "IDH2" %in% mutations
  tp53        <- "TP53" %in% mutations

  # Risk modifier
  risk_mod <- switch(cyto_risk,
    "Favorable"    = 0.6,
    "Intermediate" = 1.0,
    "Adverse"      = 1.5,
    "Very Adverse" = 2.0,
    1.0
  )
  if (tp53) risk_mod <- risk_mod * 1.6
  if (flt3_itd && regimen != "Gilteritinib") risk_mod <- risk_mod * 1.3
  if (npm1 && regimen %in% c("VEN_AZA", "7plus3")) risk_mod <- risk_mod * 0.8

  # Regimen-specific parameters
  reg_params <- list(
    "7plus3"       = list(kill_rate = 0.055, maint_kill = 0.015, cr_prob = 0.70, relapse_rate = 0.0045),
    "VEN_AZA"      = list(kill_rate = 0.040, maint_kill = 0.020, cr_prob = 0.66, relapse_rate = 0.0030),
    "Gilteritinib" = list(kill_rate = 0.035, maint_kill = 0.022, cr_prob = 0.54, relapse_rate = 0.0040),
    "Enasidenib"   = list(kill_rate = 0.028, maint_kill = 0.018, cr_prob = 0.40, relapse_rate = 0.0035),
    "ATRA_ATO"     = list(kill_rate = 0.060, maint_kill = 0.030, cr_prob = 0.95, relapse_rate = 0.0010),
    "FLAG_IDA"     = list(kill_rate = 0.065, maint_kill = 0.020, cr_prob = 0.62, relapse_rate = 0.0038),
    "LDAC_VEN"     = list(kill_rate = 0.030, maint_kill = 0.015, cr_prob = 0.48, relapse_rate = 0.0038)
  )
  rp <- reg_params[[regimen]] %||% reg_params[["VEN_AZA"]]
  kr <- rp$kill_rate / risk_mod
  mk <- rp$maint_kill / risk_mod
  cr_prob <- rp$cr_prob / sqrt(risk_mod)

  # Blast trajectory (biphasic kill + relapse model)
  blast_traj <- numeric(length(days))
  lsc_traj   <- numeric(length(days))
  mrd_traj   <- numeric(length(days))

  # LSC initial: ~1% of total leukemic burden
  lsc0 <- bm_blast0 * 0.01
  # MRD initial (log copies/PCR units)
  mrd0 <- log10(bm_blast0 * 1e4)  # ~5-6 log copies at diagnosis

  # Myelosuppression model (Friberg-type)
  circ0_wbc <- 5.0   # normal WBC baseline (10^3/µL)
  circ0_plt <- 200   # normal PLT (10^9/L)
  circ0_hgb <- 13.5  # normal Hgb (g/dL)
  circ0_anc <- 2.0   # normal ANC

  wbc_traj  <- numeric(length(days))
  anc_traj  <- numeric(length(days))
  plt_traj  <- numeric(length(days))
  hgb_traj  <- numeric(length(days))
  pb_blast_traj <- numeric(length(days))

  # Determine cycle structure
  cycle_days <- switch(regimen,
    "7plus3"    = 28,
    "VEN_AZA"   = 28,
    "FLAG_IDA"  = 28,
    "ATRA_ATO"  = 28,
    "LDAC_VEN"  = 28,
    28
  )

  cr_achieved <- FALSE
  cr_day      <- NA
  relapse_day <- NA
  relapsed    <- FALSE

  for (i in seq_along(days)) {
    d <- days[i]
    cycle <- floor(d / cycle_days) + 1

    # Drug effect: decays with time (resistance buildup)
    resist_factor <- 1 + 0.002 * d
    eff_kr <- kr / resist_factor

    if (!relapsed) {
      # Blast dynamics
      if (d == 0) {
        blast_traj[i]    <- bm_blast0
        lsc_traj[i]      <- lsc0
        mrd_traj[i]      <- mrd0
        pb_blast_traj[i] <- pb_blast0
      } else {
        prev_blast <- blast_traj[i - 1]
        prev_lsc   <- lsc_traj[i - 1]
        prev_mrd   <- mrd_traj[i - 1]
        prev_pb    <- pb_blast_traj[i - 1]

        # Net blast reduction
        growth_rate <- 0.005  # slow regrowth
        net_kill    <- eff_kr - growth_rate

        new_blast <- prev_blast * exp(-net_kill)
        new_blast <- max(new_blast, 0.01)

        # LSC follows slower kinetics (more resistant)
        lsc_kill   <- eff_kr * 0.3  # LSC 70% less sensitive
        new_lsc    <- prev_lsc * exp(-lsc_kill + 0.003)
        new_lsc    <- max(new_lsc, 0.001)

        # MRD dynamics (log scale)
        mrd_kill   <- eff_kr * 0.6
        new_mrd    <- prev_mrd - mrd_kill * 0.5
        new_mrd    <- max(new_mrd, -1)

        # PB blast tracks BM blast with ~0.6 correlation
        new_pb <- prev_pb * exp(-net_kill * 1.1)
        new_pb <- max(new_pb, 0.01)

        # CR check
        if (!cr_achieved && new_blast < 5) {
          cr_achieved <- TRUE
          cr_day      <- d
        }

        # Relapse after CR (LSC-driven)
        if (cr_achieved && d > cr_day + 30) {
          relapse_prob_daily <- rp$relapse_rate * (new_lsc / lsc0)^0.3 * risk_mod
          if (runif(1) < relapse_prob_daily) {
            relapsed    <- TRUE
            relapse_day <- d
          }
        }

        blast_traj[i]    <- new_blast
        lsc_traj[i]      <- new_lsc
        mrd_traj[i]      <- new_mrd
        pb_blast_traj[i] <- new_pb
      }
    } else {
      # Post-relapse: blast regrowth
      days_since_relapse <- d - relapse_day
      blast_traj[i]    <- min(5 * exp(0.015 * days_since_relapse), 95)
      lsc_traj[i]      <- min(lsc_traj[i - 1] * 1.02, lsc0 * 5)
      mrd_traj[i]      <- min(mrd_traj[i - 1] + 0.02, mrd0)
      pb_blast_traj[i] <- min(5 * exp(0.018 * days_since_relapse), 90)
    }

    # Myelosuppression (Friberg model approximation)
    # Nadir at ~day 14 per cycle, recovery by day 28
    day_in_cycle <- d %% cycle_days
    myelosupp_factor <- if (day_in_cycle <= 14) {
      1 - 0.5 * sin(pi * day_in_cycle / 14)
    } else {
      0.5 + 0.5 * sin(pi * (day_in_cycle - 14) / 14)
    }
    myelosupp_factor <- max(myelosupp_factor, 0.1)

    # Blast suppression also drives cytopenias
    leukemia_effect <- (1 - blast_traj[i] / 100)

    wbc_traj[i]  <- circ0_wbc  * myelosupp_factor * (0.6 + 0.4 * leukemia_effect)
    anc_traj[i]  <- circ0_anc  * myelosupp_factor * (0.5 + 0.5 * leukemia_effect)
    plt_traj[i]  <- circ0_plt  * myelosupp_factor * (0.6 + 0.4 * leukemia_effect)
    hgb_traj[i]  <- circ0_hgb  * (0.7 + 0.3 * myelosupp_factor) * (0.8 + 0.2 * leukemia_effect)
  }

  data.frame(
    day         = days,
    bm_blast    = blast_traj,
    pb_blast    = pb_blast_traj,
    lsc         = lsc_traj,
    mrd         = mrd_traj,
    wbc         = wbc_traj,
    anc         = anc_traj,
    plt         = plt_traj,
    hgb         = hgb_traj,
    cr_achieved = cr_achieved,
    cr_day      = cr_day %||% NA_real_,
    relapse_day = relapse_day %||% NA_real_
  )
}

# Null coalescing operator
`%||%` <- function(a, b) if (!is.null(a) && length(a) > 0 && !is.na(a[1])) a else b

# ============================================================
# PK Parameter Databases
# ============================================================

pk_params_db <- list(
  Venetoclax = list(
    dose = 400, unit = "mg", route = "oral",
    ka = 0.8, CL = 8.0, Vd = 256, t_half = 26,
    bioav = 0.56,
    Cmax_fasted = 1.8, Cmax_fed = 5.1,
    Tmax = 5, AUC_fasted = 32, AUC_fed = 94,
    cyp3a4_fold = 5.0,
    color = "#E74C3C",
    note = "Take with food. Strong CYP3A4 inhibitors (posaconazole) increase AUC ~5-fold — reduce dose to 10mg."
  ),
  Azacitidine = list(
    dose = 75, unit = "mg/m2", route = "SC",
    ka = 2.0, CL = 89, Vd = 76, t_half = 0.7,
    bioav = 0.89, Cmax_fasted = 750, Cmax_fed = 750,
    Tmax = 0.5, AUC_fasted = 437, AUC_fed = 437,
    cyp3a4_fold = 1.0,
    color = "#3498DB",
    note = "SC or IV administration. Short half-life; AML exposure driven by intracellular triphosphate."
  ),
  Gilteritinib = list(
    dose = 120, unit = "mg", route = "oral",
    ka = 0.5, CL = 8.6, Vd = 1092, t_half = 113,
    bioav = 0.90, Cmax_fasted = 282, Cmax_fed = 282,
    Tmax = 6, AUC_fasted = 31500, AUC_fed = 31500,
    cyp3a4_fold = 1.5,
    color = "#2ECC71",
    note = "FLT3 inhibitor. Long t1/2 allows once-daily dosing. Steady state reached ~15 days."
  ),
  Enasidenib = list(
    dose = 100, unit = "mg", route = "oral",
    ka = 0.7, CL = 0.74, Vd = 55.8, t_half = 137,
    bioav = 0.57, Cmax_fasted = 1410, Cmax_fed = 1410,
    Tmax = 4, AUC_fasted = 246000, AUC_fed = 246000,
    cyp3a4_fold = 1.2,
    color = "#9B59B6",
    note = "IDH2 inhibitor. Very long t1/2 (137h). Differentiation syndrome risk weeks 1-3."
  ),
  Ivosidenib = list(
    dose = 500, unit = "mg", route = "oral",
    ka = 0.6, CL = 3.6, Vd = 183, t_half = 93,
    bioav = 0.98, Cmax_fasted = 3430, Cmax_fed = 3430,
    Tmax = 3, AUC_fasted = 380000, AUC_fed = 380000,
    cyp3a4_fold = 1.3,
    color = "#F39C12",
    note = "IDH1 inhibitor. QTc prolongation monitoring required."
  ),
  Cytarabine = list(
    dose = 200, unit = "mg/m2", route = "IV_CI",
    ka = 10, CL = 72, Vd = 31, t_half = 2.1,
    bioav = 1.0, Cmax_fasted = 1900, Cmax_fed = 1900,
    Tmax = 1, AUC_fasted = 3800, AUC_fed = 3800,
    cyp3a4_fold = 1.0,
    color = "#1ABC9C",
    note = "IV continuous infusion 7 days. Active metabolite ara-CTP has t1/2 ~3.4h intracellularly."
  ),
  Idarubicin = list(
    dose = 12, unit = "mg/m2", route = "IV",
    ka = 10, CL = 58, Vd = 1770, t_half = 22,
    bioav = 1.0, Cmax_fasted = 97, Cmax_fed = 97,
    Tmax = 0.25, AUC_fasted = 2200, AUC_fed = 2200,
    cyp3a4_fold = 1.0,
    color = "#E67E22",
    note = "Anthracycline. Cardiotoxicity monitoring required. Active metabolite idarubicinol."
  )
)

# Simulate PK concentration-time profile
simulate_pk <- function(drug_name, doses = 7, interval = 24, time_end = NULL) {
  pk <- pk_params_db[[drug_name]]
  if (is.null(pk)) return(data.frame())

  if (is.null(time_end)) {
    time_end <- max(doses * interval + 4 * pk$t_half, 168)
  }
  times <- seq(0, time_end, by = 0.25)

  # Multi-dose 1-compartment model
  conc <- numeric(length(times))
  dose_mg <- pk$dose * pk$bioav

  for (j in seq_along(times)) {
    t <- times[j]
    c_sum <- 0
    for (k in 0:(doses - 1)) {
      t_post <- t - k * interval
      if (t_post > 0) {
        # 1-compartment oral: ka and ke
        ke <- log(2) / pk$t_half
        if (pk$route == "IV" || pk$route == "IV_CI") {
          c_sum <- c_sum + (dose_mg / pk$Vd) * exp(-ke * t_post)
        } else {
          ka <- pk$ka
          c_sum <- c_sum + (dose_mg * ka / (pk$Vd * (ka - ke))) *
            (exp(-ke * t_post) - exp(-ka * t_post))
        }
      }
    }
    conc[j] <- max(c_sum, 0)
  }

  data.frame(time = times, conc = conc, drug = drug_name)
}

# ============================================================
# Survival / OS Simulation (parametric log-normal)
# ============================================================

simulate_os <- function(regimen, n_months = 36, risk_group = "Intermediate") {
  # Median OS by regimen and risk from key trials
  median_os_table <- list(
    "7plus3"       = list(Favorable = 60, Intermediate = 18, Adverse = 9,  "Very Adverse" = 5),
    "VEN_AZA"      = list(Favorable = NA, Intermediate = 14.7, Adverse = 11, "Very Adverse" = 5),
    "Gilteritinib" = list(Favorable = NA, Intermediate = 9.3, Adverse = 7.5, "Very Adverse" = 4),
    "Enasidenib"   = list(Favorable = NA, Intermediate = 8.8, Adverse = 6.5, "Very Adverse" = 3),
    "ATRA_ATO"     = list(Favorable = NA, Intermediate = 99,  Adverse = NA,  "Very Adverse" = NA),
    "FLAG_IDA"     = list(Favorable = NA, Intermediate = 12,  Adverse = 8,   "Very Adverse" = 4),
    "LDAC_VEN"     = list(Favorable = NA, Intermediate = 10.1, Adverse = 7,  "Very Adverse" = 3.5)
  )

  med_os <- median_os_table[[regimen]][[risk_group]]
  if (is.null(med_os) || is.na(med_os)) med_os <- 8

  # Log-normal survival
  months <- seq(0, n_months, by = 0.5)
  mu  <- log(med_os)
  sig <- 0.8
  surv <- 1 - plnorm(months, meanlog = mu, sdlog = sig)

  # EFS ~ 70% of OS
  efs_med <- med_os * 0.65
  efs_surv <- 1 - plnorm(months, meanlog = log(efs_med), sdlog = 0.9)

  # Relapse incidence (cause-specific)
  relapse_incidence <- cumsum(diff(c(1, efs_surv)) * -0.7) * 100
  nrm_incidence     <- cumsum(diff(c(1, efs_surv)) * -0.3) * 100

  data.frame(
    month     = months,
    os        = pmax(surv, 0) * 100,
    efs       = pmax(efs_surv, 0) * 100,
    relapse   = c(0, relapse_incidence),
    nrm       = c(0, nrm_incidence),
    regimen   = regimen
  )
}

# ============================================================
# ELN 2022 Risk Classification
# ============================================================

classify_eln_2022 <- function(mutations, cyto_risk, subtype) {
  muts <- mutations

  # Favorable
  if (cyto_risk == "Favorable" ||
      ("NPM1" %in% muts && !("FLT3-ITD" %in% muts) && !"TP53" %in% muts)) {
    risk  <- "Favorable"
    cr    <- "85-95%"
    os2yr <- "~60%"
    mrd_neg <- "70-80%"
    rec   <- "Standard induction (7+3). Consider allo-SCT only if MRD+."
    color <- "success"
  } else if ("TP53" %in% muts || "RUNX1" %in% muts || cyto_risk == "Very Adverse") {
    risk  <- "Adverse"
    cr    <- "40-55%"
    os2yr <- "~10-15%"
    mrd_neg <- "20-30%"
    rec   <- "Clinical trial preferred. Allo-SCT in CR1. Consider decitabine/cedazuridine + VEN."
    color <- "danger"
  } else if ("FLT3-ITD" %in% muts || cyto_risk == "Adverse") {
    risk  <- "Adverse"
    cr    <- "50-65%"
    os2yr <- "~15-25%"
    mrd_neg <- "30-45%"
    rec   <- "FLT3 inhibitor (midostaurin or gilteritinib) + 7+3. Allo-SCT in CR1."
    color <- "danger"
  } else {
    risk  <- "Intermediate"
    cr    <- "65-75%"
    os2yr <- "~30-40%"
    mrd_neg <- "40-60%"
    rec   <- "Standard induction (7+3) ± targeted agent. Allo-SCT in CR1 if donor available."
    color <- "warning"
  }

  list(risk = risk, cr = cr, os2yr = os2yr, mrd_neg = mrd_neg, rec = rec, color = color)
}

# ============================================================
# UI Definition
# ============================================================

ui <- page_navbar(
  title = tags$span(
    tags$img(src = NULL, height = "30px"),
    "AML QSP Dashboard"
  ),
  theme = bs_theme(
    version = 5,
    bootswatch = "flatly",
    primary   = "#2C3E50",
    secondary = "#E74C3C",
    success   = "#27AE60",
    info      = "#2980B9",
    warning   = "#F39C12",
    danger    = "#E74C3C"
  ),
  bg = "#2C3E50",
  inverse = TRUE,

  # --------------------------------------------------------
  # TAB 1: Patient Profile
  # --------------------------------------------------------
  nav_panel(
    title = "Patient Profile",
    icon  = icon("user"),
    layout_sidebar(
      sidebar = sidebar(
        width = 320,
        bg = "#ECF0F1",
        h5("Demographics", class = "text-primary fw-bold"),
        sliderInput("age",    "Age (years)",          min = 18,  max = 90,  value = 65, step = 1),
        sliderInput("weight", "Body weight (kg)",     min = 40,  max = 150, value = 75, step = 1),
        hr(),
        h5("Presentation", class = "text-primary fw-bold"),
        numericInput("wbc_dx",     "WBC at diagnosis (×10³/µL)", value = 30,  min = 1,   max = 400),
        numericInput("bm_blast",   "BM blasts (%)",               value = 60,  min = 0,   max = 100),
        numericInput("pb_blast",   "Peripheral blasts (%)",       value = 40,  min = 0,   max = 100),
        hr(),
        h5("Classification", class = "text-primary fw-bold"),
        selectInput("aml_subtype", "AML Subtype",
          choices = c("De novo AML", "Secondary AML (sAML)", "Therapy-related AML (tAML)", "APL (t(15;17))")
        ),
        selectInput("cyto_risk", "Cytogenetics Risk Group",
          choices = c("Favorable", "Intermediate", "Adverse", "Very Adverse")
        ),
        hr(),
        h5("Molecular Mutations", class = "text-primary fw-bold"),
        checkboxGroupInput("mutations", NULL,
          choices  = c("FLT3-ITD", "FLT3-TKD", "NPM1", "DNMT3A", "IDH1", "IDH2",
                       "TP53", "CEBPA", "RUNX1", "ASXL1"),
          selected = c("FLT3-ITD", "NPM1")
        ),
        hr(),
        h5("Treatment Intent", class = "text-primary fw-bold"),
        radioButtons("tx_intent", NULL,
          choices  = c("Intensive chemotherapy", "Non-intensive (VEN-based)", "Clinical trial"),
          selected = "Intensive chemotherapy"
        ),
        hr(),
        h5("Comorbidities", class = "text-primary fw-bold"),
        checkboxGroupInput("comorbidities", NULL,
          choices = c("Cardiac disease (EF <50%)", "Prior MDS/MPN", "Prior chemotherapy/RT",
                      "Renal impairment (CrCl <45)", "Hepatic impairment")
        ),
        hr(),
        actionButton("classify_eln", "Classify ELN 2022 Risk",
          class = "btn-primary w-100", icon = icon("dna"))
      ),

      # Main Panel
      layout_columns(
        col_widths = c(12),
        card(
          card_header("ELN 2022 Risk Classification & Treatment Recommendation",
            class = "bg-primary text-white"),
          card_body(
            uiOutput("eln_result_ui")
          )
        )
      ),
      layout_columns(
        col_widths = c(3, 3, 3, 3),
        card(
          card_body(
            h6("ELN 2022 Risk Group", class = "text-muted mb-1"),
            uiOutput("eln_badge_ui"),
            class = "text-center"
          )
        ),
        card(
          card_body(
            h6("Estimated CR Rate", class = "text-muted mb-1"),
            uiOutput("cr_rate_ui"),
            class = "text-center"
          )
        ),
        card(
          card_body(
            h6("2-Year OS Estimate", class = "text-muted mb-1"),
            uiOutput("os2yr_ui"),
            class = "text-center"
          )
        ),
        card(
          card_body(
            h6("MRD-Negativity Rate", class = "text-muted mb-1"),
            uiOutput("mrd_neg_ui"),
            class = "text-center"
          )
        )
      )
    )
  ),

  # --------------------------------------------------------
  # TAB 2: Drug PK
  # --------------------------------------------------------
  nav_panel(
    title = "Drug PK",
    icon  = icon("pills"),
    layout_sidebar(
      sidebar = sidebar(
        width = 280,
        bg = "#ECF0F1",
        h5("Drug Selection", class = "text-primary fw-bold"),
        checkboxGroupInput("pk_drugs", "Select drug(s) to display:",
          choices  = names(pk_params_db),
          selected = c("Venetoclax", "Azacitidine")
        ),
        hr(),
        sliderInput("pk_doses",    "Number of doses",        min = 1, max = 28, value = 7),
        sliderInput("pk_interval", "Dosing interval (h)",    min = 8, max = 48, value = 24, step = 8),
        hr(),
        h5("Venetoclax Interaction", class = "text-primary fw-bold"),
        checkboxInput("ven_cyp_inhibitor", "Add CYP3A4 inhibitor (posaconazole)", value = FALSE),
        conditionalPanel(
          condition = "input.ven_cyp_inhibitor == true",
          div(class = "alert alert-warning p-2",
            icon("exclamation-triangle"),
            " Posaconazole increases venetoclax AUC ~5-fold.",
            br(), strong("Reduce venetoclax dose to 10 mg/day."),
            style = "font-size:0.8em"
          )
        ),
        hr(),
        h5("Venetoclax Food Effect", class = "text-primary fw-bold"),
        checkboxInput("ven_food_effect", "Show fed vs. fasted comparison", value = TRUE)
      ),
      layout_columns(
        col_widths = c(12),
        card(
          card_header("Concentration–Time Profiles (Multiple Doses)", class = "bg-info text-white"),
          card_body(plotlyOutput("pk_conc_plot", height = "400px"))
        )
      ),
      layout_columns(
        col_widths = c(6, 6),
        card(
          card_header("PK Parameter Summary Table"),
          card_body(tableOutput("pk_param_table"))
        ),
        card(
          card_header("Drug Interaction & Food Effect Notes"),
          card_body(uiOutput("pk_notes_ui"))
        )
      )
    )
  ),

  # --------------------------------------------------------
  # TAB 3: Leukemia Dynamics (PD)
  # --------------------------------------------------------
  nav_panel(
    title = "Leukemia Dynamics (PD)",
    icon  = icon("chart-line"),
    layout_sidebar(
      sidebar = sidebar(
        width = 280,
        bg = "#ECF0F1",
        h5("Treatment Regimen", class = "text-primary fw-bold"),
        selectInput("pd_regimen", "Select regimen:",
          choices = c(
            "7+3 Induction" = "7plus3",
            "Venetoclax + Azacitidine (VIALE-A)" = "VEN_AZA",
            "Gilteritinib 120 mg/day (FLT3+)" = "Gilteritinib",
            "Enasidenib 100 mg/day (IDH2+)" = "Enasidenib",
            "ATRA + ATO (APL)" = "ATRA_ATO",
            "FLAG-IDA (Salvage)" = "FLAG_IDA",
            "LDAC + Venetoclax" = "LDAC_VEN"
          ),
          selected = "VEN_AZA"
        ),
        hr(),
        sliderInput("pd_days", "Simulation duration (days)", min = 56, max = 500, value = 336, step = 28),
        hr(),
        h5("Patient (from Tab 1)", class = "text-muted"),
        p("Age, mutations, and cytogenetics are carried forward from the Patient Profile tab.",
          style = "font-size:0.8em; color:#7F8C8D"),
        hr(),
        actionButton("run_pd_sim", "Run Simulation",
          class = "btn-success w-100", icon = icon("play")),
        hr(),
        uiOutput("pd_cr_info_ui")
      ),
      layout_columns(
        col_widths = c(6, 6),
        card(
          card_header("BM Blast % Trajectory", class = "bg-danger text-white"),
          card_body(plotlyOutput("pd_blast_plot", height = "300px"))
        ),
        card(
          card_header("Leukemic Stem Cell & MRD (log10)", class = "bg-warning text-white"),
          card_body(plotlyOutput("pd_lsc_mrd_plot", height = "300px"))
        )
      ),
      layout_columns(
        col_widths = c(12),
        card(
          card_header("CBC: WBC, ANC, Plt, Hgb (Myelosuppression Model)", class = "bg-primary text-white"),
          card_body(plotlyOutput("pd_cbc_plot", height = "320px"))
        )
      )
    )
  ),

  # --------------------------------------------------------
  # TAB 4: Clinical Endpoints
  # --------------------------------------------------------
  nav_panel(
    title = "Clinical Endpoints",
    icon  = icon("hospital"),
    layout_sidebar(
      sidebar = sidebar(
        width = 280,
        bg = "#ECF0F1",
        h5("Simulation Settings", class = "text-primary fw-bold"),
        sliderInput("ep_months", "Follow-up (months)", min = 12, max = 60, value = 36, step = 6),
        selectInput("ep_risk", "Risk Group",
          choices  = c("Favorable", "Intermediate", "Adverse", "Very Adverse"),
          selected = "Intermediate"
        ),
        checkboxGroupInput("ep_regimens", "Compare regimens:",
          choices  = c(
            "7+3 Induction" = "7plus3",
            "VEN + AZA" = "VEN_AZA",
            "Gilteritinib" = "Gilteritinib",
            "Enasidenib" = "Enasidenib",
            "ATRA + ATO" = "ATRA_ATO"
          ),
          selected = c("7plus3", "VEN_AZA", "Gilteritinib")
        ),
        hr(),
        p("Based on: VIALE-A (NCT02993523), ADMIRAL (NCT02421939), QuANTUM-R, RATIFY trials.",
          style = "font-size:0.75em; color:#7F8C8D")
      ),
      layout_columns(
        col_widths = c(6, 6),
        card(
          card_header("Overall Survival", class = "bg-primary text-white"),
          card_body(plotlyOutput("ep_os_plot", height = "350px"))
        ),
        card(
          card_header("Event-Free Survival", class = "bg-info text-white"),
          card_body(plotlyOutput("ep_efs_plot", height = "350px"))
        )
      ),
      layout_columns(
        col_widths = c(6, 6),
        card(
          card_header("Cumulative Incidence: Relapse vs NRM"),
          card_body(plotlyOutput("ep_cuminc_plot", height = "300px"))
        ),
        card(
          card_header("Clinical Outcomes Summary (Key Trials)"),
          card_body(
            div(style = "overflow-x:auto;",
              tableOutput("ep_summary_table")
            )
          )
        )
      )
    )
  ),

  # --------------------------------------------------------
  # TAB 5: Scenario Comparison
  # --------------------------------------------------------
  nav_panel(
    title = "Scenario Comparison",
    icon  = icon("balance-scale"),
    layout_sidebar(
      sidebar = sidebar(
        width = 280,
        bg = "#ECF0F1",
        h5("Comparison Settings", class = "text-primary fw-bold"),
        radioButtons("sc_view", "Chart type:",
          choices  = c("Bar chart (blast reduction)",
                       "Radar chart (multi-domain)",
                       "Table (trial data)"),
          selected = "Bar chart (blast reduction)"
        ),
        hr(),
        sliderInput("sc_cycle", "Assessment at cycle:", min = 1, max = 12, value = 2),
        hr(),
        h5("7 Regimens Compared:", class = "text-muted"),
        tags$ol(
          tags$li("7+3 Induction (Ara-C + Idarubicin)"),
          tags$li("VEN + AZA (VIALE-A)"),
          tags$li("Gilteritinib 120 mg (FLT3+ R/R)"),
          tags$li("Enasidenib 100 mg (IDH2+ R/R)"),
          tags$li("LDAC + VEN (low-intensity)"),
          tags$li("FLAG-IDA (intensive salvage)"),
          tags$li("ATRA + ATO (APL standard)")
        )
      ),
      layout_columns(
        col_widths = c(12),
        card(
          card_header("Regimen Comparison — Blast Reduction & Outcomes", class = "bg-primary text-white"),
          card_body(uiOutput("sc_main_ui"))
        )
      ),
      layout_columns(
        col_widths = c(12),
        card(
          card_header("Clinical Trial Evidence Summary"),
          card_body(
            div(style = "overflow-x:auto;",
              tableOutput("sc_trial_table")
            )
          )
        )
      )
    )
  ),

  # --------------------------------------------------------
  # TAB 6: Biomarker Dashboard
  # --------------------------------------------------------
  nav_panel(
    title = "Biomarker Dashboard",
    icon  = icon("vials"),
    layout_sidebar(
      sidebar = sidebar(
        width = 280,
        bg = "#ECF0F1",
        h5("Biomarker Monitoring", class = "text-primary fw-bold"),
        selectInput("bm_regimen", "Treatment regimen:",
          choices = c(
            "VEN + AZA" = "VEN_AZA",
            "7+3" = "7plus3",
            "Gilteritinib" = "Gilteritinib",
            "Enasidenib" = "Enasidenib"
          )
        ),
        hr(),
        checkboxGroupInput("bm_show", "Show biomarkers:",
          choices  = c("NPM1-MRD (PCR)", "FLT3-ITD VAF", "IDH mutation (2-HG proxy)",
                       "WBC/ANC/PLT/Hgb", "CRP / Ferritin", "BM Cellularity"),
          selected = c("NPM1-MRD (PCR)", "FLT3-ITD VAF", "WBC/ANC/PLT/Hgb")
        ),
        hr(),
        sliderInput("bm_days", "Assessment days (current):", min = 28, max = 336, value = 84, step = 28),
        hr(),
        actionButton("run_bm", "Update Biomarkers", class = "btn-success w-100", icon = icon("sync"))
      ),
      layout_columns(
        col_widths = c(3, 3, 3, 3),
        card(
          card_body(
            class = "text-center",
            h6("MRD Status", class = "text-muted"),
            uiOutput("bm_mrd_status_box")
          )
        ),
        card(
          card_body(
            class = "text-center",
            h6("Days to MRD Negativity", class = "text-muted"),
            uiOutput("bm_days_to_mrd_box")
          )
        ),
        card(
          card_body(
            class = "text-center",
            h6("12-Month Relapse Risk", class = "text-muted"),
            uiOutput("bm_relapse_risk_box")
          )
        ),
        card(
          card_body(
            class = "text-center",
            h6("Next Assessment", class = "text-muted"),
            uiOutput("bm_next_assessment_box")
          )
        )
      ),
      layout_columns(
        col_widths = c(6, 6),
        card(
          card_header("MRD Trajectory (log10 copies)", class = "bg-danger text-white"),
          card_body(plotlyOutput("bm_mrd_plot", height = "280px"))
        ),
        card(
          card_header("FLT3-ITD Variant Allele Frequency (% VAF)", class = "bg-warning text-white"),
          card_body(plotlyOutput("bm_flt3_plot", height = "280px"))
        )
      ),
      layout_columns(
        col_widths = c(6, 6),
        card(
          card_header("CBC Trends with Normal Range Bands"),
          card_body(plotlyOutput("bm_cbc_plot", height = "280px"))
        ),
        card(
          card_header("Inflammatory Markers & BM Cellularity"),
          card_body(plotlyOutput("bm_inflam_plot", height = "280px"))
        )
      )
    )
  )
)

# ============================================================
# Server Logic
# ============================================================

server <- function(input, output, session) {

  # -----------------------------------------------------------
  # REACTIVE: AML simulation (shared across tabs)
  # -----------------------------------------------------------
  aml_sim <- reactive({
    input$run_pd_sim
    isolate({
      params <- list(
        age      = input$age,
        weight   = input$weight,
        wbc0     = input$wbc_dx,
        bm_blast0 = input$bm_blast,
        pb_blast0 = input$pb_blast,
        regimen  = input$pd_regimen,
        mutations = input$mutations,
        cyto_risk = input$cyto_risk
      )
      days <- seq(0, input$pd_days, by = 1)
      simulate_aml(params, days)
    })
  })

  bm_sim <- reactive({
    input$run_bm
    isolate({
      params <- list(
        age       = input$age,
        weight    = input$weight,
        wbc0      = input$wbc_dx,
        bm_blast0 = input$bm_blast,
        pb_blast0 = input$pb_blast,
        regimen   = input$bm_regimen,
        mutations = input$mutations,
        cyto_risk = input$cyto_risk
      )
      simulate_aml(params, seq(0, 336, by = 1))
    })
  })

  # -----------------------------------------------------------
  # TAB 1: ELN Classification
  # -----------------------------------------------------------

  eln_result <- eventReactive(input$classify_eln, {
    classify_eln_2022(input$mutations, input$cyto_risk, input$aml_subtype)
  }, ignoreNULL = FALSE)

  output$eln_result_ui <- renderUI({
    r <- eln_result()
    div(class = paste0("alert alert-", r$color),
      h5(paste("ELN 2022 Risk Group:", r$risk), class = "fw-bold"),
      p(strong("Recommended Treatment: "), r$rec),
      p(
        strong("CR Rate: "), r$cr, " | ",
        strong("2-year OS: "), r$os2yr, " | ",
        strong("MRD Negativity: "), r$mrd_neg
      ),
      hr(),
      p(style = "font-size:0.8em; margin-bottom:0",
        "Based on: Döhner H et al. Blood. 2022;140(12):1345-1377. PMID: 35021017")
    )
  })

  output$eln_badge_ui <- renderUI({
    r <- eln_result()
    span(class = paste0("badge bg-", r$color, " fs-4"), r$risk)
  })

  output$cr_rate_ui <- renderUI({
    r <- eln_result()
    h3(r$cr, class = paste0("text-", r$color))
  })

  output$os2yr_ui <- renderUI({
    r <- eln_result()
    h3(r$os2yr, class = paste0("text-", r$color))
  })

  output$mrd_neg_ui <- renderUI({
    r <- eln_result()
    h3(r$mrd_neg, class = paste0("text-", r$color))
  })

  # -----------------------------------------------------------
  # TAB 2: Drug PK
  # -----------------------------------------------------------

  output$pk_conc_plot <- renderPlotly({
    req(length(input$pk_drugs) > 0)
    all_pk <- lapply(input$pk_drugs, function(drug) {
      df <- simulate_pk(drug, doses = input$pk_doses, interval = input$pk_interval)

      # Venetoclax food effect
      if (drug == "Venetoclax" && input$ven_food_effect) {
        df_fasted <- df
        df_fasted$conc <- df_fasted$conc * (1810 / 5100)  # fasted ratio
        df_fasted$drug <- "Venetoclax (Fasted)"
        df$drug <- "Venetoclax (Fed)"
        df <- bind_rows(df, df_fasted)
      }

      # CYP3A4 inhibitor effect
      if (drug == "Venetoclax" && input$ven_cyp_inhibitor) {
        df$conc <- df$conc * 5.0
      }
      df
    })
    combined <- bind_rows(all_pk)

    drug_colors <- c(
      Venetoclax        = "#E74C3C",
      "Venetoclax (Fed)"    = "#E74C3C",
      "Venetoclax (Fasted)" = "#F1948A",
      Azacitidine       = "#3498DB",
      Gilteritinib      = "#2ECC71",
      Enasidenib        = "#9B59B6",
      Ivosidenib        = "#F39C12",
      Cytarabine        = "#1ABC9C",
      Idarubicin        = "#E67E22"
    )

    p <- plot_ly()
    for (d in unique(combined$drug)) {
      sub <- combined[combined$drug == d, ]
      col <- drug_colors[[d]] %||% "#7F8C8D"
      lty <- if (grepl("Fasted", d)) "dash" else "solid"
      p <- add_trace(p,
        data = sub, x = ~time, y = ~conc, type = "scatter", mode = "lines",
        name = d,
        line = list(color = col, width = 2, dash = lty)
      )
    }
    p %>%
      layout(
        xaxis = list(title = "Time (hours)"),
        yaxis = list(title = "Plasma Concentration (ng/mL)"),
        legend = list(orientation = "h", y = -0.2),
        hovermode = "x unified"
      )
  })

  output$pk_param_table <- renderTable({
    req(length(input$pk_drugs) > 0)
    rows <- lapply(input$pk_drugs, function(d) {
      p <- pk_params_db[[d]]
      cmax <- if (input$ven_food_effect && d == "Venetoclax") p$Cmax_fed else p$Cmax_fasted
      auc  <- if (input$ven_food_effect && d == "Venetoclax") p$AUC_fed  else p$AUC_fasted
      if (d == "Venetoclax" && input$ven_cyp_inhibitor) {
        cmax <- cmax * p$cyp3a4_fold
        auc  <- auc  * p$cyp3a4_fold
      }
      data.frame(
        Drug   = d,
        Dose   = paste(p$dose, p$unit),
        Route  = p$route,
        Cmax   = format(round(cmax, 1), big.mark = ","),
        Tmax_h = p$Tmax,
        AUC    = format(round(auc, 0), big.mark = ","),
        t_half_h = p$t_half,
        stringsAsFactors = FALSE
      )
    })
    do.call(rbind, rows) |> setNames(c("Drug", "Dose", "Route", "Cmax (ng/mL)", "Tmax (h)", "AUC (ng·h/mL)", "t½ (h)"))
  }, striped = TRUE, hover = TRUE, bordered = TRUE)

  output$pk_notes_ui <- renderUI({
    notes <- lapply(input$pk_drugs, function(d) {
      p <- pk_params_db[[d]]
      div(
        strong(d, ": "),
        span(p$note),
        br()
      )
    })
    div(class = "small", notes)
  })

  # -----------------------------------------------------------
  # TAB 3: Leukemia Dynamics (PD)
  # -----------------------------------------------------------

  output$pd_blast_plot <- renderPlotly({
    df <- aml_sim()
    cr_day <- df$cr_day[1]
    p <- plot_ly(df, x = ~day, y = ~bm_blast, type = "scatter", mode = "lines",
      line = list(color = "#E74C3C", width = 2.5), name = "BM Blasts")
    p <- add_trace(p, x = ~day, y = ~pb_blast, type = "scatter", mode = "lines",
      line = list(color = "#E67E22", width = 2, dash = "dash"), name = "PB Blasts")
    p <- add_lines(p, x = c(0, max(df$day)), y = c(5, 5),
      line = list(color = "green", dash = "dot", width = 1.5),
      name = "CR threshold (5%)", showlegend = TRUE)
    if (!is.na(cr_day)) {
      p <- add_segments(p, x = cr_day, xend = cr_day, y = 0, yend = 80,
        line = list(color = "blue", dash = "dot"), name = paste0("CR day ", cr_day))
    }
    p %>% layout(
      xaxis = list(title = "Day"),
      yaxis = list(title = "Blast % / count", range = c(0, 100)),
      legend = list(orientation = "h", y = -0.2)
    )
  })

  output$pd_lsc_mrd_plot <- renderPlotly({
    df <- aml_sim()
    p <- plot_ly(df, x = ~day, y = ~mrd, type = "scatter", mode = "lines",
      line = list(color = "#9B59B6", width = 2.5), name = "MRD (log10 copies)")
    p <- add_trace(p, x = ~day, y = ~log10(pmax(lsc, 1e-4)), type = "scatter", mode = "lines",
      line = list(color = "#2ECC71", width = 2, dash = "dash"), name = "LSC burden (log10)")
    p <- add_lines(p, x = c(0, max(df$day)), y = c(0, 0),
      line = list(color = "gray", dash = "dot"), name = "MRD negative (<1 copy)")
    p %>% layout(
      xaxis = list(title = "Day"),
      yaxis = list(title = "Log10 units"),
      legend = list(orientation = "h", y = -0.2)
    )
  })

  output$pd_cbc_plot <- renderPlotly({
    df <- aml_sim()
    p <- plot_ly()
    p <- add_trace(p, data = df, x = ~day, y = ~wbc,
      type = "scatter", mode = "lines", name = "WBC (10³/µL)",
      line = list(color = "#3498DB", width = 2))
    p <- add_trace(p, data = df, x = ~day, y = ~anc,
      type = "scatter", mode = "lines", name = "ANC (10³/µL)",
      line = list(color = "#E74C3C", width = 2))
    p <- add_trace(p, data = df, x = ~day, y = ~plt / 30,
      type = "scatter", mode = "lines", name = "PLT/30 (10⁹/L ÷30)",
      line = list(color = "#F39C12", width = 2, dash = "dash"))
    p <- add_trace(p, data = df, x = ~day, y = ~hgb,
      type = "scatter", mode = "lines", name = "Hgb (g/dL)",
      line = list(color = "#2ECC71", width = 2))
    # ANC < 0.5 warning line
    p <- add_lines(p, x = c(0, max(df$day)), y = c(0.5, 0.5),
      line = list(color = "red", dash = "dot", width = 1),
      name = "ANC 0.5 (fever threshold)")
    p %>% layout(
      xaxis = list(title = "Day"),
      yaxis = list(title = "CBC value (see legend for units)"),
      legend = list(orientation = "h", y = -0.25),
      hovermode = "x unified"
    )
  })

  output$pd_cr_info_ui <- renderUI({
    df <- aml_sim()
    cr_day <- df$cr_day[1]
    relapse_day <- df$relapse_day[1]

    cr_info <- if (!is.na(cr_day)) {
      div(class = "alert alert-success p-2",
        icon("check-circle"),
        strong(paste0(" CR achieved: Day ", cr_day, " (Cycle ", ceiling(cr_day / 28), ")")),
        br(),
        sprintf("BM blasts at CR: %.1f%%", df$bm_blast[df$day == min(df$day[df$day >= cr_day])][1])
      )
    } else {
      div(class = "alert alert-danger p-2", icon("times-circle"), " CR not achieved in simulation window.")
    }

    relapse_info <- if (!is.na(relapse_day)) {
      div(class = "alert alert-warning p-2",
        icon("exclamation-triangle"),
        strong(paste0(" Relapse detected: Day ", relapse_day))
      )
    } else {
      div(class = "alert alert-info p-2", icon("info-circle"), " No relapse detected in window.")
    }

    div(cr_info, relapse_info)
  })

  # -----------------------------------------------------------
  # TAB 4: Clinical Endpoints
  # -----------------------------------------------------------

  ep_os_data <- reactive({
    bind_rows(lapply(input$ep_regimens, simulate_os,
      n_months = input$ep_months, risk_group = input$ep_risk))
  })

  output$ep_os_plot <- renderPlotly({
    df <- ep_os_data()
    palette <- c("7plus3" = "#E74C3C", "VEN_AZA" = "#3498DB",
                 "Gilteritinib" = "#2ECC71", "Enasidenib" = "#9B59B6",
                 "ATRA_ATO" = "#F39C12")
    p <- plot_ly()
    for (r in unique(df$regimen)) {
      sub <- df[df$regimen == r, ]
      p <- add_trace(p, data = sub, x = ~month, y = ~os,
        type = "scatter", mode = "lines", name = r,
        line = list(color = palette[[r]] %||% "#7F8C8D", width = 2.5))
    }
    p %>% layout(
      xaxis = list(title = "Months"),
      yaxis = list(title = "Overall Survival (%)", range = c(0, 105)),
      legend = list(orientation = "h", y = -0.2)
    )
  })

  output$ep_efs_plot <- renderPlotly({
    df <- ep_os_data()
    palette <- c("7plus3" = "#E74C3C", "VEN_AZA" = "#3498DB",
                 "Gilteritinib" = "#2ECC71", "Enasidenib" = "#9B59B6",
                 "ATRA_ATO" = "#F39C12")
    p <- plot_ly()
    for (r in unique(df$regimen)) {
      sub <- df[df$regimen == r, ]
      p <- add_trace(p, data = sub, x = ~month, y = ~efs,
        type = "scatter", mode = "lines", name = r,
        line = list(color = palette[[r]] %||% "#7F8C8D", width = 2.5))
    }
    p %>% layout(
      xaxis = list(title = "Months"),
      yaxis = list(title = "Event-Free Survival (%)", range = c(0, 105)),
      legend = list(orientation = "h", y = -0.2)
    )
  })

  output$ep_cuminc_plot <- renderPlotly({
    df <- ep_os_data()
    reg <- unique(df$regimen)[1]
    sub <- df[df$regimen == reg, ]
    p <- plot_ly(sub, x = ~month)
    p <- add_trace(p, y = ~relapse, type = "scatter", mode = "lines",
      name = "Cumulative Relapse", line = list(color = "#E74C3C", width = 2))
    p <- add_trace(p, y = ~nrm, type = "scatter", mode = "lines",
      name = "Non-Relapse Mortality", line = list(color = "#7F8C8D", width = 2))
    p %>% layout(
      xaxis = list(title = "Months"),
      yaxis = list(title = "Cumulative Incidence (%)"),
      title = list(text = paste("First selected regimen:", reg)),
      legend = list(orientation = "h", y = -0.2)
    )
  })

  output$ep_summary_table <- renderTable({
    data.frame(
      Endpoint        = c("CR/CRi Rate", "Median OS (mos)", "1-yr OS (%)", "2-yr OS (%)",
                          "MRD-neg Rate (%)", "Grade 3-4 AE (%)"),
      "Best Supportive Care" = c("N/A", "2.7", "11", "5", "N/A", "—"),
      "VEN + AZA"     = c("66%", "14.7", "58", "37", "33", "78"),
      "Gilteritinib"  = c("54%", "9.3", "37", "21", "26", "89"),
      "7+3 Induction" = c("70%", "18.0", "62", "40", "45", "96"),
      "ATRA + ATO"    = c("95%", "NR", "98", "95", "97", "45"),
      "Source"        = c("—", "VIALE-A 2020", "ADMIRAL 2019", "RATIFY 2017", "Multiple", "Various"),
      check.names = FALSE, stringsAsFactors = FALSE
    )
  }, striped = TRUE, hover = TRUE, bordered = TRUE, digits = 1)

  # -----------------------------------------------------------
  # TAB 5: Scenario Comparison
  # -----------------------------------------------------------

  regimen_data <- reactive({
    all_regimens <- c("7plus3", "VEN_AZA", "Gilteritinib", "Enasidenib", "LDAC_VEN", "FLAG_IDA", "ATRA_ATO")
    cycle_day    <- input$sc_cycle * 28

    params_base <- list(
      age       = input$age,
      weight    = input$weight,
      wbc0      = input$wbc_dx,
      bm_blast0 = input$bm_blast,
      pb_blast0 = input$pb_blast,
      mutations = input$mutations,
      cyto_risk = input$cyto_risk
    )

    rows <- lapply(all_regimens, function(reg) {
      params_base$regimen <- reg
      df <- simulate_aml(params_base, days = seq(0, cycle_day, by = 1))
      init_blast <- df$bm_blast[1]
      end_blast  <- df$bm_blast[nrow(df)]
      pct_change <- (end_blast - init_blast) / init_blast * 100

      data.frame(
        regimen    = reg,
        end_blast  = end_blast,
        pct_change = pct_change,
        stringsAsFactors = FALSE
      )
    })
    do.call(rbind, rows)
  })

  output$sc_main_ui <- renderUI({
    view <- input$sc_view
    if (view == "Bar chart (blast reduction)") {
      plotlyOutput("sc_bar_plot", height = "420px")
    } else if (view == "Radar chart (multi-domain)") {
      plotlyOutput("sc_radar_plot", height = "500px")
    } else {
      tableOutput("sc_table_view")
    }
  })

  output$sc_bar_plot <- renderPlotly({
    df <- regimen_data()
    labels <- c(
      "7plus3"       = "7+3 Induction",
      "VEN_AZA"      = "VEN + AZA",
      "Gilteritinib" = "Gilteritinib",
      "Enasidenib"   = "Enasidenib",
      "LDAC_VEN"     = "LDAC + VEN",
      "FLAG_IDA"     = "FLAG-IDA",
      "ATRA_ATO"     = "ATRA + ATO"
    )
    colors <- c("#E74C3C", "#3498DB", "#2ECC71", "#9B59B6", "#F39C12", "#1ABC9C", "#E67E22")
    df$label <- sapply(df$regimen, function(r) labels[[r]] %||% r)
    df$color <- colors[seq_len(nrow(df))]

    plot_ly(df, x = ~label, y = ~pct_change, type = "bar",
      marker = list(color = colors[seq_len(nrow(df))]),
      text = ~paste0(round(pct_change, 1), "%"),
      textposition = "outside"
    ) %>%
    layout(
      xaxis = list(title = "Regimen"),
      yaxis = list(title = paste0("BM Blast % change at cycle ", input$sc_cycle)),
      bargap = 0.3
    )
  })

  output$sc_radar_plot <- renderPlotly({
    # Multi-domain radar: efficacy, tolerability, MRD neg rate, cost efficiency, ease
    regimen_labels <- c("7+3 Induction", "VEN + AZA", "Gilteritinib",
                        "Enasidenib", "LDAC + VEN", "FLAG-IDA", "ATRA + ATO")
    efficacy     <- c(85, 76, 64, 52, 60, 78, 95)
    tolerability <- c(30, 60, 65, 70, 75, 25, 85)
    mrd_neg      <- c(50, 38, 30, 25, 20, 45, 95)
    cost_eff     <- c(85, 40, 45, 42, 55, 60, 80)
    ease_admin   <- c(50, 65, 75, 78, 70, 30, 65)

    theta_vals <- c("Efficacy", "Tolerability", "MRD Negativity", "Cost Efficiency", "Ease of Admin", "Efficacy")
    colors_list <- c("#E74C3C", "#3498DB", "#2ECC71", "#9B59B6", "#F39C12", "#1ABC9C", "#E67E22")

    p <- plot_ly(type = "scatterpolar", fill = "toself")
    for (i in seq_along(regimen_labels)) {
      vals <- c(efficacy[i], tolerability[i], mrd_neg[i], cost_eff[i], ease_admin[i], efficacy[i])
      p <- add_trace(p,
        r     = vals,
        theta = theta_vals,
        name  = regimen_labels[i],
        line  = list(color = colors_list[i])
      )
    }
    p %>% layout(
      polar  = list(radialaxis = list(visible = TRUE, range = c(0, 100))),
      legend = list(orientation = "h", y = -0.15)
    )
  })

  output$sc_table_view <- renderTable({
    data.frame(
      Regimen        = c("7+3 Induction", "VEN + AZA (VIALE-A)", "Gilteritinib (ADMIRAL)",
                         "Enasidenib (AG221)", "LDAC + VEN (VIALE-C)", "FLAG-IDA", "ATRA + ATO"),
      "CR/CRi Rate"  = c("70%", "66%", "54%", "40%", "48%", "60%", "95%"),
      "Median OS"    = c("18 mos", "14.7 mos", "9.3 mos", "8.8 mos", "10.1 mos", "11 mos", "NR"),
      "1-yr OS"      = c("62%", "58%", "37%", "30%", "41%", "42%", "98%"),
      "Gr3-4 AE"     = c("96%", "78%", "89%", "82%", "75%", "98%", "45%"),
      "MRD Neg Rate" = c("45%", "33%", "26%", "18%", "15%", "42%", "97%"),
      "Key Trial"    = c("Standard", "NCT02993523", "NCT02421939", "NCT01915498",
                         "NCT03069352", "Retrospective", "APL0406"),
      check.names = FALSE, stringsAsFactors = FALSE
    )
  }, striped = TRUE, hover = TRUE, bordered = TRUE)

  output$sc_trial_table <- renderTable({
    data.frame(
      Trial        = c("VIALE-A", "ADMIRAL", "QuANTUM-R", "RATIFY", "VIALE-C", "APL0406", "AML-001"),
      Drug         = c("VEN+AZA", "Gilteritinib", "Quizartinib", "Midostaurin+7+3", "VEN+LDAC", "ATRA+ATO", "Azacitidine"),
      Population   = c("Unfit newly dx", "FLT3+ R/R", "FLT3+ R/R", "FLT3+ newly dx", "Unfit newly dx", "Newly dx APL", "Elderly unfit"),
      N            = c("431", "371", "367", "717", "211", "162", "488"),
      "CR/CRi"     = c("66%", "34%", "48%", "59%", "27%", "95%", "28%"),
      "Median OS (mos)" = c("14.7", "9.3", "6.2", "74.7", "7.2", "NR", "10.4"),
      PMID         = c("32515654", "31665578", "31068292", "28591536", "32271986", "23841725", "19295674"),
      check.names = FALSE, stringsAsFactors = FALSE
    )
  }, striped = TRUE, hover = TRUE, bordered = TRUE)

  # -----------------------------------------------------------
  # TAB 6: Biomarker Dashboard
  # -----------------------------------------------------------

  bm_current_day <- reactive({
    input$run_bm
    isolate(input$bm_days)
  })

  output$bm_mrd_status_box <- renderUI({
    df <- bm_sim()
    curr_day <- min(bm_current_day(), max(df$day))
    mrd_val  <- df$mrd[df$day == curr_day][1]
    if (is.na(mrd_val)) mrd_val <- df$mrd[nrow(df)]
    status <- if (mrd_val < 0) "MRD Negative" else "MRD Positive"
    color  <- if (mrd_val < 0) "#27AE60" else "#E74C3C"
    div(
      h4(status, style = paste0("color:", color, "; font-weight:bold")),
      p(sprintf("MRD: %.2f log10", mrd_val), style = "font-size:0.85em; color:#7F8C8D")
    )
  })

  output$bm_days_to_mrd_box <- renderUI({
    df <- bm_sim()
    neg_days <- df$day[df$mrd < 0]
    days_val <- if (length(neg_days) > 0) min(neg_days) else NA
    div(
      h4(
        if (is.na(days_val)) "Not reached" else paste(days_val, "days"),
        style = "color:#2980B9; font-weight:bold"
      ),
      p(if (is.na(days_val)) "MRD positivity persists" else paste0("(~Cycle ", ceiling(days_val/28), ")"),
        style = "font-size:0.85em; color:#7F8C8D")
    )
  })

  output$bm_relapse_risk_box <- renderUI({
    df <- bm_sim()
    # Approximate 12-month relapse risk from LSC burden
    lsc_at_12m <- df$lsc[df$day == min(360, max(df$day))][1]
    lsc_init   <- df$lsc[1]
    relapse_risk <- pmin(pmax((lsc_at_12m / lsc_init) * 60, 5), 95)
    color <- if (relapse_risk < 20) "#27AE60" else if (relapse_risk < 50) "#F39C12" else "#E74C3C"
    div(
      h4(sprintf("%.0f%%", relapse_risk), style = paste0("color:", color, "; font-weight:bold")),
      p("Estimated from LSC kinetics", style = "font-size:0.85em; color:#7F8C8D")
    )
  })

  output$bm_next_assessment_box <- renderUI({
    curr_day <- bm_current_day()
    next_day <- curr_day + 28
    div(
      h4(paste("Day", next_day), style = "color:#9B59B6; font-weight:bold"),
      p(paste0("(+", 28, " days / BM biopsy + MRD)"), style = "font-size:0.85em; color:#7F8C8D")
    )
  })

  output$bm_mrd_plot <- renderPlotly({
    df <- bm_sim()
    curr_day <- bm_current_day()

    p <- plot_ly(df, x = ~day, y = ~mrd, type = "scatter", mode = "lines",
      line = list(color = "#9B59B6", width = 2.5), name = "MRD (log10)")
    # NPM1-MRD (slightly noisier)
    if ("NPM1-MRD (PCR)" %in% input$bm_show) {
      npm1_mrd <- df$mrd + rnorm(nrow(df), 0, 0.15)
      p <- add_trace(p, x = df$day, y = npm1_mrd, type = "scatter", mode = "lines",
        line = list(color = "#E74C3C", width = 2, dash = "dot"), name = "NPM1-MRD")
    }
    p <- add_lines(p, x = c(0, max(df$day)), y = c(0, 0),
      line = list(color = "green", dash = "dot", width = 1.5), name = "MRD neg threshold")
    p <- add_lines(p, x = c(curr_day, curr_day), y = c(-2, max(df$mrd, na.rm = TRUE) + 0.5),
      line = list(color = "orange", dash = "dot", width = 2), name = "Current day")
    p %>% layout(
      xaxis = list(title = "Day"),
      yaxis = list(title = "MRD (log10 copies)"),
      legend = list(orientation = "h", y = -0.2)
    )
  })

  output$bm_flt3_plot <- renderPlotly({
    if (!"FLT3-ITD VAF" %in% input$bm_show) {
      return(plotly_empty() %>% layout(title = "FLT3-ITD VAF (select in sidebar)"))
    }
    df <- bm_sim()

    # FLT3 VAF tracks blast burden with some noise
    has_flt3 <- "FLT3-ITD" %in% (input$mutations %||% character(0))
    if (!has_flt3) {
      return(plotly_empty() %>% layout(title = "FLT3-ITD not selected in patient mutations"))
    }

    vaf_init <- 45  # baseline ~45% VAF
    flt3_vaf <- vaf_init * (df$bm_blast / df$bm_blast[1]) + rnorm(nrow(df), 0, 1.5)
    flt3_vaf <- pmax(pmin(flt3_vaf, 100), 0)

    plot_ly(x = df$day, y = flt3_vaf, type = "scatter", mode = "lines",
      line = list(color = "#F39C12", width = 2.5), name = "FLT3-ITD VAF (%)") %>%
    add_lines(x = c(0, max(df$day)), y = c(5, 5),
      line = list(color = "green", dash = "dot", width = 1.5), name = "VAF 5% threshold") %>%
    layout(
      xaxis = list(title = "Day"),
      yaxis = list(title = "FLT3-ITD VAF (%)", range = c(0, 55)),
      legend = list(orientation = "h", y = -0.2)
    )
  })

  output$bm_cbc_plot <- renderPlotly({
    if (!"WBC/ANC/PLT/Hgb" %in% input$bm_show) {
      return(plotly_empty() %>% layout(title = "CBC (select in sidebar)"))
    }
    df <- bm_sim()
    p <- plot_ly()
    # Normal range bands
    p <- add_ribbons(p, x = df$day, ymin = rep(4, nrow(df)), ymax = rep(11, nrow(df)),
      fillcolor = "rgba(52,152,219,0.08)", line = list(width = 0), name = "WBC normal range", showlegend = TRUE)
    p <- add_trace(p, data = df, x = ~day, y = ~wbc, type = "scatter", mode = "lines",
      name = "WBC (10³/µL)", line = list(color = "#3498DB", width = 2))
    p <- add_trace(p, data = df, x = ~day, y = ~anc, type = "scatter", mode = "lines",
      name = "ANC (10³/µL)", line = list(color = "#E74C3C", width = 2))
    p <- add_trace(p, data = df, x = ~day, y = ~hgb, type = "scatter", mode = "lines",
      name = "Hgb (g/dL)", line = list(color = "#2ECC71", width = 2))
    p <- add_trace(p, data = df, x = ~day, y = ~plt / 20, type = "scatter", mode = "lines",
      name = "PLT/20 (×10⁹/L)", line = list(color = "#F39C12", width = 2, dash = "dash"))
    p %>% layout(
      xaxis = list(title = "Day"),
      yaxis = list(title = "Value (mixed units — see legend)"),
      legend = list(orientation = "h", y = -0.3)
    )
  })

  output$bm_inflam_plot <- renderPlotly({
    df <- bm_sim()
    days <- df$day

    # CRP mirrors blast burden (inflammation proxy)
    crp_val <- 2 + 15 * (df$bm_blast / 100) + rnorm(length(days), 0, 1)
    crp_val <- pmax(crp_val, 0.5)

    # Ferritin (elevated in active disease, high after multiple transfusions)
    ferritin_val <- 300 + 800 * (df$bm_blast / 100) + cumsum(rnorm(length(days), 0, 5))
    ferritin_val <- pmax(ferritin_val, 50)

    # BM Cellularity (inverse of blast kill)
    bm_cellularity <- 80 * (0.3 + 0.7 * (1 - df$bm_blast / 100)) + rnorm(length(days), 0, 3)
    bm_cellularity <- pmax(pmin(bm_cellularity, 100), 5)

    p <- plot_ly()
    if ("CRP / Ferritin" %in% input$bm_show) {
      p <- add_trace(p, x = days, y = crp_val, type = "scatter", mode = "lines",
        name = "CRP (mg/L)", line = list(color = "#E74C3C", width = 2))
      p <- add_trace(p, x = days, y = ferritin_val / 100, type = "scatter", mode = "lines",
        name = "Ferritin/100 (µg/L)", line = list(color = "#E67E22", width = 2, dash = "dash"))
    }
    if ("BM Cellularity" %in% input$bm_show) {
      p <- add_trace(p, x = days, y = bm_cellularity, type = "scatter", mode = "lines",
        name = "BM Cellularity (%)", line = list(color = "#2ECC71", width = 2))
    }
    p %>% layout(
      xaxis = list(title = "Day"),
      yaxis = list(title = "Value (mixed units)"),
      legend = list(orientation = "h", y = -0.25)
    )
  })

}

# ============================================================
# Launch Application
# ============================================================

shinyApp(ui = ui, server = server)
