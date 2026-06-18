# ============================================================
# BPH QSP Model: Interactive Disease Simulation Shiny App
# Benign Prostatic Hyperplasia (BPH) - Quantitative Systems
# Pharmacology (QSP) Model
#
# Model covers:
#   - Drug PK for tamsulosin, finasteride, dutasteride, tadalafil
#   - Hormonal biomarkers (DHT, testosterone, PSA)
#   - Urodynamic endpoints (IPSS, Qmax, PVR, prostate volume)
#   - Treatment scenario comparisons (MTOPS, CombAT, NEPTUNE)
#   - Disease progression and AUR/TURP risk
#
# References:
#   - MTOPS trial (McConnell et al., N Engl J Med 2003)
#   - CombAT trial (Roehrborn et al., Eur Urol 2010)
#   - NEPTUNE trial (Oelke et al., Eur Urol 2012)
#   - REDUCE trial (Andriole et al., N Engl J Med 2010)
# ============================================================

library(shiny)
library(ggplot2)
library(dplyr)
library(plotly)
library(deSolve)
library(DT)

# ============================================================
# PK SIMULATION FUNCTIONS
# ============================================================

#' Simulate drug PK using a one-compartment model with first-order absorption
#' @param drug character: one of "tamsulosin", "finasteride", "dutasteride", "tadalafil"
#' @param time numeric vector: time points in hours
#' @param dose numeric: dose in mg (default uses approved dose)
#' @return data.frame with columns: time, concentration (ng/mL)
simulate_pk <- function(drug, time, dose = NULL) {
  # PK parameter database
  pk_params <- list(
    tamsulosin = list(
      dose_mg  = 0.4,    # mg
      ka       = 0.693,  # /h (Tmax ~6h)
      CL       = 1.8,    # L/h
      V1       = 15,     # L
      F        = 0.90,   # bioavailability
      t_half   = 10,     # h
      prot_bind = 94,    # %
      units    = "ng/mL",
      mw       = 408.5   # g/mol (for unit conversion if needed)
    ),
    finasteride = list(
      dose_mg  = 5.0,
      ka       = 0.50,
      CL       = 5.0,
      V1       = 76,
      F        = 0.63,
      t_half   = 6,
      prot_bind = 90,
      units    = "ng/mL",
      mw       = 372.5
    ),
    dutasteride = list(
      dose_mg  = 0.5,
      ka       = 0.08,
      CL       = 0.15,
      V1       = 5.0,
      F        = 0.60,
      t_half   = 840,   # ~5 weeks
      prot_bind = 99,
      units    = "ng/mL",
      mw       = 528.5
    ),
    tadalafil = list(
      dose_mg  = 5.0,
      ka       = 0.35,
      CL       = 2.5,
      V1       = 25,
      F        = 0.80,
      t_half   = 17,
      prot_bind = 94,
      units    = "ng/mL",
      mw       = 389.4
    )
  )

  p <- pk_params[[drug]]
  if (is.null(p)) stop("Unknown drug: ", drug)
  if (!is.null(dose)) p$dose_mg <- dose

  # Effective dose after bioavailability (convert mg to ug for ng/mL at L scale)
  dose_ug <- p$dose_mg * 1000 * p$F  # ug

  # ke = elimination rate constant
  ke <- log(2) / p$t_half  # /h

  # One-compartment model: C(t) = (F*D*ka) / (V*(ka-ke)) * (exp(-ke*t) - exp(-ka*t))
  # This gives concentration in ug/L = ng/mL
  conc <- (dose_ug * p$ka) / (p$V1 * (p$ka - ke)) * (exp(-ke * time) - exp(-ka * time))
  conc[conc < 0] <- 0

  data.frame(
    time          = time,
    concentration = conc,
    drug          = drug
  )
}

#' Simulate steady-state PK (multiple dosing at tau interval)
#' @param drug character
#' @param tau dosing interval in hours (default 24 = QD)
#' @param n_doses number of doses to simulate (default enough to reach SS)
#' @param time_points time points within one dosing interval (0 to tau)
simulate_pk_ss <- function(drug, tau = 24, time_points = seq(0, 24, by = 0.5)) {
  pk_params <- list(
    tamsulosin  = list(dose_mg = 0.4,  ka = 0.693, CL = 1.8,  V1 = 15,  F = 0.90, t_half = 10,  prot_bind = 94),
    finasteride = list(dose_mg = 5.0,  ka = 0.50,  CL = 5.0,  V1 = 76,  F = 0.63, t_half = 6,   prot_bind = 90),
    dutasteride = list(dose_mg = 0.5,  ka = 0.08,  CL = 0.15, V1 = 5.0, F = 0.60, t_half = 840, prot_bind = 99),
    tadalafil   = list(dose_mg = 5.0,  ka = 0.35,  CL = 2.5,  V1 = 25,  F = 0.80, t_half = 17,  prot_bind = 94)
  )

  p  <- pk_params[[drug]]
  ke <- log(2) / p$t_half
  ka <- p$ka
  D  <- p$dose_mg * 1000 * p$F  # ug

  # Superposition for steady state: sum over past doses (enough to converge)
  n_prev <- ceiling(5 * p$t_half / tau) + 5  # enough prior doses for SS
  conc <- numeric(length(time_points))

  for (i in 0:n_prev) {
    t_shifted <- time_points + i * tau
    c_i <- (D * ka) / (p$V1 * (ka - ke)) * (exp(-ke * t_shifted) - exp(-ka * t_shifted))
    c_i[c_i < 0] <- 0
    conc <- conc + c_i
  }

  data.frame(
    time          = time_points,
    concentration = conc,
    drug          = drug
  )
}

#' Simulate 1-year daily Cmax profile (approach to steady state)
simulate_pk_yearly <- function(drug) {
  pk_params <- list(
    tamsulosin  = list(dose_mg = 0.4,  ka = 0.693, CL = 1.8,  V1 = 15,  F = 0.90, t_half = 10,  prot_bind = 94),
    finasteride = list(dose_mg = 5.0,  ka = 0.50,  CL = 5.0,  V1 = 76,  F = 0.63, t_half = 6,   prot_bind = 90),
    dutasteride = list(dose_mg = 0.5,  ka = 0.08,  CL = 0.15, V1 = 5.0, F = 0.60, t_half = 840, prot_bind = 99),
    tadalafil   = list(dose_mg = 5.0,  ka = 0.35,  CL = 2.5,  V1 = 25,  F = 0.80, t_half = 17,  prot_bind = 94)
  )

  p   <- pk_params[[drug]]
  ke  <- log(2) / p$t_half
  ka  <- p$ka
  D   <- p$dose_mg * 1000 * p$F
  tau <- 24  # QD

  days      <- 1:365
  cmax_vals <- numeric(365)

  # Fine time grid within each day
  t_fine <- seq(0, tau, by = 0.5)
  t_peak <- log(ka / ke) / (ka - ke)
  if (t_peak < 0 || t_peak > tau) t_peak <- tau / 2

  for (day in days) {
    # Concentration at peak time within that dosing interval
    conc_day <- numeric(length(t_fine))
    for (j in 0:(day - 1)) {
      t_sh <- t_fine + j * tau
      c_j  <- (D * ka) / (p$V1 * (ka - ke)) * (exp(-ke * t_sh) - exp(-ka * t_sh))
      c_j[c_j < 0] <- 0
      conc_day <- conc_day + c_j
    }
    cmax_vals[day] <- max(conc_day)
  }

  data.frame(
    day           = days,
    cmax          = cmax_vals,
    drug          = drug
  )
}

# ============================================================
# PD SIMULATION FUNCTIONS
# ============================================================

#' Simulate PD endpoints (IPSS, Qmax, PVR, Prostate Volume) over time
#' Uses closed-form approximations derived from clinical trial data
#' @param treatment character: one of "watchful_waiting", "tamsulosin", "finasteride",
#'                  "dutasteride", "combination", "tadalafil"
#' @param baseline_pv  numeric: baseline prostate volume (cc)
#' @param baseline_ipss numeric: baseline IPSS score
#' @param baseline_qmax numeric: baseline Qmax (mL/s)
#' @param time_weeks numeric vector: weeks from baseline
#' @return data.frame with time_weeks, ipss, qmax, pvr, pv columns
simulate_pd <- function(treatment, baseline_pv, baseline_ipss, baseline_qmax,
                        time_weeks) {

  # --- Sigmoid Emax helper (for smooth approach to plateau effect) ---
  # E(t) = Emax * t^hill / (ET50^hill + t^hill)
  sigmoid_emax <- function(t, emax, et50, hill = 1.5) {
    emax * t^hill / (et50^hill + t^hill)
  }

  # --- Baseline PVR estimate from Qmax (empirical) ---
  baseline_pvr <- pmax(10, 350 - 20 * baseline_qmax)

  # --- Natural progression (watchful waiting) ---
  # IPSS worsens ~0.5 pts/yr; PV grows ~1.5 cc/yr; Qmax declines ~0.2 mL/s/yr
  nat_ipss_slope <- 0.5  / 52  # per week
  nat_pv_slope   <- 1.5  / 52  # cc per week
  nat_qmax_slope <- 0.2  / 52  # mL/s per week (decline)
  nat_pvr_slope  <- 3.0  / 52  # mL per week increase

  ipss <- numeric(length(time_weeks))
  qmax <- numeric(length(time_weeks))
  pvr  <- numeric(length(time_weeks))
  pv   <- numeric(length(time_weeks))

  for (i in seq_along(time_weeks)) {
    t <- time_weeks[i]
    yr <- t / 52

    # Natural progression component
    nat_ipss <- baseline_ipss + nat_ipss_slope * t
    nat_qmax <- baseline_qmax - nat_qmax_slope * t
    nat_pv   <- baseline_pv   + nat_pv_slope   * t
    nat_pvr  <- baseline_pvr  + nat_pvr_slope  * t

    if (treatment == "watchful_waiting") {
      ipss[i] <- nat_ipss
      qmax[i] <- max(1, nat_qmax)
      pv[i]   <- nat_pv
      pvr[i]  <- nat_pvr

    } else if (treatment == "tamsulosin") {
      # Tamsulosin: fast functional relief, max ~-4.5 IPSS at 4-6 wks; no volume change
      # Qmax improves ~1.8 mL/s; PVR decreases
      delta_ipss <- -sigmoid_emax(t, emax = 4.5, et50 = 4, hill = 2.0)
      delta_qmax <- +sigmoid_emax(t, emax = 1.8, et50 = 4, hill = 2.0)
      delta_pvr  <- -sigmoid_emax(t, emax = 30,  et50 = 4, hill = 2.0)
      ipss[i] <- nat_ipss + delta_ipss
      qmax[i] <- max(1, nat_qmax + delta_qmax)
      pv[i]   <- nat_pv  # no volume change
      pvr[i]  <- max(0, nat_pvr + delta_pvr)

    } else if (treatment == "finasteride") {
      # Finasteride: slower onset (3-6 months for DHT suppression); volume reduction
      # ΔIPSS = -3 at 1yr, -3.5 at 2yr; ΔPV = -20% at 1yr, -27% at 2yr (PLESS)
      delta_ipss <- -sigmoid_emax(t, emax = 3.5, et50 = 26, hill = 1.2)
      delta_qmax <- +sigmoid_emax(t, emax = 1.5, et50 = 26, hill = 1.2)
      pv_frac    <-  sigmoid_emax(t, emax = 0.27, et50 = 26, hill = 1.2)
      delta_pvr  <- -sigmoid_emax(t, emax = 25, et50 = 26, hill = 1.2)
      ipss[i] <- nat_ipss + delta_ipss
      qmax[i] <- max(1, nat_qmax + delta_qmax)
      pv[i]   <- nat_pv * (1 - pv_frac)
      pvr[i]  <- max(0, nat_pvr + delta_pvr)

    } else if (treatment == "dutasteride") {
      # Dutasteride: dual 5-ARI, deeper DHT suppression (94%)
      # ΔIPSS = -4.5 at 1yr, -5.5 at 2yr; ΔPV = -23% at 1yr, -28% at 2yr (CombAT)
      delta_ipss <- -sigmoid_emax(t, emax = 5.5, et50 = 24, hill = 1.2)
      delta_qmax <- +sigmoid_emax(t, emax = 2.0, et50 = 24, hill = 1.2)
      pv_frac    <-  sigmoid_emax(t, emax = 0.28, et50 = 24, hill = 1.2)
      delta_pvr  <- -sigmoid_emax(t, emax = 30, et50 = 24, hill = 1.2)
      ipss[i] <- nat_ipss + delta_ipss
      qmax[i] <- max(1, nat_qmax + delta_qmax)
      pv[i]   <- nat_pv * (1 - pv_frac)
      pvr[i]  <- max(0, nat_pvr + delta_pvr)

    } else if (treatment == "combination") {
      # Combination (Dut + Tams): CombAT data
      # ΔIPSS = -6 at 1yr, -7 at 2yr; ΔPV = -28% at 2yr; fast onset due to tamsulosin
      delta_ipss_fast <- -sigmoid_emax(t, emax = 3.5, et50 = 4,  hill = 2.0)  # tams component
      delta_ipss_slow <- -sigmoid_emax(t, emax = 3.5, et50 = 26, hill = 1.2)  # dut component
      delta_ipss <- delta_ipss_fast + delta_ipss_slow * (1 - abs(delta_ipss_fast) / 7)
      delta_ipss <- pmax(delta_ipss, -7)

      delta_qmax_fast <- +sigmoid_emax(t, emax = 2.0, et50 = 4,  hill = 2.0)
      delta_qmax_slow <- +sigmoid_emax(t, emax = 1.5, et50 = 26, hill = 1.2)
      delta_qmax <- delta_qmax_fast + delta_qmax_slow * 0.5

      pv_frac  <- sigmoid_emax(t, emax = 0.28, et50 = 24, hill = 1.2)
      delta_pvr <- -sigmoid_emax(t, emax = 45, et50 = 10, hill = 1.5)

      ipss[i] <- nat_ipss + delta_ipss
      qmax[i] <- max(1, nat_qmax + delta_qmax)
      pv[i]   <- nat_pv * (1 - pv_frac)
      pvr[i]  <- max(0, nat_pvr + delta_pvr)

    } else if (treatment == "tadalafil") {
      # Tadalafil 5mg QD: NEPTUNE trial ΔIPSS = -5.6 at 12wks; Qmax +2.4 mL/s
      # No significant effect on prostate volume
      delta_ipss <- -sigmoid_emax(t, emax = 5.6, et50 = 6, hill = 1.8)
      delta_qmax <- +sigmoid_emax(t, emax = 2.4, et50 = 6, hill = 1.8)
      delta_pvr  <- -sigmoid_emax(t, emax = 20,  et50 = 6, hill = 1.8)
      ipss[i] <- nat_ipss + delta_ipss
      qmax[i] <- max(1, nat_qmax + delta_qmax)
      pv[i]   <- nat_pv  # no volume effect
      pvr[i]  <- max(0, nat_pvr + delta_pvr)
    }
  }

  # Clamp IPSS to valid range
  ipss <- pmax(0, pmin(35, ipss))

  data.frame(
    time_weeks = time_weeks,
    ipss       = ipss,
    qmax       = qmax,
    pvr        = pvr,
    pv         = pv,
    treatment  = treatment
  )
}

# ============================================================
# BIOMARKER SIMULATION (DHT, Testosterone, PSA)
# ============================================================

#' Simulate hormonal biomarkers over time
#' @param treatment character
#' @param baseline_psa numeric (ng/mL)
#' @param time_weeks numeric vector
#' @return data.frame
simulate_biomarkers <- function(treatment, baseline_psa, time_weeks) {
  # Normal baseline values
  dht_plasma_baseline     <- 1.8   # nmol/L (mid-normal range)
  dht_prostate_baseline   <- 25.0  # nmol/g (intraprostatic, ~10x plasma)
  testosterone_baseline   <- 15.0  # nmol/L

  sigmoid_emax <- function(t, emax, et50, hill = 1.5) {
    emax * t^hill / (et50^hill + t^hill)
  }

  dht_p   <- numeric(length(time_weeks))
  dht_pr  <- numeric(length(time_weeks))
  testo   <- numeric(length(time_weeks))
  psa     <- numeric(length(time_weeks))

  for (i in seq_along(time_weeks)) {
    t <- time_weeks[i]

    if (treatment %in% c("watchful_waiting", "tamsulosin", "tadalafil")) {
      # Minimal hormonal effect
      dht_p[i]  <- dht_plasma_baseline
      dht_pr[i] <- dht_prostate_baseline
      testo[i]  <- testosterone_baseline
      psa[i]    <- baseline_psa

    } else if (treatment == "finasteride") {
      # 5AR type-II inhibition: ~70% DHT reduction; testosterone rises ~15%
      dht_frac  <- sigmoid_emax(t, emax = 0.70, et50 = 12, hill = 1.3)
      psa_frac  <- sigmoid_emax(t, emax = 0.50, et50 = 16, hill = 1.3)
      dht_p[i]  <- dht_plasma_baseline   * (1 - dht_frac)
      dht_pr[i] <- dht_prostate_baseline * (1 - 0.80 * dht_frac)  # deeper in prostate
      testo[i]  <- testosterone_baseline * (1 + 0.15 * dht_frac)  # compensatory rise
      psa[i]    <- baseline_psa * (1 - psa_frac)

    } else if (treatment == "dutasteride") {
      # Dual 5AR inhibition: ~94% DHT reduction
      dht_frac  <- sigmoid_emax(t, emax = 0.94, et50 = 12, hill = 1.3)
      psa_frac  <- sigmoid_emax(t, emax = 0.53, et50 = 16, hill = 1.3)
      dht_p[i]  <- dht_plasma_baseline   * (1 - dht_frac)
      dht_pr[i] <- dht_prostate_baseline * (1 - 0.90 * dht_frac)
      testo[i]  <- testosterone_baseline * (1 + 0.20 * dht_frac)
      psa[i]    <- baseline_psa * (1 - psa_frac)

    } else if (treatment == "combination") {
      # Same dutasteride hormonal effect as mono
      dht_frac  <- sigmoid_emax(t, emax = 0.94, et50 = 12, hill = 1.3)
      psa_frac  <- sigmoid_emax(t, emax = 0.53, et50 = 16, hill = 1.3)
      dht_p[i]  <- dht_plasma_baseline   * (1 - dht_frac)
      dht_pr[i] <- dht_prostate_baseline * (1 - 0.90 * dht_frac)
      testo[i]  <- testosterone_baseline * (1 + 0.20 * dht_frac)
      psa[i]    <- baseline_psa * (1 - psa_frac)
    }
  }

  data.frame(
    time_weeks      = time_weeks,
    dht_plasma      = dht_p,
    dht_prostate    = dht_pr,
    testosterone    = testo,
    psa             = psa,
    treatment       = treatment
  )
}

# ============================================================
# RISK COMPUTATION
# ============================================================

#' Compute cumulative AUR and TURP risk over years
#' Based on MTOPS, CombAT, and REDUCE trial data
#' @param pv   numeric: prostate volume (cc)
#' @param psa  numeric: PSA (ng/mL)
#' @param age  numeric: patient age
#' @param ipss numeric: baseline IPSS
#' @param treatment character
#' @param years numeric: number of years for projection
#' @return data.frame with year, aur_risk, turp_risk
compute_risk <- function(pv, psa, age, ipss, treatment, years) {

  time_yr <- seq(0, years, by = 0.5)

  # Baseline annual AUR risk (Jacobsen et al., Olmsted County Study)
  # Modified by PV, PSA, age
  base_aur_annual  <- 0.015 * (1 + (pv > 40) * 0.8 + (psa > 1.5) * 0.6 + (age > 70) * 0.4)
  base_turp_annual <- 0.012 * (1 + (ipss > 20) * 0.7 + (pv > 40) * 0.5)

  # Risk reduction factors from clinical trials
  rr_aur <- switch(treatment,
    watchful_waiting = 1.00,
    tamsulosin       = 0.85,   # modest: MTOPS 4.4 vs 4.9% AUR
    finasteride      = 0.43,   # MTOPS 57% reduction
    dutasteride      = 0.32,   # REDUCE/CombAT 68% reduction
    combination      = 0.21,   # CombAT 79% reduction
    tadalafil        = 0.75,   # extrapolated; limited long-term AUR data
    1.00
  )

  rr_turp <- switch(treatment,
    watchful_waiting = 1.00,
    tamsulosin       = 0.87,
    finasteride      = 0.36,   # MTOPS 64% reduction
    dutasteride      = 0.29,   # CombAT 71% reduction
    combination      = 0.20,   # CombAT 80% reduction
    tadalafil        = 0.80,
    1.00
  )

  # Cumulative risk (1 - (1 - annual_rate)^t) with RR applied
  eff_aur_annual  <- base_aur_annual  * rr_aur
  eff_turp_annual <- base_turp_annual * rr_turp

  aur_cum  <- 1 - (1 - eff_aur_annual)^time_yr
  turp_cum <- 1 - (1 - eff_turp_annual)^time_yr

  data.frame(
    year      = time_yr,
    aur_risk  = aur_cum  * 100,   # percent
    turp_risk = turp_cum * 100
  )
}

# ============================================================
# HELPER: Patient risk category
# ============================================================
get_risk_category <- function(pv, psa, ipss) {
  criteria_met <- sum(c(pv > 40, psa > 1.5, ipss > 20))
  if (criteria_met >= 2) {
    if (pv > 40 && psa > 1.5 && ipss > 20) return("High")
    return("Moderate")
  }
  return("Low")
}

# ============================================================
# TREATMENT LABEL MAPPING
# ============================================================
treatment_labels <- c(
  watchful_waiting = "Watchful Waiting",
  tamsulosin       = "Tamsulosin 0.4mg QD",
  finasteride      = "Finasteride 5mg QD",
  dutasteride      = "Dutasteride 0.5mg QD",
  combination      = "Combination (Dut+Tams)",
  tadalafil        = "Tadalafil 5mg QD"
)

treatment_colors <- c(
  watchful_waiting = "#999999",
  tamsulosin       = "#E69F00",
  finasteride      = "#56B4E9",
  dutasteride      = "#009E73",
  combination      = "#CC79A7",
  tadalafil        = "#D55E00"
)

# ============================================================
# UI
# ============================================================

ui <- fluidPage(

  tags$head(
    tags$style(HTML("
      .risk-card { border-radius: 8px; padding: 14px 18px; margin-top: 12px; font-size: 15px; }
      .risk-low    { background: #d4edda; border-left: 5px solid #28a745; }
      .risk-mod    { background: #fff3cd; border-left: 5px solid #ffc107; }
      .risk-high   { background: #f8d7da; border-left: 5px solid #dc3545; }
      .value-box   { border-radius: 6px; padding: 10px 14px; margin: 6px 0; }
      .val-primary { background: #cce5ff; border-left: 4px solid #004085; }
      .val-success { background: #d4edda; border-left: 4px solid #155724; }
      .val-warning { background: #fff3cd; border-left: 4px solid #856404; }
      .val-danger  { background: #f8d7da; border-left: 4px solid #721c24; }
      h4.tab-title { color: #2c3e50; border-bottom: 2px solid #3498db; padding-bottom: 6px; }
    "))
  ),

  titlePanel(
    div(
      h2("BPH QSP Model: Interactive Disease Simulation"),
      p(style = "color:#666; font-size:14px;",
        "Benign Prostatic Hyperplasia — Quantitative Systems Pharmacology Dashboard")
    )
  ),

  tabsetPanel(
    id = "main_tabs",

    # --------------------------------------------------------
    # TAB 1: Patient Profile
    # --------------------------------------------------------
    tabPanel(
      "Patient Profile",
      fluidRow(
        column(4,
          wellPanel(
            h4("Patient Characteristics", class = "tab-title"),
            sliderInput("age",          "Age (years)",               min = 50, max = 80, value = 65, step = 1),
            sliderInput("baseline_pv",  "Baseline Prostate Volume (cc)", min = 20, max = 80, value = 40, step = 1),
            sliderInput("baseline_ipss","Baseline IPSS Score",       min = 8,  max = 30, value = 18, step = 1),
            sliderInput("baseline_qmax","Baseline Qmax (mL/s)",      min = 5,  max = 20, value = 10, step = 0.5),
            sliderInput("baseline_psa", "Baseline PSA (ng/mL)",      min = 0.5, max = 10, value = 2.5, step = 0.1)
          )
        ),
        column(4,
          wellPanel(
            h4("Clinical Context", class = "tab-title"),
            checkboxGroupInput("comorbidities", "Comorbidities:",
              choices  = c("Type 2 Diabetes (T2DM)" = "t2dm",
                           "Hypertension"            = "htn",
                           "Erectile Dysfunction (ED)" = "ed"),
              selected = character(0)
            ),
            hr(),
            selectInput("treatment_profile", "Primary Treatment:",
              choices = c(
                "Watchful Waiting"          = "watchful_waiting",
                "Tamsulosin 0.4mg QD"       = "tamsulosin",
                "Finasteride 5mg QD"        = "finasteride",
                "Dutasteride 0.5mg QD"      = "dutasteride",
                "Combination (Dut+Tams)"    = "combination",
                "Tadalafil 5mg QD"          = "tadalafil"
              ),
              selected = "tamsulosin"
            )
          )
        ),
        column(4,
          h4("Patient Risk Summary", class = "tab-title"),
          uiOutput("risk_card"),
          hr(),
          h5("IPSS Severity Classification"),
          tableOutput("ipss_class_table"),
          hr(),
          uiOutput("profile_summary_boxes")
        )
      )
    ),

    # --------------------------------------------------------
    # TAB 2: Drug PK
    # --------------------------------------------------------
    tabPanel(
      "Drug PK",
      fluidRow(
        column(3,
          wellPanel(
            h4("PK Settings", class = "tab-title"),
            selectInput("pk_drug", "Select Drug:",
              choices = c(
                "Tamsulosin 0.4mg"  = "tamsulosin",
                "Finasteride 5mg"   = "finasteride",
                "Dutasteride 0.5mg" = "dutasteride",
                "Tadalafil 5mg"     = "tadalafil"
              ),
              selected = "tamsulosin"
            ),
            radioButtons("pk_timerange", "Time Range:",
              choices  = c("24h Steady State" = "ss24h", "1-Year Profile" = "yearly"),
              selected = "ss24h"
            ),
            hr(),
            h5("PK Parameters"),
            tableOutput("pk_param_table")
          )
        ),
        column(9,
          h4("Concentration–Time Profile", class = "tab-title"),
          plotlyOutput("pk_plot", height = "380px"),
          br(),
          p(style = "font-size: 12px; color: #666;",
            "Concentration simulated using one-compartment model with first-order absorption.",
            "Steady-state achieved by superposition of prior doses (QD dosing interval = 24h).",
            "Concentrations are total (protein-bound + free).")
        )
      )
    ),

    # --------------------------------------------------------
    # TAB 3: Hormonal Biomarkers
    # --------------------------------------------------------
    tabPanel(
      "Hormonal Biomarkers",
      fluidRow(
        column(3,
          wellPanel(
            h4("Biomarker Settings", class = "tab-title"),
            checkboxGroupInput("bm_select", "Biomarkers to Display:",
              choices = c(
                "DHT (Plasma, nmol/L)"     = "dht_plasma",
                "DHT (Prostate, nmol/g)"   = "dht_prostate",
                "Testosterone (nmol/L)"    = "testosterone",
                "PSA (ng/mL)"              = "psa"
              ),
              selected = c("dht_plasma", "psa")
            ),
            selectInput("bm_treatment", "Treatment Scenario:",
              choices = c(
                "Watchful Waiting"       = "watchful_waiting",
                "Tamsulosin 0.4mg QD"    = "tamsulosin",
                "Finasteride 5mg QD"     = "finasteride",
                "Dutasteride 0.5mg QD"   = "dutasteride",
                "Combination (Dut+Tams)" = "combination",
                "Tadalafil 5mg QD"       = "tadalafil"
              ),
              selected = "finasteride"
            ),
            sliderInput("bm_weeks", "Time Horizon (weeks):",
              min = 4, max = 104, value = 104, step = 4
            )
          )
        ),
        column(9,
          h4("Biomarker Trajectories", class = "tab-title"),
          plotlyOutput("bm_plot", height = "380px"),
          br(),
          h5("% Change from Baseline at Key Timepoints"),
          tableOutput("bm_change_table")
        )
      )
    ),

    # --------------------------------------------------------
    # TAB 4: Urodynamic Endpoints
    # --------------------------------------------------------
    tabPanel(
      "Urodynamic Endpoints",
      fluidRow(
        column(3,
          wellPanel(
            h4("Endpoint Settings", class = "tab-title"),
            checkboxGroupInput("uro_scenarios", "Display Scenarios:",
              choices = c(
                "Watchful Waiting"       = "watchful_waiting",
                "Tamsulosin 0.4mg QD"    = "tamsulosin",
                "Finasteride 5mg QD"     = "finasteride",
                "Dutasteride 0.5mg QD"   = "dutasteride",
                "Combination (Dut+Tams)" = "combination",
                "Tadalafil 5mg QD"       = "tadalafil"
              ),
              selected = c("watchful_waiting", "tamsulosin", "finasteride",
                           "dutasteride", "combination", "tadalafil")
            ),
            hr(),
            p(style = "font-size: 12px; color: #444;",
              strong("Reference values:"),
              br(), "IPSS <8: Mild symptoms",
              br(), "IPSS 8-19: Moderate symptoms",
              br(), "IPSS ≥20: Severe symptoms",
              br(), "Qmax ≥15 mL/s: Normal flow",
              br(), "PVR <50 mL: Normal",
              br(), "PVR >200 mL: Clinically significant"
            )
          )
        ),
        column(9,
          fluidRow(
            column(6,
              h5("IPSS Score (0-104 weeks)"),
              plotlyOutput("uro_ipss", height = "280px")
            ),
            column(6,
              h5("Qmax – Maximum Flow Rate (mL/s)"),
              plotlyOutput("uro_qmax", height = "280px")
            )
          ),
          fluidRow(
            column(6,
              h5("Post-Void Residual Volume (mL)"),
              plotlyOutput("uro_pvr", height = "280px")
            ),
            column(6,
              h5("Prostate Volume (cc)"),
              plotlyOutput("uro_pv", height = "280px")
            )
          )
        )
      )
    ),

    # --------------------------------------------------------
    # TAB 5: Treatment Scenarios Comparison
    # --------------------------------------------------------
    tabPanel(
      "Treatment Scenarios",
      fluidRow(
        column(3,
          wellPanel(
            h4("Comparison Settings", class = "tab-title"),
            radioButtons("scenario_endpoint", "Primary Endpoint (Bar Chart):",
              choices = c(
                "ΔIPSS at 52 weeks"    = "delta_ipss",
                "ΔQmax at 52 weeks"    = "delta_qmax",
                "ΔPV% at 52 weeks"     = "delta_pv_pct",
                "ΔPVR at 52 weeks"     = "delta_pvr",
                "PSA Change %"         = "delta_psa_pct"
              ),
              selected = "delta_ipss"
            ),
            hr(),
            p(style = "font-size: 12px; color: #555;",
              strong("Clinical Trial References:"),
              br(), em("MTOPS"), " (McConnell 2003): Combination ΔIPSS = -6.6, Tams = -4.4, Fin = -5.1 (4.5yr)",
              br(), em("CombAT"), " (Roehrborn 2010): Dut+Tams ΔIPSS ≈ -6 (2yr)",
              br(), em("NEPTUNE"), " (Oelke 2012): Tadalafil ΔIPSS = -5.6 (12wk)"
            )
          )
        ),
        column(9,
          h4("Endpoint Comparison at 52 Weeks", class = "tab-title"),
          plotlyOutput("scenario_bar", height = "360px"),
          br(),
          h5("Full Comparison Table (All Scenarios × Endpoints)"),
          DT::dataTableOutput("scenario_table")
        )
      )
    ),

    # --------------------------------------------------------
    # TAB 6: Disease Progression & Risk
    # --------------------------------------------------------
    tabPanel(
      "Disease Progression & Risk",
      fluidRow(
        column(3,
          wellPanel(
            h4("Progression Settings", class = "tab-title"),
            sliderInput("risk_years", "Simulation Horizon (years):",
              min = 1, max = 10, value = 5, step = 1
            ),
            selectInput("risk_treatment", "Treatment Strategy:",
              choices = c(
                "Watchful Waiting"       = "watchful_waiting",
                "Tamsulosin 0.4mg QD"    = "tamsulosin",
                "Finasteride 5mg QD"     = "finasteride",
                "Dutasteride 0.5mg QD"   = "dutasteride",
                "Combination (Dut+Tams)" = "combination",
                "Tadalafil 5mg QD"       = "tadalafil"
              ),
              selected = "combination"
            ),
            hr(),
            h5("5-Year Projected Outcomes"),
            uiOutput("risk_value_boxes")
          )
        ),
        column(9,
          fluidRow(
            column(6,
              h5("Cumulative AUR Risk (%)"),
              plotlyOutput("risk_aur_plot", height = "280px")
            ),
            column(6,
              h5("Cumulative TURP Risk (%)"),
              plotlyOutput("risk_turp_plot", height = "280px")
            )
          ),
          fluidRow(
            column(12,
              h5("IPSS Trajectory with Risk Stratification"),
              plotlyOutput("risk_ipss_plot", height = "280px")
            )
          ),
          br(),
          h5("Risk Stratification Summary Table"),
          tableOutput("risk_summary_table")
        )
      )
    )

  ) # end tabsetPanel
) # end fluidPage

# ============================================================
# SERVER
# ============================================================

server <- function(input, output, session) {

  # --------------------------------------------------------
  # TAB 1: Patient Profile — Reactive outputs
  # --------------------------------------------------------

  risk_category <- reactive({
    get_risk_category(input$baseline_pv, input$baseline_psa, input$baseline_ipss)
  })

  output$risk_card <- renderUI({
    rc  <- risk_category()
    cls <- switch(rc, High = "risk-high", Moderate = "risk-mod", Low = "risk-low")
    icon_txt <- switch(rc,
      High     = "High Progression Risk — Early intervention recommended",
      Moderate = "Moderate Progression Risk — Active surveillance appropriate",
      Low      = "Low Progression Risk — Watchful waiting acceptable"
    )
    criteria <- paste0(
      "PV > 40 cc: ",  ifelse(input$baseline_pv  > 40,  "Yes", "No"), " | ",
      "PSA > 1.5: ",   ifelse(input$baseline_psa > 1.5, "Yes", "No"), " | ",
      "IPSS > 20: ",   ifelse(input$baseline_ipss > 20, "Yes", "No")
    )
    div(class = paste("risk-card", cls),
      strong(paste0("Risk Category: ", rc)), br(),
      icon_txt, br(),
      tags$small(criteria)
    )
  })

  output$ipss_class_table <- renderTable({
    data.frame(
      Category  = c("Mild", "Moderate", "Severe"),
      IPSS      = c("0–7", "8–19", "20–35"),
      Action    = c("Watchful waiting", "Medication review", "Intervention")
    )
  }, striped = TRUE, bordered = TRUE, hover = TRUE)

  output$profile_summary_boxes <- renderUI({
    pvr_est <- round(pmax(10, 350 - 20 * input$baseline_qmax))
    div(
      div(class = "value-box val-primary",
        strong("Estimated PVR: "), paste0(pvr_est, " mL")
      ),
      div(class = "value-box val-warning",
        strong("Qmax/Normal ratio: "), paste0(round(input$baseline_qmax / 15 * 100), "%")
      ),
      div(class = "value-box val-danger",
        strong("PSA Density: "), paste0(round(input$baseline_psa / input$baseline_pv, 3), " ng/mL/cc")
      ),
      div(class = "value-box val-success",
        strong("Comorbidities: "),
        if (length(input$comorbidities) == 0) "None" else paste(input$comorbidities, collapse = ", ")
      )
    )
  })

  # --------------------------------------------------------
  # TAB 2: Drug PK — Reactive
  # --------------------------------------------------------

  pk_data_reactive <- reactive({
    drug <- input$pk_drug
    if (input$pk_timerange == "ss24h") {
      simulate_pk_ss(drug, tau = 24, time_points = seq(0, 24, by = 0.5))
    } else {
      simulate_pk_yearly(drug)
    }
  })

  output$pk_plot <- renderPlotly({
    drug <- input$pk_drug
    df   <- pk_data_reactive()

    if (input$pk_timerange == "ss24h") {
      p <- ggplot(df, aes(x = time, y = concentration)) +
        geom_line(color = "#2980b9", linewidth = 1.2) +
        geom_area(fill = "#2980b9", alpha = 0.15) +
        labs(
          title = paste0(tools::toTitleCase(drug), " — 24h Steady-State PK Profile"),
          x     = "Time (hours)",
          y     = "Plasma Concentration (ng/mL)"
        ) +
        theme_bw(base_size = 13) +
        theme(plot.title = element_text(face = "bold"))
    } else {
      p <- ggplot(df, aes(x = day, y = cmax)) +
        geom_line(color = "#27ae60", linewidth = 1.0) +
        geom_hline(yintercept = max(df$cmax) * 0.95,
                   linetype = "dashed", color = "red", alpha = 0.6) +
        annotate("text", x = max(df$day) * 0.8, y = max(df$cmax) * 0.97,
                 label = "95% Steady State", size = 3.5, color = "red") +
        labs(
          title = paste0(tools::toTitleCase(drug), " — 1-Year Daily Cmax Accumulation"),
          x     = "Day",
          y     = "Daily Cmax (ng/mL)"
        ) +
        theme_bw(base_size = 13) +
        theme(plot.title = element_text(face = "bold"))
    }

    ggplotly(p) |> layout(hovermode = "x unified")
  })

  output$pk_param_table <- renderTable({
    pk_table <- data.frame(
      Parameter = c("Dose", "Tmax", "t½", "Vd (L)", "CL (L/h)", "F (%)", "Prot. Bind."),
      Tamsulosin  = c("0.4 mg", "6 h",    "10 h",   "15",   "1.8",  "90", "94%"),
      Finasteride = c("5 mg",   "~2 h",   "6 h",    "76",   "5.0",  "63", "90%"),
      Dutasteride = c("0.5 mg", "~2.5 h", "~5 wk",  "5",    "0.15", "60", "99%"),
      Tadalafil   = c("5 mg",   "2 h",    "17 h",   "25",   "2.5",  "80", "94%")
    )
    # Highlight selected drug (not easy in renderTable; just return full table)
    pk_table
  }, striped = TRUE, bordered = TRUE)

  # --------------------------------------------------------
  # TAB 3: Hormonal Biomarkers — Reactive
  # --------------------------------------------------------

  bm_data <- reactive({
    weeks <- seq(0, input$bm_weeks, by = 2)
    simulate_biomarkers(input$bm_treatment, input$baseline_psa, weeks)
  })

  output$bm_plot <- renderPlotly({
    df      <- bm_data()
    markers <- input$bm_select
    if (length(markers) == 0) markers <- "dht_plasma"

    bm_labels <- c(
      dht_plasma   = "DHT Plasma (nmol/L)",
      dht_prostate = "DHT Prostate (nmol/g)",
      testosterone = "Testosterone (nmol/L)",
      psa          = "PSA (ng/mL)"
    )
    bm_colors <- c(
      dht_plasma   = "#e74c3c",
      dht_prostate = "#c0392b",
      testosterone = "#2980b9",
      psa          = "#27ae60"
    )

    fig <- plot_ly()

    for (bm in markers) {
      fig <- fig |>
        add_trace(
          data = df,
          x    = ~time_weeks,
          y    = as.formula(paste0("~", bm)),
          type = "scatter",
          mode = "lines",
          name = bm_labels[bm],
          line = list(color = bm_colors[bm], width = 2)
        )
    }

    # Add reference lines
    if ("dht_plasma" %in% markers) {
      fig <- fig |>
        add_trace(x = c(0, max(df$time_weeks)), y = c(0.9, 0.9),
                  type = "scatter", mode = "lines", name = "DHT lower normal",
                  line = list(dash = "dot", color = "#e74c3c", width = 1), showlegend = TRUE) |>
        add_trace(x = c(0, max(df$time_weeks)), y = c(2.5, 2.5),
                  type = "scatter", mode = "lines", name = "DHT upper normal",
                  line = list(dash = "dot", color = "#e74c3c", width = 1), showlegend = TRUE)
    }
    if ("psa" %in% markers) {
      fig <- fig |>
        add_trace(x = c(0, max(df$time_weeks)), y = c(4, 4),
                  type = "scatter", mode = "lines", name = "PSA threshold 4",
                  line = list(dash = "dash", color = "#27ae60", width = 1.5), showlegend = TRUE)
    }

    fig |> layout(
      title  = paste("Hormonal Biomarkers —", treatment_labels[input$bm_treatment]),
      xaxis  = list(title = "Weeks from Baseline"),
      yaxis  = list(title = "Concentration"),
      hovermode = "x unified",
      legend = list(orientation = "h", y = -0.2)
    )
  })

  output$bm_change_table <- renderTable({
    df <- bm_data()
    timepoints <- c(13, 26, 52, 104)  # weeks (3, 6, 12, 24 months)
    labels     <- c("3 months", "6 months", "12 months", "24 months")

    # Baseline (t=0)
    base <- df[df$time_weeks == 0, , drop = FALSE]
    if (nrow(base) == 0) {
      base_dht_p <- 1.8; base_dht_pr <- 25; base_testo <- 15; base_psa <- input$baseline_psa
    } else {
      base_dht_p  <- base$dht_plasma[1]
      base_dht_pr <- base$dht_prostate[1]
      base_testo  <- base$testosterone[1]
      base_psa    <- base$psa[1]
    }

    rows <- lapply(seq_along(timepoints), function(j) {
      tp   <- timepoints[j]
      near <- df[which.min(abs(df$time_weeks - tp)), ]
      data.frame(
        Timepoint        = labels[j],
        `DHT Plasma %`   = paste0(round((near$dht_plasma   / base_dht_p  - 1) * 100, 1), "%"),
        `DHT Prostate %` = paste0(round((near$dht_prostate / base_dht_pr - 1) * 100, 1), "%"),
        `Testosterone %` = paste0(round((near$testosterone / base_testo  - 1) * 100, 1), "%"),
        `PSA %`          = paste0(round((near$psa          / base_psa    - 1) * 100, 1), "%"),
        stringsAsFactors = FALSE, check.names = FALSE
      )
    })
    do.call(rbind, rows)
  }, striped = TRUE, bordered = TRUE)

  # --------------------------------------------------------
  # TAB 4: Urodynamic Endpoints — Reactive
  # --------------------------------------------------------

  uro_all_data <- reactive({
    weeks  <- seq(0, 104, by = 2)
    scens  <- input$uro_scenarios
    if (length(scens) == 0) scens <- "watchful_waiting"

    pv0   <- input$baseline_pv
    ipss0 <- input$baseline_ipss
    qmax0 <- input$baseline_qmax

    dfs <- lapply(scens, function(s) {
      simulate_pd(s, pv0, ipss0, qmax0, weeks)
    })
    bind_rows(dfs)
  })

  make_uro_plot <- function(df, y_var, y_lab, ref_lines = NULL) {
    scens_sel <- unique(df$treatment)
    colors    <- treatment_colors[scens_sel]

    fig <- plot_ly()

    for (sc in scens_sel) {
      dfs <- df[df$treatment == sc, ]
      fig <- fig |>
        add_trace(
          x    = dfs$time_weeks,
          y    = dfs[[y_var]],
          type = "scatter",
          mode = "lines",
          name = treatment_labels[sc],
          line = list(color = colors[sc], width = 2)
        )
    }

    # Add reference horizontal lines
    if (!is.null(ref_lines)) {
      for (rl in ref_lines) {
        fig <- fig |>
          add_trace(
            x    = c(0, 104),
            y    = c(rl$val, rl$val),
            type = "scatter",
            mode = "lines",
            name = rl$label,
            line = list(dash = "dot", color = rl$color, width = 1.5),
            showlegend = TRUE
          )
      }
    }

    fig |> layout(
      xaxis     = list(title = "Weeks"),
      yaxis     = list(title = y_lab),
      hovermode = "x unified",
      legend    = list(orientation = "h", y = -0.3)
    )
  }

  output$uro_ipss <- renderPlotly({
    df <- uro_all_data()
    make_uro_plot(df, "ipss", "IPSS Score",
      ref_lines = list(
        list(val = 8,  label = "Mild/Moderate threshold (8)", color = "#f39c12"),
        list(val = 20, label = "Moderate/Severe threshold (20)", color = "#e74c3c")
      )
    )
  })

  output$uro_qmax <- renderPlotly({
    df <- uro_all_data()
    make_uro_plot(df, "qmax", "Qmax (mL/s)",
      ref_lines = list(
        list(val = 15, label = "Normal Qmax (15 mL/s)", color = "#27ae60")
      )
    )
  })

  output$uro_pvr <- renderPlotly({
    df <- uro_all_data()
    make_uro_plot(df, "pvr", "Post-Void Residual (mL)",
      ref_lines = list(
        list(val = 50,  label = "Normal threshold (50 mL)", color = "#27ae60"),
        list(val = 200, label = "Clinically significant (200 mL)", color = "#e74c3c")
      )
    )
  })

  output$uro_pv <- renderPlotly({
    df <- uro_all_data()
    make_uro_plot(df, "pv", "Prostate Volume (cc)",
      ref_lines = list(
        list(val = 30, label = "Normal prostate (~30 cc)", color = "#27ae60"),
        list(val = 40, label = "Enlarged threshold (40 cc)", color = "#f39c12")
      )
    )
  })

  # --------------------------------------------------------
  # TAB 5: Treatment Scenarios Comparison — Reactive
  # --------------------------------------------------------

  scenario_comparison_data <- reactive({
    all_trt <- c("watchful_waiting", "tamsulosin", "finasteride",
                 "dutasteride", "combination", "tadalafil")
    week52  <- 52

    pv0   <- input$baseline_pv
    ipss0 <- input$baseline_ipss
    qmax0 <- input$baseline_qmax
    psa0  <- input$baseline_psa

    rows <- lapply(all_trt, function(trt) {
      pd  <- simulate_pd(trt, pv0, ipss0, qmax0, week52)
      bm  <- simulate_biomarkers(trt, psa0, week52)

      delta_ipss  <- pd$ipss[1]   - ipss0
      delta_qmax  <- pd$qmax[1]   - qmax0
      delta_pv_pct <- (pd$pv[1]   - pv0) / pv0 * 100
      delta_pvr   <- pd$pvr[1]    - (350 - 20 * qmax0)  # vs estimated baseline
      delta_psa_pct <- (bm$psa[1] - psa0) / psa0 * 100

      data.frame(
        treatment    = trt,
        delta_ipss   = delta_ipss,
        delta_qmax   = delta_qmax,
        delta_pv_pct = delta_pv_pct,
        delta_pvr    = delta_pvr,
        delta_psa_pct = delta_psa_pct,
        stringsAsFactors = FALSE
      )
    })
    bind_rows(rows)
  })

  output$scenario_bar <- renderPlotly({
    df  <- scenario_comparison_data()
    ep  <- input$scenario_endpoint

    y_labs <- c(
      delta_ipss    = "ΔIPSS (points, lower = better)",
      delta_qmax    = "ΔQmax (mL/s, higher = better)",
      delta_pv_pct  = "ΔProstate Volume (%)",
      delta_pvr     = "ΔPVR (mL)",
      delta_psa_pct = "ΔPSA (%)"
    )

    df$label <- treatment_labels[df$treatment]
    df$color <- treatment_colors[df$treatment]

    # Reference trial data points
    trial_data <- NULL
    if (ep == "delta_ipss") {
      trial_data <- data.frame(
        label = c("MTOPS: Combination (4.5yr)", "MTOPS: Finasteride (4.5yr)",
                  "MTOPS: Tamsulosin (4.5yr)", "CombAT: Dut+Tams (2yr)",
                  "NEPTUNE: Tadalafil (12wk)"),
        value = c(-6.6, -5.1, -4.4, -6.0, -5.6),
        trial = c("MTOPS", "MTOPS", "MTOPS", "CombAT", "NEPTUNE")
      )
    }

    fig <- plot_ly(df,
      x    = ~label,
      y    = as.formula(paste0("~", ep)),
      type = "bar",
      marker = list(color = ~color, line = list(width = 1, color = "white")),
      text   = ~paste0(round(get(ep), 2)),
      textposition = "outside",
      name   = "Model Prediction"
    )

    # Overlay trial reference points
    if (!is.null(trial_data)) {
      fig <- fig |>
        add_trace(
          data = trial_data,
          x    = ~label,
          y    = ~value,
          type = "scatter",
          mode = "markers+text",
          marker = list(size = 12, color = "black", symbol = "diamond"),
          text   = ~trial,
          textposition = "top right",
          name   = "Clinical Trial Data",
          inherit = FALSE
        )
    }

    fig |> layout(
      title  = paste("Treatment Comparison at 52 Weeks:", y_labs[ep]),
      xaxis  = list(title = "Treatment"),
      yaxis  = list(title = y_labs[ep]),
      barmode = "group",
      legend = list(orientation = "h", y = -0.25)
    )
  })

  output$scenario_table <- DT::renderDataTable({
    df <- scenario_comparison_data()
    df$Treatment   <- treatment_labels[df$treatment]
    df$`ΔIPSS`     <- round(df$delta_ipss, 2)
    df$`ΔQmax`     <- round(df$delta_qmax, 2)
    df$`ΔPV%`      <- paste0(round(df$delta_pv_pct, 1), "%")
    df$`ΔPVR (mL)` <- round(df$delta_pvr, 1)
    df$`ΔPSA%`     <- paste0(round(df$delta_psa_pct, 1), "%")

    df_out <- df[, c("Treatment", "ΔIPSS", "ΔQmax", "ΔPV%", "ΔPVR (mL)", "ΔPSA%")]

    DT::datatable(
      df_out,
      rownames  = FALSE,
      options   = list(pageLength = 6, dom = "t"),
      class     = "display nowrap"
    )
  })

  # --------------------------------------------------------
  # TAB 6: Disease Progression & Risk — Reactive
  # --------------------------------------------------------

  risk_data_all <- reactive({
    all_trt <- c("watchful_waiting", "tamsulosin", "finasteride",
                 "dutasteride", "combination", "tadalafil")
    yrs   <- input$risk_years
    pv0   <- input$baseline_pv
    psa0  <- input$baseline_psa
    age0  <- input$age
    ipss0 <- input$baseline_ipss

    dfs <- lapply(all_trt, function(trt) {
      rd <- compute_risk(pv0, psa0, age0, ipss0, trt, yrs)
      rd$treatment <- trt
      rd
    })
    bind_rows(dfs)
  })

  output$risk_aur_plot <- renderPlotly({
    df     <- risk_data_all()
    sel_tx <- input$risk_treatment

    fig <- plot_ly()
    for (trt in unique(df$treatment)) {
      dft <- df[df$treatment == trt, ]
      lw  <- ifelse(trt == sel_tx, 3, 1.2)
      fig <- fig |>
        add_trace(
          x    = dft$year,
          y    = dft$aur_risk,
          type = "scatter",
          mode = "lines",
          name = treatment_labels[trt],
          line = list(color = treatment_colors[trt], width = lw)
        )
    }
    fig |> layout(
      xaxis = list(title = "Years"),
      yaxis = list(title = "Cumulative AUR Risk (%)"),
      hovermode = "x unified",
      legend    = list(orientation = "h", y = -0.35)
    )
  })

  output$risk_turp_plot <- renderPlotly({
    df     <- risk_data_all()
    sel_tx <- input$risk_treatment

    fig <- plot_ly()
    for (trt in unique(df$treatment)) {
      dft <- df[df$treatment == trt, ]
      lw  <- ifelse(trt == sel_tx, 3, 1.2)
      fig <- fig |>
        add_trace(
          x    = dft$year,
          y    = dft$turp_risk,
          type = "scatter",
          mode = "lines",
          name = treatment_labels[trt],
          line = list(color = treatment_colors[trt], width = lw)
        )
    }
    fig |> layout(
      xaxis = list(title = "Years"),
      yaxis = list(title = "Cumulative TURP Risk (%)"),
      hovermode = "x unified",
      legend    = list(orientation = "h", y = -0.35)
    )
  })

  output$risk_ipss_plot <- renderPlotly({
    sel_tx <- input$risk_treatment
    yrs    <- input$risk_years
    weeks  <- seq(0, yrs * 52, by = 4)
    pv0    <- input$baseline_pv
    ipss0  <- input$baseline_ipss
    qmax0  <- input$baseline_qmax

    fig <- plot_ly()

    # Shading for severity bands
    fig <- fig |>
      add_trace(
        x    = c(0, max(weeks), max(weeks), 0),
        y    = c(0, 0, 7, 7),
        type = "scatter", mode = "none",
        fill = "toself", fillcolor = "rgba(40,167,69,0.12)",
        name = "Mild (0-7)", showlegend = TRUE
      ) |>
      add_trace(
        x    = c(0, max(weeks), max(weeks), 0),
        y    = c(8, 8, 19, 19),
        type = "scatter", mode = "none",
        fill = "toself", fillcolor = "rgba(255,193,7,0.12)",
        name = "Moderate (8-19)", showlegend = TRUE
      ) |>
      add_trace(
        x    = c(0, max(weeks), max(weeks), 0),
        y    = c(20, 20, 35, 35),
        type = "scatter", mode = "none",
        fill = "toself", fillcolor = "rgba(220,53,69,0.12)",
        name = "Severe (20-35)", showlegend = TRUE
      )

    for (trt in c("watchful_waiting", sel_tx)) {
      pd  <- simulate_pd(trt, pv0, ipss0, qmax0, weeks)
      lw  <- ifelse(trt == sel_tx, 3, 1.5)
      lst <- ifelse(trt == sel_tx, "solid", "dot")
      fig <- fig |>
        add_trace(
          x    = pd$time_weeks,
          y    = pd$ipss,
          type = "scatter",
          mode = "lines",
          name = treatment_labels[trt],
          line = list(color = treatment_colors[trt], width = lw, dash = lst)
        )
    }

    fig |> layout(
      xaxis     = list(title = "Weeks from Baseline"),
      yaxis     = list(title = "IPSS Score", range = c(0, 35)),
      hovermode = "x unified",
      legend    = list(orientation = "h", y = -0.3)
    )
  })

  output$risk_value_boxes <- renderUI({
    trt  <- input$risk_treatment
    rd5  <- compute_risk(input$baseline_pv, input$baseline_psa, input$age,
                         input$baseline_ipss, trt, 5)

    aur5  <- round(rd5$aur_risk[nrow(rd5)],  1)
    turp5 <- round(rd5$turp_risk[nrow(rd5)], 1)

    aur_class  <- ifelse(aur5 > 10, "val-danger", ifelse(aur5 > 5, "val-warning", "val-success"))
    turp_class <- ifelse(turp5 > 8, "val-danger", ifelse(turp5 > 4, "val-warning", "val-success"))

    div(
      div(class = paste("value-box", aur_class),
        strong("5-yr AUR Risk: "), paste0(aur5, "%")
      ),
      div(class = paste("value-box", turp_class),
        strong("5-yr TURP Risk: "), paste0(turp5, "%")
      ),
      div(class = "value-box val-primary",
        strong("Treatment: "), treatment_labels[trt]
      )
    )
  })

  output$risk_summary_table <- renderTable({
    all_trt <- c("watchful_waiting", "tamsulosin", "finasteride",
                 "dutasteride", "combination", "tadalafil")
    checkpoints <- c(1, 3, 5)

    rows <- lapply(all_trt, function(trt) {
      rd <- compute_risk(input$baseline_pv, input$baseline_psa, input$age,
                         input$baseline_ipss, trt, 5)

      vals_aur  <- sapply(checkpoints, function(y) {
        r <- rd[which.min(abs(rd$year - y)), ]
        paste0(round(r$aur_risk, 1), "%")
      })
      vals_turp <- sapply(checkpoints, function(y) {
        r <- rd[which.min(abs(rd$year - y)), ]
        paste0(round(r$turp_risk, 1), "%")
      })

      data.frame(
        Treatment   = treatment_labels[trt],
        `AUR yr1`   = vals_aur[1],
        `AUR yr3`   = vals_aur[2],
        `AUR yr5`   = vals_aur[3],
        `TURP yr1`  = vals_turp[1],
        `TURP yr3`  = vals_turp[2],
        `TURP yr5`  = vals_turp[3],
        stringsAsFactors = FALSE,
        check.names = FALSE
      )
    })
    do.call(rbind, rows)
  }, striped = TRUE, bordered = TRUE, hover = TRUE)

}  # end server

# ============================================================
# LAUNCH APP
# ============================================================

shinyApp(ui = ui, server = server)
