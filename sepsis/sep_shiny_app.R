# =============================================================================
# Sepsis QSP Model — Interactive Shiny Dashboard
# sep_shiny_app.R
#
# Quantitative Systems Pharmacology model of sepsis pathophysiology:
#   - Host-pathogen dynamics (bacteria, immune response)
#   - Cytokine cascade (TNF-α, IL-6, IL-10, IL-1β, HMGB1)
#   - Hemodynamics & organ dysfunction (MAP, SOFA)
#   - Antibiotic PK (2-compartment) for Pip/Taz, Vancomycin, Meropenem
#   - Vasopressor & corticosteroid effects
#   - 6-tab shinydashboard with dark theme
#
# References: Kumar 2006 (antibiotic timing), Reynolds 2006 (math model of sepsis),
#             Clermont 2004 (agent-based model), Madelain 2017 (PK/PD sepsis)
# =============================================================================

library(shiny)
library(shinydashboard)
library(ggplot2)
library(deSolve)
library(dplyr)
library(tidyr)
library(scales)
library(gridExtra)
library(purrr)

# =============================================================================
# HELPER: ODE SYSTEM
# =============================================================================

sepsis_ode <- function(t, state, parms) {
  with(as.list(c(state, parms)), {

    # --- Bacterial dynamics ---
    # Hill function: antibiotic effect
    AB_effect_beta  <- (AB1_c / (EC50_abx_beta  + AB1_c))  # pip/taz or meropenem
    AB_effect_vanco <- (AB2_c / (EC50_abx_vanco + AB2_c))  # vancomycin
    AB_kill <- max(AB_effect_beta, AB_effect_vanco) * kill_max

    # TNF + IL6 activate phagocytic killing
    immune_kill <- kphi * (TNF / (TNF + K_TNF)) * (1 + IL6 / (IL6 + K_IL6))

    dB <- r_B * B * (1 - B / B_max) - (AB_kill + immune_kill) * B

    # --- Damage signal (tissue injury)
    dD <- kD * B - kdD * D

    # --- TNF-alpha dynamics
    sTNF <- k_TNF_B * B + k_TNF_D * D
    dTNF <- sTNF / (1 + (IL10 / K_IL10_TNF)^2) - kdTNF * TNF

    # --- IL-6 dynamics
    sIL6 <- k_IL6_B * B + k_IL6_TNF * TNF
    dIL6 <- sIL6 / (1 + (IL10 / K_IL10_IL6)^2) - kdIL6 * IL6

    # --- IL-1beta
    sIL1 <- k_IL1_B * B + k_IL1_TNF * TNF
    dIL1 <- sIL1 / (1 + (IL10 / K_IL10_IL6)^2) - kdIL1 * IL1

    # --- IL-10 (anti-inflammatory)
    sIL10 <- k_IL10_TNF * TNF + k_IL10_IL6 * IL6
    dIL10 <- sIL10 - kdIL10 * IL10

    # --- HMGB1 (late mediator)
    dHMGB1 <- k_HMGB1 * D - kdHMGB1 * HMGB1

    # --- Procalcitonin (PCT) — driven by bacterial burden
    dPCT <- k_PCT * B - kdPCT * PCT

    # --- CRP — driven by IL-6
    dCRP <- k_CRP * IL6 - kdCRP * CRP

    # --- Cortisol (endogenous + exogenous)
    # Endogenous follows circadian + stress response (simplified HPA)
    cort_endo  <- cort_base * (1 + stress_amp * TNF / (TNF + K_TNF_cort))
    # Total cortisol = endogenous + hydrocortisone (administered as constant)
    cortisol_total <- cort_endo + HC_dose / 100  # scale dose

    # --- MAP dynamics
    # Vasodilation from TNF/NO, partial restoration by NE + cortisol
    vasoD <- k_vasoD * (TNF / (TNF + K_TNF_MAP) + 0.5 * IL6 / (IL6 + K_IL6_MAP))
    NE_restore <- k_NE * NE_rate * (1 + cortisol_total / K_cort_MAP)
    dMAP <- -vasoD * (MAP - MAP_min) + NE_restore * (MAP_target - MAP) -
      kdMAP * (MAP - MAP_baseline)

    # --- Lactate (marker of hypoperfusion)
    # Rises when MAP low, falls with treatment
    hypoperfusion <- max(0, (MAP_target - MAP) / MAP_target)
    dLac <- k_Lac * hypoperfusion * B / (B + 1e3) + 0.05 * hypoperfusion -
      kdLac * Lac * cortisol_total / (cortisol_total + 0.5)

    # --- SOFA sub-scores (continuous approximation 0-4 scale each)
    # Respiratory (PaO2/FiO2)
    PF_ratio <- 400 - 200 * (TNF / (TNF + 5)) - 100 * (IL6 / (IL6 + 20))
    PF_ratio <- max(80, min(400, PF_ratio))
    SOFA_resp <- ifelse(PF_ratio >= 400, 0,
                 ifelse(PF_ratio >= 300, 1,
                 ifelse(PF_ratio >= 200, 2,
                 ifelse(PF_ratio >= 100, 3, 4))))

    # Coagulation (platelet)
    Plt <- 200 - 100 * (D / (D + 2)) - 50 * (HMGB1 / (HMGB1 + 1))
    Plt <- max(10, min(200, Plt))
    SOFA_coag <- ifelse(Plt >= 150, 0,
                 ifelse(Plt >= 100, 1,
                 ifelse(Plt >= 50,  2,
                 ifelse(Plt >= 20,  3, 4))))

    # Hepatic (bilirubin surrogate)
    bili_proxy <- 0.8 + 5 * D / (D + 3)
    SOFA_hep <- ifelse(bili_proxy < 1.2, 0,
                ifelse(bili_proxy < 2,   1,
                ifelse(bili_proxy < 6,   2,
                ifelse(bili_proxy < 12,  3, 4))))

    # Cardiovascular
    SOFA_cv <- ifelse(MAP >= 70, 0,
               ifelse(MAP >= 65, 1,
               ifelse(NE_rate < 0.1, 2,
               ifelse(NE_rate < 0.25, 3, 4))))

    # Renal (creatinine)
    scr_proxy <- scr_base * (1 + 1.5 * D / (D + 2))
    SOFA_ren <- ifelse(scr_proxy < 1.2, 0,
                ifelse(scr_proxy < 2,   1,
                ifelse(scr_proxy < 3.5, 2,
                ifelse(scr_proxy < 5,   3, 4))))

    # CNS (GCS surrogate)
    GCS_proxy <- 15 - 6 * D / (D + 3)
    SOFA_cns <- ifelse(GCS_proxy >= 15, 0,
                ifelse(GCS_proxy >= 13, 1,
                ifelse(GCS_proxy >= 10, 2,
                ifelse(GCS_proxy >= 6,  3, 4))))

    dSOFA <- (SOFA_resp + SOFA_coag + SOFA_hep + SOFA_cv + SOFA_ren + SOFA_cns - SOFA) / 0.5

    # --- WBC dynamics
    dWBC <- k_WBC_B * B + k_WBC_IL6 * IL6 - kdWBC * (WBC - WBC_base)

    # --- Platelet (DIC tracking state variable)
    dPlt_var <- -k_Plt_D * D * Plt - k_Plt_HMGB1 * HMGB1 * Plt +
      k_Plt_prod * (200 - Plt) / 200

    # --- INR (coagulation)
    dINR <- k_INR * D - kdINR * (INR - 1.0)

    # --- Bilirubin
    dBili <- k_Bili * D - kdBili * Bili

    # === ANTIBIOTIC PK (2-compartment) ===

    # Drug 1: beta-lactam (pip/taz or meropenem)
    # Dosing is handled externally via event list; here just distribution/elimination
    dAB1_c <- -k12_ab1 * AB1_c + k21_ab1 * AB1_p - ke_ab1 * AB1_c
    dAB1_p <-  k12_ab1 * AB1_c - k21_ab1 * AB1_p

    # Drug 2: vancomycin
    dAB2_c <- -k12_ab2 * AB2_c + k21_ab2 * AB2_p - ke_ab2 * AB2_c
    dAB2_p <-  k12_ab2 * AB2_c - k21_ab2 * AB2_p

    # Drug 3: norepinephrine (fast PK, nearly instantaneous at infusion rate)
    # dNE_conc is negligible; NE_rate is direct parameter from UI

    list(c(dB, dD, dTNF, dIL6, dIL1, dIL10, dHMGB1, dPCT, dCRP,
           dMAP, dLac, dSOFA, dWBC, dPlt_var, dINR, dBili,
           dAB1_c, dAB1_p, dAB2_c, dAB2_p),
         PF_ratio    = PF_ratio,
         Plt_calc    = Plt,
         bili_proxy  = bili_proxy,
         SOFA_resp   = SOFA_resp,
         SOFA_coag   = SOFA_coag,
         SOFA_hep    = SOFA_hep,
         SOFA_cv     = SOFA_cv,
         SOFA_ren    = SOFA_ren,
         SOFA_cns    = SOFA_cns,
         cortisol    = cort_endo)
  })
}

# =============================================================================
# SIMULATION ENGINE
# =============================================================================

run_sepsis_sim <- function(
  age = 60, weight = 70, scr_base = 1.0,
  infection_type = "E. coli",
  comorbidity = "none",
  log_inoculum = 5,
  severity = "Sepsis",
  abx_choice = "Pip/Taz",
  abx_dose = 4.5,  # g per dose
  abx_interval = 8, # hours
  abx_delay = 0,    # hours delay to first dose
  NE_rate = 0.1,    # ug/kg/min
  HC_dose = 0,      # mg/day hydrocortisone
  toci_given = FALSE,
  t_end = 72,
  scenario_label = "Standard"
) {

  # --- Pathogen-specific growth/kill parameters
  bug_params <- switch(infection_type,
    "E. coli"         = list(r_B = 0.8, EC50_mult = 1.0, gram_pos = FALSE),
    "S. aureus"       = list(r_B = 0.7, EC50_mult = 1.5, gram_pos = TRUE),
    "K. pneumoniae"   = list(r_B = 0.9, EC50_mult = 2.0, gram_pos = FALSE),
    "P. aeruginosa"   = list(r_B = 1.0, EC50_mult = 2.5, gram_pos = FALSE),
    "Fungal"          = list(r_B = 0.5, EC50_mult = 5.0, gram_pos = FALSE),
    list(r_B = 0.8, EC50_mult = 1.0, gram_pos = FALSE)
  )

  # --- Comorbidity modifiers
  comor_mod <- switch(comorbidity,
    "diabetes"        = list(immune_mult = 0.7, scr_mult = 1.3),
    "CKD"             = list(immune_mult = 0.8, scr_mult = 2.0),
    "immunosuppressed"= list(immune_mult = 0.4, scr_mult = 1.0),
    list(immune_mult = 1.0, scr_mult = 1.0)
  )

  # --- Severity → initial state
  sev_init <- switch(severity,
    "SIRS"         = list(B0_mult = 0.3, MAP0 = 88,  Lac0 = 1.2),
    "Sepsis"       = list(B0_mult = 1.0, MAP0 = 80,  Lac0 = 2.0),
    "Severe Sepsis"= list(B0_mult = 3.0, MAP0 = 72,  Lac0 = 3.0),
    "Septic Shock" = list(B0_mult = 8.0, MAP0 = 58,  Lac0 = 5.0),
    list(B0_mult = 1.0, MAP0 = 80, Lac0 = 2.0)
  )

  # --- PK parameters per antibiotic
  # Drug 1 slot: beta-lactam (Pip/Taz or Meropenem); Drug 2: Vancomycin
  pk_beta <- switch(abx_choice,
    "Pip/Taz"   = list(Vd_c = 0.2*weight, Vd_p = 0.1*weight,
                        ke = 0.35, k12 = 0.8, k21 = 0.4, MIC = 16),
    "Meropenem" = list(Vd_c = 0.15*weight, Vd_p = 0.08*weight,
                        ke = 0.3, k12 = 0.6, k21 = 0.3, MIC = 2),
    "Vancomycin"= list(Vd_c = 0.25*weight, Vd_p = 0.4*weight,
                        ke = 0.08, k12 = 0.3, k21 = 0.15, MIC = 2),
    list(Vd_c = 0.2*weight, Vd_p = 0.1*weight,
         ke = 0.35, k12 = 0.8, k21 = 0.4, MIC = 16)
  )

  # Vancomycin always in slot 2 (separate concurrent drug)
  pk_vanco <- list(Vd_c = 0.25*weight, Vd_p = 0.4*weight,
                   ke = 0.08, k12 = 0.3, k21 = 0.15, MIC = 2)

  # --- Initial conditions
  B0    <- 10^log_inoculum * sev_init$B0_mult * bug_params$r_B
  MAP0  <- sev_init$MAP0
  Lac0  <- sev_init$Lac0
  WBC0  <- ifelse(severity == "SIRS", 14, ifelse(severity == "Sepsis", 18, 22))

  state0 <- c(
    B      = B0,
    D      = 0.1 * sev_init$B0_mult,
    TNF    = 2 * sev_init$B0_mult,
    IL6    = 5 * sev_init$B0_mult,
    IL1    = 1.5 * sev_init$B0_mult,
    IL10   = 1 * sev_init$B0_mult,
    HMGB1  = 0.1,
    PCT    = 0.5 * sev_init$B0_mult,
    CRP    = 5 * sev_init$B0_mult,
    MAP    = MAP0,
    Lac    = Lac0,
    SOFA   = ifelse(severity == "SIRS", 1,
             ifelse(severity == "Sepsis", 4,
             ifelse(severity == "Severe Sepsis", 8, 12))),
    WBC    = WBC0,
    Plt_var= 180,
    INR    = 1.1,
    Bili   = 0.8,
    AB1_c  = 0,
    AB1_p  = 0,
    AB2_c  = 0,
    AB2_p  = 0
  )

  # --- Parameters
  parms <- list(
    # Bacteria
    r_B         = bug_params$r_B,
    B_max       = 1e8,
    kill_max    = 2.0 * bug_params$EC50_mult^(-1),
    EC50_abx_beta  = 8   * bug_params$EC50_mult,
    EC50_abx_vanco = 2   * bug_params$EC50_mult,
    kphi        = 0.6 * comor_mod$immune_mult,
    K_TNF       = 3,
    K_IL6       = 10,

    # Damage
    kD  = 0.01,
    kdD = 0.05,

    # Cytokines
    k_TNF_B  = 0.2,  k_TNF_D  = 0.3, kdTNF = 0.8,
    k_IL6_B  = 0.1,  k_IL6_TNF= 0.4, kdIL6 = 0.5,
    k_IL1_B  = 0.15, k_IL1_TNF= 0.3, kdIL1 = 0.6,
    k_IL10_TNF = 0.3, k_IL10_IL6 = 0.2, kdIL10 = 0.4,
    K_IL10_TNF = 2,  K_IL10_IL6 = 3,
    k_HMGB1 = 0.02, kdHMGB1 = 0.05,
    # Tocilizumab: block IL-6 signaling
    IL6_block   = ifelse(toci_given, 0.9, 0),

    # PCT, CRP
    k_PCT = 0.5, kdPCT = 0.15,
    k_CRP = 2.0, kdCRP = 0.05,

    # Cortisol / HPA
    cort_base   = 1.0,
    stress_amp  = 3.0,
    K_TNF_cort  = 5,
    HC_dose     = HC_dose,
    K_cort_MAP  = 0.5,

    # Hemodynamics
    MAP_min     = 40,
    MAP_target  = 65,
    MAP_baseline= MAP0,
    k_vasoD     = 0.05,
    k_NE        = ifelse(NE_rate > 0, 0.3, 0),
    NE_rate     = NE_rate,
    kdMAP       = 0.1,
    K_TNF_MAP   = 4,
    K_IL6_MAP   = 15,

    # Lactate
    k_Lac  = 0.5,
    kdLac  = 0.3,

    # WBC
    k_WBC_B  = 0.1, k_WBC_IL6 = 0.02, kdWBC = 0.05,
    WBC_base = 7,

    # Platelets, INR, bili
    k_Plt_D    = 0.02, k_Plt_HMGB1 = 0.05, k_Plt_prod = 0.5,
    k_INR  = 0.1, kdINR = 0.2,
    k_Bili = 0.05, kdBili = 0.1,

    # PK drug 1 (beta-lactam or vanco as primary)
    ke_ab1  = pk_beta$ke,
    k12_ab1 = pk_beta$k12,
    k21_ab1 = pk_beta$k21,
    Vd_c_ab1 = pk_beta$Vd_c,

    # PK drug 2 (vancomycin secondary)
    ke_ab2  = pk_vanco$ke,
    k12_ab2 = pk_vanco$k12,
    k21_ab2 = pk_vanco$k21,
    Vd_c_ab2 = pk_vanco$Vd_c,

    # For SOFA renal
    scr_base = scr_base * comor_mod$scr_mult,
    # NE_rate for SOFA_cv
    NE_rate_param = NE_rate
  )

  # --- Dosing events (IV bolus added to central compartment)
  dose_times_ab1 <- seq(abx_delay, t_end, by = abx_interval)
  dose_amount_ab1 <- abx_dose * 1000 / pk_beta$Vd_c  # mg/L in central compartment

  # Vancomycin: 15-20 mg/kg q8-12h (use 15 mg/kg q12h as default)
  dose_times_ab2 <- seq(0, t_end, by = 12)
  dose_amount_ab2 <- (15 * weight) / pk_vanco$Vd_c

  events_df <- rbind(
    data.frame(var = "AB1_c", time = dose_times_ab1, value = dose_amount_ab1, method = "add"),
    data.frame(var = "AB2_c", time = dose_times_ab2, value = dose_amount_ab2, method = "add")
  )
  events_df <- events_df[order(events_df$time), ]

  # --- Time vector
  times <- seq(0, t_end, by = 0.25)

  # --- Solve
  tryCatch({
    out <- deSolve::ode(
      y     = state0,
      times = times,
      func  = sepsis_ode,
      parms = parms,
      events = list(data = events_df),
      method = "lsoda"
    )

    df <- as.data.frame(out)
    df$scenario <- scenario_label
    df$MIC_ab1  <- pk_beta$MIC

    # Add PF_ratio, SOFA components as columns from output
    df
  }, error = function(e) {
    # Return empty df with correct structure on error
    data.frame(time = seq(0, t_end, by = 0.25), scenario = scenario_label)
  })
}

# =============================================================================
# PK METRICS
# =============================================================================

calc_pk_metrics <- function(sim_out, abx_choice, MIC) {
  if (!"AB1_c" %in% names(sim_out)) return(NULL)

  # AUC via trapezoidal rule over 24h
  t24 <- sim_out[sim_out$time <= 24, ]
  AUC24 <- tryCatch(
    sum(diff(t24$time) * (head(t24$AB1_c, -1) + tail(t24$AB1_c, -1)) / 2),
    error = function(e) NA
  )

  # %fT>MIC (free fraction: f=0.5 for beta-lactams, 0.1 for vancomycin)
  f_free <- ifelse(abx_choice == "Vancomycin", 0.1, 0.5)
  pct_fT_MIC <- mean(sim_out$AB1_c * f_free > MIC, na.rm = TRUE) * 100

  list(
    AUC24     = round(AUC24, 1),
    AUC_MIC   = round(AUC24 / MIC, 1),
    pct_fT_MIC= round(pct_fT_MIC, 1),
    Cmax      = round(max(sim_out$AB1_c, na.rm = TRUE), 2),
    Ctrough   = round(min(sim_out$AB1_c[sim_out$time > 6], na.rm = TRUE), 2)
  )
}

# =============================================================================
# MORTALITY SCORE (APACHE II-inspired)
# =============================================================================

calc_mortality <- function(age, SOFA_score, Lac, MAP, severity, comorbidity) {
  # Simplified 28-day mortality estimate
  apache_score <- 0
  # Age points
  apache_score <- apache_score + ifelse(age < 45, 0,
                              ifelse(age < 55, 2,
                              ifelse(age < 65, 3,
                              ifelse(age < 75, 5, 6))))
  # SOFA weight
  apache_score <- apache_score + SOFA_score * 1.5

  # Lactate
  apache_score <- apache_score + ifelse(Lac < 2, 0,
                              ifelse(Lac < 4, 3, 7))
  # Comorbidity
  apache_score <- apache_score + switch(comorbidity,
    "immunosuppressed" = 5,
    "CKD"              = 3,
    "diabetes"         = 2, 0)

  # Severity offset
  apache_score <- apache_score + switch(severity,
    "Septic Shock" = 8,
    "Severe Sepsis" = 4,
    "Sepsis" = 1, 0)

  # Logistic conversion
  p_death <- 1 / (1 + exp(-(apache_score * 0.08 - 2.5)))
  round(min(p_death * 100, 99), 1)
}

# =============================================================================
# GGPLOT THEME
# =============================================================================

theme_sep <- function(base_size = 11) {
  theme_minimal(base_size = base_size) +
    theme(
      plot.background    = element_rect(fill = "#1a1f2e", color = NA),
      panel.background   = element_rect(fill = "#1a1f2e", color = NA),
      panel.grid.major   = element_line(color = "#2d3548", size = 0.4),
      panel.grid.minor   = element_line(color = "#252a3a", size = 0.2),
      axis.text          = element_text(color = "#9aa3b8", size = 9),
      axis.title         = element_text(color = "#c8cfe0", size = 10),
      plot.title         = element_text(color = "#e8ecf4", size = 12, face = "bold"),
      plot.subtitle      = element_text(color = "#7a849e", size = 9),
      legend.background  = element_rect(fill = "#1a1f2e", color = NA),
      legend.text        = element_text(color = "#9aa3b8"),
      legend.title       = element_text(color = "#c8cfe0"),
      strip.text         = element_text(color = "#c8cfe0", face = "bold"),
      strip.background   = element_rect(fill = "#252a3a"),
      plot.margin        = margin(10, 15, 10, 10)
    )
}

PAL <- c(
  "#4e9af1",  # blue
  "#e05c5c",  # red
  "#52c4a0",  # teal
  "#f0a354",  # amber
  "#b47cf5",  # violet
  "#f5d44e"   # yellow
)

# =============================================================================
# UI
# =============================================================================

ui <- dashboardPage(
  skin = "black",

  dashboardHeader(
    title = tags$span(
      tags$b("Sepsis QSP", style = "color:#4e9af1; font-size:17px;"),
      tags$span(" | Clinical Simulator", style = "color:#7a849e; font-size:13px;")
    ),
    titleWidth = 300
  ),

  dashboardSidebar(
    width = 240,
    tags$head(
      tags$style(HTML("
        .skin-black .main-header .logo { background-color:#0d1117 !important; }
        .skin-black .main-header .navbar { background-color:#0d1117 !important; }
        .skin-black .main-sidebar { background-color:#111622 !important; }
        .skin-black .sidebar-menu > li > a { color:#8a93a8 !important; }
        .skin-black .sidebar-menu > li.active > a { color:#4e9af1 !important; border-left:3px solid #4e9af1; }
        .content-wrapper { background-color:#141824 !important; }
        .box { background:#1a1f2e !important; border-top:none !important; }
        .box-header { background:#1a1f2e !important; color:#c8cfe0 !important; }
        .box-title { color:#c8cfe0 !important; }
        .shiny-input-container label { color:#9aa3b8 !important; }
        .irs-bar, .irs-bar-edge { background:#4e9af1 !important; border-color:#4e9af1 !important; }
        .irs-single { background:#4e9af1 !important; }
        .selectize-input { background:#252a3a !important; color:#c8cfe0 !important; border-color:#3a4158 !important; }
        .selectize-dropdown { background:#252a3a !important; color:#c8cfe0 !important; }
        .btn-primary { background-color:#4e9af1 !important; border-color:#4e9af1 !important; }
        .value-box { border-radius:6px !important; }
        table.dataTable { background:#1a1f2e !important; color:#c8cfe0 !important; }
        .nav-tabs-custom > .tab-content { background:#1a1f2e !important; }
        .checkbox label { color:#9aa3b8 !important; }
        hr { border-color:#2d3548 !important; }
        .info-box { background:#252a3a !important; color:#c8cfe0 !important; }
        .info-box-icon { font-size:28px !important; line-height:55px !important; }
      "))
    ),
    sidebarMenu(
      menuItem("Patient Profile", tabName = "tab_patient",  icon = icon("user-md")),
      menuItem("PK",              tabName = "tab_pk",       icon = icon("pills")),
      menuItem("Cytokine Dynamics",tabName = "tab_cyto",   icon = icon("biohazard")),
      menuItem("Hemodynamics",    tabName = "tab_hemo",     icon = icon("heartbeat")),
      menuItem("Scenarios",       tabName = "tab_scenario", icon = icon("chart-bar")),
      menuItem("Biomarkers",      tabName = "tab_bio",      icon = icon("flask"))
    ),
    tags$hr(),
    tags$div(style = "padding:10px 15px;",
      tags$p(style = "color:#5c6478; font-size:11px;",
        "Sepsis QSP Model v1.0", tags$br(),
        "Based on Reynolds 2006,", tags$br(),
        "Kumar 2006, Clermont 2004"
      )
    )
  ),

  dashboardBody(
    tabItems(

      # -----------------------------------------------------------------------
      # TAB 1: PATIENT PROFILE
      # -----------------------------------------------------------------------
      tabItem(tabName = "tab_patient",
        fluidRow(
          box(title = "Demographics", width = 4, solidHeader = FALSE,
            sliderInput("age", "Age (years)", 18, 90, 60, step = 1),
            sliderInput("weight", "Weight (kg)", 40, 150, 70, step = 1),
            sliderInput("scr_base", "Baseline Creatinine (mg/dL)", 0.5, 5.0, 1.0, step = 0.1)
          ),
          box(title = "Infection Characteristics", width = 4, solidHeader = FALSE,
            selectInput("infection_type", "Pathogen",
              choices = c("E. coli", "S. aureus", "K. pneumoniae", "P. aeruginosa", "Fungal"),
              selected = "E. coli"),
            selectInput("comorbidity", "Comorbidity",
              choices = c("none", "diabetes", "CKD", "immunosuppressed"),
              selected = "none"),
            sliderInput("log_inoculum", "Initial Inoculum (log₁₀ CFU/mL)", 3, 9, 5, step = 0.5),
            selectInput("severity", "Sepsis Severity",
              choices = c("SIRS", "Sepsis", "Severe Sepsis", "Septic Shock"),
              selected = "Sepsis")
          ),
          box(title = "Risk Summary", width = 4, solidHeader = FALSE,
            valueBoxOutput("risk_score_box", width = 12),
            tableOutput("patient_summary_table")
          )
        ),
        fluidRow(
          box(title = "Pathogen Profile & Expected Response", width = 12, solidHeader = FALSE,
            tableOutput("pathogen_table")
          )
        )
      ),

      # -----------------------------------------------------------------------
      # TAB 2: PK
      # -----------------------------------------------------------------------
      tabItem(tabName = "tab_pk",
        fluidRow(
          box(title = "Antibiotic Settings", width = 4, solidHeader = FALSE,
            selectInput("abx_choice", "Primary Antibiotic",
              choices = c("Pip/Taz", "Meropenem", "Vancomycin"),
              selected = "Pip/Taz"),
            sliderInput("abx_dose", "Dose per Administration (g)", 1, 10, 4.5, step = 0.5),
            sliderInput("abx_interval", "Dosing Interval (hours)", 4, 24, 8, step = 4),
            sliderInput("abx_delay", "Delay to First Dose (hours)", 0, 12, 0, step = 1),
            tags$hr(),
            sliderInput("NE_rate", "Norepinephrine Rate (μg/kg/min)", 0, 1.0, 0.1, step = 0.05),
            sliderInput("HC_dose", "Hydrocortisone (mg/day)", 0, 400, 0, step = 50)
          ),
          box(title = "PK Concentration-Time (48h)", width = 8, solidHeader = FALSE,
            plotOutput("pk_plot", height = "320px")
          )
        ),
        fluidRow(
          column(4,
            infoBoxOutput("auc_mic_box",      width = 12),
            infoBoxOutput("ft_mic_box",       width = 12),
            infoBoxOutput("cmax_box",         width = 12)
          ),
          column(8,
            box(title = "PK/PD Target Attainment", width = 12, solidHeader = FALSE,
              plotOutput("pk_target_plot", height = "220px")
            )
          )
        )
      ),

      # -----------------------------------------------------------------------
      # TAB 3: CYTOKINE DYNAMICS
      # -----------------------------------------------------------------------
      tabItem(tabName = "tab_cyto",
        fluidRow(
          box(title = "Cytokine Settings", width = 3, solidHeader = FALSE,
            checkboxInput("toci_given", "Add Tocilizumab (IL-6 blocker)", FALSE),
            tags$hr(),
            tags$p(style = "color:#7a849e; font-size:11px;",
              "Tocilizumab blocks IL-6 receptor signaling,",
              "modulating downstream inflammatory cascade."
            )
          ),
          box(title = "Pro-Inflammatory Cytokines", width = 9, solidHeader = FALSE,
            plotOutput("cyto_pro_plot", height = "300px")
          )
        ),
        fluidRow(
          box(title = "Anti-Inflammatory & Late Mediators", width = 6, solidHeader = FALSE,
            plotOutput("cyto_anti_plot", height = "260px")
          ),
          box(title = "Acute Phase Reactants", width = 6, solidHeader = FALSE,
            plotOutput("cyto_apr_plot", height = "260px")
          )
        )
      ),

      # -----------------------------------------------------------------------
      # TAB 4: HEMODYNAMICS & ORGAN FUNCTION
      # -----------------------------------------------------------------------
      tabItem(tabName = "tab_hemo",
        fluidRow(
          box(title = "Mean Arterial Pressure (72h)", width = 6, solidHeader = FALSE,
            plotOutput("map_plot", height = "260px")
          ),
          box(title = "Lactate Trajectory (72h)", width = 6, solidHeader = FALSE,
            plotOutput("lactate_plot", height = "260px")
          )
        ),
        fluidRow(
          box(title = "SOFA Score Components", width = 6, solidHeader = FALSE,
            plotOutput("sofa_plot", height = "280px")
          ),
          box(title = "Organ-Specific Markers", width = 6, solidHeader = FALSE,
            plotOutput("organ_plot", height = "280px")
          )
        ),
        fluidRow(
          box(title = "DIC Assessment", width = 6, solidHeader = FALSE,
            plotOutput("dic_plot", height = "240px")
          ),
          box(title = "ARDS Monitoring (PaO₂/FiO₂)", width = 6, solidHeader = FALSE,
            plotOutput("ards_plot", height = "240px")
          )
        )
      ),

      # -----------------------------------------------------------------------
      # TAB 5: TREATMENT SCENARIOS
      # -----------------------------------------------------------------------
      tabItem(tabName = "tab_scenario",
        fluidRow(
          box(title = "Scenario Selection", width = 4, solidHeader = FALSE,
            checkboxGroupInput("scenarios_sel", "Scenarios to Compare",
              choices = list(
                "No treatment (control)"                    = "s1",
                "Early antibiotics only"                    = "s2",
                "Antibiotics + vasopressors"                = "s3",
                "Full sepsis bundle (ABx + NE + steroids)"  = "s4",
                "Delayed antibiotics (6h)"                  = "s5",
                "Refractory shock (high NE + VP + steroids)"= "s6"
              ),
              selected = c("s2", "s4", "s5")
            ),
            actionButton("run_scenarios", "Run Comparison",
              class = "btn-primary", style = "width:100%; margin-top:10px;")
          ),
          box(title = "Bacterial Clearance Comparison", width = 8, solidHeader = FALSE,
            plotOutput("scenario_bact_plot", height = "300px")
          )
        ),
        fluidRow(
          box(title = "MAP Comparison", width = 6, solidHeader = FALSE,
            plotOutput("scenario_map_plot", height = "240px")
          ),
          box(title = "SOFA Score Comparison", width = 6, solidHeader = FALSE,
            plotOutput("scenario_sofa_plot", height = "240px")
          )
        ),
        fluidRow(
          box(title = "28-Day Mortality Estimates", width = 12, solidHeader = FALSE,
            tableOutput("mortality_table")
          )
        )
      ),

      # -----------------------------------------------------------------------
      # TAB 6: BIOMARKERS & ENDPOINTS
      # -----------------------------------------------------------------------
      tabItem(tabName = "tab_bio",
        fluidRow(
          valueBoxOutput("lactate_box",   width = 3),
          valueBoxOutput("pct_box",       width = 3),
          valueBoxOutput("crp_box",       width = 3),
          valueBoxOutput("sofa_score_box",width = 3)
        ),
        fluidRow(
          box(title = "Biomarker Dashboard (72h Trajectories)", width = 8,
            solidHeader = FALSE,
            plotOutput("biomarker_panel", height = "380px")
          ),
          box(title = "Clinical Status", width = 4, solidHeader = FALSE,
            tableOutput("sofa_timeline_table"),
            tags$hr(),
            tags$p(style = "color:#7a849e; font-size:11px;",
              tags$b(style = "color:#c8cfe0;", "De-escalation Guidance"),
              tags$br(),
              textOutput("deescalation_text")
            )
          )
        ),
        fluidRow(
          box(title = "Survival Probability Curve", width = 6, solidHeader = FALSE,
            plotOutput("survival_plot", height = "260px")
          ),
          box(title = "qSOFA & Organ Trend", width = 6, solidHeader = FALSE,
            plotOutput("qsofa_plot", height = "260px")
          )
        )
      )
    )
  )
)

# =============================================================================
# SERVER
# =============================================================================

server <- function(input, output, session) {

  # ---------------------------------------------------------------------------
  # Reactive: base simulation (single scenario from Tab 1/2/3/4 settings)
  # ---------------------------------------------------------------------------
  sim_base <- reactive({
    run_sepsis_sim(
      age            = input$age,
      weight         = input$weight,
      scr_base       = input$scr_base,
      infection_type = input$infection_type,
      comorbidity    = input$comorbidity,
      log_inoculum   = input$log_inoculum,
      severity       = input$severity,
      abx_choice     = input$abx_choice,
      abx_dose       = input$abx_dose,
      abx_interval   = input$abx_interval,
      abx_delay      = input$abx_delay,
      NE_rate        = input$NE_rate,
      HC_dose        = input$HC_dose,
      toci_given     = input$toci_given,
      t_end          = 72,
      scenario_label = "Current Settings"
    )
  })

  # ---------------------------------------------------------------------------
  # Reactive: scenario simulations (Tab 5)
  # ---------------------------------------------------------------------------
  scenario_configs <- list(
    s1 = list(abx_dose = 0, NE_rate = 0,    HC_dose = 0,   delay = 0,  label = "No Treatment"),
    s2 = list(abx_dose = 4.5, NE_rate = 0,  HC_dose = 0,   delay = 0,  label = "Early ABx"),
    s3 = list(abx_dose = 4.5, NE_rate = 0.15, HC_dose = 0, delay = 0,  label = "ABx + NE"),
    s4 = list(abx_dose = 4.5, NE_rate = 0.15, HC_dose = 200, delay = 0, label = "Full Bundle"),
    s5 = list(abx_dose = 4.5, NE_rate = 0,  HC_dose = 0,   delay = 6,  label = "Delayed ABx (6h)"),
    s6 = list(abx_dose = 4.5, NE_rate = 0.5,HC_dose = 300, delay = 0,  label = "Refractory Shock")
  )

  sim_scenarios <- eventReactive(input$run_scenarios, {
    req(input$scenarios_sel)
    purrr::map_dfr(input$scenarios_sel, function(sc) {
      cfg <- scenario_configs[[sc]]
      run_sepsis_sim(
        age            = input$age,
        weight         = input$weight,
        scr_base       = input$scr_base,
        infection_type = input$infection_type,
        comorbidity    = input$comorbidity,
        log_inoculum   = input$log_inoculum,
        severity       = input$severity,
        abx_choice     = input$abx_choice,
        abx_dose       = cfg$abx_dose,
        abx_interval   = input$abx_interval,
        abx_delay      = cfg$delay,
        NE_rate        = cfg$NE_rate,
        HC_dose        = cfg$HC_dose,
        toci_given     = FALSE,
        t_end          = 72,
        scenario_label = cfg$label
      )
    })
  }, ignoreNULL = FALSE)

  # Note: eventReactive with ignoreNULL = FALSE handles initial state automatically.

  # ---------------------------------------------------------------------------
  # TAB 1 OUTPUTS
  # ---------------------------------------------------------------------------

  output$risk_score_box <- renderValueBox({
    sim  <- sim_base()
    sofa_t0 <- sim$SOFA[1]
    lac_t0  <- sim$Lac[1]
    mort    <- calc_mortality(
      age       = input$age,
      SOFA_score= sofa_t0,
      Lac       = lac_t0,
      MAP       = sim$MAP[1],
      severity  = input$severity,
      comorbidity = input$comorbidity
    )
    color <- ifelse(mort > 50, "red", ifelse(mort > 30, "yellow", "green"))
    valueBox(
      value    = paste0(mort, "%"),
      subtitle = "Estimated 28-day Mortality Risk",
      icon     = icon("exclamation-triangle"),
      color    = color
    )
  })

  output$patient_summary_table <- renderTable({
    data.frame(
      Parameter = c("Age", "Weight", "Baseline Cr", "Pathogen", "Comorbidity",
                    "Inoculum", "Severity"),
      Value = c(
        paste(input$age, "yr"),
        paste(input$weight, "kg"),
        paste(input$scr_base, "mg/dL"),
        input$infection_type,
        input$comorbidity,
        paste0("10^", input$log_inoculum, " CFU/mL"),
        input$severity
      )
    )
  }, striped = FALSE, bordered = FALSE,
  rownames = FALSE,
  align = "lr",
  width = "100%")

  output$pathogen_table <- renderTable({
    data.frame(
      Pathogen = c("E. coli", "S. aureus", "K. pneumoniae", "P. aeruginosa", "Fungal"),
      `Gram Stain`      = c("Neg", "Pos", "Neg", "Neg", "N/A"),
      `Growth Rate`     = c("Fast", "Moderate", "Fast", "Very Fast", "Slow"),
      `Preferred ABx`   = c("Pip/Taz or Meropenem", "Vancomycin", "Meropenem", "Meropenem", "Antifungal"),
      `Resistance Risk` = c("Low", "Moderate (MRSA)", "High (KPC)", "Very High", "Moderate"),
      `Mortality Risk`  = c("Moderate", "Moderate-High", "High", "High", "Very High"),
      check.names = FALSE
    )
  }, striped = FALSE, width = "100%")

  # ---------------------------------------------------------------------------
  # TAB 2 OUTPUTS (PK)
  # ---------------------------------------------------------------------------

  output$pk_plot <- renderPlot({
    sim <- sim_base()
    req("AB1_c" %in% names(sim))

    MIC_val <- sim$MIC_ab1[1]
    pk_long <- sim[sim$time <= 48, ] %>%
      dplyr::select(time, `Central (C1)` = AB1_c, `Peripheral (C2)` = AB1_p,
                    `Vancomycin (C1)` = AB2_c) %>%
      pivot_longer(-time, names_to = "Compartment", values_to = "Conc")

    ggplot(pk_long, aes(time, Conc, color = Compartment)) +
      geom_line(size = 0.9) +
      geom_hline(yintercept = MIC_val, linetype = "dashed",
                 color = "#e05c5c", size = 0.7, alpha = 0.8) +
      annotate("text", x = 1, y = MIC_val * 1.15,
               label = paste0("MIC = ", MIC_val, " mg/L"),
               color = "#e05c5c", size = 3.2, hjust = 0) +
      scale_color_manual(values = c("#4e9af1", "#52c4a0", "#f0a354")) +
      scale_x_continuous(breaks = seq(0, 48, 8)) +
      labs(title = paste(input$abx_choice, "PK Profile"),
           subtitle = paste0("Dose: ", input$abx_dose, "g q", input$abx_interval, "h | Weight: ", input$weight, "kg"),
           x = "Time (hours)", y = "Concentration (mg/L)",
           color = NULL) +
      theme_sep()
  }, bg = "#1a1f2e")

  output$pk_target_plot <- renderPlot({
    sim <- sim_base()
    if (!"AB1_c" %in% names(sim)) return(NULL)

    # %fT>MIC across dose range
    doses <- seq(1, 10, by = 0.5)
    f_free <- ifelse(input$abx_choice == "Vancomycin", 0.1, 0.5)
    MIC_val <- sim$MIC_ab1[1]

    # Quick approximation for target attainment curves
    ta_df <- data.frame(
      dose  = doses,
      ta_50 = pmin(100, 100 * (1 - exp(-doses / MIC_val * 0.4 * f_free))),
      ta_90 = pmin(100, 100 * (1 - exp(-doses / MIC_val * 0.25 * f_free)))
    )

    current_ta <- pmin(100, 100 * (1 - exp(-input$abx_dose / MIC_val * 0.4 * f_free)))

    ggplot(ta_df, aes(x = dose)) +
      geom_ribbon(aes(ymin = ta_90, ymax = ta_50), fill = "#4e9af1", alpha = 0.15) +
      geom_line(aes(y = ta_50, color = "50% target"), size = 0.9) +
      geom_line(aes(y = ta_90, color = "90% target"), size = 0.9, linetype = "dashed") +
      geom_vline(xintercept = input$abx_dose, color = "#52c4a0",
                 linetype = "dashed", size = 0.8) +
      geom_point(data = data.frame(x = input$abx_dose, y = current_ta),
                 aes(x = x, y = y), color = "#52c4a0", size = 3) +
      scale_color_manual(values = c("#4e9af1", "#b47cf5")) +
      scale_y_continuous(limits = c(0, 100)) +
      labs(title = "Target Attainment Probability", x = "Dose (g)", y = "% PTA", color = NULL) +
      theme_sep()
  }, bg = "#1a1f2e")

  output$auc_mic_box <- renderInfoBox({
    sim  <- sim_base()
    if (!"AB1_c" %in% names(sim)) return(infoBox("AUC/MIC", "N/A", icon = icon("chart-area")))
    MIC_val <- sim$MIC_ab1[1]
    metrics <- calc_pk_metrics(sim, input$abx_choice, MIC_val)
    target  <- ifelse(input$abx_choice == "Vancomycin", "≥400", "≥125")
    infoBox(
      title = "AUC₀₋₂₄/MIC",
      value = metrics$AUC_MIC,
      subtitle = paste("Target:", target),
      icon  = icon("chart-area"),
      color = ifelse(metrics$AUC_MIC >= 125, "green", "red"),
      fill  = TRUE
    )
  })

  output$ft_mic_box <- renderInfoBox({
    sim <- sim_base()
    if (!"AB1_c" %in% names(sim)) return(infoBox("%fT>MIC", "N/A", icon = icon("clock")))
    MIC_val <- sim$MIC_ab1[1]
    metrics <- calc_pk_metrics(sim, input$abx_choice, MIC_val)
    infoBox(
      title = "%fT>MIC",
      value = paste0(metrics$pct_fT_MIC, "%"),
      subtitle = "Target ≥40–50%",
      icon  = icon("clock"),
      color = ifelse(metrics$pct_fT_MIC >= 40, "green", "red"),
      fill  = TRUE
    )
  })

  output$cmax_box <- renderInfoBox({
    sim <- sim_base()
    if (!"AB1_c" %in% names(sim)) return(infoBox("Cmax", "N/A", icon = icon("arrow-up")))
    MIC_val <- sim$MIC_ab1[1]
    metrics <- calc_pk_metrics(sim, input$abx_choice, MIC_val)
    infoBox(
      title = "Cmax",
      value = paste0(metrics$Cmax, " mg/L"),
      subtitle = paste0("Ctrough: ", metrics$Ctrough, " mg/L"),
      icon  = icon("arrow-up"),
      color = "blue",
      fill  = TRUE
    )
  })

  # ---------------------------------------------------------------------------
  # TAB 3 OUTPUTS (CYTOKINES)
  # ---------------------------------------------------------------------------

  output$cyto_pro_plot <- renderPlot({
    sim <- sim_base()
    req(all(c("TNF", "IL6", "IL1") %in% names(sim)))

    cyto_long <- sim %>%
      dplyr::select(time, `TNF-α` = TNF, `IL-6` = IL6, `IL-1β` = IL1) %>%
      pivot_longer(-time, names_to = "Cytokine", values_to = "Level")

    ggplot(cyto_long, aes(time, Level, color = Cytokine)) +
      geom_line(size = 0.9) +
      geom_area(aes(fill = Cytokine), alpha = 0.08, position = "identity") +
      scale_color_manual(values = PAL[1:3]) +
      scale_fill_manual(values  = PAL[1:3]) +
      labs(title = "Pro-Inflammatory Cytokines",
           subtitle = "Relative Units",
           x = "Time (h)", y = "Concentration (RU)", color = NULL, fill = NULL) +
      theme_sep()
  }, bg = "#1a1f2e")

  output$cyto_anti_plot <- renderPlot({
    sim <- sim_base()
    req(all(c("IL10", "HMGB1") %in% names(sim)))

    cyto_long <- sim %>%
      dplyr::select(time, `IL-10` = IL10, `HMGB1 (×10)` = HMGB1) %>%
      mutate(`HMGB1 (×10)` = `HMGB1 (×10)` * 10) %>%
      pivot_longer(-time, names_to = "Mediator", values_to = "Level")

    ggplot(cyto_long, aes(time, Level, color = Mediator)) +
      geom_line(size = 0.9) +
      scale_color_manual(values = c("#52c4a0", "#f0a354")) +
      labs(title = "Anti-Inflammatory & Late Mediators",
           x = "Time (h)", y = "Concentration (RU)", color = NULL) +
      theme_sep()
  }, bg = "#1a1f2e")

  output$cyto_apr_plot <- renderPlot({
    sim <- sim_base()
    req(all(c("PCT", "CRP") %in% names(sim)))

    apr_long <- sim %>%
      dplyr::select(time, `PCT (ng/mL)` = PCT, `CRP/10 (mg/L)` = CRP) %>%
      mutate(`CRP/10 (mg/L)` = `CRP/10 (mg/L)` / 10) %>%
      pivot_longer(-time, names_to = "Marker", values_to = "Level")

    ggplot(apr_long, aes(time, Level, color = Marker)) +
      geom_line(size = 0.9) +
      geom_hline(yintercept = 0.5, linetype = "dashed",
                 color = "#f0a354", alpha = 0.7, size = 0.6) +
      annotate("text", x = 2, y = 0.65, label = "PCT norm ≤0.5",
               color = "#f0a354", size = 3) +
      scale_color_manual(values = c("#e05c5c", "#b47cf5")) +
      labs(title = "Acute Phase Reactants",
           x = "Time (h)", y = "Concentration", color = NULL) +
      theme_sep()
  }, bg = "#1a1f2e")

  # ---------------------------------------------------------------------------
  # TAB 4 OUTPUTS (HEMODYNAMICS)
  # ---------------------------------------------------------------------------

  output$map_plot <- renderPlot({
    sim <- sim_base()
    req("MAP" %in% names(sim))

    ggplot(sim, aes(time, MAP)) +
      geom_ribbon(aes(ymin = 65, ymax = pmax(MAP, 65)), fill = "#52c4a0", alpha = 0.1) +
      geom_ribbon(aes(ymin = pmin(MAP, 65), ymax = 65), fill = "#e05c5c", alpha = 0.12) +
      geom_line(color = "#4e9af1", size = 1) +
      geom_hline(yintercept = 65, linetype = "dashed", color = "#52c4a0",
                 size = 0.8, alpha = 0.9) +
      annotate("text", x = 1, y = 67, label = "Target MAP 65 mmHg",
               color = "#52c4a0", size = 3.2, hjust = 0) +
      scale_y_continuous(limits = c(40, 110)) +
      labs(title = "Mean Arterial Pressure",
           x = "Time (h)", y = "MAP (mmHg)") +
      theme_sep()
  }, bg = "#1a1f2e")

  output$lactate_plot <- renderPlot({
    sim <- sim_base()
    req("Lac" %in% names(sim))

    # Lactate clearance at 6h
    lac_0   <- sim$Lac[1]
    lac_6   <- sim$Lac[which.min(abs(sim$time - 6))]
    clearance <- round((lac_0 - lac_6) / lac_0 * 100, 1)

    ggplot(sim, aes(time, Lac)) +
      geom_ribbon(aes(ymin = 0, ymax = Lac), fill = "#e05c5c", alpha = 0.12) +
      geom_line(color = "#f0a354", size = 1) +
      geom_hline(yintercept = 2, linetype = "dashed", color = "#f5d44e",
                 size = 0.7, alpha = 0.8) +
      geom_hline(yintercept = 4, linetype = "dashed", color = "#e05c5c",
                 size = 0.7, alpha = 0.8) +
      annotate("text", x = 1, y = 2.2, label = "Normal <2", color = "#f5d44e", size = 3) +
      annotate("text", x = 1, y = 4.2, label = "Shock >4",  color = "#e05c5c", size = 3) +
      annotate("text", x = 40, y = max(sim$Lac) * 0.9,
               label = paste0("6h clearance: ", clearance, "%"),
               color = "#52c4a0", size = 3.5) +
      labs(title = "Lactate Dynamics",
           x = "Time (h)", y = "Lactate (mmol/L)") +
      theme_sep()
  }, bg = "#1a1f2e")

  output$sofa_plot <- renderPlot({
    sim <- sim_base()
    req("SOFA" %in% names(sim) && "SOFA_resp" %in% names(sim))

    sofa_long <- sim[seq(1, nrow(sim), by = 4), ] %>%
      dplyr::select(time,
                    Respiratory = SOFA_resp, Coagulation = SOFA_coag,
                    Hepatic = SOFA_hep, Cardiovascular = SOFA_cv,
                    Renal = SOFA_ren, CNS = SOFA_cns) %>%
      pivot_longer(-time, names_to = "Component", values_to = "Score")

    ggplot(sofa_long, aes(time, Score, color = Component)) +
      geom_line(size = 0.8, alpha = 0.9) +
      scale_color_manual(values = PAL) +
      scale_y_continuous(breaks = 0:4, limits = c(0, 4.5)) +
      labs(title = "SOFA Component Scores",
           subtitle = "0=normal, 4=severe",
           x = "Time (h)", y = "Score", color = NULL) +
      theme_sep()
  }, bg = "#1a1f2e")

  output$organ_plot <- renderPlot({
    sim <- sim_base()
    req(all(c("Bili", "INR") %in% names(sim)))

    # Creatinine derived from SOFA_ren
    sim$scr_est <- input$scr_base * (1 + 1.5 * sim$D / (sim$D + 2))

    org_long <- sim[seq(1, nrow(sim), by = 4), ] %>%
      dplyr::select(time,
                    `Creatinine (mg/dL)` = scr_est,
                    `Bilirubin (mg/dL)`  = Bili,
                    `INR`                = INR) %>%
      pivot_longer(-time, names_to = "Marker", values_to = "Value")

    ggplot(org_long, aes(time, Value, color = Marker)) +
      geom_line(size = 0.9) +
      scale_color_manual(values = c("#f0a354", "#b47cf5", "#52c4a0")) +
      facet_wrap(~Marker, scales = "free_y", ncol = 1) +
      labs(title = "Organ-Specific Biochemistry",
           x = "Time (h)", y = NULL, color = NULL) +
      theme_sep() +
      theme(legend.position = "none")
  }, bg = "#1a1f2e")

  output$dic_plot <- renderPlot({
    sim <- sim_base()
    req(all(c("Plt_var", "INR", "HMGB1") %in% names(sim)))

    dic_long <- sim[seq(1, nrow(sim), by = 4), ] %>%
      dplyr::select(time,
                    `Platelet (×10³/μL)` = Plt_var,
                    `INR (×100)`         = INR,
                    `HMGB1 (×50)`        = HMGB1) %>%
      mutate(`INR (×100)`  = `INR (×100)` * 100,
             `HMGB1 (×50)` = `HMGB1 (×50)` * 50) %>%
      pivot_longer(-time, names_to = "Marker", values_to = "Value")

    ggplot(dic_long, aes(time, Value, color = Marker)) +
      geom_line(size = 0.9) +
      scale_color_manual(values = c("#f5d44e", "#e05c5c", "#b47cf5")) +
      labs(title = "DIC Markers (Scaled)",
           x = "Time (h)", y = "Scaled Value", color = NULL) +
      theme_sep()
  }, bg = "#1a1f2e")

  output$ards_plot <- renderPlot({
    sim <- sim_base()
    req("PF_ratio" %in% names(sim))

    ggplot(sim[seq(1, nrow(sim), by = 4), ], aes(time, PF_ratio)) +
      geom_ribbon(aes(ymin = 100, ymax = pmin(PF_ratio, 200)),
                  fill = "#e05c5c", alpha = 0.2) +
      geom_ribbon(aes(ymin = 200, ymax = pmin(PF_ratio, 300)),
                  fill = "#f0a354", alpha = 0.15) +
      geom_line(color = "#4e9af1", size = 1) +
      geom_hline(yintercept = 300, color = "#f5d44e",
                 linetype = "dashed", size = 0.7) +
      geom_hline(yintercept = 200, color = "#f0a354",
                 linetype = "dashed", size = 0.7) +
      geom_hline(yintercept = 100, color = "#e05c5c",
                 linetype = "dashed", size = 0.7) +
      annotate("text", x = 2, y = 310, label = "Mild ARDS <300", color = "#f5d44e", size = 3) +
      annotate("text", x = 2, y = 210, label = "Moderate <200", color = "#f0a354", size = 3) +
      annotate("text", x = 2, y = 110, label = "Severe <100",   color = "#e05c5c", size = 3) +
      scale_y_continuous(limits = c(80, 420)) +
      labs(title = "PaO₂/FiO₂ Ratio (ARDS Assessment)",
           x = "Time (h)", y = "PaO₂/FiO₂ (mmHg)") +
      theme_sep()
  }, bg = "#1a1f2e")

  # ---------------------------------------------------------------------------
  # TAB 5 OUTPUTS (SCENARIOS)
  # ---------------------------------------------------------------------------

  output$scenario_bact_plot <- renderPlot({
    sc <- sim_scenarios()
    req(!is.null(sc) && "B" %in% names(sc))

    sc_thin <- sc[seq(1, nrow(sc), by = 4), ]
    sc_thin$log_B <- log10(pmax(sc_thin$B, 1))

    ggplot(sc_thin, aes(time, log_B, color = scenario)) +
      geom_line(size = 1) +
      scale_color_manual(values = PAL) +
      geom_hline(yintercept = 3, linetype = "dashed", color = "#7a849e",
                 size = 0.6, alpha = 0.7) +
      annotate("text", x = 2, y = 3.2, label = "Clinical threshold",
               color = "#7a849e", size = 3) +
      labs(title = "Bacterial Load Over Time",
           x = "Time (h)", y = "log₁₀ CFU/mL", color = "Scenario") +
      theme_sep()
  }, bg = "#1a1f2e")

  output$scenario_map_plot <- renderPlot({
    sc <- sim_scenarios()
    req(!is.null(sc) && "MAP" %in% names(sc))
    sc_thin <- sc[seq(1, nrow(sc), by = 4), ]

    ggplot(sc_thin, aes(time, MAP, color = scenario)) +
      geom_line(size = 0.9) +
      geom_hline(yintercept = 65, linetype = "dashed", color = "#52c4a0",
                 size = 0.7, alpha = 0.8) +
      scale_color_manual(values = PAL) +
      scale_y_continuous(limits = c(40, 110)) +
      labs(title = "MAP by Scenario",
           x = "Time (h)", y = "MAP (mmHg)", color = NULL) +
      theme_sep() +
      theme(legend.position = "bottom",
            legend.key.size = unit(0.4, "cm"),
            legend.text = element_text(size = 8))
  }, bg = "#1a1f2e")

  output$scenario_sofa_plot <- renderPlot({
    sc <- sim_scenarios()
    req(!is.null(sc) && "SOFA" %in% names(sc))
    sc_thin <- sc[seq(1, nrow(sc), by = 4), ]

    ggplot(sc_thin, aes(time, SOFA, color = scenario)) +
      geom_line(size = 0.9) +
      scale_color_manual(values = PAL) +
      labs(title = "Total SOFA Score by Scenario",
           x = "Time (h)", y = "SOFA Score", color = NULL) +
      theme_sep() +
      theme(legend.position = "bottom",
            legend.key.size = unit(0.4, "cm"),
            legend.text = element_text(size = 8))
  }, bg = "#1a1f2e")

  output$mortality_table <- renderTable({
    sc <- sim_scenarios()
    req(!is.null(sc) && "SOFA" %in% names(sc))

    scenarios_unique <- unique(sc$scenario)

    mort_df <- purrr::map_dfr(scenarios_unique, function(sname) {
      sub_df <- sc[sc$scenario == sname, ]
      sofa_72 <- sub_df$SOFA[which.min(abs(sub_df$time - 72))]
      lac_72  <- sub_df$Lac[which.min(abs(sub_df$time - 72))]
      map_72  <- sub_df$MAP[which.min(abs(sub_df$time - 72))]
      b_72    <- sub_df$B[which.min(abs(sub_df$time - 72))]

      mort <- calc_mortality(
        age       = input$age,
        SOFA_score= max(0, sofa_72),
        Lac       = max(0.5, lac_72),
        MAP       = map_72,
        severity  = input$severity,
        comorbidity = input$comorbidity
      )

      data.frame(
        Scenario          = sname,
        `SOFA (72h)`      = round(sofa_72, 1),
        `MAP (72h, mmHg)` = round(map_72, 1),
        `Lactate (72h)`   = round(lac_72, 2),
        `log₁₀B (72h)`    = round(log10(max(b_72, 1)), 2),
        `28d Mortality (%)` = mort,
        check.names = FALSE
      )
    })
    mort_df[order(mort_df$`28d Mortality (%)`), ]
  }, striped = FALSE, width = "100%", digits = 1)

  # ---------------------------------------------------------------------------
  # TAB 6 OUTPUTS (BIOMARKERS)
  # ---------------------------------------------------------------------------

  output$lactate_box <- renderValueBox({
    sim <- sim_base()
    lac <- round(sim$Lac[1], 1)
    color <- ifelse(lac > 4, "red", ifelse(lac > 2, "yellow", "green"))
    valueBox(paste0(lac, " mmol/L"), "Lactate (0h)",
             icon = icon("tint"), color = color)
  })

  output$pct_box <- renderValueBox({
    sim <- sim_base()
    pct <- round(sim$PCT[1], 2)
    color <- ifelse(pct > 10, "red", ifelse(pct > 2, "yellow", "green"))
    valueBox(paste0(pct, " ng/mL"), "Procalcitonin",
             icon = icon("vial"), color = color)
  })

  output$crp_box <- renderValueBox({
    sim <- sim_base()
    crp <- round(sim$CRP[1], 1)
    color <- ifelse(crp > 100, "red", ifelse(crp > 50, "yellow", "green"))
    valueBox(paste0(crp, " mg/L"), "CRP",
             icon = icon("fire"), color = color)
  })

  output$sofa_score_box <- renderValueBox({
    sim <- sim_base()
    sofa <- round(sim$SOFA[1], 0)
    color <- ifelse(sofa > 10, "red", ifelse(sofa > 6, "yellow", "green"))
    valueBox(sofa, "SOFA Score (0h)",
             icon = icon("stethoscope"), color = color)
  })

  output$biomarker_panel <- renderPlot({
    sim <- sim_base()
    req(all(c("Lac", "PCT", "CRP", "WBC", "Plt_var", "INR", "Bili") %in% names(sim)))

    sim_thin <- sim[seq(1, nrow(sim), by = 4), ]

    # Normalize to common scale (0-100%) for dashboard view
    norm_col <- function(x) {
      rng <- range(x, na.rm = TRUE)
      if (diff(rng) == 0) return(rep(0.5, length(x)))
      (x - rng[1]) / diff(rng) * 100
    }

    bio_long <- sim_thin %>%
      dplyr::select(time, Lactate = Lac, PCT, CRP, WBC, Platelet = Plt_var,
                    INR, Bilirubin = Bili) %>%
      pivot_longer(-time, names_to = "Biomarker", values_to = "Value")

    ggplot(bio_long, aes(time, Value, color = Biomarker, fill = Biomarker)) +
      geom_line(size = 0.8) +
      geom_area(alpha = 0.06, position = "identity") +
      facet_wrap(~Biomarker, scales = "free_y", ncol = 2) +
      scale_color_manual(values = PAL) +
      scale_fill_manual(values  = PAL) +
      labs(title = "Biomarker Trajectories (72h)",
           x = "Time (h)", y = NULL) +
      theme_sep() +
      theme(legend.position = "none",
            strip.text = element_text(size = 9))
  }, bg = "#1a1f2e")

  output$sofa_timeline_table <- renderTable({
    sim <- sim_base()
    req("SOFA" %in% names(sim))

    time_points <- c(0, 6, 12, 24, 48, 72)
    rows <- purrr::map_dfr(time_points, function(tp) {
      idx <- which.min(abs(sim$time - tp))
      data.frame(
        `Time (h)` = tp,
        SOFA    = round(sim$SOFA[idx], 1),
        MAP     = round(sim$MAP[idx], 1),
        Lactate = round(sim$Lac[idx], 2),
        PCT     = round(sim$PCT[idx], 2),
        check.names = FALSE
      )
    })
    rows
  }, striped = FALSE, digits = 1, width = "100%")

  output$deescalation_text <- renderText({
    sim  <- sim_base()
    req("PCT" %in% names(sim) && "B" %in% names(sim))

    pct_24 <- sim$PCT[which.min(abs(sim$time - 24))]
    pct_0  <- sim$PCT[1]
    pct_decline <- (pct_0 - pct_24) / pct_0 * 100

    if (pct_decline > 80 && pct_24 < 1) {
      "PCT declining >80% and <1 ng/mL: Consider antibiotic de-escalation or stopping at 5-7 days per Procalcitonin Stop trial (de Jong 2016)."
    } else if (pct_decline > 50) {
      "PCT declining >50%: Treatment appears effective. Re-evaluate at 72h for potential de-escalation based on clinical response and culture results."
    } else {
      "PCT not declining adequately: Reassess antibiotic spectrum, consider source control, evaluate for resistant organisms. Repeat cultures recommended."
    }
  })

  output$survival_plot <- renderPlot({
    sim <- sim_base()
    req("SOFA" %in% names(sim))

    # Survival probability over 28 days based on evolving SOFA
    days <- seq(0, 28, by = 1)
    times_72 <- pmin(sim$time, 72)

    # Use SOFA trajectory to estimate time-varying mortality hazard
    sofa_t <- approxfun(sim$time, sim$SOFA, rule = 2)
    lac_t  <- approxfun(sim$time, sim$Lac,  rule = 2)

    surv <- numeric(length(days))
    surv[1] <- 1
    for (i in 2:length(days)) {
      day_h  <- (days[i] - 1) * 24
      day_h  <- pmin(day_h, 72)
      sofa_i <- max(0, sofa_t(day_h))
      lac_i  <- max(0.5, lac_t(day_h))
      # Daily hazard (simplified Weibull-like)
      h_daily <- (0.004 + 0.008 * sofa_i + 0.005 * max(0, lac_i - 2)) *
                 ifelse(input$comorbidity == "immunosuppressed", 1.5,
                 ifelse(input$comorbidity == "CKD", 1.2, 1.0)) *
                 (1 + (input$age - 60) / 100)
      surv[i] <- surv[i-1] * exp(-h_daily)
    }

    surv_df <- data.frame(Day = days, Survival = pmax(0.01, surv) * 100)

    ggplot(surv_df, aes(Day, Survival)) +
      geom_ribbon(aes(ymin = pmax(0, Survival - 8), ymax = pmin(100, Survival + 8)),
                  fill = "#4e9af1", alpha = 0.12) +
      geom_line(color = "#4e9af1", size = 1.1) +
      geom_hline(yintercept = 50, linetype = "dashed", color = "#7a849e",
                 size = 0.6, alpha = 0.7) +
      scale_y_continuous(limits = c(0, 100), labels = function(x) paste0(x, "%")) +
      scale_x_continuous(breaks = c(0, 7, 14, 21, 28)) +
      labs(title = "Estimated Survival Probability",
           subtitle = "Based on SOFA trajectory and patient risk factors",
           x = "Day", y = "Survival (%)") +
      theme_sep()
  }, bg = "#1a1f2e")

  output$qsofa_plot <- renderPlot({
    sim <- sim_base()
    req(all(c("SOFA", "MAP") %in% names(sim)))

    # qSOFA components: RR > 22, AMS (simplified), SBP < 100
    # MAP ~= (SBP + 2*DBP)/3; approximate SBP = MAP * 1.4
    sim_thin <- sim[seq(1, nrow(sim), by = 4), ]
    sim_thin$SBP_approx <- sim_thin$MAP * 1.4
    sim_thin$qSOFA_cv  <- as.integer(sim_thin$SBP_approx < 100)
    sim_thin$qSOFA_neuro <- as.integer(sim_thin$D > 2)  # damage proxy for AMS
    sim_thin$qSOFA_resp  <- as.integer(sim_thin$SOFA_resp >= 2)
    sim_thin$qSOFA_total <- sim_thin$qSOFA_cv + sim_thin$qSOFA_neuro + sim_thin$qSOFA_resp

    qsofa_long <- sim_thin %>%
      dplyr::select(time,
                    `Total qSOFA` = qSOFA_total,
                    `SOFA/3`      = SOFA) %>%
      mutate(`SOFA/3` = `SOFA/3` / 3) %>%
      pivot_longer(-time, names_to = "Score", values_to = "Value")

    ggplot(qsofa_long, aes(time, Value, color = Score)) +
      geom_line(size = 0.9) +
      geom_hline(yintercept = 2, linetype = "dashed", color = "#e05c5c",
                 size = 0.7, alpha = 0.7) +
      annotate("text", x = 2, y = 2.15, label = "qSOFA ≥2 = high risk",
               color = "#e05c5c", size = 3) +
      scale_color_manual(values = c("#4e9af1", "#52c4a0")) +
      scale_y_continuous(breaks = 0:7) +
      labs(title = "qSOFA & SOFA Trends",
           x = "Time (h)", y = "Score", color = NULL) +
      theme_sep()
  }, bg = "#1a1f2e")

}

# =============================================================================
# LAUNCH
# =============================================================================

shinyApp(ui = ui, server = server)
