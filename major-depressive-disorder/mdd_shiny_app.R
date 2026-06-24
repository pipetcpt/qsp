# =============================================================================
# Major Depressive Disorder (MDD) — Comprehensive QSP Shiny Dashboard
# =============================================================================
# Description : Interactive 6-tab Shiny app for MDD QSP simulation
# Tabs        : 1. Patient Profile
#               2. Drug PK (plasma PK, SERT/NET occupancy)
#               3. Neurotransmitter Dynamics (5-HT, NE, DA)
#               4. HPA Axis & Neuroinflammation (cortisol, IL-6, BDNF)
#               5. Clinical Endpoints (HDRS, response/remission)
#               6. Scenario Comparison (forest plot, radar chart)
# Date        : 2026-06-20
# =============================================================================

library(shiny)
library(shinydashboard)
library(ggplot2)
library(dplyr)
library(tidyr)
library(plotly)
library(DT)

# =============================================================================
# INLINE QSP ODE SIMULATION (Pure R — no mrgsolve dependency for Shiny)
# =============================================================================

run_mdd_sim <- function(
    # Drug doses
    dose_esc   = 0,    # Escitalopram mg/day
    dose_ven   = 0,    # Venlafaxine mg/day
    dose_ket   = 0,    # Ketamine mg/kg (single IV bolus at t=0)
    dose_bup   = 0,    # Bupropion mg/day
    # Disease severity
    mdd_5ht_def  = 0.60,  # 5-HT deficit fraction
    mdd_ne_def   = 0.70,
    mdd_da_def   = 0.65,
    mdd_bdnf_def = 0.65,
    mdd_cort_ex  = 1.40,
    mdd_il6_ex   = 1.80,
    hdrs_base    = 22.0,
    stress_level = 0.01,  # ongoing stress input
    # Augmentation
    aripiprazole = FALSE,
    # Simulation duration
    weeks_sim = 8,
    bwt = 70
) {
  # --- Parameters ---
  SS_5HT <- 1.5; SS_NE <- 1.0; SS_DA <- 0.8
  SS_BDNF <- 2.0; SS_CORT <- 15.0; SS_IL6 <- 0.625
  EMAX_SERT <- 0.95
  KI_SERT_ESC <- 1.1   # nM
  KI_SERT_VEN <- 7.5   # nM
  KI_NET_VEN  <- 2.7   # nM
  KI_NET_BUP  <- 52    # nM (bupropion NET)
  KI_DAT_BUP  <- 526   # nM (bupropion DAT)
  KI_NMDA_KET <- 3.0   # uM

  # PK constants
  KA_ESC <- 0.46; CL_ESC <- 37; V_ESC <- 1090; F_ESC <- 0.80; MW_ESC <- 324.39
  KA_VEN <- 1.2;  CL_VEN <- 90; V_VEN  <- 500;  F_VEN  <- 0.45; MW_VEN <- 277.40
  KE_KET <- 1.73; V_KET  <- 3.0 * bwt; MW_KET <- 237.73
  KA_BUP <- 0.35; CL_BUP <- 200; V_BUP <- 3000; F_BUP <- 0.85; MW_BUP <- 239.74

  # Neurotransmitter
  KSYN_5HT <- 0.12; KDEG_5HT <- 0.08; KREUP_5HT <- 0.15; K_AUTO5HT <- 0.3
  KSYN_NE  <- 0.08; KDEG_NE  <- 0.06; KREUP_NE  <- 0.12
  KSYN_DA  <- 0.06; KDEG_DA  <- 0.04; KREUP_DA  <- 0.10

  if (aripiprazole) {
    KSYN_DA  <- KSYN_DA  * 1.15
    KSYN_NE  <- KSYN_NE  * 1.15
    KSYN_5HT <- KSYN_5HT * 1.10
  }

  # HPA
  KCRH_SYNTH <- 0.10; KCRH_DEG <- 0.15
  KACTH_SYN  <- 0.20; KACTH_DEG <- 0.12
  KCORT_SYN  <- 0.30; KCORT_DEG <- 0.05; CORT_FB <- 0.04

  # IL-6
  KIL6_SYNTH <- 0.05; KIL6_DEG <- 0.08; KIL6_CORT <- 0.03

  # BDNF / mTOR / Neurogenesis
  KBDNF_SYNTH <- 0.04; KBDNF_DEG <- 0.02
  BDNF_5HT_K  <- 0.30; BDNF_CORT_K <- 0.15; BDNF_IL6_K <- 0.10; BDNF_KET_K <- 2.0; BDNF_NE_K <- 0.15
  KMTOR_ACT <- 0.05; KMTOR_DEG <- 0.03
  KNEURO_SYN <- 0.01; KNEURO_DEG <- 0.005

  # HDRS mapping
  HDRS_MAX_5HT <- 6.0; HDRS_MAX_NE <- 4.0; HDRS_MAX_BDNF <- 5.0
  HDRS_MAX_NEURO <- 3.0; HDRS_MAX_CORT <- 2.0; EC50_HDRS_5HT <- 0.5

  # --- Initial conditions (MDD state) ---
  dt    <- 2  # hours
  tmax  <- weeks_sim * 168
  times <- seq(0, tmax, by = dt)
  n     <- length(times)

  # Initial state vector
  state <- list(
    DEPOT_ESC  = 0, CENTRAL_ESC = 0,
    DEPOT_VEN  = 0, CENTRAL_VEN = 0,
    CENTRAL_KET = (dose_ket * bwt) / V_KET,  # IV bolus C0
    DEPOT_BUP  = 0, CENTRAL_BUP = 0,
    X5HT = SS_5HT * mdd_5ht_def,
    NE   = SS_NE  * mdd_ne_def,
    DA   = SS_DA  * mdd_da_def,
    CRH  = 0.67  * mdd_cort_ex,
    ACTH = 1.67  * mdd_cort_ex,
    CORT = SS_CORT * mdd_cort_ex,
    IL6  = SS_IL6  * mdd_il6_ex,
    BDNF = SS_BDNF * mdd_bdnf_def,
    MTOR = 0.5,
    NEURO = 0.75,
    HDRS = hdrs_base
  )

  # Output matrices
  out <- data.frame(
    time = times,
    ESC_nM = NA, VEN_nM = NA, KET_uM = NA, BUP_nM = NA,
    SERT_PCT = NA, NET_PCT = NA,
    X5HT = NA, NE = NA, DA = NA,
    CRH = NA, ACTH = NA, CORT = NA,
    IL6 = NA, BDNF = NA, MTOR = NA, NEURO = NA,
    HDRS = NA, PHQ9 = NA, MADRS = NA,
    RESPONSE = NA, REMISSION = NA
  )

  dose_rate_esc <- dose_esc / 24
  dose_rate_ven <- dose_ven / 24
  dose_rate_bup <- dose_bup / 24

  s <- state
  for (i in seq_along(times)) {
    # --- PK: Escitalopram ---
    dDEPOT_ESC  <- F_ESC * dose_rate_esc - KA_ESC * s$DEPOT_ESC
    dCENT_ESC   <- KA_ESC * s$DEPOT_ESC - (CL_ESC / V_ESC) * s$CENTRAL_ESC

    # --- PK: Venlafaxine ---
    dDEPOT_VEN  <- F_VEN * dose_rate_ven - KA_VEN * s$DEPOT_VEN
    dCENT_VEN   <- KA_VEN * s$DEPOT_VEN - (CL_VEN / V_VEN) * s$CENTRAL_VEN

    # --- PK: Ketamine ---
    dCENT_KET   <- -KE_KET * s$CENTRAL_KET

    # --- PK: Bupropion ---
    dDEPOT_BUP  <- F_BUP * dose_rate_bup - KA_BUP * s$DEPOT_BUP
    dCENT_BUP   <- KA_BUP * s$DEPOT_BUP - (CL_BUP / V_BUP) * s$CENTRAL_BUP

    # Concentrations in nM / uM
    Cp_ESC_nM <- s$CENTRAL_ESC / MW_ESC * 1e6
    Cp_VEN_nM <- s$CENTRAL_VEN / MW_VEN * 1e6
    Cp_KET_uM <- s$CENTRAL_KET / MW_KET * 1e3
    Cp_BUP_nM <- s$CENTRAL_BUP / MW_BUP * 1e6

    # SERT occupancy
    SERT_ESC <- EMAX_SERT * Cp_ESC_nM^1.5 / (KI_SERT_ESC^1.5 + Cp_ESC_nM^1.5)
    SERT_VEN <- EMAX_SERT * Cp_VEN_nM / (KI_SERT_VEN + Cp_VEN_nM)
    SERT_tot <- 1 - (1 - SERT_ESC) * (1 - SERT_VEN)

    # NET occupancy
    NET_VEN <- EMAX_SERT * Cp_VEN_nM / (KI_NET_VEN + Cp_VEN_nM)
    NET_BUP <- EMAX_SERT * Cp_BUP_nM / (KI_NET_BUP + Cp_BUP_nM)
    NET_tot <- 1 - (1 - NET_VEN) * (1 - NET_BUP)

    # DAT (bupropion)
    DAT_BUP <- EMAX_SERT * Cp_BUP_nM / (KI_DAT_BUP + Cp_BUP_nM)

    # NMDA block (ketamine)
    NMDA_bl <- Cp_KET_uM / (KI_NMDA_KET + Cp_KET_uM)

    # --- 5-HT dynamics ---
    reup_5HT <- KREUP_5HT * (1 - SERT_tot) * s$X5HT
    fb_5HT   <- K_AUTO5HT * s$X5HT / (s$X5HT + SS_5HT)
    synth_5HT <- KSYN_5HT * (1 - fb_5HT * 0.5)
    dX5HT <- synth_5HT - KDEG_5HT * s$X5HT - reup_5HT

    # --- NE dynamics ---
    reup_NE <- KREUP_NE * (1 - NET_tot) * s$NE
    dNE <- KSYN_NE - KDEG_NE * s$NE - reup_NE

    # --- DA dynamics ---
    reup_DA <- KREUP_DA * (1 - DAT_BUP) * s$DA
    dDA <- KSYN_DA - KDEG_DA * s$DA - reup_DA

    # --- HPA axis ---
    fb_cort <- 1 / (1 + CORT_FB * s$CORT)
    dCRH  <- (KCRH_SYNTH + stress_level) * fb_cort - KCRH_DEG * s$CRH
    dACTH <- KACTH_SYN * s$CRH - KACTH_DEG * s$ACTH
    dCORT <- KCORT_SYN * s$ACTH - KCORT_DEG * s$CORT

    # --- IL-6 ---
    il6_stress <- KIL6_CORT * (s$CORT - SS_CORT) / SS_CORT
    il6_ssri   <- 0.1 * SERT_tot
    dIL6 <- KIL6_SYNTH * (1 + il6_stress) - KIL6_DEG * s$IL6 - il6_ssri * s$IL6

    # --- BDNF ---
    stim_BDNF <- BDNF_5HT_K  * (s$X5HT - SS_5HT * mdd_5ht_def)  / SS_5HT +
                 BDNF_NE_K   * (s$NE   - SS_NE  * mdd_ne_def)    / SS_NE  +
                 BDNF_KET_K  * NMDA_bl -
                 BDNF_CORT_K * (s$CORT - SS_CORT) / SS_CORT -
                 BDNF_IL6_K  * (s$IL6  - SS_IL6) / SS_IL6
    dBDNF <- KBDNF_SYNTH * (1 + stim_BDNF) - KBDNF_DEG * s$BDNF
    dBDNF <- max(dBDNF, -0.05 * s$BDNF)  # stability guard

    # --- mTOR ---
    dMTOR <- KMTOR_ACT * (s$BDNF / SS_BDNF) + 0.5 * NMDA_bl - KMTOR_DEG * s$MTOR

    # --- Neurogenesis ---
    neuro_cort_sup <- 0.02 * max(0, (s$CORT - SS_CORT) / SS_CORT)
    dNEURO <- KNEURO_SYN * (s$BDNF / SS_BDNF) * s$MTOR - KNEURO_DEG * s$NEURO - neuro_cort_sup * s$NEURO

    # --- HDRS ---
    d5HT_d   <- s$X5HT - SS_5HT * mdd_5ht_def
    dNE_d    <- s$NE   - SS_NE  * mdd_ne_def
    dBDNF_d  <- s$BDNF - SS_BDNF * mdd_bdnf_def
    dNEURO_d <- s$NEURO - 0.75
    dCORT_d  <- s$CORT - SS_CORT * mdd_cort_ex

    hdrs_5ht   <- HDRS_MAX_5HT   * d5HT_d   / (EC50_HDRS_5HT + abs(d5HT_d))
    hdrs_ne    <- HDRS_MAX_NE    * dNE_d    / (0.3 + abs(dNE_d))
    hdrs_bdnf  <- HDRS_MAX_BDNF  * dBDNF_d  / (0.5 + abs(dBDNF_d))
    hdrs_neuro <- HDRS_MAX_NEURO * dNEURO_d / (0.2 + abs(dNEURO_d))
    hdrs_cort  <- HDRS_MAX_CORT  * (-dCORT_d) / (5.0 + abs(dCORT_d))

    hdrs_tgt <- hdrs_base - (hdrs_5ht + hdrs_ne + hdrs_bdnf + hdrs_neuro + hdrs_cort)
    hdrs_tgt <- max(0, min(52, hdrs_tgt))
    dHDRS <- 0.005 * (hdrs_tgt - s$HDRS)

    # --- Euler update ---
    s$DEPOT_ESC  <- max(0, s$DEPOT_ESC  + dt * dDEPOT_ESC)
    s$CENTRAL_ESC<- max(0, s$CENTRAL_ESC+ dt * dCENT_ESC)
    s$DEPOT_VEN  <- max(0, s$DEPOT_VEN  + dt * dDEPOT_VEN)
    s$CENTRAL_VEN<- max(0, s$CENTRAL_VEN+ dt * dCENT_VEN)
    s$CENTRAL_KET<- max(0, s$CENTRAL_KET+ dt * dCENT_KET)
    s$DEPOT_BUP  <- max(0, s$DEPOT_BUP  + dt * dDEPOT_BUP)
    s$CENTRAL_BUP<- max(0, s$CENTRAL_BUP+ dt * dCENT_BUP)
    s$X5HT       <- max(0, s$X5HT + dt * dX5HT)
    s$NE         <- max(0, s$NE   + dt * dNE)
    s$DA         <- max(0, s$DA   + dt * dDA)
    s$CRH        <- max(0, s$CRH  + dt * dCRH)
    s$ACTH       <- max(0, s$ACTH + dt * dACTH)
    s$CORT       <- max(0, s$CORT + dt * dCORT)
    s$IL6        <- max(0, s$IL6  + dt * dIL6)
    s$BDNF       <- max(0, s$BDNF + dt * dBDNF)
    s$MTOR       <- max(0, s$MTOR + dt * dMTOR)
    s$NEURO      <- max(0, s$NEURO + dt * dNEURO)
    s$HDRS       <- max(0, s$HDRS + dt * dHDRS)

    # Store outputs
    out$ESC_nM[i]   <- Cp_ESC_nM
    out$VEN_nM[i]   <- Cp_VEN_nM
    out$KET_uM[i]   <- Cp_KET_uM
    out$BUP_nM[i]   <- Cp_BUP_nM
    out$SERT_PCT[i] <- SERT_tot * 100
    out$NET_PCT[i]  <- NET_tot  * 100
    out$X5HT[i]     <- s$X5HT
    out$NE[i]       <- s$NE
    out$DA[i]       <- s$DA
    out$CRH[i]      <- s$CRH
    out$ACTH[i]     <- s$ACTH
    out$CORT[i]     <- s$CORT
    out$IL6[i]      <- s$IL6
    out$BDNF[i]     <- s$BDNF
    out$MTOR[i]     <- s$MTOR
    out$NEURO[i]    <- s$NEURO
    out$HDRS[i]     <- s$HDRS
    out$PHQ9[i]     <- s$HDRS * 0.55
    out$MADRS[i]    <- s$HDRS * 1.70
    out$RESPONSE[i] <- as.integer(s$HDRS <= hdrs_base * 0.5)
    out$REMISSION[i]<- as.integer(s$HDRS * 0.55 < 5)
  }
  out$week <- out$time / 168
  out
}

# =============================================================================
# PREDEFINED SCENARIOS
# =============================================================================

run_all_scenarios <- function() {
  list(
    "No Treatment" = run_mdd_sim(dose_esc = 0, dose_ven = 0, dose_ket = 0,
                                  stress_level = 0.02),
    "Escitalopram 10mg" = run_mdd_sim(dose_esc = 10, stress_level = 0.01),
    "Escitalopram 20mg" = run_mdd_sim(dose_esc = 20, stress_level = 0.01),
    "Venlafaxine 150mg" = run_mdd_sim(dose_ven = 150, stress_level = 0.01),
    "Ketamine IV 0.5mg/kg" = run_mdd_sim(dose_ket = 0.5, stress_level = 0.01),
    "ESC + Aripiprazole" = run_mdd_sim(dose_esc = 10, aripiprazole = TRUE,
                                        stress_level = 0.01),
    "Bupropion 300mg" = run_mdd_sim(dose_bup = 300, stress_level = 0.01),
    "TRD (High Stress)" = run_mdd_sim(dose_esc = 20, stress_level = 0.05,
                                       mdd_5ht_def = 0.50, mdd_bdnf_def = 0.50,
                                       mdd_cort_ex = 1.60, mdd_il6_ex = 2.50,
                                       hdrs_base = 26)
  )
}

# Colour palette
SCEN_COLORS <- c(
  "No Treatment"         = "#E74C3C",
  "Escitalopram 10mg"    = "#3498DB",
  "Escitalopram 20mg"    = "#1A5276",
  "Venlafaxine 150mg"    = "#2ECC71",
  "Ketamine IV 0.5mg/kg" = "#9B59B6",
  "ESC + Aripiprazole"   = "#F39C12",
  "Bupropion 300mg"      = "#16A085",
  "TRD (High Stress)"    = "#C0392B"
)

# =============================================================================
# UI
# =============================================================================

ui <- dashboardPage(
  skin = "blue",
  dashboardHeader(
    title = span(icon("brain"), "MDD QSP Dashboard"),
    titleWidth = 300
  ),

  dashboardSidebar(
    width = 300,
    sidebarMenu(
      id = "tabs",
      menuItem("Patient Profile",       tabName = "tab_patient",   icon = icon("user-md")),
      menuItem("Drug PK",               tabName = "tab_pk",        icon = icon("pills")),
      menuItem("Neurotransmitter Dyn.", tabName = "tab_neuro",     icon = icon("bolt")),
      menuItem("HPA & Neuroinflam.",    tabName = "tab_hpa",       icon = icon("fire")),
      menuItem("Clinical Endpoints",    tabName = "tab_clinical",  icon = icon("chart-line")),
      menuItem("Scenario Comparison",   tabName = "tab_compare",   icon = icon("exchange-alt"))
    ),
    hr(),
    h4("Treatment Settings", style = "padding-left:10px; color:#ECF0F1;"),
    sliderInput("dose_esc",  "Escitalopram (mg/day)", 0, 40, 10, step = 5),
    sliderInput("dose_ven",  "Venlafaxine (mg/day)",  0, 225, 0, step = 37.5),
    sliderInput("dose_bup",  "Bupropion (mg/day)",    0, 450, 0, step = 75),
    sliderInput("dose_ket",  "Ketamine IV (mg/kg)",   0, 1, 0, step = 0.1),
    checkboxInput("aripiprazole", "Aripiprazole augmentation", FALSE),
    hr(),
    h4("Disease Severity", style = "padding-left:10px; color:#ECF0F1;"),
    sliderInput("hdrs_base",  "Baseline HDRS-17", 14, 35, 22, step = 1),
    sliderInput("stress_level","Stress Level",    0, 0.1, 0.01, step = 0.005),
    sliderInput("mdd_5ht_def","5-HT Deficit (%)", 30, 90, 60, step = 5,
                post = "%"),
    sliderInput("mdd_bdnf_def","BDNF Deficit (%)",30, 90, 65, step = 5,
                post = "%"),
    hr(),
    sliderInput("weeks_sim", "Simulation Duration (weeks)", 4, 24, 8, step = 2),
    actionButton("run_sim", "Run Simulation", icon = icon("play"),
                 class = "btn-success btn-block"),
    hr(),
    p("MDD QSP Model v1.0", style = "color:#BDC3C7; font-size:11px; padding:5px;"),
    p("2026-06-20", style = "color:#95A5A6; font-size:10px; padding:5px;")
  ),

  dashboardBody(
    tags$head(
      tags$style(HTML("
        .content-wrapper, .right-side { background-color: #f8f9fa; }
        .box { border-radius: 8px; }
        .info-box { border-radius: 8px; }
      "))
    ),
    tabItems(

      # ==============================================================
      # TAB 1: Patient Profile
      # ==============================================================
      tabItem(tabName = "tab_patient",
        fluidRow(
          box(title = "Patient Demographics & MDD Profile",
              status = "primary", solidHeader = TRUE, width = 12,
            fluidRow(
              column(4,
                h4("Patient Parameters"),
                numericInput("age",      "Age (years)",              35, 18, 85, 1),
                selectInput("sex",       "Sex",                      c("Female", "Male")),
                numericInput("bmi",      "BMI (kg/m2)",              24, 15, 50, 0.5),
                selectInput("severity",  "MDD Severity",
                  c("Mild (PHQ-9: 5-9)", "Moderate (PHQ-9: 10-14)",
                    "Moderate-Severe (PHQ-9: 15-19)", "Severe (PHQ-9: ≥20)")),
                numericInput("duration_ep", "Episode Duration (months)", 6, 1, 120, 1)
              ),
              column(4,
                h4("Symptom Assessment"),
                numericInput("phq9_current", "Current PHQ-9 Score",  18,  0, 27, 1),
                numericInput("hdrs_current", "Current HDRS-17 Score",22,  0, 52, 1),
                numericInput("madrs_current","Current MADRS Score",   35,  0, 60, 1),
                checkboxGroupInput("symptoms", "Core Symptoms:",
                  choices = c("Depressed mood", "Anhedonia",
                    "Insomnia/Hypersomnia", "Fatigue", "Concentration impairment",
                    "Psychomotor changes", "Suicidal ideation",
                    "Weight/appetite change", "Worthlessness/guilt"),
                  selected = c("Depressed mood","Anhedonia","Insomnia/Hypersomnia",
                    "Fatigue","Concentration impairment")
                )
              ),
              column(4,
                h4("Treatment History & Comorbidities"),
                numericInput("prior_tx", "# Prior Antidepressant Trials", 1, 0, 10, 1),
                checkboxGroupInput("comorbidities", "Comorbidities:",
                  choices = c("Anxiety disorder","PTSD","Chronic pain",
                    "Cardiovascular disease","Type 2 diabetes",
                    "Hypothyroidism","Substance use disorder",
                    "Inflammatory condition"),
                  selected = c("Anxiety disorder")
                ),
                selectInput("metabolizer","CYP2C19/2D6 Phenotype",
                  c("Normal metabolizer", "Poor metabolizer (2D6)",
                    "Rapid metabolizer (2C19)", "Ultra-rapid metabolizer")),
                checkboxInput("trd", "Treatment-Resistant Depression (TRD)", FALSE)
              )
            )
          )
        ),
        fluidRow(
          valueBoxOutput("vb_phq9",    width = 3),
          valueBoxOutput("vb_hdrs",    width = 3),
          valueBoxOutput("vb_severity",width = 3),
          valueBoxOutput("vb_trd",     width = 3)
        ),
        fluidRow(
          box(title = "Biomarker Profile at Baseline",
              status = "info", solidHeader = TRUE, width = 6,
            plotOutput("plot_biomarker_radar", height = 300)
          ),
          box(title = "MDD Pathophysiology Summary",
              status = "warning", solidHeader = TRUE, width = 6,
            DTOutput("table_patho")
          )
        )
      ),

      # ==============================================================
      # TAB 2: Drug PK
      # ==============================================================
      tabItem(tabName = "tab_pk",
        fluidRow(
          box(title = "Plasma Concentration–Time Profile",
              status = "primary", solidHeader = TRUE, width = 8,
            plotlyOutput("plot_pk_conc", height = 350)
          ),
          box(title = "Steady-State PK Summary",
              status = "info", solidHeader = TRUE, width = 4,
            DTOutput("table_pk_ss")
          )
        ),
        fluidRow(
          box(title = "SERT Occupancy (%) Over Time",
              status = "success", solidHeader = TRUE, width = 6,
            plotlyOutput("plot_sert", height = 300)
          ),
          box(title = "NET Occupancy (%) Over Time",
              status = "warning", solidHeader = TRUE, width = 6,
            plotlyOutput("plot_net", height = 300)
          )
        ),
        fluidRow(
          box(title = "SERT Occupancy–Dose Relationship",
              status = "primary", solidHeader = TRUE, width = 6,
            plotlyOutput("plot_sert_dose", height = 280)
          ),
          box(title = "PK Parameters Table",
              status = "info", solidHeader = TRUE, width = 6,
            DTOutput("table_pk_params")
          )
        )
      ),

      # ==============================================================
      # TAB 3: Neurotransmitter Dynamics
      # ==============================================================
      tabItem(tabName = "tab_neuro",
        fluidRow(
          box(title = "Synaptic Serotonin (5-HT) Dynamics",
              status = "primary", solidHeader = TRUE, width = 6,
            plotlyOutput("plot_5ht", height = 300)
          ),
          box(title = "Synaptic Norepinephrine (NE) Dynamics",
              status = "info", solidHeader = TRUE, width = 6,
            plotlyOutput("plot_ne", height = 300)
          )
        ),
        fluidRow(
          box(title = "Synaptic Dopamine (DA) Dynamics",
              status = "success", solidHeader = TRUE, width = 6,
            plotlyOutput("plot_da", height = 300)
          ),
          box(title = "Monoamine Summary at Week 8",
              status = "warning", solidHeader = TRUE, width = 6,
            plotlyOutput("plot_mono_bar", height = 300)
          )
        ),
        fluidRow(
          box(title = "Neurotransmitter Time-Series Data",
              status = "primary", solidHeader = TRUE, width = 12,
            DTOutput("table_neuro_ts")
          )
        )
      ),

      # ==============================================================
      # TAB 4: HPA Axis & Neuroinflammation
      # ==============================================================
      tabItem(tabName = "tab_hpa",
        fluidRow(
          box(title = "Cortisol Profile (HPA Axis)",
              status = "danger", solidHeader = TRUE, width = 6,
            plotlyOutput("plot_cortisol", height = 300)
          ),
          box(title = "IL-6 (Neuroinflammation)",
              status = "warning", solidHeader = TRUE, width = 6,
            plotlyOutput("plot_il6", height = 300)
          )
        ),
        fluidRow(
          box(title = "BDNF Dynamics",
              status = "success", solidHeader = TRUE, width = 6,
            plotlyOutput("plot_bdnf", height = 300)
          ),
          box(title = "mTOR Activity & Neurogenesis",
              status = "info", solidHeader = TRUE, width = 6,
            plotlyOutput("plot_neuro_index", height = 300)
          )
        ),
        fluidRow(
          box(title = "Kynurenine/Tryptophan Ratio (KYN:TRP proxy)",
              status = "primary", solidHeader = TRUE, width = 6,
            plotlyOutput("plot_kyn", height = 250)
          ),
          box(title = "Biomarker Summary Table",
              status = "info", solidHeader = TRUE, width = 6,
            DTOutput("table_biomarker")
          )
        )
      ),

      # ==============================================================
      # TAB 5: Clinical Endpoints
      # ==============================================================
      tabItem(tabName = "tab_clinical",
        fluidRow(
          valueBoxOutput("vb_hdrs_wk8",    width = 3),
          valueBoxOutput("vb_phq9_wk8",    width = 3),
          valueBoxOutput("vb_response",    width = 3),
          valueBoxOutput("vb_remission",   width = 3)
        ),
        fluidRow(
          box(title = "HDRS-17 Score Over Time",
              status = "primary", solidHeader = TRUE, width = 8,
            plotlyOutput("plot_hdrs", height = 350)
          ),
          box(title = "Response & Remission Timeline",
              status = "success", solidHeader = TRUE, width = 4,
            plotOutput("plot_resp_timeline", height = 350)
          )
        ),
        fluidRow(
          box(title = "PHQ-9 & MADRS Trajectories",
              status = "info", solidHeader = TRUE, width = 6,
            plotlyOutput("plot_phq9_madrs", height = 280)
          ),
          box(title = "Weekly Clinical Summary",
              status = "warning", solidHeader = TRUE, width = 6,
            DTOutput("table_clinical_weekly")
          )
        )
      ),

      # ==============================================================
      # TAB 6: Scenario Comparison
      # ==============================================================
      tabItem(tabName = "tab_compare",
        fluidRow(
          box(title = "HDRS-17 Comparison: All Scenarios",
              status = "primary", solidHeader = TRUE, width = 8,
            plotlyOutput("plot_compare_hdrs", height = 350)
          ),
          box(title = "Select Scenarios to Compare",
              status = "info", solidHeader = TRUE, width = 4,
            checkboxGroupInput("compare_scenarios",
              "Treatment Arms:",
              choices = c("No Treatment","Escitalopram 10mg","Escitalopram 20mg",
                          "Venlafaxine 150mg","Ketamine IV 0.5mg/kg",
                          "ESC + Aripiprazole","Bupropion 300mg","TRD (High Stress)"),
              selected = c("No Treatment","Escitalopram 10mg",
                           "Venlafaxine 150mg","Ketamine IV 0.5mg/kg",
                           "ESC + Aripiprazole")),
            actionButton("run_compare", "Compare All", icon = icon("chart-bar"),
                         class = "btn-primary btn-block")
          )
        ),
        fluidRow(
          box(title = "Effect Size Forest Plot (Week 8 vs Baseline)",
              status = "success", solidHeader = TRUE, width = 6,
            plotOutput("plot_forest", height = 320)
          ),
          box(title = "8-Week Outcome Summary Table",
              status = "warning", solidHeader = TRUE, width = 6,
            DTOutput("table_compare_summary")
          )
        ),
        fluidRow(
          box(title = "BDNF Comparison Across Treatments",
              status = "info", solidHeader = TRUE, width = 6,
            plotlyOutput("plot_compare_bdnf", height = 280)
          ),
          box(title = "Time-to-Response (50% reduction)",
              status = "primary", solidHeader = TRUE, width = 6,
            plotOutput("plot_time_response", height = 280)
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

  # ---- Reactive simulation (user-controlled) ----
  sim_data <- eventReactive(input$run_sim, {
    withProgress(message = "Running QSP simulation...", value = 0.5, {
      run_mdd_sim(
        dose_esc    = input$dose_esc,
        dose_ven    = input$dose_ven,
        dose_bup    = input$dose_bup,
        dose_ket    = input$dose_ket,
        aripiprazole= input$aripiprazole,
        hdrs_base   = input$hdrs_base,
        stress_level= input$stress_level,
        mdd_5ht_def = input$mdd_5ht_def / 100,
        mdd_bdnf_def= input$mdd_bdnf_def / 100,
        weeks_sim   = input$weeks_sim
      )
    })
  }, ignoreNULL = FALSE)

  # ---- All scenarios (for comparison tab) ----
  all_scen_data <- eventReactive(input$run_compare, {
    withProgress(message = "Running all scenarios...", value = 0, {
      all_s <- run_all_scenarios()
      incProgress(1)
      all_s
    })
  }, ignoreNULL = FALSE)

  # Auto-run comparison on load
  observe({
    input$run_compare
  })

  # ==========================================================================
  # TAB 1: Patient Profile
  # ==========================================================================

  output$vb_phq9 <- renderValueBox({
    val <- input$phq9_current
    col <- if (val >= 20) "red" else if (val >= 15) "orange" else if (val >= 10) "yellow" else "green"
    valueBox(val, "PHQ-9 Score", icon = icon("thermometer-three-quarters"), color = col)
  })

  output$vb_hdrs <- renderValueBox({
    val <- input$hdrs_current
    col <- if (val >= 24) "red" else if (val >= 17) "orange" else if (val >= 8) "yellow" else "green"
    valueBox(val, "HDRS-17 Score", icon = icon("chart-bar"), color = col)
  })

  output$vb_severity <- renderValueBox({
    sev <- input$severity
    label <- strsplit(sev, " \\(")[[1]][1]
    col <- switch(label,
      "Mild" = "green", "Moderate" = "yellow",
      "Moderate-Severe" = "orange", "Severe" = "red", "blue")
    valueBox(label, "Severity", icon = icon("exclamation-triangle"), color = col)
  })

  output$vb_trd <- renderValueBox({
    n <- input$prior_tx
    col <- if (n >= 2) "red" else if (n == 1) "orange" else "green"
    label <- if (n >= 2) "TRD Risk" else "Not TRD"
    valueBox(paste0(n, " prior trials"), label, icon = icon("redo"), color = col)
  })

  output$plot_biomarker_radar <- renderPlot({
    # Simulated biomarker values relative to healthy (1.0 = normal)
    df <- data.frame(
      Biomarker = c("5-HT", "NE", "DA", "BDNF", "Cortisol", "IL-6", "NEURO"),
      Patient   = c(input$mdd_5ht_def/100, 0.70, 0.65,
                    input$mdd_bdnf_def/100, input$stress_level*20+1.2,
                    1.5, 0.75),
      Normal    = rep(1.0, 7)
    ) %>% pivot_longer(-Biomarker, names_to = "Group", values_to = "Value")

    ggplot(df, aes(x = Biomarker, y = Value, fill = Group, group = Group)) +
      geom_polygon(alpha = 0.3, color = NA) +
      geom_line(aes(color = Group)) +
      geom_point(aes(color = Group), size = 3) +
      coord_polar() +
      scale_fill_manual(values = c(Patient = "#E74C3C", Normal = "#2ECC71")) +
      scale_color_manual(values = c(Patient = "#E74C3C", Normal = "#2ECC71")) +
      geom_hline(yintercept = 1.0, linetype = "dashed", color = "gray50") +
      labs(title = "Biomarker Profile (Patient vs. Normal)") +
      theme_bw(base_size = 10) +
      theme(legend.position = "right")
  })

  output$table_patho <- renderDT({
    df <- data.frame(
      System       = c("Serotonergic","Noradrenergic","Dopaminergic","HPA Axis",
                       "BDNF","Neurogenesis","Kynurenine","Glutamate","Circadian"),
      `MDD Status` = c("↓ 5-HT synapse","↓ NE (LC)","↓ DA reward","↑ Cortisol (HPA dysreg.)",
                       "↓ BDNF (Val66Met)","↓ Hippocampal","↑ QUIN / IDO1",
                       "↑ Synaptic glu / ↓ AMPA","Disrupted CLOCK/PER"),
      `Drug Target` = c("SERT (SSRI/SNRI)","NET (SNRI/TCA)","DAT (Bupropion)",
                        "Antidepressant chronic","Ketamine/SSRI","Ketamine/mTOR",
                        "Anti-inflammatory","NMDA antagonist","Agomelatine"),
      check.names = FALSE
    )
    datatable(df, options = list(dom = "t", pageLength = 20),
              rownames = FALSE, class = "table-sm")
  })

  # ==========================================================================
  # TAB 2: Drug PK
  # ==========================================================================

  output$plot_pk_conc <- renderPlotly({
    d <- sim_data()
    # Downsample for display
    d_ds <- d[seq(1, nrow(d), by = 6), ]
    p <- ggplot(d_ds) +
      geom_line(aes(x = week, y = ESC_nM, color = "Escitalopram"), size = 1) +
      geom_line(aes(x = week, y = VEN_nM, color = "Venlafaxine"),  size = 1) +
      geom_line(aes(x = week, y = BUP_nM, color = "Bupropion"),    size = 1) +
      scale_color_manual(values = c(Escitalopram="#3498DB", Venlafaxine="#2ECC71",
                                     Bupropion="#E67E22")) +
      labs(x = "Time (weeks)", y = "Plasma Concentration (nM)",
           title = "Plasma Concentration–Time Profile", color = "Drug") +
      theme_bw(base_size = 11)
    ggplotly(p)
  })

  output$table_pk_ss <- renderDT({
    d <- sim_data()
    last <- tail(d, 1)
    df <- data.frame(
      Drug         = c("Escitalopram","Venlafaxine","Bupropion","Ketamine"),
      `Css (nM)`   = round(c(last$ESC_nM, last$VEN_nM, last$BUP_nM, last$KET_uM*1000), 1),
      `SERT % occ` = round(c(last$SERT_PCT, last$SERT_PCT, 0, 0), 1),
      `NET % occ`  = round(c(0, last$NET_PCT, last$NET_PCT*0.5, 0), 1),
      check.names  = FALSE
    )
    datatable(df, options = list(dom = "t"), rownames = FALSE, class = "table-sm")
  })

  output$plot_sert <- renderPlotly({
    d <- sim_data()
    d_ds <- d[seq(1, nrow(d), by = 6), ]
    p <- ggplot(d_ds, aes(x = week, y = SERT_PCT)) +
      geom_line(color = "#3498DB", size = 1.2) +
      geom_hline(yintercept = 80, linetype = "dashed", color = "gray50") +
      annotate("text", x = max(d_ds$week)*0.8, y = 81, label = "80% target", size = 3) +
      labs(x = "Time (weeks)", y = "SERT Occupancy (%)",
           title = "SERT Occupancy Over Time") +
      theme_bw()
    ggplotly(p)
  })

  output$plot_net <- renderPlotly({
    d <- sim_data()
    d_ds <- d[seq(1, nrow(d), by = 6), ]
    p <- ggplot(d_ds, aes(x = week, y = NET_PCT)) +
      geom_line(color = "#9B59B6", size = 1.2) +
      labs(x = "Time (weeks)", y = "NET Occupancy (%)",
           title = "NET Occupancy Over Time") +
      theme_bw()
    ggplotly(p)
  })

  output$plot_sert_dose <- renderPlotly({
    doses <- seq(0, 40, by = 2)
    # SERT occupancy at steady state for escitalopram
    # Using Hill equation: Css ≈ dose * F / (CL) ... simplified
    css_nM <- doses * 0.80 / 37 * 1e6 / 324.39 * (1/0.46) * 0.5  # rough Css
    sert <- 0.95 * css_nM^1.5 / (1.1^1.5 + css_nM^1.5)
    df <- data.frame(dose = doses, SERT_pct = sert * 100)
    p <- ggplot(df, aes(x = dose, y = SERT_pct)) +
      geom_line(color = "#3498DB", size = 1.2) +
      geom_point(color = "#3498DB", size = 2) +
      geom_hline(yintercept = 80, linetype = "dashed", color = "gray50") +
      labs(x = "Escitalopram Dose (mg/day)", y = "SERT Occupancy (%)",
           title = "SERT Occupancy vs Dose (Escitalopram)") +
      theme_bw()
    ggplotly(p)
  })

  output$table_pk_params <- renderDT({
    df <- data.frame(
      Drug = c("Escitalopram","Venlafaxine","Bupropion","Ketamine","Amitriptyline"),
      `t½ (h)` = c(27, 5, 21, 0.4, 16),
      `CL (L/h)` = c(37, 90, 200, 122, 62),
      `Vd (L)` = c(1090, 500, 3000, 210, 2500),
      `F (%)` = c(80, 45, 85, "IV", 48),
      `Ki SERT (nM)` = c(1.1, 7.5, ">1000", "-", 35),
      `Ki NET (nM)` = c("-", 2.7, 52, "-", 43),
      check.names = FALSE
    )
    datatable(df, options = list(dom = "t"), rownames = FALSE, class = "table-sm")
  })

  # ==========================================================================
  # TAB 3: Neurotransmitter Dynamics
  # ==========================================================================

  output$plot_5ht <- renderPlotly({
    d <- sim_data()
    d_ds <- d[seq(1, nrow(d), by = 6), ]
    p <- ggplot(d_ds, aes(x = week, y = X5HT)) +
      geom_line(color = "#3498DB", size = 1.2) +
      geom_hline(yintercept = 1.5, linetype = "dashed", color = "green4") +
      annotate("text", x = max(d_ds$week)*0.7, y = 1.55,
               label = "Healthy baseline", size = 3, color = "green4") +
      labs(x = "Weeks", y = "[5-HT] Synaptic (nM)",
           title = "Synaptic Serotonin (5-HT)") + theme_bw()
    ggplotly(p)
  })

  output$plot_ne <- renderPlotly({
    d <- sim_data()
    d_ds <- d[seq(1, nrow(d), by = 6), ]
    p <- ggplot(d_ds, aes(x = week, y = NE)) +
      geom_line(color = "#E74C3C", size = 1.2) +
      geom_hline(yintercept = 1.0, linetype = "dashed", color = "green4") +
      labs(x = "Weeks", y = "[NE] Synaptic (nM)",
           title = "Synaptic Norepinephrine (NE)") + theme_bw()
    ggplotly(p)
  })

  output$plot_da <- renderPlotly({
    d <- sim_data()
    d_ds <- d[seq(1, nrow(d), by = 6), ]
    p <- ggplot(d_ds, aes(x = week, y = DA)) +
      geom_line(color = "#9B59B6", size = 1.2) +
      geom_hline(yintercept = 0.8, linetype = "dashed", color = "green4") +
      labs(x = "Weeks", y = "[DA] Synaptic (nM)",
           title = "Synaptic Dopamine (DA)") + theme_bw()
    ggplotly(p)
  })

  output$plot_mono_bar <- renderPlotly({
    d <- sim_data()
    last <- tail(d, 1)
    df <- data.frame(
      NT    = c("5-HT","NE","DA"),
      Value = c(last$X5HT, last$NE, last$DA),
      Normal = c(1.5, 1.0, 0.8)
    ) %>% pivot_longer(-NT, names_to = "Type", values_to = "Conc")
    p <- ggplot(df, aes(x = NT, y = Conc, fill = Type)) +
      geom_bar(stat = "identity", position = "dodge") +
      scale_fill_manual(values = c(Value = "#3498DB", Normal = "#2ECC71")) +
      labs(x = "Neurotransmitter", y = "Concentration (nM)",
           title = "Monoamine Levels at Week 8") + theme_bw()
    ggplotly(p)
  })

  output$table_neuro_ts <- renderDT({
    d <- sim_data()
    df <- d %>%
      filter(week %in% c(0, 1, 2, 4, 6, 8)) %>%
      select(week, X5HT, NE, DA, SERT_PCT, NET_PCT) %>%
      mutate(across(where(is.numeric), ~round(.x, 3)))
    colnames(df) <- c("Week","5-HT (nM)","NE (nM)","DA (nM)","SERT occ (%)","NET occ (%)")
    datatable(df, options = list(dom = "t"), rownames = FALSE, class = "table-sm")
  })

  # ==========================================================================
  # TAB 4: HPA Axis & Neuroinflammation
  # ==========================================================================

  output$plot_cortisol <- renderPlotly({
    d <- sim_data()
    d_ds <- d[seq(1, nrow(d), by = 6), ]
    p <- ggplot(d_ds, aes(x = week, y = CORT)) +
      geom_line(color = "#E67E22", size = 1.2) +
      geom_hline(yintercept = 15, linetype = "dashed", color = "green4") +
      labs(x = "Weeks", y = "Cortisol (nmol/L)",
           title = "Cortisol Trajectory (HPA Axis)") + theme_bw()
    ggplotly(p)
  })

  output$plot_il6 <- renderPlotly({
    d <- sim_data()
    d_ds <- d[seq(1, nrow(d), by = 6), ]
    p <- ggplot(d_ds, aes(x = week, y = IL6)) +
      geom_line(color = "#E74C3C", size = 1.2) +
      geom_hline(yintercept = 0.625, linetype = "dashed", color = "green4") +
      labs(x = "Weeks", y = "IL-6 (pg/mL)",
           title = "IL-6 Neuroinflammation") + theme_bw()
    ggplotly(p)
  })

  output$plot_bdnf <- renderPlotly({
    d <- sim_data()
    d_ds <- d[seq(1, nrow(d), by = 6), ]
    p <- ggplot(d_ds, aes(x = week, y = BDNF)) +
      geom_line(color = "#27AE60", size = 1.2) +
      geom_hline(yintercept = 2.0, linetype = "dashed", color = "gray50") +
      labs(x = "Weeks", y = "BDNF (ng/mL)",
           title = "BDNF Dynamics") + theme_bw()
    ggplotly(p)
  })

  output$plot_neuro_index <- renderPlotly({
    d <- sim_data()
    d_ds <- d[seq(1, nrow(d), by = 6), ]
    p <- ggplot(d_ds) +
      geom_line(aes(x = week, y = MTOR, color = "mTOR Activity"), size = 1.2) +
      geom_line(aes(x = week, y = NEURO, color = "Neurogenesis Index"), size = 1.2) +
      scale_color_manual(values = c("mTOR Activity"="#9B59B6","Neurogenesis Index"="#1ABC9C")) +
      labs(x = "Weeks", y = "Index (normalized)",
           title = "mTOR & Neurogenesis", color = "") + theme_bw()
    ggplotly(p)
  })

  output$plot_kyn <- renderPlotly({
    d <- sim_data()
    d_ds <- d[seq(1, nrow(d), by = 6), ]
    # KYN:TRP proxy = inversely proportional to 5-HT, driven by IL-6
    d_ds$KYN_ratio <- (d_ds$IL6 / 0.625) * (1.5 / (d_ds$X5HT + 0.1))
    p <- ggplot(d_ds, aes(x = week, y = KYN_ratio)) +
      geom_line(color = "#C0392B", size = 1.2) +
      geom_hline(yintercept = 1.0, linetype = "dashed", color = "green4") +
      labs(x = "Weeks", y = "KYN:TRP Ratio (normalized)",
           title = "Kynurenine Pathway Activity") + theme_bw()
    ggplotly(p)
  })

  output$table_biomarker <- renderDT({
    d <- sim_data()
    last <- tail(d, 1)
    df <- data.frame(
      Biomarker    = c("Cortisol","IL-6","BDNF","mTOR Activity","Neurogenesis"),
      `Week 0`     = c(round(d$CORT[1],1), round(d$IL6[1],3),
                       round(d$BDNF[1],2), round(d$MTOR[1],3), round(d$NEURO[1],3)),
      `Week 8`     = c(round(last$CORT,1), round(last$IL6,3),
                       round(last$BDNF,2), round(last$MTOR,3), round(last$NEURO,3)),
      Units        = c("nmol/L","pg/mL","ng/mL","index","index"),
      Normal       = c("~15","~0.6","~2.0","~1.0","~1.0"),
      check.names  = FALSE
    )
    datatable(df, options = list(dom = "t"), rownames = FALSE, class = "table-sm")
  })

  # ==========================================================================
  # TAB 5: Clinical Endpoints
  # ==========================================================================

  output$vb_hdrs_wk8 <- renderValueBox({
    d <- sim_data()
    val <- round(tail(d$HDRS, 1), 1)
    col <- if (val <= 7) "green" else if (val <= 17) "yellow" else "red"
    valueBox(val, "HDRS-17 at Week 8", icon = icon("chart-line"), color = col)
  })

  output$vb_phq9_wk8 <- renderValueBox({
    d <- sim_data()
    val <- round(tail(d$PHQ9, 1), 1)
    col <- if (val < 5) "green" else if (val < 10) "yellow" else "red"
    valueBox(val, "PHQ-9 at Week 8", icon = icon("clipboard-check"), color = col)
  })

  output$vb_response <- renderValueBox({
    d <- sim_data()
    resp <- tail(d$RESPONSE, 1)
    col  <- if (resp == 1) "green" else "red"
    label <- if (resp == 1) "Yes (≥50% reduction)" else "Not Yet"
    valueBox(label, "Clinical Response", icon = icon("check-circle"), color = col)
  })

  output$vb_remission <- renderValueBox({
    d <- sim_data()
    rem  <- tail(d$REMISSION, 1)
    col  <- if (rem == 1) "green" else "orange"
    label <- if (rem == 1) "Yes (PHQ-9 < 5)" else "Not Yet"
    valueBox(label, "Remission", icon = icon("star"), color = col)
  })

  output$plot_hdrs <- renderPlotly({
    d <- sim_data()
    d_ds <- d[seq(1, nrow(d), by = 6), ]
    p <- ggplot(d_ds, aes(x = week, y = HDRS)) +
      geom_ribbon(aes(ymin = HDRS * 0.9, ymax = HDRS * 1.1),
                  fill = "#3498DB", alpha = 0.2) +
      geom_line(color = "#2980B9", size = 1.3) +
      geom_hline(yintercept = 7, linetype = "dashed", color = "green4") +
      geom_hline(yintercept = input$hdrs_base * 0.5, linetype = "dotted", color = "orange") +
      annotate("text", x = max(d_ds$week)*0.8, y = 6.0,
               label = "Remission (HDRS≤7)", size = 3, color = "green4") +
      annotate("text", x = max(d_ds$week)*0.8, y = input$hdrs_base*0.5 - 1,
               label = "Response (50%↓)", size = 3, color = "orange") +
      labs(x = "Time (weeks)", y = "HDRS-17 Score",
           title = "HDRS-17 Trajectory") +
      scale_x_continuous(breaks = seq(0, input$weeks_sim, by = 2)) +
      theme_bw(base_size = 12)
    ggplotly(p)
  })

  output$plot_resp_timeline <- renderPlot({
    d <- sim_data()
    d$status <- ifelse(d$REMISSION == 1, "Remission",
                       ifelse(d$RESPONSE == 1, "Response", "No Response"))
    # Summarise by week
    wk_df <- d %>%
      group_by(week = round(week)) %>%
      summarise(status = last(status), .groups = "drop") %>%
      filter(week <= input$weeks_sim)
    cols <- c(Remission = "#27AE60", Response = "#F39C12", "No Response" = "#E74C3C")
    ggplot(wk_df, aes(x = week, y = 1, fill = status)) +
      geom_tile(height = 0.8, color = "white") +
      scale_fill_manual(values = cols) +
      scale_x_continuous(breaks = 0:input$weeks_sim) +
      labs(x = "Week", y = "", title = "Response Status by Week", fill = "Status") +
      theme_bw(base_size = 11) +
      theme(axis.text.y = element_blank(), axis.ticks.y = element_blank())
  })

  output$plot_phq9_madrs <- renderPlotly({
    d <- sim_data()
    d_ds <- d[seq(1, nrow(d), by = 6), ]
    p <- ggplot(d_ds) +
      geom_line(aes(x = week, y = PHQ9,  color = "PHQ-9"),  size = 1.2) +
      geom_line(aes(x = week, y = MADRS/3, color = "MADRS/3"), size = 1.2) +
      scale_color_manual(values = c("PHQ-9"="#E74C3C","MADRS/3"="#9B59B6")) +
      labs(x = "Weeks", y = "Score", title = "PHQ-9 & MADRS Trajectories",
           color = "Scale") + theme_bw()
    ggplotly(p)
  })

  output$table_clinical_weekly <- renderDT({
    d <- sim_data()
    df <- d %>%
      filter(week %in% 0:min(8, input$weeks_sim)) %>%
      group_by(week = round(week)) %>%
      summarise(HDRS = round(mean(HDRS), 1),
                PHQ9 = round(mean(PHQ9), 1),
                MADRS = round(mean(MADRS), 1),
                Response = max(RESPONSE),
                Remission = max(REMISSION),
                .groups = "drop") %>%
      distinct(week, .keep_all = TRUE)
    colnames(df) <- c("Week","HDRS-17","PHQ-9","MADRS","Response","Remission")
    datatable(df, options = list(dom = "t", pageLength = 12),
              rownames = FALSE, class = "table-sm")
  })

  # ==========================================================================
  # TAB 6: Scenario Comparison
  # ==========================================================================

  compare_data <- reactive({
    all_s <- run_all_scenarios()
    # Filter to selected
    selected <- input$compare_scenarios
    bind_rows(
      lapply(selected, function(nm) {
        if (!is.null(all_s[[nm]])) {
          all_s[[nm]] %>% mutate(scenario = nm)
        }
      })
    )
  })

  output$plot_compare_hdrs <- renderPlotly({
    d <- compare_data()
    if (nrow(d) == 0) return(NULL)
    d_ds <- d %>% filter(time %% 12 == 0)
    p <- ggplot(d_ds, aes(x = week, y = HDRS, color = scenario)) +
      geom_line(size = 1.1) +
      geom_hline(yintercept = 7, linetype = "dashed", color = "gray40") +
      scale_color_manual(values = SCEN_COLORS) +
      scale_x_continuous(breaks = 0:8) +
      labs(x = "Time (weeks)", y = "HDRS-17 Score",
           title = "HDRS-17 Comparison Across Treatment Arms",
           color = "Scenario") +
      theme_bw(base_size = 11) +
      theme(legend.position = "right", legend.text = element_text(size = 8))
    ggplotly(p) %>% layout(legend = list(font = list(size = 9)))
  })

  output$plot_forest <- renderPlot({
    all_s <- run_all_scenarios()
    summary_df <- lapply(names(all_s), function(nm) {
      d <- all_s[[nm]]
      hdrs0 <- d$HDRS[1]
      hdrs8 <- tail(d$HDRS, 1)
      change <- hdrs0 - hdrs8
      pct    <- change / hdrs0 * 100
      data.frame(scenario = nm, change = change, pct = pct,
                 ci_lo = pct - 8, ci_hi = pct + 8)
    }) %>% bind_rows()

    ggplot(summary_df, aes(x = pct, y = reorder(scenario, pct))) +
      geom_vline(xintercept = 50, linetype = "dashed", color = "orange") +
      geom_vline(xintercept = 0,  linetype = "solid",  color = "gray50") +
      geom_errorbarh(aes(xmin = ci_lo, xmax = ci_hi), height = 0.3,
                     color = "gray40") +
      geom_point(aes(color = pct >= 50), size = 4) +
      scale_color_manual(values = c("TRUE"="#27AE60","FALSE"="#E74C3C"),
                         labels = c("TRUE"="Response","FALSE"="Non-Response")) +
      labs(x = "HDRS-17 % Reduction (Week 8)",
           y = NULL, title = "Forest Plot: Effect Size (Week 8 vs Baseline)",
           color = "Outcome") +
      theme_bw(base_size = 11)
  })

  output$table_compare_summary <- renderDT({
    all_s <- run_all_scenarios()
    df <- lapply(names(all_s), function(nm) {
      d <- all_s[[nm]]
      hdrs0 <- round(d$HDRS[1], 1)
      hdrs8 <- round(tail(d$HDRS, 1), 1)
      bdnf8 <- round(tail(d$BDNF, 1), 2)
      sert8 <- round(tail(d$SERT_PCT, 1), 0)
      pct   <- round((hdrs0 - hdrs8) / hdrs0 * 100, 1)
      data.frame(Scenario = nm,
                 `HDRS Wk0` = hdrs0, `HDRS Wk8` = hdrs8,
                 `% Reduction` = pct,
                 Response = ifelse(pct >= 50, "Yes", "No"),
                 `BDNF Wk8` = bdnf8,
                 `SERT% Wk8` = sert8,
                 check.names = FALSE)
    }) %>% bind_rows()
    datatable(df, options = list(dom = "t"), rownames = FALSE, class = "table-sm",
              selection = "none") %>%
      formatStyle("% Reduction",
        background = styleColorBar(c(0, 100), "#2ECC71"),
        backgroundSize = "100% 70%", backgroundRepeat = "no-repeat",
        backgroundPosition = "center")
  })

  output$plot_compare_bdnf <- renderPlotly({
    d <- compare_data()
    if (nrow(d) == 0) return(NULL)
    d_ds <- d %>% filter(time %% 12 == 0)
    p <- ggplot(d_ds, aes(x = week, y = BDNF, color = scenario)) +
      geom_line(size = 1.1) +
      geom_hline(yintercept = 2.0, linetype = "dashed", color = "gray40") +
      scale_color_manual(values = SCEN_COLORS) +
      labs(x = "Weeks", y = "BDNF (ng/mL)",
           title = "BDNF Comparison", color = "Scenario") +
      theme_bw()
    ggplotly(p)
  })

  output$plot_time_response <- renderPlot({
    all_s <- run_all_scenarios()
    ttr_df <- lapply(names(all_s), function(nm) {
      d <- all_s[[nm]]
      hdrs0 <- d$HDRS[1]
      ttr_idx <- which(d$HDRS <= hdrs0 * 0.5)[1]
      ttr_wk  <- if (!is.na(ttr_idx)) d$week[ttr_idx] else NA
      data.frame(scenario = nm, ttr_weeks = ttr_wk)
    }) %>% bind_rows()

    ggplot(ttr_df, aes(x = reorder(scenario, ttr_weeks, na.last = TRUE),
                        y = ttr_weeks, fill = !is.na(ttr_weeks))) +
      geom_bar(stat = "identity") +
      geom_hline(yintercept = 2, linetype = "dashed", color = "orange") +
      scale_fill_manual(values = c("TRUE"="#2ECC71","FALSE"="#BDC3C7"),
                        labels = c("TRUE"="Response Achieved","FALSE"="No Response")) +
      labs(x = NULL, y = "Time to Response (weeks)",
           title = "Time to ≥50% HDRS Reduction", fill = "") +
      coord_flip() +
      theme_bw(base_size = 11)
  })
}

# =============================================================================
# RUN APP
# =============================================================================

shinyApp(ui = ui, server = server)
