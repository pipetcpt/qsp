## ============================================================================
## Gastric Cancer QSP — Interactive Shiny Dashboard
## Version: 1.0  |  Date: 2026-06-23
##
## Tabs:
##   1. 환자 프로파일 (Patient Profile)
##   2. 약물 PK (Drug Pharmacokinetics)
##   3. 종양 동태 (Tumor Growth Dynamics)
##   4. 임상 엔드포인트 (Clinical Endpoints — OS/PFS/ORR)
##   5. 치료 시나리오 비교 (Treatment Scenario Comparison)
##   6. 바이오마커 분석 (Biomarker Analysis)
## ============================================================================

library(shiny)
library(shinydashboard)
library(ggplot2)
library(dplyr)
library(tidyr)
library(plotly)
library(DT)
library(survival)
library(survminer)
library(RColorBrewer)

## ============================================================================
## Helper Functions (standalone PK/PD without mrgsolve dependency)
## ============================================================================

# Two-compartment IV PK model (analytic solution)
two_cmt_iv <- function(dose, CL, V1, V2, Q, times, tinf = 0.5) {
  k10 <- CL / V1
  k12 <- Q / V1
  k21 <- Q / V2
  ksum <- k10 + k12 + k21
  alpha <- (ksum + sqrt(ksum^2 - 4 * k10 * k21)) / 2
  beta  <- (ksum - sqrt(ksum^2 - 4 * k10 * k21)) / 2
  A <- (alpha - k21) / (V1 * (alpha - beta))
  B <- (k21 - beta)  / (V1 * (alpha - beta))

  conc <- sapply(times, function(t) {
    if (t <= 0) return(0)
    if (t <= tinf) {
      # During infusion
      rate <- dose / tinf
      Ci <- rate * (A/alpha * (1 - exp(-alpha*t)) + B/beta * (1 - exp(-beta*t)))
    } else {
      # Post-infusion
      t_post <- t - tinf
      rate <- dose / tinf
      Ci_end <- rate * (A/alpha * (1 - exp(-alpha*tinf)) + B/beta * (1 - exp(-beta*tinf)))
      Ci <- Ci_end * (A/(A+B) * exp(-alpha*t_post) + B/(A+B) * exp(-beta*t_post))
    }
    max(0, Ci)
  })
  conc
}

# Multiple dose IV (superposition)
multi_dose_iv <- function(doses_df, CL, V1, V2, Q, times, tinf = 0.5) {
  conc_total <- rep(0, length(times))
  for (i in seq_len(nrow(doses_df))) {
    dose_amt  <- doses_df$amt[i]
    dose_time <- doses_df$time[i]
    shifted_times <- pmax(times - dose_time, 0)
    conc_total <- conc_total + two_cmt_iv(dose_amt, CL, V1, V2, Q, shifted_times, tinf)
  }
  conc_total
}

# Oral 1-compartment (capecitabine)
oral_1cmt <- function(dose, ka, F_oral, CL, V, times, tlag = 0) {
  k <- CL / V
  conc <- sapply(times, function(t) {
    t_adj <- t - tlag
    if (t_adj <= 0) return(0)
    F_oral * dose * ka / (V * (ka - k)) * (exp(-k * t_adj) - exp(-ka * t_adj))
  })
  pmax(0, conc)
}

# Simplified Simeoni TGI model
sim_simeoni_tgi <- function(tv0, lambda0, lambda1, k1, k2, E_total, times) {
  # ODE via Euler method
  dt     <- diff(times)
  n      <- length(times)
  TV     <- numeric(n)
  TV_D1  <- numeric(n)
  TV_D2  <- numeric(n)
  TV[1]  <- tv0
  TV_D1[1] <- 0
  TV_D2[1] <- 0

  for (i in seq_len(n - 1)) {
    dtv_total <- TV[i] + TV_D1[i] + TV_D2[i]
    if (dtv_total > 0) {
      growth  <- 2 * lambda0 * lambda1 / (lambda1 + 2 * lambda0 * TV[i]) * TV[i]
    } else {
      growth <- 0
    }
    E_t <- if(length(E_total) == n) E_total[i] else E_total
    TV[i+1]    <- TV[i]   + dt[i] * (growth - k1 * E_t * TV[i])
    TV_D1[i+1] <- TV_D1[i] + dt[i] * (k1 * E_t * TV[i] - k2 * TV_D1[i])
    TV_D2[i+1] <- TV_D2[i] + dt[i] * (k2 * TV_D1[i] - k2 * TV_D2[i])
    TV[i+1] <- max(0, TV[i+1])
  }
  TV + TV_D1 + TV_D2
}

## ============================================================================
## Simulate treatment scenarios for the app
## ============================================================================
simulate_scenario <- function(scenario_name, params, times) {
  with(params, {

    # Drug PK
    if (scenario_name == "Trastuzumab + FOLFOX (HER2+ 1L)") {
      tras_doses <- data.frame(
        time = c(0, 21, 42, 63, 84, 105, 126),
        amt  = c(8, 6, 6, 6, 6, 6, 6) * weight / 148000 * 1000 / V1_tras
      )
      conc_tras <- multi_dose_iv(tras_doses, CL_tras, V1_tras, V2_tras, Q_tras, times)
      E_tras    <- Emax_tras * conc_tras^n_tras / (EC50_tras^n_tras + conc_tras^n_tras)
      conc_chemo <- oral_1cmt(cape_dose, ka_cape, F_cape, CL_cape, V_cape, times)
      E_chemo   <- Emax_chemo * conc_chemo / (EC50_FU + conc_chemo)
      E_total   <- 1 - (1 - E_tras) * (1 - E_chemo)
      E_immune  <- rep(0, length(times))
      conc_nivo <- rep(0, length(times))
      conc_ramu <- rep(0, length(times))

    } else if (scenario_name == "Ramucirumab + Paclitaxel (2L)") {
      ramu_doses <- data.frame(
        time = seq(0, 112, by = 14),
        amt  = 8 * weight / 147000 * 1000 / V1_ramu
      )
      conc_ramu <- multi_dose_iv(ramu_doses, CL_ramu, V1_ramu, V2_ramu, Q_ramu, times)
      E_ramu    <- Emax_ramu * conc_ramu / (EC50_ramu + conc_ramu)
      conc_chemo <- oral_1cmt(cape_dose * 0.8, ka_cape, F_cape, CL_cape, V_cape, times)
      E_chemo   <- Emax_chemo * conc_chemo / (EC50_FU + conc_chemo)
      E_total   <- 1 - (1 - E_ramu) * (1 - E_chemo)
      conc_tras <- rep(0, length(times))
      E_immune  <- rep(0, length(times))
      conc_nivo <- rep(0, length(times))

    } else if (scenario_name == "Nivolumab + Chemo (CPS≥5, 1L)") {
      nivo_doses <- data.frame(
        time = seq(0, 168, by = 21),
        amt  = 360 / 146000 * 1000 / V1_nivo
      )
      conc_nivo <- multi_dose_iv(nivo_doses, CL_nivo, V1_nivo, V2_nivo, Q_nivo, times)
      PD1_occ   <- conc_nivo / (conc_nivo + EC50_nivo_CD8)
      E_immune  <- Emax_nivo * PD1_occ
      conc_chemo <- oral_1cmt(cape_dose, ka_cape, F_cape, CL_cape, V_cape, times)
      E_chemo   <- Emax_chemo * conc_chemo / (EC50_FU + conc_chemo)
      E_total   <- 1 - (1 - E_immune) * (1 - E_chemo)
      conc_tras <- rep(0, length(times))
      conc_ramu <- rep(0, length(times))

    } else if (scenario_name == "T-DXd (HER2+ 2L)") {
      tdxd_doses <- data.frame(
        time = seq(0, 168, by = 21),
        amt  = 6.4 * weight / 184000 * 1000 / V1_tdxd
      )
      conc_tdxd <- multi_dose_iv(tdxd_doses, CL_tdxd, V1_tdxd, V1_tdxd * 0.9, 0.1, times)
      E_tdxd    <- Emax_tdxd * conc_tdxd / (EC50_dxd + conc_tdxd)
      E_total   <- E_tdxd
      conc_tras <- rep(0, length(times))
      conc_ramu <- rep(0, length(times))
      conc_nivo <- rep(0, length(times))
      E_immune  <- rep(0, length(times))
      conc_chemo <- rep(0, length(times))

    } else if (scenario_name == "Zolbetuximab + mFOLFOX6 (CLDN18.2+ 1L)") {
      zolbe_doses <- data.frame(
        time = c(0, seq(21, 168, by = 21)),
        amt  = c(800, rep(600, 8)) * 1.7 / 148000 * 1000 / V1_zolbe
      )
      conc_zolbe <- multi_dose_iv(zolbe_doses, CL_zolbe, V1_zolbe, V1_zolbe * 0.8, 0.2, times)
      E_zolbe    <- Emax_zolbe * conc_zolbe / (EC50_zolbe_kill + conc_zolbe)
      conc_chemo <- oral_1cmt(cape_dose, ka_cape, F_cape, CL_cape, V_cape, times)
      E_chemo   <- Emax_chemo * conc_chemo / (EC50_FU + conc_chemo)
      E_total   <- 1 - (1 - E_zolbe) * (1 - E_chemo)
      conc_tras <- rep(0, length(times))
      conc_ramu <- rep(0, length(times))
      conc_nivo <- rep(0, length(times))
      E_immune  <- rep(0, length(times))

    } else {  # FLOT
      conc_chemo <- sapply(times, function(t) {
        cycle <- floor(t / 14)
        t_in_cycle <- t - cycle * 14
        if (cycle >= 8) return(0)
        oral_1cmt(cape_dose, ka_cape, F_cape, CL_cape, V_cape, pmax(t_in_cycle, 0))
      })
      E_chemo   <- Emax_chemo * conc_chemo / (EC50_FU + conc_chemo)
      E_total   <- E_chemo * 1.2  # FLOT multi-drug synergy
      E_total   <- pmin(E_total, 0.95)
      conc_tras <- rep(0, length(times))
      conc_ramu <- rep(0, length(times))
      conc_nivo <- rep(0, length(times))
      E_immune  <- rep(0, length(times))
    }

    tumor_vol <- sim_simeoni_tgi(tv0, lambda0, lambda1, k1_tgi, k2_tgi, E_total, times)
    cea       <- pmax(10, CEA0 * tumor_vol / tv0 * 0.8 + 5)
    cd8       <- pmax(0.1, 0.67 * (1 + E_immune * 2) * exp(-0.001 * times))

    list(
      times      = times,
      tumor_vol  = tumor_vol,
      conc_tras  = if(exists("conc_tras")) conc_tras else rep(0, length(times)),
      conc_ramu  = if(exists("conc_ramu")) conc_ramu else rep(0, length(times)),
      conc_nivo  = if(exists("conc_nivo")) conc_nivo else rep(0, length(times)),
      conc_chemo = if(exists("conc_chemo")) conc_chemo else rep(0, length(times)),
      E_total    = E_total,
      E_immune   = E_immune,
      cea        = cea,
      cd8        = cd8
    )
  })
}

## ============================================================================
## UI
## ============================================================================
ui <- dashboardPage(
  skin = "red",  # Gastric cancer theme

  dashboardHeader(
    title = "Gastric Cancer QSP Dashboard",
    titleWidth = 320
  ),

  dashboardSidebar(
    width = 280,
    sidebarMenu(
      id = "tabs",
      menuItem("환자 프로파일", tabName = "profile",
               icon = icon("user-md")),
      menuItem("약물 PK", tabName = "pk",
               icon = icon("pills")),
      menuItem("종양 동태", tabName = "tumor",
               icon = icon("chart-line")),
      menuItem("임상 엔드포인트", tabName = "endpoints",
               icon = icon("heartbeat")),
      menuItem("치료 시나리오 비교", tabName = "scenarios",
               icon = icon("balance-scale")),
      menuItem("바이오마커 분석", tabName = "biomarkers",
               icon = icon("microscope"))
    ),
    hr(),
    h5("Quick Settings", style = "color:white; padding-left:15px;"),
    sliderInput("sim_duration", "Simulation Duration (days):",
                min = 84, max = 365, value = 252, step = 14),
    selectInput("her2_status", "HER2 Status:",
                choices = c("HER2-positive (IHC 3+)" = "pos3",
                            "HER2 equivocal (IHC 2+/FISH+)" = "pos2",
                            "HER2-low (IHC 1+)" = "low",
                            "HER2-negative" = "neg"),
                selected = "pos3"),
    selectInput("line_therapy", "Line of Therapy:",
                choices = c("1st Line" = "1L",
                            "2nd Line" = "2L",
                            "3rd Line+" = "3L"),
                selected = "1L")
  ),

  dashboardBody(
    tags$head(
      tags$style(HTML("
        .skin-red .main-header .navbar { background-color: #C62828; }
        .skin-red .main-header .logo { background-color: #B71C1C; }
        .skin-red .sidebar { background-color: #212121; }
        .box-title { font-weight: bold; }
        .value-box .inner h3 { font-size: 28px; }
      "))
    ),

    tabItems(

      # ====================================================================
      # TAB 1: Patient Profile
      # ====================================================================
      tabItem(tabName = "profile",
        fluidRow(
          box(title = "Patient Characteristics", status = "danger", solidHeader = TRUE,
              width = 6,
              fluidRow(
                column(6,
                  numericInput("age", "Age (years):", value = 62, min = 18, max = 90),
                  selectInput("sex", "Sex:",
                              choices = c("Male" = "M", "Female" = "F"), selected = "M"),
                  numericInput("weight", "Weight (kg):", value = 65, min = 30, max = 150),
                  numericInput("bsa", "BSA (m²):", value = 1.70, min = 1.0, max = 3.0, step = 0.05)
                ),
                column(6,
                  selectInput("stage", "Disease Stage:",
                              choices = c("III (localized)" = "III",
                                          "IV (metastatic)" = "IV"),
                              selected = "IV"),
                  selectInput("lauren", "Lauren Classification:",
                              choices = c("Intestinal" = "int",
                                          "Diffuse" = "diff",
                                          "Mixed" = "mixed"),
                              selected = "int"),
                  selectInput("ecog", "ECOG PS:",
                              choices = c("0", "1", "2"), selected = "1"),
                  numericInput("prior_lines", "Prior Lines of Therapy:",
                               value = 0, min = 0, max = 5)
                )
              )
          ),
          box(title = "Molecular Biomarkers", status = "warning", solidHeader = TRUE,
              width = 6,
              fluidRow(
                column(6,
                  selectInput("her2_ihc", "HER2 IHC Score:",
                              choices = c("0", "1+", "2+", "3+"), selected = "3+"),
                  selectInput("fish", "HER2 FISH:",
                              choices = c("Positive (ratio≥2.0)" = "pos",
                                          "Negative" = "neg",
                                          "Not done" = "NA"),
                              selected = "pos"),
                  numericInput("cps", "PD-L1 CPS Score:", value = 8, min = 0, max = 100),
                  selectInput("msi", "MSI Status:",
                              choices = c("MSI-High" = "msi_h",
                                          "MSS / MSI-Low" = "mss"),
                              selected = "mss")
                ),
                column(6,
                  selectInput("cldn182", "CLDN18.2 IHC:",
                              choices = c("≥2+ in ≥75% cells (positive)" = "pos",
                                          "Negative / insufficient" = "neg"),
                              selected = "neg"),
                  numericInput("tmb", "TMB (mut/Mb):", value = 5, min = 0, max = 100),
                  selectInput("ebv", "EBV Status:",
                              choices = c("EBER-positive" = "pos",
                                          "EBER-negative" = "neg"),
                              selected = "neg"),
                  numericInput("cea_baseline", "Baseline CEA (ng/mL):", value = 45, min = 0)
                )
              )
          )
        ),
        fluidRow(
          box(title = "Molecular Subtype & Treatment Eligibility",
              status = "primary", solidHeader = TRUE, width = 12,
              fluidRow(
                valueBoxOutput("subtype_box", width = 3),
                valueBoxOutput("her2_box",    width = 3),
                valueBoxOutput("ici_box",     width = 3),
                valueBoxOutput("cldn_box",    width = 3)
              )
          )
        ),
        fluidRow(
          box(title = "TCGA Molecular Subtype Classifier (Lauren + Molecular)",
              status = "info", solidHeader = TRUE, width = 12,
              DTOutput("subtype_table")
          )
        )
      ),  # end Tab 1

      # ====================================================================
      # TAB 2: Drug PK
      # ====================================================================
      tabItem(tabName = "pk",
        fluidRow(
          box(title = "PK Parameters", status = "success", solidHeader = TRUE, width = 3,
              selectInput("pk_drug", "Select Drug:",
                          choices = c("Trastuzumab (8→6 mg/kg Q3W)",
                                      "Ramucirumab (8 mg/kg Q2W)",
                                      "Nivolumab (360 mg Q3W)",
                                      "Pembrolizumab (200 mg Q3W)",
                                      "T-DXd (6.4 mg/kg Q3W)",
                                      "Zolbetuximab (800→600 mg/m² Q3W)",
                                      "Capecitabine (1250 mg/m² BID D1-14)"),
                          selected = "Trastuzumab (8→6 mg/kg Q3W)"),
              numericInput("n_cycles_pk", "Number of Cycles:", value = 6, min = 1, max = 18),
              checkboxInput("show_therapeutic", "Show Therapeutic Window", TRUE),
              numericInput("c_trough", "Target Trough (nmol/L):", value = 1.0, min = 0, step = 0.1),
              actionButton("run_pk", "Simulate PK", class = "btn-success btn-block")
          ),
          box(title = "PK Profile", status = "success", solidHeader = TRUE, width = 9,
              plotlyOutput("pk_plot", height = "400px")
          )
        ),
        fluidRow(
          box(title = "PK Summary Statistics", status = "success", solidHeader = TRUE, width = 6,
              DTOutput("pk_summary_table")
          ),
          box(title = "PK Parameters Reference", status = "info", solidHeader = TRUE, width = 6,
              DTOutput("pk_params_table")
          )
        )
      ),  # end Tab 2

      # ====================================================================
      # TAB 3: Tumor Dynamics
      # ====================================================================
      tabItem(tabName = "tumor",
        fluidRow(
          box(title = "TGI Model Parameters", status = "danger", solidHeader = TRUE, width = 3,
              selectInput("tgi_scenario", "Treatment Scenario:",
                          choices = c("FLOT perioperative",
                                      "Trastuzumab + FOLFOX (HER2+ 1L)",
                                      "Ramucirumab + Paclitaxel (2L)",
                                      "Nivolumab + Chemo (CPS≥5, 1L)",
                                      "T-DXd (HER2+ 2L)",
                                      "Zolbetuximab + mFOLFOX6 (CLDN18.2+ 1L)"),
                          selected = "Trastuzumab + FOLFOX (HER2+ 1L)"),
              sliderInput("tv0", "Initial Tumor Volume (mm³):",
                          min = 100, max = 5000, value = 800, step = 100),
              sliderInput("lambda0", "Exponential Growth Rate (λ₀, 1/day):",
                          min = 0.005, max = 0.05, value = 0.0215, step = 0.001),
              sliderInput("lambda1", "Linear Growth Rate (λ₁, mm³/day):",
                          min = 0.1, max = 1.0, value = 0.385, step = 0.01),
              sliderInput("E_drug", "Drug Efficacy (Emax):",
                          min = 0.0, max = 1.0, value = 0.70, step = 0.05),
              actionButton("run_tgi", "Simulate Tumor Dynamics", class = "btn-danger btn-block")
          ),
          box(title = "Tumor Volume Over Time (Simeoni TGI Model)", status = "danger",
              solidHeader = TRUE, width = 9,
              plotlyOutput("tumor_plot", height = "400px")
          )
        ),
        fluidRow(
          box(title = "Waterfall Plot — Best % Change from Baseline",
              status = "warning", solidHeader = TRUE, width = 6,
              plotlyOutput("waterfall_plot", height = "300px")
          ),
          box(title = "Spider Plot — Individual Tumor Response",
              status = "warning", solidHeader = TRUE, width = 6,
              plotlyOutput("spider_plot", height = "300px")
          )
        )
      ),  # end Tab 3

      # ====================================================================
      # TAB 4: Clinical Endpoints
      # ====================================================================
      tabItem(tabName = "endpoints",
        fluidRow(
          box(title = "Endpoint Parameters", status = "primary", solidHeader = TRUE, width = 3,
              selectInput("endpoint_scenario", "Reference Trial:",
                          choices = c("CheckMate 649 (Nivo+Chemo vs Chemo, CPS≥5)",
                                      "ToGA (Tras+Chemo vs Chemo, HER2+)",
                                      "RAINBOW (Ramu+Pac vs Pac, 2L)",
                                      "SPOTLIGHT (Zolbe+mFOLFOX6 vs mFOLFOX6, CLDN18.2+)",
                                      "DESTINY-Gastric01 (T-DXd vs PC, HER2+ 2L)",
                                      "FLOT4 (FLOT vs ECF/ECX, periop)"),
                          selected = "CheckMate 649 (Nivo+Chemo vs Chemo, CPS≥5)"),
              numericInput("n_pts", "Virtual Patients (n):", value = 200, min = 50, max = 1000),
              sliderInput("os_hr", "OS Hazard Ratio:", min = 0.4, max = 1.0, value = 0.71, step = 0.01),
              sliderInput("pfs_hr", "PFS Hazard Ratio:", min = 0.4, max = 1.0, value = 0.68, step = 0.01),
              sliderInput("os_median_trt", "Median OS Treated (months):", min = 6, max = 30, value = 14.4),
              sliderInput("orr_trt", "ORR Treated (%):", min = 10, max = 90, value = 60),
              sliderInput("orr_ctrl", "ORR Control (%):", min = 10, max = 90, value = 45),
              actionButton("run_endpoints", "Simulate Endpoints", class = "btn-primary btn-block")
          ),
          box(title = "Kaplan-Meier — Overall Survival", status = "primary",
              solidHeader = TRUE, width = 9,
              plotOutput("km_os_plot", height = "420px")
          )
        ),
        fluidRow(
          box(title = "Kaplan-Meier — Progression-Free Survival", status = "info",
              solidHeader = TRUE, width = 6,
              plotOutput("km_pfs_plot", height = "350px")
          ),
          box(title = "Response Rates (ORR/DCR/CR)", status = "success",
              solidHeader = TRUE, width = 6,
              plotlyOutput("response_bar_plot", height = "350px")
          )
        )
      ),  # end Tab 4

      # ====================================================================
      # TAB 5: Treatment Scenario Comparison
      # ====================================================================
      tabItem(tabName = "scenarios",
        fluidRow(
          box(title = "Simulation Settings", status = "warning", solidHeader = TRUE, width = 3,
              checkboxGroupInput("selected_scenarios", "Select Scenarios:",
                choices = c("FLOT perioperative",
                            "Trastuzumab + FOLFOX (HER2+ 1L)",
                            "Ramucirumab + Paclitaxel (2L)",
                            "Nivolumab + Chemo (CPS≥5, 1L)",
                            "T-DXd (HER2+ 2L)",
                            "Zolbetuximab + mFOLFOX6 (CLDN18.2+ 1L)"),
                selected = c("Trastuzumab + FOLFOX (HER2+ 1L)",
                             "Nivolumab + Chemo (CPS≥5, 1L)",
                             "T-DXd (HER2+ 2L)")),
              sliderInput("tv0_comp", "Initial Tumor Volume (mm³):",
                          min = 200, max = 5000, value = 1000, step = 100),
              actionButton("run_compare", "Run Comparison", class = "btn-warning btn-block")
          ),
          box(title = "Tumor Volume Comparison", status = "warning",
              solidHeader = TRUE, width = 9,
              plotlyOutput("compare_tv_plot", height = "380px")
          )
        ),
        fluidRow(
          box(title = "Drug Exposure Comparison (AUC)", status = "info",
              solidHeader = TRUE, width = 6,
              plotlyOutput("compare_auc_plot", height = "300px")
          ),
          box(title = "Scenario Summary Table", status = "success",
              solidHeader = TRUE, width = 6,
              DTOutput("scenario_summary_table")
          )
        )
      ),  # end Tab 5

      # ====================================================================
      # TAB 6: Biomarker Analysis
      # ====================================================================
      tabItem(tabName = "biomarkers",
        fluidRow(
          box(title = "Biomarker Panel Settings", status = "purple", solidHeader = TRUE, width = 3,
              selectInput("bm_scenario", "Treatment Scenario:",
                          choices = c("Trastuzumab + FOLFOX (HER2+ 1L)",
                                      "Nivolumab + Chemo (CPS≥5, 1L)",
                                      "T-DXd (HER2+ 2L)",
                                      "Zolbetuximab + mFOLFOX6 (CLDN18.2+ 1L)"),
                          selected = "Nivolumab + Chemo (CPS≥5, 1L)"),
              numericInput("cea_base_bm", "Baseline CEA (ng/mL):", value = 45, min = 0),
              numericInput("ca199_base", "Baseline CA19-9 (U/mL):", value = 320, min = 0),
              numericInput("ctdna_base", "Baseline ctDNA (copies/mL):", value = 2500, min = 0),
              selectInput("her2_copy", "HER2 Copy Number:",
                          choices = c("≥6 (amplified)" = "amp",
                                      "4-6 (gain)" = "gain",
                                      "<4 (normal)" = "norm"),
                          selected = "amp"),
              actionButton("run_bm", "Analyze Biomarkers", class = "btn-block",
                           style = "background-color: #7B1FA2; color: white;")
          ),
          box(title = "Serum Biomarker Dynamics (CEA / CA19-9 / ctDNA)",
              status = "purple", solidHeader = TRUE, width = 9,
              plotlyOutput("bm_serum_plot", height = "400px")
          )
        ),
        fluidRow(
          box(title = "HER2 IHC Score Distribution", status = "danger",
              solidHeader = TRUE, width = 4,
              plotlyOutput("her2_ihc_plot", height = "280px")
          ),
          box(title = "CPS Score & Immunotherapy Eligibility", status = "warning",
              solidHeader = TRUE, width = 4,
              plotlyOutput("cps_plot", height = "280px")
          ),
          box(title = "TMB / MSI / EBV Molecular Subtype",
              status = "info", solidHeader = TRUE, width = 4,
              plotlyOutput("molecular_subtype_plot", height = "280px")
          )
        ),
        fluidRow(
          box(title = "Biomarker Correlation Matrix", status = "primary",
              solidHeader = TRUE, width = 6,
              plotlyOutput("bm_correlation_plot", height = "350px")
          ),
          box(title = "ctDNA Variant Allele Frequency (VAF) Over Time",
              status = "success", solidHeader = TRUE, width = 6,
              plotlyOutput("ctdna_vaf_plot", height = "350px")
          )
        )
      )  # end Tab 6

    )  # end tabItems
  )  # end dashboardBody
)  # end UI

## ============================================================================
## SERVER
## ============================================================================
server <- function(input, output, session) {

  ## ---- PK Parameters (reference values) ----
  pk_params_ref <- data.frame(
    Drug = c("Trastuzumab", "Ramucirumab", "Nivolumab", "Pembrolizumab", "T-DXd", "Zolbetuximab"),
    MW_kDa = c(148, 147, 146, 149, 184, 148),
    CL_L_day = c(0.218, 0.282, 0.312, 0.225, 0.198, 0.245),
    V1_L = c(2.91, 3.14, 3.75, 3.20, 3.28, 3.05),
    V2_L = c(2.45, 1.89, 2.94, 2.60, 2.80, 2.40),
    Q_L_day = c(0.742, 0.562, 0.724, 0.680, 0.550, 0.520),
    t_half_days = c(26.9, 21.4, 25.0, 22.0, 17.3, 22.5),
    Mechanism = c("HER2 (ERBB2)", "VEGFR2", "PD-1", "PD-1", "HER2 ADC", "CLDN18.2")
  )

  ## ---- Tab 1: Patient Profile ----
  output$subtype_box <- renderValueBox({
    subtype <- if (input$ebv == "pos") "EBV (9%)"
               else if (input$msi == "msi_h") "MSI-H (22%)"
               else if (input$her2_ihc == "3+" || (input$her2_ihc == "2+" && input$fish == "pos")) "CIN / HER2+"
               else if (input$lauren == "diff") "GS (Genomic Stable)"
               else "CIN (Chromosomal Instab.)"
    valueBox(value = subtype, subtitle = "TCGA Molecular Subtype",
             icon = icon("dna"), color = "purple")
  })

  output$her2_box <- renderValueBox({
    eligible <- if (input$her2_ihc == "3+" || (input$her2_ihc == "2+" && input$fish == "pos"))
      "Eligible (ToGA / KEYNOTE-811)" else if (input$her2_ihc == "1+" || input$her2_ihc == "2+")
      "Eligible for T-DXd" else "Not eligible"
    valueBox(value = paste("HER2:", input$her2_ihc), subtitle = eligible,
             icon = icon("check-circle"), color = "red")
  })

  output$ici_box <- renderValueBox({
    cps_val <- as.numeric(input$cps)
    ici_eligible <- if (cps_val >= 5) "CPS≥5: Nivolumab eligible (CMA 649)"
                    else if (cps_val >= 1) "CPS≥1: Consider pembrolizumab (KN-059)"
                    else "CPS<1: ICI benefit uncertain"
    valueBox(value = paste("CPS =", input$cps), subtitle = ici_eligible,
             icon = icon("shield-alt"), color = if(cps_val >= 5) "green" else "yellow")
  })

  output$cldn_box <- renderValueBox({
    cldn_stat <- if (input$cldn182 == "pos") "Eligible: Zolbetuximab + mFOLFOX6 (SPOTLIGHT)"
                 else "Not eligible for Zolbetuximab"
    valueBox(value = paste("CLDN18.2:", ifelse(input$cldn182 == "pos", "Positive", "Negative")),
             subtitle = cldn_stat,
             icon = icon("atom"), color = if(input$cldn182 == "pos") "blue" else "grey")
  })

  output$subtype_table <- renderDT({
    df <- data.frame(
      Subtype = c("EBV-positive (~9%)", "MSI-High (~22%)", "Genomically Stable (GS, ~20%)",
                  "Chromosomal Instability (CIN, ~50%)"),
      Key_Features = c(
        "EBER+, PIK3CA mut, high PD-L1, CDKN2A silencing, ARID1A mut",
        "MLH1 silencing, ARID1A mut, PIK3CA mut, high TMB, EBV-neg",
        "RHOA/CDH1 mut, diffuse Lauren, CLDN18.2 amp, EMT",
        "TP53 mut, RTK amp (HER2/EGFR/FGFR2/MET), intestinal Lauren"
      ),
      Prevalence_Stage_IV = c("~9%", "~22%", "~20%", "~50%"),
      Preferred_Therapy = c(
        "ICI (high PD-L1), HER2 if CIN overlap",
        "ICI highly effective (pembrolizumab, nivolumab)",
        "Zolbetuximab (CLDN18.2+), FOLFOX ± ICI",
        "Trastuzumab (HER2+), FGFR2i (BGJ398), METi"
      ),
      stringsAsFactors = FALSE
    )
    datatable(df, options = list(dom = 't', paging = FALSE),
              rownames = FALSE, class = "table-bordered table-striped")
  })

  ## ---- Tab 2: Drug PK ----
  pk_sim_data <- eventReactive(input$run_pk, {
    drug <- input$pk_drug
    times <- seq(0, input$n_cycles_pk * 21, by = 0.25)

    params <- list(
      CL = 0.218, V1 = 2.91, V2 = 2.45, Q = 0.742,
      dose = 6 * 65 / 148000 * 1e6,  # nmol dose
      interval = 21, tinf = 0.5, unit = "nmol/L",
      drug_name = "Trastuzumab",
      target = "HER2", Cmin = 0.8, Cmax = 15
    )

    if (grepl("Trastuzumab", drug)) {
      dose_loading <- 8 * 65 / 148000 * 1e3
      dose_maint   <- 6 * 65 / 148000 * 1e3
      dose_df <- data.frame(
        time = c(0, seq(21, (input$n_cycles_pk-1)*21, by=21)),
        amt  = c(dose_loading, rep(dose_maint, input$n_cycles_pk-1))
      )
      conc <- multi_dose_iv(dose_df, 0.218, 2.91, 2.45, 0.742, times)
      trough <- 0.8; peak <- 15; unit <- "nmol/L"
    } else if (grepl("Ramucirumab", drug)) {
      dose_val <- 8 * 65 / 147000 * 1e3
      dose_df <- data.frame(time = seq(0, (input$n_cycles_pk-1)*14, by=14), amt = dose_val)
      conc <- multi_dose_iv(dose_df, 0.282, 3.14, 1.89, 0.562, times)
      trough <- 0.5; peak <- 12; unit <- "nmol/L"
    } else if (grepl("Nivolumab", drug)) {
      dose_val <- 360 / 146000 * 1e3
      dose_df <- data.frame(time = seq(0, (input$n_cycles_pk-1)*21, by=21), amt = dose_val)
      conc <- multi_dose_iv(dose_df, 0.312, 3.75, 2.94, 0.724, times)
      trough <- 0.3; peak <- 8; unit <- "nmol/L"
    } else if (grepl("Pembrolizumab", drug)) {
      dose_val <- 200 / 149000 * 1e3
      dose_df <- data.frame(time = seq(0, (input$n_cycles_pk-1)*21, by=21), amt = dose_val)
      conc <- multi_dose_iv(dose_df, 0.225, 3.20, 2.60, 0.680, times)
      trough <- 0.3; peak <- 8; unit <- "nmol/L"
    } else if (grepl("T-DXd", drug)) {
      dose_val <- 6.4 * 65 / 184000 * 1e3
      dose_df <- data.frame(time = seq(0, (input$n_cycles_pk-1)*21, by=21), amt = dose_val)
      conc <- multi_dose_iv(dose_df, 0.198, 3.28, 2.80, 0.550, times)
      trough <- 0.2; peak <- 10; unit <- "nmol/L"
    } else if (grepl("Zolbetuximab", drug)) {
      dose_loading <- 800 * 1.70 / 148000 * 1e3
      dose_maint   <- 600 * 1.70 / 148000 * 1e3
      dose_df <- data.frame(
        time = c(0, seq(21, (input$n_cycles_pk-1)*21, by=21)),
        amt  = c(dose_loading, rep(dose_maint, input$n_cycles_pk-1))
      )
      conc <- multi_dose_iv(dose_df, 0.245, 3.05, 2.40, 0.520, times)
      trough <- 0.5; peak <- 14; unit <- "nmol/L"
    } else {  # Capecitabine
      cape_times <- seq(0, input$n_cycles_pk * 21, by = 0.25)
      n_cyc <- input$n_cycles_pk
      conc <- rep(0, length(times))
      for (cyc in 0:(n_cyc-1)) {
        for (d in 0:13) {
          dose_time <- cyc * 21 + d
          shifted   <- pmax(times - dose_time, 0)
          conc <- conc + oral_1cmt(3500, 1.92, 0.70, 38.4, 25.2, shifted)
        }
      }
      conc <- conc / 1000  # convert to μg/mL
      trough <- 0.1; peak <- 3; unit <- "μg/mL"
    }

    data.frame(time = times, conc = conc, trough = trough, peak = peak, unit = unit)
  })

  output$pk_plot <- renderPlotly({
    df <- pk_sim_data()
    drug_name <- input$pk_drug

    p <- plot_ly(df, x = ~time, y = ~conc, type = "scatter", mode = "lines",
                 line = list(color = "#1565C0", width = 2), name = drug_name) %>%
      layout(
        title = paste("PK Profile —", drug_name),
        xaxis = list(title = "Time (days)"),
        yaxis = list(title = paste("Concentration (", df$unit[1], ")"))
      )

    if (input$show_therapeutic) {
      p <- p %>%
        add_trace(x = range(df$time), y = c(df$trough[1], df$trough[1]),
                  mode = "lines", line = list(color = "green", dash = "dash"),
                  name = "Target Trough") %>%
        add_trace(x = range(df$time), y = c(df$peak[1], df$peak[1]),
                  mode = "lines", line = list(color = "orange", dash = "dash"),
                  name = "Max Conc Reference")
    }
    p
  })

  output$pk_summary_table <- renderDT({
    df <- pk_sim_data()
    # Compute PK stats per cycle (21-day)
    n_cyc <- min(input$n_cycles_pk, 6)
    stats <- lapply(1:n_cyc, function(cyc) {
      cyc_data <- df[df$time >= (cyc-1)*21 & df$time < cyc*21, ]
      if (nrow(cyc_data) == 0) return(NULL)
      data.frame(
        Cycle = cyc,
        Cmax  = round(max(cyc_data$conc, na.rm=TRUE), 3),
        Ctrough = round(min(cyc_data$conc[cyc_data$time > (cyc-1)*21 + 18], na.rm=TRUE), 3),
        AUC   = round(sum(diff(cyc_data$time) * cyc_data$conc[-1], na.rm=TRUE), 1)
      )
    })
    do.call(rbind, Filter(Negate(is.null), stats)) %>%
      datatable(options = list(dom = 't', paging = FALSE),
                rownames = FALSE, colnames = c("Cycle", "Cmax", "Ctrough", "AUC₀₋₂₁"))
  })

  output$pk_params_table <- renderDT({
    pk_params_ref %>%
      datatable(options = list(dom = 't', paging = FALSE, scrollX = TRUE),
                rownames = FALSE)
  })

  ## ---- Tab 3: Tumor Dynamics ----
  tgi_params_reactive <- reactive({
    list(
      weight = 65, bsa = 1.70,
      # Trastuzumab PK
      V1_tras = 2.91, CL_tras = 0.218, V2_tras = 2.45, Q_tras = 0.742,
      Emax_tras = 0.72, EC50_tras = 0.85, n_tras = 1.4,
      # Ramucirumab PK
      V1_ramu = 3.14, CL_ramu = 0.282, V2_ramu = 1.89, Q_ramu = 0.562,
      Emax_ramu = 0.55, EC50_ramu = 1.12,
      # Nivolumab PK
      V1_nivo = 3.75, CL_nivo = 0.312, V2_nivo = 2.94, Q_nivo = 0.724,
      Emax_nivo = 0.65, EC50_nivo_CD8 = 0.65,
      # T-DXd PK
      V1_tdxd = 3.28, CL_tdxd = 0.198,
      Emax_tdxd = 0.85, EC50_dxd = 0.28,
      # Zolbetuximab PK
      V1_zolbe = 3.05, CL_zolbe = 0.245,
      Emax_zolbe = 0.68, EC50_zolbe_kill = 0.92,
      # Chemo
      ka_cape = 1.92, F_cape = 0.70, CL_cape = 38.4, V_cape = 25.2,
      Emax_chemo = 0.88, EC50_FU = 0.45,
      cape_dose = 3500,
      # TGI
      tv0 = input$tv0,
      lambda0 = input$lambda0, lambda1 = input$lambda1,
      k1_tgi = 0.034 * input$E_drug / 0.70,
      k2_tgi = 0.018,
      CEA0 = 40
    )
  })

  tgi_result <- eventReactive(input$run_tgi, {
    times  <- seq(0, input$sim_duration, by = 0.5)
    params <- tgi_params_reactive()
    res    <- simulate_scenario(input$tgi_scenario, params, times)
    data.frame(time = res$times, tumor_vol = res$tumor_vol,
               cea = res$cea, cd8 = res$cd8,
               scenario = input$tgi_scenario)
  })

  output$tumor_plot <- renderPlotly({
    req(tgi_result())
    df <- tgi_result()
    p <- plot_ly(df, x = ~time, y = ~tumor_vol, type = "scatter", mode = "lines",
                 line = list(color = "#C62828", width = 2), name = input$tgi_scenario) %>%
      add_trace(x = c(0, max(df$time)),
                y = c(df$tumor_vol[1] * 0.7, df$tumor_vol[1] * 0.7),
                mode = "lines", line = list(color = "green", dash = "dash"),
                name = "-30% (PR threshold)") %>%
      add_trace(x = c(0, max(df$time)),
                y = c(df$tumor_vol[1] * 1.2, df$tumor_vol[1] * 1.2),
                mode = "lines", line = list(color = "red", dash = "dash"),
                name = "+20% (PD threshold)") %>%
      layout(
        title = paste("Tumor Volume — Simeoni TGI:", input$tgi_scenario),
        xaxis = list(title = "Time (days)"),
        yaxis = list(title = "Tumor Volume (mm³)")
      )
    p
  })

  output$waterfall_plot <- renderPlotly({
    # Simulated waterfall for n=30 virtual patients
    set.seed(42)
    n_pts <- 30
    scenarios <- c("FLOT", "Tras+FOLFOX", "Ramu+Pac", "Nivo+Chemo", "T-DXd", "Zolbe+mFOLFOX")
    colors_wf  <- c("#E53935", "#8E24AA", "#1E88E5", "#00ACC1", "#43A047", "#FB8C00")
    wf_data <- data.frame(
      patient  = 1:n_pts,
      pct_change = c(sort(rnorm(5, -55, 20)), sort(rnorm(5, -50, 25)),
                     sort(rnorm(5, -30, 30)), sort(rnorm(5, -45, 25)),
                     sort(rnorm(5, -60, 20)), sort(rnorm(5, -42, 22))),
      scenario = rep(scenarios, each = 5)
    ) %>% arrange(pct_change) %>% mutate(pt_rank = row_number())

    plot_ly(wf_data, x = ~pt_rank, y = ~pct_change, type = "bar",
            color = ~scenario, colors = colors_wf) %>%
      add_trace(x = c(0, n_pts+1), y = c(-30, -30), type = "scatter", mode = "lines",
                line = list(color = "green", dash = "dash"), name = "PR threshold (-30%)",
                showlegend = FALSE) %>%
      add_trace(x = c(0, n_pts+1), y = c(20, 20), type = "scatter", mode = "lines",
                line = list(color = "red", dash = "dash"), name = "PD threshold (+20%)",
                showlegend = FALSE) %>%
      layout(title = "Waterfall Plot — Best % Change from Baseline",
             xaxis = list(title = "Patient"), yaxis = list(title = "% Change from Baseline"),
             barmode = "group")
  })

  output$spider_plot <- renderPlotly({
    set.seed(123)
    times_sp <- c(0, 6, 12, 18, 24)
    scenarios_sp <- c("Tras+FOLFOX", "Nivo+Chemo", "T-DXd", "Zolbe+mFOLFOX")
    colors_sp    <- c("#C62828", "#1565C0", "#2E7D32", "#F57F17")

    spider_df <- do.call(rbind, lapply(seq_along(scenarios_sp), function(si) {
      do.call(rbind, lapply(1:5, function(pt) {
        traj <- cumsum(c(0, rnorm(4, mean = c(-15, -10, -5, 5) + si*2, sd = 8)))
        data.frame(time = times_sp, pct_change = traj, patient = pt, scenario = scenarios_sp[si])
      }))
    }))

    p <- plot_ly()
    for (sc in scenarios_sp) {
      sc_data <- spider_df[spider_df$scenario == sc, ]
      for (pt in unique(sc_data$patient)) {
        pt_data <- sc_data[sc_data$patient == pt, ]
        p <- p %>% add_trace(data = pt_data, x = ~time, y = ~pct_change,
                             type = "scatter", mode = "lines",
                             line = list(color = colors_sp[which(scenarios_sp == sc)],
                                         width = 0.8),
                             legendgroup = sc, showlegend = (pt == 1),
                             name = sc)
      }
    }
    p %>% layout(title = "Spider Plot — % Change from Baseline",
                 xaxis = list(title = "Time (weeks)"),
                 yaxis = list(title = "% Change from Baseline", zeroline = TRUE,
                              zerolinecolor = "black"))
  })

  ## ---- Tab 4: Clinical Endpoints ----
  km_data <- eventReactive(input$run_endpoints, {
    set.seed(2024)
    n <- input$n_pts

    # Generate OS from exponential distribution
    os_median_trt  <- input$os_median_trt
    os_median_ctrl <- os_median_trt / input$os_hr
    hr_pfs         <- input$pfs_hr
    pfs_median_trt <- os_median_trt * 0.55
    pfs_median_ctrl<- pfs_median_trt / hr_pfs

    lambda_os_trt  <- log(2) / os_median_trt
    lambda_os_ctrl <- log(2) / os_median_ctrl
    lambda_pfs_trt <- log(2) / pfs_median_trt
    lambda_pfs_ctrl<- log(2) / pfs_median_ctrl

    # OS
    os_trt  <- rexp(n/2, lambda_os_trt)
    os_ctrl <- rexp(n/2, lambda_os_ctrl)
    # Censoring at 36 months
    cens_trt  <- runif(n/2, 18, 36)
    cens_ctrl <- runif(n/2, 18, 36)

    os_data <- data.frame(
      time   = c(pmin(os_trt, cens_trt), pmin(os_ctrl, cens_ctrl)),
      status = c(as.integer(os_trt < cens_trt), as.integer(os_ctrl < cens_ctrl)),
      group  = c(rep("Treatment", n/2), rep("Control", n/2))
    )

    # PFS
    pfs_trt  <- rexp(n/2, lambda_pfs_trt)
    pfs_ctrl <- rexp(n/2, lambda_pfs_ctrl)
    pfs_data <- data.frame(
      time   = c(pmin(pfs_trt, cens_trt), pmin(pfs_ctrl, cens_ctrl)),
      status = c(as.integer(pfs_trt < cens_trt), as.integer(pfs_ctrl < cens_ctrl)),
      group  = c(rep("Treatment", n/2), rep("Control", n/2))
    )

    list(os = os_data, pfs = pfs_data)
  })

  output$km_os_plot <- renderPlot({
    req(km_data())
    df <- km_data()$os
    fit <- survfit(Surv(time, status) ~ group, data = df)
    ggsurvplot(fit, data = df,
               palette  = c("#C62828", "#1565C0"),
               conf.int = TRUE,
               pval     = TRUE,
               risk.table = TRUE,
               xlab = "Time (months)",
               ylab = "Overall Survival Probability",
               title = paste("Kaplan-Meier OS —", input$endpoint_scenario),
               legend.labs  = c("Control", "Treatment"),
               ggtheme = theme_bw(base_size = 12),
               tables.height = 0.25)
  })

  output$km_pfs_plot <- renderPlot({
    req(km_data())
    df <- km_data()$pfs
    fit <- survfit(Surv(time, status) ~ group, data = df)
    ggsurvplot(fit, data = df,
               palette  = c("#2E7D32", "#F57F17"),
               conf.int = TRUE,
               pval     = TRUE,
               xlab = "Time (months)",
               ylab = "PFS Probability",
               title = "Progression-Free Survival",
               legend.labs = c("Control", "Treatment"),
               ggtheme = theme_bw(base_size = 11))
  })

  output$response_bar_plot <- renderPlotly({
    orr_trt  <- input$orr_trt / 100
    orr_ctrl <- input$orr_ctrl / 100
    cr_trt   <- orr_trt * 0.15
    cr_ctrl  <- orr_ctrl * 0.08
    dcr_trt  <- orr_trt + 0.25
    dcr_ctrl <- orr_ctrl + 0.20

    df_resp <- data.frame(
      endpoint = rep(c("CR", "PR", "SD", "DCR", "ORR"), 2),
      value    = c(cr_trt*100, (orr_trt-cr_trt)*100, (dcr_trt-orr_trt)*100, dcr_trt*100, orr_trt*100,
                   cr_ctrl*100, (orr_ctrl-cr_ctrl)*100, (dcr_ctrl-orr_ctrl)*100, dcr_ctrl*100, orr_ctrl*100),
      group    = rep(c("Treatment", "Control"), each = 5)
    )

    plot_ly(df_resp, x = ~endpoint, y = ~value, color = ~group,
            colors = c("#C62828", "#1565C0"), type = "bar", barmode = "group") %>%
      layout(title = "Response Rates — Treatment vs Control",
             xaxis = list(title = "Endpoint"),
             yaxis = list(title = "Rate (%)", range = c(0, 100)))
  })

  ## ---- Tab 5: Scenario Comparison ----
  comparison_data <- eventReactive(input$run_compare, {
    req(length(input$selected_scenarios) >= 1)
    times  <- seq(0, input$sim_duration, by = 1)
    params <- list(
      weight = 65, bsa = 1.70,
      V1_tras = 2.91, CL_tras = 0.218, V2_tras = 2.45, Q_tras = 0.742,
      Emax_tras = 0.72, EC50_tras = 0.85, n_tras = 1.4,
      V1_ramu = 3.14, CL_ramu = 0.282, V2_ramu = 1.89, Q_ramu = 0.562,
      Emax_ramu = 0.55, EC50_ramu = 1.12,
      V1_nivo = 3.75, CL_nivo = 0.312, V2_nivo = 2.94, Q_nivo = 0.724,
      Emax_nivo = 0.65, EC50_nivo_CD8 = 0.65,
      V1_tdxd = 3.28, CL_tdxd = 0.198, Emax_tdxd = 0.85, EC50_dxd = 0.28,
      V1_zolbe = 3.05, CL_zolbe = 0.245, Emax_zolbe = 0.68, EC50_zolbe_kill = 0.92,
      ka_cape = 1.92, F_cape = 0.70, CL_cape = 38.4, V_cape = 25.2,
      Emax_chemo = 0.88, EC50_FU = 0.45, cape_dose = 3500,
      tv0 = input$tv0_comp, lambda0 = 0.0215, lambda1 = 0.385,
      k1_tgi = 0.034, k2_tgi = 0.018, CEA0 = 40
    )

    do.call(rbind, lapply(input$selected_scenarios, function(sc) {
      res <- tryCatch(
        simulate_scenario(sc, params, times),
        error = function(e) {
          list(times = times, tumor_vol = rep(params$tv0, length(times)),
               cea = rep(40, length(times)), cd8 = rep(0.67, length(times)))
        }
      )
      data.frame(time = res$times, tumor_vol = res$tumor_vol,
                 cea = res$cea, cd8 = res$cd8, scenario = sc)
    }))
  })

  output$compare_tv_plot <- renderPlotly({
    req(comparison_data())
    df <- comparison_data()
    colors_comp <- brewer.pal(max(3, length(unique(df$scenario))), "Set1")

    p <- plot_ly()
    for (i in seq_along(unique(df$scenario))) {
      sc   <- unique(df$scenario)[i]
      sc_d <- df[df$scenario == sc, ]
      p <- p %>% add_trace(data = sc_d, x = ~time, y = ~tumor_vol,
                           type = "scatter", mode = "lines",
                           name = sc,
                           line = list(color = colors_comp[i], width = 2.5))
    }
    p %>% layout(title = "Tumor Volume — Treatment Scenario Comparison",
                 xaxis = list(title = "Time (days)"),
                 yaxis = list(title = "Tumor Volume (mm³)", type = "log"))
  })

  output$compare_auc_plot <- renderPlotly({
    req(comparison_data())
    df <- comparison_data()
    # AUC of tumor suppression
    auc_data <- df %>%
      group_by(scenario) %>%
      summarize(
        AUC_TV    = sum(diff(time) * tumor_vol[-1], na.rm = TRUE) / 1e6,
        TV_nadir  = round(min(tumor_vol), 0),
        TV_final  = round(last(tumor_vol), 0),
        .groups   = "drop"
      ) %>%
      mutate(pct_suppression = round((1 - TV_final / TV_nadir[1]) * 100, 1))

    plot_ly(auc_data, x = ~scenario, y = ~AUC_TV, type = "bar",
            marker = list(color = brewer.pal(max(3, nrow(auc_data)), "Set1"))) %>%
      layout(title = "Tumor Burden AUC (×10⁶ mm³·day)",
             xaxis = list(title = ""), yaxis = list(title = "AUC (×10⁶ mm³·day)"))
  })

  output$scenario_summary_table <- renderDT({
    req(comparison_data())
    df <- comparison_data()
    summary_df <- df %>%
      group_by(scenario) %>%
      summarize(
        TV_baseline   = round(first(tumor_vol), 0),
        TV_min        = round(min(tumor_vol), 0),
        TV_final      = round(last(tumor_vol), 0),
        Pct_change    = round((TV_min - TV_baseline) / TV_baseline * 100, 1),
        CEA_final     = round(last(cea), 1),
        CD8_final     = round(last(cd8), 3),
        .groups       = "drop"
      ) %>%
      rename(`TV baseline` = TV_baseline, `TV nadir` = TV_min,
             `TV final` = TV_final, `Best % change` = Pct_change,
             `CEA final` = CEA_final, `CD8 final` = CD8_final)
    datatable(summary_df, options = list(dom = 't', paging = FALSE, scrollX = TRUE),
              rownames = FALSE)
  })

  ## ---- Tab 6: Biomarkers ----
  bm_sim_data <- eventReactive(input$run_bm, {
    times  <- seq(0, input$sim_duration, by = 1)
    params <- list(
      weight = 65, bsa = 1.70,
      V1_tras = 2.91, CL_tras = 0.218, V2_tras = 2.45, Q_tras = 0.742,
      Emax_tras = 0.72, EC50_tras = 0.85, n_tras = 1.4,
      V1_ramu = 3.14, CL_ramu = 0.282, V2_ramu = 1.89, Q_ramu = 0.562,
      Emax_ramu = 0.55, EC50_ramu = 1.12,
      V1_nivo = 3.75, CL_nivo = 0.312, V2_nivo = 2.94, Q_nivo = 0.724,
      Emax_nivo = 0.65, EC50_nivo_CD8 = 0.65,
      V1_tdxd = 3.28, CL_tdxd = 0.198, Emax_tdxd = 0.85, EC50_dxd = 0.28,
      V1_zolbe = 3.05, CL_zolbe = 0.245, Emax_zolbe = 0.68, EC50_zolbe_kill = 0.92,
      ka_cape = 1.92, F_cape = 0.70, CL_cape = 38.4, V_cape = 25.2,
      Emax_chemo = 0.88, EC50_FU = 0.45, cape_dose = 3500,
      tv0 = 1200, lambda0 = 0.0215, lambda1 = 0.385,
      k1_tgi = 0.034, k2_tgi = 0.018, CEA0 = input$cea_base_bm
    )
    res <- simulate_scenario(input$bm_scenario, params, times)

    # Simulate additional biomarkers
    tumor_ratio <- res$tumor_vol / 1200
    ca199 <- pmax(5, input$ca199_base * tumor_ratio * 0.85 + rnorm(length(times), 0, 10))
    ctdna  <- pmax(10, input$ctdna_base * tumor_ratio^1.5 * exp(-0.003 * times) +
                     rnorm(length(times), 0, 50))
    vaf    <- pmax(0, 15 * tumor_ratio * exp(-0.005 * times))

    data.frame(
      time      = times,
      tumor_vol = res$tumor_vol,
      cea       = res$cea,
      ca199     = ca199,
      ctdna     = ctdna,
      vaf       = vaf,
      cd8       = res$cd8,
      scenario  = input$bm_scenario
    )
  })

  output$bm_serum_plot <- renderPlotly({
    req(bm_sim_data())
    df <- bm_sim_data()

    p <- plot_ly() %>%
      add_trace(data = df, x = ~time, y = ~cea, type = "scatter", mode = "lines",
                name = "CEA (ng/mL)", line = list(color = "#C62828", width = 2)) %>%
      add_trace(data = df, x = ~time, y = ~ca199/10, type = "scatter", mode = "lines",
                name = "CA19-9 (×10 U/mL)", line = list(color = "#1565C0", width = 2)) %>%
      add_trace(data = df, x = ~time, y = ~ctdna/100, type = "scatter", mode = "lines",
                name = "ctDNA (×100 copies/mL)", line = list(color = "#2E7D32", width = 2)) %>%
      layout(title = paste("Serum Biomarkers —", input$bm_scenario),
             xaxis = list(title = "Time (days)"),
             yaxis = list(title = "Biomarker Level (scaled)"),
             hovermode = "x unified")
    p
  })

  output$her2_ihc_plot <- renderPlotly({
    # HER2 IHC distribution in gastric cancer (TCGA + clinical data)
    df_her2 <- data.frame(
      score = c("0", "1+", "2+/FISH+", "2+/FISH-", "3+"),
      pct   = c(48, 24, 9, 8, 11),
      color = c("#BDBDBD", "#90CAF9", "#FFCA28", "#FF7043", "#C62828")
    )
    plot_ly(df_her2, x = ~score, y = ~pct, type = "bar",
            marker = list(color = df_her2$color)) %>%
      layout(title = "HER2 IHC Score Distribution\n(Gastric/GEJ Adenocarcinoma)",
             xaxis = list(title = "HER2 IHC Score"),
             yaxis = list(title = "Percentage (%)"))
  })

  output$cps_plot <- renderPlotly({
    # CPS distribution and ICI eligibility
    cps_vals <- c(0, 1, 5, 10, 50, 100)
    labels   <- c("CPS=0", "CPS 1-4\n(pembrolizumab\nKN-059)", "CPS 5-9\n(nivolumab\nCM-649)",
                  "CPS 10-49", "CPS 50+\n(high responder)")
    pcts     <- c(18, 24, 22, 20, 16)

    plot_ly(x = labels[-length(labels)], y = pcts, type = "bar",
            marker = list(color = c("#BDBDBD", "#FFEB3B", "#FF9800", "#F44336", "#B71C1C"))) %>%
      layout(title = "PD-L1 CPS Distribution\n(Gastric Cancer)",
             xaxis = list(title = "CPS Category"),
             yaxis = list(title = "% Patients"))
  })

  output$molecular_subtype_plot <- renderPlotly({
    # TCGA molecular subtypes
    df_sub <- data.frame(
      subtype = c("EBV+", "MSI-High", "Genomic Stable", "CIN"),
      pct     = c(9, 22, 20, 50),
      color   = c("#7B1FA2", "#1565C0", "#2E7D32", "#C62828")
    )
    plot_ly(df_sub, labels = ~subtype, values = ~pct, type = "pie",
            marker = list(colors = df_sub$color),
            textinfo = "label+percent") %>%
      layout(title = "TCGA Molecular Subtypes\n(Gastric Cancer)")
  })

  output$bm_correlation_plot <- renderPlotly({
    req(bm_sim_data())
    df <- bm_sim_data()
    # Correlation matrix
    cor_df <- cor(df[, c("tumor_vol", "cea", "ca199", "ctdna", "cd8", "vaf")],
                  use = "complete.obs")
    plot_ly(z = cor_df, x = colnames(cor_df), y = rownames(cor_df),
            type = "heatmap", colorscale = "RdBu", reversescale = TRUE,
            zmin = -1, zmax = 1) %>%
      layout(title = "Biomarker Correlation Matrix")
  })

  output$ctdna_vaf_plot <- renderPlotly({
    req(bm_sim_data())
    df <- bm_sim_data()
    # Simulate multi-clone VAF dynamics
    set.seed(99)
    n_clones <- 4
    clone_names <- c("TP53 R175H", "KRAS G12D", "PIK3CA H1047R", "HER2 amp")
    clone_freqs <- matrix(0, nrow = length(df$time), ncol = n_clones)
    for (cl in 1:n_clones) {
      initial_vaf <- runif(1, 5, 25)
      decay_rate  <- runif(1, 0.003, 0.015)
      clone_freqs[, cl] <- pmax(0, initial_vaf * df$tumor_vol / max(df$tumor_vol) * exp(-decay_rate * df$time) + rnorm(nrow(df), 0, 0.3))
    }

    p <- plot_ly()
    colors_clone <- c("#C62828", "#1565C0", "#2E7D32", "#F57F17")
    for (cl in 1:n_clones) {
      p <- p %>% add_trace(x = df$time, y = clone_freqs[, cl],
                           type = "scatter", mode = "lines",
                           name = clone_names[cl],
                           line = list(color = colors_clone[cl], width = 2))
    }
    p %>% layout(title = "ctDNA VAF — Clonal Dynamics Over Time",
                 xaxis = list(title = "Time (days)"),
                 yaxis = list(title = "Variant Allele Frequency (%)"))
  })

}  # end server

## ============================================================================
## RUN APP
## ============================================================================
shinyApp(ui = ui, server = server)
