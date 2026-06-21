# =============================================================================
# Autoimmune Pancreatitis (AIP) QSP Shiny Application
# Quantitative Systems Pharmacology Interactive Dashboard
# =============================================================================
# Author: QSP Disease Model Library (CCR)
# Date: 2026-06-19
# Version: 1.0
# =============================================================================

library(shiny)
library(ggplot2)
library(dplyr)
library(tidyr)
library(shinythemes)
library(gridExtra)
library(scales)

# =============================================================================
# HELPER FUNCTIONS & MODEL EQUATIONS
# =============================================================================

# PK model: one-compartment oral / IV
pk_prednisolone <- function(time, dose, weight, ka = 1.5, Vd_L_kg = 0.7,
                             CL_L_h_kg = 0.1) {
  Vd <- Vd_L_kg * weight
  CL <- CL_L_h_kg * weight
  ke <- CL / Vd
  C <- (dose / Vd) * (ka / (ka - ke)) * (exp(-ke * time) - exp(-ka * time))
  pmax(C, 0)
}

pk_rituximab <- function(time, dose, Vd = 3.1, CL = 0.014, t_half_alpha = 20,
                          t_half_beta = 206) {
  # Two-compartment IV bolus approximation (mg/L, dose in mg)
  A <- dose / Vd * 0.6
  B <- dose / Vd * 0.4
  alpha <- log(2) / t_half_alpha
  beta  <- log(2) / t_half_beta
  C <- A * exp(-alpha * time) + B * exp(-beta * time)
  pmax(C, 0)
}

pk_azathioprine <- function(time, dose, weight, ka = 0.8, Vd_kg = 0.6,
                             CL_kg = 0.04) {
  Vd <- Vd_kg * weight
  CL <- CL_kg * weight
  ke <- CL / Vd
  C <- (dose / Vd) * (ka / (ka - ke)) * (exp(-ke * time) - exp(-ka * time))
  pmax(C, 0)
}

pk_mmf <- function(time, dose, weight, ka = 1.2, Vd_kg = 3.0, CL_kg = 0.16) {
  Vd <- Vd_kg * weight
  CL <- CL_kg * weight
  ke <- CL / Vd
  C <- (dose / Vd) * (ka / (ka - ke)) * (exp(-ke * time) - exp(-ka * time))
  pmax(C, 0)
}

# IgG4 dynamics model
# dIgG4/dt = production(B cells) - clearance - drug suppression
igg4_dynamics <- function(time_weeks, baseline_igg4, drug,
                           dose, weight, aip_type,
                           extra_pancreatic, ka = 0.8, Vd_kg = 0.7,
                           CL_kg = 0.1) {
  n  <- length(time_weeks)
  ig <- numeric(n)
  ig[1] <- baseline_igg4

  # Production and clearance rates
  k_prod <- 0.05   # baseline IgG4 production (mg/dL per week)
  k_clear <- 0.02  # natural clearance
  k_extra <- if (extra_pancreatic) 1.3 else 1.0  # modifier for extra-pancreatic

  for (i in 2:n) {
    t_h  <- time_weeks[i] * 168  # convert weeks to hours
    dt   <- (time_weeks[i] - time_weeks[i - 1]) * 7  # days

    # Drug effect on IgG4 production
    eff <- switch(drug,
      "Prednisolone" = {
        Cp <- pk_prednisolone(t_h %% 24, dose, weight)
        EC50 <- 0.15; Emax <- 0.75
        Emax * Cp / (EC50 + Cp)
      },
      "Rituximab" = {
        # Rituximab: pulses at week 0 and 24
        t_since <- min(t_h, t_h - 24 * 168)
        Cp <- pk_rituximab(max(t_h - 1, 0), dose)
        EC50 <- 0.05; Emax <- 0.92
        Emax * Cp / (EC50 + Cp)
      },
      "Azathioprine" = {
        Cp <- pk_azathioprine(t_h %% 24, dose, weight)
        EC50 <- 0.8; Emax <- 0.55
        Emax * Cp / (EC50 + Cp)
      },
      "MMF" = {
        Cp <- pk_mmf(t_h %% 24, dose, weight)
        EC50 <- 1.2; Emax <- 0.60
        Emax * Cp / (EC50 + Cp)
      },
      "Combination" = {
        Cp_p <- pk_prednisolone(t_h %% 24, dose * 0.5, weight)
        Cp_a <- pk_azathioprine(t_h %% 24, dose * 0.5, weight)
        eff_p <- 0.75 * Cp_p / (0.15 + Cp_p)
        eff_a <- 0.55 * Cp_a / (0.80 + Cp_a)
        min(eff_p + eff_a * (1 - eff_p), 0.95)
      },
      0.0
    )

    ig[i] <- ig[i - 1] + dt * (
      k_prod * ig[i - 1] * k_extra * (1 - eff) -
        k_clear * ig[i - 1]
    )
    ig[i] <- max(ig[i], 5)  # floor at 5 mg/dL
  }
  ig
}

# Exocrine function dynamics (FE-1 proxy, μg/g)
exocrine_dynamics <- function(time_weeks, baseline_fe1, drug, dose, weight,
                               igg4_vec) {
  n  <- length(time_weeks)
  fe <- numeric(n)
  fe[1] <- baseline_fe1
  normal_fe1 <- 200  # threshold

  for (i in 2:n) {
    dt <- (time_weeks[i] - time_weeks[i - 1]) * 7
    igg4_effect <- (igg4_vec[i] - 135) / 1000  # IgG4 suppression on pancreas
    igg4_effect <- max(igg4_effect, 0)

    k_recovery <- 0.008  # exocrine recovery rate
    k_damage   <- 0.003 * igg4_effect

    fe[i] <- fe[i - 1] + dt * (
      k_recovery * (normal_fe1 - fe[i - 1]) - k_damage * fe[i - 1]
    )
    fe[i] <- max(min(fe[i], 500), 10)
  }
  fe
}

# HbA1c / endocrine dynamics
hba1c_dynamics <- function(time_weeks, baseline_hba1c, drug, dose, igg4_vec) {
  n    <- length(time_weeks)
  hba1c <- numeric(n)
  hba1c[1] <- baseline_hba1c
  target <- 5.7  # near-normal

  for (i in 2:n) {
    dt <- (time_weeks[i] - time_weeks[i - 1]) * 7
    igg4_drive <- pmax(igg4_vec[i] - 135, 0) / 2000

    steroid_effect <- if (drug %in% c("Prednisolone", "Combination")) {
      0.003 * (dose / 40)
    } else { 0 }

    k_improve <- 0.005
    hba1c[i] <- hba1c[i - 1] + dt * (
      -k_improve * (hba1c[i - 1] - target) +
        steroid_effect +
        igg4_drive
    )
    hba1c[i] <- max(hba1c[i], 4.5)
  }
  hba1c
}

# B-cell depletion for rituximab
bcell_depletion <- function(time_weeks, dose) {
  # B cells return to normal in ~6-9 months
  baseline <- 100  # %
  depletion <- numeric(length(time_weeks))
  for (i in seq_along(time_weeks)) {
    t <- time_weeks[i]
    if (t <= 1) {
      depletion[i] <- baseline * exp(-3.5 * t)
    } else {
      k_recovery <- 0.015
      depletion[i] <- 5 + (baseline - 5) * (1 - exp(-k_recovery * (t - 1)))
    }
  }
  pmax(depletion, 0)
}

# Relapse risk scoring
relapse_risk <- function(igg4_6m, extra_pancreatic, steroid_dose_at_6m) {
  score <- 0
  if (igg4_6m > 135) score <- score + 2
  else if (igg4_6m > 67) score <- score + 1
  if (extra_pancreatic) score <- score + 1
  if (steroid_dose_at_6m > 5) score <- score + 1
  if (score <= 1) "Low" else if (score <= 2) "Medium" else "High"
}

# =============================================================================
# UI DEFINITION
# =============================================================================

ui <- fluidPage(
  theme = shinytheme("flatly"),

  tags$head(
    tags$style(HTML("
      .main-title {
        background: linear-gradient(135deg, #2c3e50, #3498db);
        color: white; padding: 20px 30px; margin-bottom: 20px;
        border-radius: 8px;
      }
      .main-title h2 { margin: 0; font-size: 24px; }
      .main-title p  { margin: 5px 0 0; opacity: 0.85; font-size: 13px; }
      .patient-card {
        background: #f8f9fa; border: 1px solid #dee2e6;
        border-left: 5px solid #3498db;
        border-radius: 6px; padding: 15px; margin-bottom: 10px;
      }
      .biomarker-table th { background-color: #3498db; color: white; }
      .risk-low    { color: #27ae60; font-weight: bold; }
      .risk-medium { color: #f39c12; font-weight: bold; }
      .risk-high   { color: #e74c3c; font-weight: bold; }
      .tab-label   { font-size: 13px; }
      .section-header { color: #2c3e50; border-bottom: 2px solid #3498db;
                        padding-bottom: 6px; margin-top: 15px; }
    "))
  ),

  div(class = "main-title",
    tags$h2("Autoimmune Pancreatitis (AIP) — QSP Dashboard"),
    tags$p("Quantitative Systems Pharmacology Model | IgG4-Related Disease · Pancreatic Function · Immunotherapy")
  ),

  tabsetPanel(
    id = "main_tabs",

    # =========================================================================
    # TAB 1: PATIENT PROFILE
    # =========================================================================
    tabPanel(
      title = span(class = "tab-label", icon("user"), " Patient Profile"),
      value = "tab_patient",
      br(),
      sidebarLayout(
        sidebarPanel(
          width = 4,
          h4(class = "section-header", "Demographics"),
          sliderInput("age",   "Age (years)",      min = 20,  max = 80,  value = 58, step = 1),
          sliderInput("weight","Weight (kg)",       min = 40,  max = 120, value = 68, step = 1),
          selectInput("sex",   "Sex",              choices = c("Male", "Female"), selected = "Male"),
          sliderInput("disease_duration", "Disease Duration (months)", min = 1, max = 60, value = 6),

          h4(class = "section-header", "Baseline Biomarkers"),
          numericInput("baseline_igg4", "Serum IgG4 (mg/dL)", value = 420, min = 10, max = 3000),
          numericInput("baseline_fe1",  "FE-1 (μg/g stool)",  value = 120, min = 5,  max = 500),
          numericInput("baseline_hba1c","HbA1c (%)",           value = 6.8, min = 4,  max = 14, step = 0.1),

          h4(class = "section-header", "Disease Characteristics"),
          selectInput("aip_type", "AIP Type",
                      choices = c("Type 1 (IgG4-related)" = "Type1",
                                  "Type 2 (IDCP)" = "Type2")),
          checkboxGroupInput("extra_pancreatic",
                             "Extra-Pancreatic Manifestations",
                             choices = c("Biliary stricture"     = "biliary",
                                         "Salivary gland"        = "salivary",
                                         "Retroperitoneal fibrosis" = "retro",
                                         "Renal involvement"     = "kidney"),
                             selected = c("biliary"))
        ),

        mainPanel(
          width = 8,
          fluidRow(
            column(12,
              h4(class = "section-header", "Patient Summary"),
              uiOutput("patient_card")
            )
          ),
          br(),
          fluidRow(
            column(6,
              h4(class = "section-header", "Baseline Biomarkers"),
              tableOutput("biomarker_table")
            ),
            column(6,
              h4(class = "section-header", "Disease Classification"),
              tableOutput("disease_class_table"),
              br(),
              uiOutput("aip_type_info")
            )
          )
        )
      )
    ),

    # =========================================================================
    # TAB 2: TREATMENT & PK
    # =========================================================================
    tabPanel(
      title = span(class = "tab-label", icon("pills"), " Treatment & PK"),
      value = "tab_pk",
      br(),
      sidebarLayout(
        sidebarPanel(
          width = 4,
          h4(class = "section-header", "Drug Selection"),
          selectInput("drug_pk", "Treatment",
                      choices = c("Prednisolone", "Rituximab",
                                  "Azathioprine", "MMF", "Combination")),

          h4(class = "section-header", "Dose Parameters"),
          conditionalPanel(
            condition = "input.drug_pk == 'Prednisolone' || input.drug_pk == 'Combination'",
            sliderInput("pred_dose",  "Prednisolone dose (mg/day)", min = 5, max = 60, value = 40, step = 5)
          ),
          conditionalPanel(
            condition = "input.drug_pk == 'Rituximab'",
            sliderInput("rtx_dose",  "Rituximab dose (mg infusion)", min = 375, max = 1000, value = 1000, step = 125)
          ),
          conditionalPanel(
            condition = "input.drug_pk == 'Azathioprine' || input.drug_pk == 'Combination'",
            sliderInput("aza_dose",  "Azathioprine dose (mg/day)", min = 25, max = 200, value = 100, step = 25)
          ),
          conditionalPanel(
            condition = "input.drug_pk == 'MMF'",
            sliderInput("mmf_dose",  "MMF dose (mg BID)", min = 250, max = 1500, value = 750, step = 250)
          ),

          h4(class = "section-header", "PK Time Window"),
          sliderInput("pk_time", "Observation window (hours)", min = 4, max = 96, value = 48)
        ),

        mainPanel(
          width = 8,
          h4(class = "section-header", "Concentration–Time Profile"),
          plotOutput("pk_plot", height = "320px"),
          br(),
          h4(class = "section-header", "PK Parameter Summary"),
          tableOutput("pk_params_table")
        )
      )
    ),

    # =========================================================================
    # TAB 3: IgG4 & IMMUNE RESPONSE
    # =========================================================================
    tabPanel(
      title = span(class = "tab-label", icon("vial"), " IgG4 & Immune Response"),
      value = "tab_igg4",
      br(),
      sidebarLayout(
        sidebarPanel(
          width = 4,
          h4(class = "section-header", "Treatment Settings"),
          selectInput("drug_igg4", "Treatment",
                      choices = c("Prednisolone", "Rituximab",
                                  "Azathioprine", "MMF", "Combination")),
          sliderInput("igg4_dose", "Dose (reference units)", min = 10, max = 60, value = 40),
          sliderInput("igg4_horizon", "Follow-up (weeks)", min = 12, max = 104, value = 52),

          h4(class = "section-header", "Response Thresholds"),
          numericInput("igg4_uln",     "IgG4 ULN (mg/dL)",        value = 135),
          numericInput("igg4_remission","Remission target (mg/dL)", value = 135),

          h4(class = "section-header", "Disease Modifier"),
          checkboxInput("ep_check", "Extra-pancreatic involvement", value = TRUE)
        ),

        mainPanel(
          width = 8,
          h4(class = "section-header", "Serum IgG4 Over Time"),
          plotOutput("igg4_plot", height = "300px"),
          br(),
          conditionalPanel(
            condition = "input.drug_igg4 == 'Rituximab'",
            h4(class = "section-header", "B-Cell Depletion (Rituximab)"),
            plotOutput("bcell_plot", height = "230px")
          ),
          br(),
          h4(class = "section-header", "Treatment Response Classification"),
          uiOutput("response_classification")
        )
      )
    ),

    # =========================================================================
    # TAB 4: PANCREATIC FUNCTION
    # =========================================================================
    tabPanel(
      title = span(class = "tab-label", icon("heartbeat"), " Pancreatic Function"),
      value = "tab_pancreas",
      br(),
      sidebarLayout(
        sidebarPanel(
          width = 4,
          h4(class = "section-header", "Treatment"),
          selectInput("drug_panc", "Treatment",
                      choices = c("Prednisolone", "Rituximab",
                                  "Azathioprine", "MMF", "Combination")),
          sliderInput("panc_dose",    "Dose (reference units)", min = 10, max = 60, value = 40),
          sliderInput("panc_horizon", "Follow-up (weeks)",      min = 12, max = 104, value = 52),

          h4(class = "section-header", "Pancreatic Parameters"),
          numericInput("pd_baseline_fe1",  "Baseline FE-1 (μg/g)",  value = 120, min = 5,  max = 500),
          numericInput("pd_baseline_hba1c","Baseline HbA1c (%)",     value = 6.8, min = 4,  max = 14, step = 0.1),
          numericInput("pd_duct_index",    "Duct Diameter Index (0–1)", value = 0.6, min = 0, max = 1, step = 0.05),

          h4(class = "section-header", "Thresholds"),
          numericInput("fe1_threshold",  "FE-1 normal threshold (μg/g)", value = 200),
          numericInput("hba1c_threshold","HbA1c target (%)",             value = 7.0, step = 0.1)
        ),

        mainPanel(
          width = 8,
          h4(class = "section-header", "Exocrine Function (FE-1)"),
          plotOutput("fe1_plot",   height = "230px"),
          br(),
          h4(class = "section-header", "Endocrine Function (HbA1c)"),
          plotOutput("hba1c_plot", height = "230px"),
          br(),
          h4(class = "section-header", "Combined Organ Function Score"),
          plotOutput("organ_score_plot", height = "180px")
        )
      )
    ),

    # =========================================================================
    # TAB 5: SCENARIO COMPARISON
    # =========================================================================
    tabPanel(
      title = span(class = "tab-label", icon("chart-bar"), " Scenario Comparison"),
      value = "tab_scenario",
      br(),
      sidebarLayout(
        sidebarPanel(
          width = 4,
          h4(class = "section-header", "Scenario 1"),
          selectInput("sc1_drug", "Drug",  choices = c("Prednisolone","Rituximab","Azathioprine","MMF","Combination"), selected = "Prednisolone"),
          sliderInput("sc1_dose", "Dose",  min = 10, max = 60, value = 40),

          h4(class = "section-header", "Scenario 2"),
          selectInput("sc2_drug", "Drug",  choices = c("Prednisolone","Rituximab","Azathioprine","MMF","Combination"), selected = "Rituximab"),
          sliderInput("sc2_dose", "Dose",  min = 10, max = 60, value = 40),

          h4(class = "section-header", "Scenario 3"),
          selectInput("sc3_drug", "Drug",  choices = c("Prednisolone","Rituximab","Azathioprine","MMF","Combination"), selected = "Azathioprine"),
          sliderInput("sc3_dose", "Dose",  min = 10, max = 60, value = 40),

          h4(class = "section-header", "Scenario 4"),
          selectInput("sc4_drug", "Drug",  choices = c("Prednisolone","Rituximab","Azathioprine","MMF","Combination"), selected = "Combination"),
          sliderInput("sc4_dose", "Dose",  min = 10, max = 60, value = 40),

          br(),
          actionButton("run_scenarios", "Run Comparison", class = "btn-primary btn-block")
        ),

        mainPanel(
          width = 8,
          h4(class = "section-header", "IgG4 Trajectories — All Scenarios"),
          plotOutput("scenario_igg4_plot",  height = "280px"),
          br(),
          h4(class = "section-header", "Summary Table"),
          tableOutput("scenario_summary_table"),
          br(),
          h4(class = "section-header", "Forest Plot — Relative Treatment Effects (vs. Prednisolone)"),
          plotOutput("forest_plot", height = "220px")
        )
      )
    ),

    # =========================================================================
    # TAB 6: BIOMARKERS & MONITORING
    # =========================================================================
    tabPanel(
      title = span(class = "tab-label", icon("stethoscope"), " Biomarkers & Monitoring"),
      value = "tab_biomarker",
      br(),
      sidebarLayout(
        sidebarPanel(
          width = 4,
          h4(class = "section-header", "Patient Values at 6 Months"),
          numericInput("bm_igg4_6m",  "IgG4 at 6 months (mg/dL)", value = 180, min = 5, max = 3000),
          numericInput("bm_ca199_6m", "CA 19-9 (U/mL)",            value = 25,  min = 0, max = 500),
          numericInput("bm_fe1_6m",   "FE-1 (μg/g)",               value = 175, min = 5, max = 500),
          numericInput("bm_hba1c_6m", "HbA1c (%)",                 value = 6.5, min = 4, max = 14, step = 0.1),
          checkboxInput("bm_extra",   "Extra-pancreatic involvement", value = TRUE),
          numericInput("bm_steroid_dose", "Steroid dose at assessment (mg/day)", value = 5, min = 0, max = 60),

          h4(class = "section-header", "Monitoring Interval"),
          selectInput("monitor_interval", "Monitoring frequency",
                      choices = c("Every month (active)" = 4,
                                  "Every 3 months"       = 13,
                                  "Every 6 months"       = 26),
                      selected = 13)
        ),

        mainPanel(
          width = 8,
          h4(class = "section-header", "Monitoring Schedule & Biomarker Status"),
          plotOutput("monitoring_plot", height = "300px"),
          br(),
          h4(class = "section-header", "Biomarker Threshold Reference"),
          tableOutput("biomarker_threshold_table"),
          br(),
          h4(class = "section-header", "Relapse Risk Classification"),
          uiOutput("relapse_risk_ui")
        )
      )
    )
  )
)

# =============================================================================
# SERVER DEFINITION
# =============================================================================

server <- function(input, output, session) {

  # ---------------------------------------------------------------------------
  # Reactive: base patient inputs
  # ---------------------------------------------------------------------------
  patient_data <- reactive({
    list(
      age              = input$age,
      weight           = input$weight,
      sex              = input$sex,
      duration_months  = input$disease_duration,
      igg4             = input$baseline_igg4,
      fe1              = input$baseline_fe1,
      hba1c            = input$baseline_hba1c,
      aip_type         = input$aip_type,
      extra_pancreatic = input$extra_pancreatic
    )
  })

  # ---------------------------------------------------------------------------
  # TAB 1 OUTPUTS
  # ---------------------------------------------------------------------------
  output$patient_card <- renderUI({
    p   <- patient_data()
    ep  <- if (length(p$extra_pancreatic) > 0) paste(p$extra_pancreatic, collapse = ", ") else "None"
    igg_status <- if (p$igg4 > 135) "ELEVATED" else "Normal"
    igg_col    <- if (p$igg4 > 135) "#e74c3c" else "#27ae60"
    div(class = "patient-card",
      fluidRow(
        column(4, strong("Age:"), tags$span(paste(p$age, "years")), br(),
               strong("Weight:"), tags$span(paste(p$weight, "kg")), br(),
               strong("Sex:"), tags$span(p$sex)),
        column(4, strong("Duration:"), tags$span(paste(p$duration_months, "months")), br(),
               strong("AIP Type:"), tags$span(gsub("Type", "Type ", p$aip_type)), br(),
               strong("Extra-pancreatic:"), tags$span(ep)),
        column(4, strong("IgG4 status:"),
               tags$span(style = paste("color:", igg_col, "; font-weight:bold"), igg_status), br(),
               strong("FE-1:"), tags$span(paste(p$fe1, "μg/g")), br(),
               strong("HbA1c:"), tags$span(paste(p$hba1c, "%")))
      )
    )
  })

  output$biomarker_table <- renderTable({
    p <- patient_data()
    data.frame(
      Biomarker  = c("Serum IgG4 (mg/dL)", "FE-1 (μg/g)", "HbA1c (%)"),
      Value      = c(p$igg4, p$fe1, p$hba1c),
      Normal     = c("< 135", "> 200", "< 5.7"),
      Status     = c(ifelse(p$igg4 > 135, "High", "Normal"),
                     ifelse(p$fe1  < 200, "Low",  "Normal"),
                     ifelse(p$hba1c > 6.5,"High", "Normal"))
    )
  }, striped = TRUE, hover = TRUE, bordered = TRUE)

  output$disease_class_table <- renderTable({
    p <- patient_data()
    data.frame(
      Feature = c("AIP Type", "Extra-pancreatic", "Serum IgG4", "Exocrine",   "Glycemic"),
      Status  = c(gsub("Type", "Type ", p$aip_type),
                  ifelse(length(p$extra_pancreatic) > 0, "Present", "Absent"),
                  ifelse(p$igg4 > 135, "Elevated", "Normal"),
                  ifelse(p$fe1 < 200,  "Insufficient", "Sufficient"),
                  ifelse(p$hba1c > 6.5,"Impaired", "Normal"))
    )
  }, striped = TRUE, bordered = TRUE)

  output$aip_type_info <- renderUI({
    if (input$aip_type == "Type1") {
      div(class = "patient-card",
        strong("Type 1 AIP (IgG4-related)"), br(),
        tags$ul(
          tags$li("IgG4+ plasma cell infiltration"),
          tags$li("Storiform fibrosis, obliterative phlebitis"),
          tags$li("Responds well to corticosteroids"),
          tags$li("Higher relapse rate (~30–40%)")
        )
      )
    } else {
      div(class = "patient-card",
        strong("Type 2 AIP (IDCP)"), br(),
        tags$ul(
          tags$li("Granulocytic epithelial lesions (GEL)"),
          tags$li("IBD association in ~30%"),
          tags$li("Excellent steroid response"),
          tags$li("Lower relapse rate (~10–15%)")
        )
      )
    }
  })

  # ---------------------------------------------------------------------------
  # TAB 2 OUTPUTS — PK
  # ---------------------------------------------------------------------------
  pk_time_vec <- reactive({
    seq(0, input$pk_time, by = 0.5)
  })

  pk_conc_df <- reactive({
    t   <- pk_time_vec()
    drug <- input$drug_pk
    wt   <- input$weight

    C <- switch(drug,
      "Prednisolone" = pk_prednisolone(t, input$pred_dose, wt),
      "Rituximab"    = pk_rituximab(t, input$rtx_dose),
      "Azathioprine" = pk_azathioprine(t, input$aza_dose, wt),
      "MMF"          = pk_mmf(t, input$mmf_dose, wt),
      "Combination"  = {
        Cp <- pk_prednisolone(t, input$pred_dose * 0.5, wt)
        Ca <- pk_azathioprine(t, input$aza_dose * 0.5, wt)
        Cp + Ca  # simplified combined profile
      }
    )
    data.frame(time = t, conc = C, drug = drug)
  })

  output$pk_plot <- renderPlot({
    df   <- pk_conc_df()
    drug <- input$drug_pk
    unit_label <- switch(drug,
      "Rituximab"    = "Concentration (mg/L)",
      "Combination"  = "Concentration (normalized, mg/L equiv.)",
      "Concentration (mg/L)")

    ggplot(df, aes(x = time, y = conc)) +
      geom_line(color = "#2980b9", linewidth = 1.4) +
      geom_area(fill = "#2980b9", alpha = 0.12) +
      labs(title = paste(drug, "— PK Profile"),
           x = "Time (hours)", y = unit_label) +
      scale_x_continuous(breaks = pretty_breaks(8)) +
      scale_y_continuous(labels = number_format(accuracy = 0.01)) +
      theme_minimal(base_size = 13) +
      theme(plot.title = element_text(face = "bold", color = "#2c3e50"))
  })

  output$pk_params_table <- renderTable({
    df   <- pk_conc_df()
    Cmax <- max(df$conc)
    t_max <- df$time[which.max(df$conc)]
    AUC  <- sum(diff(df$time) * (head(df$conc, -1) + tail(df$conc, -1)) / 2)
    t50  <- if (any(df$conc < Cmax / 2 & df$time > t_max)) {
      df$time[min(which(df$conc < Cmax / 2 & df$time > t_max))]
    } else { NA }
    data.frame(
      Parameter = c("Cmax", "Tmax (h)", "AUC (0-t)", "t½ (approx, h)"),
      Value     = c(round(Cmax, 3), round(t_max, 1),
                    round(AUC, 2), round(t50, 1))
    )
  }, striped = TRUE, bordered = TRUE)

  # ---------------------------------------------------------------------------
  # TAB 3 OUTPUTS — IgG4 & Immune
  # ---------------------------------------------------------------------------
  igg4_time_vec <- reactive({
    seq(0, input$igg4_horizon, by = 1)
  })

  igg4_vec_r <- reactive({
    igg4_dynamics(igg4_time_vec(), input$baseline_igg4,
                  input$drug_igg4, input$igg4_dose, input$weight,
                  input$aip_type, input$ep_check)
  })

  output$igg4_plot <- renderPlot({
    tw  <- igg4_time_vec()
    ig  <- igg4_vec_r()
    df  <- data.frame(week = tw, igg4 = ig)

    ggplot(df, aes(x = week, y = igg4)) +
      geom_ribbon(aes(ymin = pmin(igg4, input$igg4_uln),
                      ymax = igg4), fill = "#e74c3c", alpha = 0.18) +
      geom_line(color = "#c0392b", linewidth = 1.3) +
      geom_hline(yintercept = input$igg4_uln, linetype = "dashed",
                 color = "#e74c3c", linewidth = 0.9) +
      geom_hline(yintercept = input$igg4_remission, linetype = "dotted",
                 color = "#27ae60", linewidth = 0.9) +
      annotate("text", x = max(tw) * 0.95, y = input$igg4_uln + 15,
               label = "ULN (135 mg/dL)", color = "#e74c3c", size = 3.5, hjust = 1) +
      annotate("text", x = max(tw) * 0.95, y = input$igg4_remission - 15,
               label = "Remission target", color = "#27ae60", size = 3.5, hjust = 1) +
      labs(title = paste("Serum IgG4 —", input$drug_igg4),
           x = "Weeks", y = "IgG4 (mg/dL)") +
      theme_minimal(base_size = 13) +
      theme(plot.title = element_text(face = "bold", color = "#2c3e50"))
  })

  output$bcell_plot <- renderPlot({
    tw <- igg4_time_vec()
    bc <- bcell_depletion(tw, input$igg4_dose)
    df <- data.frame(week = tw, bcell = bc)
    ggplot(df, aes(x = week, y = bcell)) +
      geom_line(color = "#8e44ad", linewidth = 1.3) +
      geom_area(fill = "#8e44ad", alpha = 0.12) +
      geom_hline(yintercept = 20, linetype = "dashed", color = "#8e44ad") +
      labs(title = "B-Cell Count (% of baseline) — Rituximab Depletion",
           x = "Weeks", y = "B cells (% baseline)") +
      theme_minimal(base_size = 13)
  })

  output$response_classification <- renderUI({
    ig      <- igg4_vec_r()
    tw      <- igg4_time_vec()
    idx_6m  <- which.min(abs(tw - 26))
    idx_12m <- which.min(abs(tw - 52))
    igg4_6m  <- if (idx_6m  <= length(ig)) ig[idx_6m]  else tail(ig, 1)
    igg4_12m <- if (idx_12m <= length(ig)) ig[idx_12m] else tail(ig, 1)

    resp_6m <- if (igg4_6m < 135)       "Complete" else if (igg4_6m < input$baseline_igg4 * 0.5) "Partial" else "No response"
    resp_col <- if (resp_6m == "Complete") "#27ae60" else if (resp_6m == "Partial") "#f39c12" else "#e74c3c"

    div(
      fluidRow(
        column(4,
          div(class = "patient-card",
            strong("6-Month IgG4:"), br(),
            tags$span(style = "font-size:18px; font-weight:bold;",
                      round(igg4_6m, 0), " mg/dL")
          )
        ),
        column(4,
          div(class = "patient-card",
            strong("12-Month IgG4:"), br(),
            tags$span(style = "font-size:18px; font-weight:bold;",
                      round(igg4_12m, 0), " mg/dL")
          )
        ),
        column(4,
          div(class = "patient-card",
            strong("Treatment Response:"), br(),
            tags$span(style = paste("font-size:18px; font-weight:bold; color:", resp_col),
                      resp_6m)
          )
        )
      )
    )
  })

  # ---------------------------------------------------------------------------
  # TAB 4 OUTPUTS — Pancreatic Function
  # ---------------------------------------------------------------------------
  panc_igg4 <- reactive({
    tw <- seq(0, input$panc_horizon, by = 1)
    igg4_dynamics(tw, input$baseline_igg4, input$drug_panc,
                  input$panc_dose, input$weight, input$aip_type, FALSE)
  })

  panc_fe1 <- reactive({
    tw <- seq(0, input$panc_horizon, by = 1)
    exocrine_dynamics(tw, input$pd_baseline_fe1,
                      input$drug_panc, input$panc_dose, input$weight, panc_igg4())
  })

  panc_hba1c <- reactive({
    tw <- seq(0, input$panc_horizon, by = 1)
    hba1c_dynamics(tw, input$pd_baseline_hba1c,
                   input$drug_panc, input$panc_dose, panc_igg4())
  })

  output$fe1_plot <- renderPlot({
    tw <- seq(0, input$panc_horizon, by = 1)
    df <- data.frame(week = tw, fe1 = panc_fe1())
    ggplot(df, aes(x = week, y = fe1)) +
      geom_line(color = "#16a085", linewidth = 1.3) +
      geom_hline(yintercept = input$fe1_threshold, linetype = "dashed",
                 color = "#16a085") +
      annotate("text", x = max(tw) * 0.05, y = input$fe1_threshold + 10,
               label = "Normal threshold (200 μg/g)",
               color = "#16a085", size = 3.5, hjust = 0) +
      labs(title = "Exocrine Function (FE-1 Proxy)",
           x = "Weeks", y = "FE-1 (μg/g)") +
      theme_minimal(base_size = 13) +
      theme(plot.title = element_text(face = "bold", color = "#2c3e50"))
  })

  output$hba1c_plot <- renderPlot({
    tw <- seq(0, input$panc_horizon, by = 1)
    df <- data.frame(week = tw, hba1c = panc_hba1c())
    ggplot(df, aes(x = week, y = hba1c)) +
      geom_line(color = "#e67e22", linewidth = 1.3) +
      geom_hline(yintercept = 6.5, linetype = "dashed", color = "#e67e22") +
      geom_hline(yintercept = input$hba1c_threshold, linetype = "dotted",
                 color = "#e74c3c") +
      annotate("text", x = max(tw) * 0.05, y = 6.65,
               label = "Prediabetes (6.5%)", color = "#e67e22", size = 3.5, hjust = 0) +
      labs(title = "Endocrine Function (HbA1c)",
           x = "Weeks", y = "HbA1c (%)") +
      theme_minimal(base_size = 13) +
      theme(plot.title = element_text(face = "bold", color = "#2c3e50"))
  })

  output$organ_score_plot <- renderPlot({
    tw   <- seq(0, input$panc_horizon, by = 1)
    fe   <- panc_fe1()
    ha   <- panc_hba1c()
    ig   <- panc_igg4()

    fe_score  <- pmin(fe / 200, 1)
    ha_score  <- pmax(1 - (ha - 5.7) / 3, 0)
    igg_score <- pmax(1 - pmax(ig - 135, 0) / 1000, 0)
    overall   <- (fe_score + ha_score + igg_score) / 3

    df <- data.frame(week = tw, score = overall)
    ggplot(df, aes(x = week, y = score)) +
      geom_line(color = "#2980b9", linewidth = 1.3) +
      geom_hline(yintercept = 0.8, linetype = "dashed", color = "#27ae60",
                 linewidth = 0.8) +
      scale_y_continuous(limits = c(0, 1), labels = percent_format()) +
      annotate("text", x = max(tw) * 0.85, y = 0.83,
               label = "Target (80%)", color = "#27ae60", size = 3.5) +
      labs(title = "Combined Organ Function Score",
           x = "Weeks", y = "Score (0–1)") +
      theme_minimal(base_size = 13) +
      theme(plot.title = element_text(face = "bold", color = "#2c3e50"))
  })

  # ---------------------------------------------------------------------------
  # TAB 5 OUTPUTS — Scenario Comparison
  # ---------------------------------------------------------------------------
  scenario_data <- eventReactive(input$run_scenarios, {
    tw  <- seq(0, 52, by = 1)
    wt  <- isolate(input$weight)
    bl  <- isolate(input$baseline_igg4)
    at  <- isolate(input$aip_type)

    scenarios <- list(
      list(drug = input$sc1_drug, dose = input$sc1_dose, label = paste("Sc1:", input$sc1_drug)),
      list(drug = input$sc2_drug, dose = input$sc2_dose, label = paste("Sc2:", input$sc2_drug)),
      list(drug = input$sc3_drug, dose = input$sc3_dose, label = paste("Sc3:", input$sc3_drug)),
      list(drug = input$sc4_drug, dose = input$sc4_dose, label = paste("Sc4:", input$sc4_drug))
    )

    lapply(scenarios, function(sc) {
      ig <- igg4_dynamics(tw, bl, sc$drug, sc$dose, wt, at, FALSE)
      data.frame(week = tw, igg4 = ig, scenario = sc$label)
    })
  }, ignoreNULL = FALSE)

  output$scenario_igg4_plot <- renderPlot({
    sd_list <- scenario_data()
    df <- bind_rows(sd_list)
    colors_sc <- c("#e74c3c","#3498db","#2ecc71","#f39c12")
    ggplot(df, aes(x = week, y = igg4, color = scenario)) +
      geom_line(linewidth = 1.2) +
      geom_hline(yintercept = 135, linetype = "dashed", color = "grey40") +
      scale_color_manual(values = colors_sc) +
      labs(title = "Serum IgG4 — Scenario Comparison",
           x = "Weeks", y = "IgG4 (mg/dL)", color = "Scenario") +
      theme_minimal(base_size = 13) +
      theme(legend.position = "right",
            plot.title = element_text(face = "bold", color = "#2c3e50"))
  })

  output$scenario_summary_table <- renderTable({
    sd_list <- scenario_data()
    bl      <- isolate(input$baseline_igg4)
    tw      <- seq(0, 52, by = 1)

    rows <- lapply(sd_list, function(df) {
      ig    <- df$igg4
      igg4_6m  <- ig[which.min(abs(tw - 26))]
      igg4_12m <- ig[which.min(abs(tw - 52))]

      # Time to remission (first week IgG4 < 135)
      rem_wk <- tw[which(ig < 135)[1]]
      rem_str <- if (!is.na(rem_wk)) paste(round(rem_wk, 0), "wks") else "Not achieved"

      # Relapse rate proxy
      relapse_rate <- if (tail(ig, 1) > 135) "~35%" else "~12%"

      data.frame(
        Scenario         = df$scenario[1],
        `IgG4 at 6m`    = round(igg4_6m, 0),
        `IgG4 at 12m`   = round(igg4_12m, 0),
        `Time to Remission` = rem_str,
        `Est. Relapse Rate` = relapse_rate
      )
    })
    bind_rows(rows)
  }, striped = TRUE, bordered = TRUE)

  output$forest_plot <- renderPlot({
    sd_list <- scenario_data()
    tw      <- seq(0, 52, by = 1)
    ref_ig  <- sd_list[[1]]$igg4[which.min(abs(tw - 26))]

    df <- lapply(sd_list, function(d) {
      ig6 <- d$igg4[which.min(abs(tw - 26))]
      rr  <- if (ref_ig > 0) ig6 / ref_ig else 1
      data.frame(scenario = d$scenario[1], RR = rr,
                 lo = rr * 0.8, hi = rr * 1.2)
    }) %>% bind_rows()
    df$scenario <- factor(df$scenario, levels = rev(df$scenario))

    ggplot(df, aes(x = RR, y = scenario)) +
      geom_vline(xintercept = 1, linetype = "dashed", color = "grey60") +
      geom_errorbarh(aes(xmin = lo, xmax = hi), height = 0.25,
                     color = "#2c3e50", linewidth = 0.9) +
      geom_point(color = "#e74c3c", size = 4) +
      labs(title = "IgG4 Ratio at 6 Months vs. Reference (Scenario 1)",
           x = "Ratio (1 = same as reference)", y = NULL) +
      theme_minimal(base_size = 13) +
      theme(plot.title = element_text(face = "bold", color = "#2c3e50"))
  })

  # ---------------------------------------------------------------------------
  # TAB 6 OUTPUTS — Biomarkers & Monitoring
  # ---------------------------------------------------------------------------
  output$monitoring_plot <- renderPlot({
    int_wk <- as.integer(input$monitor_interval)
    total  <- 52
    visits <- seq(0, total, by = int_wk)

    # Generate "simulated" biomarker trajectories for visualization
    tw     <- seq(0, total, by = 1)
    ig_sim <- igg4_dynamics(tw, input$baseline_igg4, "Prednisolone",
                            40, input$weight, "Type1", input$bm_extra)

    visit_df <- data.frame(
      week  = visits,
      igg4  = ig_sim[pmin(visits + 1, length(ig_sim))],
      ca199 = pmax(input$bm_ca199_6m * exp(-0.02 * visits) + runif(length(visits), -2, 2), 0),
      fe1   = exocrine_dynamics(visits, input$baseline_fe1, "Prednisolone",
                                40, input$weight, ig_sim[pmin(visits + 1, length(ig_sim))])
    )

    df_long <- visit_df %>%
      pivot_longer(cols = c(igg4, ca199, fe1),
                   names_to = "biomarker", values_to = "value") %>%
      mutate(biomarker = recode(biomarker,
                                "igg4"  = "IgG4 (mg/dL)",
                                "ca199" = "CA 19-9 (U/mL)",
                                "fe1"   = "FE-1 (μg/g)"))

    ggplot(df_long, aes(x = week, y = value, color = biomarker, group = biomarker)) +
      geom_line(linewidth = 0.9, alpha = 0.7) +
      geom_point(size = 2.5) +
      facet_wrap(~ biomarker, scales = "free_y", ncol = 3) +
      labs(title = "Monitoring Schedule — Biomarker Time Course",
           x = "Week", y = "Value") +
      theme_minimal(base_size = 12) +
      theme(legend.position = "none",
            plot.title      = element_text(face = "bold", color = "#2c3e50"),
            strip.text      = element_text(face = "bold"))
  })

  output$biomarker_threshold_table <- renderTable({
    data.frame(
      Biomarker   = c("Serum IgG4", "CA 19-9", "FE-1", "HbA1c", "IgG4/IgG ratio"),
      Threshold   = c("135 mg/dL", "37 U/mL", "200 μg/g", "6.5%", "> 10%"),
      Direction   = c("Elevated = abnormal","Elevated = suspicious",
                       "Below = exocrine insufficiency",
                       "Above = diabetes risk",
                       "Elevated = AIP support"),
      Action      = c("Remission check","Rule out malignancy",
                       "Enzyme replacement", "Diabetes management",
                       "Diagnostic support"),
      Status      = c(
        ifelse(input$bm_igg4_6m > 135, "ABOVE", "Normal"),
        ifelse(input$bm_ca199_6m > 37, "ABOVE", "Normal"),
        ifelse(input$bm_fe1_6m < 200,  "BELOW", "Normal"),
        ifelse(input$bm_hba1c_6m > 6.5,"ABOVE", "Normal"),
        "N/A"
      )
    )
  }, striped = TRUE, bordered = TRUE)

  output$relapse_risk_ui <- renderUI({
    risk <- relapse_risk(
      igg4_6m         = input$bm_igg4_6m,
      extra_pancreatic = input$bm_extra,
      steroid_dose_at_6m = input$bm_steroid_dose
    )

    risk_col   <- switch(risk, "Low" = "#27ae60", "Medium" = "#f39c12", "High" = "#e74c3c")
    risk_class <- switch(risk, "Low" = "risk-low", "Medium" = "risk-medium", "High" = "risk-high")

    rec <- switch(risk,
      "Low"    = "Standard monitoring every 6 months. Continue maintenance dose.",
      "Medium" = "Increase monitoring to every 3 months. Consider maintenance therapy.",
      "High"   = "Close monitoring monthly. Early re-treatment if relapse confirmed."
    )

    div(
      fluidRow(
        column(4,
          div(class = "patient-card",
            h4("Relapse Risk Score"),
            tags$span(class = risk_class, style = "font-size:32px;", risk)
          )
        ),
        column(8,
          div(class = "patient-card",
            strong("Clinical Recommendation:"), br(),
            tags$span(rec), br(), br(),
            strong("Risk Factors Present:"), br(),
            tags$ul(
              if (input$bm_igg4_6m > 135) tags$li("IgG4 still elevated at 6 months (> 135 mg/dL)"),
              if (input$bm_extra)          tags$li("Extra-pancreatic involvement"),
              if (input$bm_steroid_dose > 5) tags$li("High steroid dose still required")
            )
          )
        )
      )
    )
  })
}

# =============================================================================
# LAUNCH APPLICATION
# =============================================================================

shinyApp(ui = ui, server = server)
