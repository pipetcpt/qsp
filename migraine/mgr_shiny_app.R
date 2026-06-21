##############################################################################
# Migraine QSP — Interactive Shiny Dashboard
# 6 Tabs: Patient Profile · PK · CGRP/PD · Clinical Endpoints ·
#         Scenario Comparison · Biomarkers
##############################################################################

library(shiny)
library(shinydashboard)
library(ggplot2)
library(dplyr)
library(tidyr)
library(plotly)
library(DT)

# ── Helper: simple 2-compartment PK ODE solver (Euler) ───────────────────────
pk_2cmt <- function(dose, F_abs, ka, Vc, Vp, CL, Q,
                    t_end = 48, dt = 0.1, t_dose = 0) {
  times <- seq(0, t_end, by = dt)
  depot <- Cp <- Ct <- numeric(length(times))
  depot[1] <- dose * F_abs
  for (i in seq_len(length(times) - 1)) {
    d_depot <- -ka * depot[i]
    d_Cp    <-  ka * depot[i] - (CL + Q) / Vc * Cp[i] + Q / Vp * Ct[i]
    d_Ct    <-  Q / Vc * Cp[i] - Q / Vp * Ct[i]
    depot[i+1] <- depot[i] + d_depot * dt
    Cp[i+1]    <- Cp[i]    + d_Cp    * dt
    Ct[i+1]    <- Ct[i]    + d_Ct    * dt
  }
  data.frame(time = times, Cp = pmax(Cp, 0), Ct = pmax(Ct, 0))
}

# ── Helper: simple pain ODE (2-state: CGRP, pain) ───────────────────────────
sim_pain <- function(drug_occ_fn, t_end = 24, dt = 0.05) {
  times   <- seq(0, t_end, by = dt)
  CGRP    <- numeric(length(times))
  pain    <- numeric(length(times))
  CGRP[1] <- 2.5   # elevated at attack onset
  pain[1] <- 7.0   # moderate pain at onset
  for (i in seq_len(length(times) - 1)) {
    occ  <- drug_occ_fn(times[i])
    dC   <- 0.05 - 0.15 * CGRP[i] * (1 + 2 * occ)
    dpain <- 0.30 * (CGRP[i] / (CGRP[i] + 0.5)) * (1 - occ) * 10 -
             0.05 * pain[i]
    CGRP[i+1]  <- pmax(CGRP[i]  + dC    * dt, 0)
    pain[i+1]  <- pmin(pmax(pain[i] + dpain * dt, 0), 10)
  }
  data.frame(time = times, CGRP = CGRP, pain = pain)
}

##############################################################################
# UI
##############################################################################
ui <- dashboardPage(
  skin = "purple",

  dashboardHeader(
    title = "Migraine QSP Dashboard",
    titleWidth = 280
  ),

  dashboardSidebar(
    width = 280,
    sidebarMenu(
      menuItem("Patient Profile",      tabName = "patient",   icon = icon("user")),
      menuItem("PK Profiles",          tabName = "pk",        icon = icon("chart-line")),
      menuItem("CGRP / PD Dynamics",   tabName = "pd",        icon = icon("dna")),
      menuItem("Clinical Endpoints",   tabName = "endpoints", icon = icon("hospital")),
      menuItem("Scenario Comparison",  tabName = "compare",   icon = icon("layer-group")),
      menuItem("Biomarker Dashboard",  tabName = "biomarker", icon = icon("flask"))
    ),

    hr(),
    h5("Global Settings", style = "padding-left:15px; color:#ddd"),

    sliderInput("BW",   "Body Weight (kg)",   40, 120, 70, step = 5),
    sliderInput("age",  "Age (years)",         18, 80,  35, step = 1),
    selectInput("migr_type", "Migraine Type",
                choices = c("Episodic (<15 MHD)"   = "episodic",
                            "Chronic (≥15 MHD)"    = "chronic",
                            "Hemiplegic"            = "hemiplegic"),
                selected = "episodic"),
    sliderInput("baseline_MMD", "Baseline MMD",  2, 25, 10, step = 1)
  ),

  dashboardBody(
    tags$head(tags$style(HTML("
      .content-wrapper { background-color: #f4f4f8; }
      .box { border-radius: 8px; }
      .value-box { border-radius: 8px; }
    "))),

    tabItems(

      # ── TAB 1: PATIENT PROFILE ────────────────────────────────────────────
      tabItem(tabName = "patient",
        fluidRow(
          valueBoxOutput("vb_mhd"),
          valueBoxOutput("vb_midas"),
          valueBoxOutput("vb_cgrp_risk")
        ),
        fluidRow(
          box(title = "Patient Parameters", width = 4, status = "primary",
            solidHeader = TRUE,
            selectInput("sex", "Sex", c("Female", "Male"), "Female"),
            sliderInput("attack_dur", "Typical Attack Duration (h)", 4, 72, 24),
            sliderInput("pain_peak", "Peak Pain VAS",  0, 10, 8),
            checkboxGroupInput("assoc_sx", "Associated Symptoms",
              choices = c("Nausea/Vomiting", "Photophobia", "Phonophobia",
                          "Aura (visual)", "Allodynia", "Osmophobia"),
              selected = c("Nausea/Vomiting", "Photophobia", "Phonophobia")),
            checkboxGroupInput("triggers", "Known Triggers",
              choices = c("Hormonal (menstrual)", "Stress", "Sleep disruption",
                          "Diet (alcohol/tyramine)", "Weather", "Fasting"),
              selected = c("Stress", "Sleep disruption")),
            checkboxGroupInput("comorbid", "Comorbidities",
              choices = c("Depression/Anxiety", "Cardiovascular disease",
                          "Epilepsy", "Obesity", "Hypertension",
                          "Medication overuse"),
              selected = c())
          ),
          box(title = "CGRP Pathway Risk Score", width = 4, status = "warning",
            solidHeader = TRUE,
            plotlyOutput("radar_risk", height = "320px")
          ),
          box(title = "Treatment History", width = 4, status = "info",
            solidHeader = TRUE,
            checkboxGroupInput("prev_tx", "Prior Preventive Agents",
              choices = c("Topiramate", "Propranolol", "Amitriptyline",
                          "Valproate", "Candesartan", "CGRP mAb (1st line)",
                          "CGRP mAb (2nd line)"),
              selected = c()),
            checkboxGroupInput("acute_tx", "Current Acute Treatments",
              choices = c("NSAID", "Triptan", "Gepant (rimegepant/ubrogepant)",
                          "Lasmiditan", "Ergot", "Opioid"),
              selected = c("Triptan")),
            br(),
            div(class = "well",
              h5("MIDAS Score Interpretation"),
              p("0–5: Little/No disability"),
              p("6–10: Mild disability"),
              p("11–20: Moderate disability"),
              p("≥21: Severe disability")
            )
          )
        ),
        fluidRow(
          box(title = "Mechanistic Pathway Activity", width = 12, status = "success",
            solidHeader = TRUE,
            plotlyOutput("pathway_bar", height = "250px")
          )
        )
      ),

      # ── TAB 2: PK PROFILES ───────────────────────────────────────────────
      tabItem(tabName = "pk",
        fluidRow(
          box(title = "Acute Treatment PK Settings", width = 3, status = "primary",
            solidHeader = TRUE,
            h4("Sumatriptan"),
            radioButtons("sum_route", "Route",
                         c("SC 6 mg" = "sc", "Oral 100 mg" = "oral"),
                         "sc"),
            numericInput("sum_dose", "Custom Dose (mg)", 6, 1, 100),

            h4("Lasmiditan"),
            numericInput("lam_dose", "Dose (mg)", 200, 50, 400, step = 50),

            h4("Rimegepant"),
            numericInput("rim_dose", "Dose (mg)", 75, 50, 150, step = 25),

            h4("Ubrogepant"),
            numericInput("ubr_dose", "Dose (mg)", 100, 50, 200, step = 50)
          ),
          box(title = "Acute Therapy PK — Concentration vs Time", width = 9,
            status = "primary", solidHeader = TRUE,
            plotlyOutput("pk_acute_plot", height = "420px")
          )
        ),
        fluidRow(
          box(title = "Preventive Treatment PK Settings", width = 3,
            status = "success", solidHeader = TRUE,
            h4("Erenumab"),
            radioButtons("ere_dose_mg", "Dose",
                         c("70 mg SC" = 70, "140 mg SC" = 140), 140),
            radioButtons("ere_freq", "Frequency", c("Monthly" = 720), 720),
            numericInput("ere_months", "Simulation Duration (months)", 3, 1, 12),

            h4("Fremanezumab"),
            radioButtons("frem_reg", "Regimen",
                         c("225 mg monthly" = "monthly",
                           "675 mg quarterly" = "quarterly"),
                         "monthly"),

            h4("Topiramate"),
            sliderInput("top_dose", "Daily Dose (mg/day)", 25, 200, 100, step = 25),

            h4("Propranolol"),
            sliderInput("prop_dose", "Daily Dose (mg/day)", 40, 240, 160, step = 40)
          ),
          box(title = "Preventive Therapy PK — Concentration vs Time", width = 9,
            status = "success", solidHeader = TRUE,
            plotlyOutput("pk_prev_plot", height = "420px")
          )
        )
      ),

      # ── TAB 3: CGRP / PD DYNAMICS ────────────────────────────────────────
      tabItem(tabName = "pd",
        fluidRow(
          box(title = "CGRP Model Parameters", width = 3, status = "danger",
            solidHeader = TRUE,
            sliderInput("ksyn_CGRP", "CGRP Synthesis Rate (pmol/min)", 0.01, 0.20, 0.05),
            sliderInput("kdeg_CGRP", "CGRP Degradation Rate (1/min)",  0.05, 0.50, 0.15),
            sliderInput("CSD_init",  "CSD Initiation Amplitude", 0.0, 1.0, 0.8),
            sliderInput("EC50_CGRP", "CGRP EC50 at CLR/RAMP1 (pmol/L)", 0.1, 3.0, 0.5),
            hr(),
            h5("Drug Effect Parameters"),
            sliderInput("ere_KD",   "Erenumab KD (nM)",     0.001, 0.1, 0.01),
            sliderInput("rim_Ki",   "Rimegepant Ki (nM)",    0.01,  1.0, 0.027),
            sliderInput("sum_EC50", "Sumatriptan EC50 (nM)", 1.0,  10.0, 3.2)
          ),
          box(title = "CGRP Dynamics During Attack", width = 9,
            status = "danger", solidHeader = TRUE,
            plotlyOutput("pd_cgrp_plot", height = "380px")
          )
        ),
        fluidRow(
          box(title = "CGRP Receptor Occupancy vs Drug Concentration", width = 6,
            status = "warning", solidHeader = TRUE,
            plotlyOutput("pd_occ_plot", height = "300px")
          ),
          box(title = "CSD Propagation Model", width = 6,
            status = "warning", solidHeader = TRUE,
            plotlyOutput("pd_csd_plot", height = "300px")
          )
        )
      ),

      # ── TAB 4: CLINICAL ENDPOINTS ─────────────────────────────────────────
      tabItem(tabName = "endpoints",
        fluidRow(
          box(title = "Pain Score Over Time (VAS)", width = 6, status = "danger",
            solidHeader = TRUE,
            plotlyOutput("ep_pain_plot", height = "320px")
          ),
          box(title = "Monthly Migraine Days — Preventive Response", width = 6,
            status = "info", solidHeader = TRUE,
            plotlyOutput("ep_mmd_plot", height = "320px")
          )
        ),
        fluidRow(
          box(title = "Preventive Responder Rates (≥50% MMD reduction)", width = 6,
            status = "success", solidHeader = TRUE,
            plotlyOutput("ep_responder", height = "300px")
          ),
          box(title = "Clinical Trial Benchmark Data", width = 6,
            status = "primary", solidHeader = TRUE,
            DTOutput("trial_table")
          )
        )
      ),

      # ── TAB 5: SCENARIO COMPARISON ────────────────────────────────────────
      tabItem(tabName = "compare",
        fluidRow(
          box(title = "Acute Treatment Comparison", width = 12,
            status = "primary", solidHeader = TRUE,
            column(4,
              h4("Select Treatments to Compare"),
              checkboxGroupInput("acute_compare",
                label = NULL,
                choices = c("Untreated",
                            "Sumatriptan SC 6 mg",
                            "Sumatriptan Oral 100 mg",
                            "Lasmiditan 200 mg",
                            "Rimegepant 75 mg",
                            "Ubrogepant 100 mg",
                            "NSAID (naproxen 500 mg)"),
                selected = c("Untreated", "Sumatriptan SC 6 mg",
                             "Lasmiditan 200 mg", "Rimegepant 75 mg")),
              hr(),
              h4("Key Outcome"),
              radioButtons("compare_outcome", NULL,
                           c("Pain VAS" = "pain",
                             "CGRP Level" = "cgrp",
                             "TG Activation" = "tg",
                             "CGRP-R Occupancy" = "occ"),
                           "pain")
            ),
            column(8,
              plotlyOutput("compare_plot", height = "400px")
            )
          )
        ),
        fluidRow(
          box(title = "Preventive Treatment Comparison (6-Month Simulation)", width = 12,
            status = "success", solidHeader = TRUE,
            column(4,
              checkboxGroupInput("prev_compare",
                label = "Select Preventive Agents",
                choices = c("No treatment",
                            "Erenumab 70 mg",
                            "Erenumab 140 mg",
                            "Fremanezumab 225 mg monthly",
                            "Galcanezumab 120 mg",
                            "Topiramate 100 mg/day",
                            "Propranolol 160 mg/day",
                            "Amitriptyline 50 mg/day"),
                selected = c("No treatment", "Erenumab 140 mg",
                             "Topiramate 100 mg/day")),
              sliderInput("prev_sim_months", "Simulation Duration (months)",
                          1, 12, 6)
            ),
            column(8,
              plotlyOutput("prev_compare_plot", height = "380px")
            )
          )
        )
      ),

      # ── TAB 6: BIOMARKER DASHBOARD ────────────────────────────────────────
      tabItem(tabName = "biomarker",
        fluidRow(
          valueBoxOutput("vb_cgrp_plasma"),
          valueBoxOutput("vb_5ht_platelet"),
          valueBoxOutput("vb_pge2")
        ),
        fluidRow(
          box(title = "Plasma CGRP Trajectory", width = 6, status = "danger",
            solidHeader = TRUE,
            plotlyOutput("bio_cgrp_plot", height = "300px")
          ),
          box(title = "Platelet 5-HT Dynamics", width = 6, status = "info",
            solidHeader = TRUE,
            plotlyOutput("bio_5ht_plot", height = "300px")
          )
        ),
        fluidRow(
          box(title = "PGE2 & Neuroinflammation Markers", width = 6,
            status = "warning", solidHeader = TRUE,
            plotlyOutput("bio_pge2_plot", height = "300px")
          ),
          box(title = "Biomarker Reference Ranges", width = 6,
            status = "primary", solidHeader = TRUE,
            DTOutput("bio_ref_table")
          )
        )
      )
    )
  )
)

##############################################################################
# SERVER
##############################################################################
server <- function(input, output, session) {

  # ── Reactive: attack/patient parameters ──────────────────────────────────
  patient_params <- reactive({
    list(
      BW          = input$BW,
      age         = input$age,
      migr_type   = input$migr_type,
      baseline_MMD = input$baseline_MMD,
      pain_peak   = input$pain_peak,
      attack_dur  = input$attack_dur
    )
  })

  # ── Value boxes ──────────────────────────────────────────────────────────
  output$vb_mhd <- renderValueBox({
    valueBox(
      value = input$baseline_MMD,
      subtitle = "Monthly Headache Days (MHD)",
      icon = icon("calendar"),
      color = if (input$baseline_MMD >= 15) "red" else if (input$baseline_MMD >= 8) "orange" else "green"
    )
  })

  output$vb_midas <- renderValueBox({
    midas_est <- input$baseline_MMD * input$pain_peak * 0.3
    valueBox(
      value = round(midas_est),
      subtitle = "Estimated MIDAS Score",
      icon = icon("wheelchair"),
      color = if (midas_est >= 21) "red" else if (midas_est >= 11) "orange" else "yellow"
    )
  })

  output$vb_cgrp_risk <- renderValueBox({
    risk <- min(100, round(
      input$baseline_MMD * 3 +
        (input$migr_type == "chronic") * 20 +
        length(input$triggers) * 5
    ))
    valueBox(
      value = paste0(risk, "%"),
      subtitle = "CGRP-Pathway Activation Score",
      icon = icon("dna"),
      color = if (risk >= 70) "red" else if (risk >= 40) "orange" else "yellow"
    )
  })

  output$vb_cgrp_plasma <- renderValueBox({
    valueBox(
      value = "↑ 2–4×",
      subtitle = "Ictal CGRP vs Interictal (pg/mL)",
      icon = icon("arrow-up"),
      color = "red"
    )
  })
  output$vb_5ht_platelet <- renderValueBox({
    valueBox(
      value = "↓ during attack",
      subtitle = "Platelet 5-HT (ng/mL)",
      icon = icon("arrow-down"),
      color = "blue"
    )
  })
  output$vb_pge2 <- renderValueBox({
    valueBox(
      value = "↑ COX-2",
      subtitle = "PGE2 / Neurogenic Inflammation",
      icon = icon("fire"),
      color = "orange"
    )
  })

  # ── Radar chart (risk factors) ────────────────────────────────────────────
  output$radar_risk <- renderPlotly({
    categories <- c("CGRP Activity", "CSD Risk", "Central Sensitization",
                    "TG Activation", "Serotonin Deficiency", "Hormonal Influence")
    vals <- c(
      min(100, input$baseline_MMD * 4),
      min(100, (input$migr_type == "chronic") * 60 + 20),
      min(100, (input$migr_type == "chronic") * 50 + length(input$assoc_sx) * 8),
      min(100, input$pain_peak * 8),
      min(100, 30 + (input$baseline_MMD > 10) * 30),
      min(100, ("Hormonal (menstrual)" %in% input$triggers) * 70 + 10)
    )
    plot_ly(
      type = "scatterpolar",
      r = c(vals, vals[1]),
      theta = c(categories, categories[1]),
      fill = "toself",
      fillcolor = "rgba(123,45,139,0.3)",
      line = list(color = "#7b2d8b")
    ) %>%
      layout(polar = list(radialaxis = list(visible = TRUE, range = c(0, 100))),
             margin = list(t = 20, b = 20))
  })

  # ── Pathway bar ───────────────────────────────────────────────────────────
  output$pathway_bar <- renderPlotly({
    pathways <- data.frame(
      pathway = c("CGRP/CLR-RAMP1", "5-HT1B/1D", "CSD/Cortical",
                  "PGE2/COX-2", "NO/sGC", "NK1/Substance P",
                  "TRPV1/TRPA1", "NF-κB/Cytokines"),
      activity = c(
        min(1, input$baseline_MMD / 15),
        0.6,
        (input$migr_type == "chronic") * 0.8 + 0.2,
        0.65,
        0.50,
        0.45,
        0.55,
        (input$migr_type == "chronic") * 0.7 + 0.2
      )
    )
    plot_ly(pathways, x = ~activity, y = ~reorder(pathway, activity),
            type = "bar", orientation = "h",
            marker = list(color = "#7b2d8b")) %>%
      layout(xaxis = list(title = "Relative Activity (0–1)", range = c(0, 1)),
             yaxis = list(title = ""),
             margin = list(l = 150))
  })

  # ── PK: Acute treatments ─────────────────────────────────────────────────
  output$pk_acute_plot <- renderPlotly({
    BW <- input$BW
    t_end <- 24

    # Sumatriptan
    if (input$sum_route == "sc") {
      pk_sum <- pk_2cmt(input$sum_dose, F_abs = 0.97, ka = 3.0,
                        Vc = 2.4*BW, Vp = 1.8*BW, CL = 72, Q = 20,
                        t_end = t_end)
    } else {
      pk_sum <- pk_2cmt(input$sum_dose, F_abs = 0.14, ka = 0.8,
                        Vc = 2.4*BW, Vp = 1.8*BW, CL = 72, Q = 20,
                        t_end = t_end)
    }
    pk_sum$drug <- paste("Sumatriptan", if(input$sum_route=="sc") "SC" else "Oral")
    pk_sum$Cp_nM <- pk_sum$Cp * 1000 / 413.5

    # Lasmiditan
    pk_lam <- data.frame(
      time = seq(0, t_end, by = 0.1),
      Cp = input$lam_dose * 0.38 / 140 *
             exp(-0.011 * seq(0, t_end, by = 0.1))
    )
    pk_lam$drug <- "Lasmiditan Oral"
    pk_lam$Cp_nM <- pk_lam$Cp * 1e6 / 439.4

    # Rimegepant
    pk_rim <- data.frame(
      time = seq(0, t_end, by = 0.1)
    ) %>%
      mutate(
        Cp = input$rim_dose * 0.64 / 113 *
               (exp(-0.665 * time) - exp(-0.8 * time)) * 6,
        drug = "Rimegepant Oral",
        Cp_nM = pmax(Cp * 1e6 / 534.6, 0)
      )

    # Ubrogepant
    pk_ubr <- data.frame(
      time = seq(0, t_end, by = 0.1),
      Cp = input$ubr_dose * 0.44 / 350 * exp(-0.12 * seq(0, t_end, 0.1))
    )
    pk_ubr$drug <- "Ubrogepant Oral"
    pk_ubr$Cp_nM <- pk_ubr$Cp * 1e6 / 511.6

    df_all <- bind_rows(
      select(pk_sum, time, Cp_nM, drug),
      select(pk_lam, time, Cp_nM, drug),
      select(pk_rim, time, Cp_nM, drug),
      select(pk_ubr, time, Cp_nM, drug)
    )

    plot_ly(df_all, x = ~time, y = ~Cp_nM, color = ~drug,
            type = "scatter", mode = "lines",
            line = list(width = 2)) %>%
      layout(
        xaxis = list(title = "Time (hours)"),
        yaxis = list(title = "Plasma Concentration (nM)"),
        hovermode = "x unified"
      )
  })

  # ── PK: Preventive ───────────────────────────────────────────────────────
  output$pk_prev_plot <- renderPlotly({
    n_months <- input$ere_months
    t_end_h  <- n_months * 30 * 24

    # Erenumab (simplified multi-dose 1-cmt)
    ere_dose <- as.numeric(input$ere_dose_mg)
    times_ere <- seq(0, t_end_h, by = 6)
    doses_h   <- seq(0, t_end_h, by = 720)
    Cp_ere    <- numeric(length(times_ere))
    for (td in doses_h) {
      Cp_ere <- Cp_ere + ere_dose * 0.82 / 3.86 *
        exp(-0.000193 * pmax(times_ere - td, 0))
    }
    df_ere <- data.frame(time = times_ere / 24, Cp = Cp_ere, drug = "Erenumab")

    # Topiramate (steady-state approximation)
    top_half <- input$top_dose / 2  # twice daily
    times_top <- seq(0, min(t_end_h, 720), by = 0.5)
    Cp_top <- numeric(length(times_top))
    ka_t <- 0.35; CL_t <- input$top_dose / 22  # Css approx
    for (i in seq_along(times_top)) {
      t <- times_top[i]
      n_doses <- floor(t / 12)
      for (j in 0:n_doses) {
        dt <- t - j * 12
        Cp_top[i] <- Cp_top[i] +
          top_half * 0.8 / 45.5 * exp(-0.033 * dt)
      }
    }
    df_top <- data.frame(time = times_top / 24, Cp = Cp_top, drug = "Topiramate")

    plot_ly() %>%
      add_lines(data = df_ere, x = ~time, y = ~Cp, name = "Erenumab (mg/L)",
                line = list(color = "#2c7bb6", width = 2)) %>%
      add_lines(data = df_top, x = ~time, y = ~Cp,
                name = "Topiramate (μg/mL)",
                line = list(color = "#d7191c", width = 2, dash = "dash"),
                yaxis = "y2") %>%
      layout(
        xaxis = list(title = "Time (days)"),
        yaxis  = list(title = "Erenumab (mg/L)", side = "left"),
        yaxis2 = list(title = "Topiramate (μg/mL)", overlaying = "y",
                      side = "right"),
        hovermode = "x unified",
        legend = list(x = 0.6, y = 0.9)
      )
  })

  # ── PD: CGRP dynamics ────────────────────────────────────────────────────
  output$pd_cgrp_plot <- renderPlotly({
    times <- seq(0, 24, by = 0.1)

    # CGRP dynamics for each drug scenario
    sim_cgrp <- function(drug_block) {
      CGRP <- numeric(length(times))
      CGRP[1] <- 2.5
      dt <- 0.1
      for (i in seq_len(length(times)-1)) {
        dC <- input$ksyn_CGRP +
              0.3 * input$CSD_init * exp(-0.15 * times[i]) -
              input$kdeg_CGRP * CGRP[i] * (1 + 2 * drug_block)
        CGRP[i+1] <- pmax(CGRP[i] + dC * dt, 0)
      }
      CGRP
    }

    df_cgrp <- data.frame(
      time = times,
      Untreated       = sim_cgrp(0),
      Sumatriptan_SC  = sim_cgrp(0.35),
      Erenumab_140mg  = sim_cgrp(0.85),
      Rimegepant_75mg = sim_cgrp(0.65)
    ) %>%
      pivot_longer(-time, names_to = "scenario", values_to = "CGRP")

    plot_ly(df_cgrp, x = ~time, y = ~CGRP, color = ~scenario,
            type = "scatter", mode = "lines",
            line = list(width = 2)) %>%
      layout(
        xaxis = list(title = "Time (hours)"),
        yaxis = list(title = "CGRP (pmol/L)"),
        shapes = list(
          list(type = "line", x0 = 0, x1 = 24, y0 = 0.5, y1 = 0.5,
               line = list(dash = "dot", color = "grey"))
        )
      )
  })

  # ── PD: Receptor occupancy ────────────────────────────────────────────────
  output$pd_occ_plot <- renderPlotly({
    conc_nM <- 10^seq(-3, 3, by = 0.05)
    df_occ <- data.frame(
      conc_nM  = conc_nM,
      Erenumab  = conc_nM / (conc_nM + input$ere_KD),
      Rimegepant = conc_nM / (conc_nM + input$rim_Ki),
      Sumatriptan_1D = conc_nM / (conc_nM + input$sum_EC50)
    ) %>%
      pivot_longer(-conc_nM, names_to = "drug", values_to = "occupancy")

    plot_ly(df_occ, x = ~conc_nM, y = ~occupancy, color = ~drug,
            type = "scatter", mode = "lines") %>%
      layout(
        xaxis = list(title = "Drug Concentration (nM)", type = "log"),
        yaxis = list(title = "Receptor Occupancy (fraction)", range = c(0, 1))
      )
  })

  # ── PD: CSD propagation ───────────────────────────────────────────────────
  output$pd_csd_plot <- renderPlotly({
    times <- seq(0, 60, by = 0.2)
    CSD <- numeric(length(times))
    CSD[1] <- input$CSD_init
    dt <- 0.2
    for (i in seq_len(length(times)-1)) {
      dCSD <- 0.10 * CSD[i] * (1 - CSD[i]) - 0.05 * CSD[i]
      CSD[i+1] <- pmax(pmin(CSD[i] + dCSD * dt, 1), 0)
    }
    df_csd <- data.frame(time = times, CSD_activity = CSD)

    plot_ly(df_csd, x = ~time, y = ~CSD_activity,
            type = "scatter", mode = "lines",
            line = list(color = "#c0392b", width = 2)) %>%
      layout(
        xaxis = list(title = "Time (minutes)"),
        yaxis = list(title = "CSD Activity (0–1)", range = c(0, 1))
      )
  })

  # ── Endpoints: pain over time ─────────────────────────────────────────────
  output$ep_pain_plot <- renderPlotly({
    times <- seq(0, 24, by = 0.1)
    make_pain <- function(block_frac) {
      pain <- numeric(length(times))
      pain[1] <- input$pain_peak
      dt <- 0.1
      for (i in seq_len(length(times)-1)) {
        dpain <- -0.05 * pain[i] - 0.20 * block_frac * pain[i]
        pain[i+1] <- pmax(pain[i] + dpain * dt, 0)
      }
      pain
    }

    df_pain <- data.frame(
      time = times,
      Untreated       = make_pain(0),
      Sumatriptan_SC  = make_pain(0.70),
      Lasmiditan      = make_pain(0.60),
      Rimegepant      = make_pain(0.50),
      Ubrogepant      = make_pain(0.45)
    ) %>%
      pivot_longer(-time, names_to = "treatment", values_to = "VAS")

    plot_ly(df_pain, x = ~time, y = ~VAS, color = ~treatment,
            type = "scatter", mode = "lines") %>%
      add_segments(x = 0, xend = 24, y = 0.5, yend = 0.5,
                   name = "Pain freedom threshold",
                   line = list(dash = "dot", color = "grey")) %>%
      layout(
        xaxis = list(title = "Time (hours)"),
        yaxis = list(title = "VAS Pain Score (0–10)", range = c(0, 10))
      )
  })

  # ── Endpoints: MMD reduction ──────────────────────────────────────────────
  output$ep_mmd_plot <- renderPlotly({
    months <- 1:6
    base   <- input$baseline_MMD
    df_mmd <- data.frame(
      month = months,
      "No Treatment"           = rep(base, 6),
      "Erenumab 140 mg"        = pmax(base * (1 - 0.50 * (1 - exp(-months/2))), 0),
      "Fremanezumab 225 mg"    = pmax(base * (1 - 0.45 * (1 - exp(-months/2))), 0),
      "Topiramate 100 mg"      = pmax(base * (1 - 0.40 * (1 - exp(-months/3))), 0),
      "Propranolol 160 mg"     = pmax(base * (1 - 0.38 * (1 - exp(-months/3))), 0),
      check.names = FALSE
    ) %>%
      pivot_longer(-month, names_to = "treatment", values_to = "MMD")

    plot_ly(df_mmd, x = ~month, y = ~MMD, color = ~treatment,
            type = "scatter", mode = "lines+markers") %>%
      add_segments(x = 0, xend = 7, y = 15, yend = 15,
                   name = "Chronic threshold (15 MHD)",
                   line = list(dash = "dot", color = "red")) %>%
      layout(
        xaxis = list(title = "Month of Treatment", dtick = 1),
        yaxis = list(title = "Monthly Migraine Days")
      )
  })

  # ── Endpoints: responder rates ────────────────────────────────────────────
  output$ep_responder <- renderPlotly({
    drugs <- c("Erenumab 70mg", "Erenumab 140mg", "Fremanezumab 225mg",
               "Galcanezumab 120mg", "Rimegepant 75mg QOD",
               "Topiramate 100mg", "Propranolol 160mg")
    rates_50 <- c(0.40, 0.47, 0.43, 0.52, 0.28, 0.37, 0.35)
    rates_75 <- c(0.22, 0.28, 0.26, 0.31, 0.14, 0.20, 0.18)

    plot_ly() %>%
      add_bars(x = drugs, y = rates_50 * 100,
               name = "≥50% Responder", marker = list(color = "#2c7bb6")) %>%
      add_bars(x = drugs, y = rates_75 * 100,
               name = "≥75% Responder", marker = list(color = "#abd9e9")) %>%
      layout(
        barmode = "group",
        xaxis = list(title = "", tickangle = -30),
        yaxis = list(title = "Responder Rate (%)", range = c(0, 60))
      )
  })

  # ── Trial benchmark table ─────────────────────────────────────────────────
  output$trial_table <- renderDT({
    data.frame(
      Trial = c("ARISE", "STRIVE", "HALO-EM monthly", "EVOLVE-1",
                "SAMURAI", "SPARTAN", "ARTISAN-EM", "ACHIEVE-I"),
      Drug = c("Erenumab 70mg", "Erenumab 140mg", "Fremanezumab 225mg",
               "Galcanezumab 120mg", "Lasmiditan 200mg", "Lasmiditan 200mg",
               "Rimegepant 75mg", "Ubrogepant 100mg"),
      "Primary Endpoint" = c(
        "MMD reduction -2.9", "MMD reduction -3.7", "MMD reduction -3.7",
        "MMD reduction -4.7", "Pain free 2h 32%", "Pain free 2h 39%",
        "Pain free 2h 21%", "Pain free 2h 19%"
      ),
      "50% Responder" = c("40%", "47%", "43%", "52%", "N/A", "N/A",
                           "28% (prevention)", "N/A"),
      Year = c(2017, 2017, 2017, 2018, 2019, 2019, 2021, 2019),
      stringsAsFactors = FALSE
    )
  }, options = list(pageLength = 8, scrollX = TRUE), rownames = FALSE)

  # ── Scenario comparison ───────────────────────────────────────────────────
  output$compare_plot <- renderPlotly({
    times <- seq(0, 24, by = 0.1)
    scenarios <- list(
      "Untreated"             = list(block = 0.0,  label_col = "#e63946"),
      "Sumatriptan SC 6 mg"   = list(block = 0.70, label_col = "#2c7bb6"),
      "Sumatriptan Oral 100 mg" = list(block = 0.55, label_col = "#457b9d"),
      "Lasmiditan 200 mg"     = list(block = 0.60, label_col = "#1d3557"),
      "Rimegepant 75 mg"      = list(block = 0.50, label_col = "#2a9d8f"),
      "Ubrogepant 100 mg"     = list(block = 0.45, label_col = "#264653"),
      "NSAID (naproxen 500 mg)" = list(block = 0.35, label_col = "#e9c46a")
    )

    selected <- input$acute_compare
    plt <- plot_ly()

    for (nm in selected) {
      if (!nm %in% names(scenarios)) next
      blk <- scenarios[[nm]]$block
      col <- scenarios[[nm]]$label_col

      if (input$compare_outcome == "pain") {
        vals <- sapply(seq_along(times), function(i) {
          pmax(input$pain_peak * exp(-(0.05 + 0.20 * blk) * times[i]), 0)
        })
        ylab <- "VAS Pain Score"
      } else if (input$compare_outcome == "cgrp") {
        vals <- pmax(2.5 * exp(-(0.10 + 0.15 * blk) * times), 0)
        ylab <- "CGRP (pmol/L)"
      } else if (input$compare_outcome == "tg") {
        vals <- pmax(0.8 * exp(-(0.05 + 0.25 * blk) * times), 0)
        ylab <- "TG Activation"
      } else {
        vals <- rep(blk * 100, length(times))
        ylab <- "CGRP-R Occupancy (%)"
      }

      plt <- plt %>% add_lines(x = times, y = vals, name = nm,
                                line = list(color = col, width = 2))
    }

    plt %>% layout(
      xaxis = list(title = "Time (hours)"),
      yaxis = list(title = ylab),
      hovermode = "x unified"
    )
  })

  # ── Preventive comparison ─────────────────────────────────────────────────
  output$prev_compare_plot <- renderPlotly({
    months   <- seq(0, input$prev_sim_months, by = 0.1)
    base_MMD <- input$baseline_MMD

    prev_params <- list(
      "No treatment"                = list(eff = 0.00, onset = 1),
      "Erenumab 70 mg"              = list(eff = 0.40, onset = 1.5),
      "Erenumab 140 mg"             = list(eff = 0.50, onset = 1.5),
      "Fremanezumab 225 mg monthly" = list(eff = 0.45, onset = 1.5),
      "Galcanezumab 120 mg"         = list(eff = 0.52, onset = 1.5),
      "Topiramate 100 mg/day"       = list(eff = 0.40, onset = 2.5),
      "Propranolol 160 mg/day"      = list(eff = 0.38, onset = 3.0),
      "Amitriptyline 50 mg/day"     = list(eff = 0.35, onset = 2.5)
    )

    selected <- input$prev_compare
    plt <- plot_ly()
    colors <- c("#e63946","#2c7bb6","#457b9d","#2a9d8f","#264653","#e9c46a","#f4a261","#a8dadc")

    for (i in seq_along(selected)) {
      nm <- selected[i]
      if (!nm %in% names(prev_params)) next
      p  <- prev_params[[nm]]
      mmd_vals <- base_MMD * (1 - p$eff * (1 - exp(-months / p$onset)))
      plt <- plt %>%
        add_lines(x = months, y = mmd_vals, name = nm,
                  line = list(color = colors[i %% length(colors) + 1], width = 2))
    }

    plt %>%
      add_segments(x = 0, xend = input$prev_sim_months, y = 15, yend = 15,
                   name = "Chronic threshold", inherit = FALSE,
                   line = list(dash = "dot", color = "red")) %>%
      layout(
        xaxis = list(title = "Month"),
        yaxis = list(title = "Monthly Migraine Days (MMD)"),
        hovermode = "x unified"
      )
  })

  # ── Biomarker: CGRP trajectory ────────────────────────────────────────────
  output$bio_cgrp_plot <- renderPlotly({
    phases <- c("Prodrome", "Aura", "Headache onset", "Headache peak",
                "Resolution", "Postdrome", "Interictal")
    cgrp_ictal <- c(60, 80, 150, 200, 130, 90, 55)  # pg/mL (approximate)

    plot_ly(x = phases, y = cgrp_ictal, type = "scatter", mode = "lines+markers",
            line = list(color = "#c0392b", width = 2),
            marker = list(size = 8)) %>%
      add_segments(x = "Prodrome", xend = "Interictal",
                   y = 50, yend = 50, name = "Normal range upper",
                   line = list(dash = "dot", color = "grey")) %>%
      layout(
        xaxis = list(title = "Migraine Phase"),
        yaxis = list(title = "Plasma CGRP (pg/mL)")
      )
  })

  # ── Biomarker: 5-HT platelet ──────────────────────────────────────────────
  output$bio_5ht_plot <- renderPlotly({
    phases <- c("Interictal", "Prodrome", "Aura", "Headache", "Resolution", "Postdrome")
    ht5_levels <- c(0.90, 0.85, 0.70, 0.45, 0.55, 0.75)

    plot_ly(x = phases, y = ht5_levels * 100, type = "bar",
            marker = list(color = "#1b4f72")) %>%
      layout(
        xaxis = list(title = "Migraine Phase"),
        yaxis = list(title = "Relative Platelet 5-HT (% of interictal)", range = c(0, 110))
      )
  })

  # ── Biomarker: PGE2 ──────────────────────────────────────────────────────
  output$bio_pge2_plot <- renderPlotly({
    times <- seq(0, 24, by = 0.5)
    pge2_vals <- 20 + 60 * exp(-0.08 * times)
    il6_vals  <- 5 + 15 * exp(-0.12 * times)

    plot_ly() %>%
      add_lines(x = times, y = pge2_vals, name = "PGE2 (pg/mL)",
                line = list(color = "#e74c3c", width = 2)) %>%
      add_lines(x = times, y = il6_vals, name = "IL-6 (pg/mL)",
                line = list(color = "#e67e22", width = 2, dash = "dash"),
                yaxis = "y2") %>%
      layout(
        xaxis = list(title = "Time from attack onset (h)"),
        yaxis  = list(title = "PGE2 (pg/mL)", side = "left"),
        yaxis2 = list(title = "IL-6 (pg/mL)", overlaying = "y", side = "right")
      )
  })

  # ── Biomarker reference table ──────────────────────────────────────────────
  output$bio_ref_table <- renderDT({
    data.frame(
      Biomarker = c("Plasma CGRP", "Jugular venous CGRP", "Platelet 5-HT",
                    "PGE2 (CSF)", "IL-1β (plasma)", "TNF-α (plasma)",
                    "NO metabolites", "BDNF (serum)"),
      Interictal = c("~40–70 pg/mL", "~35–60 pg/mL", "~0.8–1.2 ng/mL",
                     "<1 pg/mL", "<2 pg/mL", "<10 pg/mL",
                     "~15 μM", "~15 ng/mL"),
      Ictal = c("↑ 100–400 pg/mL", "↑ 100–350 pg/mL", "↓ 0.3–0.6 ng/mL",
                "↑ 2–5 pg/mL", "↑ 4–8 pg/mL", "↑ 15–40 pg/mL",
                "↑ 20–30 μM", "↑ 20–30 ng/mL"),
      Clinical_significance = c("Primary CGRP mAb target", "Validates CGRP release",
                                 "5-HT deficiency hypothesis", "Neuroinflammation",
                                 "Central sensitization", "Mast cell activation",
                                 "NTG model / NO pathway", "Neurotrophin sensitization"),
      stringsAsFactors = FALSE
    )
  }, options = list(pageLength = 8), rownames = FALSE)

}

# ── Run ────────────────────────────────────────────────────────────────────
shinyApp(ui, server)
