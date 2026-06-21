## ============================================================
## Obstructive Sleep Apnea (OSA) QSP — Interactive Shiny App
## 폐쇄성 수면 무호흡증 정량적 시스템 약리학 대시보드
##
## Tabs:
##   1. Patient Profile (환자 프로파일)
##   2. PK Profiles (약동학 프로파일)
##   3. Airway Mechanics (기도 역학)
##   4. Cardiovascular Effects (심혈관 효과)
##   5. Metabolic & Inflammatory (대사·염증 지표)
##   6. Clinical Endpoints (임상 종말점)
##   7. Scenario Comparison (시나리오 비교)
##   8. Biomarker Heatmap (바이오마커 히트맵)
## ============================================================

library(shiny)
library(shinydashboard)
library(ggplot2)
library(dplyr)
library(tidyr)
library(scales)
library(plotly)
library(DT)

## ── Embedded ODE solver (no mrgsolve required for Shiny) ──────
## Uses simple Euler integration for deployment portability

osa_ode_euler <- function(params, init, times) {
  dt <- times[2] - times[1]
  n  <- length(times)

  # Unpack parameters
  Pcrit_base     <- params$Pcrit_base
  Pcrit_BMI_coef <- params$Pcrit_BMI_coef
  LG_base        <- params$LG_base
  AT_base        <- params$AT_base
  AHI_base       <- params$AHI_base
  SpO2_base      <- params$SpO2_base
  SpO2_max       <- 98.5
  CPAP_on        <- params$CPAP_on
  CPAP_pressure  <- params$CPAP_pressure
  MODA_on        <- params$MODA_on
  SOLRI_on       <- params$SOLRI_on
  TIRZ_on        <- params$TIRZ_on
  ESZOP_on       <- params$ESZOP_on
  ACETZ_on       <- params$ACETZ_on

  # State vector: Pcrit, LG, AT, AHI, SpO2, SNA, SBP, HR, CRP, HOMA, ESS, Weight
  state <- unlist(init)
  names(state) <- names(init)

  out <- matrix(NA_real_, nrow = n, ncol = length(state) + 1)
  colnames(out) <- c("time", names(state))
  out[1, ] <- c(times[1], state)

  for (i in 2:n) {
    Pcrit  <- state["Pcrit"]
    LG     <- state["LG"]
    AT     <- state["AT"]
    AHI    <- state["AHI"]
    SpO2   <- state["SpO2"]
    SNA    <- state["SNA"]
    SBP    <- state["SBP"]
    HR     <- state["HR"]
    CRP    <- state["CRP"]
    HOMA   <- state["HOMA"]
    ESS    <- state["ESS"]
    Weight <- state["Weight"]

    BMI <- Weight / 1.75^2

    # ---- Drug PD effects (steady-state approximation) ----
    CPAP_Pcrit_red <- if (CPAP_on) (-CPAP_pressure - 1.5) else 0
    Tirz_wt_eff    <- if (TIRZ_on) pmin(0.20, 0.20 * (1 - exp(-times[i] / (52 * 168 * 0.5)))) else 0
    Tirz_Pcrit_eff <- -Tirz_wt_eff * 0.10 * (BMI - 25)
    ESZOP_AT_eff   <- if (ESZOP_on) 0.45 else 0
    ACETZ_LG_eff   <- if (ACETZ_on) 0.35 else 0
    MODA_ESS_eff   <- if (MODA_on) 0.60 else 0
    SOLRI_ESS_eff  <- if (SOLRI_on) 0.65 else 0

    # ---- Pcrit ----
    Pcrit_target <- Pcrit_base + Pcrit_BMI_coef * (BMI - 25) +
                    CPAP_Pcrit_red + Tirz_Pcrit_eff
    dPcrit <- 0.02 * (Pcrit_target - Pcrit)

    # ---- LG ----
    LG_hyp_eff   <- 0.1 * max(0, 95 - SpO2) / 10
    LG_target    <- (LG_base + LG_hyp_eff) * (1 - ACETZ_LG_eff)
    dLG          <- 0.05 * (LG_target - LG)

    # ---- AT ----
    AT_target <- AT_base * (1 + ESZOP_AT_eff)
    dAT       <- 0.05 * (AT_target - AT)

    # ---- AHI ----
    AHI_t_raw <- 4.5 * (Pcrit - (-2)) + 25 * (LG - 0.70) - 0.8 * (AT - 14)
    AHI_target <- max(1, AHI_t_raw)
    dAHI       <- 0.3 * (AHI_target - AHI)

    # ---- SpO2 ----
    SpO2_target <- SpO2_base + (-0.2) * AHI +
                   ifelse(CPAP_on, (SpO2_max - SpO2_base) * 0.9, 0)
    SpO2_target <- max(70, min(SpO2_max, SpO2_target))
    dSpO2       <- 2.0 * (SpO2_target - SpO2)

    # ---- SNA ----
    SNA_target <- 1.0 + 0.6 * (AHI / AHI_base - 1)
    SNA_target <- max(0.2, SNA_target)
    dSNA       <- 0.8 * (SNA_target - SNA)

    # ---- SBP ----
    SBP_target <- 138 + 12 * (SNA - 1)
    dSBP       <- 0.3 * (SBP_target - SBP)

    # ---- HR ----
    HR_target <- 72 + 10 * (SNA - 1)
    dHR       <- 1.0 * (HR_target - HR)

    # ---- CRP ----
    CRP_target <- 3.5 + 0.05 * AHI + 0.02 * (SBP - 130)
    dCRP       <- 0.02 * (CRP_target - CRP)

    # ---- HOMA ----
    HOMA_target <- max(1, 3.2 + 0.15 * (CRP - 1) + 0.08 * (SNA - 1) -
                         Tirz_wt_eff * 0.08)
    dHOMA       <- 0.01 * (HOMA_target - HOMA)

    # ---- ESS ----
    drug_eff    <- min(0.8, MODA_ESS_eff + SOLRI_ESS_eff)
    ESS_target  <- pmax(0, pmin(24, ESS * (1 - drug_eff) +
                                 0.08 * (AHI - AHI_base) +
                                 (-0.15) * (SpO2 - SpO2_base)))
    dESS        <- 0.05 * (ESS_target - ESS)

    # ---- Weight ----
    dWeight <- -0.001 * Tirz_wt_eff * Weight

    state["Pcrit"]  <- Pcrit + dt * dPcrit
    state["LG"]     <- max(0.1, LG  + dt * dLG)
    state["AT"]     <- AT     + dt * dAT
    state["AHI"]    <- max(1,   AHI + dt * dAHI)
    state["SpO2"]   <- max(70,  min(99, SpO2 + dt * dSpO2))
    state["SNA"]    <- max(0.2, SNA  + dt * dSNA)
    state["SBP"]    <- SBP    + dt * dSBP
    state["HR"]     <- HR     + dt * dHR
    state["CRP"]    <- max(0.1, CRP  + dt * dCRP)
    state["HOMA"]   <- HOMA   + dt * dHOMA
    state["ESS"]    <- max(0, min(24, ESS + dt * dESS))
    state["Weight"] <- max(50, Weight + dt * dWeight)

    out[i, ] <- c(times[i], state)
  }

  as.data.frame(out) %>%
    mutate(week = time / 168,
           BMI  = Weight / 1.75^2,
           SpO2_class = case_when(
             SpO2 >= 95 ~ "Normal (≥95%)",
             SpO2 >= 90 ~ "Mild (90-94%)",
             SpO2 >= 85 ~ "Moderate (85-89%)",
             TRUE       ~ "Severe (<85%)"
           ),
           AHI_class = case_when(
             AHI < 5   ~ "Normal",
             AHI < 15  ~ "Mild",
             AHI < 30  ~ "Moderate",
             TRUE       ~ "Severe"
           ))
}

## ── PK profile (1-compartment oral, Bateman equation) ─────────
pk_profile <- function(dose, CL, V, ka, F_oral = 1, hours = 48, label = "") {
  t  <- seq(0, hours, by = 0.25)
  kel <- CL / V
  Cp <- (dose * F_oral * ka) / (V * (ka - kel)) *
        (exp(-kel * t) - exp(-ka * t))
  Cp <- pmax(0, Cp)
  tibble(time_h = t, Cp = Cp, drug = label)
}

## ══════════════════════════════════════════════════════════════
## UI
## ══════════════════════════════════════════════════════════════
ui <- dashboardPage(
  skin = "blue",

  dashboardHeader(
    title = span(icon("lungs"), "OSA QSP Model",
                 style = "font-size:16px; font-weight:bold;"),
    titleWidth = 260
  ),

  dashboardSidebar(
    width = 260,
    sidebarMenu(
      id = "tabs",
      menuItem("Patient Profile",   tabName = "tab_profile",  icon = icon("user-md")),
      menuItem("PK Profiles",       tabName = "tab_pk",       icon = icon("pills")),
      menuItem("Airway Mechanics",  tabName = "tab_airway",   icon = icon("wind")),
      menuItem("Cardiovascular",    tabName = "tab_cv",       icon = icon("heartbeat")),
      menuItem("Metabolic Effects", tabName = "tab_metabolic",icon = icon("dna")),
      menuItem("Clinical Endpoints",tabName = "tab_clinical", icon = icon("chart-line")),
      menuItem("Scenario Comparison",tabName ="tab_compare",  icon = icon("layer-group")),
      menuItem("Biomarker Heatmap", tabName = "tab_heatmap",  icon = icon("th"))
    ),
    hr(),
    h5("  Patient Parameters", style = "color:#ccc; padding-left:12px;"),
    sliderInput("bmi",   "BMI (kg/m²)",    min = 20, max = 55, value = 33, step = 0.5),
    sliderInput("ahi_init", "Initial AHI (ev/hr)", min = 5, max = 80, value = 35),
    sliderInput("age",   "Age (years)",     min = 20, max = 80, value = 48),
    selectInput("sex",   "Sex",             choices = c("Male","Female"), selected = "Male"),
    hr(),
    h5("  Treatment", style = "color:#ccc; padding-left:12px;"),
    checkboxInput("cpap_on",  "CPAP", FALSE),
    conditionalPanel("input.cpap_on",
      sliderInput("cpap_p", "CPAP Pressure (cmH2O)", 4, 20, 10, step = 0.5)
    ),
    checkboxInput("tirz_on",  "Tirzepatide", FALSE),
    checkboxInput("moda_on",  "Modafinil", FALSE),
    checkboxInput("solri_on", "Solriamfetol", FALSE),
    checkboxInput("eszop_on", "Eszopiclone", FALSE),
    checkboxInput("acetz_on", "Acetazolamide", FALSE),
    hr(),
    sliderInput("sim_weeks", "Simulation (weeks)", 4, 104, 52, step = 4),
    actionButton("run_btn", "Run Simulation", icon = icon("play"),
                 class = "btn-success", style = "margin:10px 12px; width:220px;")
  ),

  dashboardBody(
    tags$head(tags$style(HTML("
      .content-wrapper { background-color: #f8f9fa; }
      .box { border-top: 3px solid #3c8dbc; }
      .box-title { font-weight: bold; }
      .ahi-badge { padding: 4px 10px; border-radius: 12px;
                   color: white; font-weight: bold; font-size: 13px; }
    "))),

    tabItems(

      ## ─── TAB 1: Patient Profile ──────────────────────────
      tabItem(tabName = "tab_profile",
        fluidRow(
          valueBoxOutput("vbox_ahi",    width = 3),
          valueBoxOutput("vbox_spo2",   width = 3),
          valueBoxOutput("vbox_sbp",    width = 3),
          valueBoxOutput("vbox_ess",    width = 3)
        ),
        fluidRow(
          box(title = "Patient Overview", status = "primary", solidHeader = TRUE,
              width = 4,
              tableOutput("patient_table")
          ),
          box(title = "OSA Endotype Classification", status = "warning",
              solidHeader = TRUE, width = 4,
              plotOutput("endotype_radar", height = "280px"),
              p("Endotypes based on Eckert et al. Sleep 2013", style = "font-size:11px; color:#666;")
          ),
          box(title = "OSA Severity Classification", status = "danger",
              solidHeader = TRUE, width = 4,
              plotOutput("severity_gauge", height = "280px")
          )
        ),
        fluidRow(
          box(title = "Baseline Biomarker Profile", status = "info",
              solidHeader = TRUE, width = 12,
              plotlyOutput("baseline_bar", height = "260px"))
        )
      ),

      ## ─── TAB 2: PK Profiles ──────────────────────────────
      tabItem(tabName = "tab_pk",
        fluidRow(
          box(title = "Wake-Promoting Agent PK (Single Dose)", status = "primary",
              solidHeader = TRUE, width = 6,
              plotlyOutput("pk_wake", height = "320px")
          ),
          box(title = "Arousal/Loop-Gain Modulator PK", status = "warning",
              solidHeader = TRUE, width = 6,
              plotlyOutput("pk_arousal", height = "320px")
          )
        ),
        fluidRow(
          box(title = "Tirzepatide PK (Multi-dose, 10mg QW sc)", status = "success",
              solidHeader = TRUE, width = 6,
              plotlyOutput("pk_tirz", height = "320px")
          ),
          box(title = "PK Parameter Summary", status = "info",
              solidHeader = TRUE, width = 6,
              DTOutput("pk_table"))
        )
      ),

      ## ─── TAB 3: Airway Mechanics ─────────────────────────
      tabItem(tabName = "tab_airway",
        fluidRow(
          box(title = "Critical Closing Pressure (Pcrit) Over Time",
              status = "danger", solidHeader = TRUE, width = 7,
              plotlyOutput("pcrit_plot", height = "350px"),
              p("Pcrit > 0 cmH2O = passive collapse; Pcrit 0 to -2 = partial obstruction;
                 Pcrit < -2 = normal", style = "font-size:11px; color:#555;")
          ),
          box(title = "Loop Gain & Arousal Threshold", status = "warning",
              solidHeader = TRUE, width = 5,
              plotlyOutput("lg_at_plot", height = "350px")
          )
        ),
        fluidRow(
          box(title = "AHI Trajectory", status = "primary", solidHeader = TRUE, width = 8,
              plotlyOutput("ahi_time_plot", height = "320px")
          ),
          box(title = "Airway Mechanics Key Concepts", status = "info",
              solidHeader = TRUE, width = 4,
              h5("Starling Resistor Model"),
              p("Airway collapse occurs when intraluminal pressure falls below
                the critical closing pressure (Pcrit)."),
              h5("Four OSA Endotypes"),
              tags$ul(
                tags$li("High pharyngeal collapsibility (Pcrit > 0)"),
                tags$li("Elevated loop gain (ventilatory instability)"),
                tags$li("Low arousal threshold"),
                tags$li("Inadequate muscle responsiveness")
              ),
              h5("CPAP Mechanism"),
              p("Positive airway pressure stents the pharynx, lowering effective
                Pcrit to below atmospheric pressure, eliminating obstructions.")
          )
        )
      ),

      ## ─── TAB 4: Cardiovascular ───────────────────────────
      tabItem(tabName = "tab_cv",
        fluidRow(
          box(title = "Blood Pressure (SBP/DBP) Over Time", status = "danger",
              solidHeader = TRUE, width = 6,
              plotlyOutput("bp_plot", height = "320px")
          ),
          box(title = "Heart Rate & SNS Activity", status = "warning",
              solidHeader = TRUE, width = 6,
              plotlyOutput("hr_sna_plot", height = "320px")
          )
        ),
        fluidRow(
          box(title = "Cardiovascular Risk Pathway Activation",
              status = "primary", solidHeader = TRUE, width = 12,
              plotlyOutput("cv_risk_plot", height = "280px"),
              p("Pathway: AHI → Intermittent hypoxia → SNS activation → Hypertension → LVH → MACE",
                style = "font-size:11px; color:#555; margin-top:5px;")
          )
        )
      ),

      ## ─── TAB 5: Metabolic Effects ────────────────────────
      tabItem(tabName = "tab_metabolic",
        fluidRow(
          box(title = "Insulin Resistance (HOMA-IR)", status = "warning",
              solidHeader = TRUE, width = 6,
              plotlyOutput("homa_plot", height = "300px")
          ),
          box(title = "Inflammatory Marker (hsCRP)", status = "danger",
              solidHeader = TRUE, width = 6,
              plotlyOutput("crp_plot", height = "300px")
          )
        ),
        fluidRow(
          box(title = "Weight / BMI Trajectory (Tirzepatide)", status = "success",
              solidHeader = TRUE, width = 6,
              plotlyOutput("weight_plot", height = "300px")
          ),
          box(title = "Metabolic Pathway Summary", status = "info",
              solidHeader = TRUE, width = 6,
              h5("Intermittent Hypoxia → Metabolic Dysregulation"),
              tags$ul(
                tags$li("HIF-1α → lipolysis → ↑ FFA → ↑ TG, LDL"),
                tags$li("SNS → HPA axis → ↑ cortisol → insulin resistance"),
                tags$li("Adipose inflammation → TNF-α, IL-6 → ↑ CRP"),
                tags$li("Leptin resistance → appetite dysregulation"),
                tags$li("Ghrelin ↑ with sleep fragmentation → obesity")
              ),
              h5("SURMOUNT-OSA Trial (Wolk et al. NEJM 2024)"),
              p("Tirzepatide 10-15mg QW reduced AHI by ~55 events/hr in
                non-CPAP group (BMI ~38 kg/m²) over 52 weeks.")
          )
        )
      ),

      ## ─── TAB 6: Clinical Endpoints ───────────────────────
      tabItem(tabName = "tab_clinical",
        fluidRow(
          box(title = "Epworth Sleepiness Scale (ESS) Over Time",
              status = "primary", solidHeader = TRUE, width = 6,
              plotlyOutput("ess_plot", height = "320px")
          ),
          box(title = "SpO2 Dynamics", status = "info",
              solidHeader = TRUE, width = 6,
              plotlyOutput("spo2_plot", height = "320px")
          )
        ),
        fluidRow(
          box(title = "Multi-Endpoint Trajectory (Normalized)", status = "success",
              solidHeader = TRUE, width = 8,
              plotlyOutput("multi_endpoint", height = "320px")
          ),
          box(title = "Endpoint Legend & Thresholds", status = "warning",
              solidHeader = TRUE, width = 4,
              h5("Clinical Thresholds"),
              tableOutput("threshold_table")
          )
        )
      ),

      ## ─── TAB 7: Scenario Comparison ──────────────────────
      tabItem(tabName = "tab_compare",
        fluidRow(
          box(title = "Predefined Scenario Comparison: AHI", status = "primary",
              solidHeader = TRUE, width = 12,
              plotlyOutput("compare_ahi", height = "380px")
          )
        ),
        fluidRow(
          box(title = "Scenario: ESS Improvement", status = "warning",
              solidHeader = TRUE, width = 6,
              plotlyOutput("compare_ess", height = "320px")
          ),
          box(title = "Scenario: SBP Effect", status = "danger",
              solidHeader = TRUE, width = 6,
              plotlyOutput("compare_sbp", height = "320px")
          )
        )
      ),

      ## ─── TAB 8: Biomarker Heatmap ────────────────────────
      tabItem(tabName = "tab_heatmap",
        fluidRow(
          box(title = "Biomarker Trajectory Heatmap (Normalized)",
              status = "primary", solidHeader = TRUE, width = 12,
              plotlyOutput("bm_heatmap", height = "450px"),
              p("Values normalized to baseline (1.0). Color scale: blue = improved, red = worsened.",
                style = "font-size:11px; color:#555;")
          )
        ),
        fluidRow(
          box(title = "Simulation Data Export", status = "info",
              solidHeader = TRUE, width = 12,
              DTOutput("export_table"),
              downloadButton("dl_csv", "Download CSV", class = "btn-primary")
          )
        )
      )
    )
  )
)

## ══════════════════════════════════════════════════════════════
## SERVER
## ══════════════════════════════════════════════════════════════
server <- function(input, output, session) {

  ## ── Reactive simulation ─────────────────────────────────────
  sim_data <- eventReactive(input$run_btn, {
    withProgress(message = "Running OSA simulation...", value = 0.5, {
      times <- seq(0, input$sim_weeks * 168, by = 2)

      # BMI → Pcrit offset
      pcrit_baseline <- -0.5 + 0.25 * (input$bmi - 25)  # adjust for BMI
      sex_pcrit <- if (input$sex == "Male") 0.5 else 0
      pcrit_baseline <- pcrit_baseline + sex_pcrit

      params <- list(
        Pcrit_base     = pcrit_baseline,
        Pcrit_BMI_coef = 0.25,
        LG_base        = 0.80,
        AT_base        = 14.0,
        AHI_base       = input$ahi_init,
        SpO2_base      = pmax(82, 98 - input$ahi_init * 0.17),
        CPAP_on        = as.integer(input$cpap_on),
        CPAP_pressure  = if (input$cpap_on) input$cpap_p else 10,
        MODA_on        = as.integer(input$moda_on),
        SOLRI_on       = as.integer(input$solri_on),
        TIRZ_on        = as.integer(input$tirz_on),
        ESZOP_on       = as.integer(input$eszop_on),
        ACETZ_on       = as.integer(input$acetz_on)
      )
      init <- list(
        Pcrit  = pcrit_baseline,
        LG     = 0.80 + 0.05 * (input$ahi_init / 35 - 1),
        AT     = 14.0,
        AHI    = input$ahi_init,
        SpO2   = pmax(82, 98 - input$ahi_init * 0.17),
        SNA    = 1.0,
        SBP    = 128 + input$ahi_init * 0.25,
        HR     = 72,
        CRP    = 1.0 + input$ahi_init * 0.07,
        HOMA   = 1.5 + input$ahi_init * 0.05,
        ESS    = 4 + input$ahi_init * 0.22,
        Weight = input$bmi * 1.75^2
      )
      osa_ode_euler(params, init, times)
    })
  }, ignoreNULL = FALSE)

  # Multi-scenario data (fixed scenarios for comparison tab)
  multi_sim <- reactive({
    base_init <- list(
      Pcrit = 1.5, LG = 0.82, AT = 14, AHI = 35,
      SpO2 = 92.5, SNA = 1, SBP = 138, HR = 72,
      CRP = 3.5, HOMA = 3.2, ESS = 12, Weight = 100
    )
    times <- seq(0, 52 * 168, by = 4)

    scenarios <- list(
      list(name = "Untreated",               CPAP_on=0, MODA_on=0, TIRZ_on=0, SOLRI_on=0, ESZOP_on=0, ACETZ_on=0),
      list(name = "CPAP 10 cmH2O",           CPAP_on=1, MODA_on=0, TIRZ_on=0, SOLRI_on=0, ESZOP_on=0, ACETZ_on=0),
      list(name = "CPAP + Modafinil",         CPAP_on=1, MODA_on=1, TIRZ_on=0, SOLRI_on=0, ESZOP_on=0, ACETZ_on=0),
      list(name = "CPAP + Solriamfetol",      CPAP_on=1, MODA_on=0, TIRZ_on=0, SOLRI_on=1, ESZOP_on=0, ACETZ_on=0),
      list(name = "Tirzepatide (No CPAP)",    CPAP_on=0, MODA_on=0, TIRZ_on=1, SOLRI_on=0, ESZOP_on=0, ACETZ_on=0),
      list(name = "Acetazolamide+Eszopiclone",CPAP_on=0, MODA_on=0, TIRZ_on=0, SOLRI_on=0, ESZOP_on=1, ACETZ_on=1)
    )

    lapply(scenarios, function(s) {
      p <- list(Pcrit_base=1.5, Pcrit_BMI_coef=0.3, LG_base=0.82, AT_base=14,
                AHI_base=35, SpO2_base=92.5, CPAP_on=s$CPAP_on,
                CPAP_pressure=10, MODA_on=s$MODA_on, SOLRI_on=s$SOLRI_on,
                TIRZ_on=s$TIRZ_on, ESZOP_on=s$ESZOP_on, ACETZ_on=s$ACETZ_on)
      osa_ode_euler(p, base_init, times) %>% mutate(scenario = s$name)
    }) %>% bind_rows()
  })

  ## ── Value Boxes ────────────────────────────────────────────
  get_final <- function(col) {
    d <- sim_data()
    tail(d[[col]], 1)
  }

  ahi_color <- function(ahi) {
    if (ahi < 5) "green" else if (ahi < 15) "yellow" else if (ahi < 30) "orange" else "red"
  }

  output$vbox_ahi <- renderValueBox({
    ahi_v <- round(get_final("AHI"), 1)
    cat <- if (ahi_v < 5) "Normal" else if (ahi_v < 15) "Mild" else if (ahi_v < 30) "Moderate" else "Severe"
    valueBox(paste0(ahi_v, " ev/hr"), paste("AHI —", cat),
             icon = icon("wind"), color = ahi_color(ahi_v))
  })
  output$vbox_spo2 <- renderValueBox({
    v <- round(get_final("SpO2"), 1)
    clr <- if (v >= 95) "green" else if (v >= 90) "yellow" else "red"
    valueBox(paste0(v, "%"), "Mean SpO2", icon = icon("lungs"), color = clr)
  })
  output$vbox_sbp <- renderValueBox({
    v <- round(get_final("SBP"))
    clr <- if (v < 130) "green" else if (v < 140) "yellow" else "red"
    valueBox(paste0(v, " mmHg"), "Systolic BP", icon = icon("heartbeat"), color = clr)
  })
  output$vbox_ess <- renderValueBox({
    v <- round(get_final("ESS"), 1)
    clr <- if (v <= 10) "green" else if (v <= 15) "yellow" else "red"
    valueBox(v, "ESS Score", icon = icon("bed"), color = clr)
  })

  ## ── Patient Table ──────────────────────────────────────────
  output$patient_table <- renderTable({
    d <- sim_data()
    init_row <- d[1, ]
    tibble(
      Parameter = c("Age", "Sex", "BMI", "AHI (baseline)", "SpO2 (baseline)",
                    "SBP (baseline)", "hsCRP", "HOMA-IR", "ESS"),
      Value = c(input$age, input$sex, round(input$bmi, 1),
                round(init_row$AHI, 1), round(init_row$SpO2, 1),
                round(init_row$SBP, 1), round(init_row$CRP, 2),
                round(init_row$HOMA, 2), round(init_row$ESS, 1)),
      Unit = c("years", "", "kg/m²", "events/hr", "%",
               "mmHg", "mg/L", "", "0-24")
    )
  }, striped = TRUE, hover = TRUE, bordered = TRUE)

  ## ── Severity gauge (simple bar) ───────────────────────────
  output$severity_gauge <- renderPlot({
    ahi_val <- sim_data()$AHI[1]
    df <- tibble(
      Category = factor(c("Normal","Mild","Moderate","Severe"),
                        levels = c("Normal","Mild","Moderate","Severe")),
      Upper = c(5, 15, 30, 80),
      color = c("#4CAF50","#FFC107","#FF9800","#F44336")
    )
    ggplot(df, aes(x = Category, y = Upper, fill = color)) +
      geom_col(width = 0.6, alpha = 0.6) +
      geom_hline(yintercept = ahi_val, color = "black", linewidth = 2.5, linetype = "dashed") +
      annotate("text", x = 0.6, y = ahi_val + 2, label = paste("Current AHI:", round(ahi_val, 1)),
               hjust = 0, fontface = "bold", size = 4.5) +
      scale_fill_identity() +
      labs(title = "AHI Severity Bands", y = "AHI (events/hr)", x = NULL) +
      theme_minimal(base_size = 11) +
      theme(legend.position = "none")
  })

  ## ── Endotype radar ─────────────────────────────────────────
  output$endotype_radar <- renderPlot({
    d <- sim_data()
    vals <- c(
      "Collapsibility" = pmin(1, (d$Pcrit[1] + 2) / 8),
      "Loop Gain"      = pmin(1, d$LG[1]),
      "Low AT"         = 1 - pmin(1, d$AT[1] / 20),
      "Muscle deficit" = 0.5,
      "Obesity"        = pmin(1, (d$BMI[1] - 20) / 20)
    )
    df <- tibble(axis = names(vals), value = unname(vals))
    angles <- seq(0, 2*pi, length.out = nrow(df) + 1)[-(nrow(df)+1)]
    df$x <- df$value * cos(angles)
    df$y <- df$value * sin(angles)
    df2 <- rbind(df, df[1,])

    ggplot(df2, aes(x = x, y = y)) +
      geom_polygon(fill = "#3c8dbc", alpha = 0.4, color = "#3c8dbc") +
      geom_point(data = df, aes(x = x, y = y), size = 3, color = "#3c8dbc") +
      geom_text(data = df, aes(x = x * 1.25, y = y * 1.25, label = axis), size = 3.5) +
      coord_equal() + xlim(-1.5, 1.5) + ylim(-1.5, 1.5) +
      theme_void() + ggtitle("OSA Endotype Radar")
  })

  ## ── Baseline Biomarker Bar ─────────────────────────────────
  output$baseline_bar <- renderPlotly({
    d <- sim_data()[1, ]
    bm <- tibble(
      Biomarker = c("AHI\n(ev/hr)", "SpO2\n(%)", "SBP\n(mmHg)", "HR\n(bpm)",
                    "CRP\n(mg/L)", "HOMA-IR", "ESS", "BMI\n(kg/m²)"),
      Value     = c(d$AHI, d$SpO2, d$SBP, d$HR, d$CRP, d$HOMA, d$ESS, d$BMI),
      Normal    = c(5, 97, 120, 70, 1.0, 2.0, 8, 25),
      Unit      = c("ev/hr", "%", "mmHg", "bpm", "mg/L", "", "", "kg/m²")
    ) %>%
      mutate(Ratio = Value / Normal,
             Color = ifelse(Biomarker %in% c("SpO2\n(%)"), Ratio < 1, Ratio > 1))

    plot_ly(bm, x = ~Biomarker, y = ~Value, type = "bar",
            marker = list(color = ifelse(bm$Color, "#F44336", "#4CAF50"),
                          line  = list(color = "white", width = 1)),
            text = ~paste0(round(Value, 1), " ", Unit),
            hoverinfo = "text") %>%
      layout(title = "Patient Baseline Biomarker Profile",
             yaxis = list(title = "Value"),
             xaxis = list(title = ""))
  })

  ## ── PK wake agents ─────────────────────────────────────────
  output$pk_wake <- renderPlotly({
    pk_d <- bind_rows(
      pk_profile(200, 35, 35, 1.0, F_oral=0.40, hours=48, label="Modafinil 200mg (μg/mL)"),
      pk_profile(150, 12.5, 118, 1.5, F_oral=0.95, hours=48, label="Solriamfetol 150mg (μg/mL)")
    )
    plot_ly(pk_d, x = ~time_h, y = ~Cp, color = ~drug, type = "scatter", mode = "lines",
            line = list(width = 2.5)) %>%
      layout(title = "Wake-Promoting Agents: PK Profile",
             xaxis = list(title = "Time (h)"),
             yaxis = list(title = "Plasma Concentration (μg/mL)"),
             legend = list(orientation = "h"))
  })

  output$pk_arousal <- renderPlotly({
    pk_d <- bind_rows(
      pk_profile(3, 21, 80, 1.2, F_oral=0.80, hours=24, label="Eszopiclone 3mg (ng/mL)") %>%
        mutate(Cp = Cp * 1000),
      pk_profile(250, 3.5, 30, 1.0, F_oral=0.90, hours=24, label="Acetazolamide 250mg (μg/mL)")
    )
    plot_ly(pk_d, x = ~time_h, y = ~Cp, color = ~drug, type = "scatter", mode = "lines",
            line = list(width = 2.5)) %>%
      layout(title = "Arousal/LG Modulators: PK Profile",
             xaxis = list(title = "Time (h)"),
             yaxis = list(title = "Plasma Concentration"),
             legend = list(orientation = "h"))
  })

  output$pk_tirz <- renderPlotly({
    t <- seq(0, 52*7*24, by = 12)  # 52 weeks, 12-hourly
    kel <- 0.036 / 6.5
    ka  <- 0.02
    # Accumulation model (weekly dosing)
    dose_mg <- 10
    Cp_acc <- numeric(length(t))
    dose_times <- seq(0, 52*7*24 - 1, by = 7*24)
    for (td in dose_times) {
      idx <- t >= td
      tt <- t[idx] - td
      Cp_acc[idx] <- Cp_acc[idx] +
        (dose_mg * ka / (6.5 * (ka - kel))) * (exp(-kel*tt) - exp(-ka*tt)) * 1e3
    }
    plot_ly(x = t/168, y = Cp_acc, type = "scatter", mode = "lines",
            line = list(color = "#4CAF50", width = 2)) %>%
      layout(title = "Tirzepatide 10mg QW: Steady-State Accumulation",
             xaxis = list(title = "Weeks"),
             yaxis = list(title = "Plasma Concentration (ng/mL)"))
  })

  output$pk_table <- renderDT({
    tibble(
      Drug           = c("Modafinil","Solriamfetol","Eszopiclone","Acetazolamide","Tirzepatide"),
      Dose           = c("200mg QD","75-150mg QD","3mg QHS","250mg BID","10mg QW sc"),
      `F (%)` = c(40,95,80,90,"~80"),
      `t½`           = c("15h","7.1h","6h","10h","5 days"),
      `CL/F`         = c("35L/h","12.5L/h","21L/h","3.5L/h","0.036L/h"),
      Mechanism      = c("Orexin/Hist","DAT/NET inhibit","GABA-A PAM","CA inhibit","GLP-1/GIP"),
      Indication     = c("Residual EDS","Residual EDS","↑ Arousal threshold","↓ Loop gain","Weight/AHI")
    ) %>% datatable(options = list(pageLength = 10, dom = "t"), rownames = FALSE)
  })

  ## ── Airway plots ───────────────────────────────────────────
  output$pcrit_plot <- renderPlotly({
    d <- sim_data()
    plot_ly(d, x = ~week, y = ~Pcrit, type = "scatter", mode = "lines",
            line = list(color = "#F44336", width = 2)) %>%
      add_segments(x = 0, xend = max(d$week), y = 0, yend = 0,
                   line = list(dash = "dash", color = "black")) %>%
      add_segments(x = 0, xend = max(d$week), y = -2, yend = -2,
                   line = list(dash = "dot", color = "green")) %>%
      layout(title = "Effective Critical Closing Pressure (Pcrit)",
             xaxis = list(title = "Weeks"),
             yaxis = list(title = "Pcrit (cmH2O)"))
  })

  output$lg_at_plot <- renderPlotly({
    d <- sim_data()
    plot_ly(d, x = ~week) %>%
      add_trace(y = ~LG, type = "scatter", mode = "lines", name = "Loop Gain",
                line = list(color = "steelblue", width = 2)) %>%
      add_trace(y = ~AT / 20, type = "scatter", mode = "lines", name = "AT/20 (scaled)",
                line = list(color = "orange", width = 2, dash = "dash")) %>%
      add_segments(x = 0, xend = max(d$week), y = 1, yend = 1,
                   line = list(dash = "dot", color = "red"), name = "LG=1 threshold") %>%
      layout(title = "Loop Gain & Arousal Threshold",
             xaxis = list(title = "Weeks"),
             yaxis = list(title = "Value"),
             legend = list(orientation = "h"))
  })

  output$ahi_time_plot <- renderPlotly({
    d <- sim_data()
    plot_ly(d, x = ~week, y = ~AHI, type = "scatter", mode = "lines",
            line = list(color = "#FF5722", width = 2.5),
            fill = "tozeroy", fillcolor = "rgba(255,87,34,0.1)") %>%
      add_segments(x=0, xend=max(d$week), y=30, yend=30,
                   line=list(dash="dash",color="red")) %>%
      add_segments(x=0, xend=max(d$week), y=15, yend=15,
                   line=list(dash="dash",color="orange")) %>%
      add_segments(x=0, xend=max(d$week), y=5, yend=5,
                   line=list(dash="dash",color="green")) %>%
      layout(title = "AHI Trajectory",
             xaxis = list(title = "Weeks"),
             yaxis = list(title = "AHI (events/hr)"))
  })

  ## ── Cardiovascular plots ───────────────────────────────────
  output$bp_plot <- renderPlotly({
    d <- sim_data()
    plot_ly(d, x = ~week) %>%
      add_trace(y = ~SBP, type="scatter", mode="lines", name="SBP",
                line=list(color="#F44336", width=2)) %>%
      add_trace(y = ~(SBP - 40), type="scatter", mode="lines", name="DBP (approx)",
                line=list(color="#FF9800", width=2, dash="dash")) %>%
      add_segments(x=0, xend=max(d$week), y=130, yend=130,
                   line=list(dash="dot", color="black"), name="130 threshold") %>%
      layout(title="Blood Pressure Over Time",
             xaxis=list(title="Weeks"), yaxis=list(title="BP (mmHg)"),
             legend=list(orientation="h"))
  })

  output$hr_sna_plot <- renderPlotly({
    d <- sim_data()
    plot_ly(d, x = ~week) %>%
      add_trace(y = ~HR, type="scatter", mode="lines", name="Heart Rate (bpm)",
                line=list(color="steelblue", width=2)) %>%
      add_trace(y = ~(SNA * 70), type="scatter", mode="lines",
                name="SNS Activity (scaled)", yaxis="y2",
                line=list(color="purple", width=2, dash="dash")) %>%
      layout(title="Heart Rate & Sympathetic Activity",
             xaxis=list(title="Weeks"),
             yaxis=list(title="Heart Rate (bpm)"),
             yaxis2=list(title="SNS (normalized×70)", overlaying="y", side="right"),
             legend=list(orientation="h"))
  })

  output$cv_risk_plot <- renderPlotly({
    d <- sim_data() %>% filter(week %in% seq(0, input$sim_weeks, by = 2))
    plot_ly(d, x = ~week) %>%
      add_trace(y = ~(AHI/50), name="AHI/50", type="scatter", mode="lines") %>%
      add_trace(y = ~(SNA-0.5), name="SNS-0.5", type="scatter", mode="lines") %>%
      add_trace(y = ~((SBP-110)/60), name="(SBP-110)/60", type="scatter", mode="lines") %>%
      add_trace(y = ~(CRP/10), name="CRP/10", type="scatter", mode="lines") %>%
      layout(title="CV Risk Pathway (Normalized Markers)",
             xaxis=list(title="Weeks"), yaxis=list(title="Normalized Value"),
             legend=list(orientation="h"))
  })

  ## ── Metabolic plots ────────────────────────────────────────
  output$homa_plot <- renderPlotly({
    d <- sim_data()
    plot_ly(d, x=~week, y=~HOMA, type="scatter", mode="lines",
            line=list(color="darkorange", width=2)) %>%
      add_segments(x=0, xend=max(d$week), y=2.5, yend=2.5,
                   line=list(dash="dash", color="forestgreen"), name="Normal <2.5") %>%
      layout(title="HOMA-IR (Insulin Resistance)",
             xaxis=list(title="Weeks"), yaxis=list(title="HOMA-IR"))
  })

  output$crp_plot <- renderPlotly({
    d <- sim_data()
    plot_ly(d, x=~week, y=~CRP, type="scatter", mode="lines",
            line=list(color="#F44336", width=2)) %>%
      add_segments(x=0, xend=max(d$week), y=3.0, yend=3.0,
                   line=list(dash="dash", color="orange"), name="Elevated >3") %>%
      layout(title="hsCRP (mg/L)",
             xaxis=list(title="Weeks"), yaxis=list(title="hsCRP (mg/L)"))
  })

  output$weight_plot <- renderPlotly({
    d <- sim_data()
    plot_ly(d, x=~week) %>%
      add_trace(y=~Weight, type="scatter", mode="lines", name="Weight (kg)",
                line=list(color="steelblue", width=2)) %>%
      add_trace(y=~BMI, type="scatter", mode="lines", name="BMI", yaxis="y2",
                line=list(color="tomato", width=2, dash="dash")) %>%
      layout(title="Weight / BMI Trajectory",
             xaxis=list(title="Weeks"),
             yaxis=list(title="Weight (kg)"),
             yaxis2=list(title="BMI (kg/m²)", overlaying="y", side="right"),
             legend=list(orientation="h"))
  })

  ## ── Clinical Endpoint plots ────────────────────────────────
  output$ess_plot <- renderPlotly({
    d <- sim_data()
    plot_ly(d, x=~week, y=~ESS, type="scatter", mode="lines",
            line=list(color="#9C27B0", width=2),
            fill="tozeroy", fillcolor="rgba(156,39,176,0.1)") %>%
      add_segments(x=0, xend=max(d$week), y=10, yend=10,
                   line=list(dash="dash", color="orange"), name="EDS threshold") %>%
      layout(title="Epworth Sleepiness Scale",
             xaxis=list(title="Weeks"), yaxis=list(title="ESS (0-24)", range=c(0,24)))
  })

  output$spo2_plot <- renderPlotly({
    d <- sim_data()
    plot_ly(d, x=~week, y=~SpO2, type="scatter", mode="lines",
            line=list(color="#2196F3", width=2)) %>%
      add_segments(x=0, xend=max(d$week), y=95, yend=95,
                   line=list(dash="dash", color="green"), name="Normal 95%") %>%
      layout(title="Mean Overnight SpO2",
             xaxis=list(title="Weeks"),
             yaxis=list(title="SpO2 (%)", range=c(75,100)))
  })

  output$multi_endpoint <- renderPlotly({
    d <- sim_data()
    # Normalize each to baseline
    d_norm <- d %>%
      mutate(
        AHI_n  = AHI / AHI[1],
        SpO2_n = SpO2 / SpO2[1],
        SBP_n  = SBP / SBP[1],
        ESS_n  = ESS / ESS[1],
        HOMA_n = HOMA / HOMA[1],
        CRP_n  = CRP / CRP[1]
      )
    plot_ly(d_norm, x=~week) %>%
      add_trace(y=~AHI_n,  name="AHI",  type="scatter", mode="lines") %>%
      add_trace(y=~SpO2_n, name="SpO2", type="scatter", mode="lines") %>%
      add_trace(y=~SBP_n,  name="SBP",  type="scatter", mode="lines") %>%
      add_trace(y=~ESS_n,  name="ESS",  type="scatter", mode="lines") %>%
      add_trace(y=~HOMA_n, name="HOMA-IR", type="scatter", mode="lines") %>%
      add_trace(y=~CRP_n,  name="CRP",  type="scatter", mode="lines") %>%
      add_segments(x=0, xend=max(d$week), y=1, yend=1,
                   line=list(dash="dot", color="black"), name="Baseline") %>%
      layout(title="Multi-Endpoint Trajectory (Normalized to Baseline=1)",
             xaxis=list(title="Weeks"), yaxis=list(title="Ratio to Baseline"),
             legend=list(orientation="h"))
  })

  output$threshold_table <- renderTable({
    tibble(
      Endpoint = c("AHI", "SpO2", "SBP", "HR", "ESS", "HOMA-IR", "hsCRP"),
      `Normal/Target` = c("<5 ev/hr", "≥95%", "<130 mmHg", "<80 bpm", "≤10", "<2.5", "<1.0 mg/L"),
      `Concern` = c("≥30", "<90%", "≥140", ">100", ">16", ">3.5", ">3.0 mg/L")
    )
  }, striped = TRUE, bordered = TRUE)

  ## ── Scenario Comparison plots ──────────────────────────────
  output$compare_ahi <- renderPlotly({
    d <- multi_sim()
    plot_ly(d %>% filter(week %in% seq(0, 52, by = 1)),
            x=~week, y=~AHI, color=~scenario, type="scatter", mode="lines",
            line=list(width=2)) %>%
      add_segments(x=0,xend=52,y=30,yend=30,line=list(dash="dash",color="red"),
                   name="Severe threshold") %>%
      add_segments(x=0,xend=52,y=15,yend=15,line=list(dash="dash",color="orange"),
                   name="Moderate threshold") %>%
      add_segments(x=0,xend=52,y=5,yend=5,line=list(dash="dash",color="green"),
                   name="Normal threshold") %>%
      layout(title="AHI by Treatment Scenario",
             xaxis=list(title="Weeks"),
             yaxis=list(title="AHI (events/hr)"),
             legend=list(orientation="h"))
  })

  output$compare_ess <- renderPlotly({
    d <- multi_sim()
    plot_ly(d %>% filter(week %in% seq(0, 52, by = 1)),
            x=~week, y=~ESS, color=~scenario, type="scatter", mode="lines") %>%
      add_segments(x=0,xend=52,y=10,yend=10,
                   line=list(dash="dash",color="orange"), name="EDS threshold") %>%
      layout(title="ESS by Scenario",
             xaxis=list(title="Weeks"), yaxis=list(title="ESS"))
  })

  output$compare_sbp <- renderPlotly({
    d <- multi_sim()
    plot_ly(d %>% filter(week %in% seq(0, 52, by = 1)),
            x=~week, y=~SBP, color=~scenario, type="scatter", mode="lines") %>%
      add_segments(x=0,xend=52,y=130,yend=130,
                   line=list(dash="dash",color="black"), name="130 mmHg") %>%
      layout(title="SBP by Scenario",
             xaxis=list(title="Weeks"), yaxis=list(title="SBP (mmHg)"))
  })

  ## ── Biomarker Heatmap ──────────────────────────────────────
  output$bm_heatmap <- renderPlotly({
    d <- sim_data()
    weeks_sel <- seq(0, input$sim_weeks, by = max(1, floor(input$sim_weeks/20)))
    d_sel <- d %>%
      filter(round(week) %in% weeks_sel) %>%
      group_by(round(week)) %>% slice_tail(n=1) %>% ungroup()

    bm_names <- c("AHI","SpO2","SBP","HR","CRP","HOMA","ESS","Weight")
    init_vals <- unlist(d[1, bm_names])
    mat <- d_sel %>%
      select(all_of(bm_names)) %>%
      sweep(2, init_vals, "/") %>%
      as.matrix()
    rownames(mat) <- round(d_sel$week, 1)

    # For SpO2: higher is better (invert for color scale)
    mat[, "SpO2"] <- 2 - mat[, "SpO2"]

    plot_ly(z = t(mat),
            x = rownames(mat),
            y = bm_names,
            type = "heatmap",
            colorscale = list(c(0,"#1565C0"), c(0.5,"#FFFFFF"), c(1,"#B71C1C")),
            zmid = 1.0) %>%
      layout(title = "Biomarker Heatmap (Ratio to Baseline)",
             xaxis = list(title = "Week"),
             yaxis = list(title = "Biomarker"))
  })

  ## ── Export table ───────────────────────────────────────────
  output$export_table <- renderDT({
    d <- sim_data() %>%
      filter(week %in% seq(0, input$sim_weeks, by = 4)) %>%
      select(week, AHI, SpO2, SBP, HR, CRP, HOMA, ESS, Weight, BMI, LG, AT, Pcrit) %>%
      mutate(across(where(is.numeric), ~round(.x, 2))) %>%
      rename(Week=week)
    datatable(d, options=list(pageLength=15, scrollX=TRUE), rownames=FALSE)
  })

  output$dl_csv <- downloadHandler(
    filename = function() paste0("osa_qsp_", Sys.Date(), ".csv"),
    content  = function(file) write.csv(sim_data(), file, row.names = FALSE)
  )
}

## ── Launch ────────────────────────────────────────────────────
shinyApp(ui = ui, server = server)
