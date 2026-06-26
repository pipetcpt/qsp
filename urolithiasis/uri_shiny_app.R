################################################################################
# Chronic Recurrent Urolithiasis — Interactive QSP Shiny Dashboard
# 6 Tabs: Patient Profile · PK Profiles · Urine Chemistry · Stone Risk ·
#         Scenario Comparison · Biomarker Panel
################################################################################

library(shiny)
library(shinydashboard)
library(shinyWidgets)
library(ggplot2)
library(dplyr)
library(tidyr)
library(DT)
library(plotly)

# ─────────────────────────────────────────────────────────────────────────────
# Core ODE model functions (lightweight, no mrgsolve dependency for Shiny)
# ─────────────────────────────────────────────────────────────────────────────

compute_urine_chem <- function(
    Ca_diet = 1000, Ox_diet = 200, Pro_diet = 96, Na_diet = 150,
    FluIn = 2.0, VitC = 500, PH1 = FALSE, MetSyn = FALSE, IBD = FALSE,
    HPT = FALSE, PTH_0 = 40,
    Dose_HCTZ = 0, Dose_KCit = 0, Dose_Allo = 0, Dose_Tam = 0,
    Lumasiran = FALSE
) {
  # Urine volume
  UV <- max(0.5, 1.5 + 0.85 * (FluIn - 1.5))

  # PTH
  PTH <- PTH_0 * (1 + HPT * 2)

  # Calcium absorption & urine Ca
  Ca_abs_frac <- 0.30 * (1 + 0.20 * PTH / (60 + PTH))
  uCa <- Ca_diet * Ca_abs_frac + 50 + Na_diet * 0.6
  uCa <- uCa * (1 - min(0.45, 0.45 * Dose_HCTZ / (25 + Dose_HCTZ)))

  # Oxalate
  Ox_hepatic <- 15 * (1 + PH1 * 4) * (1 - Lumasiran * 0.53)
  Ox_gut_abs <- (Ox_diet * 0.10 * (1 + IBD * 1.0) + VitC * 0.005)
  uOx <- Ox_hepatic + Ox_gut_abs

  # Uric Acid
  UA_prod <- 800 + Pro_diet * 3 + MetSyn * 200
  XO_inhib <- min(0.9, 0.9 * Dose_Allo / (300 + Dose_Allo))
  uUA <- UA_prod * (1 - XO_inhib)

  # Citrate
  KCit_eff <- min(0.6, 0.6 * Dose_KCit / (60 + Dose_KCit))
  uCit <- 600 * (1 - 0.3 * (PTH / PTH_0 - 1)) * (1 + KCit_eff)
  uCit <- max(50, uCit)

  # Urine pH
  urinepH <- 5.8 - MetSyn * 0.5 + KCit_eff * 1.2

  # Magnesium (simplified)
  uMg <- 80 + 0.4 * Ca_diet * 0.05

  # Supersaturation (CaOx, CaP, UA)
  Ca_mM  <- (uCa / 400.07) / UV
  Ox_mM  <- (uOx / 88.02)  / UV
  UA_mM  <- (uUA / 168.11) / UV
  Cit_mM <- (uCit / 192.12) / UV
  Mg_mM  <- (uMg / 24.31)  / UV

  freeCa_mM <- Ca_mM * (1 - 0.35 * Cit_mM / (Cit_mM + 2) - 0.25 * Mg_mM / (Mg_mM + 0.5))

  Ksp_CaOx <- 2.57e-9
  Ksp_CaP  <- 2.35e-7
  Ksp_UA   <- 9.8e-5

  SS_CaOx <- (freeCa_mM * 1e-3) * (Ox_mM * 1e-3) / Ksp_CaOx
  SS_CaP  <- (freeCa_mM * 1e-3) * (0.8e-3) / Ksp_CaP
  SS_UA   <- (UA_mM * 1e-3) / (Ksp_UA * 10^(urinepH - 5.35))

  list(
    uCa = round(uCa, 1), uOx = round(uOx, 1), uUA = round(uUA, 1),
    uCit = round(uCit, 1), UV = round(UV, 2), urinepH = round(urinepH, 2),
    uMg = round(uMg, 1),
    SS_CaOx = round(max(0, SS_CaOx), 2),
    SS_CaP  = round(max(0, SS_CaP),  2),
    SS_UA   = round(max(0, SS_UA),   2),
    Hypercalciuria  = uCa > 300,
    Hyperoxaluria   = uOx > 45,
    Hyperuricosuria = uUA > 800,
    Hypocitraturia  = uCit < 320,
    LowUV           = UV < 2.0,
    Stone_risk      = min(10, SS_CaOx * 3 + (uCa > 300) * 1.5 +
                         (uOx > 45) * 2 + (uCit < 320) * 1 +
                         (UV < 2.0) * 1.5)
  )
}

# Stone size trajectory (simplified ODE)
compute_stone_trajectory <- function(
    SS_CaOx, SS_UA, GFR_0 = 90, duration_yr = 5, stone0 = 2,
    Dose_Tam = 0
) {
  times <- seq(0, duration_yr, length.out = 300)
  dt    <- diff(times)[1]
  stone <- numeric(length(times))
  gfr   <- numeric(length(times))
  infl  <- numeric(length(times))
  stone[1] <- stone0
  gfr[1]   <- GFR_0
  infl[1]  <- 0.05

  k_grow  <- 0.5   # mm/yr per SS unit above 1
  n_grow  <- 1.5
  k_pass  <- 0.08 * (1 + 0.28 * Dose_Tam / (0.4 + Dose_Tam))
  k_infl  <- 0.10
  k_decay <- 0.30
  k_GFR   <- 0.5

  for (i in 2:length(times)) {
    SS_eff     <- max(0, SS_CaOx - 1)^n_grow
    crystal_l  <- max(0, SS_CaOx - 1)
    d_stone    <- (k_grow * SS_eff - k_pass * stone[i-1]) * dt
    d_infl     <- (k_infl * crystal_l - k_decay * infl[i-1]) * dt
    d_gfr      <- -k_GFR * infl[i-1] * dt
    stone[i]   <- max(0, stone[i-1] + d_stone)
    infl[i]    <- max(0, infl[i-1] + d_infl)
    gfr[i]     <- max(15, gfr[i-1] + d_gfr)
  }

  data.frame(Time_yr = times, Stone_mm = stone, GFR = gfr, Inflammation = infl)
}

# PK concentration curve (1-compartment oral)
pk_curve <- function(dose, ka, CL, V, F, times_h) {
  if (dose == 0) return(rep(0, length(times_h)))
  k_elim <- CL / V
  dose_abs <- dose * F
  Cp <- (dose_abs * ka) / (V * (ka - k_elim)) *
    (exp(-k_elim * times_h) - exp(-ka * times_h))
  pmax(0, Cp)
}

# ─────────────────────────────────────────────────────────────────────────────
# UI
# ─────────────────────────────────────────────────────────────────────────────
ui <- dashboardPage(
  skin = "blue",
  dashboardHeader(title = "Urolithiasis QSP Dashboard", titleWidth = 320),

  dashboardSidebar(
    width = 270,
    sidebarMenu(
      id = "tabs",
      menuItem("Patient Profile",    tabName = "tab_patient",   icon = icon("user")),
      menuItem("PK Profiles",        tabName = "tab_pk",        icon = icon("chart-line")),
      menuItem("Urine Chemistry",    tabName = "tab_urine",     icon = icon("flask")),
      menuItem("Stone Risk",         tabName = "tab_risk",      icon = icon("exclamation-triangle")),
      menuItem("Scenario Comparison",tabName = "tab_scenario",  icon = icon("balance-scale")),
      menuItem("Biomarker Panel",    tabName = "tab_biomarker", icon = icon("dna"))
    ),
    hr(),
    tags$div(style = "padding: 8px 15px; font-size:11px; color:#aaa;",
             "Model calibrated to ILLUMINATE-A,",
             "Ettinger NEJM 1986, Pak 1985,",
             "Pearle AUA Guidelines 2014")
  ),

  dashboardBody(
    tags$head(tags$style(HTML("
      .content-wrapper { background-color: #f4f7f9; }
      .box { border-radius: 8px; }
      .info-box { border-radius: 6px; }
      .risk-badge { padding: 4px 10px; border-radius: 12px;
                    font-weight:bold; font-size:13px; }
    "))),

    tabItems(

      # ── TAB 1: PATIENT PROFILE ─────────────────────────────────────────────
      tabItem(tabName = "tab_patient",
        fluidRow(
          box(title = "Demographics & Comorbidities", width = 4, status = "primary",
              solidHeader = TRUE,
              selectInput("sex", "Sex", c("Male"=1, "Female"=0), selected=1),
              sliderInput("age",  "Age (years)", 20, 80, 45, step=1),
              sliderInput("bmi",  "BMI (kg/m²)",  18, 45, 28, step=1),
              sliderInput("GFR0", "Baseline GFR (mL/min)", 15, 120, 90, step=5),
              hr(),
              checkboxInput("HPT",    "Primary Hyperparathyroidism", FALSE),
              checkboxInput("MetSyn", "Metabolic Syndrome / T2DM",    FALSE),
              checkboxInput("IBD",    "IBD / Bariatric Surgery",       FALSE),
              checkboxInput("PH1",    "Primary Hyperoxaluria Type 1",  FALSE),
              checkboxInput("Cystin", "Cystinuria",                    FALSE)
          ),
          box(title = "Dietary & Lifestyle Inputs", width = 4, status = "warning",
              solidHeader = TRUE,
              sliderInput("Ca_diet",   "Dietary Ca (mg/day)", 200, 2000, 1000, step=50),
              sliderInput("Ox_diet",   "Dietary OX (mg/day)", 10,  800,  200,  step=10),
              sliderInput("Pro_diet",  "Animal Protein (g/day)", 20, 250, 96,  step=5),
              sliderInput("Na_diet",   "Dietary Na (mEq/day)", 40, 300, 150, step=10),
              sliderInput("FluIn",     "Fluid Intake (L/day)", 0.5, 5.0, 2.0, step=0.1),
              sliderInput("VitC",      "Vitamin C (mg/day)",   0, 2000, 500, step=50),
              sliderInput("PTH_0",     "Baseline PTH (pg/mL)", 10, 200, 40,  step=5)
          ),
          box(title = "Drug Doses", width = 4, status = "success",
              solidHeader = TRUE,
              h5("Thiazide Diuretic"),
              sliderInput("Dose_HCTZ", "HCTZ (mg/day)", 0, 50, 0, step=12.5),
              h5("Allopurinol / Febuxostat"),
              sliderInput("Dose_Allo", "Allopurinol (mg/day)", 0, 600, 0, step=100),
              h5("Potassium Citrate"),
              sliderInput("Dose_KCit", "K-Citrate (mEq/day)", 0, 120, 0, step=10),
              h5("Medical Expulsive Therapy"),
              sliderInput("Dose_Tam", "Tamsulosin (mg/day)", 0, 0.8, 0, step=0.4),
              checkboxInput("Lumasiran", "Lumasiran (PH1 siRNA therapy)", FALSE),
              hr(),
              sliderInput("duration_yr", "Simulation Duration (years)", 1, 10, 5, step=1),
              sliderInput("stone0_mm",   "Initial Stone Size (mm)", 0.5, 20, 2, step=0.5)
          )
        ),
        fluidRow(
          box(title = "Patient Summary", width = 12, status = "info",
              tableOutput("patient_summary_tbl"))
        )
      ),

      # ── TAB 2: PK PROFILES ────────────────────────────────────────────────
      tabItem(tabName = "tab_pk",
        fluidRow(
          box(title = "HCTZ Plasma Concentration (24h)", width = 6, status = "primary",
              plotlyOutput("pk_HCTZ", height = 280)),
          box(title = "Allopurinol & Oxypurinol (24h)", width = 6, status = "primary",
              plotlyOutput("pk_Allo", height = 280))
        ),
        fluidRow(
          box(title = "K-Citrate Plasma (24h)", width = 6, status = "primary",
              plotlyOutput("pk_KCit", height = 280)),
          box(title = "Tamsulosin Plasma (24h)", width = 6, status = "primary",
              plotlyOutput("pk_Tam", height = 280))
        ),
        fluidRow(
          box(title = "PK Parameters Summary", width = 12, status = "info",
              DTOutput("pk_params_tbl"))
        )
      ),

      # ── TAB 3: URINE CHEMISTRY ─────────────────────────────────────────────
      tabItem(tabName = "tab_urine",
        fluidRow(
          infoBoxOutput("box_uCa",  width = 3),
          infoBoxOutput("box_uOx",  width = 3),
          infoBoxOutput("box_uUA",  width = 3),
          infoBoxOutput("box_uCit", width = 3)
        ),
        fluidRow(
          infoBoxOutput("box_UV",      width = 3),
          infoBoxOutput("box_pH",      width = 3),
          infoBoxOutput("box_SS_CaOx", width = 3),
          infoBoxOutput("box_SS_UA",   width = 3)
        ),
        fluidRow(
          box(title = "24h Urine Abnormalities", width = 6, status = "warning",
              solidHeader = TRUE,
              tableOutput("urine_flags_tbl")),
          box(title = "Stone Risk Score", width = 6, status = "danger",
              solidHeader = TRUE,
              plotlyOutput("risk_gauge", height = 250))
        ),
        fluidRow(
          box(title = "Supersaturation Comparison (CaOx / CaP / UA)", width = 12,
              status = "info", solidHeader = TRUE,
              plotlyOutput("ss_bar", height = 280))
        )
      ),

      # ── TAB 4: STONE RISK ─────────────────────────────────────────────────
      tabItem(tabName = "tab_risk",
        fluidRow(
          box(title = "Stone Size Trajectory", width = 7, status = "danger",
              solidHeader = TRUE,
              plotlyOutput("stone_traj", height = 350)),
          box(title = "GFR Trajectory", width = 5, status = "warning",
              solidHeader = TRUE,
              plotlyOutput("gfr_traj", height = 350))
        ),
        fluidRow(
          box(title = "Renal Inflammation Index", width = 6, status = "info",
              plotlyOutput("infl_traj", height = 280)),
          box(title = "Stone Type Probability", width = 6, status = "primary",
              plotlyOutput("stone_type_pie", height = 280))
        ),
        fluidRow(
          box(title = "Passage Probability by Size", width = 12, status = "success",
              plotlyOutput("passage_plot", height = 200))
        )
      ),

      # ── TAB 5: SCENARIO COMPARISON ────────────────────────────────────────
      tabItem(tabName = "tab_scenario",
        fluidRow(
          box(title = "Define Comparison Scenarios", width = 4, status = "primary",
              solidHeader = TRUE,
              h5("Scenario A (current settings from Tab 1)"),
              hr(),
              h5("Scenario B — Override Settings"),
              sliderInput("scB_HCTZ", "B: HCTZ (mg/day)",  0, 50, 25, step=12.5),
              sliderInput("scB_KCit", "B: K-Citrate (mEq)", 0, 120, 0, step=10),
              sliderInput("scB_Allo", "B: Allopurinol (mg)", 0, 600, 0, step=100),
              sliderInput("scB_Flu",  "B: Fluid (L/day)",   0.5, 5, 2.5, step=0.5),
              hr(),
              h5("Scenario C — Intensive Lifestyle"),
              sliderInput("scC_Ca",   "C: Dietary Ca (mg)", 500, 2000, 1000, step=100),
              sliderInput("scC_Flu",  "C: Fluid (L/day)",    0.5, 5, 3.0, step=0.5),
              sliderInput("scC_Na",   "C: Na (mEq/day)",     40, 300, 80, step=10),
              sliderInput("scC_Prot", "C: Protein (g/day)",  20, 250, 70, step=5)
          ),
          box(title = "Stone Size: A vs B vs C", width = 8, status = "success",
              solidHeader = TRUE,
              plotlyOutput("scen_stone", height = 320),
              plotlyOutput("scen_SS",    height = 220))
        ),
        fluidRow(
          box(title = "Urine Chemistry Comparison (Bar)", width = 12, status = "info",
              plotlyOutput("scen_urine_bar", height = 300))
        )
      ),

      # ── TAB 6: BIOMARKER PANEL ────────────────────────────────────────────
      tabItem(tabName = "tab_biomarker",
        fluidRow(
          box(title = "24h Urine Profile (Radar)", width = 6, status = "primary",
              solidHeader = TRUE,
              plotlyOutput("radar_urine", height = 400)),
          box(title = "Target vs. Achieved Biomarkers", width = 6, status = "success",
              solidHeader = TRUE,
              DTOutput("biomarker_table", height = 380))
        ),
        fluidRow(
          box(title = "Treatment Response: ΔBiomarkers", width = 6, status = "warning",
              plotlyOutput("delta_biomarker", height = 300)),
          box(title = "Calcium Balance", width = 6, status = "info",
              plotlyOutput("ca_balance_plot", height = 300))
        ),
        fluidRow(
          box(title = "Model Information & Evidence Base", width = 12,
              status = "primary",
              tags$ul(
                tags$li(strong("CaOx SS Model:"), " EQUIL2 approximation; calibrated to Lingeman et al. 2003"),
                tags$li(strong("Thiazide PK:"), " HCTZ F=0.65, t½=6-15h, 2-cpt (Beermann & Groschinsky-Grind 1977)"),
                tags$li(strong("Allopurinol/Oxypurinol:"), " Mechanism-based XO inhibition (IC₅₀=0.7µM, Pacher 2006)"),
                tags$li(strong("K-Citrate:"), " Pak et al. NEJM 2014; 50-60% reduction in CaOx recurrence"),
                tags$li(strong("Lumasiran (PH1):"), " ILLUMINATE-A (NEJM 2021): 53% uOx reduction"),
                tags$li(strong("Stone growth:"), " Kok et al. 1990; Moe 2006; Preminger 2009"),
                tags$li(strong("GFR impact:"), " Alexander et al. JASN 2012 — stones → 1.3× CKD risk"),
                tags$li(strong("Randall's plaque:"), " Evan et al. 2003; Matlaga 2007")
              ))
        )
      )
    )
  )
)

# ─────────────────────────────────────────────────────────────────────────────
# SERVER
# ─────────────────────────────────────────────────────────────────────────────
server <- function(input, output, session) {

  # Reactive: compute urine chemistry for current inputs (Tab 1 settings)
  urine_chem <- reactive({
    compute_urine_chem(
      Ca_diet = input$Ca_diet, Ox_diet = input$Ox_diet,
      Pro_diet = input$Pro_diet, Na_diet = input$Na_diet,
      FluIn = input$FluIn, VitC = input$VitC,
      PH1 = input$PH1, MetSyn = input$MetSyn, IBD = input$IBD,
      HPT = input$HPT, PTH_0 = input$PTH_0,
      Dose_HCTZ = input$Dose_HCTZ, Dose_KCit = input$Dose_KCit,
      Dose_Allo = input$Dose_Allo, Dose_Tam = input$Dose_Tam,
      Lumasiran = input$Lumasiran
    )
  })

  # Reactive: stone trajectory for current settings
  stone_traj_data <- reactive({
    uc <- urine_chem()
    compute_stone_trajectory(
      SS_CaOx = uc$SS_CaOx, SS_UA = uc$SS_UA,
      GFR_0 = input$GFR0, duration_yr = input$duration_yr,
      stone0 = input$stone0_mm, Dose_Tam = input$Dose_Tam
    )
  })

  # ── TAB 1: Patient Summary ──────────────────────────────────────────────
  output$patient_summary_tbl <- renderTable({
    uc <- urine_chem()
    data.frame(
      Parameter = c("Sex", "Age", "BMI", "Baseline GFR",
                    "Daily Fluid Intake", "Urine Volume",
                    "PTH", "HPT", "Metabolic Syndrome", "PH1"),
      Value = c(
        ifelse(input$sex == 1, "Male", "Female"),
        paste(input$age, "years"),
        paste(input$bmi, "kg/m²"),
        paste(input$GFR0, "mL/min/1.73m²"),
        paste(input$FluIn, "L/day"),
        paste(round(uc$UV, 2), "L/day"),
        paste(input$PTH_0, "pg/mL"),
        ifelse(input$HPT, "Yes", "No"),
        ifelse(input$MetSyn, "Yes", "No"),
        ifelse(input$PH1, "Yes", "No")
      )
    )
  })

  # ── TAB 2: PK Profiles ──────────────────────────────────────────────────
  times_h <- seq(0, 24, length.out = 200)

  output$pk_HCTZ <- renderPlotly({
    Cp <- pk_curve(input$Dose_HCTZ, ka=1.5, CL=18, V=40, F=0.65, times_h)
    plot_ly(x=~times_h, y=~Cp, type="scatter", mode="lines",
            line=list(color="steelblue", width=2.5)) %>%
      layout(xaxis=list(title="Time (h)"),
             yaxis=list(title="Concentration (mg/L)"),
             title=list(text="HCTZ Plasma Conc."))
  })

  output$pk_Allo <- renderPlotly({
    Cp_allo <- pk_curve(input$Dose_Allo, ka=2.0, CL=10, V=35, F=0.9, times_h)
    Cp_oxp  <- pk_curve(input$Dose_Allo, ka=0.15, CL=2.5, V=100, F=0.9, times_h)
    plot_ly(x=~times_h) %>%
      add_trace(y=~Cp_allo, name="Allopurinol", mode="lines",
                line=list(color="royalblue", width=2)) %>%
      add_trace(y=~Cp_oxp,  name="Oxypurinol", mode="lines",
                line=list(color="coral", width=2, dash="dash")) %>%
      layout(xaxis=list(title="Time (h)"),
             yaxis=list(title="Concentration (mg/L)"),
             legend=list(orientation="h"))
  })

  output$pk_KCit <- renderPlotly({
    Cp <- pk_curve(input$Dose_KCit, ka=1.2, CL=15, V=25, F=0.95, times_h)
    plot_ly(x=~times_h, y=~Cp, type="scatter", mode="lines",
            line=list(color="darkorange", width=2.5)) %>%
      layout(xaxis=list(title="Time (h)"),
             yaxis=list(title="Concentration (mEq/L)"))
  })

  output$pk_Tam <- renderPlotly({
    Cp <- pk_curve(input$Dose_Tam, ka=0.8, CL=3.5, V=65, F=0.9, times_h)
    plot_ly(x=~times_h, y=~Cp*1000, type="scatter", mode="lines",
            line=list(color="purple", width=2.5)) %>%
      layout(xaxis=list(title="Time (h)"),
             yaxis=list(title="Concentration (µg/L)"))
  })

  output$pk_params_tbl <- renderDT({
    df <- data.frame(
      Drug = c("HCTZ", "Allopurinol", "Oxypurinol", "K-Citrate", "Tamsulosin"),
      `Dose (current)` = c(input$Dose_HCTZ, input$Dose_Allo, "—",
                            input$Dose_KCit, input$Dose_Tam),
      F = c("0.65", "0.90", "—", "0.95", "0.90"),
      `t½ (h)` = c("6-15", "0.5-1.5", "18-30", "1-2", "9-16"),
      `Vd (L)` = c("40+80", "35", "100", "25", "65"),
      `CL (L/h)` = c("18", "10", "2.5", "15", "3.5"),
      `Primary effect` = c("↓uCa 30-45%", "↓XO activity", "↓XO (active)",
                           "↑uCit, ↑pH", "↑Stone passage"),
      check.names = FALSE
    )
    datatable(df, options = list(dom = "t"), rownames = FALSE)
  })

  # ── TAB 3: Urine Chemistry info boxes ──────────────────────────────────
  mk_infobox <- function(val, label, threshold_lo, threshold_hi,
                          unit, better = "low") {
    ok <- if (better == "low") val <= threshold_hi else val >= threshold_lo
    color <- if (ok) "green" else "red"
    icon_nm <- if (ok) "check-circle" else "exclamation-circle"
    infoBox(label, paste0(val, " ", unit), icon = icon(icon_nm),
            color = color, fill = TRUE)
  }

  output$box_uCa  <- renderInfoBox({
    uc <- urine_chem()
    thr <- if (input$sex == 1) 300 else 250
    color <- if (uc$uCa <= thr) "green" else "red"
    infoBox("Urinary Ca", paste0(uc$uCa, " mg/day"),
            subtitle = paste0("Target: <", thr, " mg/day"),
            icon = icon(if(uc$uCa<=thr)"check-circle" else "exclamation-circle"),
            color = color, fill = TRUE)
  })

  output$box_uOx  <- renderInfoBox({
    uc <- urine_chem()
    color <- if (uc$uOx <= 45) "green" else "red"
    infoBox("Urinary Oxalate", paste0(uc$uOx, " mg/day"),
            subtitle = "Target: <45 mg/day",
            icon = icon(if(uc$uOx<=45)"check-circle" else "exclamation-circle"),
            color = color, fill = TRUE)
  })

  output$box_uUA  <- renderInfoBox({
    uc <- urine_chem()
    thr <- if (input$sex == 1) 800 else 750
    color <- if (uc$uUA <= thr) "green" else "red"
    infoBox("Urinary UA", paste0(uc$uUA, " mg/day"),
            subtitle = paste0("Target: <", thr, " mg/day"),
            icon = icon(if(uc$uUA<=thr)"check-circle" else "exclamation-circle"),
            color = color, fill = TRUE)
  })

  output$box_uCit <- renderInfoBox({
    uc <- urine_chem()
    color <- if (uc$uCit >= 320) "green" else "red"
    infoBox("Urinary Citrate", paste0(uc$uCit, " mg/day"),
            subtitle = "Target: ≥320 mg/day",
            icon = icon(if(uc$uCit>=320)"check-circle" else "exclamation-circle"),
            color = color, fill = TRUE)
  })

  output$box_UV <- renderInfoBox({
    uc <- urine_chem()
    color <- if (uc$UV >= 2.0) "green" else "red"
    infoBox("Urine Volume", paste0(uc$UV, " L/day"),
            subtitle = "Target: ≥2.0 L/day",
            icon = icon(if(uc$UV>=2)"check-circle" else "exclamation-circle"),
            color = color, fill = TRUE)
  })

  output$box_pH <- renderInfoBox({
    uc <- urine_chem()
    ok <- uc$urinepH >= 6.0 && uc$urinepH <= 7.0
    color <- if (ok) "green" else "orange"
    infoBox("Urine pH", uc$urinepH,
            subtitle = "Optimal: 6.0-7.0",
            icon = icon("tint"), color = color, fill = TRUE)
  })

  output$box_SS_CaOx <- renderInfoBox({
    uc <- urine_chem()
    color <- if (uc$SS_CaOx < 1) "green" else if (uc$SS_CaOx < 3) "yellow" else "red"
    infoBox("CaOx SS", uc$SS_CaOx,
            subtitle = "Target: <1.0 (no crystallization)",
            icon = icon("atom"), color = color, fill = TRUE)
  })

  output$box_SS_UA <- renderInfoBox({
    uc <- urine_chem()
    color <- if (uc$SS_UA < 1) "green" else "red"
    infoBox("UA SS", uc$SS_UA,
            subtitle = "Target: <1.0",
            icon = icon("atom"), color = color, fill = TRUE)
  })

  output$urine_flags_tbl <- renderTable({
    uc <- urine_chem()
    data.frame(
      Abnormality = c("Hypercalciuria", "Hyperoxaluria", "Hyperuricosuria",
                      "Hypocitraturia", "Low Urine Volume"),
      Present = c(uc$Hypercalciuria, uc$Hyperoxaluria, uc$Hyperuricosuria,
                  uc$Hypocitraturia, uc$LowUV),
      `Threshold` = c("M>300, F>250 mg/day", ">45 mg/day", "M>800, F>750 mg/day",
                      "<320 mg/day", "<2.0 L/day"),
      `Risk Contribution` = c("High", "Very High", "Moderate", "Moderate", "High"),
      check.names = FALSE
    ) %>% mutate(Present = ifelse(Present, "YES ⚠", "No ✓"))
  })

  output$risk_gauge <- renderPlotly({
    uc <- urine_chem()
    risk <- round(min(10, uc$Stone_risk), 1)
    color <- if (risk < 3) "#27ae60" else if (risk < 6) "#f39c12" else "#e74c3c"

    plot_ly(type = "indicator", mode = "gauge+number",
            value = risk,
            title = list(text = "Stone Risk Score (0-10)"),
            gauge = list(
              axis = list(range = list(0, 10)),
              bar = list(color = color),
              steps = list(
                list(range=c(0,3),   color="#d5f5e3"),
                list(range=c(3,6),   color="#fef9e7"),
                list(range=c(6,10),  color="#fadbd8")
              ),
              threshold = list(
                line = list(color="red", width=4),
                thickness = 0.75, value = 7
              )
            )) %>%
      layout(margin = list(t=50, b=10))
  })

  output$ss_bar <- renderPlotly({
    uc <- urine_chem()
    df_ss <- data.frame(
      SS_type = c("CaOx SS", "CaP SS", "UA SS"),
      Value   = c(uc$SS_CaOx, uc$SS_CaP, uc$SS_UA),
      Color   = c("#e74c3c", "#e67e22", "#8e44ad")
    )
    plot_ly(df_ss, x=~SS_type, y=~Value, type="bar", color=~SS_type,
            colors = df_ss$Color) %>%
      add_segments(x=0.5, xend=3.5, y=1, yend=1,
                   line=list(dash="dash", color="red", width=2),
                   name="Crystallization threshold") %>%
      layout(xaxis=list(title=""),
             yaxis=list(title="Supersaturation (relative to Ksp)"),
             showlegend=FALSE)
  })

  # ── TAB 4: Stone Risk ───────────────────────────────────────────────────
  output$stone_traj <- renderPlotly({
    df <- stone_traj_data()
    plot_ly(df, x=~Time_yr, y=~Stone_mm, type="scatter", mode="lines",
            line=list(color="#e74c3c", width=2.5), name="Stone size") %>%
      add_segments(x=0, xend=max(df$Time_yr), y=5, yend=5,
                   line=list(dash="dash", color="orange", width=1.5),
                   name="5mm (URS threshold)") %>%
      add_segments(x=0, xend=max(df$Time_yr), y=10, yend=10,
                   line=list(dash="dash", color="red", width=1.5),
                   name="10mm (PCNL threshold)") %>%
      layout(xaxis=list(title="Time (years)"),
             yaxis=list(title="Stone Size (mm)"),
             legend=list(orientation="h", y=-0.2))
  })

  output$gfr_traj <- renderPlotly({
    df <- stone_traj_data()
    plot_ly(df, x=~Time_yr, y=~GFR, type="scatter", mode="lines",
            line=list(color="steelblue", width=2.5)) %>%
      add_segments(x=0, xend=max(df$Time_yr), y=60, yend=60,
                   line=list(dash="dash", color="#f1c40f", width=1.5),
                   name="CKD G3a") %>%
      add_segments(x=0, xend=max(df$Time_yr), y=30, yend=30,
                   line=list(dash="dash", color="orange", width=1.5),
                   name="CKD G4") %>%
      layout(xaxis=list(title="Time (years)"),
             yaxis=list(title="eGFR (mL/min/1.73m²)"),
             legend=list(orientation="h", y=-0.2))
  })

  output$infl_traj <- renderPlotly({
    df <- stone_traj_data()
    plot_ly(df, x=~Time_yr, y=~Inflammation, type="scatter", mode="lines",
            fill="tozeroy", fillcolor="rgba(231,76,60,0.2)",
            line=list(color="#e74c3c", width=2)) %>%
      layout(xaxis=list(title="Time (years)"),
             yaxis=list(title="Inflammation Index (AU)"))
  })

  output$stone_type_pie <- renderPlotly({
    uc <- urine_chem()
    # Rough probability of each stone type based on risk factors
    p_caox <- 0.70 - 0.3 * uc$Hyperuricosuria - 0.1 * input$IBD + 0.1 * uc$Hypercalciuria
    p_caox <- max(0.1, min(0.85, p_caox))
    p_ua   <- 0.10 + 0.2 * uc$Hyperuricosuria + 0.15 * input$MetSyn
    p_ua   <- max(0.02, min(0.6, p_ua))
    p_cap  <- 0.08 + 0.1 * (input$PTH_0 > 80)
    p_strv <- 0.05
    p_cys  <- if (input$Cystin) 0.8 else 0.02

    tot <- p_caox + p_ua + p_cap + p_strv + p_cys
    df_pie <- data.frame(
      Type = c("CaOx", "UA", "CaP", "Struvite", "Cystine"),
      Prob = c(p_caox, p_ua, p_cap, p_strv, p_cys) / tot
    )
    plot_ly(df_pie, labels=~Type, values=~Prob, type="pie",
            textinfo="label+percent",
            marker=list(colors=c("#e74c3c","#8e44ad","#e67e22","#27ae60","#2980b9")))
  })

  output$passage_plot <- renderPlotly({
    sizes <- c(1,2,3,4,5,6,7,8,9,10)
    base_pass <- c(0.87, 0.80, 0.73, 0.62, 0.48, 0.35, 0.22, 0.14, 0.08, 0.05)
    tam_boost  <- pmin(1, base_pass * (1 + 0.28 * input$Dose_Tam / (0.4 + input$Dose_Tam)))

    plot_ly(x=sizes) %>%
      add_trace(y=base_pass, name="No MET", mode="lines+markers",
                line=list(color="steelblue"), marker=list(size=6)) %>%
      add_trace(y=tam_boost, name=paste0("Tamsulosin ", input$Dose_Tam, "mg"),
                mode="lines+markers", line=list(color="darkorange", dash="dash"),
                marker=list(size=6)) %>%
      layout(xaxis=list(title="Stone Size (mm)"),
             yaxis=list(title="Spontaneous Passage Rate", range=c(0,1),
                        tickformat=".0%"),
             legend=list(orientation="h", y=-0.3))
  })

  # ── TAB 5: Scenario Comparison ─────────────────────────────────────────
  scenarios_reactive <- reactive({
    # Scenario A: current settings
    ucA <- urine_chem()
    trajA <- compute_stone_trajectory(SS_CaOx = ucA$SS_CaOx, SS_UA = ucA$SS_UA,
                                       GFR_0 = input$GFR0,
                                       duration_yr = input$duration_yr,
                                       stone0 = input$stone0_mm) %>%
      mutate(Scenario = "A: Current")

    # Scenario B
    ucB <- compute_urine_chem(
      Ca_diet = input$Ca_diet, Ox_diet = input$Ox_diet,
      Pro_diet = input$Pro_diet, Na_diet = input$Na_diet,
      FluIn = input$scB_Flu, PH1 = input$PH1, MetSyn = input$MetSyn,
      HPT = input$HPT, PTH_0 = input$PTH_0,
      Dose_HCTZ = input$scB_HCTZ, Dose_KCit = input$scB_KCit,
      Dose_Allo = input$scB_Allo
    )
    trajB <- compute_stone_trajectory(SS_CaOx = ucB$SS_CaOx, SS_UA = ucB$SS_UA,
                                       GFR_0 = input$GFR0,
                                       duration_yr = input$duration_yr,
                                       stone0 = input$stone0_mm) %>%
      mutate(Scenario = "B: Drug Rx")

    # Scenario C: Lifestyle
    ucC <- compute_urine_chem(
      Ca_diet = input$scC_Ca, Ox_diet = 100, Pro_diet = input$scC_Prot,
      Na_diet = input$scC_Na, FluIn = input$scC_Flu,
      PH1 = input$PH1, MetSyn = input$MetSyn, HPT = input$HPT,
      PTH_0 = input$PTH_0
    )
    trajC <- compute_stone_trajectory(SS_CaOx = ucC$SS_CaOx, SS_UA = ucC$SS_UA,
                                       GFR_0 = input$GFR0,
                                       duration_yr = input$duration_yr,
                                       stone0 = input$stone0_mm) %>%
      mutate(Scenario = "C: Lifestyle")

    list(traj = bind_rows(trajA, trajB, trajC),
         ucA = ucA, ucB = ucB, ucC = ucC)
  })

  output$scen_stone <- renderPlotly({
    sc <- scenarios_reactive()
    plot_ly(sc$traj, x=~Time_yr, y=~Stone_mm, color=~Scenario,
            type="scatter", mode="lines",
            colors=c("#e74c3c","#2980b9","#27ae60")) %>%
      layout(xaxis=list(title="Time (years)"),
             yaxis=list(title="Stone Size (mm)"),
             legend=list(orientation="h"))
  })

  output$scen_SS <- renderPlotly({
    sc <- scenarios_reactive()
    df_ss <- data.frame(
      Scenario = c("A: Current","B: Drug Rx","C: Lifestyle"),
      SS_CaOx  = c(sc$ucA$SS_CaOx, sc$ucB$SS_CaOx, sc$ucC$SS_CaOx)
    )
    plot_ly(df_ss, x=~Scenario, y=~SS_CaOx, type="bar",
            marker=list(color=c("#e74c3c","#2980b9","#27ae60"))) %>%
      add_segments(x=0.5, xend=3.5, y=1, yend=1,
                   line=list(dash="dash", color="red"),
                   name="SS=1 threshold") %>%
      layout(xaxis=list(title=""), yaxis=list(title="CaOx SS"),
             showlegend=FALSE)
  })

  output$scen_urine_bar <- renderPlotly({
    sc <- scenarios_reactive()
    vars <- c("uCa","uOx","uUA","uCit","UV")
    labs <- c("Urinary Ca (mg/d)","Urinary OX (mg/d)","Urinary UA (mg/d)",
              "Urinary Citrate (mg/d)","Urine Volume (L/d × 100)")

    df_bar <- data.frame(
      Variable = rep(labs, 3),
      Value    = c(
        c(sc$ucA$uCa, sc$ucA$uOx, sc$ucA$uUA, sc$ucA$uCit, sc$ucA$UV * 100),
        c(sc$ucB$uCa, sc$ucB$uOx, sc$ucB$uUA, sc$ucB$uCit, sc$ucB$UV * 100),
        c(sc$ucC$uCa, sc$ucC$uOx, sc$ucC$uUA, sc$ucC$uCit, sc$ucC$UV * 100)
      ),
      Scenario = rep(c("A: Current","B: Drug Rx","C: Lifestyle"), each=5)
    )

    plot_ly(df_bar, x=~Variable, y=~Value, color=~Scenario, type="bar",
            colors=c("#e74c3c","#2980b9","#27ae60")) %>%
      layout(barmode="group", xaxis=list(title=""),
             yaxis=list(title="Value"), legend=list(orientation="h"))
  })

  # ── TAB 6: Biomarker Panel ─────────────────────────────────────────────
  output$radar_urine <- renderPlotly({
    uc <- urine_chem()
    # Normalize to target (1.0 = at target)
    cats <- c("Urine Vol","Citrate","Urine pH","↓Ca","↓Oxalate","↓UA")
    actual <- c(
      min(2, uc$UV) / 2.0,
      min(800, uc$uCit) / 800,
      min(7.5, max(4, uc$urinepH)) / 7.5,
      1 - min(1, uc$uCa / 400),
      1 - min(1, uc$uOx / 100),
      1 - min(1, uc$uUA / 1000)
    )
    target <- rep(1.0, 6)

    plot_ly(type="scatterpolar", mode="lines+markers") %>%
      add_trace(r=c(target, target[1]), theta=c(cats, cats[1]),
                name="Target", line=list(color="green", dash="dash")) %>%
      add_trace(r=c(actual, actual[1]), theta=c(cats, cats[1]),
                name="Current", fill="toself", fillcolor="rgba(41,128,185,0.3)",
                line=list(color="steelblue")) %>%
      layout(polar=list(radialaxis=list(range=c(0,1.2))),
             showlegend=TRUE, legend=list(x=0.8, y=1.1))
  })

  output$biomarker_table <- renderDT({
    uc <- urine_chem()
    df_bm <- data.frame(
      Biomarker = c("24h Urine Ca (mg/day)", "24h Urine Oxalate (mg/day)",
                    "24h Urine UA (mg/day)", "24h Urine Citrate (mg/day)",
                    "Urine Volume (L/day)", "Urine pH",
                    "CaOx Supersaturation", "UA Supersaturation"),
      Current   = c(uc$uCa, uc$uOx, uc$uUA, uc$uCit,
                    uc$UV, uc$urinepH, uc$SS_CaOx, uc$SS_UA),
      Target    = c(
        ifelse(input$sex == 1, "<300", "<250"),
        "<45", ifelse(input$sex == 1, "<800", "<750"),
        ">320", ">2.0", "6.0-7.0", "<1.0", "<1.0"
      ),
      Status = c(
        ifelse(uc$Hypercalciuria, "ABNORMAL", "NORMAL"),
        ifelse(uc$Hyperoxaluria,  "ABNORMAL", "NORMAL"),
        ifelse(uc$Hyperuricosuria,"ABNORMAL", "NORMAL"),
        ifelse(uc$Hypocitraturia, "ABNORMAL", "NORMAL"),
        ifelse(uc$LowUV,          "ABNORMAL", "NORMAL"),
        ifelse(uc$urinepH < 6.0 | uc$urinepH > 7.0, "REVIEW", "NORMAL"),
        ifelse(uc$SS_CaOx > 1,   "HIGH RISK", "OK"),
        ifelse(uc$SS_UA   > 1,   "HIGH RISK", "OK")
      )
    )
    datatable(df_bm, options = list(dom="t", pageLength=10), rownames=FALSE) %>%
      formatStyle("Status",
                  backgroundColor = styleEqual(
                    c("ABNORMAL","HIGH RISK","REVIEW","NORMAL","OK"),
                    c("#fadbd8","#f1948a","#fef9e7","#d5f5e3","#d5f5e3")
                  ))
  })

  output$delta_biomarker <- renderPlotly({
    # Baseline (no treatment)
    uc_base <- compute_urine_chem(
      Ca_diet = input$Ca_diet, Ox_diet = input$Ox_diet,
      Pro_diet = input$Pro_diet, Na_diet = input$Na_diet,
      FluIn = input$FluIn, VitC = input$VitC,
      PH1 = input$PH1, MetSyn = input$MetSyn, IBD = input$IBD,
      HPT = input$HPT, PTH_0 = input$PTH_0,
      Dose_HCTZ = 0, Dose_KCit = 0, Dose_Allo = 0
    )
    uc_rx <- urine_chem()

    delta_Ca  <- round((uc_rx$uCa  - uc_base$uCa)  / uc_base$uCa  * 100, 1)
    delta_Ox  <- round((uc_rx$uOx  - uc_base$uOx)  / uc_base$uOx  * 100, 1)
    delta_UA  <- round((uc_rx$uUA  - uc_base$uUA)  / uc_base$uUA  * 100, 1)
    delta_Cit <- round((uc_rx$uCit - uc_base$uCit) / uc_base$uCit * 100, 1)
    delta_SS  <- round((uc_rx$SS_CaOx - uc_base$SS_CaOx) / max(0.01, uc_base$SS_CaOx) * 100, 1)

    df_delta <- data.frame(
      Param = c("Urinary Ca","Urinary OX","Urinary UA","Urinary Citrate","CaOx SS"),
      Delta = c(delta_Ca, delta_Ox, delta_UA, delta_Cit, delta_SS)
    )
    df_delta$color <- ifelse(df_delta$Delta < 0, "#27ae60", "#e74c3c")
    # For Citrate, positive is good
    df_delta$color[df_delta$Param == "Urinary Citrate"] <-
      ifelse(df_delta$Delta[df_delta$Param == "Urinary Citrate"] > 0, "#27ae60", "#e74c3c")

    plot_ly(df_delta, x=~Param, y=~Delta, type="bar",
            marker=list(color=df_delta$color)) %>%
      add_segments(x=0.5, xend=5.5, y=0, yend=0,
                   line=list(color="black", width=1)) %>%
      layout(xaxis=list(title=""),
             yaxis=list(title="% Change vs. No Treatment"),
             showlegend=FALSE)
  })

  output$ca_balance_plot <- renderPlotly({
    uc <- urine_chem()
    Ca_absorbed <- input$Ca_diet * 0.30 * (1 + 0.2 * input$PTH_0 / (60 + input$PTH_0))
    Ca_bone_resorp <- 50
    Ca_urine <- uc$uCa
    Ca_balance <- Ca_absorbed + Ca_bone_resorp - Ca_urine - 200  # ~200 fecal

    df_ca <- data.frame(
      Component = c("GI Absorbed","Bone Resorption","Urinary Loss","Fecal/Other","Net Balance"),
      Value     = c(Ca_absorbed, Ca_bone_resorp, -Ca_urine, -200, Ca_balance)
    )
    df_ca$color <- ifelse(df_ca$Value >= 0, "#27ae60", "#e74c3c")

    plot_ly(df_ca, x=~Component, y=~Value, type="bar",
            marker=list(color=df_ca$color)) %>%
      add_segments(x=0.5, xend=5.5, y=0, yend=0,
                   line=list(color="black", width=1.5)) %>%
      layout(xaxis=list(title=""),
             yaxis=list(title="Calcium (mg/day)"),
             showlegend=FALSE,
             title=list(text="Daily Calcium Balance"))
  })
}

# ─────────────────────────────────────────────────────────────────────────────
shinyApp(ui, server)
