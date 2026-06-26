## ============================================================
## AIHA QSP Shiny Dashboard
## Autoimmune Hemolytic Anemia — Interactive Simulator
## Tabs: 1.Patient Profile · 2.Drug PK · 3.PD Markers
##       4.Clinical Endpoints · 5.Scenario Comparison · 6.Biomarkers
## ============================================================

library(shiny)
library(shinydashboard)
library(ggplot2)
library(dplyr)
library(tidyr)
library(DT)
library(plotly)
library(scales)

# ============================================================
# ODE SIMULATION FUNCTION (pure-R fallback)
# ============================================================
simulate_aiha <- function(params) {
  p <- params
  dt <- 0.5
  times <- seq(0, p$sim_days, by = dt)
  n <- length(times)

  # State vectors
  BCELL  <- numeric(n); BCELL[1]  <- p$BCELL_init
  PLASMA <- numeric(n); PLASMA[1] <- p$PLASMA_init
  AB     <- numeric(n); AB[1]     <- p$AB_init
  C3b    <- numeric(n); C3b[1]    <- p$C3b_init
  RBC    <- numeric(n); RBC[1]    <- p$RBC_init
  RETI   <- numeric(n); RETI[1]   <- p$RETI_init
  EPO    <- numeric(n); EPO[1]    <- 80
  LDH    <- numeric(n); LDH[1]    <- p$LDH_init
  HAP    <- numeric(n); HAP[1]    <- p$HAP_init
  BILIR  <- numeric(n); BILIR[1]  <- p$BILIR_init

  # Drug PK (single compartment)
  PRED   <- numeric(n); PRED[1]   <- 0
  RTX    <- numeric(n); RTX[1]    <- 0
  SUTI   <- numeric(n); SUTI[1]   <- 0
  FOSTA  <- numeric(n); FOSTA[1]  <- 0
  MMF    <- numeric(n); MMF[1]    <- 0
  IVIG   <- numeric(n); IVIG[1]   <- 0

  # Helper: dose event (returns dose amount at time t)
  dose_at <- function(t, times_doses, amt, dur = 0.25) {
    sum(amt[(t >= times_doses) & (t < times_doses + dur)])
  }

  # Drug PK parameters
  CL_PRED <- 15; Vd_PRED <- 45; ka_PRED <- 2.88
  CL_RTX  <- 0.33; Vd_RTX <- 4.4
  CL_SUTI <- 0.42; Vd_SUTI <- 5.0
  CL_FOSTA <- 35; Vd_FOSTA <- 250; ka_FOSTA <- 2.88
  CL_MMF  <- 25; Vd_MMF <- 100; ka_MMF <- 3.0
  CL_IVIG <- 0.20; Vd_IVIG <- 3.5

  # Dose schedules
  pred_times <- if (p$use_pred) seq(0, min(p$pred_duration, p$sim_days), by = 1) else c()
  rtx_times  <- if (p$use_rtx)  c(0, 7, 14, 21) else c()
  suti_times <- if (p$use_suti) seq(0, p$sim_days, by = 14) else c()
  fosta_times<- if (p$use_fosta)seq(0, p$sim_days, by = 0.5) else c()
  mmf_times  <- if (p$use_mmf)  seq(0, p$sim_days, by = 0.5) else c()
  ivig_times <- if (p$use_ivig) c(0, 1) else c()

  pred_dose  <- p$pred_dose
  rtx_dose   <- 700 / Vd_RTX  # nM equivalent
  suti_dose  <- p$suti_dose / Vd_SUTI / 150  # simplified unit
  fosta_dose <- p$fosta_dose
  mmf_dose   <- 1000
  ivig_dose  <- p$weight * 1.0  # g = 1g/kg

  for (i in seq_len(n - 1)) {
    t <- times[i]
    d <- dt

    # Drug infusion/absorption rates
    inPRED  <- if (length(pred_times) > 0 && any(abs(t - pred_times) < d)) pred_dose else 0
    inRTX   <- if (length(rtx_times)  > 0 && any(abs(t - rtx_times)  < d)) rtx_dose else 0
    inSUTI  <- if (length(suti_times) > 0 && any(abs(t - suti_times)  < d)) suti_dose else 0
    inFOSTA <- if (length(fosta_times)> 0 && any(abs(t - fosta_times) < d)) fosta_dose else 0
    inMMF   <- if (length(mmf_times)  > 0 && any(abs(t - mmf_times)   < d)) mmf_dose else 0
    inIVIG  <- if (length(ivig_times) > 0 && any(abs(t - ivig_times)  < d)) ivig_dose else 0

    # Drug PK (Euler)
    PRED[i+1]  <- max(0, PRED[i]  + d * (inPRED  * ka_PRED / Vd_PRED   - CL_PRED / Vd_PRED  * PRED[i]))
    RTX[i+1]   <- max(0, RTX[i]   + d * (inRTX                           - CL_RTX  / Vd_RTX   * RTX[i]))
    SUTI[i+1]  <- max(0, SUTI[i]  + d * (inSUTI                          - CL_SUTI / Vd_SUTI  * SUTI[i]))
    FOSTA[i+1] <- max(0, FOSTA[i] + d * (inFOSTA * ka_FOSTA / Vd_FOSTA  - CL_FOSTA/ Vd_FOSTA * FOSTA[i]))
    MMF[i+1]   <- max(0, MMF[i]   + d * (inMMF   * ka_MMF   / Vd_MMF    - CL_MMF  / Vd_MMF   * MMF[i]))
    IVIG[i+1]  <- max(0, IVIG[i]  + d * (inIVIG                          - CL_IVIG / Vd_IVIG  * IVIG[i]))

    # Drug PD effects
    EC50_PRED <- 0.05
    eff_PRED_FcR <- 0.75 * PRED[i] / (EC50_PRED + PRED[i])
    eff_PRED_Ab  <- 0.55 * PRED[i] / (EC50_PRED + PRED[i])
    eff_RTX_Bkill <- 0.95 * RTX[i] / (2.0 + RTX[i])
    eff_FOSTA_phago <- 0.80 * FOSTA[i] / (0.041 + FOSTA[i])
    eff_SUTI_C3b <- 0.90 * SUTI[i] / (0.10 + SUTI[i])
    eff_MMF_B  <- 0.70 * MMF[i] / (0.50 + MMF[i])
    eff_IVIG_FcR <- 0.85 * IVIG[i] / (4.0 + IVIG[i])

    GC_FcR <- max(eff_PRED_FcR, eff_PRED_FcR * 0.9)
    phago_mod <- (1 - GC_FcR) * (1 - eff_FOSTA_phago) * (1 - eff_IVIG_FcR)
    Ab_mod    <- (1 - eff_PRED_Ab) * (1 - eff_MMF_B)

    # Disease ODE
    # B cells
    dBCELL <- 0.018 * BCELL[i] * (1 - eff_MMF_B) - 0.020 * BCELL[i] - eff_RTX_Bkill * BCELL[i]
    BCELL[i+1] <- max(0, BCELL[i] + d * dBCELL)

    # Plasma cells
    dPLASMA <- 0.003 * BCELL[i]^2 / (BCELL[i] + 50) - 0.015 * PLASMA[i] - 0.3 * eff_RTX_Bkill * PLASMA[i]
    PLASMA[i+1] <- max(0, PLASMA[i] + d * dPLASMA)

    # Autoantibody
    dAB <- 0.025 * PLASMA[i] * Ab_mod - 0.004 * AB[i]
    AB[i+1] <- max(0, AB[i] + d * dAB)

    # C3b
    C3b_drive <- if (p$subtype == 2) 1.2 * (1 - eff_SUTI_C3b) * AB[i]/(AB[i]+2) else 0.3 * AB[i]/(AB[i]+5)
    dC3b <- C3b_drive - 0.30 * C3b[i]
    C3b[i+1] <- max(0, C3b[i] + d * dC3b)

    # Opsonization
    OpsonWarm <- AB[i] / (AB[i] + 1.0)
    OpsonCold <- C3b[i] / (C3b[i] + 5.0)
    Opson <- if (p$subtype == 1) OpsonWarm else 0.7*OpsonCold + 0.3*OpsonWarm

    # RBC
    Hb_i <- RBC[i] * 0.002
    EPO_stim <- 15 * (1 + 2.5 * max(0, 10 - Hb_i))
    EPO[i+1] <- max(15, EPO[i] + d * (EPO_stim - 0.10 * EPO[i]))
    dRETI <- 0.08 * EPO[i] / (EPO[i] + 15) * 500 - 1.0 * RETI[i]
    RETI[i+1] <- max(0, RETI[i] + d * dRETI)

    Rate_phago <- 0.015 * RBC[i] * Opson * phago_mod
    Rate_lysis <- if (p$subtype == 2) 0.02 * RBC[i] * C3b[i]/(C3b[i]+5) * (1-eff_SUTI_C3b) else 0.001 * RBC[i]

    dRBC <- 1.0 * RETI[i] + 2.0 - 0.0083 * RBC[i] - Rate_phago - Rate_lysis
    RBC[i+1] <- max(0, RBC[i] + d * dRBC)

    # Biomarkers
    Hemolysis <- Rate_phago + Rate_lysis
    LDH[i+1]   <- max(180, LDH[i]   + d * (80 * Hemolysis - 0.50 * (LDH[i] - 180)))
    HAP[i+1]   <- max(0,   HAP[i]   + d * (0.60 * (1.5 - HAP[i]) - 0.05 * Rate_lysis))
    BILIR[i+1] <- max(0.5, BILIR[i] + d * (0.015 * Rate_phago - 2.0 * (BILIR[i] - 0.5)))
  }

  Hemoglobin <- RBC * 0.002
  data.frame(
    time = times,
    Hemoglobin, Hematocrit = Hemoglobin * 3,
    RBC, RETI,
    Reticulocyte_pct = RETI / (RBC + RETI + 0.001) * 100,
    AB, BCELL, PLASMA, C3b, EPO,
    LDH, Haptoglobin = HAP, Bilirubin = BILIR,
    PRED_conc = PRED, RTX_conc = RTX, SUTI_conc = SUTI,
    FOSTA_conc = FOSTA, MMF_conc = MMF, IVIG_conc = IVIG
  )
}

# ============================================================
# UI
# ============================================================
ui <- dashboardPage(
  skin = "red",
  dashboardHeader(title = "AIHA QSP Simulator"),

  dashboardSidebar(
    sidebarMenu(
      menuItem("1. Patient Profile",    tabName = "tab_patient",   icon = icon("user")),
      menuItem("2. Drug PK",            tabName = "tab_pk",        icon = icon("pills")),
      menuItem("3. PD Core Markers",    tabName = "tab_pd",        icon = icon("microscope")),
      menuItem("4. Clinical Endpoints", tabName = "tab_clinical",  icon = icon("chart-line")),
      menuItem("5. Scenario Comparison",tabName = "tab_compare",   icon = icon("balance-scale")),
      menuItem("6. Biomarkers",         tabName = "tab_biomarker", icon = icon("flask")),
      menuItem("About",                 tabName = "tab_about",     icon = icon("info-circle"))
    )
  ),

  dashboardBody(
    tags$head(tags$style(HTML("
      .skin-red .main-header .logo { background-color: #B71C1C; }
      .skin-red .main-header .navbar { background-color: #C62828; }
      .content-wrapper { background-color: #F5F5F5; }
      .box-header { background-color: #FFCDD2 !important; }
    "))),

    tabItems(
      # ======================================================
      # TAB 1: Patient Profile
      # ======================================================
      tabItem(tabName = "tab_patient",
        fluidRow(
          box(title = "Patient Demographics", status = "danger", solidHeader = TRUE, width = 4,
            numericInput("weight", "Body Weight (kg)", value = 68, min = 40, max = 120),
            sliderInput("age", "Age (years)", min = 18, max = 85, value = 45),
            selectInput("sex", "Sex", choices = c("Female", "Male"), selected = "Female"),
            selectInput("subtype", "AIHA Subtype",
                        choices = c("Warm AIHA (IgG)" = 1, "Cold Agglutinin Disease (IgM)" = 2),
                        selected = 1),
            selectInput("severity", "Disease Severity at Presentation",
                        choices = c("Mild (Hb 9-10)" = "mild",
                                    "Moderate (Hb 7-9)" = "moderate",
                                    "Severe (Hb <7)" = "severe"),
                        selected = "moderate"),
            numericInput("sim_days", "Simulation Duration (days)", value = 180, min = 30, max = 365)
          ),
          box(title = "Baseline Disease Characteristics", status = "danger", solidHeader = TRUE, width = 4,
            sliderInput("AB_init", "Initial Anti-RBC Ab (AU/mL)", min = 1, max = 20, value = 8, step = 0.5),
            sliderInput("RBC_init", "Baseline RBC (×10³/μL)", min = 1500, max = 4500, value = 3500, step = 100),
            sliderInput("LDH_init", "Baseline LDH (U/L)", min = 200, max = 1200, value = 450, step = 10),
            sliderInput("HAP_init", "Baseline Haptoglobin (g/L)", min = 0.0, max = 1.5, value = 0.2, step = 0.1),
            sliderInput("BILIR_init", "Baseline Unconj. Bilirubin (mg/dL)", min = 1.0, max = 6.0, value = 2.8, step = 0.1)
          ),
          box(title = "Disease Summary at Baseline", status = "warning", solidHeader = TRUE, width = 4,
            h4("Calculated Baseline Parameters"),
            tableOutput("baseline_table"),
            hr(),
            h4("AIHA Subtype Description"),
            uiOutput("subtype_description")
          )
        ),
        fluidRow(
          box(title = "Treatment Selection", status = "danger", solidHeader = TRUE, width = 6,
            checkboxInput("use_pred",  "Prednisolone", value = TRUE),
            conditionalPanel("input.use_pred",
              sliderInput("pred_dose", "Prednisolone dose (mg/day)", 10, 140, 70, step = 5),
              sliderInput("pred_duration", "Prednisolone duration (days)", 14, 84, 28)
            ),
            checkboxInput("use_rtx",   "Rituximab (375 mg/m² ×4 weekly)", value = FALSE),
            checkboxInput("use_suti",  "Sutimlimab (CAD only)", value = FALSE),
            conditionalPanel("input.use_suti",
              selectInput("suti_dose", "Sutimlimab dose",
                          choices = c("6500 mg (<75 kg)" = 6500, "7500 mg (≥75 kg)" = 7500))
            ),
            checkboxInput("use_fosta", "Fostamatinib 150 mg BID (refractory)", value = FALSE),
            checkboxInput("use_mmf",   "Mycophenolate mofetil 1000 mg BID", value = FALSE),
            checkboxInput("use_ivig",  "IVIG 1 g/kg ×2d (acute rescue)", value = FALSE),
            actionButton("run_sim", "Run Simulation", class = "btn-danger btn-lg", icon = icon("play"))
          ),
          box(title = "Treatment Notes", status = "info", solidHeader = TRUE, width = 6,
            tags$ul(
              tags$li(strong("Prednisolone:"), " First-line for Warm AIHA. Response in 70-80% at 3 weeks."),
              tags$li(strong("Rituximab:"), " Anti-CD20, second-line or combination. CR ~65% at 1 year."),
              tags$li(strong("Sutimlimab (Enjaymo):"), " FDA-approved for CAD (Feb 2022). Inhibits C1s."),
              tags$li(strong("Fostamatinib (Tavalia):"), " Syk inhibitor, FDA-approved for ITP, emerging evidence in AIHA."),
              tags$li(strong("MMF:"), " Steroid-sparing, used for maintenance in refractory AIHA."),
              tags$li(strong("IVIG:"), " Short-term FcγR blockade for acute hemolytic crisis.")
            )
          )
        )
      ),

      # ======================================================
      # TAB 2: Drug PK
      # ======================================================
      tabItem(tabName = "tab_pk",
        fluidRow(
          box(title = "Drug Plasma Concentration — Time Profiles", status = "primary",
              solidHeader = TRUE, width = 12,
              plotlyOutput("pk_plot", height = "450px")
          )
        ),
        fluidRow(
          box(title = "PK Parameter Summary", status = "primary", solidHeader = TRUE, width = 7,
            DTOutput("pk_param_table")
          ),
          box(title = "Drug PK Display Options", status = "info", solidHeader = TRUE, width = 5,
            checkboxGroupInput("pk_drugs", "Select drugs to display:",
              choices = c("Prednisolone" = "PRED_conc", "Rituximab" = "RTX_conc",
                          "Sutimlimab" = "SUTI_conc", "Fostamatinib" = "FOSTA_conc",
                          "MMF/MPA" = "MMF_conc", "IVIG" = "IVIG_conc"),
              selected = c("PRED_conc", "RTX_conc")),
            sliderInput("pk_timerange", "Time range (days)", min = 0, max = 365, value = c(0, 90))
          )
        )
      ),

      # ======================================================
      # TAB 3: PD Core Markers
      # ======================================================
      tabItem(tabName = "tab_pd",
        fluidRow(
          box(title = "B Cell & Autoantibody Dynamics", status = "warning",
              solidHeader = TRUE, width = 6,
              plotlyOutput("pd_bcell_plot", height = "350px")
          ),
          box(title = "Complement (C3b) & Erythropoiesis", status = "warning",
              solidHeader = TRUE, width = 6,
              plotlyOutput("pd_complement_plot", height = "350px")
          )
        ),
        fluidRow(
          box(title = "EPO & Reticulocyte Response", status = "success",
              solidHeader = TRUE, width = 6,
              plotlyOutput("pd_epo_plot", height = "350px")
          ),
          box(title = "PD Mechanism Summary", status = "info", solidHeader = TRUE, width = 6,
            tableOutput("pd_summary_table")
          )
        )
      ),

      # ======================================================
      # TAB 4: Clinical Endpoints
      # ======================================================
      tabItem(tabName = "tab_clinical",
        fluidRow(
          box(title = "Hemoglobin Trajectory", status = "danger",
              solidHeader = TRUE, width = 8,
              plotlyOutput("hb_plot", height = "380px")
          ),
          box(title = "Response Assessment", status = "success",
              solidHeader = TRUE, width = 4,
              h4("Clinical Response by Day"),
              tableOutput("response_table"),
              hr(),
              h4("Key Endpoints at Day 90"),
              uiOutput("endpoint_summary")
          )
        ),
        fluidRow(
          box(title = "Reticulocyte Count (%)", status = "primary",
              solidHeader = TRUE, width = 6,
              plotlyOutput("reti_plot", height = "300px")
          ),
          box(title = "Hematocrit (%)", status = "primary",
              solidHeader = TRUE, width = 6,
              plotlyOutput("hct_plot", height = "300px")
          )
        )
      ),

      # ======================================================
      # TAB 5: Scenario Comparison
      # ======================================================
      tabItem(tabName = "tab_compare",
        fluidRow(
          box(title = "Multi-Scenario Setup", status = "danger", solidHeader = TRUE, width = 12,
            p("Compare up to 4 treatment scenarios simultaneously on the same patient."),
            fluidRow(
              column(3,
                strong("Scenario A (Base)"),
                checkboxInput("scA_pred", "Prednisolone", value = TRUE),
                checkboxInput("scA_rtx",  "Rituximab", value = FALSE),
                checkboxInput("scA_suti", "Sutimlimab", value = FALSE),
                checkboxInput("scA_fosta","Fostamatinib", value = FALSE)
              ),
              column(3,
                strong("Scenario B"),
                checkboxInput("scB_pred", "Prednisolone", value = TRUE),
                checkboxInput("scB_rtx",  "Rituximab", value = TRUE),
                checkboxInput("scB_suti", "Sutimlimab", value = FALSE),
                checkboxInput("scB_fosta","Fostamatinib", value = FALSE)
              ),
              column(3,
                strong("Scenario C"),
                checkboxInput("scC_pred", "Prednisolone", value = FALSE),
                checkboxInput("scC_rtx",  "Rituximab", value = TRUE),
                checkboxInput("scC_suti", "Sutimlimab", value = FALSE),
                checkboxInput("scC_fosta","Fostamatinib", value = FALSE)
              ),
              column(3,
                strong("Scenario D"),
                checkboxInput("scD_pred", "Prednisolone", value = FALSE),
                checkboxInput("scD_rtx",  "Rituximab", value = FALSE),
                checkboxInput("scD_suti", "Sutimlimab", value = TRUE),
                checkboxInput("scD_fosta","Fostamatinib", value = TRUE)
              )
            ),
            actionButton("run_compare", "Compare All Scenarios", class = "btn-danger")
          )
        ),
        fluidRow(
          box(title = "Hemoglobin Comparison (All Scenarios)", solidHeader = TRUE,
              status = "primary", width = 8,
              plotlyOutput("compare_hb_plot", height = "400px")
          ),
          box(title = "Outcome Table at Day 90/180", solidHeader = TRUE,
              status = "info", width = 4,
              DTOutput("compare_table")
          )
        ),
        fluidRow(
          box(title = "LDH Comparison", solidHeader = TRUE, status = "warning", width = 6,
            plotlyOutput("compare_ldh_plot", height = "300px")
          ),
          box(title = "Autoantibody Comparison", solidHeader = TRUE, status = "warning", width = 6,
            plotlyOutput("compare_ab_plot", height = "300px")
          )
        )
      ),

      # ======================================================
      # TAB 6: Biomarkers
      # ======================================================
      tabItem(tabName = "tab_biomarker",
        fluidRow(
          box(title = "Serum LDH (Hemolysis Marker)", status = "danger",
              solidHeader = TRUE, width = 6,
              plotlyOutput("ldh_plot", height = "300px")
          ),
          box(title = "Haptoglobin Level", status = "warning",
              solidHeader = TRUE, width = 6,
              plotlyOutput("hap_plot", height = "300px")
          )
        ),
        fluidRow(
          box(title = "Unconjugated Bilirubin", status = "warning",
              solidHeader = TRUE, width = 6,
              plotlyOutput("bilir_plot", height = "300px")
          ),
          box(title = "Direct Antiglobulin Test (DAT) Proxy", status = "info",
              solidHeader = TRUE, width = 6,
              plotlyOutput("dat_plot", height = "300px")
          )
        ),
        fluidRow(
          box(title = "Biomarker Reference Ranges", solidHeader = TRUE,
              status = "primary", width = 12,
            DTOutput("biomarker_ref_table")
          )
        )
      ),

      # ======================================================
      # TAB 7: About
      # ======================================================
      tabItem(tabName = "tab_about",
        fluidRow(
          box(title = "About This Model", status = "danger", solidHeader = TRUE, width = 12,
            h3("Autoimmune Hemolytic Anemia (AIHA) — QSP Dashboard"),
            p("This interactive dashboard implements a Quantitative Systems Pharmacology (QSP)
              model for Autoimmune Hemolytic Anemia, covering both:"),
            tags$ul(
              tags$li(strong("Warm AIHA (wAIHA):"), " IgG-mediated extravascular hemolysis via splenic macrophage FcγR."),
              tags$li(strong("Cold Agglutinin Disease (CAD):"), " IgM-mediated complement classical pathway activation,
                C3b opsonization, and intravascular MAC formation.")
            ),
            h4("Model Structure (26 ODE Compartments)"),
            tags$ol(
              tags$li("B Cell / Plasma Cell dynamics (GC, SHM, CSR)"),
              tags$li("Anti-RBC autoantibody (IgG warm, IgM cold)"),
              tags$li("Complement cascade (C1s → C3b → C5b-9 MAC)"),
              tags$li("RBC destruction (extravascular + intravascular)"),
              tags$li("Compensatory erythropoiesis (EPO → Reticulocyte → RBC)"),
              tags$li("Biomarkers (LDH, Haptoglobin, Bilirubin, DAT)"),
              tags$li("Drug PK: Prednisolone, Rituximab, Sutimlimab, Fostamatinib, MMF, IVIG")
            ),
            h4("Key Calibration References"),
            tags$ul(
              tags$li("Röth A et al. NEJM 2021;384:1535 — CADENZA trial (sutimlimab CAD)"),
              tags$li("Barcellini W et al. Blood 2018;131:1534 — Rituximab warm AIHA"),
              tags$li("Giaimo ME et al. Am J Hematol 2020;95:E28 — Fostamatinib AIHA"),
              tags$li("Lechner K & Jäger U. Blood 2010;116:1831 — Prednisone treatment"),
              tags$li("Berentsen S et al. Haematologica 2020;105:1308 — CAD management")
            ),
            h4("Model Limitations"),
            tags$ul(
              tags$li("Simplified pure-R ODE solver (Euler method) — use mrgsolve version for clinical research"),
              tags$li("Inter-individual variability not implemented in Shiny version"),
              tags$li("TMDD model for rituximab simplified to 1-compartment for speed"),
              tags$li("Secondary AIHA (lymphoma, infections) not modeled explicitly")
            )
          )
        )
      )
    )
  )
)

# ============================================================
# SERVER
# ============================================================
server <- function(input, output, session) {

  # ---- Make baseline parameters from inputs ----
  make_params <- function(use_pred, use_rtx, use_suti, use_fosta, use_mmf, use_ivig) {
    list(
      subtype    = as.integer(input$subtype),
      weight     = input$weight,
      sim_days   = input$sim_days,
      AB_init    = input$AB_init,
      C3b_init   = ifelse(as.integer(input$subtype) == 2, 3.0, 0.5),
      RBC_init   = input$RBC_init,
      RETI_init  = 200,
      BCELL_init = 300,
      PLASMA_init= 50,
      LDH_init   = input$LDH_init,
      HAP_init   = input$HAP_init,
      BILIR_init = input$BILIR_init,
      use_pred   = use_pred,  pred_dose = input$pred_dose, pred_duration = input$pred_duration,
      use_rtx    = use_rtx,
      use_suti   = use_suti,  suti_dose = as.numeric(input$suti_dose),
      use_fosta  = use_fosta, fosta_dose = 150,
      use_mmf    = use_mmf,
      use_ivig   = use_ivig
    )
  }

  # Reactive simulation result
  sim_result <- eventReactive(input$run_sim, {
    withProgress(message = "Running simulation...", {
      p <- make_params(input$use_pred, input$use_rtx, input$use_suti,
                       input$use_fosta, input$use_mmf, input$use_ivig)
      simulate_aiha(p)
    })
  }, ignoreNULL = FALSE)

  # Scenario comparison results
  compare_results <- eventReactive(input$run_compare, {
    withProgress(message = "Running 4 scenarios...", {
      scenarios <- list(
        A = make_params(input$scA_pred, input$scA_rtx, input$scA_suti, input$scA_fosta, FALSE, FALSE),
        B = make_params(input$scB_pred, input$scB_rtx, input$scB_suti, input$scB_fosta, FALSE, FALSE),
        C = make_params(input$scC_pred, input$scC_rtx, input$scC_suti, input$scC_fosta, FALSE, FALSE),
        D = make_params(input$scD_pred, input$scD_rtx, input$scD_suti, input$scD_fosta, FALSE, FALSE)
      )
      lapply(scenarios, simulate_aiha)
    })
  })

  # ---- TAB 1: Baseline Table ----
  output$baseline_table <- renderTable({
    Hb_init <- input$RBC_init * 0.002
    data.frame(
      Parameter = c("Hemoglobin", "Hematocrit", "RBC", "LDH", "Haptoglobin", "Bilirubin (Unc.)"),
      Value     = c(paste0(round(Hb_init, 1), " g/dL"),
                    paste0(round(Hb_init * 3, 1), "%"),
                    paste0(input$RBC_init, " ×10³/μL"),
                    paste0(input$LDH_init, " U/L"),
                    paste0(input$HAP_init, " g/L"),
                    paste0(input$BILIR_init, " mg/dL")),
      Normal_Range = c("12-16 g/dL", "36-48%", "4000-5500 ×10³/μL",
                       "<250 U/L", "0.5-2.5 g/L", "0.2-1.2 mg/dL")
    )
  })

  output$subtype_description <- renderUI({
    if (input$subtype == 1) {
      tagList(
        tags$b("Warm AIHA (IgG-mediated)"),
        tags$p("Anti-RBC IgG autoantibodies opsonize RBCs and trigger FcγR-mediated
               phagocytosis by splenic macrophages (extravascular hemolysis).
               DAT: IgG+ (±C3d). Treatment: Steroids → Rituximab.")
      )
    } else {
      tagList(
        tags$b("Cold Agglutinin Disease (IgM-mediated)"),
        tags$p("Monoclonal IgM cold agglutinins bind I/i antigen on RBCs at low temperatures,
               activating the classical complement pathway (C1q → C3b → MAC).
               Hemolysis: intravascular + extravascular (C3b opsonization).
               DAT: C3d+ (IgG-). Treatment: Sutimlimab (Enjaymo), Rituximab.")
      )
    }
  })

  # ---- TAB 2: PK Plot ----
  output$pk_plot <- renderPlotly({
    df <- sim_result()
    if (is.null(df)) return(NULL)

    drugs_sel <- input$pk_drugs
    if (length(drugs_sel) == 0) return(NULL)

    tr <- input$pk_timerange
    df_filt <- df %>% filter(time >= tr[1], time <= tr[2])

    drug_labels <- c(PRED_conc = "Prednisolone (mg/L)", RTX_conc = "Rituximab (nM)",
                     SUTI_conc = "Sutimlimab (nM)",    FOSTA_conc = "Fostamatinib/R406 (µM)",
                     MMF_conc = "MMF/MPA (µg/mL)",     IVIG_conc = "IVIG (g/L)")

    df_long <- df_filt %>%
      select(time, all_of(drugs_sel)) %>%
      pivot_longer(-time, names_to = "drug", values_to = "conc")

    p <- ggplot(df_long, aes(x = time, y = conc, color = drug)) +
      geom_line(linewidth = 1.2) +
      facet_wrap(~drug, scales = "free_y", ncol = 2) +
      labs(x = "Time (days)", y = "Concentration", title = "Drug PK Profiles") +
      theme_bw(base_size = 11) +
      theme(legend.position = "none")
    ggplotly(p)
  })

  output$pk_param_table <- renderDT({
    df <- data.frame(
      Drug        = c("Prednisolone", "Rituximab", "Sutimlimab", "Fostamatinib (R406)", "MMF (MPA)", "IVIG"),
      Admin       = c("PO", "IV", "IV", "PO", "PO", "IV"),
      Half_life   = c("2.5h", "22 days", "20 days", "14h", "18h", "21 days"),
      Vd          = c("0.5-1.0 L/kg", "4.4 L", "5.0 L", "250 L", "100 L", "3.5 L"),
      CL          = c("15 L/h", "0.33 L/day", "0.42 L/day", "35 L/day", "25 L/day", "0.20 L/day"),
      Bioavail    = c("80%", "IV (100%)", "IV (100%)", "~55% (R406)", "94% (MPA)", "IV (100%)"),
      Mechanism   = c("GR agonist", "Anti-CD20", "Anti-C1s (CP block)", "Syk inhibitor", "IMPDH inhibitor", "FcRn/FcgR block")
    )
    datatable(df, options = list(pageLength = 10, dom = 't'), rownames = FALSE)
  })

  # ---- TAB 3: PD Plots ----
  output$pd_bcell_plot <- renderPlotly({
    df <- sim_result()
    if (is.null(df)) return(NULL)
    p <- ggplot(df, aes(x = time)) +
      geom_line(aes(y = BCELL, color = "B Cells (×10⁶/L)"), linewidth = 1) +
      geom_line(aes(y = PLASMA * 5, color = "Plasma Cells (×10⁶/L ×5)"), linewidth = 1) +
      geom_line(aes(y = AB * 30, color = "Anti-RBC Ab (AU×30)"), linewidth = 1) +
      labs(x = "Time (days)", y = "Cell/Ab Level",
           title = "B Cell & Autoantibody Dynamics", color = "") +
      theme_bw(base_size = 11)
    ggplotly(p)
  })

  output$pd_complement_plot <- renderPlotly({
    df <- sim_result()
    if (is.null(df)) return(NULL)
    p <- ggplot(df, aes(x = time)) +
      geom_line(aes(y = C3b * 50, color = "C3b on RBC (AU×50)"), linewidth = 1) +
      geom_line(aes(y = RBC / 50, color = "RBC/50 (×10³/μL)"), linewidth = 1) +
      labs(x = "Time (days)", y = "Level (normalized)",
           title = "Complement C3b & RBC", color = "") +
      theme_bw(base_size = 11)
    ggplotly(p)
  })

  output$pd_epo_plot <- renderPlotly({
    df <- sim_result()
    if (is.null(df)) return(NULL)
    p <- ggplot(df, aes(x = time)) +
      geom_line(aes(y = EPO, color = "EPO (mIU/mL)"), linewidth = 1) +
      geom_line(aes(y = RETI / 5, color = "Reticulocytes/5 (×10³/μL)"), linewidth = 1) +
      geom_line(aes(y = Reticulocyte_pct * 10, color = "Reti% ×10"), linewidth = 1) +
      labs(x = "Time (days)", y = "Level", title = "EPO & Reticulocyte Response", color = "") +
      theme_bw(base_size = 11)
    ggplotly(p)
  })

  output$pd_summary_table <- renderTable({
    df <- sim_result()
    if (is.null(df)) return(NULL)
    snap <- df %>% filter(abs(time - 90) < 1) %>% slice(1)
    data.frame(
      Parameter = c("B Cells (×10⁶/L)", "Plasma Cells (×10⁶/L)",
                    "Anti-RBC Ab (AU/mL)", "C3b on RBC (AU)",
                    "EPO (mIU/mL)", "Reticulocytes (×10³/μL)"),
      Day_0     = c(300, 50, input$AB_init, ifelse(input$subtype==2, 3.0, 0.5), 80, 200),
      Day_90    = round(c(snap$BCELL, snap$PLASMA, snap$AB, snap$C3b,
                          snap$EPO, snap$RETI), 2)
    )
  })

  # ---- TAB 4: Clinical Endpoints ----
  output$hb_plot <- renderPlotly({
    df <- sim_result()
    if (is.null(df)) return(NULL)
    p <- ggplot(df, aes(x = time, y = Hemoglobin)) +
      geom_line(color = "#C62828", linewidth = 1.5) +
      geom_hline(yintercept = 10, linetype = "dashed", color = "orange") +
      geom_hline(yintercept = 7, linetype = "dotted", color = "red") +
      annotate("text", x = max(df$time)*0.7, y = 10.3, label = "CR threshold (10 g/dL)", size = 3.5) +
      annotate("text", x = max(df$time)*0.7, y = 7.3, label = "Transfusion threshold (7 g/dL)", size = 3.5, color = "red") +
      labs(x = "Time (days)", y = "Hemoglobin (g/dL)",
           title = "Hemoglobin Trajectory",
           subtitle = "Complete Response = Hb ≥ 10 g/dL") +
      scale_y_continuous(limits = c(5, 16)) +
      theme_bw(base_size = 12)
    ggplotly(p)
  })

  output$reti_plot <- renderPlotly({
    df <- sim_result()
    if (is.null(df)) return(NULL)
    p <- ggplot(df, aes(x = time, y = Reticulocyte_pct)) +
      geom_line(color = "#1565C0", linewidth = 1) +
      geom_hline(yintercept = 2.0, linetype = "dashed", color = "gray50") +
      labs(x = "Time (days)", y = "Reticulocyte (%)",
           title = "Reticulocyte Percentage") +
      theme_bw(base_size = 11)
    ggplotly(p)
  })

  output$hct_plot <- renderPlotly({
    df <- sim_result()
    if (is.null(df)) return(NULL)
    p <- ggplot(df, aes(x = time, y = Hematocrit)) +
      geom_line(color = "#4A148C", linewidth = 1) +
      geom_hline(yintercept = 36, linetype = "dashed", color = "gray50") +
      labs(x = "Time (days)", y = "Hematocrit (%)",
           title = "Hematocrit Trajectory") +
      theme_bw(base_size = 11)
    ggplotly(p)
  })

  output$response_table <- renderTable({
    df <- sim_result()
    if (is.null(df)) return(NULL)
    check_days <- c(14, 30, 60, 90, 120, 180)
    check_days <- check_days[check_days <= max(df$time)]
    res <- lapply(check_days, function(d) {
      snap <- df %>% filter(abs(time - d) < 1) %>% slice(1)
      data.frame(
        Day = d,
        Hb  = round(snap$Hemoglobin, 1),
        Response = case_when(
          snap$Hemoglobin >= 10 ~ "CR",
          snap$Hemoglobin >= 8  ~ "PR",
          snap$Hemoglobin >= 6  ~ "SD",
          TRUE ~ "Refractory"
        )
      )
    })
    bind_rows(res)
  }, striped = TRUE)

  output$endpoint_summary <- renderUI({
    df <- sim_result()
    if (is.null(df)) return(NULL)
    snap90 <- df %>% filter(abs(time - 90) < 1) %>% slice(1)
    tagList(
      tags$p(strong("Hb at Day 90: "), round(snap90$Hemoglobin, 1), " g/dL"),
      tags$p(strong("LDH at Day 90: "), round(snap90$LDH, 0), " U/L"),
      tags$p(strong("Reticulocytes: "), round(snap90$Reticulocyte_pct, 1), "%"),
      tags$p(strong("Response: "),
             if (snap90$Hemoglobin >= 10) tags$span("Complete Response", style="color:green;font-weight:bold")
             else if (snap90$Hemoglobin >= 8) tags$span("Partial Response", style="color:orange;font-weight:bold")
             else tags$span("Refractory / Insufficient Response", style="color:red;font-weight:bold"))
    )
  })

  # ---- TAB 5: Scenario Comparison ----
  output$compare_hb_plot <- renderPlotly({
    cr <- compare_results()
    if (is.null(cr)) return(NULL)
    colors <- c("A"="#C62828","B"="#1565C0","C"="#2E7D32","D"="#E65100")
    p <- plot_ly()
    for (sc in names(cr)) {
      p <- add_lines(p, data = cr[[sc]], x = ~time, y = ~Hemoglobin,
                     name = paste("Scenario", sc), line = list(color = colors[sc], width = 2))
    }
    p %>%
      add_segments(x=0, xend=max(cr$A$time), y=10, yend=10,
                   line=list(dash="dash",color="gray",width=1), showlegend=FALSE) %>%
      layout(xaxis = list(title = "Time (days)"),
             yaxis = list(title = "Hemoglobin (g/dL)", range = c(4, 16)),
             title = "Hemoglobin Comparison (All Scenarios)")
  })

  output$compare_table <- renderDT({
    cr <- compare_results()
    if (is.null(cr)) return(NULL)
    lapply(names(cr), function(sc) {
      df <- cr[[sc]]
      d90 <- df %>% filter(abs(time - 90) < 1) %>% slice(1)
      d180 <- df %>% filter(abs(time - 180) < 1) %>% slice(1)
      data.frame(
        Scenario  = sc,
        Hb_Day90  = round(d90$Hemoglobin, 1),
        LDH_Day90 = round(d90$LDH, 0),
        Hb_Day180 = round(d180$Hemoglobin, 1),
        CR_Day90  = ifelse(d90$Hemoglobin >= 10, "Yes", "No")
      )
    }) %>% bind_rows() %>%
      datatable(options = list(dom = 't', pageLength = 4), rownames = FALSE)
  })

  output$compare_ldh_plot <- renderPlotly({
    cr <- compare_results()
    if (is.null(cr)) return(NULL)
    colors <- c("A"="#C62828","B"="#1565C0","C"="#2E7D32","D"="#E65100")
    p <- plot_ly()
    for (sc in names(cr)) {
      p <- add_lines(p, data = cr[[sc]], x = ~time, y = ~LDH,
                     name = paste("Scenario", sc), line = list(color = colors[sc]))
    }
    p %>% layout(xaxis = list(title="Time (days)"), yaxis = list(title="LDH (U/L)"),
                 title = "LDH Comparison")
  })

  output$compare_ab_plot <- renderPlotly({
    cr <- compare_results()
    if (is.null(cr)) return(NULL)
    colors <- c("A"="#C62828","B"="#1565C0","C"="#2E7D32","D"="#E65100")
    p <- plot_ly()
    for (sc in names(cr)) {
      p <- add_lines(p, data = cr[[sc]], x = ~time, y = ~AB,
                     name = paste("Scenario", sc), line = list(color = colors[sc]))
    }
    p %>% layout(xaxis = list(title="Time (days)"), yaxis = list(title="Anti-RBC Ab (AU/mL)"),
                 title = "Autoantibody Dynamics")
  })

  # ---- TAB 6: Biomarkers ----
  output$ldh_plot <- renderPlotly({
    df <- sim_result()
    if (is.null(df)) return(NULL)
    p <- ggplot(df, aes(x = time, y = LDH)) +
      geom_line(color = "#C62828", linewidth = 1) +
      geom_hline(yintercept = 250, linetype = "dashed", color = "gray50") +
      annotate("text", x = max(df$time)*0.6, y = 260, label = "ULN 250 U/L", size = 3) +
      labs(x = "Time (days)", y = "LDH (U/L)", title = "Serum LDH") +
      theme_bw(base_size = 11)
    ggplotly(p)
  })

  output$hap_plot <- renderPlotly({
    df <- sim_result()
    if (is.null(df)) return(NULL)
    p <- ggplot(df, aes(x = time, y = Haptoglobin)) +
      geom_line(color = "#1565C0", linewidth = 1) +
      geom_hline(yintercept = 0.5, linetype = "dashed", color = "orange") +
      labs(x = "Time (days)", y = "Haptoglobin (g/L)", title = "Serum Haptoglobin") +
      theme_bw(base_size = 11)
    ggplotly(p)
  })

  output$bilir_plot <- renderPlotly({
    df <- sim_result()
    if (is.null(df)) return(NULL)
    p <- ggplot(df, aes(x = time, y = Bilirubin)) +
      geom_line(color = "#F9A825", linewidth = 1) +
      geom_hline(yintercept = 1.2, linetype = "dashed", color = "gray50") +
      labs(x = "Time (days)", y = "Bilirubin (mg/dL)", title = "Unconjugated Bilirubin") +
      theme_bw(base_size = 11)
    ggplotly(p)
  })

  output$dat_plot <- renderPlotly({
    df <- sim_result()
    if (is.null(df)) return(NULL)
    df <- df %>% mutate(DAT_proxy = AB / (AB + 1.0))
    p <- ggplot(df, aes(x = time, y = DAT_proxy)) +
      geom_line(color = "#6A1B9A", linewidth = 1) +
      geom_hline(yintercept = 0.5, linetype = "dashed", color = "red") +
      scale_y_continuous(limits = c(0, 1), labels = percent) +
      labs(x = "Time (days)", y = "DAT Positivity Score (0-1)",
           title = "Direct Antiglobulin Test (DAT) Proxy\n(1 = strongly positive)") +
      theme_bw(base_size = 11)
    ggplotly(p)
  })

  output$biomarker_ref_table <- renderDT({
    df <- data.frame(
      Biomarker     = c("Hemoglobin", "Hematocrit", "Reticulocytes",
                        "LDH", "Haptoglobin", "Indirect Bilirubin",
                        "Direct Antiglobulin Test (DAT)", "Cold Agglutinin Titer",
                        "C3b on RBC (DAT C3d)"),
      Normal_Range  = c("12-16 g/dL (F), 13.5-17.5 (M)", "36-48% (F), 41-53% (M)",
                        "0.5-2.5%", "<250 U/L",
                        "0.5-2.5 g/L", "0.2-1.2 mg/dL",
                        "Negative", "<1:64",
                        "Negative"),
      AIHA_Value    = c("<10 g/dL", "<30%",
                        ">3-4% (compensatory)", ">300 U/L (often >500)",
                        "<0.5 g/L (often undetectable)", ">2 mg/dL",
                        "Positive IgG (warm), C3d (cold)", ">1:64 (CAD: often >1:512)",
                        "Positive (CAD)"),
      Clinical_Role = c("CR ≥10 g/dL", "Oxygen delivery",
                        "Compensatory erythropoiesis", "Hemolysis intensity",
                        "Hemolysis (intravascular)", "Extravascular hemolysis",
                        "Diagnosis", "CAD diagnosis",
                        "CAD diagnosis")
    )
    datatable(df, options = list(pageLength = 10, dom = 't'), rownames = FALSE)
  })
}

# ============================================================
# LAUNCH APP
# ============================================================
shinyApp(ui = ui, server = server)
