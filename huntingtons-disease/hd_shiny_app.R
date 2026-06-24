################################################################################
# Huntington's Disease QSP — Interactive Shiny Dashboard
# 6 Tabs:
#   Tab 1: Patient Profile & HD Staging
#   Tab 2: Drug Pharmacokinetics
#   Tab 3: PD Key Indices (mHTT, BDNF, Dopamine)
#   Tab 4: Clinical Endpoints (TMS, TFC, cUHDRS)
#   Tab 5: Scenario Comparison (7 treatment arms)
#   Tab 6: Biomarkers & Target Engagement
################################################################################

library(shiny)
library(shinydashboard)
library(ggplot2)
library(dplyr)
library(tidyr)
library(plotly)
library(DT)

# ============================================================================
# SIMPLIFIED ODE SOLVER (Euler for Shiny responsiveness)
# ============================================================================

simulate_hd <- function(
    CAG = 42,
    duration_yrs = 5,
    dt_hr = 24,
    # Treatment flags
    use_TBZ    = FALSE,  dose_TBZ_mgd  = 25,
    use_DTBZ   = FALSE,  dose_DTBZ_mgd = 30,
    use_VBZ    = FALSE,  dose_VBZ_mgd  = 80,
    use_tominersen = FALSE,
    use_branaplam  = FALSE,
    # Disease parameters
    mHTT_agg_rate = 0.0002,
    BDNF_suppression = 0.70,
    neuroinflam_baseline = 1.6
) {
  n_steps <- ceiling(duration_yrs * 365 * 24 / dt_hr)
  times   <- seq(0, duration_yrs * 365 * 24, by = dt_hr)

  # State variables
  mHTT_prot  <- 150.0   # nM
  mHTT_oligo <- 35.0    # nM
  BDNF       <- 4.2     # ng/mL
  dopamine   <- 0.25    # nmol
  MSN        <- 85.0    # % survival
  oxStress   <- 1.8
  inflam     <- 1.6
  TMS        <- 18.0
  TFC        <- 10.5

  # PK states (simplified one-compartment)
  Cp_VMAT2   <- 0.0     # effective VMAT2 inhibitor concentration

  # Parameters
  kprod_mHTT <- 0.008 * (1 + 0.015 * (CAG - 36))
  kdeg_mRNA  <- 0.045; mHTT_mRNA <- 1.0
  ktrans     <- 0.012;  kdeg_protein <- 0.003
  k_agg      <- mHTT_agg_rate
  k_disagg   <- 0.015;  k_UPS <- 0.025; k_UPS_sat <- 120
  k_auto     <- 0.008
  kprod_BDNF <- 1.2;    kdeg_BDNF <- 0.18
  EC50_BDNF_mHTT <- 200; Emax_BDNF_mHTT <- BDNF_suppression
  kprod_DA   <- 0.45;   kdeg_DA <- 1.8
  kdeath_MSN <- 0.0015; EC50_mHTT_death <- 80
  kprod_ROS  <- 0.002;  kdeg_ROS <- 0.08
  kprod_IL1b <- 0.15;   kdeg_IL1b <- 0.04
  kprog_TMS  <- 0.0003; kprog_TFC  <- 0.00008

  # ASO / Branaplam efficacy
  aso_eff  <- ifelse(use_tominersen, 0.74, 0.0)  # 74% mRNA reduction (GENERATION-HD1)
  bran_eff <- ifelse(use_branaplam,  0.50, 0.0)
  total_mRNA_inh <- aso_eff + bran_eff - aso_eff * bran_eff

  # VMAT2 inhibition (simplified: assume steady-state reached within hours)
  VMAT2_inh <- 0.0
  if (use_TBZ)  VMAT2_inh <- max(VMAT2_inh, 0.65 * (dose_TBZ_mgd / 25))  # reference calibration
  if (use_DTBZ) VMAT2_inh <- max(VMAT2_inh, 0.68 * (dose_DTBZ_mgd / 30))
  if (use_VBZ)  VMAT2_inh <- max(VMAT2_inh, 0.75 * (dose_VBZ_mgd  / 80))
  VMAT2_inh <- min(VMAT2_inh, 0.92)  # cap at 92%

  # Pre-allocate output
  out <- data.frame(
    time      = times,
    year      = times / (365 * 24),
    mHTT_prot = NA, mHTT_oligo = NA,
    BDNF      = NA, dopamine   = NA,
    MSN       = NA, oxStress   = NA,
    inflam    = NA, TMS        = NA,
    TFC       = NA,
    mHTT_total = NA,
    chorea_red_pct = 100 * VMAT2_inh,
    cUHDRS    = NA
  )

  for (i in seq_along(times)) {
    out$mHTT_prot[i]  <- mHTT_prot
    out$mHTT_oligo[i] <- mHTT_oligo
    out$BDNF[i]       <- BDNF
    out$dopamine[i]   <- dopamine
    out$MSN[i]        <- MSN
    out$oxStress[i]   <- oxStress
    out$inflam[i]     <- inflam
    out$TMS[i]        <- max(0, TMS)
    out$TFC[i]        <- max(0, TFC)
    out$mHTT_total[i] <- mHTT_prot + mHTT_oligo
    out$cUHDRS[i]     <- max(0, TFC + (25 - TMS/5)/2)

    if (i == length(times)) break

    # Derivatives
    mHTT_mRNA_new <- mHTT_mRNA + dt_hr * (kprod_mHTT * (1 - total_mRNA_inh) - kdeg_mRNA * mHTT_mRNA)
    mHTT_mRNA <- max(0, mHTT_mRNA_new)

    UPS_eff  <- k_UPS * k_UPS_sat / (k_UPS_sat + mHTT_prot)
    dmHTT_p  <- ktrans * mHTT_mRNA - (UPS_eff + k_auto) * mHTT_prot -
                k_agg * mHTT_prot^2 + k_disagg * mHTT_oligo
    dmHTT_o  <- k_agg * mHTT_prot^2 - k_disagg * mHTT_oligo - k_auto * mHTT_oligo

    BDNF_mHTT_sup <- Emax_BDNF_mHTT * mHTT_oligo / (EC50_BDNF_mHTT + mHTT_oligo)
    dBDNF    <- kprod_BDNF * (1 - BDNF_mHTT_sup) + 0.15 * BDNF * (MSN/100) - kdeg_BDNF * BDNF

    VMAT2_act <- (1 - VMAT2_inh)
    dDA      <- kprod_DA * VMAT2_act * (MSN/100) - kdeg_DA * dopamine

    BDNF_prot   <- BDNF / (3.0 + BDNF)
    mHTT_kill   <- kdeath_MSN * mHTT_oligo / (EC50_mHTT_death + mHTT_oligo)
    net_death   <- (mHTT_kill + 0.0008 * oxStress + 0.0005 * inflam) * (1 - 0.6 * BDNF_prot)
    dMSN        <- -net_death * MSN

    dROS   <- kprod_ROS * (mHTT_oligo/30 + mHTT_prot/150)/2 - kdeg_ROS * oxStress
    dIL1b  <- kprod_IL1b * (mHTT_oligo/50) * (1 - MSN/100 + 0.2) - kdeg_IL1b * inflam

    DA_excess   <- dopamine / (0.25 * (MSN/100 + 0.1))
    chorea_drive <- kprog_TMS * DA_excess * (100 - MSN) / 10
    chorea_tx    <- VMAT2_inh * 0.55 * (TMS / 50)
    dmt_TMS      <- 0.0
    if (use_tominersen || use_branaplam) dmt_TMS <- 0.3 * kprog_TMS
    dTMS  <- chorea_drive - chorea_tx - dmt_TMS

    TFC_decline  <- kprog_TFC * (100 - MSN) / 20
    TFC_stabilize <- ifelse(use_tominersen || use_branaplam, 0.3 * TFC_decline, 0)
    dTFC  <- -TFC_decline + TFC_stabilize

    # Update states
    mHTT_prot  <- max(0, mHTT_prot  + dt_hr * dmHTT_p)
    mHTT_oligo <- max(0, mHTT_oligo + dt_hr * dmHTT_o)
    BDNF       <- max(0, BDNF       + dt_hr * dBDNF)
    dopamine   <- max(0, dopamine   + dt_hr * dDA)
    MSN        <- max(0, min(100, MSN + dt_hr * dMSN))
    oxStress   <- max(0, oxStress   + dt_hr * dROS)
    inflam     <- max(0, inflam     + dt_hr * dIL1b)
    TMS        <- max(0, TMS        + dt_hr * dTMS)
    TFC        <- max(0, TFC        + dt_hr * dTFC)
  }

  return(out)
}

# ============================================================================
# HD STAGING FUNCTION
# ============================================================================

hd_stage <- function(TFC) {
  dplyr::case_when(
    TFC >= 13             ~ "Pre-manifest",
    TFC >= 11 & TFC < 13 ~ "Stage I (mild)",
    TFC >= 7  & TFC < 11 ~ "Stage II (moderate)",
    TFC >= 3  & TFC < 7  ~ "Stage III (moderate-severe)",
    TFC >= 1  & TFC < 3  ~ "Stage IV (severe)",
    TRUE                  ~ "Stage V (end-stage)"
  )
}

cap_score <- function(CAG, age) { pmax(0, (CAG - 33.66) * age) }

# ============================================================================
# UI
# ============================================================================

ui <- dashboardPage(
  skin = "blue",
  dashboardHeader(title = "Huntington's Disease QSP"),

  dashboardSidebar(
    sidebarMenu(
      menuItem("Patient Profile",     tabName = "patient",  icon = icon("user-md")),
      menuItem("Pharmacokinetics",    tabName = "pk",       icon = icon("flask")),
      menuItem("PD Key Indices",      tabName = "pd",       icon = icon("dna")),
      menuItem("Clinical Endpoints",  tabName = "clinical", icon = icon("chart-line")),
      menuItem("Scenario Comparison", tabName = "scenario", icon = icon("balance-scale")),
      menuItem("Biomarkers",          tabName = "biomarker",icon = icon("vial"))
    )
  ),

  dashboardBody(
    tags$head(tags$style(HTML("
      .box { border-radius: 8px; }
      .info-box { border-radius: 8px; }
      .skin-blue .main-header .logo { background-color: #1a237e; }
      .skin-blue .main-header .navbar { background-color: #283593; }
    "))),

    tabItems(

      # ==================================================================
      # TAB 1: PATIENT PROFILE & HD STAGING
      # ==================================================================
      tabItem(tabName = "patient",
        fluidRow(
          box(title = "Patient Parameters", status = "primary", solidHeader = TRUE, width = 4,
            sliderInput("CAG", "CAG Repeat Length", min = 36, max = 70, value = 42, step = 1),
            sliderInput("age", "Patient Age (years)", min = 20, max = 70, value = 42, step = 1),
            numericInput("TMS_baseline", "Baseline UHDRS-TMS", value = 18, min = 0, max = 100),
            numericInput("TFC_baseline", "Baseline TFC", value = 10.5, min = 0, max = 13),
            selectInput("sex", "Sex", choices = c("Male", "Female")),
            hr(),
            h4("Inheritance Pattern"),
            selectInput("inheritance", "Family History",
                        choices = c("Paternal", "Maternal", "Sporadic")),
            actionButton("update_patient", "Update Profile", class = "btn-primary btn-block")
          ),
          box(title = "HD Staging & Risk Assessment", status = "warning",
              solidHeader = TRUE, width = 8,
            fluidRow(
              infoBoxOutput("stage_box", width = 4),
              infoBoxOutput("cap_box",   width = 4),
              infoBoxOutput("onset_box", width = 4)
            ),
            hr(),
            h4("CAP Score Interpretation"),
            tableOutput("cap_table"),
            hr(),
            h4("Disease Natural History — 10-Year Projection"),
            plotlyOutput("natural_history_plot", height = "320px")
          )
        ),
        fluidRow(
          box(title = "Pathophysiology Overview", status = "info",
              solidHeader = TRUE, width = 12,
            column(6,
              h4("Core Mechanism"),
              tags$ul(
                tags$li(strong("CAG Repeat"), " expansion in HTT gene (Ch 4p16.3) → polyglutamine (polyQ) tract"),
                tags$li(strong("mHTT Protein"), " misfolding → oligomerization → inclusion body formation"),
                tags$li(strong("BDNF Deficit"), " (REST/NRSF sequestr.) → loss of TrkB-mediated MSN survival"),
                tags$li(strong("Striatal MSN degeneration"), " (D2>D1) → chorea (early), bradykinesia (late)"),
                tags$li(strong("Excitotoxicity"), ": eNMDAR sensitization → Ca²⁺ overload → mitochondrial failure")
              )
            ),
            column(6,
              h4("Key Clinical Scales"),
              tags$ul(
                tags$li(strong("UHDRS-TMS:"), " 0–124; chorea/rigidity/dystonia subscores"),
                tags$li(strong("TFC:"), " 13→0; work, finances, domestic, ADL, self-care"),
                tags$li(strong("cUHDRS:"), " composite of TFC + TMS + VFC + SDMT (TRACK-HD)"),
                tags$li(strong("CAP Score:"), " (CAG−33.66)×Age; predicts disease burden"),
                tags$li(strong("CSF NfL:"), " neurofilament light chain — neurodegeneration rate"),
                tags$li(strong("CSF mHTT:"), " target engagement biomarker for HTT-lowering therapies")
              )
            )
          )
        )
      ),

      # ==================================================================
      # TAB 2: PHARMACOKINETICS
      # ==================================================================
      tabItem(tabName = "pk",
        fluidRow(
          box(title = "PK Parameters", status = "primary", solidHeader = TRUE, width = 3,
            h4("VMAT2 Inhibitors"),
            selectInput("pk_drug", "Select Drug",
                        choices = c("Tetrabenazine (TBZ)", "Deutetrabenazine (DTBZ)",
                                    "Valbenazine (VBZ)")),
            numericInput("pk_dose", "Dose (mg)", value = 25, min = 6, max = 150),
            selectInput("pk_freq", "Frequency",
                        choices = c("Once daily (QD)", "Twice daily (BID)", "Three times daily (TID)")),
            numericInput("pk_dur_days", "Duration (days)", value = 28, min = 1, max = 365),
            hr(),
            h4("ASO / Splicing Modifiers"),
            checkboxInput("show_aso_pk", "Show Tominersen PK (IT)", FALSE),
            numericInput("tominersen_dose_pk", "Tominersen Dose (mg IT)", 120, min = 30, max = 200),
            actionButton("run_pk", "Run PK Simulation", class = "btn-info btn-block")
          ),
          box(title = "PK Simulation", status = "info", solidHeader = TRUE, width = 9,
            plotlyOutput("pk_plasma_plot", height = "280px"),
            hr(),
            plotlyOutput("pk_brain_plot", height = "280px")
          )
        ),
        fluidRow(
          box(title = "PK Parameter Summary", status = "success", solidHeader = TRUE, width = 12,
            DT::dataTableOutput("pk_param_table")
          )
        )
      ),

      # ==================================================================
      # TAB 3: PD KEY INDICES
      # ==================================================================
      tabItem(tabName = "pd",
        fluidRow(
          box(title = "Disease & Treatment Settings", status = "primary",
              solidHeader = TRUE, width = 3,
            sliderInput("pd_CAG", "CAG Repeat Length", 36, 70, 42, step = 1),
            sliderInput("pd_dur", "Duration (years)", 1, 10, 5, step = 0.5),
            hr(),
            checkboxInput("pd_tominersen", "Tominersen (ASO, IT)", FALSE),
            checkboxInput("pd_branaplam", "Branaplam (oral)", FALSE),
            hr(),
            sliderInput("pd_agg_rate", "mHTT Aggregation Rate",
                        min = 0.0001, max = 0.0005, value = 0.0002, step = 0.00005),
            sliderInput("pd_BDNF_supp", "BDNF Suppression by mHTT",
                        min = 0.3, max = 0.9, value = 0.7, step = 0.05),
            actionButton("run_pd", "Run PD Simulation", class = "btn-warning btn-block")
          ),
          tabBox(title = "PD Dynamics", width = 9,
            tabPanel("mHTT Cascade",
              plotlyOutput("mHTT_plot", height = "320px"),
              plotlyOutput("mHTT_oligo_plot", height = "280px")
            ),
            tabPanel("BDNF & Survival",
              plotlyOutput("BDNF_plot", height = "280px"),
              plotlyOutput("MSN_plot", height = "280px")
            ),
            tabPanel("Dopamine & Oxidative",
              plotlyOutput("DA_plot", height = "280px"),
              plotlyOutput("ROS_inflam_plot", height = "280px")
            )
          )
        )
      ),

      # ==================================================================
      # TAB 4: CLINICAL ENDPOINTS
      # ==================================================================
      tabItem(tabName = "clinical",
        fluidRow(
          box(title = "Simulation Settings", status = "primary",
              solidHeader = TRUE, width = 3,
            sliderInput("clin_CAG", "CAG Repeat Length", 36, 70, 42, step = 1),
            sliderInput("clin_dur", "Duration (years)", 1, 10, 5, step = 0.5),
            hr(),
            h4("Treatment"),
            checkboxInput("clin_TBZ",   "Tetrabenazine 25 mg/d",  FALSE),
            checkboxInput("clin_DTBZ",  "Deutetrabenazine 30 mg/d", FALSE),
            checkboxInput("clin_VBZ",   "Valbenazine 80 mg/d",    FALSE),
            checkboxInput("clin_ASO",   "Tominersen (IT, Q8W)",   FALSE),
            checkboxInput("clin_bran",  "Branaplam (QW)",          FALSE),
            actionButton("run_clin", "Run Simulation", class = "btn-success btn-block")
          ),
          box(title = "Clinical Endpoint Trajectories", status = "success",
              solidHeader = TRUE, width = 9,
            plotlyOutput("TMS_plot",  height = "270px"),
            hr(),
            plotlyOutput("TFC_plot",  height = "270px")
          )
        ),
        fluidRow(
          box(title = "cUHDRS Composite Score", status = "info", solidHeader = TRUE, width = 6,
            plotlyOutput("cUHDRS_plot", height = "280px")
          ),
          box(title = "Clinical Summary at Endpoint", status = "warning",
              solidHeader = TRUE, width = 6,
            tableOutput("clin_summary_table")
          )
        )
      ),

      # ==================================================================
      # TAB 5: SCENARIO COMPARISON
      # ==================================================================
      tabItem(tabName = "scenario",
        fluidRow(
          box(title = "Scenario Configuration", status = "primary",
              solidHeader = TRUE, width = 3,
            sliderInput("sc_CAG", "CAG Length", 36, 70, 42, step = 1),
            sliderInput("sc_dur", "Duration (years)", 1, 10, 5, step = 0.5),
            hr(),
            checkboxGroupInput("sc_scenarios", "Scenarios to Include",
              choices = c(
                "Natural History"          = "NH",
                "TBZ 25 mg/d"             = "TBZ",
                "DTBZ 30 mg/d"            = "DTBZ",
                "VBZ 80 mg/d"             = "VBZ",
                "Tominersen Q8W"          = "ASO",
                "Branaplam QW"            = "Bran",
                "DTBZ + Tominersen"       = "Combo"
              ),
              selected = c("NH","TBZ","DTBZ","VBZ","ASO","Bran","Combo")
            ),
            actionButton("run_scenario", "Compare Scenarios", class = "btn-danger btn-block")
          ),
          box(title = "Scenario Comparison — TMS & TFC", status = "danger",
              solidHeader = TRUE, width = 9,
            plotlyOutput("sc_TMS_plot", height = "280px"),
            hr(),
            plotlyOutput("sc_TFC_plot", height = "280px")
          )
        ),
        fluidRow(
          box(title = "5-Year Endpoint Summary Table", status = "success",
              solidHeader = TRUE, width = 12,
            DT::dataTableOutput("sc_summary_table")
          )
        )
      ),

      # ==================================================================
      # TAB 6: BIOMARKERS
      # ==================================================================
      tabItem(tabName = "biomarker",
        fluidRow(
          box(title = "Biomarker Settings", status = "primary",
              solidHeader = TRUE, width = 3,
            sliderInput("bio_CAG",   "CAG Length",       36, 70, 42, step = 1),
            sliderInput("bio_dur",   "Duration (years)",  1, 10, 5, step = 0.5),
            checkboxInput("bio_ASO",  "Tominersen (IT)", FALSE),
            checkboxInput("bio_bran", "Branaplam",        FALSE),
            hr(),
            h4("CSF Biomarker Reference Ranges"),
            tags$ul(
              tags$li("NfL (CSF): HD manifest ~3000–5000 pg/mL vs. healthy <500"),
              tags$li("mHTT (CSF): 0 in healthy; elevated in manifest HD"),
              tags$li("BDNF (CSF): healthy ~6.5; HD ~3–4 ng/mL"),
              tags$li("Oligomer ratio: >25% signals high toxicity burden")
            ),
            actionButton("run_bio", "Run Biomarker Simulation", class = "btn-info btn-block")
          ),
          box(title = "Biomarker Trajectory Panel", status = "info",
              solidHeader = TRUE, width = 9,
            plotlyOutput("bio_NfL_plot",    height = "240px"),
            plotlyOutput("bio_mHTT_plot",   height = "240px"),
            plotlyOutput("bio_oligo_plot",  height = "200px")
          )
        ),
        fluidRow(
          box(title = "Target Engagement (tominersen / branaplam)", status = "warning",
              solidHeader = TRUE, width = 6,
            plotlyOutput("bio_TE_plot", height = "300px")
          ),
          box(title = "Biomarker-Clinical Correlation", status = "success",
              solidHeader = TRUE, width = 6,
            plotlyOutput("bio_corr_plot", height = "300px"),
            p("Correlation between CSF NfL and UHDRS-TMS (simulated).")
          )
        )
      )
    )
  )
)

# ============================================================================
# SERVER
# ============================================================================

server <- function(input, output, session) {

  # ── Reactive: Natural history simulation (Tab 1) ──
  patient_sim <- reactive({
    input$update_patient
    isolate({
      simulate_hd(CAG = input$CAG, duration_yrs = 10, dt_hr = 48)
    })
  })

  output$stage_box <- renderInfoBox({
    df  <- patient_sim()
    stg <- hd_stage(tail(df$TFC, 1))
    infoBox("Current Stage", stg, icon = icon("heartbeat"),
            color = "red", fill = TRUE)
  })

  output$cap_box <- renderInfoBox({
    cap <- round(cap_score(input$CAG, input$age), 0)
    col <- ifelse(cap < 200, "green", ifelse(cap < 400, "yellow", "red"))
    infoBox("CAP Score", cap,
            subtitle = ifelse(cap < 200, "Low risk", ifelse(cap < 400, "Moderate risk", "High risk")),
            icon = icon("calculator"), color = col, fill = TRUE)
  })

  output$onset_box <- renderInfoBox({
    # Approximate age at onset: (CAG - 34.98) × age relationship
    onset_age <- round(21.54 + 9.556e7 / (input$CAG - 35.55)^4.1, 0)
    infoBox("Est. Motor Onset", paste(onset_age, "yrs"),
            icon = icon("clock"), color = "orange", fill = TRUE)
  })

  output$cap_table <- renderTable({
    data.frame(
      `CAP Range` = c("< 200", "200–400", "400–600", "> 600"),
      `HD Stage`  = c("Pre-manifest / Low risk", "Stage I–II onset likely",
                      "Manifest HD, moderate", "Manifest HD, severe"),
      `5-yr TMS change` = c("+1–3", "+3–8", "+8–18", "+15–30")
    )
  })

  output$natural_history_plot <- renderPlotly({
    df <- patient_sim() %>%
      select(year, TMS, TFC, MSN) %>%
      pivot_longer(-year, names_to = "variable", values_to = "value") %>%
      mutate(variable = recode(variable,
        "TMS" = "UHDRS-TMS (0–124)",
        "TFC" = "TFC (0–13)",
        "MSN" = "MSN Survival (%)"
      ))
    p <- ggplot(df, aes(x = year, y = value, color = variable)) +
      geom_line(linewidth = 1.0) +
      facet_wrap(~variable, scales = "free_y", ncol = 3) +
      theme_bw(base_size = 11) +
      labs(x = "Year", y = NULL, color = NULL) +
      theme(legend.position = "none")
    ggplotly(p) %>% layout(showlegend = FALSE)
  })

  # ── Tab 2: PK ──
  pk_data <- reactive({
    input$run_pk
    isolate({
      drug <- input$pk_drug
      dose <- input$pk_dose
      n_days <- input$pk_dur_days
      times <- seq(0, n_days * 24, by = 0.5)

      # Simplified one-compartment PK
      pk_params <- list(
        "Tetrabenazine (TBZ)"    = list(ka=0.8, F=0.20, CL=58, Vc=84, Kp=0.35, ii=8,  name="TBZ"),
        "Deutetrabenazine (DTBZ)"= list(ka=0.65,F=0.82, CL=14, Vc=180,Kp=0.42, ii=12, name="DTBZ"),
        "Valbenazine (VBZ)"      = list(ka=0.45,F=0.49, CL=7.8,Vc=280,Kp=0.68, ii=24, name="VBZ")
      )
      pp <- pk_params[[drug]]
      ii <- pp$ii

      dose_times <- seq(0, max(times) - 0.01, by = ii)
      Cp <- numeric(length(times))
      for (dt in dose_times) {
        tpost <- pmax(0, times - dt)
        Cp <- Cp + (dose * pp$F * pp$ka) / (pp$Vc * (pp$ka - pp$CL/pp$Vc)) *
          (exp(-pp$CL/pp$Vc * tpost) - exp(-pp$ka * tpost))
      }
      Cp_brain <- Cp * pp$Kp

      data.frame(time_hr = times, day = times/24,
                 Cp_plasma = pmax(0, Cp),
                 Cp_brain  = pmax(0, Cp_brain),
                 drug = pp$name)
    })
  })

  output$pk_plasma_plot <- renderPlotly({
    df <- pk_data()
    p <- ggplot(df, aes(x = day, y = Cp_plasma)) +
      geom_line(color = "#2980b9", linewidth = 1.1) +
      labs(title = paste("Plasma Concentration —", unique(df$drug)),
           x = "Day", y = "Cp (mg/L)") + theme_bw()
    ggplotly(p)
  })

  output$pk_brain_plot <- renderPlotly({
    df <- pk_data()
    p <- ggplot(df, aes(x = day, y = Cp_brain)) +
      geom_line(color = "#8e44ad", linewidth = 1.1) +
      labs(title = "Brain Concentration (Kp corrected)",
           x = "Day", y = "Brain Conc. (mg/L)") + theme_bw()
    ggplotly(p)
  })

  output$pk_param_table <- DT::renderDataTable({
    df <- data.frame(
      Drug          = c("Tetrabenazine", "Deutetrabenazine", "Valbenazine",
                        "Tominersen (IT)", "Riluzole", "Memantine"),
      Dose          = c("25–100 mg/d (TID)", "12–48 mg/d (BID)", "40–80 mg/d (QD)",
                        "120 mg Q8W (IT)", "50 mg BID", "10–20 mg BID"),
      Bioavailability = c("~20%", "~82%", "~49%", "N/A (IT)", "~60%", "~100%"),
      `T½ (hr)`    = c("2–4 (HTBZ 4–8)", "9–10", "15–22", "~2000 hr", "12", "60–80"),
      `Mechanism`  = c("VMAT2 inhibitor", "VMAT2 inhib (d-KIE)", "VMAT2 inhib (QD)",
                        "ASO → RNase H1", "Na+ chan/Glu block", "NMDAR antagonist"),
      `Clinical trial` = c("TETRA-HD (NEJM 2008)", "FIRST-HD (NEJM 2016)",
                            "KINECT-HD (NEJM 2023)", "GENERATION-HD1 (NEJM 2022)",
                            "Small studies", "Observational")
    )
    DT::datatable(df, options = list(pageLength = 10, dom = 't'), rownames = FALSE)
  })

  # ── Tab 3: PD Indices ──
  pd_sim <- reactive({
    input$run_pd
    isolate({
      simulate_hd(
        CAG = input$pd_CAG, duration_yrs = input$pd_dur, dt_hr = 24,
        use_tominersen = input$pd_tominersen, use_branaplam = input$pd_branaplam,
        mHTT_agg_rate  = input$pd_agg_rate,
        BDNF_suppression = input$pd_BDNF_supp
      )
    })
  })

  output$mHTT_plot <- renderPlotly({
    df <- pd_sim() %>% select(year, mHTT_prot, mHTT_total) %>%
      pivot_longer(-year, names_to = "var", values_to = "nM") %>%
      mutate(var = recode(var, mHTT_prot = "mHTT Soluble", mHTT_total = "mHTT Total"))
    p <- ggplot(df, aes(x=year, y=nM, color=var)) + geom_line(linewidth=1) +
      labs(title="mHTT Protein Dynamics", x="Year", y="Concentration (nM)", color=NULL) +
      theme_bw()
    ggplotly(p)
  })

  output$mHTT_oligo_plot <- renderPlotly({
    df <- pd_sim()
    p <- ggplot(df, aes(x=year, y=mHTT_oligo)) + geom_line(color="#c0392b", linewidth=1) +
      labs(title="mHTT Oligomers (toxic species)", x="Year", y="Oligomer (nM)") + theme_bw()
    ggplotly(p)
  })

  output$BDNF_plot <- renderPlotly({
    df <- pd_sim()
    p <- ggplot(df, aes(x=year, y=BDNF)) + geom_line(color="#2980b9", linewidth=1) +
      geom_hline(yintercept=6.5, linetype="dashed", color="darkgreen") +
      labs(title="BDNF Level (CSF proxy)", x="Year", y="BDNF (ng/mL)") + theme_bw()
    ggplotly(p)
  })

  output$MSN_plot <- renderPlotly({
    df <- pd_sim()
    p <- ggplot(df, aes(x=year, y=MSN)) + geom_line(color="#27ae60", linewidth=1) +
      labs(title="MSN Survival (%)", x="Year", y="MSN Survival (%)") + theme_bw()
    ggplotly(p)
  })

  output$DA_plot <- renderPlotly({
    df <- pd_sim()
    p <- ggplot(df, aes(x=year, y=dopamine)) + geom_line(color="#8e44ad", linewidth=1) +
      labs(title="Synaptic Dopamine (nmol)", x="Year", y="Dopamine (nmol)") + theme_bw()
    ggplotly(p)
  })

  output$ROS_inflam_plot <- renderPlotly({
    df <- pd_sim() %>% select(year, oxStress, inflam) %>%
      pivot_longer(-year, names_to="var", values_to="value") %>%
      mutate(var = recode(var, oxStress="Oxidative Stress Index", inflam="Neuroinflammation Index"))
    p <- ggplot(df, aes(x=year, y=value, color=var)) + geom_line(linewidth=1) +
      labs(title="Oxidative Stress & Neuroinflammation", x="Year", y="Index (AU)", color=NULL) +
      theme_bw()
    ggplotly(p)
  })

  # ── Tab 4: Clinical ──
  clin_sim <- reactive({
    input$run_clin
    isolate({
      simulate_hd(
        CAG = input$clin_CAG, duration_yrs = input$clin_dur,
        use_TBZ  = input$clin_TBZ,  dose_TBZ_mgd = 25,
        use_DTBZ = input$clin_DTBZ, dose_DTBZ_mgd = 30,
        use_VBZ  = input$clin_VBZ,  dose_VBZ_mgd = 80,
        use_tominersen = input$clin_ASO,
        use_branaplam  = input$clin_bran
      )
    })
  })

  output$TMS_plot <- renderPlotly({
    df <- clin_sim()
    p <- ggplot(df, aes(x=year, y=TMS)) + geom_line(color="#e74c3c", linewidth=1.2) +
      labs(title="UHDRS Total Motor Score (TMS)", x="Year", y="TMS (0–124)") + theme_bw()
    ggplotly(p)
  })

  output$TFC_plot <- renderPlotly({
    df <- clin_sim()
    p <- ggplot(df, aes(x=year, y=TFC)) + geom_line(color="#2980b9", linewidth=1.2) +
      labs(title="Total Functional Capacity (TFC)", x="Year", y="TFC (0–13)") + theme_bw()
    ggplotly(p)
  })

  output$cUHDRS_plot <- renderPlotly({
    df <- clin_sim()
    p <- ggplot(df, aes(x=year, y=cUHDRS)) + geom_line(color="#27ae60", linewidth=1.2) +
      labs(title="cUHDRS Composite Score", x="Year", y="cUHDRS") + theme_bw()
    ggplotly(p)
  })

  output$clin_summary_table <- renderTable({
    df <- clin_sim()
    n  <- nrow(df)
    data.frame(
      Endpoint    = c("UHDRS-TMS", "TFC", "MSN Survival (%)", "BDNF (ng/mL)",
                      "mHTT Total (nM)", "cUHDRS", "HD Stage"),
      Baseline    = c(round(df$TMS[1],1), round(df$TFC[1],1), round(df$MSN[1],1),
                      round(df$BDNF[1],2), round(df$mHTT_total[1],0),
                      round(df$cUHDRS[1],1), hd_stage(df$TFC[1])),
      Endpoint_   = c(round(df$TMS[n],1), round(df$TFC[n],1), round(df$MSN[n],1),
                      round(df$BDNF[n],2), round(df$mHTT_total[n],0),
                      round(df$cUHDRS[n],1), hd_stage(df$TFC[n])),
      Change      = c(round(df$TMS[n]-df$TMS[1],1), round(df$TFC[n]-df$TFC[1],1),
                      round(df$MSN[n]-df$MSN[1],1), round(df$BDNF[n]-df$BDNF[1],2),
                      round(df$mHTT_total[n]-df$mHTT_total[1],0),
                      round(df$cUHDRS[n]-df$cUHDRS[1],1), "")
    )
  }, colnames = TRUE)

  # ── Tab 5: Scenarios ──
  sc_sims <- reactive({
    input$run_scenario
    isolate({
      sc_map <- list(
        NH    = list(use_TBZ=F, use_DTBZ=F, use_VBZ=F, use_tominersen=F, use_branaplam=F),
        TBZ   = list(use_TBZ=T, use_DTBZ=F, use_VBZ=F, use_tominersen=F, use_branaplam=F),
        DTBZ  = list(use_TBZ=F, use_DTBZ=T, use_VBZ=F, use_tominersen=F, use_branaplam=F),
        VBZ   = list(use_TBZ=F, use_DTBZ=F, use_VBZ=T, use_tominersen=F, use_branaplam=F),
        ASO   = list(use_TBZ=F, use_DTBZ=F, use_VBZ=F, use_tominersen=T, use_branaplam=F),
        Bran  = list(use_TBZ=F, use_DTBZ=F, use_VBZ=F, use_tominersen=F, use_branaplam=T),
        Combo = list(use_TBZ=F, use_DTBZ=T, use_VBZ=F, use_tominersen=T, use_branaplam=F)
      )
      sc_labels <- c(NH="Natural History", TBZ="TBZ 25 mg/d",
                     DTBZ="DTBZ 30 mg/d", VBZ="VBZ 80 mg/d",
                     ASO="Tominersen Q8W", Bran="Branaplam QW",
                     Combo="DTBZ+Tominersen")

      selected <- input$sc_scenarios
      results  <- lapply(selected, function(sc) {
        args <- c(list(CAG=input$sc_CAG, duration_yrs=input$sc_dur), sc_map[[sc]])
        df <- do.call(simulate_hd, args)
        df$scenario <- sc_labels[[sc]]
        df
      })
      bind_rows(results)
    })
  })

  sc_colors <- c(
    "Natural History"   = "#666666",
    "TBZ 25 mg/d"       = "#e74c3c",
    "DTBZ 30 mg/d"      = "#e67e22",
    "VBZ 80 mg/d"       = "#f39c12",
    "Tominersen Q8W"    = "#2980b9",
    "Branaplam QW"      = "#8e44ad",
    "DTBZ+Tominersen"   = "#27ae60"
  )

  output$sc_TMS_plot <- renderPlotly({
    df <- sc_sims()
    p <- ggplot(df, aes(x=year, y=TMS, color=scenario)) + geom_line(linewidth=1) +
      scale_color_manual(values=sc_colors) +
      labs(title="UHDRS-TMS — Treatment Scenario Comparison",
           x="Year", y="TMS", color=NULL) + theme_bw()
    ggplotly(p) %>% layout(legend=list(orientation="h", y=-0.2))
  })

  output$sc_TFC_plot <- renderPlotly({
    df <- sc_sims()
    p <- ggplot(df, aes(x=year, y=TFC, color=scenario)) + geom_line(linewidth=1) +
      scale_color_manual(values=sc_colors) +
      labs(title="TFC — Treatment Scenario Comparison",
           x="Year", y="TFC (0–13)", color=NULL) + theme_bw()
    ggplotly(p) %>% layout(legend=list(orientation="h", y=-0.2))
  })

  output$sc_summary_table <- DT::renderDataTable({
    df <- sc_sims()
    summary_df <- df %>% group_by(scenario) %>%
      summarise(
        TMS_start = round(first(TMS), 1),
        TMS_end   = round(last(TMS), 1),
        TMS_delta = round(last(TMS) - first(TMS), 1),
        TFC_start = round(first(TFC), 1),
        TFC_end   = round(last(TFC), 1),
        TFC_delta = round(last(TFC) - first(TFC), 1),
        MSN_pct   = round(last(MSN), 1),
        mHTT_end  = round(last(mHTT_total), 0),
        Chorea_red = round(mean(chorea_red_pct), 0),
        .groups   = "drop"
      ) %>%
      arrange(TMS_delta)
    DT::datatable(summary_df, options=list(dom="t", pageLength=10), rownames=FALSE) %>%
      DT::formatStyle("TMS_delta",
        backgroundColor = DT::styleInterval(c(0, 5, 15),
                                             c("#c8e6c9","#fff9c4","#ffcc80","#ef9a9a")))
  })

  # ── Tab 6: Biomarkers ──
  bio_sim <- reactive({
    input$run_bio
    isolate({
      simulate_hd(CAG=input$bio_CAG, duration_yrs=input$bio_dur,
                  use_tominersen=input$bio_ASO, use_branaplam=input$bio_bran)
    })
  })

  output$bio_NfL_plot <- renderPlotly({
    df <- bio_sim() %>%
      mutate(NfL = 2000 + 80 * (100 - MSN))  # proxy: NfL ~ MSN loss
    p <- ggplot(df, aes(x=year, y=NfL)) + geom_line(color="#c0392b", linewidth=1) +
      geom_hline(yintercept=500, linetype="dashed", color="darkgreen") +
      annotate("text", x=0.2, y=550, label="Healthy NfL threshold", size=3, color="darkgreen") +
      labs(title="CSF Neurofilament Light Chain (NfL proxy)", x="Year", y="NfL (pg/mL)") +
      theme_bw()
    ggplotly(p)
  })

  output$bio_mHTT_plot <- renderPlotly({
    df <- bio_sim()
    p <- ggplot(df, aes(x=year, y=mHTT_total)) + geom_line(color="#e74c3c", linewidth=1) +
      labs(title="CSF mHTT (Target Engagement Biomarker)", x="Year", y="mHTT (nM)") +
      theme_bw()
    ggplotly(p)
  })

  output$bio_oligo_plot <- renderPlotly({
    df <- bio_sim() %>%
      mutate(oligo_pct = 100 * mHTT_oligo / (mHTT_total + 0.01))
    p <- ggplot(df, aes(x=year, y=oligo_pct)) + geom_line(color="#8e44ad", linewidth=1) +
      geom_hline(yintercept=25, linetype="dashed", color="orange") +
      annotate("text", x=0.3, y=26.5, label="High toxicity threshold (25%)", size=3, color="orange") +
      labs(title="mHTT Oligomer Fraction (%)", x="Year", y="Oligomer %") + theme_bw()
    ggplotly(p)
  })

  output$bio_TE_plot <- renderPlotly({
    df <- bio_sim()
    df_NH <- simulate_hd(CAG=input$bio_CAG, duration_yrs=input$bio_dur,
                          use_tominersen=FALSE, use_branaplam=FALSE)
    combined <- bind_rows(
      df    %>% mutate(scenario="With Treatment"),
      df_NH %>% mutate(scenario="Natural History")
    )
    reduction <- combined %>%
      group_by(scenario) %>%
      mutate(mHTT_pct_of_baseline = 100 * mHTT_total / first(df_NH$mHTT_total))
    p <- ggplot(reduction, aes(x=year, y=mHTT_pct_of_baseline, color=scenario)) +
      geom_line(linewidth=1) +
      scale_color_manual(values=c("Natural History"="#666666","With Treatment"="#2980b9")) +
      labs(title="mHTT Reduction vs Natural History (%)",
           x="Year", y="mHTT (% of NH baseline)", color=NULL) + theme_bw()
    ggplotly(p)
  })

  output$bio_corr_plot <- renderPlotly({
    df <- bio_sim() %>%
      mutate(NfL = 2000 + 80 * (100 - MSN)) %>%
      filter(year > 0.05)
    p <- ggplot(df, aes(x=NfL, y=TMS)) +
      geom_point(alpha=0.4, color="#e74c3c", size=1.5) +
      geom_smooth(method="lm", color="#2980b9", se=TRUE) +
      labs(title="CSF NfL vs UHDRS-TMS (Pearson r~0.85 in literature)",
           x="CSF NfL (pg/mL)", y="UHDRS-TMS") + theme_bw()
    ggplotly(p)
  })
}

# ============================================================================
# RUN APP
# ============================================================================
shinyApp(ui = ui, server = server)
