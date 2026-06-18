## ============================================================
## Hypertrophic Cardiomyopathy (HCM) QSP — Interactive Shiny App
## 7 Tabs: Patient Profile · Drug PK · Cardiac Mechanics ·
##         Hypertrophy & Fibrosis · Clinical Endpoints ·
##         Scenario Comparison · Risk Assessment
## ============================================================

library(shiny)
library(shinydashboard)
library(ggplot2)
library(dplyr)
library(tidyr)
library(scales)
library(DT)

# ── ODE solver (deSolve-based, independent of mrgsolve for Shiny) ──────────
library(deSolve)

hcm_odes <- function(t, state, parms) {
  with(as.list(c(state, parms)), {

    # Mavacamten PK
    Cm <- A_c_m / V1_m  # nM
    E_mava <- Emax_mava * Cm^Hill_mava / (EC50_mava^Hill_mava + Cm^Hill_mava)
    E_mava <- min(E_mava, 0.80)

    dA_gut_m <- -ka_m * A_gut_m
    dA_c_m   <-  ka_m * A_gut_m - (CL_m_adj/V1_m)*A_c_m - (Q_m/V1_m)*A_c_m + (Q_m/V2_m)*A_p_m
    dA_p_m   <-  (Q_m/V1_m)*A_c_m - (Q_m/V2_m)*A_p_m

    # Beta-blocker PK
    Cbb     <- A_bb / V_bb
    E_bb_hr <- Emax_bb * Cbb / (EC50_bb + Cbb)
    E_bb_c  <- 0.60 * E_bb_hr

    dA_bb <- -(CL_bb/V_bb)*A_bb

    # Duty ratio
    DR_eff  <- DR_HCM * (1 - E_mava) * (1 - 0.5*E_bb_c)
    DR_eff  <- max(DR_eff, DR_normal*0.5)
    DR_rat  <- DR_eff / DR_normal

    # Calcium
    Ca_infl  <- 60 * DR_rat * (1 + 0.3*(Ca_SR/420 - 1))
    dCa_cyt  <- Ca_infl - 0.80*Ca_cyt - 0.20*Ca_cyt
    dCa_SR   <- 0.80*Ca_cyt*(200/420) - 0.05*(Ca_SR/420)*(Ca_cyt/200)*(420/200)

    # Calcineurin-NFAT
    CaN_t <- CaN_base + (1-CaN_base)*(Ca_cyt/200-1)/(1+abs(Ca_cyt/200-1))
    CaN_t <- max(0.05, min(1.0, CaN_t))
    dCalcn <- 0.50*(CaN_t - Calcn)

    NFAT_t <- NFAT_base + 0.80*Calcn
    NFAT_t <- min(NFAT_t, 1.0)
    dNFAT  <- 0.30*(NFAT_t - NFAT)

    # ERK
    LVOT_n  <- LVOT / LVOT_base
    TGFb_n  <- TGFb / TGFb_base
    ERK_t   <- ERK_base + 0.40*(LVOT_n-1)/(1+abs(LVOT_n-1)) + 0.30*(TGFb_n-1)/(1+abs(TGFb_n-1))
    ERK_t   <- max(0, min(1, ERK_t))
    dERK    <- 0.30*(ERK_t - ERK)

    # Hypertrophy
    hyp_sig <- NFAT*(1 + 0.4*ERK)
    dIVS    <- 2.5e-4*hyp_sig*(IVS_max - IVS) - 5e-5*IVS*(1-hyp_sig)
    dLVmass <- 0.001*(1 + 0.60*(IVS-11)/11 - LVmass)

    # TGFβ / Fibrosis
    TGFb_prod <- 0.015*(1 + 0.5*(LVOT_n-1) + 0.3*(DR_rat-1) + 0.2*(Col-1))
    dTGFb <- TGFb_prod - 0.020*TGFb
    dCol  <- 0.005*TGFb/TGFb_base - 0.002*Col

    # LVOT gradient
    LVOT_t <- max(LVOT_min, 2.5*(IVS-11) + 40*(DR_rat-1) + 0.25*(HR_st-70) + LVOT_min)
    dLVOT  <- 0.50*(LVOT_t - LVOT)

    # LVEDP
    C_eff   <- 1/(1+0.5*(Col-1))
    LVEDP_t <- max(5, 8 + 12/C_eff + 0.3*LVOT - 2*E_bb_hr)
    dLVEDP  <- 0.20*(LVEDP_t - LVEDP)

    # HR
    HR_t <- max(45, min(130, HR_base*(1-E_bb_hr) + 0.15*LVOT))
    dHR  <- 0.10*(HR_t - HR_st)

    # NT-proBNP
    NTpBNP_t <- 50*exp(0.14*LVEDP)*(1+0.20*(LVmass-1))
    dNTpBNP  <- 0.05*(NTpBNP_t - NTpBNP)

    # Troponin
    MVO2_idx  <- DR_rat*(1+0.20*LVOT/50)
    isch      <- max(0, (MVO2_idx-1.3)/1.3)
    Trop_t    <- TropI_base*(1+3*isch)
    dTropI    <- 0.02*(Trop_t - TropI) - 0.04*(TropI - TropI_base)

    # AF hazard
    LA_p <- LVEDP + 3 + 0.20*LVOT
    dAF  <- AF_haz_base*(LA_p/15)/(24*365)

    list(c(dA_gut_m, dA_c_m, dA_p_m, dA_bb,
           dCa_cyt, dCa_SR, dCalcn, dNFAT, dERK,
           dIVS, dLVmass, dTGFb, dCol,
           dLVOT, dLVEDP, dHR,
           dNTpBNP, dTropI, dAF),
         Cm_nM = Cm, E_mava_pct = 100*E_mava, DR_eff = DR_eff,
         Cbb_ng = Cbb, E_bb_pct = 100*E_bb_hr,
         LVEF = max(40, 72 - 2*(Col-1.3) - 1.5*(LVmass-1.6)),
         CO = max(2, (SV_base*(1-0.5*(LVOT-LVOT_base)/LVOT_base))*HR_st/1000),
         ECV = 20 + 8*(Col-1),
         NYHA = as.integer(max(1, min(4,
                               ifelse((30 - 0.15*(LVOT-45) - 0.5*(LVEDP-10)) > 20, 1,
                               ifelse((30 - 0.15*(LVOT-45) - 0.5*(LVEDP-10)) > 14, 2,
                               ifelse((30 - 0.15*(LVOT-45) - 0.5*(LVEDP-10)) > 10, 3, 4)))))),
         peakVO2 = max(8, 30 - 0.15*(LVOT-45) - 0.5*(LVEDP-10)))
  })
}

run_sim <- function(params, mava_mg, bb_mg, duration_days = 365,
                    cyp_factor = 1.0, srt = FALSE) {
  dt    <- 24  # hourly to daily steps
  times <- seq(0, duration_days*24, by = dt)
  MW_m  <- 471.0
  dose_mava_nmol <- mava_mg * 1e6 / MW_m * 0.93

  state0 <- c(
    A_gut_m = 0, A_c_m = 0, A_p_m = 0, A_bb = 0,
    Ca_cyt = 200, Ca_SR = 420,
    Calcn = 0.30, NFAT = 0.20, ERK = 0.15,
    IVS = if(srt) 11.0 else params$IVS_0,
    LVmass = params$LVmass_0,
    TGFb = params$TGFb_base,
    Col = params$Col_base,
    LVOT = if(srt) 10.0 else params$LVOT_base,
    LVEDP = params$LVEDP_base,
    HR_st = params$HR_base,
    NTpBNP = params$NTpBNP_0,
    TropI = params$TropI_base,
    AF = 0
  )

  # Build forcing functions (dose events)
  dose_events_mava <- NULL
  dose_events_bb   <- NULL
  if (mava_mg > 0 && !srt) {
    dose_times <- seq(0, duration_days*24, by = 24)
    dose_events_mava <- data.frame(
      var   = "A_gut_m",
      time  = dose_times,
      value = dose_mava_nmol,
      method= "add"
    )
  }
  if (bb_mg > 0) {
    dose_times_bb <- seq(0, duration_days*24, by = 24)
    dose_events_bb <- data.frame(
      var   = "A_bb",
      time  = dose_times_bb,
      value = bb_mg * 1000,
      method= "add"
    )
  }
  dose_events <- rbind(dose_events_mava, dose_events_bb)

  p_run <- c(params, CL_m_adj = params$CL_m * cyp_factor)

  out <- ode(y = state0, times = times, func = hcm_odes,
             parms = p_run, events = list(data = dose_events),
             method = "lsoda")
  as.data.frame(out) %>% mutate(day = time/24)
}

default_params <- list(
  ka_m = 0.693, CL_m = 1.80, V1_m = 42.0, Q_m = 0.90, V2_m = 125.0,
  EC50_mava = 85.0, Emax_mava = 0.65, Hill_mava = 1.2,
  ka_bb = 0.50, CL_bb = 18.0, V_bb = 180.0,
  EC50_bb = 25.0, Emax_bb = 0.35,
  DR_normal = 0.04, DR_HCM = 0.09,
  CaN_base = 0.30, NFAT_base = 0.20, ERK_base = 0.15,
  IVS_0 = 18.0, IVS_max = 32.0, LVmass_0 = 1.60, LVOT_base = 45.0, LVOT_min = 8.0,
  LVEDP_base = 20.0, HR_base = 80.0,
  TGFb_base = 1.0, Col_base = 1.30,
  NTpBNP_0 = 400.0, TropI_base = 25.0, AF_haz_base = 4.5,
  SV_base = 72.0
)

# ═══════════════════════════════════════════════════════════════
# UI
# ═══════════════════════════════════════════════════════════════
ui <- dashboardPage(
  skin = "red",

  dashboardHeader(
    title = span(icon("heartbeat"), "HCM QSP Dashboard"),
    titleWidth = 280
  ),

  dashboardSidebar(
    width = 280,
    sidebarMenu(
      id = "tabs",
      menuItem("Patient Profile",    tabName = "profile",   icon = icon("user-md")),
      menuItem("Drug PK",            tabName = "pk",        icon = icon("pills")),
      menuItem("Cardiac Mechanics",  tabName = "mechanics", icon = icon("heartbeat")),
      menuItem("Hypertrophy & Fibrosis", tabName = "remodel", icon = icon("chart-area")),
      menuItem("Clinical Endpoints", tabName = "clinical",  icon = icon("stethoscope")),
      menuItem("Scenario Comparison",tabName = "scenarios", icon = icon("exchange-alt")),
      menuItem("Risk Assessment",    tabName = "risk",      icon = icon("exclamation-triangle"))
    ),
    hr(),
    h5("  Global Parameters", style = "color:#ccc; padding-left:10px;"),
    sliderInput("sim_days", "Simulation Duration (days)",
                min=30, max=730, value=365, step=30),
    sliderInput("IVS_init", "Baseline IVS (mm)",
                min=13, max=30, value=18, step=1),
    sliderInput("LVOT_init", "Baseline LVOT Gradient (mmHg)",
                min=10, max=100, value=45, step=5),
    sliderInput("LVEDP_init", "Baseline LVEDP (mmHg)",
                min=10, max=35, value=20, step=1),
    selectInput("cyp_genotype", "CYP2C19 Genotype",
                choices = c("EM (extensive)" = 1.0,
                            "PM (poor)"      = 2.5,
                            "UM (ultra-rapid)"= 0.4),
                selected = 1.0)
  ),

  dashboardBody(
    tags$head(tags$style(HTML("
      .content-wrapper { background-color: #f4f6f9; }
      .box { border-top: 3px solid #c0392b; }
      .value-box .inner h3 { font-size:22px; }
    "))),

    tabItems(

      # ─────────────────────────────────────────────────────────
      # TAB 1: PATIENT PROFILE
      # ─────────────────────────────────────────────────────────
      tabItem(tabName = "profile",
        fluidRow(
          box(title = "Disease Overview — Hypertrophic Cardiomyopathy", width = 12,
              status = "danger", solidHeader = TRUE,
              p("Hypertrophic cardiomyopathy (HCM) is the most common inherited
                cardiac disease (~1:500 prevalence). It is caused by pathogenic
                variants in sarcomere genes (most commonly MYH7 and MYBPC3),
                leading to hypercontractility, asymmetric hypertrophy, diastolic
                dysfunction, LVOT obstruction, and increased risk of SCD and AF."),
              p(strong("Mavacamten (CAMZYOS®)"," — first-in-class cardiac myosin inhibitor
                (FDA-approved 2022) — stabilises the super-relaxed (SRX) state of myosin,
                reducing duty ratio, LVOT gradient, and myocardial energy demand.")))
        ),
        fluidRow(
          column(4,
            box(title = "Patient Parameters", width = NULL, status = "danger",
                solidHeader = TRUE,
                numericInput("age", "Age (years)", 45, 18, 80),
                selectInput("sex", "Sex", c("Male","Female")),
                selectInput("genotype", "Causal Variant",
                            c("MYH7" = "MYH7", "MYBPC3" = "MYBPC3",
                              "TNNT2" = "TNNT2", "Unknown" = "Unknown")),
                selectInput("obstr_type", "HCM Phenotype",
                            c("Obstructive (LVOTO)" = "obstr",
                              "Non-obstructive"     = "nonobstr",
                              "Labile obstructive"  = "labile")),
                sliderInput("LVEF_base", "Baseline LVEF (%)", 50, 85, 72, 1),
                sliderInput("HR_base_ui", "Baseline HR (bpm)", 50, 110, 80, 5)
            )
          ),
          column(4,
            box(title = "Current Medications", width = NULL, status = "warning",
                solidHeader = TRUE,
                numericInput("mava_dose", "Mavacamten Dose (mg QD)\n[0 = none]",
                             value = 5, min = 0, max = 15, step = 2.5),
                numericInput("bb_dose", "Beta-Blocker Dose (mg QD)\n[0 = none]",
                             value = 100, min = 0, max = 400, step = 25),
                checkboxInput("on_SRT", "Post-Septal Reduction Therapy (SRT)", FALSE),
                br(),
                actionButton("run_sim", "Run Simulation",
                             class = "btn-danger btn-lg", width = "100%",
                             icon = icon("play"))
            )
          ),
          column(4,
            box(title = "Baseline Status", width = NULL, status = "info",
                solidHeader = TRUE,
                valueBoxOutput("vb_IVS",   width = 12),
                valueBoxOutput("vb_LVOT",  width = 12),
                valueBoxOutput("vb_NTpBNP",width = 12),
                valueBoxOutput("vb_NYHA",  width = 12)
            )
          )
        ),
        fluidRow(
          box(title = "Disease Mechanistic Map (Overview)", width = 12,
              status = "danger", solidHeader = TRUE,
              img(src = "hcm_qsp_model.png",
                  style = "max-width:100%;border:1px solid #ddd;border-radius:4px;",
                  alt = "HCM QSP Mechanistic Map"))
        )
      ),

      # ─────────────────────────────────────────────────────────
      # TAB 2: DRUG PK
      # ─────────────────────────────────────────────────────────
      tabItem(tabName = "pk",
        fluidRow(
          box(title = "Mavacamten PK Parameters", width = 4,
              status = "primary", solidHeader = TRUE,
              sliderInput("pk_mava_dose", "Mavacamten Dose (mg QD)", 0, 15, 5, 2.5),
              sliderInput("pk_EC50", "EC50 (nM)", 40, 200, 85, 5),
              sliderInput("pk_Emax", "Emax (fraction)", 0.3, 0.9, 0.65, 0.05),
              selectInput("pk_cyp", "CYP2C19",
                          c("EM (standard)" = 1.0, "PM (x2.5 expo)" = 2.5, "UM (x0.4 expo)" = 0.4)),
              hr(),
              h5("Beta-Blocker"),
              sliderInput("pk_bb_dose", "BB Dose (mg QD)", 0, 400, 100, 25)
          ),
          column(8,
            fluidRow(
              box(title = "Mavacamten Plasma Concentration (nM)", width = 12,
                  status = "primary", solidHeader = TRUE,
                  plotOutput("pk_plot_mava", height = "280px"))
            ),
            fluidRow(
              box(title = "PD Effect: % Duty Ratio Suppression", width = 6,
                  status = "info", solidHeader = TRUE,
                  plotOutput("pk_plot_effect", height = "230px")),
              box(title = "Beta-Blocker Plasma Conc (ng/mL)", width = 6,
                  status = "info", solidHeader = TRUE,
                  plotOutput("pk_plot_bb", height = "230px"))
            )
          )
        ),
        fluidRow(
          box(title = "CYP2C19 Genotype Comparison — Mavacamten Exposure (30-day steady state)",
              width = 12, status = "primary", solidHeader = TRUE,
              plotOutput("pk_cyp_compare", height = "300px"))
        )
      ),

      # ─────────────────────────────────────────────────────────
      # TAB 3: CARDIAC MECHANICS
      # ─────────────────────────────────────────────────────────
      tabItem(tabName = "mechanics",
        fluidRow(
          valueBoxOutput("mb_LVOT",  width = 3),
          valueBoxOutput("mb_LVEDP", width = 3),
          valueBoxOutput("mb_LVEF",  width = 3),
          valueBoxOutput("mb_CO",    width = 3)
        ),
        fluidRow(
          box(title = "LVOT Gradient Over Time", width = 6,
              status = "danger", solidHeader = TRUE,
              plotOutput("mech_lvot", height = "280px")),
          box(title = "LVEDP (Diastolic Filling Pressure)", width = 6,
              status = "danger", solidHeader = TRUE,
              plotOutput("mech_lvedp", height = "280px"))
        ),
        fluidRow(
          box(title = "LVEF (%)", width = 4, status = "info", solidHeader = TRUE,
              plotOutput("mech_lvef", height = "250px")),
          box(title = "Cardiac Output (L/min)", width = 4, status = "info", solidHeader = TRUE,
              plotOutput("mech_co", height = "250px")),
          box(title = "Heart Rate (bpm)", width = 4, status = "info", solidHeader = TRUE,
              plotOutput("mech_hr", height = "250px"))
        )
      ),

      # ─────────────────────────────────────────────────────────
      # TAB 4: HYPERTROPHY & FIBROSIS
      # ─────────────────────────────────────────────────────────
      tabItem(tabName = "remodel",
        fluidRow(
          box(title = "IVS Thickness (mm) Over Time", width = 6,
              status = "warning", solidHeader = TRUE,
              plotOutput("remo_ivs", height = "280px")),
          box(title = "Relative Collagen Content (Fibrosis)", width = 6,
              status = "warning", solidHeader = TRUE,
              plotOutput("remo_col", height = "280px"))
        ),
        fluidRow(
          box(title = "TGF-β1 (Fibrosis Driver)", width = 4,
              status = "info", solidHeader = TRUE,
              plotOutput("remo_tgfb", height = "250px")),
          box(title = "NFAT Nuclear Fraction (Hypertrophy Signal)", width = 4,
              status = "info", solidHeader = TRUE,
              plotOutput("remo_nfat", height = "250px")),
          box(title = "Estimated CMR ECV (%)", width = 4,
              status = "info", solidHeader = TRUE,
              plotOutput("remo_ecv", height = "250px"))
        ),
        fluidRow(
          box(title = "Calcineurin–NFAT Hypertrophy Pathway",
              width = 12, status = "warning", solidHeader = TRUE,
              plotOutput("remo_signal_path", height = "260px"))
        )
      ),

      # ─────────────────────────────────────────────────────────
      # TAB 5: CLINICAL ENDPOINTS
      # ─────────────────────────────────────────────────────────
      tabItem(tabName = "clinical",
        fluidRow(
          valueBoxOutput("clin_NYHA",    width = 3),
          valueBoxOutput("clin_VO2",     width = 3),
          valueBoxOutput("clin_NTpBNP",  width = 3),
          valueBoxOutput("clin_TropI",   width = 3)
        ),
        fluidRow(
          box(title = "NT-proBNP (pg/mL) Over Time", width = 6,
              status = "success", solidHeader = TRUE,
              plotOutput("clin_ntpbnp_plot", height = "280px")),
          box(title = "Cardiac Troponin I (ng/L) Over Time", width = 6,
              status = "success", solidHeader = TRUE,
              plotOutput("clin_trop_plot", height = "280px"))
        ),
        fluidRow(
          box(title = "Estimated NYHA Class Over Time", width = 6,
              status = "success", solidHeader = TRUE,
              plotOutput("clin_nyha_plot", height = "250px")),
          box(title = "Estimated Peak VO₂ (mL/kg/min)", width = 6,
              status = "success", solidHeader = TRUE,
              plotOutput("clin_vo2_plot", height = "250px"))
        ),
        fluidRow(
          box(title = "Results Table (Monthly Snapshots)", width = 12,
              status = "success", solidHeader = TRUE,
              DTOutput("clin_table"))
        )
      ),

      # ─────────────────────────────────────────────────────────
      # TAB 6: SCENARIO COMPARISON
      # ─────────────────────────────────────────────────────────
      tabItem(tabName = "scenarios",
        fluidRow(
          box(title = "Scenario Definitions", width = 4,
              status = "primary", solidHeader = TRUE,
              checkboxGroupInput("scen_sel", "Select Scenarios to Compare",
                choices = c("Untreated HCM"              = "s1",
                            "Mavacamten 5mg (EM)"        = "s2",
                            "Mavacamten 10mg (EM)"       = "s3",
                            "Mavacamten 2.5mg (PM)"      = "s4",
                            "Beta-Blocker 200mg"         = "s5",
                            "Mava 5mg + BB 100mg"        = "s6",
                            "Post-SRT (Myectomy)"        = "s7"),
                selected = c("s1","s2","s3","s5","s6","s7")),
              hr(),
              selectInput("scen_var", "Primary Outcome",
                          choices = c("LVOT" = "LVOT",
                                      "LVEDP" = "LVEDP",
                                      "NTpBNP" = "NTpBNP",
                                      "IVS" = "IVS",
                                      "Col" = "Col",
                                      "TropI" = "TropI"),
                          selected = "LVOT"),
              actionButton("run_scen", "Run All Scenarios",
                           class = "btn-primary btn-block",
                           icon = icon("play"))
          ),
          column(8,
            fluidRow(
              box(title = "Primary Outcome Comparison Over Time", width = 12,
                  status = "primary", solidHeader = TRUE,
                  plotOutput("scen_main_plot", height = "320px"))
            )
          )
        ),
        fluidRow(
          box(title = "LVOT Gradient", width = 4, status = "info", solidHeader = TRUE,
              plotOutput("scen_lvot", height = "230px")),
          box(title = "NT-proBNP", width = 4, status = "info", solidHeader = TRUE,
              plotOutput("scen_bnp", height = "230px")),
          box(title = "IVS Thickness", width = 4, status = "info", solidHeader = TRUE,
              plotOutput("scen_ivs", height = "230px"))
        ),
        fluidRow(
          box(title = "52-Week Outcomes Summary Table", width = 12,
              status = "primary", solidHeader = TRUE,
              DTOutput("scen_table"))
        )
      ),

      # ─────────────────────────────────────────────────────────
      # TAB 7: RISK ASSESSMENT
      # ─────────────────────────────────────────────────────────
      tabItem(tabName = "risk",
        fluidRow(
          box(title = "HCM Risk-SCD Score Calculator", width = 4,
              status = "danger", solidHeader = TRUE,
              p("Based on HCM Risk-SCD model (O'Mahony et al., 2014) —
                ESC-endorsed 5-year SCD risk calculator"),
              numericInput("risk_age",    "Age (years)",      45, 16, 80),
              sliderInput("risk_IVS",     "Max IVS (mm)",     13, 35, 18, 1),
              sliderInput("risk_LA",      "LA Diameter (mm)", 28, 60, 42, 1),
              sliderInput("risk_LVOT",    "LVOT Gradient (mmHg)", 0, 100, 45, 5),
              checkboxInput("risk_FHx",   "Family Hx of SCD", FALSE),
              checkboxInput("risk_NSVT",  "NSVT on Holter",   FALSE),
              checkboxInput("risk_syncope","Unexplained Syncope", FALSE),
              actionButton("calc_risk", "Calculate SCD Risk",
                           class = "btn-danger", icon = icon("calculator"))
          ),
          column(8,
            fluidRow(
              valueBoxOutput("risk_scd5yr",  width = 4),
              valueBoxOutput("risk_af_prob", width = 4),
              valueBoxOutput("risk_icd_rec", width = 4)
            ),
            fluidRow(
              box(title = "SCD Risk Gauge", width = 6,
                  status = "danger", solidHeader = TRUE,
                  plotOutput("risk_gauge", height = "280px")),
              box(title = "Risk Factor Contribution", width = 6,
                  status = "danger", solidHeader = TRUE,
                  plotOutput("risk_factors", height = "280px"))
            )
          )
        ),
        fluidRow(
          box(title = "AF Cumulative Incidence Over Time", width = 6,
              status = "warning", solidHeader = TRUE,
              plotOutput("risk_af_plot", height = "280px")),
          box(title = "ICD Indication Threshold Tracking", width = 6,
              status = "warning", solidHeader = TRUE,
              plotOutput("risk_icd_plot", height = "280px"))
        ),
        fluidRow(
          box(title = "HCM Risk-SCD — Clinical Decision Support", width = 12,
              status = "danger", solidHeader = TRUE,
              tableOutput("risk_decision_table"))
        )
      )
    )
  )
)

# ═══════════════════════════════════════════════════════════════
# SERVER
# ═══════════════════════════════════════════════════════════════
server <- function(input, output, session) {

  pal7 <- c("#E53935","#1E88E5","#43A047","#F4511E",
            "#8E24AA","#039BE5","#795548")

  # ── Reactive: main simulation ──────────────────────────────
  sim_data <- eventReactive(input$run_sim, {
    params_run <- modifyList(default_params, list(
      IVS_0      = input$IVS_init,
      LVOT_base  = input$LVOT_init,
      LVEDP_base = input$LVEDP_init,
      HR_base    = input$HR_base_ui
    ))
    withProgress(message = "Running HCM simulation...", value = 0.5, {
      run_sim(params_run,
              mava_mg  = input$mava_dose,
              bb_mg    = input$bb_dose,
              duration_days = input$sim_days,
              cyp_factor    = as.numeric(input$cyp_genotype),
              srt           = input$on_SRT)
    })
  }, ignoreNULL = FALSE)

  last_row <- reactive({
    d <- sim_data(); d[nrow(d), ]
  })

  # ── Scenario simulations ────────────────────────────────────
  scen_data <- eventReactive(input$run_scen, {
    p <- default_params
    p$IVS_0 <- input$IVS_init; p$LVOT_base <- input$LVOT_init
    p$LVEDP_base <- input$LVEDP_init; p$HR_base <- input$HR_base_ui
    days <- input$sim_days

    defs <- list(
      s1 = list(mava=0,    bb=0,   srt=FALSE, cyp=1.0, name="Untreated HCM"),
      s2 = list(mava=5,    bb=0,   srt=FALSE, cyp=1.0, name="Mavacamten 5mg (EM)"),
      s3 = list(mava=10,   bb=0,   srt=FALSE, cyp=1.0, name="Mavacamten 10mg (EM)"),
      s4 = list(mava=2.5,  bb=0,   srt=FALSE, cyp=2.5, name="Mavacamten 2.5mg (PM)"),
      s5 = list(mava=0,    bb=200, srt=FALSE, cyp=1.0, name="Beta-Blocker 200mg"),
      s6 = list(mava=5,    bb=100, srt=FALSE, cyp=1.0, name="Mava 5mg + BB 100mg"),
      s7 = list(mava=0,    bb=0,   srt=TRUE,  cyp=1.0, name="Post-SRT (Myectomy)")
    )
    sel <- intersect(input$scen_sel, names(defs))
    withProgress(message = "Running scenarios...", value = 0, {
      result <- lapply(seq_along(sel), function(i) {
        sc <- defs[[sel[i]]]
        setProgress(i/length(sel))
        run_sim(p, sc$mava, sc$bb, days, sc$cyp, sc$srt) %>%
          mutate(Scenario = sc$name)
      })
      bind_rows(result)
    })
  }, ignoreNULL = FALSE)

  # ── Profile tab value boxes ──────────────────────────────────
  output$vb_IVS    <- renderValueBox(
    valueBox(paste0(input$IVS_init, " mm"), "IVS Thickness",
             icon = icon("ruler"), color = "red"))
  output$vb_LVOT   <- renderValueBox(
    valueBox(paste0(input$LVOT_init, " mmHg"), "LVOT Gradient",
             icon = icon("heartbeat"), color = "red"))
  output$vb_NTpBNP <- renderValueBox(
    valueBox("400 pg/mL", "NT-proBNP (baseline)",
             icon = icon("flask"), color = "orange"))
  output$vb_NYHA   <- renderValueBox(
    valueBox("II–III", "NYHA Class (baseline)",
             icon = icon("user"), color = "yellow"))

  # ── PK plots ─────────────────────────────────────────────────
  pk_sim <- reactive({
    p <- modifyList(default_params, list(EC50_mava = input$pk_EC50,
                                         Emax_mava = input$pk_Emax))
    run_sim(p, input$pk_mava_dose, input$pk_bb_dose, 90,
            as.numeric(input$pk_cyp), FALSE)
  })

  output$pk_plot_mava <- renderPlot({
    d <- pk_sim()
    ggplot(d, aes(day, Cm_nM)) +
      geom_line(color = "#1E88E5", linewidth = 1.2) +
      geom_hline(yintercept = input$pk_EC50, linetype="dashed", color="red") +
      labs(x="Day", y="Mavacamten (nM)",
           caption=paste0("Dashed: EC50 = ", input$pk_EC50, " nM")) +
      theme_bw(base_size = 13)
  })

  output$pk_plot_effect <- renderPlot({
    d <- pk_sim()
    ggplot(d, aes(day, E_mava_pct)) +
      geom_line(color = "#E53935", linewidth = 1.2) +
      geom_hline(yintercept = 50, linetype="dashed", color="grey50") +
      scale_y_continuous(limits = c(0,85)) +
      labs(x="Day", y="Duty Ratio Suppression (%)") +
      theme_bw(base_size = 13)
  })

  output$pk_plot_bb <- renderPlot({
    d <- pk_sim()
    ggplot(d, aes(day, Cbb_ng)) +
      geom_line(color = "#8E24AA", linewidth = 1.2) +
      labs(x="Day", y="Beta-Blocker (ng/mL)") +
      theme_bw(base_size = 13)
  })

  output$pk_cyp_compare <- renderPlot({
    p <- default_params
    cyps <- c("EM (1.0x)" = 1.0, "PM (2.5x)" = 2.5, "UM (0.4x)" = 0.4)
    dd <- lapply(names(cyps), function(nm) {
      run_sim(p, input$pk_mava_dose, 0, 30, cyps[nm], FALSE) %>%
        mutate(CYP = nm)
    }) %>% bind_rows()
    ggplot(dd, aes(day, Cm_nM, color = CYP)) +
      geom_line(linewidth = 1.2) +
      geom_hline(yintercept = 85, linetype="dashed", color="grey40") +
      scale_color_manual(values = c("#1E88E5","#E53935","#43A047")) +
      labs(x="Day", y="Mavacamten (nM)", color="CYP2C19",
           caption="Dashed: EC50 = 85 nM") +
      theme_bw(base_size = 13) + theme(legend.position = "bottom")
  })

  # ── Mechanics value boxes ────────────────────────────────────
  output$mb_LVOT  <- renderValueBox({
    d <- last_row()
    valueBox(sprintf("%.0f mmHg", d$LVOT), "Final LVOT Gradient",
             icon=icon("heartbeat"),
             color=if(d$LVOT<30)"green" else if(d$LVOT<50)"yellow" else "red")
  })
  output$mb_LVEDP <- renderValueBox({
    d <- last_row()
    valueBox(sprintf("%.0f mmHg", d$LVEDP), "Final LVEDP",
             icon=icon("tachometer-alt"),
             color=if(d$LVEDP<20)"green" else "yellow")
  })
  output$mb_LVEF <- renderValueBox({
    d <- last_row()
    valueBox(sprintf("%.0f %%", d$LVEF), "Final LVEF",
             icon=icon("heart"), color="blue")
  })
  output$mb_CO <- renderValueBox({
    d <- last_row()
    valueBox(sprintf("%.1f L/min", d$CO), "Cardiac Output",
             icon=icon("water"), color="purple")
  })

  # ── Mechanics plots ──────────────────────────────────────────
  output$mech_lvot  <- renderPlot({
    d <- sim_data()
    ggplot(d, aes(day, LVOT)) + geom_line(color="#E53935", size=1.1) +
      geom_hline(yintercept=30, linetype="dashed", color="grey50") +
      labs(x="Day", y="LVOT Gradient (mmHg)",
           caption="30 mmHg = obstruction threshold") + theme_bw(base_size=13)
  })
  output$mech_lvedp <- renderPlot({
    d <- sim_data()
    ggplot(d, aes(day, LVEDP)) + geom_line(color="#1E88E5", size=1.1) +
      geom_hline(yintercept=16, linetype="dashed", color="grey50") +
      labs(x="Day", y="LVEDP (mmHg)", caption="16 mmHg = elevated threshold") +
      theme_bw(base_size=13)
  })
  output$mech_lvef  <- renderPlot({
    d <- sim_data()
    ggplot(d, aes(day, LVEF)) + geom_line(color="#43A047", size=1.1) +
      scale_y_continuous(limits=c(40,90)) +
      labs(x="Day", y="LVEF (%)") + theme_bw(base_size=13)
  })
  output$mech_co    <- renderPlot({
    d <- sim_data()
    ggplot(d, aes(day, CO)) + geom_line(color="#8E24AA", size=1.1) +
      labs(x="Day", y="CO (L/min)") + theme_bw(base_size=13)
  })
  output$mech_hr    <- renderPlot({
    d <- sim_data()
    ggplot(d, aes(day, HR_st)) + geom_line(color="#F4511E", size=1.1) +
      labs(x="Day", y="Heart Rate (bpm)") + theme_bw(base_size=13)
  })

  # ── Remodeling plots ─────────────────────────────────────────
  output$remo_ivs  <- renderPlot({
    d <- sim_data()
    ggplot(d, aes(day, IVS)) + geom_line(color="#E53935", size=1.1) +
      geom_hline(yintercept=15, linetype="dashed") +
      labs(x="Day", y="IVS Thickness (mm)") + theme_bw(base_size=13)
  })
  output$remo_col  <- renderPlot({
    d <- sim_data()
    ggplot(d, aes(day, Col)) + geom_line(color="#F4511E", size=1.1) +
      geom_hline(yintercept=1.0, linetype="dashed") +
      labs(x="Day", y="Relative Collagen Content") + theme_bw(base_size=13)
  })
  output$remo_tgfb <- renderPlot({
    d <- sim_data()
    ggplot(d, aes(day, TGFb)) + geom_line(color="#795548", size=1.1) +
      labs(x="Day", y="TGF-β1 (ng/mL)") + theme_bw(base_size=13)
  })
  output$remo_nfat <- renderPlot({
    d <- sim_data()
    ggplot(d, aes(day, NFAT)) + geom_line(color="#AD1457", size=1.1) +
      scale_y_continuous(limits=c(0,1)) +
      labs(x="Day", y="Nuclear NFAT Fraction") + theme_bw(base_size=13)
  })
  output$remo_ecv  <- renderPlot({
    d <- sim_data()
    ggplot(d, aes(day, ECV)) + geom_line(color="#00838F", size=1.1) +
      geom_hline(yintercept=27, linetype="dashed") +
      labs(x="Day", y="Estimated CMR ECV (%)", caption="27% = normal upper limit") +
      theme_bw(base_size=13)
  })
  output$remo_signal_path <- renderPlot({
    d <- sim_data()
    d2 <- d %>% select(day, Calcn, NFAT, ERK) %>%
      pivot_longer(-day, names_to="Signal", values_to="Value")
    ggplot(d2, aes(day, Value, color=Signal)) +
      geom_line(size=1.1) +
      scale_color_manual(values=c(Calcn="#1E88E5", NFAT="#E53935", ERK="#43A047"),
                         labels=c("Calcineurin","NFAT (nuclear)","ERK1/2")) +
      scale_y_continuous(limits=c(0,1)) +
      labs(x="Day", y="Normalized Activity", color="Pathway") +
      theme_bw(base_size=13) + theme(legend.position="bottom")
  })

  # ── Clinical endpoint value boxes ────────────────────────────
  output$clin_NYHA   <- renderValueBox({
    d <- last_row()
    valueBox(paste("NYHA", d$NYHA), "Functional Class",
             icon=icon("user-md"),
             color=c("green","yellow","orange","red")[max(1,min(4,d$NYHA))])
  })
  output$clin_VO2    <- renderValueBox({
    d <- last_row()
    valueBox(sprintf("%.1f mL/kg/min", d$peakVO2), "Peak VO₂",
             icon=icon("running"),
             color=if(d$peakVO2>20)"green" else if(d$peakVO2>14)"yellow" else "red")
  })
  output$clin_NTpBNP <- renderValueBox({
    d <- last_row()
    valueBox(sprintf("%.0f pg/mL", d$NTpBNP), "NT-proBNP",
             icon=icon("flask"),
             color=if(d$NTpBNP<125)"green" else if(d$NTpBNP<400)"yellow" else "red")
  })
  output$clin_TropI  <- renderValueBox({
    d <- last_row()
    valueBox(sprintf("%.1f ng/L", d$TropI), "hs-Troponin I",
             icon=icon("vial"),
             color=if(d$TropI<30)"green" else "yellow")
  })

  output$clin_ntpbnp_plot <- renderPlot({
    d <- sim_data()
    ggplot(d, aes(day, NTpBNP)) + geom_line(color="#43A047",size=1.1) +
      geom_hline(yintercept=125,linetype="dashed",color="grey40") +
      geom_hline(yintercept=400,linetype="dotted",color="red") +
      labs(x="Day", y="NT-proBNP (pg/mL)",
           caption="125 pg/mL = normal; 400 pg/mL = elevated in HCM") +
      theme_bw(base_size=13)
  })
  output$clin_trop_plot   <- renderPlot({
    d <- sim_data()
    ggplot(d, aes(day, TropI)) + geom_line(color="#1E88E5",size=1.1) +
      geom_hline(yintercept=53,linetype="dashed",color="grey40") +
      labs(x="Day", y="Troponin I (ng/L)", caption="53 ng/L = 99th percentile URL") +
      theme_bw(base_size=13)
  })
  output$clin_nyha_plot   <- renderPlot({
    d <- sim_data()
    ggplot(d, aes(day, NYHA)) + geom_step(color="#E53935",size=1.1) +
      scale_y_continuous(breaks=1:4, labels=paste("NYHA",1:4)) +
      labs(x="Day", y="NYHA Class") + theme_bw(base_size=13)
  })
  output$clin_vo2_plot    <- renderPlot({
    d <- sim_data()
    ggplot(d, aes(day, peakVO2)) + geom_line(color="#8E24AA",size=1.1) +
      geom_hline(yintercept=20,linetype="dashed",color="grey40") +
      labs(x="Day", y="Peak VO₂ (mL/kg/min)", caption="20 mL/kg/min = NYHA II/III boundary") +
      theme_bw(base_size=13)
  })

  output$clin_table <- renderDT({
    d <- sim_data()
    d_monthly <- d %>% filter(day %% 30 == 0 | day == 0) %>%
      select(day, LVOT, LVEDP, NTpBNP, TropI, IVS, Col, NYHA, peakVO2, LVEF, ECV) %>%
      mutate(across(where(is.numeric), ~round(.,1)))
    datatable(d_monthly, options=list(pageLength=15, scrollX=TRUE),
              colnames=c("Day","LVOT(mmHg)","LVEDP(mmHg)","NT-proBNP(pg/mL)",
                         "TropI(ng/L)","IVS(mm)","Collagen","NYHA",
                         "Peak VO₂","LVEF(%)","ECV(%)"))
  })

  # ── Scenario comparison ──────────────────────────────────────
  output$scen_main_plot <- renderPlot({
    d <- scen_data()
    yvar <- input$scen_var
    ggplot(d, aes_string("day", yvar, color="Scenario")) +
      geom_line(size=1.0) +
      scale_color_manual(values=pal7) +
      labs(x="Day", y=yvar) +
      theme_bw(base_size=13) + theme(legend.position="bottom")
  })
  output$scen_lvot <- renderPlot({
    d <- scen_data()
    ggplot(d, aes(day, LVOT, color=Scenario)) +
      geom_line(size=0.9) + scale_color_manual(values=pal7) +
      geom_hline(yintercept=30, linetype="dashed") +
      labs(x="Day", y="LVOT (mmHg)") + theme_bw(base_size=11) +
      theme(legend.position="none")
  })
  output$scen_bnp  <- renderPlot({
    d <- scen_data()
    ggplot(d, aes(day, NTpBNP, color=Scenario)) +
      geom_line(size=0.9) + scale_color_manual(values=pal7) +
      labs(x="Day", y="NT-proBNP (pg/mL)") + theme_bw(base_size=11) +
      theme(legend.position="none")
  })
  output$scen_ivs  <- renderPlot({
    d <- scen_data()
    ggplot(d, aes(day, IVS, color=Scenario)) +
      geom_line(size=0.9) + scale_color_manual(values=pal7) +
      labs(x="Day", y="IVS (mm)") + theme_bw(base_size=11) +
      theme(legend.position="none")
  })

  output$scen_table <- renderDT({
    d <- scen_data()
    end_day <- max(d$day)
    tbl <- d %>%
      filter(abs(day - end_day) < 1) %>%
      select(Scenario, LVOT, LVEDP, NTpBNP, TropI, IVS, NYHA, peakVO2, LVEF) %>%
      mutate(across(where(is.numeric), ~round(.,1)))
    datatable(tbl, options=list(dom="t", pageLength=10))
  })

  # ── Risk Assessment ──────────────────────────────────────────
  risk_calc <- eventReactive(input$calc_risk, {
    # Simplified HCM Risk-SCD formula (continuous model approximation)
    score <- -0.998 * exp(
      0.15939858 * input$risk_IVS - 0.00294271 * input$risk_IVS^2
      + 0.0259082 * input$risk_LA
      + 0.00446131 * input$risk_LVOT
      + 0.4583082 * as.numeric(input$risk_FHx)
      + 0.82639195 * as.numeric(input$risk_NSVT)
      + 0.71650361 * as.numeric(input$risk_syncope)
      - 0.01799934 * input$risk_age
    )
    scd_5yr <- round((1 - 0.998^exp(score)) * 100, 2)
    list(scd_5yr = scd_5yr,
         af_annual = round(4.5 * (input$risk_LA / 42)^2, 1),
         icd = scd_5yr >= 6.0)
  })

  output$risk_scd5yr <- renderValueBox({
    r <- risk_calc()
    valueBox(paste0(r$scd_5yr, "%"), "5-Year SCD Risk",
             icon = icon("exclamation-triangle"),
             color = if(r$scd_5yr < 4)"green" else if(r$scd_5yr < 6)"yellow" else "red")
  })
  output$risk_af_prob <- renderValueBox({
    r <- risk_calc()
    valueBox(paste0(r$af_annual, "%/yr"), "Annual AF Incidence",
             icon = icon("heartbeat"), color = "orange")
  })
  output$risk_icd_rec <- renderValueBox({
    r <- risk_calc()
    valueBox(if(r$icd)"ICD Recommended" else "ICD Not Required",
             "ICD Decision",
             icon = icon("bolt"),
             color = if(r$icd)"red" else "green")
  })

  output$risk_gauge <- renderPlot({
    req(input$calc_risk)
    r <- risk_calc()
    df_g <- data.frame(
      category = c("Low\n(<4%)", "Intermediate\n(4-6%)", "High\n(≥6%)"),
      width = c(4, 2, 20),
      col = c("#43A047","#FDD835","#E53935")
    )
    ggplot(df_g, aes(x=0, y=width, fill=col)) +
      geom_bar(stat="identity") +
      geom_hline(yintercept=r$scd_5yr, color="black", size=2) +
      geom_text(aes(x=0, y=r$scd_5yr, label=paste0(r$scd_5yr,"%")),
                hjust=-0.2, size=6, fontface="bold") +
      scale_fill_identity() +
      coord_flip() +
      labs(title=paste("5-Year SCD Risk:", r$scd_5yr, "%"),
           y="Risk (%)", x="") +
      theme_bw(base_size=14)
  })

  output$risk_factors <- renderPlot({
    req(input$calc_risk)
    df_rf <- data.frame(
      Factor = c("IVS thickness","LA diameter","LVOT gradient",
                 "Family Hx SCD","NSVT on Holter","Unexplained syncope","Age"),
      Contribution = c(
        0.159 * input$risk_IVS - 0.00294 * input$risk_IVS^2,
        0.026 * input$risk_LA,
        0.0045 * input$risk_LVOT,
        0.458 * as.numeric(input$risk_FHx),
        0.826 * as.numeric(input$risk_NSVT),
        0.717 * as.numeric(input$risk_syncope),
        -0.018 * input$risk_age
      )
    )
    ggplot(df_rf, aes(reorder(Factor, Contribution), Contribution,
                      fill = Contribution > 0)) +
      geom_col() +
      scale_fill_manual(values = c("#1E88E5","#E53935"), guide = FALSE) +
      coord_flip() +
      labs(x = "", y = "Log-scale Contribution to SCD Risk",
           title = "Risk Factor Contributions") +
      theme_bw(base_size = 12)
  })

  output$risk_af_plot <- renderPlot({
    d <- sim_data()
    ggplot(d, aes(day, AF)) + geom_line(color="#F4511E", size=1.1) +
      labs(x="Day", y="Cumulative AF Hazard (a.u.)",
           title="AF Cumulative Incidence Over Time") + theme_bw(base_size=13)
  })

  output$risk_icd_plot <- renderPlot({
    d <- sim_data()
    ggplot(d, aes(day, IVS)) + geom_line(color="#E53935", size=1.1) +
      geom_hline(yintercept=30, linetype="dashed", color="red",
                 linewidth=1.2) +
      annotate("text", x=5, y=31, label="IVS ≥30mm → Massive hypertrophy",
               hjust=0, size=3.5, color="red") +
      labs(x="Day", y="Max IVS (mm)", title="IVS Tracking vs ICD Risk Threshold") +
      theme_bw(base_size=13)
  })

  output$risk_decision_table <- renderTable({
    data.frame(
      `5-yr SCD Risk` = c("<4%", "4–<6%", "≥6%"),
      `ICD Recommendation` = c("ICD not indicated", "ICD may be considered",
                                "ICD recommended (Class IIa)"),
      `Annual SCD Rate` = c("<0.8%/yr", "0.8–1.2%/yr", "≥1.2%/yr"),
      `ESC Guideline (2023)` = c("Conservative", "Shared decision", "Implant ICD"),
      check.names = FALSE
    )
  }, bordered = TRUE, striped = TRUE)
}

shinyApp(ui = ui, server = server)
