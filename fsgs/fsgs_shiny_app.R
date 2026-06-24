# =============================================================================
# FSGS QSP Shiny Dashboard
# Focal Segmental Glomerulosclerosis - Quantitative Systems Pharmacology
# deSolve-based ODE simulation, shinydashboard, plotly, DT
# =============================================================================

library(shiny)
library(shinydashboard)
library(plotly)
library(deSolve)
library(DT)

# =============================================================================
# ODE SYSTEM
# =============================================================================

fsgs_ode <- function(t, y, parms) {
  with(as.list(c(y, parms)), {

    t_day <- t / 24  # convert hours to days for disease ODEs

    # --- Drug concentrations (volume-normalized to ng/mL or ug/L equivalents) ---
    # Prednisolone: central (PRED1 in mg), convert to ng/mL using V1_pred=30L
    PRED_conc <- PRED1 / V1_pred * 1000  # ng/mL
    # Tacrolimus: central (TAC1 in ug), convert to ng/mL using V1_tac=85L
    TAC_conc  <- TAC1  / V1_tac  * 1000  # ng/mL
    # Rituximab: central (RTX1 in mg), convert to mg/L
    RTX_conc  <- RTX1  / V1_rtx         # mg/L
    # Sparsentan: central (SPARS1 in mg), convert to ng/mL
    SPARS_conc <- SPARS1 / V1_spars * 1000  # ng/mL

    # --- Emax drug effects ---
    E_pred  <- PRED_conc  / (EC50_pred  + PRED_conc  + 1e-9)
    E_tac   <- TAC_conc   / (EC50_tac   + TAC_conc   + 1e-9)
    E_rtx   <- RTX_conc   / (EC50_rtx   + RTX_conc   + 1e-9)
    E_spars <- SPARS_conc / (EC50_spars + SPARS_conc + 1e-9)

    # --- Disease state dynamics (per day, scaled to /hour for deSolve) ---
    k_scale <- 1/24  # /day -> /hour

    # Circulating permeability factor
    # Primary FSGS: B-cell dependent production; secondary: constitutive
    CLCF_prod <- (DIS_TYPE * BCELL * 0.8 + (1 - DIS_TYPE) * 0.2) * k_CLCF_prod
    CLCF_clear <- k_CLCF_cl * CLCF
    CLCF_rtx  <- 0.7 * E_rtx * CLCF  # RTX reduces CLCF via B-cell depletion
    dCLCF <- (CLCF_prod - CLCF_clear - CLCF_rtx) * k_scale

    # Podocyte fraction (0-1)
    CLCF_excess <- max(CLCF - 1.0, 0)
    POD_renewal  <- k_pod_renew * (1 - POD)
    POD_loss     <- (k_pod_CLCF * CLCF_excess + k_pod_COMP * (COMP - 1) + k_pod_RAAS * (RAAS - 1)) * POD
    dPOD <- (POD_renewal - POD_loss) * k_scale

    # Foot process effacement (0-1)
    FPE_form    <- k_FPE_form * CLCF_excess + k_FPE_COMP * max(COMP - 1, 0)
    FPE_repair  <- k_FPE_repair * (1 + 1.5 * E_pred + 3.0 * E_tac) * FPE
    dFPE <- (FPE_form - FPE_repair) * k_scale

    # Proteinuria (g/day)
    PROT_drive  <- k_prot_FPE * FPE + k_prot_POD * (1 - POD) + k_prot_RAAS * max(RAAS - 1, 0)
    PROT_target <- PROT_base + PROT_drive
    PROT_spars  <- 0.65 * E_spars
    dPROT <- k_prot_adj * (PROT_target * (1 - PROT_spars) - PROT) * k_scale

    # eGFR (mL/min/1.73m²)
    GFR_scar_loss  <- k_GFR_scar  * SCAR * GFR_c
    GFR_RAAS_loss  <- k_GFR_RAAS  * max(RAAS - 1, 0) * GFR_c
    GFR_spars_rec  <- k_GFR_spars * E_spars * max(GFR0 - GFR_c, 0)
    dGFR_c <- (-GFR_scar_loss - GFR_RAAS_loss + GFR_spars_rec) * k_scale

    # Glomerular sclerosis (0-1, irreversible)
    SCAR_acc <- k_SCAR_POD * max(1 - POD - 0.2, 0) + k_SCAR_TGFb * max(TGFb - 1, 0)
    dSCAR <- SCAR_acc * k_scale

    # TGF-beta (normalized, baseline=1)
    TGFb_base_prod <- k_TGFb_base
    TGFb_RAAS_drive <- k_TGFb_RAAS * max(RAAS - 1, 0)
    TGFb_POD_drive  <- k_TGFb_POD  * max(1 - POD, 0)
    TGFb_pred_supp  <- 0.50 * E_pred * TGFb
    TGFb_clear      <- k_TGFb_cl * (TGFb - 1)
    dTGFb <- (TGFb_base_prod + TGFb_RAAS_drive + TGFb_POD_drive - TGFb_pred_supp - TGFb_clear) * k_scale

    # Complement activity (normalized, baseline=1)
    COMP_prod  <- k_COMP_prod * INFLAM
    COMP_pred  <- 0.30 * E_pred * COMP
    COMP_clear <- k_COMP_cl * (COMP - 1)
    dCOMP <- (COMP_prod - COMP_pred - COMP_clear) * k_scale

    # Inflammation index (normalized, baseline=1)
    INFLAM_COMP  <- k_INFLAM_COMP * max(COMP - 1, 0)
    INFLAM_pred  <- 0.60 * E_pred * INFLAM
    INFLAM_clear <- k_INFLAM_cl * (INFLAM - 1)
    dINFLAM <- (INFLAM_COMP - INFLAM_pred - INFLAM_clear) * k_scale

    # RAAS activity (normalized, baseline=1)
    RAAS_drive  <- k_RAAS_drive * max(1 - GFR_c / GFR0, 0)
    RAAS_spars  <- 0.70 * E_spars * (RAAS - 1)
    RAAS_clear  <- k_RAAS_cl * (RAAS - 1)
    dRAAS <- (RAAS_drive - RAAS_spars - RAAS_clear) * k_scale

    # B-cell fraction (0-1)
    BCELL_renew <- k_BCELL_renew * (1 - BCELL)
    BCELL_rtx   <- k_BCELL_rtx * E_rtx * BCELL
    dBCELL <- (BCELL_renew - BCELL_rtx) * k_scale

    # ==========================================================================
    # Drug PK (2-CMT, per hour, oral absorption for pred/tac/spars, IV for rtx)
    # ==========================================================================

    # --- Prednisolone ---
    dPRED1 <- -( CL_pred/V1_pred + Q_pred/V1_pred ) * PRED1 + (Q_pred/V2_pred) * PRED2
    dPRED2 <-  ( Q_pred/V1_pred ) * PRED1 - (Q_pred/V2_pred) * PRED2

    # --- Tacrolimus ---
    dTAC1 <- -( CL_tac/V1_tac + Q_tac/V1_tac ) * TAC1 + (Q_tac/V2_tac) * TAC2
    dTAC2 <-  ( Q_tac/V1_tac ) * TAC1 - (Q_tac/V2_tac) * TAC2

    # --- Rituximab ---
    dRTX1 <- -( CL_rtx/V1_rtx + Q_rtx/V1_rtx ) * RTX1 + (Q_rtx/V2_rtx) * RTX2
    dRTX2 <-  ( Q_rtx/V1_rtx ) * RTX1 - (Q_rtx/V2_rtx) * RTX2

    # --- Sparsentan ---
    dSPARS1 <- -( CL_spars/V1_spars + Q_spars/V1_spars ) * SPARS1 + (Q_spars/V2_spars) * SPARS2
    dSPARS2 <-  ( Q_spars/V1_spars ) * SPARS1 - (Q_spars/V2_spars) * SPARS2

    list(c(dCLCF, dPOD, dFPE, dPROT, dGFR_c, dSCAR, dTGFb, dCOMP, dINFLAM, dRAAS, dBCELL,
           dPRED1, dPRED2, dTAC1, dTAC2, dRTX1, dRTX2, dSPARS1, dSPARS2))
  })
}

# =============================================================================
# DEFAULT PARAMETERS
# =============================================================================

default_parms <- function(age = 40, wt = 70, base_gfr = 90,
                          base_prot = 5, dis_type = 1) {

  bsa <- 0.007184 * wt^0.425 * (age * 0.5 + 165)^0.725 / 10000  # rough BSA estimate

  list(
    # Disease type: 1 = primary (B-cell dependent), 0 = secondary
    DIS_TYPE     = dis_type,
    GFR0         = base_gfr,
    PROT_base    = 0.15,

    # Disease kinetics
    k_CLCF_prod  = 0.12,   k_CLCF_cl  = 0.12,
    k_pod_renew  = 0.002,  k_pod_CLCF = 0.05, k_pod_COMP = 0.02, k_pod_RAAS = 0.015,
    k_FPE_form   = 0.08,   k_FPE_COMP = 0.03, k_FPE_repair = 0.06,
    k_prot_FPE   = 8.0,    k_prot_POD = 3.0,  k_prot_RAAS  = 1.5,  k_prot_adj = 0.3,
    k_GFR_scar   = 0.002,  k_GFR_RAAS = 0.01, k_GFR_spars  = 0.005,
    k_SCAR_POD   = 0.003,  k_SCAR_TGFb = 0.002,
    k_TGFb_base  = 0.05,   k_TGFb_RAAS = 0.08, k_TGFb_POD = 0.10, k_TGFb_cl = 0.15,
    k_COMP_prod  = 0.10,   k_COMP_cl   = 0.15,
    k_INFLAM_COMP = 0.08,  k_INFLAM_cl = 0.12,
    k_RAAS_drive  = 0.10,  k_RAAS_cl   = 0.15,
    k_BCELL_renew = 0.004, k_BCELL_rtx = 0.5,

    # Prednisolone PK (L, L/h)
    CL_pred  = 12,   V1_pred  = 30,   V2_pred  = 45,   Q_pred  = 8,    ka_pred = 2.0,
    EC50_pred = 150,

    # Tacrolimus PK
    CL_tac   = 2.5,  V1_tac   = 85,   V2_tac   = 150,  Q_tac   = 15,   ka_tac  = 0.4,
    EC50_tac  = 8,

    # Rituximab PK
    CL_rtx   = 0.014, V1_rtx  = 3.1,  V2_rtx  = 1.7,  Q_rtx   = 0.012,
    EC50_rtx  = 0.05,

    # Sparsentan PK
    CL_spars = 8.5,  V1_spars = 60,   V2_spars = 80,   Q_spars = 5,    ka_spars = 1.2,
    EC50_spars = 200,

    # BSA for RTX dosing
    BSA = max(bsa, 1.5)
  )
}

# =============================================================================
# INITIAL CONDITIONS
# =============================================================================

make_inits <- function(base_gfr = 90, base_prot = 5, dis_type = 1) {
  CLCF_init <- ifelse(dis_type == 1, 2.5, 1.4)
  c(
    CLCF=CLCF_init, POD=0.70, FPE=0.55, PROT=base_prot, GFR_c=base_gfr,
    SCAR=0.05, TGFb=1.3, COMP=1.2, INFLAM=1.2, RAAS=1.3, BCELL=1.0,
    PRED1=0, PRED2=0, TAC1=0, TAC2=0, RTX1=0, RTX2=0, SPARS1=0, SPARS2=0
  )
}

# =============================================================================
# DOSING EVENTS
# =============================================================================

make_events <- function(scenario, parms, dose_pred = 60, dose_tac_trough = 6,
                        dose_spars = 800, wt = 70) {

  bsa  <- parms$BSA
  events <- NULL

  if (scenario %in% c("PRED", "PRED_TAC", "PRED_TAC_RTX", "COMBO")) {
    # Prednisolone: 60mg/day wk1-8, 30mg/day wk9-16, 20mg/day wk17-24, 10mg/day thereafter
    pred_dose <- dose_pred  # mg
    pred_times_phase1 <- seq(0,    56*24 - 24, by=24)  # 8 weeks
    pred_times_phase2 <- seq(56*24, 112*24 - 24, by=24)  # wk 9-16
    pred_times_phase3 <- seq(112*24, 168*24 - 24, by=24)  # wk17-24
    pred_times_phase4 <- seq(168*24, 363*24, by=24)

    if (scenario == "COMBO") pred_dose <- 30

    d1 <- data.frame(var="PRED1", time=pred_times_phase1, value=pred_dose * 0.82, method="add")
    d2 <- data.frame(var="PRED1", time=pred_times_phase2, value=pred_dose * 0.5  * 0.82, method="add")
    d3 <- data.frame(var="PRED1", time=pred_times_phase3, value=pred_dose * 0.33 * 0.82, method="add")
    d4 <- data.frame(var="PRED1", time=pred_times_phase4, value=pred_dose * 0.17 * 0.82, method="add")
    events <- rbind(events, d1, d2, d3, d4)
  }

  if (scenario %in% c("PRED_TAC", "PRED_TAC_RTX", "COMBO")) {
    # Tacrolimus: 0.05 mg/kg/day -> target trough 5-8 ng/mL
    tac_dose_mg <- wt * 0.05  # mg/day
    tac_times   <- seq(0, 363*24, by=12)  # BID
    dt <- data.frame(var="TAC1", time=tac_times, value=(tac_dose_mg/2) * 0.25, method="add")
    events <- rbind(events, dt)
  }

  if (scenario == "PRED_TAC_RTX") {
    # Rituximab 375 mg/m² IV weekly x4
    rtx_dose <- 375 * bsa
    rtx_times <- c(0, 7, 14, 21) * 24
    dr <- data.frame(var="RTX1", time=rtx_times, value=rtx_dose, method="add")
    events <- rbind(events, dr)
  }

  if (scenario %in% c("SPARS", "COMBO")) {
    # Sparsentan 800 mg QD
    spars_times <- seq(0, 363*24, by=24)
    ds <- data.frame(var="SPARS1", time=spars_times, value=dose_spars * 0.85, method="add")
    events <- rbind(events, ds)
  }

  if (!is.null(events)) {
    events <- events[order(events$time), ]
    rownames(events) <- NULL
  }
  events
}

# =============================================================================
# SIMULATION RUNNER
# =============================================================================

run_sim <- function(scenario, parms, inits, dose_pred=60, dose_tac_trough=6,
                    dose_spars=800, wt=70) {

  times <- seq(0, 364*24, by=6)  # every 6 hours for 52 weeks
  ev    <- make_events(scenario, parms, dose_pred, dose_tac_trough, dose_spars, wt)

  if (!is.null(ev) && nrow(ev) > 0) {
    out <- tryCatch(
      ode(y=inits, times=times, func=fsgs_ode, parms=parms,
          method="lsoda", events=list(data=ev)),
      error = function(e) NULL
    )
  } else {
    out <- tryCatch(
      ode(y=inits, times=times, func=fsgs_ode, parms=parms, method="lsoda"),
      error = function(e) NULL
    )
  }

  if (is.null(out)) return(NULL)
  df <- as.data.frame(out)
  df$time_day  <- df$time / 24
  df$time_week <- df$time / 168
  df$scenario  <- scenario
  # Clamp physically meaningful states
  df$POD   <- pmax(0, pmin(1, df$POD))
  df$FPE   <- pmax(0, pmin(1, df$FPE))
  df$SCAR  <- pmax(0, pmin(1, df$SCAR))
  df$BCELL <- pmax(0, pmin(1, df$BCELL))
  df$PROT  <- pmax(0, df$PROT)
  df$GFR_c <- pmax(0, df$GFR_c)
  df$CLCF  <- pmax(0, df$CLCF)
  # Simulated biomarkers
  df$suPAR    <- 2000 + df$CLCF * 1500 + (1 - df$POD) * 800
  df$uPodocalyxin <- (1 - df$POD) * 500 + df$FPE * 300   # pg/mL
  df$NephrinShed  <- df$FPE * 250 + (1 - df$POD) * 180    # ng/mL
  # Drug concentrations
  df$PRED_conc  <- df$PRED1 / 30  * 1000   # ng/mL
  df$TAC_conc   <- df$TAC1  / 85  * 1000   # ng/mL
  df$RTX_conc   <- df$RTX1  / 3.1         # mg/L
  df$SPARS_conc <- df$SPARS1/ 60  * 1000   # ng/mL
  df
}

# =============================================================================
# SCENARIO CONFIG
# =============================================================================
scenario_labels <- c(
  NatHist       = "1: Natural History",
  PRED          = "2: Prednisolone",
  PRED_TAC      = "3: Pred + Tacrolimus",
  PRED_TAC_RTX  = "4: Pred + Tac + Rituximab",
  SPARS         = "5: Sparsentan",
  COMBO         = "6: Pred + Tac + Sparsentan"
)
scenario_colors <- c(
  NatHist       = "#E63946",
  PRED          = "#F4A261",
  PRED_TAC      = "#2A9D8F",
  PRED_TAC_RTX  = "#264653",
  SPARS         = "#8338EC",
  COMBO         = "#3A86FF"
)

# =============================================================================
# UI
# =============================================================================

ui <- dashboardPage(
  skin = "blue",
  dashboardHeader(title = "FSGS QSP Dashboard"),

  dashboardSidebar(
    sidebarMenu(
      menuItem("Patient Profile",         tabName = "tab1", icon = icon("user")),
      menuItem("Proteinuria & Antibodies",tabName = "tab2", icon = icon("flask")),
      menuItem("Glomerular Pathophysiology", tabName="tab3", icon=icon("microscope")),
      menuItem("Clinical Endpoints",      tabName = "tab4", icon = icon("heartbeat")),
      menuItem("Scenario Comparison",     tabName = "tab5", icon = icon("chart-line")),
      menuItem("Drug PK & Biomarkers",    tabName = "tab6", icon = icon("pills"))
    )
  ),

  dashboardBody(
    tags$head(tags$style(HTML("
      .content-wrapper { background-color: #f4f6f9; }
      .box { border-radius: 6px; }
      .box-header { font-weight: 600; }
      .nav-tabs-custom>.nav-tabs>li.active { border-top-color:#3c8dbc; }
      .small-box { border-radius: 6px; }
    "))),

    tabItems(

      # -----------------------------------------------------------------------
      # TAB 1: Patient Profile
      # -----------------------------------------------------------------------
      tabItem(tabName="tab1",
        fluidRow(
          box(title="Patient Parameters", status="primary", solidHeader=TRUE, width=4,
            sliderInput("age",   "Age (years)",            min=10, max=80, value=40, step=1),
            sliderInput("wt",    "Weight (kg)",            min=30, max=120, value=70, step=1),
            sliderInput("base_gfr",  "Baseline eGFR (mL/min/1.73m²)", min=20, max=130, value=85, step=1),
            sliderInput("base_prot", "Baseline Proteinuria (g/day)",   min=0.5, max=20, value=6, step=0.1),
            radioButtons("dis_type", "Disease Type",
              choices=c("Primary FSGS"=1, "Secondary FSGS"=0, "Genetic FSGS"=0),
              selected=1)
          ),
          box(title="Treatment Selection", status="info", solidHeader=TRUE, width=4,
            radioButtons("scenario", "Treatment Scenario",
              choices=setNames(names(scenario_labels), scenario_labels),
              selected="PRED_TAC"),
            hr(),
            numericInput("dose_pred",  "Prednisolone initial dose (mg/day)", value=60, min=5, max=120),
            numericInput("dose_tac",   "Tacrolimus target trough (ng/mL)",   value=6,  min=2, max=15),
            numericInput("dose_spars", "Sparsentan dose (mg/day)",           value=800, min=100, max=1200)
          ),
          box(title="Patient Summary", status="success", solidHeader=TRUE, width=4,
            DTOutput("patientSummary")
          )
        ),
        fluidRow(
          valueBoxOutput("vb_gfr",  width=3),
          valueBoxOutput("vb_prot", width=3),
          valueBoxOutput("vb_stage",width=3),
          valueBoxOutput("vb_risk", width=3)
        ),
        fluidRow(
          box(title="Baseline Disease State Overview", status="primary", solidHeader=TRUE, width=12,
            plotlyOutput("baselineRadar", height="300px")
          )
        )
      ),

      # -----------------------------------------------------------------------
      # TAB 2: Proteinuria & Antibody Dynamics
      # -----------------------------------------------------------------------
      tabItem(tabName="tab2",
        fluidRow(
          box(title="Proteinuria Over Time", status="primary", solidHeader=TRUE, width=8,
            plotlyOutput("protPlot", height="320px")
          ),
          box(title="Remission Summary", status="info", solidHeader=TRUE, width=4,
            DTOutput("remissionTable")
          )
        ),
        fluidRow(
          box(title="Circulating Permeability Factor (CLCF)", status="warning", solidHeader=TRUE, width=6,
            plotlyOutput("clcfPlot", height="280px")
          ),
          box(title="B-Cell Dynamics", status="danger", solidHeader=TRUE, width=6,
            plotlyOutput("bcellPlot", height="280px")
          )
        )
      ),

      # -----------------------------------------------------------------------
      # TAB 3: Glomerular Pathophysiology
      # -----------------------------------------------------------------------
      tabItem(tabName="tab3",
        fluidRow(
          box(title="Podocyte Fraction", status="danger", solidHeader=TRUE, width=6,
            plotlyOutput("podPlot", height="280px")
          ),
          box(title="Foot Process Effacement (%)", status="warning", solidHeader=TRUE, width=6,
            plotlyOutput("fpePlot", height="280px")
          )
        ),
        fluidRow(
          box(title="Complement Activity", status="info", solidHeader=TRUE, width=4,
            plotlyOutput("compPlot", height="260px")
          ),
          box(title="Inflammation Index", status="warning", solidHeader=TRUE, width=4,
            plotlyOutput("inflamPlot", height="260px")
          ),
          box(title="TGF-β & Glomerular Sclerosis", status="danger", solidHeader=TRUE, width=4,
            plotlyOutput("tgfScarPlot", height="260px")
          )
        )
      ),

      # -----------------------------------------------------------------------
      # TAB 4: Clinical Endpoints
      # -----------------------------------------------------------------------
      tabItem(tabName="tab4",
        fluidRow(
          box(title="eGFR Trajectory (52 weeks)", status="primary", solidHeader=TRUE, width=7,
            plotlyOutput("gfrPlot", height="320px")
          ),
          box(title="Clinical Outcome Table", status="info", solidHeader=TRUE, width=5,
            DTOutput("outcomeTable")
          )
        ),
        fluidRow(
          box(title="Remission Probability Over Time", status="success", solidHeader=TRUE, width=6,
            plotlyOutput("remProbPlot", height="280px")
          ),
          box(title="ESRD Risk (GFR Slope Estimate)", status="danger", solidHeader=TRUE, width=6,
            plotlyOutput("esrdPlot", height="280px")
          )
        )
      ),

      # -----------------------------------------------------------------------
      # TAB 5: Scenario Comparison
      # -----------------------------------------------------------------------
      tabItem(tabName="tab5",
        fluidRow(
          box(title="All Scenarios: 2×2 Panel (PROT · GFR · POD · SCAR)",
              status="primary", solidHeader=TRUE, width=12,
            plotlyOutput("comparisonPanel", height="520px")
          )
        ),
        fluidRow(
          box(title="Complete Remission Rate at 6 Months", status="success", solidHeader=TRUE, width=6,
            plotlyOutput("crBarPlot", height="300px")
          ),
          box(title="eGFR Change from Baseline at 52 Weeks", status="info", solidHeader=TRUE, width=6,
            plotlyOutput("gfrChangeBar", height="300px")
          )
        )
      ),

      # -----------------------------------------------------------------------
      # TAB 6: Drug PK & Biomarkers
      # -----------------------------------------------------------------------
      tabItem(tabName="tab6",
        fluidRow(
          box(title="Prednisolone PK (ng/mL)", status="primary", solidHeader=TRUE, width=6,
            plotlyOutput("predPK", height="260px")
          ),
          box(title="Tacrolimus PK (ng/mL) with Therapeutic Window", status="info", solidHeader=TRUE, width=6,
            plotlyOutput("tacPK", height="260px")
          )
        ),
        fluidRow(
          box(title="Rituximab PK (mg/L)", status="warning", solidHeader=TRUE, width=6,
            plotlyOutput("rtxPK", height="260px")
          ),
          box(title="Sparsentan PK (ng/mL)", status="success", solidHeader=TRUE, width=6,
            plotlyOutput("sparsPK", height="260px")
          )
        ),
        fluidRow(
          box(title="Biomarker Panel", status="danger", solidHeader=TRUE, width=8,
            plotlyOutput("biomarkerPlot", height="300px")
          ),
          box(title="Tacrolimus Target Attainment", status="info", solidHeader=TRUE, width=4,
            gaugeOutput("tacGauge"),
            verbatimTextOutput("tacTA")
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

  # --- Reactive: parameters ---
  params_r <- reactive({
    default_parms(age=input$age, wt=input$wt, base_gfr=input$base_gfr,
                  base_prot=input$base_prot, dis_type=as.numeric(input$dis_type))
  })

  inits_r <- reactive({
    make_inits(input$base_gfr, input$base_prot, as.numeric(input$dis_type))
  })

  # --- Reactive: single scenario simulation ---
  sim_r <- reactive({
    parms <- params_r(); inits <- inits_r()
    withProgress(message="Running simulation...", value=0.5, {
      run_sim(input$scenario, parms, inits,
              dose_pred=input$dose_pred, dose_tac_trough=input$dose_tac,
              dose_spars=input$dose_spars, wt=input$wt)
    })
  })

  # --- Reactive: all 6 scenarios ---
  all_scen_r <- reactive({
    parms <- params_r(); inits <- inits_r()
    withProgress(message="Running all scenarios...", value=0, {
      out <- lapply(names(scenario_labels), function(sc) {
        incProgress(1/6, detail=sc)
        run_sim(sc, parms, inits, dose_pred=input$dose_pred,
                dose_tac_trough=input$dose_tac, dose_spars=input$dose_spars,
                wt=input$wt)
      })
    })
    names(out) <- names(scenario_labels)
    out
  })

  # Helper: weekly snapshot
  weekly_snap <- function(df, col, wk) {
    idx <- which.min(abs(df$time_week - wk))
    df[[col]][idx]
  }

  # ==========================================================================
  # TAB 1 OUTPUTS
  # ==========================================================================
  output$patientSummary <- renderDT({
    gfr <- input$base_gfr; prot <- input$base_prot
    ckd <- ifelse(gfr>=60, "G2", ifelse(gfr>=45, "G3a", ifelse(gfr>=30, "G3b", ifelse(gfr>=15,"G4","G5"))))
    neph <- ifelse(prot>=3.5,"Nephrotic","Sub-nephrotic")
    df <- data.frame(
      Parameter = c("Age","Weight","Baseline eGFR","CKD Stage","Proteinuria","Category","Disease Type"),
      Value     = c(paste0(input$age," yr"), paste0(input$wt," kg"),
                    paste0(gfr," mL/min/1.73m²"), ckd,
                    paste0(prot," g/day"), neph,
                    ifelse(input$dis_type==1,"Primary","Secondary/Genetic"))
    )
    datatable(df, rownames=FALSE, options=list(dom="t", pageLength=10))
  })

  output$vb_gfr <- renderValueBox({
    valueBox(paste0(input$base_gfr," mL/min"), "Baseline eGFR", icon=icon("tint"), color="blue")
  })
  output$vb_prot <- renderValueBox({
    valueBox(paste0(input$base_prot," g/day"), "Proteinuria", icon=icon("flask"), color="orange")
  })
  output$vb_stage <- renderValueBox({
    gfr <- input$base_gfr
    ckd <- ifelse(gfr>=60,"G2",ifelse(gfr>=45,"G3a",ifelse(gfr>=30,"G3b",ifelse(gfr>=15,"G4","G5"))))
    valueBox(ckd, "CKD Stage", icon=icon("heartbeat"), color="red")
  })
  output$vb_risk <- renderValueBox({
    risk <- ifelse(input$base_prot>8,"Very High",ifelse(input$base_prot>3.5,"High","Moderate"))
    valueBox(risk, "ESRD Risk", icon=icon("exclamation-triangle"), color="yellow")
  })

  output$baselineRadar <- renderPlotly({
    prot_n <- min(input$base_prot/15, 1)
    gfr_n  <- 1 - min(input$base_gfr/130, 1)
    theta <- c("Proteinuria","GFR loss","CLCF","Inflammation","COMP","TGF-β","Proteinuria")
    vals  <- c(prot_n, gfr_n, 0.65, 0.45, 0.42, 0.50, prot_n)
    plot_ly(type="scatterpolar", r=vals, theta=theta, fill="toself",
            line=list(color="#3c8dbc"), fillcolor="rgba(60,141,188,0.3)") %>%
      layout(polar=list(radialaxis=list(visible=TRUE,range=c(0,1))),
             paper_bgcolor="#ffffff", plot_bgcolor="#ffffff",
             margin=list(l=40,r=40,t=20,b=20))
  })

  # ==========================================================================
  # TAB 2 OUTPUTS
  # ==========================================================================
  output$protPlot <- renderPlotly({
    df <- sim_r(); req(!is.null(df))
    init_prot <- df$PROT[1]
    cr_thresh <- 0.3; pr_thresh <- init_prot * 0.5
    plot_ly(df, x=~time_week, y=~PROT, type="scatter", mode="lines",
            line=list(color=scenario_colors[input$scenario], width=2.5),
            name="Proteinuria") %>%
      add_lines(x=~time_week, y=rep(cr_thresh, nrow(df)), line=list(color="green",dash="dash",width=1.5), name="CR (<0.3 g/day)") %>%
      add_lines(x=~time_week, y=rep(pr_thresh, nrow(df)), line=list(color="orange",dash="dot",width=1.5), name="PR threshold") %>%
      layout(xaxis=list(title="Time (weeks)"), yaxis=list(title="Proteinuria (g/day)"),
             paper_bgcolor="#ffffff", plot_bgcolor="#fafafa",
             legend=list(orientation="h",y=-0.2))
  })

  output$clcfPlot <- renderPlotly({
    df <- sim_r(); req(!is.null(df))
    plot_ly(df, x=~time_week, y=~CLCF, type="scatter", mode="lines",
            line=list(color="#E63946", width=2), name="CLCF") %>%
      add_lines(x=~time_week, y=rep(1, nrow(df)), line=list(color="gray",dash="dash"), name="Normal") %>%
      layout(xaxis=list(title="Week"), yaxis=list(title="CLCF (normalized)"),
             paper_bgcolor="#ffffff", plot_bgcolor="#fafafa")
  })

  output$bcellPlot <- renderPlotly({
    df <- sim_r(); req(!is.null(df))
    plot_ly(df, x=~time_week, y=~BCELL*100, type="scatter", mode="lines",
            line=list(color="#264653", width=2), name="B cells") %>%
      layout(xaxis=list(title="Week"), yaxis=list(title="B-Cell Fraction (%)"),
             paper_bgcolor="#ffffff", plot_bgcolor="#fafafa")
  })

  output$remissionTable <- renderDT({
    df <- sim_r(); req(!is.null(df))
    init_prot <- df$PROT[1]
    get_remission <- function(wk) {
      p <- weekly_snap(df, "PROT", wk)
      if (p < 0.3) "Complete" else if (p < init_prot*0.5) "Partial" else "None"
    }
    tab <- data.frame(
      Timepoint = c("Week 4","Week 12","Week 24","Week 52"),
      Proteinuria = round(c(weekly_snap(df,"PROT",4), weekly_snap(df,"PROT",12),
                            weekly_snap(df,"PROT",24), weekly_snap(df,"PROT",52)), 2),
      Remission = c(get_remission(4), get_remission(12), get_remission(24), get_remission(52))
    )
    datatable(tab, rownames=FALSE, options=list(dom="t"))
  })

  # ==========================================================================
  # TAB 3 OUTPUTS
  # ==========================================================================
  output$podPlot <- renderPlotly({
    df <- sim_r(); req(!is.null(df))
    plot_ly(df, x=~time_week, y=~POD*100, type="scatter", mode="lines",
            line=list(color="#E63946", width=2.5)) %>%
      add_lines(x=~time_week, y=rep(80, nrow(df)), line=list(color="orange",dash="dash"), name="Critical threshold") %>%
      layout(xaxis=list(title="Week"), yaxis=list(title="Podocyte Fraction (%)"),
             paper_bgcolor="#ffffff", plot_bgcolor="#fafafa")
  })

  output$fpePlot <- renderPlotly({
    df <- sim_r(); req(!is.null(df))
    plot_ly(df, x=~time_week, y=~FPE*100, type="scatter", mode="lines",
            line=list(color="#F4A261", width=2.5)) %>%
      layout(xaxis=list(title="Week"), yaxis=list(title="Foot Process Effacement (%)"),
             paper_bgcolor="#ffffff", plot_bgcolor="#fafafa")
  })

  output$compPlot <- renderPlotly({
    df <- sim_r(); req(!is.null(df))
    plot_ly(df, x=~time_week, y=~COMP, type="scatter", mode="lines",
            line=list(color="#3A86FF", width=2)) %>%
      layout(xaxis=list(title="Week"), yaxis=list(title="Complement (normalized)"),
             paper_bgcolor="#ffffff", plot_bgcolor="#fafafa")
  })

  output$inflamPlot <- renderPlotly({
    df <- sim_r(); req(!is.null(df))
    plot_ly(df, x=~time_week, y=~INFLAM, type="scatter", mode="lines",
            line=list(color="#8338EC", width=2)) %>%
      layout(xaxis=list(title="Week"), yaxis=list(title="Inflammation Index"),
             paper_bgcolor="#ffffff", plot_bgcolor="#fafafa")
  })

  output$tgfScarPlot <- renderPlotly({
    df <- sim_r(); req(!is.null(df))
    plot_ly(df, x=~time_week, y=~TGFb, name="TGF-β", type="scatter", mode="lines",
            line=list(color="#E63946", width=2)) %>%
      add_trace(y=~SCAR*5, name="Sclerosis (×5)", line=list(color="#264653",dash="dot",width=2)) %>%
      layout(xaxis=list(title="Week"), yaxis=list(title="Index / Sclerosis×5"),
             paper_bgcolor="#ffffff", plot_bgcolor="#fafafa",
             legend=list(orientation="h",y=-0.25))
  })

  # ==========================================================================
  # TAB 4 OUTPUTS
  # ==========================================================================
  output$gfrPlot <- renderPlotly({
    df <- sim_r(); req(!is.null(df))
    plot_ly(df, x=~time_week, y=~GFR_c, type="scatter", mode="lines",
            line=list(color=scenario_colors[input$scenario], width=2.5)) %>%
      add_lines(x=~time_week, y=rep(60, nrow(df)), line=list(color="orange",dash="dash"), name="CKD G2/G3 boundary") %>%
      add_lines(x=~time_week, y=rep(30, nrow(df)), line=list(color="red",dash="dot"), name="CKD G3b/G4 boundary") %>%
      layout(xaxis=list(title="Week"), yaxis=list(title="eGFR (mL/min/1.73m²)"),
             paper_bgcolor="#ffffff", plot_bgcolor="#fafafa",
             legend=list(orientation="h",y=-0.25))
  })

  output$remProbPlot <- renderPlotly({
    df <- sim_r(); req(!is.null(df))
    init_prot <- df$PROT[1]
    # P(CR) estimated from PROT normalized
    prob_cr <- pmax(0, pmin(1, (init_prot - df$PROT) / init_prot * 1.2 - 0.05))
    plot_ly(df, x=~time_week, y=prob_cr*100, type="scatter", mode="lines",
            line=list(color="#2A9D8F", width=2.5)) %>%
      layout(xaxis=list(title="Week"), yaxis=list(title="P(Complete Remission) %"),
             paper_bgcolor="#ffffff", plot_bgcolor="#fafafa")
  })

  output$esrdPlot <- renderPlotly({
    df <- sim_r(); req(!is.null(df))
    gfr_slope <- diff(df$GFR_c[c(1, nrow(df))]) / (364/365)  # mL/min/yr
    times_to_esrd <- ifelse(gfr_slope < 0, (df$GFR_c - 15) / abs(gfr_slope), NA)
    times_to_esrd <- pmax(0, pmin(times_to_esrd, 30))
    esrd_risk <- pmax(0, 1 - times_to_esrd/30)
    plot_ly(df, x=~time_week, y=esrd_risk*100, type="scatter", mode="lines",
            line=list(color="#E63946", width=2)) %>%
      layout(xaxis=list(title="Week"), yaxis=list(title="Estimated ESRD Risk (%)"),
             paper_bgcolor="#ffffff", plot_bgcolor="#fafafa")
  })

  output$outcomeTable <- renderDT({
    df <- sim_r(); req(!is.null(df))
    init_prot <- df$PROT[1]
    wks <- c(12, 24, 52)
    tab <- do.call(rbind, lapply(wks, function(wk) {
      p <- weekly_snap(df,"PROT",wk)
      g <- weekly_snap(df,"GFR_c",wk)
      s <- weekly_snap(df,"SCAR",wk)
      rem <- if(p<0.3) "CR" else if(p<init_prot*0.5) "PR" else "NR"
      data.frame(Week=wk, eGFR=round(g,1), Proteinuria=round(p,2),
                 Sclerosis=paste0(round(s*100,1),"%"), Remission=rem)
    }))
    datatable(tab, rownames=FALSE, options=list(dom="t"))
  })

  # ==========================================================================
  # TAB 5 OUTPUTS
  # ==========================================================================
  output$comparisonPanel <- renderPlotly({
    all_s <- all_scen_r()
    make_trace <- function(lst, yvar, label_fn=identity) {
      traces <- lapply(names(scenario_labels), function(sc) {
        df <- lst[[sc]]; req(!is.null(df))
        list(x=df$time_week, y=label_fn(df[[yvar]]),
             name=sc, line=list(color=scenario_colors[sc], width=1.8))
      })
      traces
    }
    # Build subplot
    p1 <- plot_ly(); p2 <- plot_ly(); p3 <- plot_ly(); p4 <- plot_ly()
    for (sc in names(scenario_labels)) {
      df <- all_s[[sc]]; if (is.null(df)) next
      p1 <- p1 %>% add_trace(x=df$time_week, y=df$PROT, name=sc, type="scatter",
                              mode="lines", line=list(color=scenario_colors[sc], width=1.8),
                              legendgroup=sc, showlegend=TRUE)
      p2 <- p2 %>% add_trace(x=df$time_week, y=df$GFR_c, name=sc, type="scatter",
                              mode="lines", line=list(color=scenario_colors[sc], width=1.8),
                              legendgroup=sc, showlegend=FALSE)
      p3 <- p3 %>% add_trace(x=df$time_week, y=df$POD*100, name=sc, type="scatter",
                              mode="lines", line=list(color=scenario_colors[sc], width=1.8),
                              legendgroup=sc, showlegend=FALSE)
      p4 <- p4 %>% add_trace(x=df$time_week, y=df$SCAR*100, name=sc, type="scatter",
                              mode="lines", line=list(color=scenario_colors[sc], width=1.8),
                              legendgroup=sc, showlegend=FALSE)
    }
    p1 <- p1 %>% layout(yaxis=list(title="Proteinuria (g/day)"), xaxis=list(title="Week"))
    p2 <- p2 %>% layout(yaxis=list(title="eGFR (mL/min/1.73m²)"), xaxis=list(title="Week"))
    p3 <- p3 %>% layout(yaxis=list(title="Podocyte Fraction (%)"), xaxis=list(title="Week"))
    p4 <- p4 %>% layout(yaxis=list(title="Sclerosis (%)"), xaxis=list(title="Week"))
    subplot(p1, p2, p3, p4, nrows=2, shareX=FALSE, titleY=TRUE, margin=0.07) %>%
      layout(paper_bgcolor="#ffffff", plot_bgcolor="#fafafa",
             legend=list(orientation="h", x=0, y=-0.08))
  })

  output$crBarPlot <- renderPlotly({
    all_s <- all_scen_r()
    cr_rates <- sapply(names(scenario_labels), function(sc) {
      df <- all_s[[sc]]
      if (is.null(df)) return(0)
      p <- weekly_snap(df, "PROT", 24)
      init_p <- df$PROT[1]
      ifelse(p < 0.3, 1, ifelse(p < init_p * 0.5, 0.5, 0)) * 100
    })
    sc_labels <- unname(scenario_labels)
    plot_ly(x=sc_labels, y=cr_rates, type="bar",
            marker=list(color=unname(scenario_colors)),
            text=paste0(round(cr_rates,0),"%"), textposition="outside") %>%
      layout(xaxis=list(title="", tickangle=-30),
             yaxis=list(title="Remission Score (%) at 6mo", range=c(0,120)),
             paper_bgcolor="#ffffff", plot_bgcolor="#fafafa")
  })

  output$gfrChangeBar <- renderPlotly({
    all_s <- all_scen_r()
    gfr_change <- sapply(names(scenario_labels), function(sc) {
      df <- all_s[[sc]]
      if (is.null(df)) return(0)
      weekly_snap(df,"GFR_c",52) - df$GFR_c[1]
    })
    bar_colors <- ifelse(gfr_change >= 0, "#2A9D8F", "#E63946")
    sc_labels <- unname(scenario_labels)
    plot_ly(x=sc_labels, y=gfr_change, type="bar",
            marker=list(color=bar_colors),
            text=paste0(round(gfr_change,1)," mL/min"), textposition="outside") %>%
      layout(xaxis=list(title="", tickangle=-30),
             yaxis=list(title="ΔeGFR at 52 weeks (mL/min)"),
             paper_bgcolor="#ffffff", plot_bgcolor="#fafafa")
  })

  # ==========================================================================
  # TAB 6 OUTPUTS
  # ==========================================================================
  output$predPK <- renderPlotly({
    df <- sim_r(); req(!is.null(df))
    df_sub <- df[df$time_day <= 56, ]  # first 8 weeks for clarity
    plot_ly(df_sub, x=~time_day, y=~PRED_conc, type="scatter", mode="lines",
            line=list(color="#F4A261", width=2)) %>%
      layout(xaxis=list(title="Days"), yaxis=list(title="Prednisolone (ng/mL)"),
             paper_bgcolor="#ffffff", plot_bgcolor="#fafafa")
  })

  output$tacPK <- renderPlotly({
    df <- sim_r(); req(!is.null(df))
    df_sub <- df[df$time_day <= 84, ]
    plot_ly(df_sub, x=~time_day, y=~TAC_conc, type="scatter", mode="lines",
            line=list(color="#2A9D8F", width=2), name="Tac level") %>%
      add_lines(x=df_sub$time_day, y=rep(5, nrow(df_sub)), line=list(color="green",dash="dash"), name="Lower TW (5 ng/mL)") %>%
      add_lines(x=df_sub$time_day, y=rep(10, nrow(df_sub)), line=list(color="red",dash="dot"), name="Upper TW (10 ng/mL)") %>%
      layout(xaxis=list(title="Days"), yaxis=list(title="Tacrolimus (ng/mL)"),
             paper_bgcolor="#ffffff", plot_bgcolor="#fafafa",
             legend=list(orientation="h",y=-0.3))
  })

  output$rtxPK <- renderPlotly({
    df <- sim_r(); req(!is.null(df))
    df_sub <- df[df$time_day <= 90, ]
    plot_ly(df_sub, x=~time_day, y=~RTX_conc, type="scatter", mode="lines",
            line=list(color="#264653", width=2)) %>%
      layout(xaxis=list(title="Days"), yaxis=list(title="Rituximab (mg/L)"),
             paper_bgcolor="#ffffff", plot_bgcolor="#fafafa")
  })

  output$sparsPK <- renderPlotly({
    df <- sim_r(); req(!is.null(df))
    df_sub <- df[df$time_day <= 14, ]
    plot_ly(df_sub, x=~time_day, y=~SPARS_conc, type="scatter", mode="lines",
            line=list(color="#8338EC", width=2)) %>%
      add_lines(x=df_sub$time_day, y=rep(200, nrow(df_sub)), line=list(color="green",dash="dash"), name="EC50") %>%
      layout(xaxis=list(title="Days (first 2 weeks)"), yaxis=list(title="Sparsentan (ng/mL)"),
             paper_bgcolor="#ffffff", plot_bgcolor="#fafafa")
  })

  output$biomarkerPlot <- renderPlotly({
    df <- sim_r(); req(!is.null(df))
    plot_ly(df, x=~time_week, y=~suPAR, name="suPAR (pg/mL)", type="scatter",
            mode="lines", line=list(color="#E63946", width=2)) %>%
      add_trace(y=~uPodocalyxin, name="u-Podocalyxin (pg/mL)",
                line=list(color="#3A86FF", width=2)) %>%
      add_trace(y=~NephrinShed*10, name="Nephrin Shedding ×10 (ng/mL)",
                line=list(color="#F4A261", width=2)) %>%
      layout(xaxis=list(title="Week"), yaxis=list(title="Biomarker Level"),
             paper_bgcolor="#ffffff", plot_bgcolor="#fafafa",
             legend=list(orientation="h", y=-0.25))
  })

  output$tacTA <- renderText({
    df <- sim_r(); req(!is.null(df))
    in_window <- mean(df$TAC_conc >= 5 & df$TAC_conc <= 10) * 100
    paste0("Tacrolimus time in therapeutic window (5-10 ng/mL):\n",
           round(in_window, 1), "%\n\n",
           ifelse(in_window > 50, "Target attainment: ADEQUATE",
                  ifelse(in_window > 20, "Target attainment: MARGINAL", "Target attainment: INSUFFICIENT")))
  })

  # Gauge output requires flexdashboard or custom - use verbatim fallback
  output$tacGauge <- renderPlotly({
    df <- sim_r(); req(!is.null(df))
    in_window <- mean(df$TAC_conc >= 5 & df$TAC_conc <= 10) * 100
    plot_ly(
      type = "indicator",
      mode = "gauge+number",
      value = round(in_window, 1),
      number = list(suffix = "%"),
      title = list(text = "Time in Window"),
      gauge = list(
        axis = list(range = list(0, 100)),
        bar = list(color = "#2A9D8F"),
        steps = list(
          list(range=c(0,30),  color="#FFCCCC"),
          list(range=c(30,60), color="#FFEECC"),
          list(range=c(60,100),color="#CCFFCC")
        ),
        threshold = list(line=list(color="red",width=3), thickness=0.8, value=50)
      )
    ) %>%
      layout(paper_bgcolor="#ffffff", margin=list(l=20,r=20,t=40,b=20))
  })
}

# =============================================================================
shinyApp(ui = ui, server = server)
