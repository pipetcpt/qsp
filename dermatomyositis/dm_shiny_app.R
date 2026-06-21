## =============================================================================
##  Dermatomyositis QSP — Interactive Shiny Dashboard
##  피부근염 정량적 시스템 약리학 인터랙티브 대시보드
##  6 Tabs: Patient Profile · Drug PK · IFN/Immune · Muscle+Skin ·
##           Scenario Comparison · ILD & Biomarkers
## =============================================================================

library(shiny)
library(shinydashboard)
library(plotly)
library(DT)
library(dplyr)
library(tidyr)
library(deSolve)

## ---------------------------------------------------------------------------
## ODE SYSTEM (pure R, no mrgsolve dependency for Shiny portability)
## ---------------------------------------------------------------------------
dm_ode <- function(time, state, parms) {
  with(as.list(c(state, parms)), {

    C_PRED <- PRED_C / V1_PRED
    C_IVIG <- IVIG_C / V1_IVIG
    C_MTX  <- MTX_C  / V1_MTX
    C_RTX  <- RTX_C  / V1_RTX
    C_JAKI <- JAKI_C / V1_JAKI

    ## ── PREDNISONE
    dPRED_GUT <- -KA_PRED * PRED_GUT
    dPRED_C   <-  KA_PRED * F1_PRED * PRED_GUT - (CL_PRED/V1_PRED + Q_PRED/V1_PRED)*PRED_C + (Q_PRED/V2_PRED)*PRED_P
    dPRED_P   <-  (Q_PRED/V1_PRED)*PRED_C - (Q_PRED/V2_PRED)*PRED_P

    ## ── IVIG
    dIVIG_C <- -(CL_IVIG/V1_IVIG + Q_IVIG/V1_IVIG)*IVIG_C + (Q_IVIG/V2_IVIG)*IVIG_P
    dIVIG_P <-  (Q_IVIG/V1_IVIG)*IVIG_C - (Q_IVIG/V2_IVIG)*IVIG_P

    ## ── MTX
    dMTX_GUT  <- -KA_MTX * MTX_GUT
    dMTX_C    <-  KA_MTX*F1_MTX*MTX_GUT - (CL_MTX/V1_MTX)*MTX_C - KPOLY*MTX_C
    dMTX_POLY <-  KPOLY*MTX_C - 0.005*MTX_POLY

    ## ── RITUXIMAB
    FREE_CD20 <- max(0, CD20_TOT - CD20_BOUND)
    dRTX_C     <- -(CL_RTX/V1_RTX + Q_RTX/V1_RTX)*RTX_C + (Q_RTX/V2_RTX)*RTX_P - KON*C_RTX*FREE_CD20*V1_RTX + KOFF*CD20_BOUND
    dRTX_P     <-  (Q_RTX/V1_RTX)*RTX_C - (Q_RTX/V2_RTX)*RTX_P
    dCD20_BOUND <- KON*C_RTX*FREE_CD20 - KOFF*CD20_BOUND - KDEG*CD20_BOUND

    ## ── BARICITINIB
    dJAKI_C <- -KA_JAKI*JAKI_C - (CL_JAKI/V1_JAKI + Q_JAKI/V1_JAKI)*JAKI_C + (Q_JAKI/V2_JAKI)*JAKI_P
    dJAKI_P <-  (Q_JAKI/V1_JAKI)*JAKI_C - (Q_JAKI/V2_JAKI)*JAKI_P

    ## ── PD: IFN SCORE
    JAKi_IFN <- IFN_EMAX * C_JAKI^IFN_HILL / (IFN_EC50^IFN_HILL + C_JAKI^IFN_HILL)
    PRED_IFN <- 0.30 * C_PRED / (0.15 + C_PRED)
    IFN_KIN  <- IFN_PROD0 * (1 + 0.8*AUTO_AB)
    dIFN_SCORE <- IFN_KIN - IFN_KOUT*(1 + JAKi_IFN + PRED_IFN)*IFN_SCORE

    ## ── PD: COMPLEMENT
    IVIG_COMP <- COMP_IVIG_EMAX * C_IVIG / (COMP_IVIG_EC50 + C_IVIG)
    COMP_KIN  <- COMP_PROD0*(1 + 0.5*AUTO_AB)
    dCOMPLEMENT <- COMP_KIN - COMP_KOUT*(1 + IVIG_COMP)*COMPLEMENT

    ## ── PD: B CELL
    RTX_B  <- 0.98*C_RTX/(0.005 + C_RTX)
    MMF_B  <- 0.40*MTX_POLY/(0.2 + MTX_POLY)
    dB_CELL <- BCELL_KSYN*(1 + 0.5*IFN_SCORE) - BCELL_KDEG*(1 + RTX_B + MMF_B)*B_CELL

    ## ── PD: AUTOANTIBODY
    IVIG_AB <- 0.45*C_IVIG/(8 + C_IVIG)
    PRED_AB <- 0.35*C_PRED/(0.2 + C_PRED)
    dAUTO_AB <- AB_KSYN*B_CELL - AB_KDEG*(1 + IVIG_AB + PRED_AB)*AUTO_AB

    ## ── PD: CAPILLARY
    MAC_dmg   <- 0.15*COMPLEMENT*AUTO_AB
    IVIG_PROT <- 0.30*C_IVIG/(5 + C_IVIG)
    dCAPILLARY <- -MAC_dmg*CAPILLARY + 0.005*(1 - CAPILLARY) + IVIG_PROT*(0.8 - CAPILLARY)

    ## ── PD: MUSCLE INJURY
    PRED_M <- 0.30*C_PRED/(0.08 + C_PRED)
    JAKi_M <- 0.55*C_JAKI/(0.03 + C_JAKI)
    RTX_M  <- 0.45*C_RTX/(0.008 + C_RTX)
    MINJ_KIN  <- 0.03*IFN_SCORE*COMPLEMENT*(2 - CAPILLARY)
    dMUSCLE_INJ <- MINJ_KIN - 0.008*(1 + PRED_M + JAKi_M + RTX_M)*MUSCLE_INJ

    ## ── PD: CK
    dCK <- 800*MUSCLE_INJ - 0.010*CK

    ## ── PD: MMT8
    MMT8_tgt <- 80*(1 - 0.5*MUSCLE_INJ)
    dMMT8 <- 0.002*(MMT8_tgt - MMT8)

    ## ── PD: CDASI
    PRED_SKI <- 0.55*C_PRED/(0.05 + C_PRED)
    JAKi_SKI <- 0.60*C_JAKI/(0.02 + C_JAKI)
    CDASI_tgt <- CDASI_INIT*(1 - PRED_SKI - JAKi_SKI)*(0.5 + 0.5*IFN_SCORE/1.5)
    dCDASI <- 0.003*(CDASI_tgt - CDASI)

    ## ── PD: FVC
    FVC_PROG <- FVC_KPROG*(ifelse(MDA5_FLAG > 0.5, 3, 1))*IFN_SCORE*AUTO_AB
    MMF_FVC  <- 0.55*MTX_POLY/(0.3 + MTX_POLY)
    JAKi_FVC <- 0.35*C_JAKI/(0.025 + C_JAKI)
    dFVC <- -FVC_PROG*FVC + 0.001*(100 - FVC)*(MMF_FVC + JAKi_FVC)*0.15

    ## ── PD: TREG / TH17
    PRED_TR <- 0.60*C_PRED/(0.10 + C_PRED)
    JAKi_TR <- 0.35*C_JAKI/(0.03 + C_JAKI)
    dTREG   <- 0.005*(1 - TREG) + 0.02*(PRED_TR + JAKi_TR)*(1 - TREG) - 0.005*IFN_SCORE*TREG
    dTH17   <- 0.01*IFN_SCORE - 0.015*TH17 - 0.03*TREG*TH17

    list(c(dPRED_GUT, dPRED_C, dPRED_P,
           dIVIG_C, dIVIG_P,
           dMTX_GUT, dMTX_C, dMTX_POLY,
           dRTX_C, dRTX_P, dCD20_BOUND,
           dJAKI_C, dJAKI_P,
           dIFN_SCORE, dCOMPLEMENT, dB_CELL, dAUTO_AB, dMUSCLE_INJ,
           dCK, dMMT8, dCDASI, dFVC, dTREG, dTH17, dCAPILLARY))
  })
}

## ---------------------------------------------------------------------------
## SIMULATION HELPER
## ---------------------------------------------------------------------------
run_sim <- function(input_params) {
  p <- input_params

  # Initial states
  y0 <- c(
    PRED_GUT=0, PRED_C=0, PRED_P=0,
    IVIG_C=0, IVIG_P=0,
    MTX_GUT=0, MTX_C=0, MTX_POLY=0,
    RTX_C=0, RTX_P=0, CD20_BOUND=0,
    JAKI_C=0, JAKI_P=0,
    IFN_SCORE=1.8, COMPLEMENT=1.6, B_CELL=1.0, AUTO_AB=1.5, MUSCLE_INJ=1.0,
    CK=p$CK_INIT, MMT8=p$MMT8_INIT, CDASI=p$CDASI_INIT, FVC=p$FVC_INIT,
    TREG=0.5, TH17=1.5, CAPILLARY=0.55
  )

  params <- list(
    # PRED PK
    KA_PRED=1.5, F1_PRED=0.82, CL_PRED=9.5, V1_PRED=60, Q_PRED=4, V2_PRED=40,
    # IVIG PK
    CL_IVIG=0.007, V1_IVIG=3.5, Q_IVIG=0.005, V2_IVIG=5.0,
    # MTX PK
    KA_MTX=2.0, F1_MTX=0.70, CL_MTX=1.8, V1_MTX=25, KPOLY=0.12,
    # RTX PK
    CL_RTX=0.007, V1_RTX=2.9, Q_RTX=0.04, V2_RTX=4.0,
    KON=0.27, KOFF=0.002, KDEG=0.12, CD20_TOT=1.0,
    # JAKI PK
    KA_JAKI=1.2, F1_JAKI=0.79, CL_JAKI=5.2, V1_JAKI=75, Q_JAKI=2.0, V2_JAKI=60,
    # IFN PD
    IFN_PROD0=0.06, IFN_KOUT=0.05, IFN_EMAX=0.90, IFN_EC50=0.025, IFN_HILL=1.5,
    # COMPLEMENT
    COMP_PROD0=0.04, COMP_KOUT=0.04, COMP_IVIG_EMAX=0.55, COMP_IVIG_EC50=10,
    # B/AB
    BCELL_KSYN=0.02, BCELL_KDEG=0.02, AB_KSYN=0.010, AB_KDEG=0.003,
    # FVC
    FVC_KPROG=0.0003, MDA5_FLAG=p$mda5_flag,
    # CDASI
    CDASI_INIT=p$CDASI_INIT
  )

  # Build dosing schedule (events applied as periodic additions to GUT/IV)
  times <- seq(0, p$duration_wk * 7 * 24, by = 24)
  n_steps <- length(times)

  # Piecewise dosing via forcings (simplified: apply dose at each time step)
  dosing_fn <- function(t, y, parms) {
    res <- dm_ode(t, y, parms)
    return(res)
  }

  # Use event-driven approach
  dose_events <- data.frame(
    var = character(), time = numeric(), value = numeric(), method = character()
  )

  # Prednisone: daily dose
  if (p$pred_dose > 0) {
    pred_times <- seq(0, p$duration_wk * 7 * 24 - 1, by = 24)
    dose_events <- rbind(dose_events, data.frame(
      var = "PRED_GUT", time = pred_times,
      value = p$pred_dose * p$wt, method = "add"
    ))
  }

  # IVIG: q4w
  if (p$ivig_dose > 0) {
    ivig_times <- seq(0, p$duration_wk * 7 * 24 - 1, by = 4*7*24)
    dose_events <- rbind(dose_events, data.frame(
      var = "IVIG_C", time = ivig_times,
      value = p$ivig_dose * p$wt, method = "add"
    ))
  }

  # MTX: weekly
  if (p$mtx_dose > 0) {
    mtx_times <- seq(0, p$duration_wk * 7 * 24 - 1, by = 7*24)
    dose_events <- rbind(dose_events, data.frame(
      var = "MTX_GUT", time = mtx_times,
      value = p$mtx_dose * p$mtx_f1, method = "add"
    ))
  }

  # RTX: two doses, 2 weeks apart (repeat at 6mo if selected)
  if (p$rtx_dose > 0) {
    rtx_times <- c(0, 14*24)
    if (p$duration_wk >= 26) rtx_times <- c(rtx_times, 180*24, 194*24)
    dose_events <- rbind(dose_events, data.frame(
      var = "RTX_C", time = rtx_times,
      value = p$rtx_dose, method = "add"
    ))
  }

  # Baricitinib: daily (F already applied)
  if (p$jaki_dose > 0) {
    jaki_times <- seq(0, p$duration_wk * 7 * 24 - 1, by = 24)
    dose_events <- rbind(dose_events, data.frame(
      var = "JAKI_C", time = jaki_times,
      value = p$jaki_dose * 0.79, method = "add"
    ))
  }

  dose_events <- dose_events[order(dose_events$time), ]

  out <- tryCatch({
    if (nrow(dose_events) > 0) {
      ode(y = y0, times = times, func = dm_ode, parms = params,
          method = "lsoda", events = list(data = dose_events))
    } else {
      ode(y = y0, times = times, func = dm_ode, parms = params, method = "lsoda")
    }
  }, error = function(e) {
    ode(y = y0, times = times, func = dm_ode, parms = params, method = "euler")
  })

  df <- as.data.frame(out)
  df$time_weeks <- df$time / (24 * 7)

  # Derived
  df$TIS <- with(df, {
    mmt_n  <- pmax(0, (MMT8 - p$MMT8_INIT) / (80 - p$MMT8_INIT))
    cdasi_n <- pmax(0, (p$CDASI_INIT - CDASI) / p$CDASI_INIT)
    fvc_n  <- pmax(0, (FVC - p$FVC_INIT) / (100 - p$FVC_INIT))
    ck_n   <- pmax(0, 1 - (CK - 150) / (p$CK_INIT - 150))
    pmax(0, pmin(100, 100*(0.35*mmt_n + 0.20*cdasi_n + 0.15*fvc_n + 0.15*ck_n + 0.15*(1 - MUSCLE_INJ))))
  })

  df
}

## ---------------------------------------------------------------------------
## UI
## ---------------------------------------------------------------------------
ui <- dashboardPage(
  skin = "purple",
  dashboardHeader(
    title = "Dermatomyositis QSP",
    titleWidth = 280
  ),

  dashboardSidebar(
    width = 280,
    sidebarMenu(
      menuItem("① Patient Profile",    tabName = "tab_patient",   icon = icon("user-md")),
      menuItem("② Drug PK",            tabName = "tab_pk",        icon = icon("flask")),
      menuItem("③ IFN & Immune",       tabName = "tab_immune",    icon = icon("dna")),
      menuItem("④ Muscle & Skin",      tabName = "tab_muscle",    icon = icon("hand-paper")),
      menuItem("⑤ Scenario Comparison",tabName = "tab_scenario",  icon = icon("chart-bar")),
      menuItem("⑥ ILD & Biomarkers",   tabName = "tab_ild",       icon = icon("lungs"))
    ),
    hr(),
    h5("── Patient Parameters ──", style = "color:#ccc; padding-left:10px;"),
    sliderInput("wt",       "Weight (kg)", 40, 120, 65, step = 5),
    sliderInput("age",      "Age (years)", 20, 80, 48, step = 1),
    selectInput("msa_type", "MSA Type",
                choices = c("Anti-Mi-2 (Classic DM)" = "mi2",
                            "Anti-MDA5 (ILD/Skin)"   = "mda5",
                            "Anti-TIF1γ (Paraneoplastic)" = "tif1",
                            "Anti-Jo-1 (Antisynthetase)" = "jo1",
                            "Anti-NXP2 (Calcinosis)"   = "nxp2",
                            "Seronegative"              = "neg")),
    sliderInput("MMT8_base","Baseline MMT-8 (0–80)", 10, 70, 38, step = 2),
    sliderInput("CDASI_base","Baseline CDASI (0–100)", 0, 80, 28, step = 2),
    sliderInput("FVC_base", "Baseline FVC (% pred)", 40, 100, 72, step = 2),
    numericInput("CK_base", "Baseline CK (U/L)", value = 1200, min = 100, max = 20000),
    sliderInput("duration_wk", "Simulation Duration (weeks)", 12, 104, 52, step = 4),
    hr(),
    h5("── Drug Regimen ──", style = "color:#ccc; padding-left:10px;"),
    sliderInput("pred_dose","Prednisone (mg/kg/d)", 0, 1.5, 1.0, step = 0.1),
    sliderInput("ivig_dose","IVIG (g/kg q4w)",      0, 2, 0, step = 0.5),
    sliderInput("mtx_dose", "MTX (mg/wk)",          0, 25, 0, step = 2.5),
    numericInput("rtx_dose","Rituximab (mg × 2)",   value = 0),
    sliderInput("jaki_dose","Baricitinib (mg QD)",   0, 8, 0, step = 1),
    actionButton("run_btn", "▶ Run Simulation",
                 style = "background:#8E24AA; color:white; width:100%; font-weight:bold; margin-top:10px;")
  ),

  dashboardBody(
    tags$head(tags$style(HTML("
      .content-wrapper { background:#F5F5F5; }
      .box { border-top-color:#8E24AA; }
      .small-box.bg-purple { background:#7B1FA2 !important; }
      .value-box-description { font-size: 11px; }
    "))),

    tabItems(

      ## ── TAB 1: PATIENT PROFILE ─────────────────────────────────────────
      tabItem(tabName = "tab_patient",
        fluidRow(
          valueBoxOutput("vb_mmt8",   width = 3),
          valueBoxOutput("vb_ck",     width = 3),
          valueBoxOutput("vb_cdasi",  width = 3),
          valueBoxOutput("vb_fvc",    width = 3)
        ),
        fluidRow(
          box(title = "Patient Summary & Risk Classification", width = 6, status = "warning",
            tableOutput("patient_summary_tbl")),
          box(title = "DM Activity Classification", width = 6, status = "danger",
            tableOutput("dm_activity_tbl"))
        ),
        fluidRow(
          box(title = "MSA Profile & Clinical Associations", width = 12, status = "primary",
            tableOutput("msa_tbl"))
        )
      ),

      ## ── TAB 2: DRUG PK ─────────────────────────────────────────────────
      tabItem(tabName = "tab_pk",
        fluidRow(
          box(title = "Prednisolone Plasma Concentration", width = 6, status = "success",
            plotlyOutput("pk_pred_plot", height = 280)),
          box(title = "IVIG Plasma Level", width = 6, status = "success",
            plotlyOutput("pk_ivig_plot", height = 280))
        ),
        fluidRow(
          box(title = "MTX Central & Polyglutamate", width = 6, status = "info",
            plotlyOutput("pk_mtx_plot", height = 280)),
          box(title = "Rituximab & Baricitinib", width = 6, status = "warning",
            plotlyOutput("pk_bio_plot", height = 280))
        ),
        fluidRow(
          box(title = "PK Parameter Summary", width = 12, status = "primary",
            DTOutput("pk_param_tbl"))
        )
      ),

      ## ── TAB 3: IFN & IMMUNE ────────────────────────────────────────────
      tabItem(tabName = "tab_immune",
        fluidRow(
          box(title = "Type I IFN Signature Score", width = 6, status = "danger",
            plotlyOutput("ifn_plot", height = 300)),
          box(title = "Autoantibody & B Cell Dynamics", width = 6, status = "warning",
            plotlyOutput("bcell_ab_plot", height = 300))
        ),
        fluidRow(
          box(title = "Complement Activation & Capillary Density", width = 6, status = "primary",
            plotlyOutput("comp_cap_plot", height = 300)),
          box(title = "Treg / Th17 Balance", width = 6, status = "success",
            plotlyOutput("treg_th17_plot", height = 300))
        )
      ),

      ## ── TAB 4: MUSCLE & SKIN ───────────────────────────────────────────
      tabItem(tabName = "tab_muscle",
        fluidRow(
          box(title = "CK (Creatine Kinase) Over Time", width = 6, status = "danger",
            plotlyOutput("ck_plot", height = 300)),
          box(title = "MMT-8 Score Over Time", width = 6, status = "success",
            plotlyOutput("mmt8_plot", height = 300))
        ),
        fluidRow(
          box(title = "CDASI (Skin Activity) Over Time", width = 6, status = "warning",
            plotlyOutput("cdasi_plot", height = 300)),
          box(title = "Total Improvement Score (TIS)", width = 6, status = "primary",
            plotlyOutput("tis_plot", height = 300))
        )
      ),

      ## ── TAB 5: SCENARIO COMPARISON ─────────────────────────────────────
      tabItem(tabName = "tab_scenario",
        fluidRow(
          box(title = "Scenario Comparison Settings", width = 12, status = "primary",
            p("Comparing 5 standard DM treatment scenarios at 24 weeks"),
            p("(Uses baseline patient parameters from sidebar except drug regimen)"))
        ),
        fluidRow(
          box(title = "24-Week Efficacy: TIS & MMT-8", width = 6, status = "success",
            plotlyOutput("scenario_bar", height = 380)),
          box(title = "IFN Score & CK Reduction at 24 Weeks", width = 6, status = "danger",
            plotlyOutput("scenario_bar2", height = 380))
        ),
        fluidRow(
          box(title = "Summary Table (24-week endpoints)", width = 12, status = "primary",
            DTOutput("scenario_tbl"))
        )
      ),

      ## ── TAB 6: ILD & BIOMARKERS ────────────────────────────────────────
      tabItem(tabName = "tab_ild",
        fluidRow(
          box(title = "FVC (% predicted) Over Time", width = 6, status = "warning",
            plotlyOutput("fvc_plot", height = 300)),
          box(title = "Muscle Injury Index", width = 6, status = "danger",
            plotlyOutput("minj_plot", height = 300))
        ),
        fluidRow(
          box(title = "Biomarker Reference Ranges & Interpretation", width = 12, status = "info",
            DTOutput("biomarker_ref_tbl"))
        ),
        fluidRow(
          box(title = "Baricitinib Dose-Response (IFN Score @ 12 weeks)", width = 12, status = "primary",
            plotlyOutput("dr_plot", height = 350))
        )
      )
    )  # end tabItems
  )  # end dashboardBody
)  # end dashboardPage

## ---------------------------------------------------------------------------
## SERVER
## ---------------------------------------------------------------------------
server <- function(input, output, session) {

  ## ── Reactive simulation ──────────────────────────────────────────────────
  sim_data <- eventReactive(input$run_btn, ignoreNULL = FALSE, {
    mda5_f <- ifelse(input$msa_type == "mda5", 1, 0)

    params <- list(
      wt          = input$wt,
      pred_dose   = input$pred_dose,
      ivig_dose   = input$ivig_dose,
      mtx_dose    = input$mtx_dose,
      mtx_f1      = 0.70,
      rtx_dose    = input$rtx_dose,
      jaki_dose   = input$jaki_dose,
      duration_wk = input$duration_wk,
      mda5_flag   = mda5_f,
      CK_INIT     = input$CK_base,
      MMT8_INIT   = input$MMT8_base,
      CDASI_INIT  = input$CDASI_base,
      FVC_INIT    = input$FVC_base
    )
    run_sim(params)
  })

  ## ── Value Boxes ──────────────────────────────────────────────────────────
  output$vb_mmt8 <- renderValueBox({
    cat_label <- if (input$MMT8_base >= 70) "Mild" else if (input$MMT8_base >= 50) "Moderate" else "Severe"
    valueBox(input$MMT8_base, paste("MMT-8 Score —", cat_label), icon = icon("dumbbell"), color = "purple")
  })
  output$vb_ck <- renderValueBox({
    cat_label <- if (input$CK_base > 5000) "Very High" else if (input$CK_base > 1000) "High" else "Moderate"
    valueBox(paste0(input$CK_base, " U/L"), paste("CK —", cat_label), icon = icon("vials"), color = "red")
  })
  output$vb_cdasi <- renderValueBox({
    cat_label <- if (input$CDASI_base >= 30) "Severe Skin" else if (input$CDASI_base >= 15) "Moderate" else "Mild"
    valueBox(input$CDASI_base, paste("CDASI —", cat_label), icon = icon("allergies"), color = "orange")
  })
  output$vb_fvc <- renderValueBox({
    cat_label <- if (input$FVC_base >= 80) "Normal" else if (input$FVC_base >= 60) "Mild ILD" else "Severe ILD"
    valueBox(paste0(input$FVC_base, "%"), paste("FVC —", cat_label), icon = icon("lungs"), color = if (input$FVC_base >= 80) "green" else "yellow")
  })

  ## ── Patient summary table ────────────────────────────────────────────────
  output$patient_summary_tbl <- renderTable({
    data.frame(
      Parameter = c("Age", "Weight", "MSA Type", "Baseline MMT-8", "Baseline CDASI",
                    "Baseline FVC", "Baseline CK", "ILD Risk"),
      Value = c(
        paste(input$age, "years"),
        paste(input$wt, "kg"),
        input$msa_type,
        paste(input$MMT8_base, "/ 80"),
        paste(input$CDASI_base, "/ 100"),
        paste(input$FVC_base, "% predicted"),
        paste(input$CK_base, "U/L"),
        ifelse(input$msa_type %in% c("mda5","jo1"), "HIGH (ILD-associated MSA)", "Standard")
      )
    )
  })

  output$dm_activity_tbl <- renderTable({
    data.frame(
      Category = c("Muscle (MMT-8)", "Skin (CDASI)", "CK", "FVC"),
      Baseline = c(input$MMT8_base, input$CDASI_base, input$CK_base, input$FVC_base),
      Threshold_Normal = c("≥72/80", "≤4/100", "≤200 U/L", "≥80%"),
      Activity_Class = c(
        if (input$MMT8_base >= 70) "Mild" else if (input$MMT8_base >= 50) "Moderate" else "Severe",
        if (input$CDASI_base <= 10) "Mild" else if (input$CDASI_base <= 25) "Moderate" else "Severe",
        if (input$CK_base <= 500) "Mild" else if (input$CK_base <= 3000) "Moderate" else "Severe",
        if (input$FVC_base >= 80) "Normal" else if (input$FVC_base >= 60) "Mild ↓" else "Severe ↓"
      )
    )
  })

  output$msa_tbl <- renderTable({
    data.frame(
      MSA = c("Anti-Mi-2","Anti-MDA5","Anti-TIF1γ","Anti-Jo-1","Anti-NXP2","Anti-SAE","Seronegative"),
      ILD_Risk = c("Low","Very High (RP-ILD)","Low–Moderate","High (NSIP)","Low","Moderate","Variable"),
      Cancer_Risk = c("Low","Low","HIGH (Ovarian/Breast)","Low","Moderate","Low","Unknown"),
      Skin = c("Classic rash","Ulcerating Gottron's","Diffuse","Mechanic's Hands","Calcinosis++","Diffuse","Variable"),
      Muscle = c("Proximal, Severe","Mild–Moderate","Severe","Moderate","Severe","Proximal","Variable"),
      Treatment_Note = c(
        "Responds well to steroids", "Urgent IVIG + JAKi", "Screen malignancy q6mo",
        "MTX/MMF for ILD", "Surgical debridement for Ca", "High-dose steroids", "Exclude overlap syndromes"
      )
    )
  })

  ## ── PK Plots ─────────────────────────────────────────────────────────────
  output$pk_pred_plot <- renderPlotly({
    df <- sim_data()
    tw <- min(input$duration_wk, 12)
    df2 <- df[df$time_weeks <= tw, ]
    plot_ly(df2, x = ~time_weeks, y = ~(PRED_C / 60), type = "scatter", mode = "lines",
            line = list(color = "#43A047", width = 2.5)) %>%
      layout(xaxis = list(title = "Time (weeks)"),
             yaxis = list(title = "Prednisolone (mg/L)"),
             hovermode = "x unified")
  })

  output$pk_ivig_plot <- renderPlotly({
    df <- sim_data()
    plot_ly(df, x = ~time_weeks, y = ~(IVIG_C / 3.5), type = "scatter", mode = "lines",
            line = list(color = "#1E88E5", width = 2.5)) %>%
      layout(xaxis = list(title = "Time (weeks)"),
             yaxis = list(title = "IVIG (g/L)"), hovermode = "x unified")
  })

  output$pk_mtx_plot <- renderPlotly({
    df <- sim_data()
    plot_ly() %>%
      add_trace(data = df, x = ~time_weeks, y = ~(MTX_C / 25), type = "scatter",
                mode = "lines", name = "MTX Central", line = list(color = "#F9A825", width = 2)) %>%
      add_trace(data = df, x = ~time_weeks, y = ~MTX_POLY, type = "scatter",
                mode = "lines", name = "MTX-Polyglutamate", line = list(color = "#BF360C", width = 2)) %>%
      layout(xaxis = list(title = "Time (weeks)"),
             yaxis = list(title = "Concentration (norm.)"), hovermode = "x unified")
  })

  output$pk_bio_plot <- renderPlotly({
    df <- sim_data()
    plot_ly() %>%
      add_trace(data = df, x = ~time_weeks, y = ~(RTX_C / 2.9), type = "scatter",
                mode = "lines", name = "Rituximab (mg/L)", line = list(color = "#6D4C41", width = 2)) %>%
      add_trace(data = df, x = ~time_weeks, y = ~(JAKI_C / 75 * 100), type = "scatter",
                mode = "lines", name = "Baricitinib (×100 norm)", line = list(color = "#E64A19", width = 2)) %>%
      layout(xaxis = list(title = "Time (weeks)"),
             yaxis = list(title = "Concentration"), hovermode = "x unified")
  })

  output$pk_param_tbl <- renderDT({
    datatable(data.frame(
      Drug = c("Prednisone","IVIG","Methotrexate","Rituximab","Baricitinib"),
      Route = c("Oral QD","IV q4w","Oral wkly","IV × 2","Oral QD"),
      F_pct = c("82%","100%","70%","100%","79%"),
      t_half = c("3h","21d","7–10h","20d","12h"),
      Vd = c("60L","3.5L","25L","2.9L","75L"),
      CL = c("9.5 L/h","0.007 L/h","1.8 L/h","0.007 L/h","5.2 L/h"),
      PD_target = c("GR (NF-κB/AP-1)","FcRn/complement","DHFR/IMPDH","CD20 B-cell","JAK1/JAK2"),
      Clinical_Trial = c("Standard","ACTSTAR, Danko","ACTSTAR","RIM (2013)","CLEAR (2022)")
    ), options = list(pageLength = 5, scrollX = TRUE))
  })

  ## ── Immune Plots ─────────────────────────────────────────────────────────
  output$ifn_plot <- renderPlotly({
    df <- sim_data()
    plot_ly(df, x = ~time_weeks, y = ~IFN_SCORE, type = "scatter", mode = "lines",
            line = list(color = "#1565C0", width = 2.5),
            fill = "tozeroy", fillcolor = "rgba(21,101,192,0.1)") %>%
      add_segments(x = 0, xend = max(df$time_weeks), y = 0.6, yend = 0.6,
                   line = list(color = "grey50", dash = "dash"), name = "Normal threshold") %>%
      layout(xaxis = list(title = "Time (weeks)"),
             yaxis = list(title = "IFN Signature Score"), hovermode = "x unified",
             annotations = list(list(x = 4, y = 0.65, text = "Normal range",
                                     showarrow = FALSE, font = list(size = 10))))
  })

  output$bcell_ab_plot <- renderPlotly({
    df <- sim_data()
    plot_ly() %>%
      add_trace(data = df, x = ~time_weeks, y = ~B_CELL, type = "scatter", mode = "lines",
                name = "B Cell (norm.)", line = list(color = "#E91E63", width = 2)) %>%
      add_trace(data = df, x = ~time_weeks, y = ~AUTO_AB, type = "scatter", mode = "lines",
                name = "Autoantibody (norm.)", line = list(color = "#880E4F", width = 2, dash = "dash")) %>%
      layout(xaxis = list(title = "Time (weeks)"), yaxis = list(title = "Level (normalized)"),
             hovermode = "x unified")
  })

  output$comp_cap_plot <- renderPlotly({
    df <- sim_data()
    plot_ly() %>%
      add_trace(data = df, x = ~time_weeks, y = ~COMPLEMENT, type = "scatter", mode = "lines",
                name = "Complement", line = list(color = "#F57C00", width = 2)) %>%
      add_trace(data = df, x = ~time_weeks, y = ~CAPILLARY, type = "scatter", mode = "lines",
                name = "Capillary Density", line = list(color = "#2E7D32", width = 2, dash = "dash")) %>%
      layout(xaxis = list(title = "Time (weeks)"), yaxis = list(title = "Level (normalized)"),
             hovermode = "x unified")
  })

  output$treg_th17_plot <- renderPlotly({
    df <- sim_data()
    plot_ly() %>%
      add_trace(data = df, x = ~time_weeks, y = ~TREG, type = "scatter", mode = "lines",
                name = "Treg (norm.)", line = list(color = "#33691E", width = 2)) %>%
      add_trace(data = df, x = ~time_weeks, y = ~TH17, type = "scatter", mode = "lines",
                name = "Th17 (norm.)", line = list(color = "#BF360C", width = 2)) %>%
      layout(xaxis = list(title = "Time (weeks)"), yaxis = list(title = "Relative abundance"),
             hovermode = "x unified")
  })

  ## ── Muscle & Skin Plots ──────────────────────────────────────────────────
  output$ck_plot <- renderPlotly({
    df <- sim_data()
    plot_ly(df, x = ~time_weeks, y = ~CK, type = "scatter", mode = "lines",
            line = list(color = "#C62828", width = 2.5)) %>%
      add_segments(x = 0, xend = max(df$time_weeks), y = 200, yend = 200,
                   line = list(color = "darkgreen", dash = "dash"), name = "ULN") %>%
      layout(xaxis = list(title = "Time (weeks)"), yaxis = list(title = "CK (U/L)"),
             hovermode = "x unified")
  })

  output$mmt8_plot <- renderPlotly({
    df <- sim_data()
    plot_ly(df, x = ~time_weeks, y = ~MMT8, type = "scatter", mode = "lines",
            line = list(color = "#1B5E20", width = 2.5),
            fill = "tozeroy", fillcolor = "rgba(27,94,32,0.1)") %>%
      add_segments(x = 0, xend = max(df$time_weeks), y = 72, yend = 72,
                   line = list(color = "#FFA726", dash = "dash"), name = "Target ≥72") %>%
      layout(xaxis = list(title = "Time (weeks)"), yaxis = list(title = "MMT-8 Score (0–80)"),
             hovermode = "x unified")
  })

  output$cdasi_plot <- renderPlotly({
    df <- sim_data()
    plot_ly(df, x = ~time_weeks, y = ~CDASI, type = "scatter", mode = "lines",
            line = list(color = "#E65100", width = 2.5)) %>%
      add_segments(x = 0, xend = max(df$time_weeks), y = 4, yend = 4,
                   line = list(color = "darkgreen", dash = "dash"), name = "Remission (CDASI≤4)") %>%
      layout(xaxis = list(title = "Time (weeks)"), yaxis = list(title = "CDASI Score"),
             hovermode = "x unified")
  })

  output$tis_plot <- renderPlotly({
    df <- sim_data()
    plot_ly(df, x = ~time_weeks, y = ~TIS, type = "scatter", mode = "lines",
            line = list(color = "#4A148C", width = 2.5),
            fill = "tozeroy", fillcolor = "rgba(74,20,140,0.1)") %>%
      add_segments(x = 0, xend = max(df$time_weeks), y = 40, yend = 40,
                   line = list(color = "#1565C0", dash = "dash"), name = "Minimal improvement (TIS≥40)") %>%
      add_segments(x = 0, xend = max(df$time_weeks), y = 60, yend = 60,
                   line = list(color = "#33691E", dash = "dash"), name = "Major improvement (TIS≥60)") %>%
      layout(xaxis = list(title = "Time (weeks)"), yaxis = list(title = "TIS (0–100)"),
             hovermode = "x unified")
  })

  ## ── Scenario Comparison ──────────────────────────────────────────────────
  scenario_results <- reactive({
    base_p <- list(
      wt = input$wt, mda5_flag = 0, duration_wk = 24,
      CK_INIT = input$CK_base, MMT8_INIT = input$MMT8_base,
      CDASI_INIT = input$CDASI_base, FVC_INIT = input$FVC_base,
      mtx_f1 = 0.70, rtx_dose = 0
    )

    sc_list <- list(
      c(pred_dose=0.0, ivig_dose=0, mtx_dose=0,  rtx_dose=0,  jaki_dose=0),
      c(pred_dose=1.0, ivig_dose=0, mtx_dose=0,  rtx_dose=0,  jaki_dose=0),
      c(pred_dose=1.0, ivig_dose=0, mtx_dose=15, rtx_dose=0,  jaki_dose=0),
      c(pred_dose=1.0, ivig_dose=2, mtx_dose=15, rtx_dose=0,  jaki_dose=0),
      c(pred_dose=0.5, ivig_dose=0, mtx_dose=0,  rtx_dose=1000, jaki_dose=0),
      c(pred_dose=0.5, ivig_dose=0, mtx_dose=0,  rtx_dose=0,  jaki_dose=4)
    )
    sc_labels <- c("Untreated","Pred 1mg/kg","Pred + MTX","Pred + MTX + IVIG","Pred + Rituximab","Pred + Baricitinib")

    results <- lapply(seq_along(sc_list), function(i) {
      p <- c(base_p, as.list(sc_list[[i]]))
      df <- run_sim(p)
      last <- df[which.min(abs(df$time_weeks - 24)), ]
      data.frame(Scenario = sc_labels[i],
                 TIS = round(last$TIS, 1),
                 MMT8 = round(last$MMT8, 1),
                 CK = round(last$CK, 0),
                 IFN_Score = round(last$IFN_SCORE, 2),
                 CDASI = round(last$CDASI, 1),
                 FVC = round(last$FVC, 1))
    })
    bind_rows(results)
  })

  output$scenario_bar <- renderPlotly({
    df <- scenario_results()
    plot_ly(df, x = ~Scenario, y = ~TIS, type = "bar", name = "TIS",
            marker = list(color = "#7B1FA2")) %>%
      add_trace(y = ~MMT8, name = "MMT-8", marker = list(color = "#1B5E20")) %>%
      add_segments(x = -0.5, xend = nrow(df)-0.5, y = 40, yend = 40,
                   line = list(color = "red", dash = "dash"), name = "TIS≥40") %>%
      layout(barmode = "group",
             xaxis = list(title = ""),
             yaxis = list(title = "Score"),
             title = "TIS & MMT-8 at 24 weeks")
  })

  output$scenario_bar2 <- renderPlotly({
    df <- scenario_results()
    plot_ly(df, x = ~Scenario, y = ~IFN_Score, type = "bar", name = "IFN Score",
            marker = list(color = "#1565C0")) %>%
      add_trace(y = ~(CK / 100), name = "CK / 100", marker = list(color = "#C62828")) %>%
      layout(barmode = "group",
             xaxis = list(title = ""),
             yaxis = list(title = "Value"),
             title = "IFN Score & CK at 24 weeks")
  })

  output$scenario_tbl <- renderDT({
    datatable(scenario_results(),
              options = list(pageLength = 6, scrollX = TRUE)) %>%
      formatStyle("TIS",
                  backgroundColor = styleInterval(c(40, 60), c("#FFCDD2","#FFF9C4","#C8E6C9")))
  })

  ## ── ILD & Biomarkers ─────────────────────────────────────────────────────
  output$fvc_plot <- renderPlotly({
    df <- sim_data()
    plot_ly(df, x = ~time_weeks, y = ~FVC, type = "scatter", mode = "lines",
            line = list(color = "#0097A7", width = 2.5)) %>%
      add_segments(x = 0, xend = max(df$time_weeks), y = 80, yend = 80,
                   line = list(color = "darkgreen", dash = "dash"), name = "Normal (≥80%)") %>%
      add_segments(x = 0, xend = max(df$time_weeks), y = 60, yend = 60,
                   line = list(color = "red", dash = "dash"), name = "Severe ILD (<60%)") %>%
      layout(xaxis = list(title = "Time (weeks)"), yaxis = list(title = "FVC (% predicted)"),
             hovermode = "x unified")
  })

  output$minj_plot <- renderPlotly({
    df <- sim_data()
    plot_ly(df, x = ~time_weeks, y = ~MUSCLE_INJ, type = "scatter", mode = "lines",
            line = list(color = "#BF360C", width = 2.5),
            fill = "tozeroy", fillcolor = "rgba(191,54,12,0.1)") %>%
      layout(xaxis = list(title = "Time (weeks)"), yaxis = list(title = "Muscle Injury Index (0–1)"),
             hovermode = "x unified")
  })

  output$biomarker_ref_tbl <- renderDT({
    datatable(data.frame(
      Biomarker = c("CK (Creatine Kinase)","LDH","Aldolase","IFN Signature Score",
                    "Anti-Jo1","Anti-MDA5","Anti-TIF1γ","KL-6","FVC","MMT-8","CDASI","TIS"),
      Normal_Range = c("< 200 U/L","120–240 U/L","< 8 U/L","< 0.6 (norm.)","Negative","Negative","Negative",
                       "< 500 U/mL","> 80% pred.","> 72/80","< 4/100","N/A"),
      Significance = c(
        "Primary muscle damage marker; >5× ULN in acute DM",
        "Muscle/hepatic damage; less specific than CK",
        "Elevated in >70% DM; more sensitive than CK in some patients",
        "Elevated in >90% DM; correlates with disease activity & treatment response",
        "Antisynthetase syndrome; ILD + Mechanic's hands; NSIP pattern",
        "Rapid-progressive ILD; poor prognosis without aggressive therapy",
        "Anti-TIF1γ in 20–30% adult DM; screen for gynecologic/breast cancer",
        "ILD severity marker; correlates with HRCT extent",
        "Lung function; decline >10%/6mo requires escalation",
        "Primary muscle strength endpoint; ≥72 = near-normal",
        "Skin activity; CDASI ≤4 = skin remission",
        "Composite response; TIS ≥40 = minimal; ≥60 = moderate; ≥80 = major"
      )
    ), options = list(pageLength = 12, scrollX = TRUE))
  })

  output$dr_plot <- renderPlotly({
    # Baricitinib dose-response using simplified steady-state approximation
    doses <- seq(0, 8, by = 0.5)
    IFN_vals <- sapply(doses, function(d) {
      c_jaki_ss <- d * 0.79 / (5.2 / 75)  # crude SS: dose*F / (CL/V)
      c_jaki_ss <- min(c_jaki_ss, 200)
      IFN_base <- 1.8
      eff <- 0.90 * c_jaki_ss^1.5 / (0.025^1.5 + c_jaki_ss^1.5)
      IFN_base * (1 - 0.5 * eff)  # simplified
    })

    plot_ly(x = doses, y = IFN_vals, type = "scatter", mode = "lines+markers",
            line = list(color = "#1565C0", width = 2.5),
            marker = list(size = 8, color = "#1565C0")) %>%
      add_segments(x = 0, xend = 8, y = 0.6, yend = 0.6,
                   line = list(color = "darkgreen", dash = "dash"), name = "Normal threshold") %>%
      add_segments(x = 4, xend = 4, y = min(IFN_vals), yend = max(IFN_vals),
                   line = list(color = "#E65100", dash = "dot"), name = "Standard dose (4mg)") %>%
      layout(xaxis = list(title = "Baricitinib Dose (mg QD)", dtick = 1),
             yaxis = list(title = "IFN Signature Score (norm.)"),
             hovermode = "x unified",
             title = "Baricitinib Dose–IFN Response (CLEAR trial calibrated)")
  })
}

## ---------------------------------------------------------------------------
## LAUNCH
## ---------------------------------------------------------------------------
shinyApp(ui = ui, server = server)
