## =============================================================================
## MCTD QSP Shiny Interactive Dashboard
## Mixed Connective Tissue Disease — Quantitative Systems Pharmacology
## =============================================================================
## Tabs:
##   1. Patient Profile & Disease Characterization
##   2. Drug PK — Concentration-Time Profiles
##   3. Immunological PD — Cytokines, Autoantibodies, Immune Cells
##   4. Organ Function — Lung (ILD/PAH), Muscle, Joints
##   5. Treatment Scenario Comparison
##   6. Biomarker Heatmap & Risk Stratification
## =============================================================================

library(shiny)
library(shinydashboard)
library(ggplot2)
library(dplyr)
library(tidyr)
library(DT)
library(plotly)

## ─── Simplified ODE solver (without mrgsolve for portability) ────────────────
solve_mctd <- function(
    dose_hcq    = 0,
    dose_pred   = 0,
    dose_mmf    = 0,
    dose_rtx    = 0,
    dose_bos    = 0,
    duration    = 365,
    wt          = 60,
    # Disease severity modifiers
    Ab0         = 300,
    CK0         = 800,
    FVC0        = 75,
    RVSP0       = 45,
    PVR0        = 4.0
) {
  dt  <- 1     # days
  N   <- duration

  # PK parameters
  HCQ_CL <- 19.5; HCQ_Vd <- 800; HCQ_Ka <- 8.4   # Ka in /day
  MPA_CL <- 11.6 * 24; MPA_Vd <- 3.6; MPA_Ka <- 26.4
  PRED_CL<- 12.0 * 24; PRED_Vd<- 0.55; PRED_Ka<- 20.4
  RTX_CL <- 0.24; RTX_Vd <- 3.5
  BOS_CL <- 2.0 * 24; BOS_Vd <- 0.18; BOS_Ka <- 19.2

  # State vectors
  t_seq <- seq(0, N, by = dt)
  n     <- length(t_seq)

  # Initialize
  vars <- list(
    HCQ_gut=0, HCQ_cp=0,
    MPA_gut=0, MPA_cp=0,
    PRED_gut=0, PRED_cp=0,
    RTX_cp=0,
    BOS_gut=0, BOS_cp=0,
    Th1=150, Th17=60, Treg=80,
    Bnai=200, PC=20,
    AbU1=Ab0, TNF=25, IL6=15, IL17=18,
    IFNg=8, TGFb=5, IFNa=20,
    ET1=2.5, PVR=PVR0,
    FVC=FVC0, DLCO=65, Collagen=0.3,
    CK=CK0, MMT=70, SJC=6, C3=80, C4=15
  )

  results <- data.frame(time = numeric(n))
  for (nm in names(vars)) results[[nm]] <- NA_real_
  results$time <- t_seq
  for (nm in names(vars)) results[1, nm] <- vars[[nm]]

  # Drug effect function
  emax_fn <- function(Cp, EC50, Emax, hill = 1.2) {
    Cp_h <- max(Cp, 0)^hill
    Emax * Cp_h / (EC50^hill + Cp_h)
  }

  # RTX peak: give at day 0 and day 14 (2 doses)
  RTX_dose_each <- dose_rtx  # mg per dose

  for (i in 2:n) {
    t <- t_seq[i]
    v <- vars

    # ── Drug dosing (daily, converted to per-day) ──────────────────────────
    # HCQ: 400mg/day oral
    HCQ_input  <- ifelse(dose_hcq > 0,  dose_hcq * 0.79, 0)
    PRED_input <- ifelse(dose_pred > 0, dose_pred * 0.92, 0)
    MPA_input  <- ifelse(dose_mmf > 0,  dose_mmf * 0.94, 0)
    BOS_input  <- ifelse(dose_bos > 0,  dose_bos * 0.50, 0)

    # RTX IV bolus at day 0 and 14 (single day spike)
    RTX_input  <- ifelse(dose_rtx > 0 & (t == 1 | t == 15),
                         RTX_dose_each, 0)

    # ── PK (daily Cp approximation) ────────────────────────────────────────
    v$HCQ_cp  <- v$HCQ_cp  + HCQ_input / (HCQ_Vd * wt) -
                 (HCQ_CL / (HCQ_Vd * wt)) * dt * v$HCQ_cp
    v$MPA_cp  <- v$MPA_cp  + MPA_input / (MPA_Vd * wt) -
                 (MPA_CL / (MPA_Vd * wt)) * dt * v$MPA_cp
    v$PRED_cp <- v$PRED_cp + PRED_input / (PRED_Vd * wt) -
                 (PRED_CL / (PRED_Vd * wt)) * dt * v$PRED_cp
    v$RTX_cp  <- v$RTX_cp  + RTX_input / RTX_Vd -
                 (RTX_CL / RTX_Vd) * dt * v$RTX_cp
    v$BOS_cp  <- v$BOS_cp  + BOS_input / (BOS_Vd * wt) -
                 (BOS_CL / (BOS_Vd * wt)) * dt * v$BOS_cp

    # WBC concentration (HCQ concentrates ~2000x)
    HCQ_wbc <- v$HCQ_cp * 200

    # Drug effects
    hcq_e  <- emax_fn(HCQ_wbc,  0.8,   0.75)
    mpa_e  <- emax_fn(v$MPA_cp, 1.5,   0.85)
    pred_e <- emax_fn(v$PRED_cp, 0.05, 0.90)
    rtx_e  <- emax_fn(v$RTX_cp, 5.0,   0.99, 2.0)
    bos_e  <- emax_fn(v$BOS_cp, 0.15,  0.80)

    Th_inflam <- (v$Th1 + v$Th17) / (150 + 60)
    Ab_drive  <- v$AbU1 / Ab0

    # ── T Cells ────────────────────────────────────────────────────────────
    Th1_prod  <- 0.15 * 150 * (1 + 0.5 * v$IFNa / 20) *
                 (1 - pred_e * 0.7) * (1 - hcq_e * 0.4)
    v$Th1 <- max(0, v$Th1  + (Th1_prod  - 0.10 * v$Th1)  * dt)

    Th17_prod <- 0.12 * 60  * (1 + 0.3 * v$IL6 / 15) *
                 (1 - pred_e * 0.8) * (1 - mpa_e * 0.5)
    v$Th17 <- max(0, v$Th17 + (Th17_prod - 0.10 * v$Th17) * dt)

    Treg_prod <- 0.08 * 80  * (1 + hcq_e * 0.2)
    v$Treg <- max(0, v$Treg + (Treg_prod - 0.05 * v$Treg) * dt)

    # ── B Cells & Autoantibody ─────────────────────────────────────────────
    Bnai_loss <- 0.20 * v$Bnai * (1 - rtx_e) * (1 - mpa_e * 0.8)
    v$Bnai <- max(0, v$Bnai + (0.20 * 200 * (1 - rtx_e) -
                                Bnai_loss - 0.08 * v$Bnai) * dt)
    v$PC   <- max(0, v$PC + (0.10 * v$Bnai * (1 - rtx_e) -
                              0.08 * v$PC) * dt)

    Ab_prod <- 0.005 * v$PC * (1 - hcq_e * 0.4)
    v$AbU1 <- max(0, v$AbU1 + (Ab_prod - 0.0693 * v$AbU1) * dt)

    # ── Cytokines ──────────────────────────────────────────────────────────
    v$TNF  <- max(0, v$TNF + (2.0 * Th_inflam * (1 - pred_e * 0.85) -
                               0.50 * v$TNF) * dt)
    v$IL6  <- max(0, v$IL6 + (3.0 * Th_inflam * (1 - pred_e * 0.80) -
                               0.50 * v$IL6) * dt)
    v$IL17 <- max(0, v$IL17 + (1.5 * v$Th17 / 60 * (1 - pred_e * 0.75) -
                                0.50 * v$IL17) * dt)
    v$IFNg <- max(0, v$IFNg + (1.0 * v$Th1 / 150 * (1 - pred_e * 0.70) -
                                0.50 * v$IFNg) * dt)
    v$TGFb <- max(0, v$TGFb + (0.5 * (1 + 0.3 * v$IL6 / 15) *
                                 (1 - pred_e * 0.50) - 0.50 * v$TGFb) * dt)
    v$IFNa <- max(0, v$IFNa + (0.15 * Ab_drive * (1 - hcq_e * 0.70) -
                                0.40 * v$IFNa) * dt)

    # ── Vascular ───────────────────────────────────────────────────────────
    v$ET1 <- max(0, v$ET1 + (0.05 * (1 + 0.3 * v$IL6 / 15) *
                               (1 - bos_e) - 0.30 * v$ET1) * dt)
    v$PVR <- max(1.5, v$PVR + (0.08 * v$ET1 * (1 - bos_e * 0.6) -
                                0.05 * (v$PVR - 1.5)) * dt)

    # ── Lung ───────────────────────────────────────────────────────────────
    v$Collagen <- min(1, max(0, v$Collagen + (0.002 * v$TGFb / 5 *
                      (1 - pred_e * 0.4) * (1 - mpa_e * 0.5) -
                      0.001 * v$Collagen) * dt))
    v$FVC  <- max(20, v$FVC + (-0.0003 * (1 + 5 * v$Collagen) * v$FVC +
                                0.01 * pred_e * (FVC0 - v$FVC)) * dt)
    v$DLCO <- max(20, v$DLCO + (-0.00024 * (1 + 3 * v$PVR / PVR0) * v$DLCO +
                                  0.008 * pred_e * (65 - v$DLCO)) * dt)

    # ── Muscle / CK ────────────────────────────────────────────────────────
    v$CK  <- max(50, v$CK + (5.0 * v$IFNg / 8 * (1 - pred_e * 0.70) *
                               (1 - mpa_e * 0.50) - 0.15 * v$CK) * dt)
    v$MMT <- min(80, v$MMT + (0.005 * (80 - v$MMT) * pred_e -
                               0.002 * v$CK / CK0 * (80 - 70)) * dt)

    # ── Joints ─────────────────────────────────────────────────────────────
    v$SJC <- max(0, v$SJC + (0.03 * Th_inflam * (1 - pred_e * 0.75) *
                               (1 - hcq_e * 0.50) - 0.05 * v$SJC) * dt)

    # ── Complement ─────────────────────────────────────────────────────────
    v$C3 <- max(10, v$C3 + (0.20 * 80 - 0.002 * Ab_drive * v$C3 -
                              0.10 * v$C3) * dt)
    v$C4 <- max(5,  v$C4 + (0.10 * 15 - 0.003 * Ab_drive * v$C4 -
                              0.12 * v$C4) * dt)

    vars <- v
    for (nm in names(v)) results[i, nm] <- v[[nm]]
  }

  # Derived endpoints
  results <- results %>%
    mutate(
      RVSP    = PVR * 5 + 8,
      SixMWT  = pmax(50, 450 - 80 * (PVR - 2.0)),
      DAS28   = 0.56 * sqrt(pmax(0, SJC)) + 0.28 * sqrt(pmax(0, SJC * 1.3)) +
                0.70 * log(pmax(1, 40 * IL6 / 15)) + 0.014 * 50,
      ESR     = 30 * IL6 / 15,
      CRP     = 10 * IL6 / 15,
      WHO_FC  = case_when(RVSP < 35 ~ "I", RVSP < 50 ~ "II",
                          RVSP < 70 ~ "III", TRUE ~ "IV"),
      MCTD_AI = (as.integer(AbU1 > Ab0 * 1.5)) +
                (as.integer(CK > 1000) * 2 + as.integer(CK > 500)) +
                (as.integer(RVSP > 40) * 2 + as.integer(RVSP > 35)) +
                (as.integer(FVC < 70) * 2 + as.integer(FVC < 80)) +
                as.integer(SJC > 8) + as.integer(ESR > 60),
      HCQ_wbc = HCQ_cp * 200
    )

  return(results)
}

## ─── Scenario labels ─────────────────────────────────────────────────────────
scenarios <- list(
  "No Treatment"        = list(0, 0, 0, 0, 0),
  "HCQ Monotherapy"     = list(400, 0, 0, 0, 0),
  "HCQ + Pred 20mg"     = list(400, 20, 0, 0, 0),
  "HCQ + Pred + MMF"    = list(400, 20, 2000, 0, 0),
  "HCQ + Pred + MMF + RTX" = list(400, 20, 2000, 1000, 0),
  "HCQ + Pred + MMF + BOS" = list(400, 20, 2000, 0, 250)
)

## ─── UI ──────────────────────────────────────────────────────────────────────
ui <- dashboardPage(
  skin = "purple",
  header = dashboardHeader(
    title = "MCTD QSP Dashboard",
    titleWidth = 280
  ),

  sidebar = dashboardSidebar(
    width = 280,
    sidebarMenu(
      menuItem("Patient Profile",   tabName = "tab1", icon = icon("user-circle")),
      menuItem("Drug PK",           tabName = "tab2", icon = icon("pills")),
      menuItem("Immunological PD",  tabName = "tab3", icon = icon("dna")),
      menuItem("Organ Function",    tabName = "tab4", icon = icon("lungs")),
      menuItem("Scenario Comparison",tabName= "tab5", icon = icon("chart-bar")),
      menuItem("Biomarker & Risk",  tabName = "tab6", icon = icon("heartbeat"))
    ),
    hr(),
    h5("  Treatment Parameters", style="color:#ccc; margin-left:15px"),

    sliderInput("dose_hcq",  "HCQ (mg/day)",  0, 400,  400, step = 200),
    sliderInput("dose_pred", "Prednisone (mg/day)", 0, 60, 20, step = 5),
    sliderInput("dose_mmf",  "MMF (mg/day)",  0, 3000, 2000, step = 500),
    sliderInput("dose_rtx",  "Rituximab (mg, ×2 doses)", 0, 1000, 0, step = 500),
    sliderInput("dose_bos",  "Bosentan (mg/day, ERA)", 0, 250, 0, step = 125),
    hr(),
    h5("  Patient Parameters", style="color:#ccc; margin-left:15px"),
    sliderInput("wt",     "Body weight (kg)",     40, 100, 60),
    sliderInput("Ab0",    "Initial Anti-U1-RNP (AU/mL)", 50, 1000, 300, step = 50),
    sliderInput("CK0",    "Initial CK (U/L)",    100, 5000, 800, step = 100),
    sliderInput("FVC0",   "Initial FVC %pred",    40, 95, 75),
    sliderInput("RVSP0",  "Initial RVSP (mmHg)",  20, 80, 45),
    sliderInput("duration", "Simulation duration (days)", 90, 730, 365, step = 90)
  ),

  body = dashboardBody(
    tags$head(tags$style(HTML("
      .content-wrapper { background: #f4f6f9; }
      .box { border-top: 3px solid #6C3483; }
      .value-box { border-radius: 8px; }
    "))),
    tabItems(

      ## ── Tab 1: Patient Profile ─────────────────────────────────────────────
      tabItem("tab1",
        h2("MCTD Patient Profile & Disease Overview"),
        fluidRow(
          valueBoxOutput("vb_AbU1", width = 3),
          valueBoxOutput("vb_CK",   width = 3),
          valueBoxOutput("vb_FVC",  width = 3),
          valueBoxOutput("vb_RVSP", width = 3)
        ),
        fluidRow(
          box(title = "Disease Overview — What is MCTD?",
              width = 6, solidHeader = TRUE, status = "purple",
              HTML("
              <h4>Mixed Connective Tissue Disease (MCTD)</h4>
              <p><b>Hallmark:</b> High-titer anti-U1-RNP autoantibodies</p>
              <p><b>Overlap features:</b></p>
              <ul>
                <li><b>SLE:</b> Arthritis, serositis, lymphopenia, low complement</li>
                <li><b>Systemic Sclerosis:</b> Raynaud's, esophageal dysmotility, puffy hands, ILD</li>
                <li><b>Polymyositis:</b> Proximal muscle weakness, elevated CK, myositis</li>
              </ul>
              <p><b>Key pathophysiology:</b></p>
              <ul>
                <li>U1-snRNP complex exposure → TLR7/8 activation → IFN-α surge</li>
                <li>T cell (Th1/Th17) and B cell activation → anti-U1-RNP IgG</li>
                <li>Vascular endothelial injury → Raynaud's, PAH</li>
                <li>TGF-β driven fibrosis → ILD</li>
              </ul>
              <p><b>Prevalence:</b> ~2-4 per 100,000; F:M ratio ≈ 9:1</p>
              <p><b>Treatment:</b> HCQ (backbone) + corticosteroids + MMF/AZA ±
              rituximab; ERA/PDE5i for PAH</p>
              ")),
          box(title = "Baseline Clinical Summary", width = 6,
              solidHeader = TRUE, status = "purple",
              tableOutput("tbl_baseline"))
        ),
        fluidRow(
          box(title = "Diagnostic Criteria (Alarcon-Segovia & Kasukawa)",
              width = 12, solidHeader = TRUE,
              DTOutput("dt_criteria"))
        )
      ),

      ## ── Tab 2: Drug PK ────────────────────────────────────────────────────
      tabItem("tab2",
        h2("Drug PK — Concentration-Time Profiles"),
        fluidRow(
          box(title = "HCQ & Prednisolone Plasma Concentrations",
              width = 6, plotlyOutput("pk_hcq_pred")),
          box(title = "MPA & Bosentan Plasma Concentrations",
              width = 6, plotlyOutput("pk_mpa_bos"))
        ),
        fluidRow(
          box(title = "Rituximab — Target-Mediated PK & B Cell Depletion",
              width = 6, plotlyOutput("pk_rtx")),
          box(title = "HCQ WBC Concentration (Key Pharmacological Compartment)",
              width = 6, plotlyOutput("pk_hcq_wbc"))
        ),
        fluidRow(
          box(title = "PK Parameters Summary", width = 12,
              DTOutput("dt_pk_params"))
        )
      ),

      ## ── Tab 3: Immunological PD ──────────────────────────────────────────
      tabItem("tab3",
        h2("Immunological PD — Cytokines, Autoantibodies & Immune Cells"),
        fluidRow(
          box(title = "Anti-U1-RNP Antibody Titer (MCTD Biomarker)",
              width = 6, plotlyOutput("pd_antibody")),
          box(title = "Cytokine Dynamics",
              width = 6, plotlyOutput("pd_cytokines"))
        ),
        fluidRow(
          box(title = "T Cell Subsets (Th1 / Th17 / Treg)",
              width = 6, plotlyOutput("pd_tcell")),
          box(title = "B Cells & Complement",
              width = 6, plotlyOutput("pd_bcell"))
        ),
        fluidRow(
          box(title = "IFN-α Signature (Type I Interferon)",
              width = 6, plotlyOutput("pd_ifna")),
          box(title = "Immune Summary Table at Key Timepoints",
              width = 6, DTOutput("dt_immune"))
        )
      ),

      ## ── Tab 4: Organ Function ────────────────────────────────────────────
      tabItem("tab4",
        h2("Organ Function — Lung (ILD/PAH), Muscle, Joints"),
        fluidRow(
          box(title = "Lung Function — FVC & DLCO (ILD)",
              width = 6, plotlyOutput("org_lung")),
          box(title = "Pulmonary Vascular — RVSP & PVR (PAH)",
              width = 6, plotlyOutput("org_pah"))
        ),
        fluidRow(
          box(title = "6-Minute Walk Test & WHO Functional Class",
              width = 6, plotlyOutput("org_6mwt")),
          box(title = "Muscle Function — CK & MMT-8 Score",
              width = 6, plotlyOutput("org_muscle"))
        ),
        fluidRow(
          box(title = "Joint Inflammation — SJC & DAS28",
              width = 6, plotlyOutput("org_joint")),
          box(title = "Vascular — Endothelin-1 & Raynaud's Severity",
              width = 6, plotlyOutput("org_vascular"))
        )
      ),

      ## ── Tab 5: Scenario Comparison ─────────────────────────────────────
      tabItem("tab5",
        h2("Treatment Scenario Comparison"),
        fluidRow(
          box(title = "Select Scenarios to Compare", width = 12,
              checkboxGroupInput("sc_select",
                label = "Choose scenarios:",
                choices = names(scenarios),
                selected = names(scenarios),
                inline = TRUE))
        ),
        fluidRow(
          box(title = "MCTD Activity Index — All Scenarios",
              width = 6, plotlyOutput("sc_mctd_ai")),
          box(title = "Anti-U1-RNP Titer — All Scenarios",
              width = 6, plotlyOutput("sc_antibody"))
        ),
        fluidRow(
          box(title = "FVC % Predicted — Lung ILD",
              width = 6, plotlyOutput("sc_fvc")),
          box(title = "RVSP — PAH Endpoint",
              width = 6, plotlyOutput("sc_rvsp"))
        ),
        fluidRow(
          box(title = "Endpoint Summary Table (Day 90 / 180 / 365)",
              width = 12, DTOutput("dt_sc_summary"))
        )
      ),

      ## ── Tab 6: Biomarker & Risk ─────────────────────────────────────────
      tabItem("tab6",
        h2("Biomarker Heatmap & Risk Stratification"),
        fluidRow(
          box(title = "Risk Classification at End of Simulation",
              width = 4, status = "danger", solidHeader = TRUE,
              htmlOutput("risk_class")),
          box(title = "Key Biomarker Trend (Baseline vs. End)",
              width = 8, plotlyOutput("bm_comparison"))
        ),
        fluidRow(
          box(title = "Biomarker Radar Chart",
              width = 6, plotlyOutput("bm_radar")),
          box(title = "PAH Risk Assessment (ESC/ERS 2022 Criteria)",
              width = 6, htmlOutput("pah_risk"))
        ),
        fluidRow(
          box(title = "Complete Biomarker Table", width = 12,
              DTOutput("dt_biomarkers"))
        )
      )
    )
  )
)

## ─── Server ──────────────────────────────────────────────────────────────────
server <- function(input, output, session) {

  ## Reactive simulation (single scenario from sidebar)
  sim_data <- reactive({
    solve_mctd(
      dose_hcq  = input$dose_hcq,
      dose_pred = input$dose_pred,
      dose_mmf  = input$dose_mmf,
      dose_rtx  = input$dose_rtx,
      dose_bos  = input$dose_bos,
      duration  = input$duration,
      wt        = input$wt,
      Ab0       = input$Ab0,
      CK0       = input$CK0,
      FVC0      = input$FVC0,
      RVSP0     = input$RVSP0,
      PVR0      = input$RVSP0 / 5 - 1.6
    )
  })

  ## Reactive multi-scenario data
  all_scenario_data <- reactive({
    req(input$sc_select)
    sc_list <- lapply(input$sc_select, function(nm) {
      args <- scenarios[[nm]]
      df   <- solve_mctd(
        dose_hcq  = args[[1]], dose_pred = args[[2]],
        dose_mmf  = args[[3]], dose_rtx  = args[[4]],
        dose_bos  = args[[5]],
        duration  = input$duration, wt = input$wt,
        Ab0 = input$Ab0, CK0 = input$CK0,
        FVC0 = input$FVC0, RVSP0 = input$RVSP0
      )
      df$Scenario <- nm
      df
    })
    bind_rows(sc_list)
  })

  ## ── Value Boxes ──────────────────────────────────────────────────────────
  output$vb_AbU1 <- renderValueBox({
    d <- tail(sim_data(), 1)
    valueBox(round(d$AbU1, 0), "Anti-U1-RNP (AU/mL)",
             icon = icon("vials"),
             color = ifelse(d$AbU1 > 500, "red", ifelse(d$AbU1 > 200, "yellow", "green")))
  })
  output$vb_CK <- renderValueBox({
    d <- tail(sim_data(), 1)
    valueBox(round(d$CK, 0), "CK (U/L)",
             icon = icon("dumbbell"),
             color = ifelse(d$CK > 1000, "red", ifelse(d$CK > 400, "yellow", "green")))
  })
  output$vb_FVC <- renderValueBox({
    d <- tail(sim_data(), 1)
    valueBox(paste0(round(d$FVC, 1), "%"), "FVC %pred",
             icon = icon("lungs"),
             color = ifelse(d$FVC < 60, "red", ifelse(d$FVC < 75, "yellow", "green")))
  })
  output$vb_RVSP <- renderValueBox({
    d <- tail(sim_data(), 1)
    valueBox(round(d$RVSP, 0), "RVSP (mmHg)",
             icon = icon("heartbeat"),
             color = ifelse(d$RVSP > 60, "red", ifelse(d$RVSP > 40, "yellow", "green")))
  })

  ## ── Tab 1: Baseline Table ────────────────────────────────────────────────
  output$tbl_baseline <- renderTable({
    d0 <- sim_data()[1, ]
    data.frame(
      Parameter = c("Anti-U1-RNP (AU/mL)", "CK (U/L)", "FVC (%pred)",
                    "RVSP (mmHg)", "ESR (mm/hr)", "CRP (mg/L)",
                    "SJC", "DAS28", "MCTD Activity Index"),
      Value = c(round(d0$AbU1, 0), round(d0$CK, 0), round(d0$FVC, 1),
                round(d0$RVSP, 0), round(d0$ESR, 0), round(d0$CRP, 1),
                round(d0$SJC, 1), round(d0$DAS28, 2), round(d0$MCTD_AI, 0)),
      Normal_Range = c("<25", "<200", ">80", "<35", "<20", "<5",
                       "0", "<2.6", "0")
    )
  })

  output$dt_criteria <- renderDT({
    data.frame(
      Criterion = c("Clinical 1: Raynaud's Phenomenon",
                    "Clinical 2: Synovitis / Arthritis",
                    "Clinical 3: Myositis (proximal weakness)",
                    "Clinical 4: Swollen hands / Sclerodactyly",
                    "Clinical 5: Acrosclerosis / SSc overlap",
                    "Serological: Anti-U1-RNP ≥ 1:320",
                    "Confirmatory: High anti-U1-RNP (>200 AU/mL)"),
      Kasukawa = c("Required", "≥1 of 2 groups", "≥1 of 2 groups",
                   "Required", "—", "Required", "—"),
      Alarcon_Segovia = c("Required", "✓", "✓", "✓", "✓",
                          "Required", "Required"),
      Notes = c("Digital vasospasm on cold/stress",
                "Non-erosive, RA-like",
                "Proximal weakness; ↑CK; ↑LDH",
                "Puffy fingers, sausage digits",
                "Lung/esophageal involvement",
                "Speckled ANA, anti-snRNP",
                "Key distinguishing feature")
    )
  }, options = list(pageLength = 7, scrollX = TRUE))

  ## ── Tab 2: PK Plots ──────────────────────────────────────────────────────
  output$pk_hcq_pred <- renderPlotly({
    d <- sim_data()
    p <- ggplot(d, aes(x = time)) +
      geom_line(aes(y = HCQ_cp, color = "HCQ (μg/mL)"), size = 1.2) +
      geom_line(aes(y = PRED_cp * 100, color = "Prednisolone (×100, mg/L)"),
                size = 1.2, linetype = "dashed") +
      scale_color_manual(values = c("HCQ (μg/mL)" = "#3498DB",
                                    "Prednisolone (×100, mg/L)" = "#E74C3C")) +
      labs(x = "Time (days)", y = "Concentration", color = NULL) +
      theme_minimal()
    ggplotly(p)
  })

  output$pk_mpa_bos <- renderPlotly({
    d <- sim_data()
    p <- ggplot(d, aes(x = time)) +
      geom_line(aes(y = MPA_cp, color = "MPA (mg/L)"), size = 1.2) +
      geom_line(aes(y = BOS_cp * 10, color = "Bosentan (×10, mg/L)"),
                size = 1.2, linetype = "dashed") +
      scale_color_manual(values = c("MPA (mg/L)" = "#2ECC71",
                                    "Bosentan (×10, mg/L)" = "#E67E22")) +
      labs(x = "Time (days)", y = "Concentration", color = NULL) +
      theme_minimal()
    ggplotly(p)
  })

  output$pk_rtx <- renderPlotly({
    d <- sim_data()
    p <- ggplot(d, aes(x = time)) +
      geom_line(aes(y = RTX_cp, color = "Rituximab (μg/mL)"),
                size = 1.2, color = "#9B59B6") +
      geom_line(aes(y = Bnai / 200 * 10, color = "B cells (norm.)"),
                size = 1.2, color = "#E74C3C", linetype = "dashed") +
      labs(x = "Time (days)", y = "RTX conc. / B cells (normalized)",
           title = "Rituximab PK & B Cell Depletion") +
      theme_minimal()
    ggplotly(p)
  })

  output$pk_hcq_wbc <- renderPlotly({
    d <- sim_data()
    p <- ggplot(d, aes(x = time, y = HCQ_wbc)) +
      geom_line(size = 1.2, color = "#3498DB") +
      geom_hline(yintercept = 0.8, color = "red", linetype = "dashed") +
      annotate("text", x = 10, y = 0.9,
               label = "IC50 = 0.8 μg/mL (TLR7/8 inhibition)") +
      labs(x = "Time (days)", y = "HCQ WBC Concentration (μg/mL)") +
      theme_minimal()
    ggplotly(p)
  })

  output$dt_pk_params <- renderDT({
    data.frame(
      Drug = c("HCQ", "MPA (from MMF)", "Prednisolone", "Rituximab", "Bosentan"),
      Dose = c("400 mg/day", "1-3 g/day", "5-60 mg/day", "1000 mg ×2", "62.5-125 mg BID"),
      F_pct = c("79%", "94%", "92%", "100% (IV)", "50%"),
      t_half = c("40-56 days", "17h", "2-4h", "22 days", "5h"),
      Vd = c("800 L/kg", "3.6 L/kg", "0.55 L/kg", "3.5 L", "0.18 L/kg"),
      PK_target = c("WBC 200:1 plasma", "AUC 30-60 mg·h/L", "GR occupancy >50%", "CD20 B cell saturation", "ETA receptor >90%"),
      Key_PD = c("TLR7/8 inhibition, ↓IFN-α", "IMPDH inhibition, ↓B/T cells", "NF-κB transrepression", "B cell depletion >99%", "ET-1 blockade, ↓PVR")
    )
  }, options = list(scrollX = TRUE))

  ## ── Tab 3: PD Plots ──────────────────────────────────────────────────────
  output$pd_antibody <- renderPlotly({
    d <- sim_data()
    p <- ggplot(d, aes(x = time, y = AbU1)) +
      geom_line(size = 1.2, color = "#E74C3C") +
      geom_hline(yintercept = c(100, 300, 500), linetype = "dashed",
                 color = c("green", "orange", "red"), alpha = 0.6) +
      labs(x = "Time (days)", y = "Anti-U1-RNP (AU/mL)",
           title = "MCTD Hallmark Autoantibody") +
      theme_minimal()
    ggplotly(p)
  })

  output$pd_cytokines <- renderPlotly({
    d <- sim_data() %>%
      select(time, TNF, IL6, IL17, IFNg, TGFb) %>%
      pivot_longer(-time, names_to = "Cytokine", values_to = "Level")
    p <- ggplot(d, aes(x = time, y = Level, color = Cytokine)) +
      geom_line(size = 1.1) +
      scale_color_brewer(palette = "Set1") +
      labs(x = "Time (days)", y = "Cytokine (relative pg/mL)") +
      theme_minimal()
    ggplotly(p)
  })

  output$pd_tcell <- renderPlotly({
    d <- sim_data() %>%
      select(time, Th1, Th17, Treg) %>%
      pivot_longer(-time, names_to = "Cell", values_to = "Count")
    p <- ggplot(d, aes(x = time, y = Count, color = Cell)) +
      geom_line(size = 1.2) +
      scale_color_manual(values = c(Th1="#E74C3C", Th17="#9B59B6", Treg="#2ECC71")) +
      labs(x = "Time (days)", y = "Cells/μL") +
      theme_minimal()
    ggplotly(p)
  })

  output$pd_bcell <- renderPlotly({
    d <- sim_data() %>%
      select(time, Bnai, PC, C3, C4) %>%
      pivot_longer(-time, names_to = "Marker", values_to = "Value")
    p <- ggplot(d, aes(x = time, y = Value, color = Marker)) +
      geom_line(size = 1.1) +
      scale_color_manual(values = c(Bnai="#3498DB", PC="#E67E22",
                                    C3="#27AE60", C4="#8E44AD")) +
      labs(x = "Time (days)", y = "Cells/μL or mg/dL") +
      theme_minimal()
    ggplotly(p)
  })

  output$pd_ifna <- renderPlotly({
    d <- sim_data()
    p <- ggplot(d, aes(x = time, y = IFNa)) +
      geom_line(size = 1.2, color = "#8E44AD") +
      geom_hline(yintercept = 10, linetype = "dashed", color = "red") +
      annotate("text", x = 10, y = 11, label = "Normal upper limit") +
      labs(x = "Time (days)", y = "IFN-α (relative units)",
           title = "Type I Interferon Signature") +
      theme_minimal()
    ggplotly(p)
  })

  output$dt_immune <- renderDT({
    d <- sim_data() %>%
      filter(time %in% c(0, 30, 90, 180, 365)) %>%
      select(time, AbU1, Th1, Th17, Treg, Bnai, TNF, IL6, IL17, IFNa, C3, C4) %>%
      mutate(across(where(is.numeric), ~round(.x, 1)))
    datatable(d, options = list(scrollX = TRUE, pageLength = 5))
  })

  ## ── Tab 4: Organ Function ────────────────────────────────────────────────
  output$org_lung <- renderPlotly({
    d <- sim_data()
    p <- ggplot(d, aes(x = time)) +
      geom_line(aes(y = FVC, color = "FVC %pred"), size = 1.2) +
      geom_line(aes(y = DLCO, color = "DLCO %pred"), size = 1.2, linetype = "dashed") +
      geom_hline(yintercept = c(70, 80), linetype = "dotted",
                 color = c("red", "orange"), alpha = 0.7) +
      scale_color_manual(values = c("FVC %pred" = "#2E86C1",
                                    "DLCO %pred" = "#E67E22")) +
      labs(x = "Time (days)", y = "% Predicted", color = NULL) +
      theme_minimal()
    ggplotly(p)
  })

  output$org_pah <- renderPlotly({
    d <- sim_data()
    p <- ggplot(d, aes(x = time)) +
      geom_line(aes(y = RVSP, color = "RVSP (mmHg)"), size = 1.2) +
      geom_line(aes(y = PVR * 10, color = "PVR ×10 (WU)"), size = 1.2,
                linetype = "dashed") +
      geom_hline(yintercept = c(35, 50, 70), linetype = "dashed",
                 color = c("yellow", "orange", "red"), alpha = 0.5) +
      scale_color_manual(values = c("RVSP (mmHg)" = "#C0392B",
                                    "PVR ×10 (WU)" = "#8E44AD")) +
      labs(x = "Time (days)", y = "RVSP (mmHg) / PVR×10 (WU)", color = NULL) +
      theme_minimal()
    ggplotly(p)
  })

  output$org_6mwt <- renderPlotly({
    d <- sim_data()
    p <- ggplot(d, aes(x = time, y = SixMWT)) +
      geom_line(size = 1.2, color = "#16A085") +
      geom_hline(yintercept = c(300, 400), linetype = "dashed",
                 color = c("red", "orange"), alpha = 0.7) +
      annotate("text", x = 20, y = 300, label = "WHO FC III threshold", size = 3) +
      labs(x = "Time (days)", y = "6-Minute Walk Distance (m)") +
      theme_minimal()
    ggplotly(p)
  })

  output$org_muscle <- renderPlotly({
    d <- sim_data()
    p <- ggplot(d, aes(x = time)) +
      geom_line(aes(y = CK / 50, color = "CK/50 (U/L)"), size = 1.2) +
      geom_line(aes(y = MMT, color = "MMT-8 Score"), size = 1.2, linetype = "dashed") +
      scale_color_manual(values = c("CK/50 (U/L)" = "#E74C3C",
                                    "MMT-8 Score" = "#27AE60")) +
      labs(x = "Time (days)", y = "CK/50 or MMT-8 Score", color = NULL) +
      theme_minimal()
    ggplotly(p)
  })

  output$org_joint <- renderPlotly({
    d <- sim_data()
    p <- ggplot(d, aes(x = time)) +
      geom_line(aes(y = SJC, color = "Swollen Joint Count"), size = 1.2) +
      geom_line(aes(y = DAS28, color = "DAS28"), size = 1.2, linetype = "dashed") +
      geom_hline(yintercept = c(2.6, 3.2), linetype = "dashed",
                 color = c("green", "orange"), alpha = 0.5) +
      scale_color_manual(values = c("Swollen Joint Count" = "#E74C3C",
                                    "DAS28" = "#9B59B6")) +
      labs(x = "Time (days)", y = "SJC / DAS28", color = NULL) +
      theme_minimal()
    ggplotly(p)
  })

  output$org_vascular <- renderPlotly({
    d <- sim_data()
    p <- ggplot(d, aes(x = time, y = ET1)) +
      geom_line(size = 1.2, color = "#E74C3C") +
      geom_hline(yintercept = 2, linetype = "dashed", color = "orange") +
      labs(x = "Time (days)", y = "Endothelin-1 (pg/mL)",
           title = "Endothelin-1 — Vascular Biomarker") +
      theme_minimal()
    ggplotly(p)
  })

  ## ── Tab 5: Scenario Comparison ──────────────────────────────────────────
  sc_colors <- c(
    "No Treatment"            = "#E74C3C",
    "HCQ Monotherapy"         = "#3498DB",
    "HCQ + Pred 20mg"         = "#F39C12",
    "HCQ + Pred + MMF"        = "#2ECC71",
    "HCQ + Pred + MMF + RTX"  = "#9B59B6",
    "HCQ + Pred + MMF + BOS"  = "#1ABC9C"
  )

  make_sc_plot <- function(var, ylab, threshold = NULL, thresh_color = "red") {
    d <- all_scenario_data()
    p <- ggplot(d, aes(x = time, y = .data[[var]], color = Scenario)) +
      geom_line(size = 1.1) +
      scale_color_manual(values = sc_colors) +
      labs(x = "Time (days)", y = ylab, color = NULL) +
      theme_minimal() +
      theme(legend.position = "bottom")
    if (!is.null(threshold)) {
      p <- p + geom_hline(yintercept = threshold, linetype = "dashed",
                          color = thresh_color, alpha = 0.6)
    }
    ggplotly(p) %>% layout(legend = list(orientation = "h", y = -0.2))
  }

  output$sc_mctd_ai  <- renderPlotly(make_sc_plot("MCTD_AI", "MCTD Activity Index", 3))
  output$sc_antibody <- renderPlotly(make_sc_plot("AbU1",    "Anti-U1-RNP (AU/mL)", 200))
  output$sc_fvc      <- renderPlotly(make_sc_plot("FVC",     "FVC % Predicted", 70, "red"))
  output$sc_rvsp     <- renderPlotly(make_sc_plot("RVSP",    "RVSP (mmHg)", 40, "orange"))

  output$dt_sc_summary <- renderDT({
    d <- all_scenario_data() %>%
      filter(time %in% c(0, 90, 180, 365)) %>%
      select(Scenario, time, AbU1, CK, FVC, RVSP, DAS28, MCTD_AI, SixMWT) %>%
      mutate(across(where(is.numeric), ~round(.x, 1)))
    datatable(d, filter = "top",
              options = list(pageLength = 10, scrollX = TRUE))
  })

  ## ── Tab 6: Biomarker & Risk ──────────────────────────────────────────────
  output$risk_class <- renderUI({
    d <- tail(sim_data(), 1)
    ai <- d$MCTD_AI
    color <- ifelse(ai >= 6, "red", ifelse(ai >= 3, "orange", "green"))
    label <- ifelse(ai >= 6, "High Activity",
                    ifelse(ai >= 3, "Moderate Activity", "Low/Remission"))
    tags$div(
      tags$h3(style = paste0("color:", color), label),
      tags$p(paste("MCTD Activity Index:", round(ai, 0))),
      tags$hr(),
      tags$p(paste("Anti-U1-RNP:", round(d$AbU1, 0), "AU/mL")),
      tags$p(paste("CK:", round(d$CK, 0), "U/L")),
      tags$p(paste("FVC:", round(d$FVC, 1), "% pred")),
      tags$p(paste("RVSP:", round(d$RVSP, 0), "mmHg")),
      tags$p(paste("WHO Functional Class:", d$WHO_FC)),
      tags$p(paste("DAS28:", round(d$DAS28, 2))),
      tags$p(paste("ESR:", round(d$ESR, 0), "mm/hr"))
    )
  })

  output$bm_comparison <- renderPlotly({
    d     <- sim_data()
    d_bm  <- data.frame(
      Biomarker = c("Anti-U1-RNP/10", "CK/100", "FVC%", "RVSP",
                    "DAS28×10", "ESR", "CRP×5", "C3%"),
      Baseline  = c(d$AbU1[1]/10, d$CK[1]/100, d$FVC[1], d$RVSP[1],
                    d$DAS28[1]*10, d$ESR[1], d$CRP[1]*5, d$C3[1]),
      End       = c(tail(d$AbU1, 1)/10, tail(d$CK, 1)/100, tail(d$FVC, 1),
                    tail(d$RVSP, 1), tail(d$DAS28, 1)*10, tail(d$ESR, 1),
                    tail(d$CRP, 1)*5, tail(d$C3, 1))
    ) %>% pivot_longer(-Biomarker, names_to = "Timepoint", values_to = "Value")

    p <- ggplot(d_bm, aes(x = Biomarker, y = Value, fill = Timepoint)) +
      geom_col(position = "dodge") +
      scale_fill_manual(values = c(Baseline = "#E74C3C", End = "#2ECC71")) +
      coord_flip() +
      labs(y = "Normalized Value", fill = NULL) +
      theme_minimal()
    ggplotly(p)
  })

  output$bm_radar <- renderPlotly({
    d <- tail(sim_data(), 1)
    cats <- c("Anti-U1-RNP", "CK", "FVC", "RVSP", "DAS28", "MCTD-AI")
    vals <- c(
      min(d$AbU1 / 500, 1),
      min(d$CK / 2000, 1),
      1 - min(d$FVC / 90, 1),
      min(d$RVSP / 70, 1),
      min(d$DAS28 / 6, 1),
      min(d$MCTD_AI / 8, 1)
    )
    plot_ly(type = "scatterpolar", fill = "toself") %>%
      add_trace(r = c(vals, vals[1]),
                theta = c(cats, cats[1]),
                name = "Current State",
                fillcolor = "rgba(231,76,60,0.3)",
                line = list(color = "#E74C3C")) %>%
      layout(polar = list(radialaxis = list(visible = TRUE, range = c(0, 1))),
             title = "Disease Activity Radar")
  })

  output$pah_risk <- renderUI({
    d <- tail(sim_data(), 1)
    risk <- ifelse(d$RVSP > 60, "HIGH", ifelse(d$RVSP > 45, "INTERMEDIATE", "LOW"))
    color <- ifelse(risk == "HIGH", "red", ifelse(risk == "INTERMEDIATE", "orange", "green"))
    tags$div(
      tags$h3(style = paste0("color:", color), paste("PAH Risk:", risk)),
      tags$table(class = "table table-striped",
        tags$tbody(
          tags$tr(tags$th("Parameter"), tags$th("Value"), tags$th("Risk Cut-off")),
          tags$tr(tags$td("RVSP"), tags$td(round(d$RVSP, 0), " mmHg"),
                  tags$td("Low:<40 / High:>55")),
          tags$tr(tags$td("PVR"), tags$td(round(d$PVR, 1), " WU"),
                  tags$td("Low:<4 / High:>7")),
          tags$tr(tags$td("6MWT"), tags$td(round(d$SixMWT, 0), " m"),
                  tags$td("Low:>440 / High:<380")),
          tags$tr(tags$td("WHO FC"), tags$td(d$WHO_FC),
                  tags$td("Low:I-II / High:III-IV"))
        )
      )
    )
  })

  output$dt_biomarkers <- renderDT({
    d <- sim_data() %>%
      filter(time %in% c(0, 30, 90, 180, 365)) %>%
      transmute(
        "Time (days)" = time,
        "Anti-U1-RNP (AU/mL)" = round(AbU1, 0),
        "Serum CK (U/L)" = round(CK, 0),
        "FVC (%pred)" = round(FVC, 1),
        "DLCO (%pred)" = round(DLCO, 1),
        "RVSP (mmHg)" = round(RVSP, 0),
        "PVR (WU)" = round(PVR, 2),
        "6MWT (m)" = round(SixMWT, 0),
        "WHO FC" = WHO_FC,
        "DAS28" = round(DAS28, 2),
        "SJC" = round(SJC, 1),
        "ESR (mm/hr)" = round(ESR, 0),
        "CRP (mg/L)" = round(CRP, 1),
        "C3 (mg/dL)" = round(C3, 1),
        "C4 (mg/dL)" = round(C4, 1),
        "MCTD AI" = round(MCTD_AI, 0),
        "ET-1 (pg/mL)" = round(ET1, 2)
      )
    datatable(d, options = list(scrollX = TRUE, pageLength = 5))
  })
}

## ─── Launch ──────────────────────────────────────────────────────────────────
shinyApp(ui, server)
