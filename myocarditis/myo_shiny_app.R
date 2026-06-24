## ============================================================
## Myocarditis QSP Shiny Dashboard
## Interactive Simulation of Viral/Autoimmune Myocarditis
## Tabs: Overview | PK | Viral & Immune | PD Biomarkers |
##       Clinical Endpoints | Scenario Comparison | Fibrosis
## ============================================================

library(shiny)
library(shinydashboard)
library(ggplot2)
library(dplyr)
library(tidyr)
library(plotly)
library(DT)

# ============================================================
# Simplified inline ODE solver (Euler method for responsiveness)
# ============================================================
run_myocarditis_sim <- function(params, events, tend = 365, dt = 0.5) {

  p <- params

  # State vector initialization
  state <- list(
    H = 5e9, I = 0, D = 0, V = 0,
    NK = 200, M1 = 500, M2 = 300,
    IFNb = 0, IFNg = 5, TNFa = 5, IL6 = 2, IL1b = 1, TGFb = 3, IL10 = 5,
    Tn4 = 600, Th1 = 50, Th17 = 20, Treg = 30,
    Tn8 = 400, CTL = 80,
    Bn = 200, PC = 30, Ab = 0,
    CFib = 50, MFib = 5, Col = 0,
    Tn = 0.01, BNP = 50, EF = 60,
    IVIG_C = 0, PRED_C = 0, AZA_C = 0, CSA_C = 0, COLC_C = 0,
    MP6_C = 0
  )

  times <- seq(0, tend, by = dt)
  out <- vector("list", length(times))
  out[[1]] <- c(time = 0, unlist(state))

  for (idx in seq_along(times)[-1]) {
    t <- times[idx]

    # Apply dosing events
    for (ev in events) {
      if (abs(t - ev$time) < dt / 2) {
        if (ev$cmt == "V")      state$V      <- state$V      + ev$amt
        if (ev$cmt == "IVIG_C") state$IVIG_C <- state$IVIG_C + ev$amt
        if (ev$cmt == "PRED_C") state$PRED_C <- state$PRED_C + ev$amt
        if (ev$cmt == "AZA_C")  state$AZA_C  <- state$AZA_C  + ev$amt
        if (ev$cmt == "CSA_C")  state$CSA_C  <- state$CSA_C  + ev$amt
        if (ev$cmt == "COLC_C") state$COLC_C <- state$COLC_C + ev$amt
      }
    }

    s <- state  # current state

    # Drug effects
    E_IVIG <- p$Emax_IVIG * s$IVIG_C / (p$EC50_IVIG + s$IVIG_C + 1e-9)
    E_PRED <- p$Emax_PRED * s$PRED_C / (p$EC50_PRED + s$PRED_C + 1e-9)
    E_AZA  <- p$Emax_AZA  * s$MP6_C  / (p$EC50_AZA  + s$MP6_C  + 1e-9)
    E_CSA  <- p$Emax_CSA  * s$CSA_C  / (p$EC50_CSA   + s$CSA_C  + 1e-9)
    E_COLC <- p$Emax_COLC * s$COLC_C / (p$EC50_COLC  + s$COLC_C + 1e-9)
    E_IS   <- min(0.95, E_PRED + E_AZA * (1 - E_PRED) + E_CSA * (1 - E_PRED) * (1 - E_AZA))

    # Viral dynamics
    dV <- p$p_V * s$I - p$c_V * s$V * (1 + p$kIFN * s$IFNb / (500 + s$IFNb))

    # Cardiomyocytes
    IS_sig <- (s$TNFa / (50 + s$TNFa) + s$IL1b / (20 + s$IL1b)) * (1 - E_IS) * (1 - E_IVIG)
    dH <- p$r_H * s$H * (1 - (s$H + s$I) / 5e9) - p$d_H * s$H - p$beta_V * s$V * s$H
    dI <- p$beta_V * s$V * s$H - p$delta_I * s$I - (0.006 * s$NK + 0.008 * s$CTL) * s$I * (1 - E_IS)
    dD <- p$delta_I * s$I + (0.006 * s$NK + 0.008 * s$CTL) * s$I * (1 - E_IS) +
          0.04 * IS_sig * s$H - 0.5 * s$D

    # Innate immunity
    dNK  <- 50 + 0.8 * s$IFNb / (100 + s$IFNb) * s$NK - 0.1 * s$NK
    dM1  <- 200 * (s$IFNg / (10 + s$IFNg)) * (1 - s$M1 / 5000) * (1 - E_IS) - 0.12 * s$M1
    dM2  <- 100 * s$IL10 / (10 + s$IL10) * (1 - s$M2 / 2000) - 0.12 * s$M2

    # Cytokines
    DAMP <- s$D + 0.1 * s$I
    dIFNb <- 80  * s$I / (500 + s$I)   - 4   * s$IFNb
    dIFNg <- 0.5 * (s$Th1 + s$NK) / (200 + s$Th1 + s$NK) - 2 * s$IFNg
    dTNFa <- 0.4 * s$M1 * (1 - E_IS) * (1 - E_PRED) - 2   * s$TNFa
    dIL6  <- 0.6 * s$M1 * (1 - E_IS) * (1 - E_PRED) - 3   * s$IL6
    NLRP3 <- DAMP / (500 + DAMP) * (1 - E_COLC)
    dIL1b <- 0.3 * s$M1 * NLRP3 * (1 - E_IS) - 2.5 * s$IL1b
    dTGFb <- 0.2 * (s$M2 + s$Treg) / (100 + s$M2 + s$Treg) - 1.5 * s$TGFb
    dIL10 <- 0.25 * (s$M2 + s$Treg) / (100 + s$M2 + s$Treg) - 2 * s$IL10

    # Adaptive immunity
    Ag_sig <- s$I / (1e6 + s$I) + s$D / (1e6 + s$D)
    Th1_diff  <- 0.4  * s$Tn4 * s$IFNg / (10 + s$IFNg)  * Ag_sig * (1 - E_IS) * (1 - E_CSA)
    Th17_diff <- 0.25 * s$Tn4 * s$IL6  / (15 + s$IL6)   * Ag_sig * (1 - E_IS)
    Treg_diff <- 0.15 * s$Tn4 * s$TGFb / (10 + s$TGFb)
    dTn4  <- 100 - Th1_diff - Th17_diff - Treg_diff - 0.05 * s$Tn4
    dTh1  <- Th1_diff  - 0.08 * s$Th1
    dTh17 <- Th17_diff - 0.09 * s$Th17
    dTreg <- Treg_diff - 0.06 * s$Treg

    CTL_diff <- 0.35 * s$Tn8 * s$Th1 / (200 + s$Th1) * Ag_sig * (1 - E_IS) * (1 - E_CSA)
    dTn8 <- 80  - CTL_diff - 0.05 * s$Tn8
    dCTL <- CTL_diff - 0.07 * s$CTL

    B_act <- 0.15 * s$Bn * s$Th1 / (200 + s$Th1) * Ag_sig * (1 - E_IVIG) * (1 - E_IS)
    dBn   <- 60  - B_act - 0.04 * s$Bn
    dPC   <- B_act - 0.03 * s$PC
    dAb   <- 20  * s$PC * (1 - E_IVIG) - 0.02 * s$Ab

    # Remodeling
    Fib_act <- 0.3 * s$TGFb / (5 + s$TGFb) + 0.1 * s$TNFa / (20 + s$TNFa)
    dCFib <- 50 - Fib_act * s$CFib - 0.01 * s$CFib
    dMFib <- Fib_act * s$CFib - 0.05 * s$MFib
    dCol  <- 0.4 * s$MFib - 0.008 * s$Col

    # Biomarkers
    DeadRate <- p$delta_I * s$I + 0.04 * IS_sig * s$H
    WallStress <- (1 - s$H / 5e9) * 100
    dTn  <- 0.15 * DeadRate * (1 - 0.3 * E_IVIG) - 0.3 * s$Tn
    dBNP <- 0.002 * WallStress - 1.4 * s$BNP

    EF_actual <- max(20, 60 * (s$H / 5e9) * (1 - 0.3 * s$Col / (100 + s$Col)))
    dEF <- 0.03 * (EF_actual - s$EF)

    # Drug PK (simplified 1-compartment)
    dIVIG_C <- -0.0033 * s$IVIG_C
    dPRED_C <- -0.3    * s$PRED_C  # simplified decay
    dAZA_C  <- -0.2    * s$AZA_C
    dMP6_C  <- 0.5 * s$AZA_C - 0.15 * s$MP6_C
    dCSA_C  <- -0.05   * s$CSA_C
    dCOLC_C <- -0.08   * s$COLC_C

    # Euler update
    state$H     <- max(0, s$H     + dH     * dt)
    state$I     <- max(0, s$I     + dI     * dt)
    state$D     <- max(0, s$D     + dD     * dt)
    state$V     <- max(0, s$V     + dV     * dt)
    state$NK    <- max(0, s$NK    + dNK    * dt)
    state$M1    <- max(0, s$M1    + dM1    * dt)
    state$M2    <- max(0, s$M2    + dM2    * dt)
    state$IFNb  <- max(0, s$IFNb  + dIFNb  * dt)
    state$IFNg  <- max(0, s$IFNg  + dIFNg  * dt)
    state$TNFa  <- max(0, s$TNFa  + dTNFa  * dt)
    state$IL6   <- max(0, s$IL6   + dIL6   * dt)
    state$IL1b  <- max(0, s$IL1b  + dIL1b  * dt)
    state$TGFb  <- max(0, s$TGFb  + dTGFb  * dt)
    state$IL10  <- max(0, s$IL10  + dIL10  * dt)
    state$Tn4   <- max(0, s$Tn4   + dTn4   * dt)
    state$Th1   <- max(0, s$Th1   + dTh1   * dt)
    state$Th17  <- max(0, s$Th17  + dTh17  * dt)
    state$Treg  <- max(0, s$Treg  + dTreg  * dt)
    state$Tn8   <- max(0, s$Tn8   + dTn8   * dt)
    state$CTL   <- max(0, s$CTL   + dCTL   * dt)
    state$Bn    <- max(0, s$Bn    + dBn    * dt)
    state$PC    <- max(0, s$PC    + dPC    * dt)
    state$Ab    <- max(0, s$Ab    + dAb    * dt)
    state$CFib  <- max(0, s$CFib  + dCFib  * dt)
    state$MFib  <- max(0, s$MFib  + dMFib  * dt)
    state$Col   <- max(0, s$Col   + dCol   * dt)
    state$Tn    <- max(0, s$Tn    + dTn    * dt)
    state$BNP   <- max(0, s$BNP   + dBNP   * dt)
    state$EF    <- min(75, max(10, s$EF + dEF * dt))
    state$IVIG_C <- max(0, s$IVIG_C + dIVIG_C * dt)
    state$PRED_C <- max(0, s$PRED_C + dPRED_C * dt)
    state$AZA_C  <- max(0, s$AZA_C  + dAZA_C  * dt)
    state$MP6_C  <- max(0, s$MP6_C  + dMP6_C  * dt)
    state$CSA_C  <- max(0, s$CSA_C  + dCSA_C  * dt)
    state$COLC_C <- max(0, s$COLC_C + dCOLC_C * dt)

    out[[idx]] <- c(time = t, unlist(state))
  }

  as.data.frame(do.call(rbind, out))
}

# ============================================================
# Default parameters
# ============================================================
default_params <- list(
  beta_V = 2e-4, p_V = 50, c_V = 8, delta_I = 1.2,
  r_H = 0.02, d_H = 0.0003, kIFN = 0.3,
  Emax_IVIG = 0.6, EC50_IVIG = 5,
  Emax_PRED = 0.8, EC50_PRED = 2,
  Emax_AZA  = 0.65, EC50_AZA = 0.5,
  Emax_CSA  = 0.75, EC50_CSA = 300,
  Emax_COLC = 0.55, EC50_COLC = 10
)

# ============================================================
# UI
# ============================================================
ui <- dashboardPage(
  skin = "red",
  dashboardHeader(title = "Myocarditis QSP Dashboard", titleWidth = 300),

  dashboardSidebar(
    width = 270,
    sidebarMenu(
      menuItem("Overview",            tabName = "overview",    icon = icon("heart-pulse")),
      menuItem("PK Profiles",         tabName = "pk",          icon = icon("chart-line")),
      menuItem("Viral & Innate",      tabName = "viral",       icon = icon("virus")),
      menuItem("PD Biomarkers",       tabName = "pd",          icon = icon("vial")),
      menuItem("Clinical Endpoints",  tabName = "clinical",    icon = icon("stethoscope")),
      menuItem("Scenario Comparison", tabName = "scenarios",   icon = icon("sliders")),
      menuItem("Fibrosis & Remodel",  tabName = "fibrosis",    icon = icon("dna"))
    ),
    hr(),
    h5("Patient Parameters", style = "padding-left:15px; color:#aaa;"),
    sliderInput("body_weight", "Body Weight (kg)", 40, 120, 70, step = 5),
    sliderInput("viral_inoculum", "Viral Inoculum (log copies/mL)", 3, 8, 5, step = 0.5),
    sliderInput("tend", "Simulation Duration (days)", 90, 730, 365, step = 30),
    hr(),
    h5("Disease Variant", style = "padding-left:15px; color:#aaa;"),
    selectInput("variant", "Myocarditis Type",
                choices = c("Viral (CVB3/COVID-19)" = "viral",
                            "Autoimmune (Lymphocytic)" = "autoimmune",
                            "Giant Cell (Fulminant)" = "gcm",
                            "Eosinophilic" = "eos"),
                selected = "viral")
  ),

  dashboardBody(
    tags$style(HTML("
      .content-wrapper { background-color: #f4f6f9; }
      .box { border-radius: 8px; }
    ")),
    tabItems(

      # ---- Tab 1: Overview ----
      tabItem(tabName = "overview",
        fluidRow(
          valueBoxOutput("vbox_tn",   width = 3),
          valueBoxOutput("vbox_bnp",  width = 3),
          valueBoxOutput("vbox_ef",   width = 3),
          valueBoxOutput("vbox_recovery", width = 3)
        ),
        fluidRow(
          box(title = "Disease Overview — Myocarditis", status = "danger",
              solidHeader = TRUE, width = 12, collapsible = TRUE,
              HTML("
              <b>Myocarditis</b> is inflammation of the myocardium (heart muscle) caused by viral infection,
              autoimmune activation, or toxic exposures. <br><br>
              <b>Key mechanisms:</b>
              <ul>
                <li>Viral phase: CVB3 / SARS-CoV-2 infects cardiomyocytes via CAR/ACE2 receptors</li>
                <li>Innate response: IFN-α/β, NK cells, macrophage M1 activation</li>
                <li>Adaptive response: Th1/CTL-mediated cardiomyocyte killing</li>
                <li>Autoimmunity: molecular mimicry → anti-cardiac antibodies (anti-myosin, anti-β1AR)</li>
                <li>Remodeling: TGF-β → myofibroblast → collagen → fibrosis → DCM</li>
              </ul>
              <b>Drug targets:</b> IVIG (Fc-R blockade), Corticosteroids (NF-κB), Azathioprine (purine synthesis),
              Cyclosporine (calcineurin/IL-2), Colchicine (NLRP3 inflammasome)
              "))
        ),
        fluidRow(
          box(title = "LVEF Trajectory", status = "primary", solidHeader = TRUE,
              width = 6, plotlyOutput("overview_ef_plot", height = 300)),
          box(title = "Troponin Kinetics", status = "warning", solidHeader = TRUE,
              width = 6, plotlyOutput("overview_tn_plot", height = 300))
        )
      ),

      # ---- Tab 2: Drug PK Profiles ----
      tabItem(tabName = "pk",
        fluidRow(
          box(title = "Drug Dosing Setup", status = "info",
              solidHeader = TRUE, width = 4,
              h5("IVIG"),
              checkboxInput("use_IVIG", "Administer IVIG (2 g/kg IV)", FALSE),
              numericInput("t_IVIG", "Day of IVIG infusion", 5, min = 0, max = 30),
              hr(),
              h5("Prednisone"),
              numericInput("dose_PRED", "Prednisone dose (mg/day)", 50, min = 0, max = 150),
              numericInput("t_PRED", "Start day", 5, min = 0, max = 30),
              numericInput("dur_PRED", "Duration (days)", 180, min = 7, max = 365),
              hr(),
              h5("Azathioprine"),
              numericInput("dose_AZA", "Azathioprine dose (mg/day)", 150, min = 0, max = 300),
              numericInput("t_AZA", "Start day", 5, min = 0, max = 30),
              numericInput("dur_AZA", "Duration (days)", 180, min = 7, max = 365),
              hr(),
              h5("Cyclosporine"),
              numericInput("dose_CSA", "Cyclosporine dose (mg/day)", 200, min = 0, max = 600),
              numericInput("t_CSA", "Start day", 5, min = 0, max = 30),
              numericInput("dur_CSA", "Duration (days)", 180, min = 7, max = 365),
              hr(),
              h5("Colchicine"),
              numericInput("dose_COLC", "Colchicine dose (mg BID)", 0.5, min = 0, max = 1.5, step = 0.5),
              numericInput("t_COLC", "Start day", 5, min = 0, max = 30),
              numericInput("dur_COLC", "Duration (days)", 90, min = 7, max = 365),
              br(),
              actionButton("run_sim", "Run Simulation", class = "btn-danger btn-lg", width = "100%")
          ),
          box(title = "PK Concentration-Time Profiles", status = "primary",
              solidHeader = TRUE, width = 8,
              plotlyOutput("pk_plot", height = 500))
        )
      ),

      # ---- Tab 3: Viral & Innate Immunity ----
      tabItem(tabName = "viral",
        fluidRow(
          box(title = "Viral Load Kinetics", status = "danger",
              solidHeader = TRUE, width = 6, plotlyOutput("viral_plot", height = 300)),
          box(title = "Cardiomyocyte Dynamics", status = "warning",
              solidHeader = TRUE, width = 6, plotlyOutput("cmc_plot", height = 300))
        ),
        fluidRow(
          box(title = "Innate Immune Cells", status = "success",
              solidHeader = TRUE, width = 6, plotlyOutput("innate_cells_plot", height = 300)),
          box(title = "Innate Cytokines", status = "info",
              solidHeader = TRUE, width = 6, plotlyOutput("innate_cyt_plot", height = 300))
        )
      ),

      # ---- Tab 4: PD Biomarkers ----
      tabItem(tabName = "pd",
        fluidRow(
          box(title = "Troponin I (Primary Cardiac Damage Marker)", status = "danger",
              solidHeader = TRUE, width = 6, plotlyOutput("tn_plot", height = 300)),
          box(title = "BNP (Heart Failure Marker)", status = "warning",
              solidHeader = TRUE, width = 6, plotlyOutput("bnp_plot", height = 300))
        ),
        fluidRow(
          box(title = "Drug Effect (PD) Profiles", status = "primary",
              solidHeader = TRUE, width = 6, plotlyOutput("pd_effect_plot", height = 300)),
          box(title = "Anti-cardiac Antibodies (Autoimmune)", status = "info",
              solidHeader = TRUE, width = 6, plotlyOutput("ab_plot", height = 300))
        )
      ),

      # ---- Tab 5: Clinical Endpoints ----
      tabItem(tabName = "clinical",
        fluidRow(
          box(title = "LVEF Trajectory", status = "primary",
              solidHeader = TRUE, width = 8, plotlyOutput("ef_plot", height = 400)),
          box(title = "Clinical Outcome Prediction", status = "success",
              solidHeader = TRUE, width = 4,
              verbatimTextOutput("outcome_text"),
              br(),
              h5("Reference Thresholds:"),
              tags$ul(
                tags$li("EF ≥ 50%: Complete recovery"),
                tags$li("EF 35–50%: Partial recovery"),
                tags$li("EF < 35%: Dilated CMP risk"),
                tags$li("Troponin > 10 ng/mL: Severe injury"),
                tags$li("BNP > 400 pg/mL: Decompensated HF")
              ))
        ),
        fluidRow(
          box(title = "Adaptive Immune Response", status = "info",
              solidHeader = TRUE, width = 6, plotlyOutput("adaptive_plot", height = 300)),
          box(title = "Biomarker Summary Table", status = "warning",
              solidHeader = TRUE, width = 6, DTOutput("biomarker_table"))
        )
      ),

      # ---- Tab 6: Scenario Comparison ----
      tabItem(tabName = "scenarios",
        fluidRow(
          box(title = "Treatment Scenario Configuration", status = "info",
              solidHeader = TRUE, width = 12,
              column(4,
                h5("Scenario A (Blue)"),
                selectInput("scA", NULL, c("No Treatment", "IVIG only",
                  "Prednisone+Aza", "Triple IS (Pred+Aza+CsA)", "IVIG+Colchicine"),
                  selected = "No Treatment")
              ),
              column(4,
                h5("Scenario B (Green)"),
                selectInput("scB", NULL, c("No Treatment", "IVIG only",
                  "Prednisone+Aza", "Triple IS (Pred+Aza+CsA)", "IVIG+Colchicine"),
                  selected = "Prednisone+Aza")
              ),
              column(4,
                h5("Scenario C (Red)"),
                selectInput("scC", NULL, c("No Treatment", "IVIG only",
                  "Prednisone+Aza", "Triple IS (Pred+Aza+CsA)", "IVIG+Colchicine"),
                  selected = "Triple IS (Pred+Aza+CsA)")
              )
          )
        ),
        fluidRow(
          box(title = "LVEF Comparison", status = "primary",
              solidHeader = TRUE, width = 6, plotlyOutput("comp_ef_plot", height = 350)),
          box(title = "Troponin Comparison", status = "danger",
              solidHeader = TRUE, width = 6, plotlyOutput("comp_tn_plot", height = 350))
        ),
        fluidRow(
          box(title = "Fibrosis Comparison", status = "warning",
              solidHeader = TRUE, width = 6, plotlyOutput("comp_fib_plot", height = 350)),
          box(title = "Comparison Summary Table", status = "success",
              solidHeader = TRUE, width = 6, DTOutput("comp_table"))
        )
      ),

      # ---- Tab 7: Fibrosis & Remodeling ----
      tabItem(tabName = "fibrosis",
        fluidRow(
          box(title = "Fibrosis Progression Model", status = "warning",
              solidHeader = TRUE, width = 12,
              HTML("<p>Myocardial fibrosis in myocarditis is driven by:</p>
                   <ul>
                   <li><b>TGF-β:</b> Key driver of cardiac fibroblast → myofibroblast transition</li>
                   <li><b>TNF-α / IL-1β:</b> Pro-inflammatory activation of fibroblasts</li>
                   <li><b>Angiotensin II:</b> RAAS-mediated profibrotic signaling</li>
                   <li><b>MMP/TIMP imbalance:</b> Net extracellular matrix accumulation</li>
                   </ul>
                   <p>Late gadolinium enhancement (LGE) on CMR correlates with fibrosis and
                   predicts arrhythmia risk and DCM progression.</p>")
          )
        ),
        fluidRow(
          box(title = "Myofibroblast & Collagen Dynamics", status = "warning",
              solidHeader = TRUE, width = 6, plotlyOutput("fibrosis_plot", height = 350)),
          box(title = "Cytokine Drivers of Fibrosis", status = "danger",
              solidHeader = TRUE, width = 6, plotlyOutput("fibrosis_cyto_plot", height = 350))
        ),
        fluidRow(
          box(title = "TGF-β Kinetics", status = "info",
              solidHeader = TRUE, width = 6, plotlyOutput("tgfb_plot", height = 300)),
          box(title = "EF vs Fibrosis Correlation", status = "primary",
              solidHeader = TRUE, width = 6, plotlyOutput("ef_fib_plot", height = 300))
        )
      )
    )
  )
)

# ============================================================
# Server
# ============================================================
server <- function(input, output, session) {

  # Reactive simulation
  sim_data <- reactiveVal(NULL)

  # Build events from inputs
  make_events <- function() {
    events <- list(list(cmt = "V", amt = 10^input$viral_inoculum, time = 0))

    if (isTRUE(input$use_IVIG)) {
      ivig_conc <- 2000 * input$body_weight / (3.5 * input$body_weight)  # mg/mL
      events <- c(events, list(list(cmt = "IVIG_C", amt = ivig_conc, time = input$t_IVIG)))
    }
    # Prednisone (daily bolus approximation)
    if (input$dose_PRED > 0) {
      for (day in seq(input$t_PRED, input$t_PRED + input$dur_PRED - 1, by = 1)) {
        events <- c(events, list(list(cmt = "PRED_C", amt = input$dose_PRED * 1000 /
                                        (0.97 * input$body_weight), time = day)))
      }
    }
    # Azathioprine
    if (input$dose_AZA > 0) {
      for (day in seq(input$t_AZA, input$t_AZA + input$dur_AZA - 1, by = 1)) {
        events <- c(events, list(list(cmt = "AZA_C", amt = input$dose_AZA /
                                        (0.8 * input$body_weight), time = day)))
      }
    }
    # Cyclosporine
    if (input$dose_CSA > 0) {
      for (day in seq(input$t_CSA, input$t_CSA + input$dur_CSA - 1, by = 1)) {
        events <- c(events, list(list(cmt = "CSA_C", amt = input$dose_CSA * 1000 /
                                        (4 * input$body_weight), time = day)))
      }
    }
    # Colchicine
    if (input$dose_COLC > 0) {
      for (day in seq(input$t_COLC, input$t_COLC + input$dur_COLC - 1, by = 0.5)) {
        events <- c(events, list(list(cmt = "COLC_C", amt = input$dose_COLC * 1000 /
                                        (250 * input$body_weight), time = day)))
      }
    }
    events
  }

  # Run simulation on button press or on initial load
  observeEvent(input$run_sim, {
    withProgress(message = "Running simulation...", value = 0, {
      events <- make_events()
      result <- run_myocarditis_sim(default_params, events, tend = input$tend, dt = 0.5)
      sim_data(result)
      incProgress(1)
    })
  })

  # Auto-run on load
  observe({
    if (is.null(sim_data())) {
      events <- list(list(cmt = "V", amt = 1e5, time = 0))
      result <- run_myocarditis_sim(default_params, events, tend = 365, dt = 0.5)
      sim_data(result)
    }
  })

  # ---- Value Boxes ----
  output$vbox_tn <- renderValueBox({
    df <- sim_data()
    req(df)
    peak_tn <- max(df$Tn, na.rm = TRUE)
    valueBox(sprintf("%.1f ng/mL", peak_tn), "Peak Troponin I",
             icon = icon("heart"), color = ifelse(peak_tn > 10, "red", "yellow"))
  })
  output$vbox_bnp <- renderValueBox({
    df <- sim_data()
    req(df)
    peak_bnp <- max(df$BNP, na.rm = TRUE)
    valueBox(sprintf("%.0f pg/mL", peak_bnp), "Peak BNP",
             icon = icon("lungs"), color = ifelse(peak_bnp > 400, "red", "orange"))
  })
  output$vbox_ef <- renderValueBox({
    df <- sim_data()
    req(df)
    nadir_ef <- min(df$EF, na.rm = TRUE)
    valueBox(sprintf("%.0f%%", nadir_ef), "Nadir LVEF",
             icon = icon("stethoscope"), color = ifelse(nadir_ef < 35, "red",
                                                         ifelse(nadir_ef < 50, "yellow", "green")))
  })
  output$vbox_recovery <- renderValueBox({
    df <- sim_data()
    req(df)
    final_ef <- tail(df$EF, 1)
    recovery <- ifelse(final_ef >= 50, "Recovery", ifelse(final_ef >= 35, "Partial", "DCM Risk"))
    valueBox(recovery, "6-Month Outcome",
             icon = icon("chart-line"),
             color = ifelse(recovery == "Recovery", "green",
                            ifelse(recovery == "Partial", "yellow", "red")))
  })

  # ---- Overview Plots ----
  output$overview_ef_plot <- renderPlotly({
    df <- sim_data()
    req(df)
    plot_ly(df, x = ~time, y = ~EF, type = "scatter", mode = "lines",
            line = list(color = "steelblue", width = 2)) %>%
      add_segments(x = 0, xend = max(df$time), y = 50, yend = 50,
                   line = list(color = "gray", dash = "dash")) %>%
      layout(xaxis = list(title = "Time (days)"),
             yaxis = list(title = "LVEF (%)"),
             showlegend = FALSE)
  })

  output$overview_tn_plot <- renderPlotly({
    df <- sim_data()
    req(df)
    plot_ly(df, x = ~time, y = ~Tn, type = "scatter", mode = "lines",
            line = list(color = "crimson", width = 2)) %>%
      layout(xaxis = list(title = "Time (days)"),
             yaxis = list(title = "Troponin I (ng/mL)", type = "log"),
             showlegend = FALSE)
  })

  # ---- PK Plot ----
  output$pk_plot <- renderPlotly({
    df <- sim_data()
    req(df)
    df_pk <- df %>%
      select(time, IVIG_C, PRED_C, MP6_C, CSA_C, COLC_C) %>%
      pivot_longer(-time, names_to = "Drug", values_to = "Conc")

    plot_ly(df_pk, x = ~time, y = ~Conc, color = ~Drug,
            type = "scatter", mode = "lines") %>%
      layout(xaxis = list(title = "Time (days)"),
             yaxis = list(title = "Concentration (various units)"),
             title = "Drug Concentration-Time Profiles")
  })

  # ---- Viral & Innate Plots ----
  output$viral_plot <- renderPlotly({
    df <- sim_data()
    req(df)
    plot_ly(df, x = ~time, y = ~pmax(V, 0.01), type = "scatter", mode = "lines",
            line = list(color = "crimson")) %>%
      layout(yaxis = list(type = "log", title = "Viral Load (copies/mL)"),
             xaxis = list(title = "Time (days)"))
  })

  output$cmc_plot <- renderPlotly({
    df <- sim_data()
    req(df)
    plot_ly(df, x = ~time, y = ~H, name = "Healthy CMC", type = "scatter", mode = "lines",
            line = list(color = "green")) %>%
      add_trace(y = ~I, name = "Infected CMC", line = list(color = "orange")) %>%
      add_trace(y = ~D, name = "Dead CMC", line = list(color = "red")) %>%
      layout(xaxis = list(title = "Time (days)"),
             yaxis = list(title = "Cardiomyocytes (cells)"))
  })

  output$innate_cells_plot <- renderPlotly({
    df <- sim_data()
    req(df)
    plot_ly(df, x = ~time, y = ~NK, name = "NK Cells", type = "scatter", mode = "lines",
            line = list(color = "navy")) %>%
      add_trace(y = ~M1, name = "M1 Macro", line = list(color = "firebrick")) %>%
      add_trace(y = ~M2, name = "M2 Macro", line = list(color = "seagreen")) %>%
      layout(xaxis = list(title = "Time (days)"),
             yaxis = list(title = "Cells (cells/uL)"))
  })

  output$innate_cyt_plot <- renderPlotly({
    df <- sim_data()
    req(df)
    plot_ly(df, x = ~time, y = ~TNFa, name = "TNF-α", type = "scatter", mode = "lines",
            line = list(color = "red")) %>%
      add_trace(y = ~IL6, name = "IL-6", line = list(color = "purple")) %>%
      add_trace(y = ~IL1b, name = "IL-1β", line = list(color = "orange")) %>%
      add_trace(y = ~IFNg, name = "IFN-γ", line = list(color = "blue")) %>%
      layout(xaxis = list(title = "Time (days)"),
             yaxis = list(title = "Cytokine (pg/mL)"))
  })

  # ---- PD Plots ----
  output$tn_plot <- renderPlotly({
    df <- sim_data()
    req(df)
    plot_ly(df, x = ~time, y = ~Tn, type = "scatter", mode = "lines",
            line = list(color = "crimson", width = 2)) %>%
      add_segments(x = 0, xend = max(df$time), y = 0.04, yend = 0.04,
                   line = list(color = "gray", dash = "dash"),
                   name = "ULN (0.04 ng/mL)") %>%
      layout(xaxis = list(title = "Time (days)"),
             yaxis = list(title = "Troponin I (ng/mL)", type = "log"))
  })

  output$bnp_plot <- renderPlotly({
    df <- sim_data()
    req(df)
    plot_ly(df, x = ~time, y = ~BNP, type = "scatter", mode = "lines",
            line = list(color = "darkorange", width = 2)) %>%
      add_segments(x = 0, xend = max(df$time), y = 100, yend = 100,
                   line = list(color = "gray", dash = "dash")) %>%
      layout(xaxis = list(title = "Time (days)"),
             yaxis = list(title = "BNP (pg/mL)"))
  })

  output$pd_effect_plot <- renderPlotly({
    df <- sim_data()
    req(df)
    # Compute drug effects
    p <- default_params
    df_eff <- df %>%
      mutate(
        E_PRED = p$Emax_PRED * PRED_C / (p$EC50_PRED + PRED_C + 1e-9),
        E_AZA  = p$Emax_AZA  * MP6_C  / (p$EC50_AZA  + MP6_C  + 1e-9),
        E_CSA  = p$Emax_CSA  * CSA_C  / (p$EC50_CSA   + CSA_C  + 1e-9),
        E_IVIG = p$Emax_IVIG * IVIG_C / (p$EC50_IVIG  + IVIG_C + 1e-9)
      )
    plot_ly(df_eff, x = ~time, y = ~E_PRED, name = "Prednisone PD",
            type = "scatter", mode = "lines", line = list(color = "green")) %>%
      add_trace(y = ~E_AZA, name = "Azathioprine PD", line = list(color = "blue")) %>%
      add_trace(y = ~E_CSA, name = "Cyclosporine PD", line = list(color = "purple")) %>%
      add_trace(y = ~E_IVIG, name = "IVIG PD", line = list(color = "orange")) %>%
      layout(xaxis = list(title = "Time (days)"),
             yaxis = list(title = "Drug Effect (fraction)", range = c(0, 1)))
  })

  output$ab_plot <- renderPlotly({
    df <- sim_data()
    req(df)
    plot_ly(df, x = ~time, y = ~Ab, type = "scatter", mode = "lines",
            line = list(color = "purple", width = 2)) %>%
      layout(xaxis = list(title = "Time (days)"),
             yaxis = list(title = "Anti-cardiac Antibodies (AU/mL)"))
  })

  # ---- Clinical Endpoint Plots ----
  output$ef_plot <- renderPlotly({
    df <- sim_data()
    req(df)
    plot_ly(df, x = ~time, y = ~EF, type = "scatter", mode = "lines",
            line = list(color = "steelblue", width = 2.5)) %>%
      add_segments(x = 0, xend = max(df$time), y = 50, yend = 50,
                   line = list(color = "gray", dash = "dash"), name = "50% threshold") %>%
      add_segments(x = 0, xend = max(df$time), y = 35, yend = 35,
                   line = list(color = "red", dash = "dot"), name = "35% DCM threshold") %>%
      layout(xaxis = list(title = "Time (days)"),
             yaxis = list(title = "LVEF (%)"),
             title = "LVEF Trajectory — Treatment Response")
  })

  output$outcome_text <- renderText({
    df <- sim_data()
    req(df)
    final_ef  <- round(tail(df$EF, 1), 1)
    nadir_ef  <- round(min(df$EF), 1)
    peak_tn   <- round(max(df$Tn), 2)
    peak_bnp  <- round(max(df$BNP), 0)
    final_col <- round(tail(df$Col, 1), 2)
    recovery  <- ifelse(final_ef >= 50, "COMPLETE RECOVERY",
                        ifelse(final_ef >= 35, "PARTIAL RECOVERY", "RISK OF DCM"))
    sprintf(
      "OUTCOME: %s\n\nLVEF:\n  Nadir: %.1f%%\n  Final: %.1f%%\n\nBiomarkers:\n  Peak Troponin: %.2f ng/mL\n  Peak BNP: %.0f pg/mL\n\nFibrosis:\n  Final Collagen: %.2f AU\n\n%s",
      recovery, nadir_ef, final_ef, peak_tn, peak_bnp, final_col,
      ifelse(final_ef >= 50, "Patient predicted to achieve complete myocardial recovery.",
             ifelse(final_ef >= 35, "Partial recovery expected. Close follow-up recommended.",
                    "High risk of progression to dilated cardiomyopathy."))
    )
  })

  output$adaptive_plot <- renderPlotly({
    df <- sim_data()
    req(df)
    plot_ly(df, x = ~time, y = ~Th1, name = "Th1 Cells", type = "scatter", mode = "lines",
            line = list(color = "steelblue")) %>%
      add_trace(y = ~CTL, name = "CD8+ CTL", line = list(color = "darkred")) %>%
      add_trace(y = ~Treg, name = "Treg", line = list(color = "green")) %>%
      add_trace(y = ~Th17, name = "Th17", line = list(color = "orange")) %>%
      layout(xaxis = list(title = "Time (days)"),
             yaxis = list(title = "T Cells (cells/uL)"))
  })

  output$biomarker_table <- renderDT({
    df <- sim_data()
    req(df)
    key_times <- c(7, 14, 30, 90, 180, 365)
    df_sum <- df %>%
      mutate(time_round = round(time)) %>%
      filter(time_round %in% key_times) %>%
      group_by(time_round) %>%
      summarise(
        `Troponin I (ng/mL)` = round(mean(Tn), 3),
        `BNP (pg/mL)` = round(mean(BNP), 0),
        `LVEF (%)` = round(mean(EF), 1),
        `Viral Load` = sprintf("%.0e", mean(V)),
        `Th1 (cells/uL)` = round(mean(Th1), 0),
        `Fibrosis (AU)` = round(mean(Col), 2),
        .groups = "drop"
      ) %>%
      rename(`Day` = time_round)
    datatable(df_sum, options = list(pageLength = 10), rownames = FALSE)
  })

  # ---- Scenario comparison (simplified with preset scenarios) ----
  get_scenario_data <- function(sc_name) {
    base_events <- list(list(cmt = "V", amt = 1e5, time = 0))
    extra <- switch(sc_name,
      "No Treatment" = list(),
      "IVIG only" = list(list(cmt = "IVIG_C", amt = 571, time = 5)),
      "Prednisone+Aza" = {
        evs <- list()
        for (d in 5:184) {
          evs <- c(evs, list(list(cmt = "PRED_C", amt = 50000 / (0.97 * 70), time = d)))
          evs <- c(evs, list(list(cmt = "AZA_C", amt = 150 / (0.8 * 70), time = d)))
        }
        evs
      },
      "Triple IS (Pred+Aza+CsA)" = {
        evs <- list(list(cmt = "IVIG_C", amt = 571, time = 5))
        for (d in 5:364) {
          evs <- c(evs, list(list(cmt = "PRED_C", amt = 50000 / (0.97 * 70), time = d)))
          evs <- c(evs, list(list(cmt = "AZA_C", amt = 150 / (0.8 * 70), time = d)))
          evs <- c(evs, list(list(cmt = "CSA_C", amt = 200000 / (4 * 70), time = d)))
        }
        evs
      },
      "IVIG+Colchicine" = {
        evs <- list(list(cmt = "IVIG_C", amt = 571, time = 5))
        for (d in seq(5, 94, 0.5)) {
          evs <- c(evs, list(list(cmt = "COLC_C", amt = 500 / (250 * 70), time = d)))
        }
        evs
      }
    )
    events <- c(base_events, extra)
    run_myocarditis_sim(default_params, events, tend = 365, dt = 1)
  }

  sc_results <- reactive({
    list(
      A = get_scenario_data(input$scA),
      B = get_scenario_data(input$scB),
      C = get_scenario_data(input$scC)
    )
  })

  output$comp_ef_plot <- renderPlotly({
    res <- sc_results()
    plot_ly() %>%
      add_trace(data = res$A, x = ~time, y = ~EF, name = input$scA,
                type = "scatter", mode = "lines", line = list(color = "blue")) %>%
      add_trace(data = res$B, x = ~time, y = ~EF, name = input$scB,
                type = "scatter", mode = "lines", line = list(color = "green")) %>%
      add_trace(data = res$C, x = ~time, y = ~EF, name = input$scC,
                type = "scatter", mode = "lines", line = list(color = "red")) %>%
      add_segments(x = 0, xend = 365, y = 50, yend = 50,
                   line = list(color = "gray", dash = "dash")) %>%
      layout(xaxis = list(title = "Time (days)"),
             yaxis = list(title = "LVEF (%)"))
  })

  output$comp_tn_plot <- renderPlotly({
    res <- sc_results()
    plot_ly() %>%
      add_trace(data = res$A, x = ~time, y = ~pmax(Tn, 0.01), name = input$scA,
                type = "scatter", mode = "lines", line = list(color = "blue")) %>%
      add_trace(data = res$B, x = ~time, y = ~pmax(Tn, 0.01), name = input$scB,
                type = "scatter", mode = "lines", line = list(color = "green")) %>%
      add_trace(data = res$C, x = ~time, y = ~pmax(Tn, 0.01), name = input$scC,
                type = "scatter", mode = "lines", line = list(color = "red")) %>%
      layout(xaxis = list(title = "Time (days)"),
             yaxis = list(title = "Troponin I (ng/mL)", type = "log"))
  })

  output$comp_fib_plot <- renderPlotly({
    res <- sc_results()
    plot_ly() %>%
      add_trace(data = res$A, x = ~time, y = ~Col, name = input$scA,
                type = "scatter", mode = "lines", line = list(color = "blue")) %>%
      add_trace(data = res$B, x = ~time, y = ~Col, name = input$scB,
                type = "scatter", mode = "lines", line = list(color = "green")) %>%
      add_trace(data = res$C, x = ~time, y = ~Col, name = input$scC,
                type = "scatter", mode = "lines", line = list(color = "red")) %>%
      layout(xaxis = list(title = "Time (days)"),
             yaxis = list(title = "Myocardial Collagen (AU)"))
  })

  output$comp_table <- renderDT({
    res <- sc_results()
    make_row <- function(df, name) {
      df_180 <- df %>% filter(abs(time - 180) < 1) %>% tail(1)
      data.frame(
        Scenario    = name,
        `EF at 6mo` = sprintf("%.1f%%", df_180$EF),
        `Tn peak`   = sprintf("%.2f", max(df$Tn)),
        `BNP peak`  = sprintf("%.0f", max(df$BNP)),
        `Col 6mo`   = sprintf("%.2f", df_180$Col),
        Outcome     = ifelse(df_180$EF >= 50, "Recovery",
                             ifelse(df_180$EF >= 35, "Partial", "DCM"))
      )
    }
    tab <- rbind(
      make_row(res$A, input$scA),
      make_row(res$B, input$scB),
      make_row(res$C, input$scC)
    )
    datatable(tab, rownames = FALSE, options = list(dom = "t"))
  })

  # ---- Fibrosis Tab ----
  output$fibrosis_plot <- renderPlotly({
    df <- sim_data()
    req(df)
    plot_ly(df, x = ~time, y = ~MFib, name = "Myofibroblasts",
            type = "scatter", mode = "lines", line = list(color = "orange")) %>%
      add_trace(y = ~Col, name = "Collagen (AU)",
                line = list(color = "brown", dash = "dash")) %>%
      layout(xaxis = list(title = "Time (days)"),
             yaxis = list(title = "AU or cells/uL"))
  })

  output$fibrosis_cyto_plot <- renderPlotly({
    df <- sim_data()
    req(df)
    plot_ly(df, x = ~time, y = ~TGFb, name = "TGF-β",
            type = "scatter", mode = "lines", line = list(color = "darkgreen")) %>%
      add_trace(y = ~TNFa, name = "TNF-α", line = list(color = "red")) %>%
      add_trace(y = ~IL6, name = "IL-6", line = list(color = "purple")) %>%
      layout(xaxis = list(title = "Time (days)"),
             yaxis = list(title = "Cytokine (pg/mL)"))
  })

  output$tgfb_plot <- renderPlotly({
    df <- sim_data()
    req(df)
    plot_ly(df, x = ~time, y = ~TGFb, type = "scatter", mode = "lines",
            line = list(color = "forestgreen", width = 2)) %>%
      layout(xaxis = list(title = "Time (days)"),
             yaxis = list(title = "TGF-β (pg/mL)"))
  })

  output$ef_fib_plot <- renderPlotly({
    df <- sim_data()
    req(df)
    plot_ly(df, x = ~Col, y = ~EF, type = "scatter", mode = "markers",
            marker = list(color = ~time, colorscale = "Viridis", size = 5,
                          showscale = TRUE,
                          colorbar = list(title = "Time (days)"))) %>%
      layout(xaxis = list(title = "Myocardial Collagen (AU)"),
             yaxis = list(title = "LVEF (%)"))
  })
}

# ============================================================
# Launch
# ============================================================
shinyApp(ui = ui, server = server)
