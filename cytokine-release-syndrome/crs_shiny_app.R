## =============================================================================
## CRS QSP Shiny Dashboard — Cytokine Release Syndrome
## =============================================================================
## 6 tabs: Patient Profile · Drug PK · Cytokine Storm · Clinical Endpoints ·
##         Scenario Comparison · Biomarkers & CRS Grading
## =============================================================================

library(shiny)
library(shinydashboard)
library(ggplot2)
library(dplyr)
library(tidyr)
library(plotly)

# ==============================================================================
# Euler ODE solver (simplified CRS model for real-time Shiny interaction)
# ==============================================================================

run_crs_sim <- function(
    sim_days       = 28,
    dt             = 0.1,
    cart_type      = "CD19 CAR-T (CD28)",
    tumor_burden   = 1.0,
    infusion_dose  = 1.0,
    use_toci       = FALSE,
    use_dex        = FALSE,
    use_rux        = FALSE,
    use_ana        = FALSE,
    use_silt       = FALSE,
    toci_day       = 3,
    toci_dose_mg   = 560,
    dex_day        = 2,
    dex_dose_mg    = 10,
    rux_day        = 2,
    rux_dose_mg    = 10,
    ana_day        = 2,
    ana_dose_mg    = 100,
    silt_day       = 3,
    silt_dose_mg   = 770
) {

  # ---- Parameter set (from mrgsolve model) ----
  p <- list(
    kCART_act=0.6, kCART_prolif=1.8, kCART_death=0.25, kCART_exh=0.12,
    kIL2_prolif=20, kTumor_kill=0.4, kTumor_grow=0.05,
    kIFNg_prod=0.8, kIL2_prod=0.5, kGMCSF_prod=0.4, kTNFaT_prod=0.3,
    dIFNg=1.5, dIL2=6.0, dGMCSF=2.4, dTNFaT=3.0,
    kMAC_act=0.5, kMAC_deact=0.3, EC50IFNgMAC=50, EC50GMCSF=20,
    kIL6_prod=2.0, kIL1B_prod=0.6, kTNFaM_prod=0.8,
    dIL6=3.0, dIL1B=4.0, dTNFaM=3.0, IL6_feedback=0.3,
    kSTAT3_act=0.8, kSTAT3_deact=1.2, EC50IL6STAT3=30,
    kCRP_prod=0.4, dCRP=0.7, kFerr_prod=0.25, dFerr=0.15,
    kENDO_act=0.4, kENDO_deact=0.6, kVASC_leak=0.5, dVASC=0.8, EC50IL6ENDO=80,
    kCRS_onset=0.6, kCRS_resolve=0.4, kORGAN_dmg=0.08, dORGAN=0.05,
    # Tocilizumab PK
    CL_TOCI=0.22, V1_TOCI=3.1, Q_TOCI=0.29, V2_TOCI=2.9,
    Emax_TOCI=0.92, EC50_TOCI=1.5,
    # Siltuximab PK
    CL_SILT=0.015, V1_SILT=5.5, Emax_SILT=0.90, EC50_SILT=2.0,
    # Dexamethasone PK
    CL_DEX=288, V1_DEX=25, Q_DEX=8, V2_DEX=65,
    Emax_DEX=0.75, EC50_DEX=0.3,
    # Ruxolitinib PK
    ka_RUX=6, CL_RUX=432, V1_RUX=72, Q_RUX=12, V2_RUX=144,
    Emax_RUX=0.85, EC50_RUX=80,
    # Anakinra PK
    ka_ANA=4, CL_ANA=35, V1_ANA=18, Emax_ANA=0.90, EC50_ANA=30
  )

  # ---- CAR type adjustments ----
  if (cart_type == "CD19 CAR-T (CD28)") {
    p$kCART_prolif <- 2.0; p$kCART_act <- 0.7   # faster, more CRS
  } else if (cart_type == "CD19 CAR-T (4-1BB)") {
    p$kCART_prolif <- 1.5; p$kCART_act <- 0.5   # slower expansion, milder CRS
  } else if (cart_type == "BCMA CAR-T") {
    p$kCART_prolif <- 1.8; p$kIFNg_prod <- 0.7
  } else if (cart_type == "CD20 Bispecific Ab") {
    p$kCART_prolif <- 1.2; p$kCART_act <- 0.8; p$kCART_death <- 0.3
  } else if (cart_type == "BCMA Bispecific Ab") {
    p$kCART_prolif <- 1.1; p$kCART_act <- 0.75; p$kCART_death <- 0.35
  }

  # Scale by tumor burden and dose
  p$kCART_act <- p$kCART_act * tumor_burden * infusion_dose

  # ---- Initialize state variables ----
  s <- list(
    CART = 0.1 * infusion_dose, CART_exh = 0,
    TUMOR = tumor_burden,
    IFNg = 5, IL2 = 8, GMCSF = 2, TNFaT = 3,
    MAC = 0.05, IL6 = 10, IL1B = 3, TNFaM = 5,
    STAT3 = 0.1, CRP = 3, FERRITIN = 150,
    ENDO = 0.05, VASC = 0.1,
    CRS = 0, ORGAN = 0,
    TOCI_C1=0, TOCI_C2=0, SILT_C1=0,
    DEX_C1=0, DEX_C2=0,
    RUX_GUT=0, RUX_C1=0, RUX_C2=0,
    ANA_SC=0, ANA_C1=0
  )

  # Pre-compute dose concentrations
  TOCI_bolus <- toci_dose_mg * 1000 / 3100    # μg/mL
  SILT_bolus <- silt_dose_mg * 1000 / 5500    # μg/mL
  DEX_bolus  <- dex_dose_mg  * 1e6  / 25000   # ng/mL
  RUX_bolus  <- rux_dose_mg  * 1e6  / 72000   # ng/mL depot
  ANA_bolus  <- ana_dose_mg  * 1e6  / 18000   # ng/mL depot
  toci_given <- FALSE; toci2_given <- FALSE
  silt_given <- FALSE
  dex_doses_given <- 0
  rux_last_dose_day <- -1
  ana_last_dose_day <- -1

  n_steps <- round(sim_days / dt)
  times <- seq(0, sim_days, by=dt)
  n_out <- length(times)

  # Output storage
  out <- data.frame(
    time=times,
    IFNg=NA, IL6=NA, IL1B=NA, CRP=NA, FERRITIN=NA,
    TNFa=NA, CART=NA, TUMOR=NA, GMCSF=NA,
    MAC=NA, STAT3=NA, ENDO=NA, VASC=NA,
    CRS_SEV=NA, ORGAN=NA, Temp=NA, Survival=NA,
    TOCI_C1=NA, DEX_C1=NA, RUX_C1=NA, ANA_C1=NA, SILT_C1=NA,
    Fibrinogen=NA, CRS_Grade=NA, ICANS=NA, CRS_Index=NA
  )

  for (i in seq_along(times)) {
    t <- times[i]
    # ---- Dosing events ----
    if (use_toci && !toci_given && t >= toci_day) {
      s$TOCI_C1 <- s$TOCI_C1 + TOCI_bolus; toci_given <- TRUE
    }
    if (use_toci && toci_given && !toci2_given && t >= toci_day + 2) {
      s$TOCI_C1 <- s$TOCI_C1 + TOCI_bolus * 0.5; toci2_given <- TRUE
    }
    if (use_silt && !silt_given && t >= silt_day) {
      s$SILT_C1 <- s$SILT_C1 + SILT_bolus; silt_given <- TRUE
    }
    if (use_dex && t >= dex_day && dex_doses_given < 7) {
      day_floor <- floor(t)
      if (day_floor >= (dex_day + dex_doses_given) && (i==1 || floor(times[i-1]) < day_floor)) {
        s$DEX_C1 <- s$DEX_C1 + DEX_bolus * (0.9^dex_doses_given)
        dex_doses_given <- dex_doses_given + 1
      }
    }
    if (use_rux && t >= rux_day) {
      day_floor <- floor(t * 2) / 2  # BID dosing every 0.5 days
      if (day_floor > rux_last_dose_day) {
        s$RUX_GUT <- s$RUX_GUT + RUX_bolus; rux_last_dose_day <- day_floor
      }
    }
    if (use_ana && t >= ana_day) {
      day_floor <- floor(t)
      if (day_floor > ana_last_dose_day) {
        s$ANA_SC <- s$ANA_SC + ANA_bolus; ana_last_dose_day <- day_floor
      }
    }

    # ---- Drug PK ODEs ----
    ke_TOCI <- p$CL_TOCI / p$V1_TOCI
    k12_TOCI <- p$Q_TOCI / p$V1_TOCI; k21_TOCI <- p$Q_TOCI / p$V2_TOCI
    dTOCI_C1 <- -(ke_TOCI + k12_TOCI) * s$TOCI_C1 + k21_TOCI * s$TOCI_C2
    dTOCI_C2 <- k12_TOCI * s$TOCI_C1 - k21_TOCI * s$TOCI_C2

    ke_SILT <- p$CL_SILT / p$V1_SILT
    dSILT_C1 <- -ke_SILT * s$SILT_C1

    ke_DEX <- p$CL_DEX / p$V1_DEX; k12_DEX <- p$Q_DEX / p$V1_DEX; k21_DEX <- p$Q_DEX / p$V2_DEX
    dDEX_C1 <- -(ke_DEX + k12_DEX) * s$DEX_C1 + k21_DEX * s$DEX_C2
    dDEX_C2 <- k12_DEX * s$DEX_C1 - k21_DEX * s$DEX_C2

    ke_RUX <- p$CL_RUX / p$V1_RUX; k12_RUX <- p$Q_RUX / p$V1_RUX; k21_RUX <- p$Q_RUX / p$V2_RUX
    dRUX_GUT <- -p$ka_RUX * s$RUX_GUT
    dRUX_C1 <- p$ka_RUX * s$RUX_GUT - (ke_RUX + k12_RUX) * s$RUX_C1 + k21_RUX * s$RUX_C2
    dRUX_C2 <- k12_RUX * s$RUX_C1 - k21_RUX * s$RUX_C2

    ke_ANA <- p$CL_ANA / p$V1_ANA
    dANA_SC <- -p$ka_ANA * s$ANA_SC
    dANA_C1 <- p$ka_ANA * s$ANA_SC - ke_ANA * s$ANA_C1

    # ---- Drug PD ----
    Inh_TOCI <- p$Emax_TOCI * s$TOCI_C1 / (p$EC50_TOCI + s$TOCI_C1 + 1e-9)
    Inh_SILT <- p$Emax_SILT * s$SILT_C1 / (p$EC50_SILT + s$SILT_C1 + 1e-9)
    Inh_DEX_IL6  <- p$Emax_DEX * s$DEX_C1 / (p$EC50_DEX + s$DEX_C1 + 1e-9)
    Inh_DEX_IFNg <- p$Emax_DEX * 0.8 * s$DEX_C1 / (p$EC50_DEX + s$DEX_C1 + 1e-9)
    Inh_RUX <- p$Emax_RUX * s$RUX_C1 / (p$EC50_RUX + s$RUX_C1 + 1e-9)
    Inh_ANA <- p$Emax_ANA * s$ANA_C1 / (p$EC50_ANA + s$ANA_C1 + 1e-9)
    Inh_IL6_sig <- 1 - (1 - Inh_TOCI) * (1 - Inh_SILT) * (1 - Inh_DEX_IL6) * (1 - Inh_RUX)

    # ---- CAR-T dynamics ----
    IL2_prolif <- s$IL2 / (p$kIL2_prolif + s$IL2)
    tumor_stim <- s$TUMOR / (0.5 + s$TUMOR)
    dCART <- p$kCART_act * tumor_stim + p$kCART_prolif * IL2_prolif * s$CART -
             (p$kCART_death + p$kCART_exh * s$CART) * s$CART
    dTUMOR <- p$kTumor_grow * s$TUMOR - p$kTumor_kill * s$CART * s$TUMOR

    # ---- Cytokine ODEs ----
    dIFNg <- p$kIFNg_prod * s$CART * (1 - Inh_DEX_IFNg) * (1 - Inh_RUX * 0.5) + 2 - p$dIFNg * s$IFNg
    dIL2  <- p$kIL2_prod * s$CART + 2 - p$dIL2 * s$IL2
    dGMCSF <- p$kGMCSF_prod * s$CART + 0.5 - p$dGMCSF * s$GMCSF
    dTNFaT <- p$kTNFaT_prod * s$CART * (1 - Inh_DEX_IL6 * 0.7) + 1 - p$dTNFaT * s$TNFaT

    # ---- Macrophage ----
    IFNg_stim  <- s$IFNg  / (p$EC50IFNgMAC + s$IFNg)
    GMCSF_stim <- s$GMCSF / (p$EC50GMCSF   + s$GMCSF)
    MAC_stim   <- 0.6 * IFNg_stim + 0.4 * GMCSF_stim
    dMAC <- p$kMAC_act * MAC_stim * (1 - s$MAC) - p$kMAC_deact * s$MAC

    IL6_amp  <- 1 + p$IL6_feedback * s$MAC
    IL6_prod <- p$kIL6_prod * s$MAC * IL6_amp * (1 - Inh_DEX_IL6 * 0.8)
    dIL6 <- IL6_prod + 3 - p$dIL6 * s$IL6

    NLRP3 <- (s$IFNg / (s$IFNg + 100)) * (s$TNFaT / (s$TNFaT + 20))
    dIL1B <- p$kIL1B_prod * s$MAC * (1 + NLRP3) * (1 - Inh_ANA * 0.8) + 0.5 - p$dIL1B * s$IL1B
    dTNFaM <- p$kTNFaM_prod * s$MAC * (1 - Inh_DEX_IL6 * 0.7) + 1 - p$dTNFaM * s$TNFaM

    # ---- IL-6 downstream ----
    IL6_sig <- (s$IL6 / (p$EC50IL6STAT3 + s$IL6)) * (1 - Inh_IL6_sig)
    dSTAT3 <- p$kSTAT3_act * IL6_sig * (1 - s$STAT3) - p$kSTAT3_deact * s$STAT3
    dCRP <- p$kCRP_prod * s$STAT3 * 50 - p$dCRP * s$CRP
    Ferr_stim <- 1 + (s$IL6 / (s$IL6 + 100)) + (s$TNFaM / (s$TNFaM + 80))
    dFERRITIN <- p$kFerr_prod * s$MAC * Ferr_stim * 500 - p$dFerr * s$FERRITIN

    # ---- Endothelial ----
    ENDO_stim <- (s$IL6 / (p$EC50IL6ENDO + s$IL6)) + 0.5 * (s$TNFaM / (50 + s$TNFaM))
    dENDO <- p$kENDO_act * ENDO_stim * (1 - s$ENDO) - p$kENDO_deact * s$ENDO
    dVASC <- p$kVASC_leak * s$ENDO - p$dVASC * s$VASC

    # ---- CRS severity / organ damage ----
    CRS_driver <- (s$IFNg / (s$IFNg + 200)) * 2 +
                  (s$IL6  / (s$IL6  + 100)) * 1.5 +
                  s$VASC * 1.5
    dCRS <- p$kCRS_onset * CRS_driver * (4 - s$CRS) / 4 -
            p$kCRS_resolve * s$CRS * (1 - Inh_IL6_sig * 0.8)
    dmg_driver <- if (s$CRS > 2.5) (s$CRS - 2.5) * 0.4 else 0
    dORGAN <- p$kORGAN_dmg * dmg_driver * (1 - s$ORGAN) - p$dORGAN * s$ORGAN

    # ---- Euler update ----
    s$CART     <- max(0, s$CART     + dCART     * dt)
    s$TUMOR    <- max(0, s$TUMOR    + dTUMOR    * dt)
    s$IFNg     <- max(0, s$IFNg     + dIFNg     * dt)
    s$IL2      <- max(0, s$IL2      + dIL2      * dt)
    s$GMCSF    <- max(0, s$GMCSF    + dGMCSF    * dt)
    s$TNFaT    <- max(0, s$TNFaT    + dTNFaT    * dt)
    s$MAC      <- min(1, max(0, s$MAC + dMAC * dt))
    s$IL6      <- max(0, s$IL6      + dIL6      * dt)
    s$IL1B     <- max(0, s$IL1B     + dIL1B     * dt)
    s$TNFaM    <- max(0, s$TNFaM    + dTNFaM    * dt)
    s$STAT3    <- min(1, max(0, s$STAT3 + dSTAT3 * dt))
    s$CRP      <- max(0, s$CRP      + dCRP      * dt)
    s$FERRITIN <- max(0, s$FERRITIN + dFERRITIN * dt)
    s$ENDO     <- min(1, max(0, s$ENDO + dENDO * dt))
    s$VASC     <- max(0, s$VASC     + dVASC     * dt)
    s$CRS      <- min(4, max(0, s$CRS + dCRS * dt))
    s$ORGAN    <- min(1, max(0, s$ORGAN + dORGAN * dt))
    s$TOCI_C1  <- max(0, s$TOCI_C1  + dTOCI_C1  * dt)
    s$TOCI_C2  <- max(0, s$TOCI_C2  + dTOCI_C2  * dt)
    s$SILT_C1  <- max(0, s$SILT_C1  + dSILT_C1  * dt)
    s$DEX_C1   <- max(0, s$DEX_C1   + dDEX_C1   * dt)
    s$DEX_C2   <- max(0, s$DEX_C2   + dDEX_C2   * dt)
    s$RUX_GUT  <- max(0, s$RUX_GUT  + dRUX_GUT  * dt)
    s$RUX_C1   <- max(0, s$RUX_C1   + dRUX_C1   * dt)
    s$RUX_C2   <- max(0, s$RUX_C2   + dRUX_C2   * dt)
    s$ANA_SC   <- max(0, s$ANA_SC   + dANA_SC   * dt)
    s$ANA_C1   <- max(0, s$ANA_C1   + dANA_C1   * dt)

    # ---- Store output ----
    out$IFNg[i]  <- s$IFNg
    out$IL6[i]   <- s$IL6
    out$IL1B[i]  <- s$IL1B
    out$CRP[i]   <- s$CRP
    out$FERRITIN[i] <- s$FERRITIN
    out$TNFa[i]  <- s$TNFaT + s$TNFaM
    out$CART[i]  <- s$CART * 100
    out$TUMOR[i] <- s$TUMOR
    out$GMCSF[i] <- s$GMCSF
    out$MAC[i]   <- s$MAC
    out$STAT3[i] <- s$STAT3
    out$ENDO[i]  <- s$ENDO
    out$VASC[i]  <- s$VASC
    out$CRS_SEV[i]  <- s$CRS
    out$ORGAN[i] <- s$ORGAN
    out$Temp[i]  <- 37 + 2.5 * (s$IL1B / (s$IL1B + 20)) + 1.5 * (s$IL6 / (s$IL6 + 100))
    out$Survival[i] <- exp(-2.0 * s$ORGAN)
    out$TOCI_C1[i]  <- s$TOCI_C1
    out$DEX_C1[i]   <- s$DEX_C1
    out$RUX_C1[i]   <- s$RUX_C1
    out$ANA_C1[i]   <- s$ANA_C1
    out$SILT_C1[i]  <- s$SILT_C1
    out$Fibrinogen[i] <- max(0.5, 4.0 - 2.5 * (s$CRS / 4))
    out$CRS_Grade[i] <- floor(s$CRS + 0.5)
    out$ICANS[i]  <- pmin(4, (s$IFNg / (s$IFNg + 500)) * 4 + (s$IL6 / (s$IL6 + 200)) * 3)
    out$CRS_Index[i] <- (s$IFNg / 500) + (s$IL6 / 200) + (s$CRP / 100) + (s$FERRITIN / 5000)
  }
  return(out)
}

# ==============================================================================
# UI
# ==============================================================================

ui <- dashboardPage(
  skin = "red",
  dashboardHeader(title = "CRS QSP Dashboard"),
  dashboardSidebar(
    sidebarMenu(
      menuItem("Patient Profile",     tabName="profile",  icon=icon("user-md")),
      menuItem("Drug PK",             tabName="pk",       icon=icon("pills")),
      menuItem("Cytokine Storm",      tabName="cyto",     icon=icon("fire")),
      menuItem("Clinical Endpoints",  tabName="endpoints",icon=icon("heartbeat")),
      menuItem("Scenario Comparison", tabName="compare",  icon=icon("chart-bar")),
      menuItem("Biomarkers & Grading",tabName="biomark",  icon=icon("vials"))
    ),
    hr(),
    h4("Simulation Settings", style="color:#fff;padding-left:15px;"),
    selectInput("cart_type","CAR-T / Therapy Type:",
                choices=c("CD19 CAR-T (CD28)","CD19 CAR-T (4-1BB)",
                          "BCMA CAR-T","CD20 Bispecific Ab","BCMA Bispecific Ab"),
                selected="CD19 CAR-T (CD28)"),
    sliderInput("tumor_burden","Tumor Burden:",min=0.1,max=3.0,value=1.0,step=0.1),
    sliderInput("infusion_dose","Infusion Dose (relative):",min=0.3,max=3.0,value=1.0,step=0.1),
    sliderInput("sim_days","Simulation (days):",min=14,max=42,value=28,step=1)
  ),

  dashboardBody(
    tags$head(tags$style(HTML("
      .content-wrapper { background-color: #f5f5f5; }
      .box { border-radius: 8px; box-shadow: 0 2px 6px rgba(0,0,0,0.1); }
      .value-box { border-radius: 8px; }
    "))),
    tabItems(

      # =====================================================================
      # TAB 1: Patient Profile
      # =====================================================================
      tabItem(tabName="profile",
        fluidRow(
          box(title="Treatment Selection", width=4, status="danger",
            checkboxInput("use_toci","Tocilizumab (IL-6R blocker)",FALSE),
            conditionalPanel("input.use_toci",
              numericInput("toci_day","  Day to give:",3,min=1,max=14),
              numericInput("toci_dose","  Dose (mg):",560,min=100,max=1200)
            ),
            checkboxInput("use_silt","Siltuximab (IL-6 blocker)",FALSE),
            conditionalPanel("input.use_silt",
              numericInput("silt_day","  Day to give:",3,min=1,max=14),
              numericInput("silt_dose","  Dose (mg):",770,min=100,max=1500)
            ),
            checkboxInput("use_dex","Dexamethasone",FALSE),
            conditionalPanel("input.use_dex",
              numericInput("dex_day","  Day to start:",2,min=1,max=14),
              numericInput("dex_dose","  Dose (mg/day):",10,min=1,max=40)
            ),
            checkboxInput("use_rux","Ruxolitinib (JAK1/2i)",FALSE),
            conditionalPanel("input.use_rux",
              numericInput("rux_day","  Day to start:",2,min=1,max=14),
              numericInput("rux_dose","  Dose (mg BID):",10,min=5,max=20)
            ),
            checkboxInput("use_ana","Anakinra (IL-1Ra)",FALSE),
            conditionalPanel("input.use_ana",
              numericInput("ana_day","  Day to start:",2,min=1,max=14),
              numericInput("ana_dose","  Dose (mg/day):",100,min=50,max=400)
            )
          ),
          box(title="Peak Biomarkers (Day 7)", width=4, status="danger",
            valueBoxOutput("vb_ifng",width=12),
            valueBoxOutput("vb_il6",width=12),
            valueBoxOutput("vb_ferritin",width=12)
          ),
          box(title="CRS Summary", width=4, status="danger",
            valueBoxOutput("vb_grade",width=12),
            valueBoxOutput("vb_survival",width=12),
            valueBoxOutput("vb_icans",width=12)
          )
        ),
        fluidRow(
          box(title="CRS Grade Trajectory", width=12, status="danger",
            plotlyOutput("p_crs_grade", height="300px"))
        ),
        fluidRow(
          box(title="CRS Pathophysiology", width=12,
            HTML('<pre style="font-family:monospace;font-size:12px;line-height:1.5;background:#1a1a2e;color:#e0e0e0;padding:15px;border-radius:8px;">
  CAR-T Infusion (Day 0)
         |
         v
  Antigen Recognition → T Cell Activation → IL-2, IFN-γ, GM-CSF ↑↑
         |
         v
  IFN-γ + GM-CSF → Macrophage M1 Activation (JAK-STAT1)
         |
         v
  IL-6 ↑↑↑ (key driver) + IL-1β + TNF-α
         |
         ├─→ IL-6 → Endothelial activation → NO↑ → Vasodilation → HYPOTENSION
         ├─→ IL-1β + IL-6 → Hypothalamus → PGE2 → FEVER
         ├─→ VEGF + ANG2 → Vascular leak → HYPOXIA
         └─→ BBB disruption → CNS → ICANS

  TREATMENT:
  Tocilizumab → blocks IL-6R (gp130 + sIL-6R trans-signaling)
  Siltuximab  → neutralizes free IL-6
  Dexamethasone → NF-κB inhibition → ↓IL-6, TNF-α, IFN-γ
  Ruxolitinib → JAK1/2 inhibition → ↓STAT3/STAT1 → ↓IL-6 feedback
  Anakinra    → blocks IL-1R → ↓NLRP3 inflammasome pathway
            </pre>')
          )
        )
      ),

      # =====================================================================
      # TAB 2: Drug PK
      # =====================================================================
      tabItem(tabName="pk",
        fluidRow(
          box(title="Drug Concentration-Time Profiles", width=8, status="warning",
            plotlyOutput("p_pk_all",height="380px")),
          box(title="Drug Summary", width=4, status="warning",
            tableOutput("tbl_pk_params"))
        ),
        fluidRow(
          box(title="Mechanism of Action Summary", width=12,
            tableOutput("tbl_moa"))
        )
      ),

      # =====================================================================
      # TAB 3: Cytokine Storm
      # =====================================================================
      tabItem(tabName="cyto",
        fluidRow(
          box(title="IFN-γ & IL-6 (Key CRS Drivers)", width=6, status="danger",
            plotlyOutput("p_ifng_il6", height="320px")),
          box(title="IL-1β & TNF-α", width=6, status="warning",
            plotlyOutput("p_il1b_tnfa", height="320px"))
        ),
        fluidRow(
          box(title="Macrophage Activation & GM-CSF", width=6, status="danger",
            plotlyOutput("p_mac", height="300px")),
          box(title="CAR-T Expansion & Tumor Burden", width=6, status="success",
            plotlyOutput("p_cart", height="300px"))
        ),
        fluidRow(
          box(title="IL-6 Downstream: STAT3 & Endothelial Activation", width=12, status="info",
            plotlyOutput("p_il6_downstream", height="300px"))
        )
      ),

      # =====================================================================
      # TAB 4: Clinical Endpoints
      # =====================================================================
      tabItem(tabName="endpoints",
        fluidRow(
          box(title="CRS Severity Score (0-4)", width=6, status="danger",
            plotlyOutput("p_crs_sev", height="300px")),
          box(title="Survival Probability", width=6, status="success",
            plotlyOutput("p_survival", height="300px"))
        ),
        fluidRow(
          box(title="Organ Damage & Vascular Permeability", width=6, status="warning",
            plotlyOutput("p_organ", height="300px")),
          box(title="Temperature & Fibrinogen", width=6, status="info",
            plotlyOutput("p_temp_fibrin", height="300px"))
        ),
        fluidRow(
          box(title="ICANS Score", width=12, status="danger",
            plotlyOutput("p_icans", height="280px"))
        )
      ),

      # =====================================================================
      # TAB 5: Scenario Comparison
      # =====================================================================
      tabItem(tabName="compare",
        fluidRow(
          box(title="IFN-γ — All 5 Scenarios", width=6, status="danger",
            plotlyOutput("pc_ifng", height="280px")),
          box(title="IL-6 — All 5 Scenarios", width=6, status="danger",
            plotlyOutput("pc_il6", height="280px"))
        ),
        fluidRow(
          box(title="CRS Severity — All 5 Scenarios", width=6, status="warning",
            plotlyOutput("pc_crs", height="280px")),
          box(title="Survival — All 5 Scenarios", width=6, status="success",
            plotlyOutput("pc_surv", height="280px"))
        ),
        fluidRow(
          box(title="Scenario Summary Table (Day 14)", width=12,
            tableOutput("tbl_scenario_compare"))
        )
      ),

      # =====================================================================
      # TAB 6: Biomarkers & CRS Grading
      # =====================================================================
      tabItem(tabName="biomark",
        fluidRow(
          box(title="CRS Biomarker Index", width=6, status="danger",
            plotlyOutput("p_crs_index", height="280px")),
          box(title="Ferritin (log scale)", width=6, status="warning",
            plotlyOutput("p_ferritin", height="280px"))
        ),
        fluidRow(
          box(title="CRP (mg/L)", width=6, status="info",
            plotlyOutput("p_crp", height="280px")),
          box(title="ICANS Trajectory", width=6, status="danger",
            plotlyOutput("p_icans2", height="280px"))
        ),
        fluidRow(
          box(title="ASTCT CRS Grading Criteria", width=12,
            tableOutput("tbl_crs_criteria"))
        )
      )
    )
  )
)

# ==============================================================================
# Server
# ==============================================================================

server <- function(input, output, session) {

  # ---- Reactive: run simulation ----
  sim_data <- reactive({
    run_crs_sim(
      sim_days     = input$sim_days,
      cart_type    = input$cart_type,
      tumor_burden = input$tumor_burden,
      infusion_dose= input$infusion_dose,
      use_toci     = input$use_toci,
      use_dex      = input$use_dex,
      use_rux      = input$use_rux,
      use_ana      = input$use_ana,
      use_silt     = input$use_silt,
      toci_day     = input$toci_day,
      toci_dose_mg = input$toci_dose,
      dex_day      = input$dex_day,
      dex_dose_mg  = input$dex_dose,
      rux_day      = input$rux_day,
      rux_dose_mg  = input$rux_dose,
      ana_day      = input$ana_day,
      ana_dose_mg  = input$ana_dose,
      silt_day     = input$silt_day,
      silt_dose_mg = input$silt_dose
    )
  })

  # ---- Comparison: run all 5 standard scenarios ----
  compare_data <- reactive({
    scenarios_def <- list(
      list(name="Untreated",               use_toci=F,use_dex=F,use_rux=F,use_ana=F,use_silt=F),
      list(name="Tocilizumab (Grade 2)",   use_toci=T,use_dex=F,use_rux=F,use_ana=F,use_silt=F),
      list(name="Toci + DEX (Grade 3)",    use_toci=T,use_dex=T,use_rux=F,use_ana=F,use_silt=F),
      list(name="Ruxolitinib + DEX (G4)",  use_toci=F,use_dex=T,use_rux=T,use_ana=F,use_silt=F),
      list(name="Anakinra + DEX (MAS)",    use_toci=F,use_dex=T,use_rux=F,use_ana=T,use_silt=F)
    )
    bind_rows(lapply(scenarios_def, function(sc) {
      d <- run_crs_sim(sim_days=input$sim_days, cart_type=input$cart_type,
                       tumor_burden=input$tumor_burden, infusion_dose=input$infusion_dose,
                       use_toci=sc$use_toci, use_dex=sc$use_dex,
                       use_rux=sc$use_rux, use_ana=sc$use_ana, use_silt=sc$use_silt,
                       toci_day=3, dex_day=2, rux_day=2, ana_day=2)
      d$Scenario <- sc$name; d
    }))
  })

  sce_colors <- c("Untreated"="#d32f2f","Tocilizumab (Grade 2)"="#1976d2",
                  "Toci + DEX (Grade 3)"="#388e3c","Ruxolitinib + DEX (G4)"="#f57c00",
                  "Anakinra + DEX (MAS)"="#7b1fa2")

  # ---- Value boxes ----
  output$vb_ifng <- renderValueBox({
    d <- sim_data()
    peak <- round(max(d$IFNg[d$time <= 7], na.rm=TRUE))
    valueBox(paste0(peak," pg/mL"),"Peak IFN-γ (Day 7)",icon=icon("fire"),
             color=if(peak>500)"red" else if(peak>100)"orange" else "green")
  })
  output$vb_il6 <- renderValueBox({
    d <- sim_data()
    peak <- round(max(d$IL6[d$time <= 7], na.rm=TRUE))
    valueBox(paste0(peak," pg/mL"),"Peak IL-6 (Day 7)",icon=icon("flask"),
             color=if(peak>200)"red" else if(peak>80)"orange" else "green")
  })
  output$vb_ferritin <- renderValueBox({
    d <- sim_data()
    peak <- round(max(d$FERRITIN[d$time <= 14], na.rm=TRUE))
    valueBox(paste0(peak," ng/mL"),"Peak Ferritin",icon=icon("thermometer"),
             color=if(peak>2000)"red" else if(peak>500)"orange" else "green")
  })
  output$vb_grade <- renderValueBox({
    d <- sim_data()
    g <- max(d$CRS_Grade, na.rm=TRUE)
    valueBox(paste0("Grade ",g),"Max CRS Grade",icon=icon("exclamation-triangle"),
             color=c("green","yellow","orange","red","red")[min(g+1,5)])
  })
  output$vb_survival <- renderValueBox({
    d <- sim_data()
    surv <- round(tail(d$Survival,1)*100,1)
    valueBox(paste0(surv,"%"),"28-Day Survival",icon=icon("heartbeat"),
             color=if(surv>90)"green" else if(surv>70)"yellow" else "red")
  })
  output$vb_icans <- renderValueBox({
    d <- sim_data()
    ic <- round(max(d$ICANS, na.rm=TRUE),1)
    valueBox(ic,"Max ICANS Score",icon=icon("brain"),
             color=if(ic<1)"green" else if(ic<2)"yellow" else "red")
  })

  # ---- Plots ----
  make_plotly <- function(df, x_col, y_cols, y_labels, title, ylab, colors=NULL) {
    p <- plot_ly()
    for (j in seq_along(y_cols)) {
      col <- if (!is.null(colors)) colors[j] else NULL
      p <- add_lines(p, data=df, x=~get(x_col), y=~get(y_cols[j]),
                     name=y_labels[j], line=list(color=col, width=2.5))
    }
    layout(p, title=title, xaxis=list(title="Days"), yaxis=list(title=ylab),
           legend=list(orientation="h"), hovermode="x unified")
  }

  output$p_crs_grade <- renderPlotly({
    d <- sim_data()
    plot_ly(d, x=~time, y=~CRS_SEV, type="scatter", mode="lines",
            line=list(color="#d32f2f",width=3), fill="tozeroy",
            fillcolor="rgba(211,47,47,0.15)") %>%
      layout(title="CRS Severity Over Time",
             xaxis=list(title="Days"),
             yaxis=list(title="CRS Severity (0=none, 4=life-threatening)",range=c(0,4)),
             shapes=list(
               list(type="line",x0=0,x1=max(d$time),y0=0.5,y1=0.5,
                    line=list(color="#a5d6a7",dash="dash",width=1.5)),
               list(type="line",x0=0,x1=max(d$time),y0=1.5,y1=1.5,
                    line=list(color="#fff176",dash="dash",width=1.5)),
               list(type="line",x0=0,x1=max(d$time),y0=2.5,y1=2.5,
                    line=list(color="#ffa726",dash="dash",width=1.5)),
               list(type="line",x0=0,x1=max(d$time),y0=3.5,y1=3.5,
                    line=list(color="#f44336",dash="dash",width=1.5))
             ))
  })

  output$p_pk_all <- renderPlotly({
    d <- sim_data()
    p <- plot_ly()
    if (any(d$TOCI_C1 > 0)) p <- add_lines(p, data=d, x=~time, y=~TOCI_C1,
        name="Tocilizumab (μg/mL)", line=list(color="#4caf50",width=2.5))
    if (any(d$SILT_C1 > 0)) p <- add_lines(p, data=d, x=~time, y=~SILT_C1,
        name="Siltuximab (μg/mL)", line=list(color="#26a69a",width=2.5))
    if (any(d$DEX_C1 > 0)) p <- add_lines(p, data=d, x=~time, y=~DEX_C1/1000,
        name="Dexamethasone (μg/mL)", line=list(color="#2196f3",width=2.5))
    if (any(d$RUX_C1 > 0)) p <- add_lines(p, data=d, x=~time, y=~RUX_C1/1000,
        name="Ruxolitinib (μg/mL)", line=list(color="#ff9800",width=2.5))
    if (any(d$ANA_C1 > 0)) p <- add_lines(p, data=d, x=~time, y=~ANA_C1/1000,
        name="Anakinra (μg/mL)", line=list(color="#9c27b0",width=2.5))
    layout(p, title="Drug PK Profiles",
           xaxis=list(title="Days"),
           yaxis=list(title="Concentration (μg/mL; DEX/RUX/ANA scaled)"),
           hovermode="x unified", legend=list(orientation="h"))
  })

  output$tbl_pk_params <- renderTable({
    data.frame(
      Drug=c("Tocilizumab","Siltuximab","Dexamethasone","Ruxolitinib","Anakinra"),
      Class=c("Anti-IL-6R mAb","Anti-IL-6 mAb","Corticosteroid","JAK1/2i","IL-1Ra"),
      `t½`=c("11-13d","21d","3.5h","3h","4-6h"),
      Route=c("IV","IV","IV/PO","PO BID","SC QD"),
      `Primary Target`=c("IL-6R (gp130)","Free IL-6","NF-κB/GR","JAK1/2 (pSTAT3)","IL-1R"),
      check.names=FALSE
    )
  })

  output$tbl_moa <- renderTable({
    data.frame(
      Drug=c("Tocilizumab","Siltuximab","Dexamethasone","Ruxolitinib","Anakinra"),
      `Mechanism in CRS`=c(
        "Competitive IL-6R blockade; prevents cis & trans IL-6 signaling; ↓STAT3; FDA-approved for CRS",
        "Direct IL-6 neutralization; does NOT block IL-6R; alternative when TOCI unavailable",
        "GR agonist → NF-κB ↓; suppresses IL-6, TNF-α, IFN-γ; key for ICANS (stabilizes BBB)",
        "JAK1/2 inhibition → ↓pSTAT3 + ↓pSTAT1; attenuates macrophage IL-6 feedback loop",
        "IL-1R competitive antagonist; blocks NLRP3 inflammasome downstream; preferred in HLH/MAS overlap"
      ),
      `FDA/EMA Status`=c(
        "FDA-approved CRS (2017, 2021 adult extension)",
        "Off-label CRS; FDA-approved Castleman disease",
        "Standard of care grade 2+ CRS + all ICANS grades",
        "Clinical trial / compassionate use (refractory CRS/HLH)",
        "Clinical trial (pediatric HLH/MAS); off-label CRS"
      ),
      check.names=FALSE
    )
  })

  output$p_ifng_il6 <- renderPlotly({
    d <- sim_data()
    plot_ly(d) %>%
      add_lines(x=~time, y=~IFNg, name="IFN-γ (pg/mL)",
                line=list(color="#ff7043",width=2.5)) %>%
      add_lines(x=~time, y=~IL6, name="IL-6 (pg/mL)",
                line=list(color="#9c27b0",width=2.5)) %>%
      layout(title="IFN-γ & IL-6", xaxis=list(title="Days"),
             yaxis=list(title="pg/mL"), hovermode="x unified",
             legend=list(orientation="h"))
  })

  output$p_il1b_tnfa <- renderPlotly({
    d <- sim_data()
    plot_ly(d) %>%
      add_lines(x=~time, y=~IL1B, name="IL-1β (pg/mL)",
                line=list(color="#e91e63",width=2.5)) %>%
      add_lines(x=~time, y=~TNFa, name="TNF-α (pg/mL)",
                line=list(color="#ff5722",width=2.5)) %>%
      layout(title="IL-1β & TNF-α", xaxis=list(title="Days"),
             yaxis=list(title="pg/mL"), hovermode="x unified",
             legend=list(orientation="h"))
  })

  output$p_mac <- renderPlotly({
    d <- sim_data()
    plot_ly(d) %>%
      add_lines(x=~time, y=~MAC, name="Macrophage Activation (0-1)",
                line=list(color="#c2185b",width=2.5)) %>%
      add_lines(x=~time, y=~GMCSF/max(d$GMCSF,1), name="GM-CSF (norm.)",
                line=list(color="#ff9800",width=2.5,dash="dash")) %>%
      layout(title="Macrophage & GM-CSF", xaxis=list(title="Days"),
             yaxis=list(title="Activation (0-1) / normalized"), hovermode="x unified",
             legend=list(orientation="h"))
  })

  output$p_cart <- renderPlotly({
    d <- sim_data()
    plot_ly(d) %>%
      add_lines(x=~time, y=~CART, name="CAR-T cells (cells/μL equiv.)",
                line=list(color="#43a047",width=2.5)) %>%
      add_lines(x=~time, y=~TUMOR*100, name="Tumor Burden (×100)",
                line=list(color="#795548",width=2.5,dash="dot")) %>%
      layout(title="CAR-T Expansion & Tumor Burden", xaxis=list(title="Days"),
             yaxis=list(title="Cells/μL or ×100"), hovermode="x unified",
             legend=list(orientation="h"))
  })

  output$p_il6_downstream <- renderPlotly({
    d <- sim_data()
    plot_ly(d) %>%
      add_lines(x=~time, y=~STAT3, name="pSTAT3 (rel.)",
                line=list(color="#7b1fa2",width=2.5)) %>%
      add_lines(x=~time, y=~ENDO, name="Endothelial Activation (0-1)",
                line=list(color="#00897b",width=2.5)) %>%
      add_lines(x=~time, y=~VASC, name="Vascular Permeability",
                line=list(color="#0288d1",width=2.5)) %>%
      layout(title="IL-6 Downstream Signaling", xaxis=list(title="Days"),
             yaxis=list(title="Relative units"), hovermode="x unified",
             legend=list(orientation="h"))
  })

  output$p_crs_sev <- renderPlotly({
    d <- sim_data()
    plot_ly(d, x=~time, y=~CRS_SEV, type="scatter", mode="lines",
            line=list(color="#d32f2f",width=3), fill="tozeroy",
            fillcolor="rgba(211,47,47,0.15)") %>%
      layout(title="CRS Severity (0=none, 4=life-threatening)",
             xaxis=list(title="Days"), yaxis=list(title="Score",range=c(0,4.2)))
  })

  output$p_survival <- renderPlotly({
    d <- sim_data()
    plot_ly(d, x=~time, y=~Survival*100, type="scatter", mode="lines",
            line=list(color="#4caf50",width=3), fill="tozeroy",
            fillcolor="rgba(76,175,80,0.15)") %>%
      layout(title="Survival Probability (%)",
             xaxis=list(title="Days"), yaxis=list(title="%",range=c(0,105)))
  })

  output$p_organ <- renderPlotly({
    d <- sim_data()
    plot_ly(d) %>%
      add_lines(x=~time, y=~ORGAN, name="Organ Damage (0-1)",
                line=list(color="#f44336",width=2.5)) %>%
      add_lines(x=~time, y=~VASC, name="Vascular Permeability",
                line=list(color="#0288d1",width=2.5)) %>%
      layout(title="Organ Damage & Vascular Permeability",
             xaxis=list(title="Days"), yaxis=list(title="0-1 scale"),
             hovermode="x unified", legend=list(orientation="h"))
  })

  output$p_temp_fibrin <- renderPlotly({
    d <- sim_data()
    plot_ly(d) %>%
      add_lines(x=~time, y=~Temp, name="Body Temperature (°C)",
                line=list(color="#ff7043",width=2.5)) %>%
      add_lines(x=~time, y=~Fibrinogen, name="Fibrinogen (g/L)",
                line=list(color="#1976d2",width=2.5)) %>%
      layout(title="Temperature & Fibrinogen",
             xaxis=list(title="Days"), yaxis=list(title="°C / g/L"),
             hovermode="x unified", legend=list(orientation="h"))
  })

  output$p_icans <- renderPlotly({
    d <- sim_data()
    plot_ly(d, x=~time, y=~ICANS, type="scatter", mode="lines",
            line=list(color="#3f51b5",width=3), fill="tozeroy",
            fillcolor="rgba(63,81,181,0.15)") %>%
      layout(title="ICANS Score (0=none, 4=severe)",
             xaxis=list(title="Days"), yaxis=list(title="ICANS score",range=c(0,4.5)))
  })

  # Scenario comparison plots
  output$pc_ifng <- renderPlotly({
    d <- compare_data()
    p <- plot_ly()
    for (sc in unique(d$Scenario)) {
      ds <- d[d$Scenario==sc,]
      p <- add_lines(p, data=ds, x=~time, y=~IFNg, name=sc,
                     line=list(color=sce_colors[sc],width=2))
    }
    layout(p, title="IFN-γ", xaxis=list(title="Days"),
           yaxis=list(title="pg/mL"), hovermode="x unified",
           legend=list(orientation="h"))
  })

  output$pc_il6 <- renderPlotly({
    d <- compare_data()
    p <- plot_ly()
    for (sc in unique(d$Scenario)) {
      ds <- d[d$Scenario==sc,]
      p <- add_lines(p, data=ds, x=~time, y=~IL6, name=sc,
                     line=list(color=sce_colors[sc],width=2))
    }
    layout(p, title="IL-6", xaxis=list(title="Days"),
           yaxis=list(title="pg/mL"), hovermode="x unified",
           legend=list(orientation="h"))
  })

  output$pc_crs <- renderPlotly({
    d <- compare_data()
    p <- plot_ly()
    for (sc in unique(d$Scenario)) {
      ds <- d[d$Scenario==sc,]
      p <- add_lines(p, data=ds, x=~time, y=~CRS_SEV, name=sc,
                     line=list(color=sce_colors[sc],width=2))
    }
    layout(p, title="CRS Severity", xaxis=list(title="Days"),
           yaxis=list(title="Score (0-4)",range=c(0,4.2)), hovermode="x unified",
           legend=list(orientation="h"))
  })

  output$pc_surv <- renderPlotly({
    d <- compare_data()
    p <- plot_ly()
    for (sc in unique(d$Scenario)) {
      ds <- d[d$Scenario==sc,]
      p <- add_lines(p, data=ds, x=~time, y=~Survival*100, name=sc,
                     line=list(color=sce_colors[sc],width=2))
    }
    layout(p, title="Survival Probability", xaxis=list(title="Days"),
           yaxis=list(title="%",range=c(0,105)), hovermode="x unified",
           legend=list(orientation="h"))
  })

  output$tbl_scenario_compare <- renderTable({
    d <- compare_data()
    d14 <- d[abs(d$time - 14) < 0.11,]
    d14 %>%
      group_by(Scenario) %>%
      summarise(
        `CRS Grade` = round(mean(CRS_SEV),1),
        `IFN-γ (pg/mL)` = round(mean(IFNg)),
        `IL-6 (pg/mL)` = round(mean(IL6)),
        `Ferritin (ng/mL)` = round(mean(FERRITIN)),
        `CRP (mg/L)` = round(mean(CRP)),
        `Survival (%)` = round(mean(Survival)*100,1),
        .groups="drop"
      )
  })

  output$p_crs_index <- renderPlotly({
    d <- sim_data()
    plot_ly(d, x=~time, y=~CRS_Index, type="scatter", mode="lines",
            line=list(color="#d32f2f",width=3)) %>%
      layout(title="CRS Biomarker Index (IFN-γ+IL-6+CRP+Ferritin composite)",
             xaxis=list(title="Days"), yaxis=list(title="Index (composite)"))
  })

  output$p_ferritin <- renderPlotly({
    d <- sim_data()
    plot_ly(d, x=~time, y=~FERRITIN, type="scatter", mode="lines",
            line=list(color="#ff6f00",width=3)) %>%
      layout(title="Ferritin (ng/mL)",
             xaxis=list(title="Days"),
             yaxis=list(title="ng/mL", type="log",
                        tickvals=c(100,500,1000,5000,10000,50000)))
  })

  output$p_crp <- renderPlotly({
    d <- sim_data()
    plot_ly(d, x=~time, y=~CRP, type="scatter", mode="lines",
            line=list(color="#1565c0",width=3)) %>%
      layout(title="CRP (mg/L)",
             xaxis=list(title="Days"), yaxis=list(title="mg/L"))
  })

  output$p_icans2 <- renderPlotly({
    d <- sim_data()
    plot_ly(d, x=~time, y=~ICANS, type="scatter", mode="lines",
            line=list(color="#3f51b5",width=3), fill="tozeroy",
            fillcolor="rgba(63,81,181,0.15)") %>%
      layout(title="ICANS Score",
             xaxis=list(title="Days"), yaxis=list(title="Grade (0-4)",range=c(0,4.5)))
  })

  output$tbl_crs_criteria <- renderTable({
    d <- sim_data()
    last <- d[nrow(d),]
    data.frame(
      `ASTCT Criterion`=c("Fever (≥38°C)","Hypotension","Hypoxia",
                          "CRS Grade 1","CRS Grade 2","CRS Grade 3","CRS Grade 4",
                          "ICANS","Ferritin > 500 ng/mL","IL-6 > 50 pg/mL"),
      `Simulated Value`=c(
        paste0(round(last$Temp,1),"°C"),
        paste0("VASC=",round(last$VASC,2)),
        paste0("VASC=",round(last$VASC,2)),
        "Fever only",
        "Fever + IV fluids / O2",
        "Vasopressors / high-flow O2",
        "Vent / ECMO",
        round(last$ICANS,1),
        round(last$FERRITIN),
        round(last$IL6,1)
      ),
      `Status`=c(
        if(last$Temp>=38)"PRESENT" else "absent",
        if(last$VASC>0.3)"PRESENT" else "absent",
        if(last$VASC>0.5)"PRESENT" else "absent",
        if(last$CRS_SEV>=0.5 & last$CRS_SEV<1.5)"YES" else "-",
        if(last$CRS_SEV>=1.5 & last$CRS_SEV<2.5)"YES" else "-",
        if(last$CRS_SEV>=2.5 & last$CRS_SEV<3.5)"YES" else "-",
        if(last$CRS_SEV>=3.5)"YES" else "-",
        if(last$ICANS>0.5) paste0("Grade ",floor(last$ICANS+0.5)) else "absent",
        if(last$FERRITIN>500)"ELEVATED" else "normal",
        if(last$IL6>50)"ELEVATED" else "normal"
      ),
      check.names=FALSE
    )
  })
}

shinyApp(ui, server)
