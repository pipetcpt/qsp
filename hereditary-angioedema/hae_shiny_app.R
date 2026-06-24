##############################################################################
# Hereditary Angioedema (HAE) — Interactive QSP Shiny Dashboard
# 유전성 혈관부종 인터랙티브 QSP 대시보드
#
# 6 Tabs:
#   1. Patient Profile — HAE type, C1-INH level, attack history
#   2. PK Profiles — Icatibant, Berotralstat, Lanadelumab, C1-INH conc.
#   3. KKS Biology — FXII activation, Kallikrein, Bradykinin dynamics
#   4. Clinical Endpoints — Attack severity, VP, swelling score
#   5. Scenario Comparison — Side-by-side treatment comparison
#   6. Biomarkers — C4, C1-INH %, cleaved HMWK, attack frequency
##############################################################################

library(shiny)
library(ggplot2)
library(dplyr)
library(tidyr)
library(shinydashboard)

# ─────────────────────────────────────────────────────────────────────────────
# Simulation Engine (simplified ODE solver using Euler integration)
# ─────────────────────────────────────────────────────────────────────────────

simulate_hae <- function(
    HAE_type       = 1,
    C1INH_base     = 0.35,
    attack_trigger_time = 12,
    attack_duration = 8,
    dose_ICA       = 0,   t_ICA   = 14,
    dose_C1INH_IV  = 0,   t_C1INH_IV = 14,
    dose_BER_daily = 0,
    dose_LAN       = 0,   t_LAN_Q2W = TRUE,
    dose_C1INH_SC  = 0,   SC_freq_days = 3.5,
    sim_hours      = 168,
    dt             = 0.5,
    BW             = 70
) {

  # Parameters
  k_C1INH_syn   <- 0.012;  k_C1INH_deg  <- 0.0086; k_C1INH_cons <- 0.15
  k_FXII_act    <- 0.08;   k_FXIIa_inh  <- 0.25;   k_FXIIa_deg  <- 0.10
  k_Kal_form    <- 0.45;   k_Kal_inh    <- 0.35;   k_Kal_deg    <- 0.08
  k_BK_syn      <- 0.80;   k_BK_deg_ACE <- 2.50;   BK_base      <- 0.15
  kon_BK        <- 2.20;   koff_BK      <- 0.30;   kint_B2R     <- 0.15; krecyc_B2R <- 0.04
  B2R_total     <- 1.0;    EC50_BK_VP   <- 0.80;   Emax_BK_VP   <- 4.5;  Hill_VP    <- 1.8
  k_VP_decay    <- 0.18;   VP_base      <- 1.0;    SW_threshold <- 1.5
  k_SW_form     <- 0.25;   k_SW_res     <- 0.12

  # Icatibant PK
  ka_ICA  <- 0.74; CL_ICA <- 15.5; Vc_ICA <- 29.0; Vp_ICA <- 18.0; Q_ICA <- 8.0; F_ICA <- 0.97
  Ki_ICA  <- 0.47; MW_ICA <- 1304.0

  # C1-INH IV
  CL_C1INH_IV <- 0.051; Vd_C1INH_IV <- 3.3

  # C1-INH SC
  ka_SC <- 0.025; F_SC <- 0.43; CL_SC <- 0.051; Vd_SC <- 3.3

  # Berotralstat
  ka_BER <- 0.35; F_BER <- 0.57; CL_BER <- 2.0; Vd_BER <- 268.0
  IC50_BER <- 3.7e-3; Emax_BER <- 0.92; Hill_BER <- 1.5

  # Lanadelumab
  ka_LAN <- 0.0087; F_LAN <- 0.61; CL_LAN <- 0.0139; Vc_LAN <- 6.4; Vp_LAN <- 4.8; Q_LAN <- 0.025
  KD_LAN <- 0.1e-3; Emax_LAN <- 0.93

  # Time vector
  times <- seq(0, sim_hours, by = dt)
  N <- length(times)

  # State variables
  C1INH <- rep(0, N);  FXIIa <- rep(0, N);  Kal <- rep(0, N)
  BK    <- rep(0, N);  B2Rf  <- rep(0, N);  B2Rb <- rep(0, N)
  VP    <- rep(0, N);  SW    <- rep(0, N)

  # PK states
  ICA_dep <- rep(0,N); ICA_C <- rep(0,N); ICA_P <- rep(0,N)
  C1I_IV  <- rep(0,N)
  C1I_SC_dep <- rep(0,N); C1I_SC_C <- rep(0,N)
  BER_gut <- rep(0,N); BER_C <- rep(0,N)
  LAN_dep <- rep(0,N); LAN_C <- rep(0,N); LAN_P <- rep(0,N)

  # Init
  C1INH[1] <- C1INH_base; FXIIa[1] <- 0.01; Kal[1] <- 0.05
  BK[1]    <- BK_base;    B2Rf[1]  <- B2R_total; B2Rb[1] <- 0.0
  VP[1]    <- VP_base;    SW[1]    <- 0.0

  # Dose schedule
  ica_doses <- if(dose_ICA > 0) {
    dose_nmol <- dose_ICA * 1000 / MW_ICA
    c(t_ICA) } else c()

  c1inh_iv_doses <- if(dose_C1INH_IV > 0) c(t_C1INH_IV) else c()

  c1inh_sc_doses <- if(dose_C1INH_SC > 0) {
    seq(0, sim_hours, by = SC_freq_days * 24)
  } else c()

  ber_doses <- if(dose_BER_daily > 0) seq(0, sim_hours, by=24) else c()

  lan_doses <- if(dose_LAN > 0) {
    if(t_LAN_Q2W) seq(0, sim_hours, by=336) else c(0)  # 336h = 14 days
  } else c()

  # Euler integration
  for (i in 2:N) {
    t <- times[i]
    h <- dt

    # Check dosing
    if(any(abs(times[i-1] - ica_doses) < h/2)) {
      dn <- dose_ICA * 1000 / MW_ICA  # nmol
      ICA_dep[i-1] <- ICA_dep[i-1] + dn
    }
    if(any(abs(times[i-1] - c1inh_iv_doses) < h/2)) {
      C1I_IV[i-1] <- C1I_IV[i-1] + dose_C1INH_IV
    }
    if(any(abs(times[i-1] - c1inh_sc_doses) < h/2)) {
      C1I_SC_dep[i-1] <- C1I_SC_dep[i-1] + dose_C1INH_SC
    }
    if(any(abs(times[i-1] - ber_doses) < h/2)) {
      BER_gut[i-1] <- BER_gut[i-1] + dose_BER_daily * 1000  # ug
    }
    if(any(abs(times[i-1] - lan_doses) < h/2)) {
      LAN_dep[i-1] <- LAN_dep[i-1] + dose_LAN * 1000  # ug
    }

    # Drug concentrations
    C_ICA_nM <- ICA_C[i-1] / Vc_ICA
    C_C1INH_exo <- C1I_IV[i-1] / Vd_C1INH_IV + C1I_SC_C[i-1] / Vd_SC
    C_BER <- BER_C[i-1] / Vd_BER
    C_LAN <- LAN_C[i-1] / Vc_LAN

    # Drug effects
    Kd_B2R_nM <- koff_BK / kon_BK * 1000
    E_ICA <- (C_ICA_nM/Ki_ICA) / (1 + C_ICA_nM/Ki_ICA + BK[i-1]/Kd_B2R_nM)
    E_ICA <- max(0, min(1, E_ICA))

    BER_h <- C_BER^Hill_BER; IC50_h <- IC50_BER^Hill_BER
    E_BER <- Emax_BER * BER_h / (IC50_h + BER_h + 1e-10)
    E_BER <- max(0, min(Emax_BER, E_BER))

    # Lanadelumab: KD in ug/mL ≈ 0.1e-3 * 150000/1000 = 0.015 ug/mL
    KD_LAN_ugmL <- KD_LAN * 150000 / 1000
    E_LAN <- Emax_LAN * C_LAN / (KD_LAN_ugmL + C_LAN + 1e-10)
    E_LAN <- max(0, min(Emax_LAN, E_LAN))

    C1INH_total <- C1INH[i-1] + C_C1INH_exo * 0.01

    # Attack trigger
    attack_on <- (t >= attack_trigger_time && t <= attack_trigger_time + attack_duration)
    trigger <- if(attack_on) 1.0 else 0.0

    # C1-INH
    cons_rate <- k_C1INH_cons * (FXIIa[i-1] + Kal[i-1]) * C1INH[i-1]
    dC1INH <- k_C1INH_syn - k_C1INH_deg * C1INH[i-1] - cons_rate

    # FXIIa
    FXII_trig <- trigger + 0.5 * Kal[i-1]
    dFXIIa <- k_FXII_act * FXII_trig - k_FXIIa_inh * C1INH_total * FXIIa[i-1] - k_FXIIa_deg * FXIIa[i-1]

    # Kallikrein
    dKal <- k_Kal_form * FXIIa[i-1] * (1-E_LAN) - k_Kal_inh * C1INH_total * Kal[i-1] - k_Kal_deg * Kal[i-1] * (1-E_BER) + 0.05*0.001

    # BK
    dBK <- k_BK_syn * Kal[i-1] - k_BK_deg_ACE * max(0, BK[i-1]-BK_base) + BK_base*0.01

    # B2R
    eff_bind <- kon_BK * BK[i-1] * B2Rf[i-1] * (1-E_ICA)
    dissoc   <- koff_BK * B2Rb[i-1]
    intern   <- kint_B2R * B2Rb[i-1]
    recycle  <- krecyc_B2R * (B2R_total - B2Rf[i-1] - B2Rb[i-1])
    dB2Rf    <- -eff_bind + dissoc + recycle
    dB2Rb    <-  eff_bind - dissoc - intern

    # VP
    B2R_frac <- B2Rb[i-1] / B2R_total
    VP_stim  <- Emax_BK_VP * B2R_frac^Hill_VP / (EC50_BK_VP^Hill_VP/B2R_total^Hill_VP + B2R_frac^Hill_VP + 1e-10)
    VP_target <- VP_base + VP_stim
    dVP      <- k_VP_decay * (VP_target - VP[i-1])

    # SW
    SW_drive <- if(VP[i-1] > SW_threshold) k_SW_form * (VP[i-1] - SW_threshold) else 0
    dSW <- SW_drive - k_SW_res * SW[i-1]

    # Integrate
    C1INH[i] <- max(0, C1INH[i-1] + h*dC1INH)
    FXIIa[i] <- max(0, FXIIa[i-1] + h*dFXIIa)
    Kal[i]   <- max(0, Kal[i-1]   + h*dKal)
    BK[i]    <- max(0, BK[i-1]    + h*dBK)
    B2Rf[i]  <- max(0, min(B2R_total, B2Rf[i-1] + h*dB2Rf))
    B2Rb[i]  <- max(0, B2Rb[i-1] + h*dB2Rb)
    VP[i]    <- max(0, VP[i-1]   + h*dVP)
    SW[i]    <- max(0, SW[i-1]   + h*dSW)

    # PK integration
    ICA_dep[i] <- max(0, ICA_dep[i-1] + h*(-ka_ICA*ICA_dep[i-1]))
    dICA_C <- ka_ICA*ICA_dep[i-1]*F_ICA - (CL_ICA/Vc_ICA)*ICA_C[i-1] - (Q_ICA/Vc_ICA)*ICA_C[i-1] + (Q_ICA/Vp_ICA)*ICA_P[i-1]
    dICA_P <- (Q_ICA/Vc_ICA)*ICA_C[i-1] - (Q_ICA/Vp_ICA)*ICA_P[i-1]
    ICA_C[i] <- max(0, ICA_C[i-1]+h*dICA_C)
    ICA_P[i] <- max(0, ICA_P[i-1]+h*dICA_P)

    C1I_IV[i] <- max(0, C1I_IV[i-1] + h*(-(CL_C1INH_IV/Vd_C1INH_IV)*C1I_IV[i-1]))

    C1I_SC_dep[i] <- max(0, C1I_SC_dep[i-1] + h*(-ka_SC*C1I_SC_dep[i-1]))
    dC1I_SC_C <- ka_SC*C1I_SC_dep[i-1]*F_SC - (CL_SC/Vd_SC)*C1I_SC_C[i-1]
    C1I_SC_C[i] <- max(0, C1I_SC_C[i-1]+h*dC1I_SC_C)

    BER_gut[i] <- max(0, BER_gut[i-1] + h*(-ka_BER*BER_gut[i-1]))
    dBER_C <- ka_BER*BER_gut[i-1]*F_BER - (CL_BER/Vd_BER)*BER_C[i-1]
    BER_C[i] <- max(0, BER_C[i-1]+h*dBER_C)

    LAN_dep[i] <- max(0, LAN_dep[i-1] + h*(-ka_LAN*LAN_dep[i-1]))
    dLAN_C <- ka_LAN*LAN_dep[i-1]*F_LAN - (CL_LAN/Vc_LAN)*LAN_C[i-1] - (Q_LAN/Vc_LAN)*LAN_C[i-1] + (Q_LAN/Vp_LAN)*LAN_P[i-1]
    dLAN_P <- (Q_LAN/Vc_LAN)*LAN_C[i-1] - (Q_LAN/Vp_LAN)*LAN_P[i-1]
    LAN_C[i] <- max(0, LAN_C[i-1]+h*dLAN_C)
    LAN_P[i] <- max(0, LAN_P[i-1]+h*dLAN_P)
  }

  data.frame(
    time      = times,
    C1INH_pct = C1INH * 100,
    C4_proxy  = 100 * exp(-1.5 * FXIIa),
    FXIIa     = FXIIa,
    Kal_act   = Kal,
    BK        = BK,
    BK_fold   = BK / BK_base,
    B2R_occ   = B2Rb / B2R_total * 100,
    VP        = VP,
    SW_score  = SW,
    ICA_nM    = ICA_C / Vc_ICA,
    C1INH_IV_ugmL = C1I_IV / Vd_C1INH_IV * 50,  # rough scaling to ug/mL
    C1INH_SC_ugmL = C1I_SC_C / Vd_SC * 50,
    BER_ugmL  = BER_C / Vd_BER,
    LAN_ugmL  = LAN_C / Vc_LAN
  )
}

# ─────────────────────────────────────────────────────────────────────────────
# UI
# ─────────────────────────────────────────────────────────────────────────────

ui <- dashboardPage(
  skin = "purple",

  dashboardHeader(
    title = span("HAE QSP Dashboard", style="font-size:15px; font-weight:bold"),
    titleWidth = 280
  ),

  dashboardSidebar(
    width = 280,
    sidebarMenu(
      menuItem("① Patient Profile",      tabName="tab_patient",   icon=icon("user")),
      menuItem("② PK Profiles",          tabName="tab_pk",        icon=icon("pills")),
      menuItem("③ KKS Biology",          tabName="tab_biology",   icon=icon("dna")),
      menuItem("④ Clinical Endpoints",   tabName="tab_clinical",  icon=icon("chart-line")),
      menuItem("⑤ Scenario Comparison",  tabName="tab_scenarios", icon=icon("balance-scale")),
      menuItem("⑥ Biomarkers",           tabName="tab_biomarkers",icon=icon("flask"))
    ),
    hr(),
    h5("Global Parameters", style="margin-left:15px;color:#ddd"),

    selectInput("hae_type",  "HAE Type",
                choices = c("Type I (C1INH deficient)"=1,
                            "Type II (C1INH dysfunctional)"=2,
                            "Type III (FXII gain-of-function)"=3),
                selected = 1),

    sliderInput("c1inh_base", "Baseline C1-INH (%)",
                min=5, max=100, value=35, step=5,
                post="%"),

    sliderInput("attack_time", "Attack onset (h)",
                min=0, max=72, value=12, step=2),

    sliderInput("attack_dur", "Attack duration (h)",
                min=2, max=24, value=8, step=1),

    sliderInput("bw", "Body weight (kg)",
                min=40, max=120, value=70, step=5),

    sliderInput("sim_hours", "Simulation duration (h)",
                min=24, max=720, value=168, step=24),

    hr(),
    h5("Acute Treatment", style="margin-left:15px;color:#ddd"),
    checkboxInput("use_ICA",   "Icatibant 30mg SC", value=FALSE),
    conditionalPanel("input.use_ICA",
      sliderInput("t_ICA", "Give at (h post-attack):", min=0, max=24, value=2, step=0.5)),

    checkboxInput("use_C1INH_IV", "C1-INH IV (Berinert 20 IU/kg)", value=FALSE),
    conditionalPanel("input.use_C1INH_IV",
      sliderInput("t_C1INH_IV", "Give at (h):", min=0, max=24, value=2, step=0.5)),

    hr(),
    h5("Prophylaxis", style="margin-left:15px;color:#ddd"),
    checkboxInput("use_BER",   "Berotralstat 150mg QD", value=FALSE),
    checkboxInput("use_LAN",   "Lanadelumab 300mg Q2W", value=FALSE),
    checkboxInput("use_SC",    "C1-INH SC (Haegarda 60 IU/kg)", value=FALSE),

    actionButton("run_sim", "Run Simulation",
                 icon=icon("play"), class="btn-primary btn-block",
                 style="margin:10px 15px; width:250px")
  ),

  dashboardBody(
    tags$head(
      tags$style(HTML("
        .content-wrapper { background-color: #f8f9fa; }
        .box { border-radius: 8px; }
        .info-box { border-radius: 8px; }
      "))
    ),

    tabItems(

      # ══════════════════════════════════════════════════════════════
      # TAB 1: Patient Profile
      # ══════════════════════════════════════════════════════════════
      tabItem(tabName = "tab_patient",
        fluidRow(
          box(title="HAE Pathophysiology Summary", status="primary", solidHeader=TRUE, width=12,
            fluidRow(
              infoBox("HAE Type", uiOutput("hae_type_label"), icon=icon("dna"),
                      color="purple", width=4),
              infoBox("Baseline C1-INH", uiOutput("c1inh_label"),
                      icon=icon("vial"), color="red", width=4),
              infoBox("Simulation Duration", uiOutput("sim_label"),
                      icon=icon("clock"), color="blue", width=4)
            )
          )
        ),
        fluidRow(
          box(title="HAE Diagnostic Criteria", status="info", solidHeader=TRUE, width=6,
            tableOutput("dx_table")
          ),
          box(title="Treatment Options Overview", status="success", solidHeader=TRUE, width=6,
            tableOutput("tx_table")
          )
        ),
        fluidRow(
          box(title="HAE Attack Pattern (Simulated)", status="warning", solidHeader=TRUE, width=12,
            plotOutput("attack_pattern_plot", height="250px")
          )
        )
      ),

      # ══════════════════════════════════════════════════════════════
      # TAB 2: PK Profiles
      # ══════════════════════════════════════════════════════════════
      tabItem(tabName = "tab_pk",
        fluidRow(
          box(title="Icatibant Plasma Concentration (nM)", status="primary", solidHeader=TRUE, width=6,
            plotOutput("pk_ica_plot", height="280px")
          ),
          box(title="Berotralstat Concentration (ug/mL)", status="success", solidHeader=TRUE, width=6,
            plotOutput("pk_ber_plot", height="280px")
          )
        ),
        fluidRow(
          box(title="Lanadelumab Concentration (ug/mL)", status="warning", solidHeader=TRUE, width=6,
            plotOutput("pk_lan_plot", height="280px")
          ),
          box(title="C1-INH Exogenous Level (a.u.)", status="info", solidHeader=TRUE, width=6,
            plotOutput("pk_c1inh_plot", height="280px")
          )
        )
      ),

      # ══════════════════════════════════════════════════════════════
      # TAB 3: KKS Biology
      # ══════════════════════════════════════════════════════════════
      tabItem(tabName = "tab_biology",
        fluidRow(
          box(title="FXIIa Activation Dynamics", status="danger", solidHeader=TRUE, width=6,
            plotOutput("fxiia_plot", height="280px")
          ),
          box(title="Plasma Kallikrein Activity", status="warning", solidHeader=TRUE, width=6,
            plotOutput("kal_plot", height="280px")
          )
        ),
        fluidRow(
          box(title="Bradykinin Plasma Level", status="primary", solidHeader=TRUE, width=6,
            plotOutput("bk_plot", height="280px")
          ),
          box(title="B2R Occupancy (%)", status="info", solidHeader=TRUE, width=6,
            plotOutput("b2r_plot", height="280px")
          )
        )
      ),

      # ══════════════════════════════════════════════════════════════
      # TAB 4: Clinical Endpoints
      # ══════════════════════════════════════════════════════════════
      tabItem(tabName = "tab_clinical",
        fluidRow(
          infoBoxOutput("max_bk_box",    width=3),
          infoBoxOutput("max_sw_box",    width=3),
          infoBoxOutput("max_vp_box",    width=3),
          infoBoxOutput("attack_dur_box",width=3)
        ),
        fluidRow(
          box(title="Vascular Permeability Index", status="danger", solidHeader=TRUE, width=6,
            plotOutput("vp_plot", height="280px")
          ),
          box(title="Swelling Score (0-10)", status="warning", solidHeader=TRUE, width=6,
            plotOutput("sw_plot", height="280px")
          )
        ),
        fluidRow(
          box(title="C1-INH Level (% of normal)", status="success", solidHeader=TRUE, width=6,
            plotOutput("c1inh_plot", height="280px")
          ),
          box(title="Combined Clinical Cascade", status="primary", solidHeader=TRUE, width=6,
            plotOutput("cascade_plot", height="280px")
          )
        )
      ),

      # ══════════════════════════════════════════════════════════════
      # TAB 5: Scenario Comparison
      # ══════════════════════════════════════════════════════════════
      tabItem(tabName = "tab_scenarios",
        fluidRow(
          box(title="Run All Treatment Scenarios", status="primary", solidHeader=TRUE, width=12,
            fluidRow(
              column(3, actionButton("run_all_scenarios", "Compare All 6 Scenarios",
                                     icon=icon("chart-bar"), class="btn-warning btn-block")),
              column(9, p("Compares: (1) Untreated, (2) Icatibant, (3) C1-INH IV,
                           (4) Berotralstat QD, (5) Lanadelumab Q2W, (6) C1-INH SC",
                          style="padding-top:8px"))
            )
          )
        ),
        fluidRow(
          box(title="BK Level — All Scenarios", status="danger", solidHeader=TRUE, width=6,
            plotOutput("comp_bk_plot", height="280px")
          ),
          box(title="Swelling Score — All Scenarios", status="warning", solidHeader=TRUE, width=6,
            plotOutput("comp_sw_plot", height="280px")
          )
        ),
        fluidRow(
          box(title="Scenario Summary Table", status="info", solidHeader=TRUE, width=12,
            tableOutput("scenario_table")
          )
        )
      ),

      # ══════════════════════════════════════════════════════════════
      # TAB 6: Biomarkers
      # ══════════════════════════════════════════════════════════════
      tabItem(tabName = "tab_biomarkers",
        fluidRow(
          box(title="C4 Level Proxy (%)", status="info", solidHeader=TRUE, width=6,
            plotOutput("c4_plot", height="280px")
          ),
          box(title="C1-INH % of Normal", status="success", solidHeader=TRUE, width=6,
            plotOutput("c1inh_bio_plot", height="280px")
          )
        ),
        fluidRow(
          box(title="Drug Effect Summary — Kallikrein Inhibition", status="warning", solidHeader=TRUE, width=6,
            plotOutput("kal_inhib_plot", height="280px")
          ),
          box(title="Biomarker Reference Table", status="primary", solidHeader=TRUE, width=6,
            tableOutput("biomarker_table")
          )
        )
      )
    )
  )
)

# ─────────────────────────────────────────────────────────────────────────────
# SERVER
# ─────────────────────────────────────────────────────────────────────────────

server <- function(input, output, session) {

  # Reactive simulation result
  sim_result <- eventReactive(input$run_sim, {
    dose_ica  <- if(input$use_ICA)        30  else 0  # mg
    dose_civ  <- if(input$use_C1INH_IV)  20 * input$bw else 0  # IU
    dose_ber  <- if(input$use_BER)        150 else 0  # mg
    dose_lan  <- if(input$use_LAN)        300 else 0  # mg
    dose_sc   <- if(input$use_SC)         60 * input$bw else 0  # IU

    simulate_hae(
      HAE_type       = as.integer(input$hae_type),
      C1INH_base     = input$c1inh_base / 100,
      attack_trigger_time = input$attack_time,
      attack_duration     = input$attack_dur,
      dose_ICA       = dose_ica,
      t_ICA          = input$attack_time + ifelse(input$use_ICA, input$t_ICA, 0),
      dose_C1INH_IV  = dose_civ,
      t_C1INH_IV     = input$attack_time + ifelse(input$use_C1INH_IV, input$t_C1INH_IV, 0),
      dose_BER_daily = dose_ber,
      dose_LAN       = dose_lan,
      dose_C1INH_SC  = dose_sc,
      sim_hours      = input$sim_hours,
      BW             = input$bw
    )
  }, ignoreNULL = FALSE)

  # All scenarios (comparison)
  all_scenarios <- eventReactive(input$run_all_scenarios, {
    withProgress(message="Running scenarios...", value=0, {
      sim_h <- min(input$sim_hours, 168)

      s1 <- simulate_hae(C1INH_base=0.35, attack_trigger_time=12, attack_duration=8,
                          sim_hours=sim_h); s1$scenario <- "S1: Untreated"

      s2 <- simulate_hae(C1INH_base=0.35, attack_trigger_time=12, attack_duration=8,
                          dose_ICA=30, t_ICA=14, sim_hours=sim_h)
      s2$scenario <- "S2: Icatibant 30mg SC"

      s3 <- simulate_hae(C1INH_base=0.35, attack_trigger_time=12, attack_duration=8,
                          dose_C1INH_IV=1400, t_C1INH_IV=14, sim_hours=sim_h)
      s3$scenario <- "S3: C1-INH IV (Berinert)"

      s4 <- simulate_hae(C1INH_base=0.35, attack_trigger_time=12, attack_duration=8,
                          dose_BER_daily=150, sim_hours=sim_h)
      s4$scenario <- "S4: Berotralstat QD"

      s5 <- simulate_hae(C1INH_base=0.35, attack_trigger_time=12, attack_duration=8,
                          dose_LAN=300, t_LAN_Q2W=FALSE, sim_hours=sim_h)
      s5$scenario <- "S5: Lanadelumab 300mg"

      s6 <- simulate_hae(C1INH_base=0.35, attack_trigger_time=12, attack_duration=8,
                          dose_C1INH_SC=4200, SC_freq_days=3.5, sim_hours=sim_h)
      s6$scenario <- "S6: C1-INH SC (Haegarda)"

      bind_rows(s1, s2, s3, s4, s5, s6)
    })
  })

  # ─────────────────────────────────────────────────────
  # TAB 1: Info boxes and static tables
  # ─────────────────────────────────────────────────────
  output$hae_type_label <- renderUI({
    types <- c("1"="Type I (C1INH Low)", "2"="Type II (C1INH Dysfunc.)", "3"="Type III (FXII Mut)")
    h4(types[input$hae_type], style="color:white")
  })
  output$c1inh_label <- renderUI({
    h4(paste0(input$c1inh_base, "% of normal"), style="color:white")
  })
  output$sim_label <- renderUI({
    h4(paste0(input$sim_hours, " hours"), style="color:white")
  })

  output$dx_table <- renderTable({
    data.frame(
      `HAE Type` = c("Type I","Type II","Type III"),
      `C1-INH Antigen` = c("<50%","Normal",">Normal"),
      `C1-INH Function` = c("Low","Low","Normal"),
      `C4 Level` = c("Low","Low","Normal"),
      `C1q Level` = c("Normal","Normal","Normal"),
      stringsAsFactors = FALSE
    )
  }, striped=TRUE, hover=TRUE, bordered=TRUE)

  output$tx_table <- renderTable({
    data.frame(
      Drug = c("Icatibant","C1-INH IV","Ecallantide","Berotralstat","Lanadelumab","C1-INH SC"),
      Mechanism = c("B2R antagonist","C1-INH replacement","Kal. inhibitor",
                    "Oral Kal. inhib.","Anti-preKal mAb","C1-INH replacement"),
      Use = c("Acute","Acute","Acute","Prophylaxis","Prophylaxis","Prophylaxis"),
      Efficacy = c("~4h relief","~1-2h relief","~4h relief","-44% attacks","-87% attacks","-95% attacks"),
      stringsAsFactors = FALSE
    )
  }, striped=TRUE, hover=TRUE, bordered=TRUE)

  output$attack_pattern_plot <- renderPlot({
    d <- sim_result()
    ggplot(d, aes(x=time)) +
      geom_line(aes(y=SW_score, color="Swelling Score"), linewidth=1.2) +
      geom_line(aes(y=BK_fold*0.5, color="BK Fold Change (scaled)"), linewidth=1, linetype="dashed") +
      geom_vline(xintercept=input$attack_time, linetype="dotted", color="red", linewidth=0.8) +
      annotate("text", x=input$attack_time+1, y=max(d$SW_score)*0.9,
               label="Attack\nonset", hjust=0, size=3.5, color="red") +
      scale_color_manual(values=c("Swelling Score"="#c62828","BK Fold Change (scaled)"="#1976D2")) +
      labs(title="Attack Timeline", x="Time (h)", y="Score / fold-change", color="") +
      theme_bw(base_size=11) + theme(legend.position="bottom")
  })

  # ─────────────────────────────────────────────────────
  # TAB 2: PK Profiles
  # ─────────────────────────────────────────────────────
  output$pk_ica_plot <- renderPlot({
    d <- sim_result()
    ggplot(d, aes(x=time, y=ICA_nM)) +
      geom_line(color="#1565C0", linewidth=1.3) +
      geom_hline(yintercept=0.47, linetype="dashed", color="orange") +
      annotate("text", x=5, y=0.6, label="Ki(B2R) = 0.47 nM", size=3.5, color="orange") +
      labs(x="Time (h)", y="Icatibant (nM)") + theme_bw(base_size=11)
  })

  output$pk_ber_plot <- renderPlot({
    d <- sim_result()
    ggplot(d, aes(x=time, y=BER_ugmL)) +
      geom_line(color="#2E7D32", linewidth=1.3) +
      geom_hline(yintercept=3.7e-3, linetype="dashed", color="orange") +
      annotate("text", x=10, y=4e-3, label="IC50 = 3.7e-3 ug/mL", size=3.5, color="orange") +
      labs(x="Time (h)", y="Berotralstat (ug/mL)") + theme_bw(base_size=11)
  })

  output$pk_lan_plot <- renderPlot({
    d <- sim_result()
    ggplot(d, aes(x=time, y=LAN_ugmL)) +
      geom_line(color="#1B5E20", linewidth=1.3) +
      labs(x="Time (h)", y="Lanadelumab (ug/mL)") + theme_bw(base_size=11)
  })

  output$pk_c1inh_plot <- renderPlot({
    d <- sim_result()
    ggplot(d, aes(x=time, y=C1INH_IV_ugmL + C1INH_SC_ugmL)) +
      geom_line(aes(y=C1INH_IV_ugmL, color="IV C1-INH"), linewidth=1.2) +
      geom_line(aes(y=C1INH_SC_ugmL, color="SC C1-INH"), linewidth=1.2) +
      scale_color_manual(values=c("IV C1-INH"="#0D47A1","SC C1-INH"="#1B5E20")) +
      labs(x="Time (h)", y="C1-INH (a.u.)", color="") + theme_bw(base_size=11) +
      theme(legend.position="bottom")
  })

  # ─────────────────────────────────────────────────────
  # TAB 3: KKS Biology
  # ─────────────────────────────────────────────────────
  output$fxiia_plot <- renderPlot({
    d <- sim_result()
    ggplot(d, aes(x=time, y=FXIIa)) +
      geom_line(color="#E65100", linewidth=1.3) +
      geom_vline(xintercept=input$attack_time, linetype="dotted", color="red") +
      labs(x="Time (h)", y="FXIIa (normalized)") + theme_bw(base_size=11)
  })

  output$kal_plot <- renderPlot({
    d <- sim_result()
    ggplot(d, aes(x=time, y=Kal_act)) +
      geom_line(color="#E91E63", linewidth=1.3) +
      labs(x="Time (h)", y="Kallikrein (normalized)") + theme_bw(base_size=11)
  })

  output$bk_plot <- renderPlot({
    d <- sim_result()
    ggplot(d, aes(x=time, y=BK)) +
      geom_line(color="#880E4F", linewidth=1.3) +
      geom_hline(yintercept=0.15, linetype="dashed", color="gray", alpha=0.7) +
      annotate("text", x=5, y=0.20, label="Baseline BK", size=3.5, color="gray40") +
      labs(x="Time (h)", y="Bradykinin (ng/mL)") + theme_bw(base_size=11)
  })

  output$b2r_plot <- renderPlot({
    d <- sim_result()
    ggplot(d, aes(x=time, y=B2R_occ)) +
      geom_line(color="#9C27B0", linewidth=1.3) +
      geom_hline(yintercept=50, linetype="dashed", color="orange") +
      annotate("text", x=5, y=52, label="50% occupancy", size=3.5, color="orange") +
      labs(x="Time (h)", y="B2R Occupancy (%)") +
      ylim(0,100) + theme_bw(base_size=11)
  })

  # ─────────────────────────────────────────────────────
  # TAB 4: Clinical Endpoints
  # ─────────────────────────────────────────────────────
  output$max_bk_box <- renderInfoBox({
    d <- sim_result()
    infoBox("Max BK Fold", sprintf("%.1fx", max(d$BK_fold, na.rm=TRUE)),
            icon=icon("arrow-up"), color="red")
  })
  output$max_sw_box <- renderInfoBox({
    d <- sim_result()
    infoBox("Max Swelling", sprintf("%.1f/10", max(d$SW_score, na.rm=TRUE)),
            icon=icon("person"), color="orange")
  })
  output$max_vp_box <- renderInfoBox({
    d <- sim_result()
    infoBox("Max VP Index", sprintf("%.2f", max(d$VP, na.rm=TRUE)),
            icon=icon("wave-square"), color="purple")
  })
  output$attack_dur_box <- renderInfoBox({
    d <- sim_result()
    # Duration above threshold
    above <- d %>% filter(SW_score > 0.5)
    dur_h <- if(nrow(above)>0) diff(range(above$time)) else 0
    infoBox("Attack Duration", sprintf("%.0f h", dur_h),
            icon=icon("clock"), color="blue")
  })

  output$vp_plot <- renderPlot({
    d <- sim_result()
    ggplot(d, aes(x=time, y=VP)) +
      geom_line(color="#00695C", linewidth=1.3) +
      geom_hline(yintercept=1.5, linetype="dashed", color="red") +
      annotate("text", x=5, y=1.55, label="Edema threshold", size=3.5, color="red") +
      geom_ribbon(data=d %>% filter(VP>1.5), aes(ymin=1.5, ymax=VP), fill="red", alpha=0.2) +
      labs(x="Time (h)", y="VP Index") + theme_bw(base_size=11)
  })

  output$sw_plot <- renderPlot({
    d <- sim_result()
    ggplot(d, aes(x=time, y=SW_score)) +
      geom_line(color="#c62828", linewidth=1.3) +
      geom_area(fill="#c62828", alpha=0.15) +
      labs(x="Time (h)", y="Swelling Score (0-10)") +
      ylim(0, max(10, max(d$SW_score, na.rm=TRUE)*1.1)) +
      theme_bw(base_size=11)
  })

  output$c1inh_plot <- renderPlot({
    d <- sim_result()
    ggplot(d, aes(x=time, y=C1INH_pct)) +
      geom_line(color="#2E7D32", linewidth=1.3) +
      geom_hline(yintercept=50, linetype="dashed", color="red") +
      annotate("text", x=5, y=52, label="Diagnostic threshold (50%)", size=3.5, color="red") +
      labs(x="Time (h)", y="C1-INH (% of normal)") +
      ylim(0, 110) + theme_bw(base_size=11)
  })

  output$cascade_plot <- renderPlot({
    d <- sim_result() %>%
      select(time, BK_fold, B2R_occ, VP, SW_score) %>%
      mutate(VP_scaled = VP/5*100, B2R_occ=B2R_occ) %>%
      gather(key="variable", value="value",
             BK_fold, VP_scaled, SW_score) %>%
      mutate(variable = factor(variable,
               levels=c("BK_fold","VP_scaled","SW_score"),
               labels=c("BK Fold Change","VP (scaled %)","Swelling Score")))

    ggplot(d, aes(x=time, y=value, color=variable)) +
      geom_line(linewidth=1.2) +
      scale_color_manual(values=c("#880E4F","#00695C","#c62828")) +
      labs(title="KKS → VP → Edema Cascade", x="Time (h)", y="Value", color="") +
      theme_bw(base_size=11) + theme(legend.position="bottom")
  })

  # ─────────────────────────────────────────────────────
  # TAB 5: Scenario Comparison
  # ─────────────────────────────────────────────────────
  output$comp_bk_plot <- renderPlot({
    req(all_scenarios())
    d <- all_scenarios()
    ggplot(d, aes(x=time, y=BK_fold, color=scenario)) +
      geom_line(linewidth=1.1) +
      scale_color_brewer(palette="Dark2") +
      labs(x="Time (h)", y="BK Fold Change", color="") +
      theme_bw(base_size=11) + theme(legend.position="bottom",
                                      legend.text=element_text(size=8))
  })

  output$comp_sw_plot <- renderPlot({
    req(all_scenarios())
    d <- all_scenarios()
    ggplot(d, aes(x=time, y=SW_score, color=scenario)) +
      geom_line(linewidth=1.1) +
      scale_color_brewer(palette="Dark2") +
      labs(x="Time (h)", y="Swelling Score", color="") +
      theme_bw(base_size=11) + theme(legend.position="bottom",
                                      legend.text=element_text(size=8))
  })

  output$scenario_table <- renderTable({
    req(all_scenarios())
    d <- all_scenarios()
    d %>%
      group_by(scenario) %>%
      summarise(
        `Max BK Fold` = sprintf("%.2f", max(BK_fold, na.rm=TRUE)),
        `Max SW Score` = sprintf("%.2f", max(SW_score, na.rm=TRUE)),
        `Max VP Index` = sprintf("%.2f", max(VP, na.rm=TRUE)),
        `Min C1INH %`  = sprintf("%.1f%%", min(C1INH_pct, na.rm=TRUE)),
        `AUC(BK)` = sprintf("%.1f", max(cumsum(BK_fold)*0.5, na.rm=TRUE)),
        .groups="drop"
      ) %>%
      rename(Scenario=scenario)
  }, striped=TRUE, hover=TRUE, bordered=TRUE)

  # ─────────────────────────────────────────────────────
  # TAB 6: Biomarkers
  # ─────────────────────────────────────────────────────
  output$c4_plot <- renderPlot({
    d <- sim_result()
    ggplot(d, aes(x=time, y=C4_proxy)) +
      geom_line(color="#1565C0", linewidth=1.3) +
      geom_hline(yintercept=80, linetype="dashed", color="orange") +
      annotate("text", x=5, y=82, label="Lower normal", size=3.5, color="orange") +
      labs(x="Time (h)", y="C4 Level (% proxy)") +
      ylim(0, 110) + theme_bw(base_size=11)
  })

  output$c1inh_bio_plot <- renderPlot({
    d <- sim_result()
    ggplot(d, aes(x=time, y=C1INH_pct)) +
      geom_line(color="#4CAF50", linewidth=1.3) +
      geom_hline(yintercept=50, linetype="dashed", color="red") +
      annotate("text", x=5, y=52, label="HAE I diagnostic (<50%)", size=3.5, color="red") +
      labs(x="Time (h)", y="C1-INH (% of normal)") +
      ylim(0, 130) + theme_bw(base_size=11)
  })

  output$kal_inhib_plot <- renderPlot({
    d <- sim_result()
    ggplot(d, aes(x=time)) +
      geom_line(aes(y=(1-Kal_act/max(Kal_act+0.001))*100, color="Kallikrein suppression (%)"), linewidth=1.2) +
      scale_color_manual(values=c("Kallikrein suppression (%)"="#E91E63")) +
      labs(x="Time (h)", y="Kallikrein Suppression (%)", color="") +
      ylim(0,100) + theme_bw(base_size=11) + theme(legend.position="bottom")
  })

  output$biomarker_table <- renderTable({
    data.frame(
      Biomarker  = c("C4", "C1-INH antigen", "C1-INH function",
                     "C1q", "Cleaved HMWK", "Plasma BK",
                     "Kallikrein activity"),
      `Normal value` = c("90-140 mg/L", "150-350 mg/L", ">67%",
                          "50-250 mg/L", "<20% cleaved", "<1 ng/mL",
                          "Low/undetectable"),
      `HAE (interattack)` = c("<90 mg/L","<50%","<50%",
                               "Normal","Elevated","Slightly elevated","Elevated"),
      `HAE (attack)` = c("Very low","<50%","<50%",
                          "Normal",">50% cleaved",">2-10x baseline","High"),
      stringsAsFactors = FALSE
    )
  }, striped=TRUE, hover=TRUE, bordered=TRUE)

}

# ─────────────────────────────────────────────────────────────────────────────
# Run
# ─────────────────────────────────────────────────────────────────────────────
shinyApp(ui = ui, server = server)
