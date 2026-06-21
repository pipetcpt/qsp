## =============================================================================
## CIDP Interactive QSP Shiny Dashboard
## 6 Tabs: Patient Profile · Drug PK · Immune Pathways ·
##         Nerve Pathology & Electrophysiology · Scenario Comparison · Biomarkers
## =============================================================================

library(shiny)
library(shinydashboard)
library(plotly)
library(DT)
library(dplyr)
library(tidyr)
library(deSolve)

# ─────────────────────────────────────────────────────────────────
# ODE system (deSolve; no mrgsolve dependency for Shiny portability)
# ─────────────────────────────────────────────────────────────────
cidp_ode <- function(time, state, parms) {
  with(as.list(c(state, parms)), {

    # ── Drug effects ─────────────────────────────────────────────
    # FcRn-mediated catabolism multiplier (elevated at high IgG)
    FcRn_CL_mod <- 1 + (FcRn_CL_mult - 1) * IVIG_C1 / (IVIG_C1 + FcRn_Km)

    # IVIG mechanism effects
    IVIG_eff_Ab    <- if ((IVIG_FLAG + SCIG_FLAG) > 0) IVIG_C1 / (IVIG_C1 + 3) else 0
    IVIG_treg_stim <- 0.3 * IVIG_eff_Ab
    IVIG_comp_sca  <- 0.5 * IVIG_eff_Ab

    # CS effects
    CS_effect   <- CS_PLASMA / (CS_PLASMA + 0.4)
    CS_Th_inh   <- if (CS_FLAG > 0) CS_effect * 0.70 else 0
    CS_Mac_inh  <- if (CS_FLAG > 0) CS_effect * 0.50 else 0

    # PLEX
    PLEX_IgG_rm <- if (PLEX_FLAG > 0) 0.6 * PLEX_COUP else 0

    # RTX B-cell depletion
    RTX_Bdepl    <- RTX_CD20 / (RTX_CD20 + 10)
    B_RTX_factor <- if (RTX_FLAG > 0) max(0.05, 1 - RTX_Bdepl) else 1

    # EFC FcRn blockade
    EFC_FcRn_inh <- if (EFC_FLAG > 0) EFC_C / (EFC_C + EC50_FcRn) else 0
    EFC_IgG_accCL <- 1 + 4 * EFC_FcRn_inh

    # ── PK ODEs ──────────────────────────────────────────────────
    IgG_CL_eff <- CL_IgG * FcRn_CL_mod * EFC_IgG_accCL
    IgG_dist   <- Q_IgG * (IVIG_C1 - IVIG_C2)

    dIVIG_C1   <- -IgG_CL_eff * IVIG_C1 / Vd_IgG1 - IgG_dist - PLEX_IgG_rm * IVIG_C1
    dIVIG_C2   <- IgG_dist - Q_IgG * (IVIG_C2 - IVIG_C1) / Vd_IgG2

    dCS_GUT    <- -ka_CS * CS_GUT
    dCS_PLASMA <- ka_CS * CS_GUT * F_CS / Vd_CS - CL_CS / 24 * CS_PLASMA

    free_CD20  <- max(0, CD20_tot - RTX_CD20)
    dRTX_C    <- -(CL_RTX / Vd_RTX) * RTX_C - kon_RTX * RTX_C * free_CD20 + koff_RTX * RTX_CD20
    dRTX_CD20 <- kon_RTX * RTX_C * free_CD20 - koff_RTX * RTX_CD20 - 0.01 * RTX_CD20

    dEFC_C    <- -(CL_EFC / Vd_EFC) * EFC_C
    dPLEX_COUP <- -0.3 * PLEX_COUP

    # ── Immune cell ODEs ─────────────────────────────────────────
    Th1_stim  <- 1 + 0.5 * (Bc / 1)
    Th1_inhib <- (1 + Treg) * (1 + CS_Th_inh)
    dTh1  <- kprod_Th1  * Th1_stim  / Th1_inhib  - kdeg_Th1  * Th1

    Th17_stim  <- 1 + 0.3 * Mac
    Th17_inhib <- (1 + 0.8 * Treg) * (1 + CS_Th_inh)
    dTh17 <- kprod_Th17 * Th17_stim / Th17_inhib - kdeg_Th17 * Th17

    dTreg <- kprod_Treg * (1 + IVIG_treg_stim) - kdeg_Treg * Treg - 0.1 * Th17 * Treg

    Bc_prod_mod <- B_RTX_factor * (1 + 0.2 * Th1)
    dBc   <- kprod_Bc * Bc_prod_mod - kdeg_Bc * Bc
    dPC   <- kprod_PC * Bc - kdeg_PC * PC

    Mac_activ <- 0.3 * Comp + 0.2 * Th1
    dMac  <- kprod_Mac * Mac_activ * (1 - CS_Mac_inh) - kdeg_Mac * Mac * (1 + 0.5 * Treg)

    Ab_Comp_drive <- Ab_path * Mac * 0.5
    dComp <- ksynth_Comp * Ab_Comp_drive - kdeg_Comp * Comp -
             IVIG_comp_sca * Comp - PLEX_Comp_rm * Comp * 0.4

    # ── Antibody ODE ─────────────────────────────────────────────
    Ab_synth <- ksynth_Ab * PC
    Ab_deg   <- kdeg_Ab * Ab_path * FcRn_CL_mod * EFC_IgG_accCL
    dAb_path <- Ab_synth - Ab_deg - PLEX_IgG_rm * Ab_path * 0.8 - 0.02 * IVIG_eff_Ab * Ab_path

    # ── Nerve pathology ODEs ─────────────────────────────────────
    inflam_drive <- Th1 * Th17 * CIDP_SEVERITY
    dem_rate <- kdem_Ab * Ab_path * (1 + NODAL_SUBTYPE * 0.3) +
                kdem_Mac * Mac + kdem_Comp * Comp
    remy_rate <- kremy * (1 - Demyelin) * Remyel * max(0, 1 - 0.5 * inflam_drive)
    dDemyelin  <- dem_rate * (1 - Demyelin) - remy_rate
    dDemyelin  <- max(-Demyelin, min(dDemyelin, 1 - Demyelin))

    dRemyel   <- 0.02 * (1 - Remyel) - 0.01 * Mac * Remyel + 0.005 * Treg

    axon_loss  <- kaxon_dem * Demyelin + kaxon_base
    axon_regen <- kregen_axon * (1 - Axon_dens) * NCV_norm
    dAxon_dens <- axon_regen - axon_loss * Axon_dens

    dNfL <- kNfL_prod * (1 - Axon_dens) * (1 + Demyelin) - kNfL_clear * NfL + 1

    NCV_loss_rate <- kNCV_dem * Demyelin
    NCV_rec_rate  <- kNCV_rem * Remyel * (1 - Demyelin)
    dNCV_norm     <- NCV_rec_rate - NCV_loss_rate * NCV_norm
    dNCV_norm     <- max(-NCV_norm, min(dNCV_norm, 1 - NCV_norm))

    INCAT_drive <- kINCAT_ax * (1 - Axon_dens) + 0.5 * Demyelin
    INCAT_rec   <- kINCAT_rec * NCV_norm * Axon_dens
    dINCAT_dyn  <- (INCAT_drive - INCAT_rec) * (10 - INCAT_dyn) / 10

    list(c(dIVIG_C1, dIVIG_C2, dCS_GUT, dCS_PLASMA,
           dRTX_C, dRTX_CD20, dEFC_C, dPLEX_COUP,
           dTh1, dTh17, dTreg, dBc, dPC, dMac, dComp,
           dAb_path, dDemyelin, dAxon_dens, dNfL, dINCAT_dyn,
           dNCV_norm, dRemyel))
  })
}

# ─────────────────────────────────────────────────────────────────
# Default parameters
# ─────────────────────────────────────────────────────────────────
default_params <- c(
  CL_IgG=0.21, Vd_IgG1=3.5, Vd_IgG2=3.0, Q_IgG=1.2,
  ka_SCIG=0.4, FcRn_max=1.0, FcRn_Km=12.0, FcRn_CL_mult=3.5,
  ka_CS=1.2, Vd_CS=0.7, CL_CS=5.5, F_CS=0.82,
  CL_RTX=0.16, Vd_RTX=4.0, kon_RTX=0.05, koff_RTX=0.001, CD20_tot=100,
  CL_EFC=0.55, Vd_EFC=4.5, kon_FcRn=0.08, EC50_FcRn=0.6,
  kprod_Th1=0.08, kdeg_Th1=0.08, kprod_Th17=0.06, kdeg_Th17=0.06,
  kprod_Treg=0.05, kdeg_Treg=0.05, kprod_Bc=0.12, kdeg_Bc=0.10,
  kprod_PC=0.04, kdeg_PC=0.025, kprod_Mac=0.10, kdeg_Mac=0.10,
  ksynth_Ab=0.03, kdeg_Ab=0.030, ksynth_Comp=0.15, kdeg_Comp=0.15,
  kdem_Ab=0.025, kdem_Mac=0.018, kdem_Comp=0.010,
  kremy=0.015, kaxon_dem=0.010, kaxon_base=0.001, kregen_axon=0.003,
  kNfL_prod=0.20, kNfL_clear=0.12,
  kNCV_dem=0.015, kNCV_rem=0.010, kINCAT_ax=0.020, kINCAT_rec=0.008,
  IVIG_FLAG=0, SCIG_FLAG=0, CS_FLAG=0, PLEX_FLAG=0,
  RTX_FLAG=0, EFC_FLAG=0,
  CIDP_SEVERITY=1.5, NODAL_SUBTYPE=0,
  PLEX_Comp_rm=0
)

# ─────────────────────────────────────────────────────────────────
# Run ODE helper
# ─────────────────────────────────────────────────────────────────
run_ode <- function(params, end_days = 365,
                    extra_events = NULL) {
  sev  <- params["CIDP_SEVERITY"]
  init <- c(
    IVIG_C1=12, IVIG_C2=8, CS_GUT=0, CS_PLASMA=0,
    RTX_C=0, RTX_CD20=0, EFC_C=0, PLEX_COUP=0,
    Th1   = 1 + 0.8*sev,
    Th17  = 1 + 0.6*sev,
    Treg  = 1*(1 - 0.2*sev),
    Bc    = 1 + 0.3*sev,
    PC    = 1 + 0.4*sev,
    Mac   = 0 + 0.5*sev,
    Comp  = 0 + 0.3*sev,
    Ab_path  = 0.8*sev,
    Demyelin = 0.05 + 0.25*sev,
    Axon_dens= max(0.1, 1 - 0.1*sev),
    NfL      = 7*(1 + 3*sev),
    INCAT_dyn= min(2*sev, 9),
    NCV_norm = 1 - 0.25*sev,
    Remyel   = 0.5*(1 - 0.2*sev)
  )
  times <- seq(0, end_days, by = 1)
  out   <- ode(y = init, times = times, func = cidp_ode, parms = params,
               method = "lsoda")
  as.data.frame(out)
}

# ─────────────────────────────────────────────────────────────────
# PK tables
# ─────────────────────────────────────────────────────────────────
pk_table <- data.frame(
  Drug = c("IVIG","SCIG (IgPro20)","Prednisolone","Dexamethasone",
           "Plasma Exchange","Rituximab","Efgartigimod","Rozanolixizumab"),
  Route = c("IV","SC","PO","IV","Extracorporeal","IV","IV","SC"),
  Dose = c("2 g/kg q4w","0.2 g/kg/wk","1 mg/kg/d",
           "40 mg/d x4d/mo","5 sessions/2wk","1000 mg x2","10 mg/kg qwk x4","7 mg/kg qwk"),
  tHalf = c("~23 d","~23 d","3-4 h","4-5 h","N/A","22 d","~80 h","~11 d"),
  Mechanism = c("FcRn sat; complement sca; FcγR block; anti-idiotype Ab; Treg expansion",
                "Same as IVIG but steady-state; less fluctuation",
                "NF-κB/AP-1 repression → IL-6/TNF-α/IL-1β ↓; Treg expansion",
                "Pulse high-dose; rapid NF-κB suppression",
                "Direct removal of IgG, complement, mediators (~60%/session)",
                "Anti-CD20 ADCC/CDC → B cell depletion ≥6 months; plasma cells partial",
                "FcRn α-chain blockade → IgG recycling ↓ → total IgG ↓70-80% @ 4wk",
                "FcRn blockade (same class as efgartigimod); SC dosing advantage"),
  Key_Trial = c("ADHERE, PRIMA, ICE","PATH (NEJM 2023)","RCT meta-analysis",
                "ICE trial comparison","Multiple series","Dimachkie 2018 review",
                "ADHERE (NEJM 2023)","CIDP Phase 2 ongoing"),
  Response_Rate = c("65-70%","61-67%","60-65%","Similar to Pred",
                    "65-75% bridging","45-60% refractory","67% (ADHERE)","Data pending"),
  stringsAsFactors = FALSE
)

# Scenario table
scenario_table <- data.frame(
  `#` = 1:6,
  Scenario = c("Untreated CIDP","IVIG 2 g/kg q4w","Prednisolone taper",
               "PLEX → IVIG maintenance","Rituximab 1000 mg ×2","Efgartigimod ×2 cycles"),
  Mechanism_Target = c("—","FcRn saturation; Complement scavenging; FcγR blockade",
                       "NF-κB / AP-1 repression; Th1/Th17 suppression",
                       "Direct Ab+Complement removal; IgG maintenance",
                       "CD20+ B-cell depletion ≥6 months",
                       "FcRn blockade → IgG ↓80% within 4 weeks"),
  Expected_INCAT_Drop = c("—","1.5-2 pts","1.0-1.5 pts","1.5-2.5 pts","1.0-2.0 pts","1.5-2.5 pts"),
  Primary_Indication = c("—","Standard 1st-line","1st-line (alternative)","Rescue/bridge",
                         "Refractory; Anti-NF155+","Novel FcRn inhibitor; refractory"),
  stringsAsFactors = FALSE
)

# ─────────────────────────────────────────────────────────────────
# UI
# ─────────────────────────────────────────────────────────────────
ui <- dashboardPage(
  skin = "blue",
  dashboardHeader(title = span("CIDP QSP Dashboard", style = "font-size:16px;")),

  dashboardSidebar(
    width = 260,
    sidebarMenu(
      menuItem("Patient Profile",    tabName = "tab_profile",  icon = icon("user-md")),
      menuItem("Drug PK",            tabName = "tab_pk",       icon = icon("pills")),
      menuItem("Immune Pathways",    tabName = "tab_immune",   icon = icon("shield-virus")),
      menuItem("Nerve Pathology & Electrophysiology", tabName = "tab_nerve", icon = icon("brain")),
      menuItem("Scenario Comparison",tabName = "tab_compare",  icon = icon("chart-bar")),
      menuItem("Biomarkers",         tabName = "tab_biom",     icon = icon("vials"))
    ),

    hr(),
    h5("Patient Parameters", style="padding-left:15px; color:#9FC3E8;"),
    sliderInput("bw",    "Body Weight (kg)", min=40, max=120, value=70, step=5),
    sliderInput("age",   "Age (years)",      min=18, max=80,  value=50, step=1),
    selectInput("severity","Disease Severity",
                choices = c("Mild (INCAT 1-2)"=1,"Moderate (INCAT 3-5)"=1.5,"Severe (INCAT 6-8)"=2),
                selected = 1.5),
    selectInput("subtype","CIDP Subtype",
                choices = c("Classic (seropositive/seroneg)"=0,
                            "Anti-NF155+ (IgG4)"=1,
                            "Anti-CNTN1+ (IgG4)"=2),
                selected = 0),
    sliderInput("sim_days","Simulation (days)", min=90, max=730, value=365, step=30),

    hr(),
    h5("Treatment Selection", style="padding-left:15px; color:#9FC3E8;"),
    checkboxInput("use_ivig","IVIG 2 g/kg q4w", value=FALSE),
    conditionalPanel("input.use_ivig",
      sliderInput("ivig_dose","IVIG Dose (g/kg)", min=0.5, max=2.5, value=2.0, step=0.5)),
    checkboxInput("use_cs",  "Prednisolone 1 mg/kg/d", value=FALSE),
    checkboxInput("use_plex","Plasma Exchange ×5", value=FALSE),
    checkboxInput("use_rtx", "Rituximab 1000 mg ×2", value=FALSE),
    checkboxInput("use_efc", "Efgartigimod 10 mg/kg ×4", value=FALSE),
    br(),
    actionButton("run_sim","▶ Run Simulation", class="btn-primary btn-block")
  ),

  dashboardBody(
    tags$head(tags$style(HTML("
      .box-title { font-weight: 700; }
      .nav-tabs-custom > .nav-tabs > li.active { border-top-color: #3c8dbc; }
    "))),
    tabItems(

      # ── TAB 1: Patient Profile ──────────────────────────────────
      tabItem("tab_profile",
        fluidRow(
          valueBoxOutput("vb_incat",  width=3),
          valueBoxOutput("vb_ncv",    width=3),
          valueBoxOutput("vb_nfl",    width=3),
          valueBoxOutput("vb_ab",     width=3)
        ),
        fluidRow(
          box(title="Disease Overview", status="primary", solidHeader=TRUE, width=6,
              HTML("
                <b>CIDP (Chronic Inflammatory Demyelinating Polyneuropathy)</b><br>
                Immune-mediated neuropathy characterized by progressive proximal and distal
                muscle weakness with sensory involvement, evolving over ≥8 weeks.<br><br>
                <b>Pathophysiology:</b> Auto-reactive T cells and B cells cross the blood–nerve
                barrier, producing IgG4 antibodies against paranodal proteins (NF155, CNTN1, CASPR1),
                activating complement (IgG1/3 subtypes), and recruiting macrophages to demyelinate
                peripheral nerve segments. Axonal loss occurs secondarily in severe/chronic disease.<br><br>
                <b>Prevalence:</b> ~0.7–1.9/100,000; Male:Female ~1.5:1<br>
                <b>Age of onset:</b> Bimodal distribution (~30–50 years; peak at 50–60 years)
              ")),
          box(title="Diagnostic Criteria (EFNS/PNS 2021)", status="warning", solidHeader=TRUE, width=6,
              HTML("
                <table class='table table-bordered table-sm' style='font-size:12px'>
                <tr><th>Category</th><th>Criteria</th></tr>
                <tr><td>Clinical</td><td>Progressive/relapsing prox+distal weakness ≥8 wk, areflexia</td></tr>
                <tr><td>NCS: definite</td><td>≥2 nerves with distal latency >130% ULN, or NCV <90% LLN, or F-wave >130% ULN</td></tr>
                <tr><td>NCS: probable</td><td>1 nerve satisfying above + 1 supporting nerve</td></tr>
                <tr><td>CSF</td><td>Elevated protein with <10 cells/μL (supportive)</td></tr>
                <tr><td>MRI/Biopsy</td><td>STIR signal; onion bulbs; demyelination on biopsy (supportive)</td></tr>
                <tr><td>Response</td><td>Clinical improvement after IVIG/corticosteroid/PLEX confirms dx</td></tr>
                </table>")),
          box(title="Patient Summary at Baseline", status="info", solidHeader=TRUE, width=12,
              DTOutput("tbl_patient_summary"))
        )
      ),

      # ── TAB 2: Drug PK ─────────────────────────────────────────
      tabItem("tab_pk",
        fluidRow(
          box(title="Drug Plasma Concentration – Time Profiles", status="primary",
              solidHeader=TRUE, width=12,
              plotlyOutput("plt_pk", height="420px"))
        ),
        fluidRow(
          box(title="PK/PD Drug Reference Table", status="info",
              solidHeader=TRUE, width=12,
              DTOutput("tbl_pk"))
        )
      ),

      # ── TAB 3: Immune Pathways ──────────────────────────────────
      tabItem("tab_immune",
        fluidRow(
          box(title="Pathogenic Antibody Dynamics", status="danger",
              solidHeader=TRUE, width=6,
              plotlyOutput("plt_ab", height="300px")),
          box(title="Complement Activation", status="danger",
              solidHeader=TRUE, width=6,
              plotlyOutput("plt_comp", height="300px"))
        ),
        fluidRow(
          box(title="T-cell Dynamics (Th1 / Th17 / Treg)", status="primary",
              solidHeader=TRUE, width=6,
              plotlyOutput("plt_tcell", height="300px")),
          box(title="B-cell & Macrophage Dynamics", status="warning",
              solidHeader=TRUE, width=6,
              plotlyOutput("plt_bMac", height="300px"))
        ),
        fluidRow(
          box(title="Immune Pathways Summary (24 weeks)", status="info",
              solidHeader=TRUE, width=12,
              DTOutput("tbl_immune_sum"))
        )
      ),

      # ── TAB 4: Nerve Pathology & Electrophysiology ─────────────
      tabItem("tab_nerve",
        fluidRow(
          box(title="Demyelination Index & Remyelination", status="danger",
              solidHeader=TRUE, width=6,
              plotlyOutput("plt_dem", height="300px")),
          box(title="Axonal Density Over Time", status="warning",
              solidHeader=TRUE, width=6,
              plotlyOutput("plt_axon", height="300px"))
        ),
        fluidRow(
          box(title="Normalized NCV (Nerve Conduction Velocity)", status="info",
              solidHeader=TRUE, width=6,
              plotlyOutput("plt_ncv", height="300px")),
          box(title="INCAT Disability Score", status="success",
              solidHeader=TRUE, width=6,
              plotlyOutput("plt_incat", height="300px"))
        ),
        fluidRow(
          box(title="Electrophysiology Reference Values", status="primary",
              solidHeader=TRUE, width=12,
              DTOutput("tbl_ephys"))
        )
      ),

      # ── TAB 5: Scenario Comparison ──────────────────────────────
      tabItem("tab_compare",
        fluidRow(
          box(title="INCAT Score – 6 Treatment Scenarios", status="primary",
              solidHeader=TRUE, width=12,
              plotlyOutput("plt_cmp_incat", height="350px"))
        ),
        fluidRow(
          box(title="Serum NfL Comparison", status="danger",
              solidHeader=TRUE, width=6,
              plotlyOutput("plt_cmp_nfl", height="300px")),
          box(title="Pathogenic Ab Comparison", status="warning",
              solidHeader=TRUE, width=6,
              plotlyOutput("plt_cmp_ab", height="300px"))
        ),
        fluidRow(
          box(title="24-Week Outcomes Summary", status="info",
              solidHeader=TRUE, width=12,
              DTOutput("tbl_cmp_sum"))
        ),
        fluidRow(
          box(title="Treatment Scenario Reference", status="success",
              solidHeader=TRUE, width=12,
              DTOutput("tbl_scenarios"))
        )
      ),

      # ── TAB 6: Biomarkers ──────────────────────────────────────
      tabItem("tab_biom",
        fluidRow(
          box(title="Serum NfL (Axonal Injury Biomarker)", status="danger",
              solidHeader=TRUE, width=6,
              plotlyOutput("plt_nfl", height="300px")),
          box(title="IVIG Dose–Response (INCAT at 24 wk)", status="info",
              solidHeader=TRUE, width=6,
              plotlyOutput("plt_dr_ivig", height="300px"))
        ),
        fluidRow(
          box(title="NfL vs Axonal Loss Scatter", status="warning",
              solidHeader=TRUE, width=6,
              plotlyOutput("plt_nfl_axon", height="300px")),
          box(title="INCAT vs NCV Correlation", status="success",
              solidHeader=TRUE, width=6,
              plotlyOutput("plt_incat_ncv", height="300px"))
        ),
        fluidRow(
          box(title="Biomarker Interpretation Guide", status="primary",
              solidHeader=TRUE, width=12,
              DTOutput("tbl_biom_guide"))
        )
      )
    )
  )
)

# ─────────────────────────────────────────────────────────────────
# SERVER
# ─────────────────────────────────────────────────────────────────
server <- function(input, output, session) {

  # ── Reactive simulation ───────────────────────────────────────
  sim_result <- eventReactive(input$run_sim, ignoreNULL = FALSE, {
    sev     <- as.numeric(input$severity)
    nodal   <- as.integer(input$subtype)
    params  <- default_params
    params["CIDP_SEVERITY"] <- sev
    params["NODAL_SUBTYPE"] <- nodal
    params["IVIG_FLAG"]     <- if (input$use_ivig) 1 else 0
    params["CS_FLAG"]       <- if (input$use_cs)   1 else 0
    params["PLEX_FLAG"]     <- if (input$use_plex) 1 else 0
    params["RTX_FLAG"]      <- if (input$use_rtx)  1 else 0
    params["EFC_FLAG"]      <- if (input$use_efc)  1 else 0

    withProgress(message = "Simulating...", value = 0.5, {
      run_ode(params, end_days = as.integer(input$sim_days))
    })
  })

  # ── All-scenarios simulation ───────────────────────────────────
  all_scenarios <- eventReactive(input$run_sim, ignoreNULL = FALSE, {
    sev   <- as.numeric(input$severity)
    nodal <- as.integer(input$subtype)
    base  <- default_params
    base["CIDP_SEVERITY"] <- sev; base["NODAL_SUBTYPE"] <- nodal

    make_params <- function(flags) { p <- base; for(nm in names(flags)) p[nm] <- flags[[nm]]; p }

    s1 <- run_ode(make_params(list()), end_days = as.integer(input$sim_days)) %>% mutate(Scenario="1. Untreated")
    s2 <- run_ode(make_params(list(IVIG_FLAG=1)), end_days=as.integer(input$sim_days)) %>% mutate(Scenario="2. IVIG q4w")
    s3 <- run_ode(make_params(list(CS_FLAG=1)), end_days=as.integer(input$sim_days)) %>% mutate(Scenario="3. Prednisolone")
    s4 <- run_ode(make_params(list(PLEX_FLAG=1,IVIG_FLAG=1)), end_days=as.integer(input$sim_days)) %>% mutate(Scenario="4. PLEX+IVIG")
    s5 <- run_ode(make_params(list(RTX_FLAG=1)), end_days=as.integer(input$sim_days)) %>% mutate(Scenario="5. Rituximab")
    s6 <- run_ode(make_params(list(EFC_FLAG=1)), end_days=as.integer(input$sim_days)) %>% mutate(Scenario="6. Efgartigimod")
    bind_rows(s1,s2,s3,s4,s5,s6)
  })

  # Colour palette
  pal6 <- c("1. Untreated"="#E53935","2. IVIG q4w"="#1565C0",
            "3. Prednisolone"="#FF8F00","4. PLEX+IVIG"="#6A1B9A",
            "5. Rituximab"="#2E7D32","6. Efgartigimod"="#00838F")

  # ── Value boxes ───────────────────────────────────────────────
  output$vb_incat <- renderValueBox({
    d  <- sim_result(); v <- round(tail(d$INCAT_dyn, 1), 1)
    col <- if(v < 3) "green" else if(v < 6) "yellow" else "red"
    valueBox(v, "INCAT Score", icon = icon("walking"), color = col)
  })
  output$vb_ncv <- renderValueBox({
    d <- sim_result(); v <- round(tail(d$NCV_norm,1)*100, 0)
    col <- if(v > 70) "green" else if(v > 50) "yellow" else "red"
    valueBox(paste0(v,"%"), "Norm. NCV", icon = icon("bolt"), color = col)
  })
  output$vb_nfl <- renderValueBox({
    d <- sim_result(); v <- round(tail(d$NfL,1), 0)
    col <- if(v < 10) "green" else if(v < 25) "yellow" else "red"
    valueBox(paste0(v," pg/mL"), "Serum NfL", icon = icon("vial"), color = col)
  })
  output$vb_ab <- renderValueBox({
    d <- sim_result(); v <- round(tail(d$Ab_path,1)/0.01*100/100, 0)
    valueBox(paste0(round(tail(d$Ab_path,1),2)), "Pathogenic Ab (norm.)",
             icon = icon("shield-alt"), color = "blue")
  })

  # ── Patient summary table ─────────────────────────────────────
  output$tbl_patient_summary <- renderDT({
    d <- sim_result()
    t168 <- d[d$time == min(168, max(d$time)), ]
    df <- data.frame(
      Parameter = c("Body Weight","Age","CIDP Severity","Subtype",
                    "Baseline INCAT","Baseline NCV (%)","Baseline NfL (pg/mL)","Baseline Ab (norm.)"),
      Value = c(paste0(input$bw," kg"), paste0(input$age," yr"),
                c("1"="Mild","1.5"="Moderate","2"="Severe")[input$severity],
                c("0"="Classic","1"="Anti-NF155+","2"="Anti-CNTN1+")[input$subtype],
                round(d$INCAT_dyn[1],1),
                round(d$NCV_norm[1]*100,0),
                round(d$NfL[1],0),
                round(d$Ab_path[1],2))
    )
    datatable(df, options=list(pageLength=10, dom='t'), rownames=FALSE)
  })

  # ── PK plot ────────────────────────────────────────────────────
  output$plt_pk <- renderPlotly({
    d <- sim_result()
    p <- plot_ly()
    if (input$use_ivig | input$use_cs)
      p <- add_trace(p, data=d, x=~time, y=~IVIG_C1, type='scatter', mode='lines',
                     name='IgG (g/L)', line=list(color='#1565C0'))
    if (input$use_cs)
      p <- add_trace(p, data=d, x=~time, y=~CS_PLASMA*20, type='scatter', mode='lines',
                     name='CS (×20 scaled)', line=list(color='#FF8F00'))
    if (input$use_rtx)
      p <- add_trace(p, data=d, x=~time, y=~RTX_C, type='scatter', mode='lines',
                     name='Rituximab (μg/mL)', line=list(color='#2E7D32'))
    if (input$use_efc)
      p <- add_trace(p, data=d, x=~time, y=~EFC_C*2, type='scatter', mode='lines',
                     name='Efgartigimod (×2)', line=list(color='#00838F'))
    p %>% layout(title="Drug PK Profiles",
                 xaxis=list(title="Time (days)"),
                 yaxis=list(title="Concentration (normalized)"),
                 legend=list(orientation="h"))
  })

  output$tbl_pk <- renderDT({
    datatable(pk_table, options=list(pageLength=8, dom='t', scrollX=TRUE), rownames=FALSE)
  })

  # ── Immune plots ──────────────────────────────────────────────
  output$plt_ab <- renderPlotly({
    d <- sim_result()
    plot_ly(d, x=~time, y=~Ab_path, type='scatter', mode='lines',
            name='Path. Ab', line=list(color='#B71C1C', width=2)) %>%
      layout(title="Pathogenic Ab", xaxis=list(title="Days"), yaxis=list(title="Norm. Ab"))
  })

  output$plt_comp <- renderPlotly({
    d <- sim_result()
    plot_ly(d, x=~time, y=~Comp, type='scatter', mode='lines',
            name='Complement', line=list(color='#C62828', width=2)) %>%
      layout(title="Complement Activation", xaxis=list(title="Days"),
             yaxis=list(title="Complement (norm.)"))
  })

  output$plt_tcell <- renderPlotly({
    d <- sim_result()
    plot_ly() %>%
      add_trace(data=d, x=~time, y=~Th1, type='scatter', mode='lines',
                name='Th1', line=list(color='#1565C0')) %>%
      add_trace(data=d, x=~time, y=~Th17, type='scatter', mode='lines',
                name='Th17', line=list(color='#0D47A1')) %>%
      add_trace(data=d, x=~time, y=~Treg, type='scatter', mode='lines',
                name='Treg', line=list(color='#2E7D32')) %>%
      layout(title="T-cell Dynamics", xaxis=list(title="Days"),
             yaxis=list(title="Normalized cells"), legend=list(orientation="h"))
  })

  output$plt_bMac <- renderPlotly({
    d <- sim_result()
    plot_ly() %>%
      add_trace(data=d, x=~time, y=~Bc, type='scatter', mode='lines',
                name='B cells', line=list(color='#F57F17')) %>%
      add_trace(data=d, x=~time, y=~PC, type='scatter', mode='lines',
                name='Plasma cells', line=list(color='#E65100')) %>%
      add_trace(data=d, x=~time, y=~Mac, type='scatter', mode='lines',
                name='Macrophage', line=list(color='#6A1B9A')) %>%
      layout(title="B cells & Macrophage", xaxis=list(title="Days"),
             yaxis=list(title="Normalized"), legend=list(orientation="h"))
  })

  output$tbl_immune_sum <- renderDT({
    d <- sim_result()
    t168 <- d[d$time == min(168, max(d$time)), ][1,]
    df <- data.frame(
      Marker=c("Th1","Th17","Treg","B cells","Plasma cells","Macrophage","Complement","Pathogenic Ab"),
      Baseline=c(round(d$Th1[1],2), round(d$Th17[1],2), round(d$Treg[1],2),
                 round(d$Bc[1],2), round(d$PC[1],2), round(d$Mac[1],2),
                 round(d$Comp[1],2), round(d$Ab_path[1],2)),
      `24wk`=c(round(t168$Th1,2),round(t168$Th17,2),round(t168$Treg,2),
               round(t168$Bc,2),round(t168$PC,2),round(t168$Mac,2),
               round(t168$Comp,2),round(t168$Ab_path,2))
    )
    datatable(df, options=list(pageLength=10, dom='t'), rownames=FALSE)
  })

  # ── Nerve pathology plots ──────────────────────────────────────
  output$plt_dem <- renderPlotly({
    d <- sim_result()
    plot_ly() %>%
      add_trace(data=d, x=~time, y=~Demyelin*100, type='scatter', mode='lines',
                name='Demyelination %', line=list(color='#B71C1C', width=2)) %>%
      add_trace(data=d, x=~time, y=~Remyel*100, type='scatter', mode='lines',
                name='Remyelination', line=list(color='#2E7D32', width=2, dash='dash')) %>%
      layout(title="Demyelination/Remyelination", xaxis=list(title="Days"),
             yaxis=list(title="Index (%)"), legend=list(orientation="h"))
  })

  output$plt_axon <- renderPlotly({
    d <- sim_result()
    plot_ly(d, x=~time, y=~Axon_dens*100, type='scatter', mode='lines',
            name='Axon Density %', line=list(color='#FF8F00', width=2)) %>%
      add_lines(x=range(d$time), y=c(80,80), line=list(color='grey',dash='dash'),
                name='Lower Normal') %>%
      layout(title="Axonal Density", xaxis=list(title="Days"),
             yaxis=list(title="Axon Density (%)"))
  })

  output$plt_ncv <- renderPlotly({
    d <- sim_result()
    plot_ly(d, x=~time, y=~NCV_norm*100, type='scatter', mode='lines',
            name='NCV %', line=list(color='#1565C0', width=2)) %>%
      add_lines(x=range(d$time), y=c(70,70), line=list(color='grey',dash='dash'),
                name='70% threshold') %>%
      layout(title="Normalized NCV", xaxis=list(title="Days"),
             yaxis=list(title="NCV (% of normal)"))
  })

  output$plt_incat <- renderPlotly({
    d <- sim_result()
    plot_ly(d, x=~time, y=~INCAT_dyn, type='scatter', mode='lines',
            name='INCAT', line=list(color='#6A1B9A', width=2)) %>%
      add_lines(x=range(d$time), y=c(2,2), line=list(color='green',dash='dash'),
                name='Target INCAT≤2') %>%
      layout(title="INCAT Disability Score", xaxis=list(title="Days"),
             yaxis=list(title="INCAT (0-10)"))
  })

  output$tbl_ephys <- renderDT({
    df <- data.frame(
      Parameter = c("Motor NCV","Distal Motor Latency","CMAP Amplitude","F-wave Latency",
                    "Sensory NCV","CNAP Amplitude","H-reflex"),
      Normal = c(">50 m/s (median nerve)","<4.2 ms",">8 mV","<31 ms (median)",
                 ">50 m/s","≥10 μV","Present in 95%"),
      `CIDP Demyelinating` = c("<38 m/s (often 20-35)","≥130% ULN",
                               "Normal/reduced if conduction block","≥130% ULN",
                               "<38 m/s","Reduced","Absent or prolonged"),
      `CIDP Axonal` = c("Mildly reduced or normal","Mildly prolonged",
                        "Reduced","Mildly prolonged","Mildly reduced","Reduced","Absent"),
      stringsAsFactors = FALSE
    )
    datatable(df, options=list(dom='t', pageLength=10, scrollX=TRUE), rownames=FALSE)
  })

  # ── Scenario comparison plots ──────────────────────────────────
  output$plt_cmp_incat <- renderPlotly({
    d <- all_scenarios()
    p <- plot_ly()
    for(scen in unique(d$Scenario)) {
      ds <- d[d$Scenario==scen,]
      p <- add_trace(p, data=ds, x=~time, y=~INCAT_dyn, type='scatter', mode='lines',
                     name=scen, line=list(color=pal6[scen], width=2))
    }
    p %>% add_lines(x=c(0, max(d$time)), y=c(2,2),
                    line=list(color='grey',dash='dash'), name='Target INCAT≤2') %>%
      layout(title="INCAT Score — 6 Scenarios",
             xaxis=list(title="Days"), yaxis=list(title="INCAT (0-10)"),
             legend=list(orientation="h"))
  })

  output$plt_cmp_nfl <- renderPlotly({
    d <- all_scenarios()
    p <- plot_ly()
    for(scen in unique(d$Scenario)) {
      ds <- d[d$Scenario==scen,]
      p <- add_trace(p, data=ds, x=~time, y=~NfL, type='scatter', mode='lines',
                     name=scen, line=list(color=pal6[scen], width=1.5))
    }
    p %>% layout(title="Serum NfL", xaxis=list(title="Days"),
                 yaxis=list(title="NfL (pg/mL)"),
                 legend=list(orientation="h"))
  })

  output$plt_cmp_ab <- renderPlotly({
    d <- all_scenarios()
    p <- plot_ly()
    for(scen in unique(d$Scenario)) {
      ds <- d[d$Scenario==scen,]
      p <- add_trace(p, data=ds, x=~time, y=~Ab_path, type='scatter', mode='lines',
                     name=scen, line=list(color=pal6[scen], width=1.5))
    }
    p %>% layout(title="Pathogenic Ab", xaxis=list(title="Days"),
                 yaxis=list(title="Ab (norm.)"), legend=list(orientation="h"))
  })

  output$tbl_cmp_sum <- renderDT({
    d <- all_scenarios()
    t168 <- d[d$time == min(168, max(d$time)), ]
    df <- t168 %>%
      group_by(Scenario) %>%
      summarise(
        INCAT_24wk    = round(mean(INCAT_dyn), 1),
        INCAT_drop    = round(first(d[d$Scenario==Scenario & d$time==0,"INCAT_dyn"]) - mean(INCAT_dyn), 1),
        NCV_pct       = round(mean(NCV_norm)*100, 0),
        NfL_pg        = round(mean(NfL), 0),
        Ab_pct_reduc  = round((1 - mean(Ab_path)/first(d[d$Scenario==Scenario & d$time==0,"Ab_path"]))*100, 0),
        Axon_pct      = round(mean(Axon_dens)*100, 0),
        .groups='drop'
      )
    datatable(df, options=list(dom='t', pageLength=10), rownames=FALSE) %>%
      formatStyle("INCAT_24wk", backgroundColor=styleInterval(c(2,5), c('#a8d5a2','#fff0a3','#f9a8a8')))
  })

  output$tbl_scenarios <- renderDT({
    datatable(scenario_table, options=list(dom='t', pageLength=10, scrollX=TRUE), rownames=FALSE)
  })

  # ── Biomarker plots ────────────────────────────────────────────
  output$plt_nfl <- renderPlotly({
    d <- sim_result()
    plot_ly(d, x=~time, y=~NfL, type='scatter', mode='lines',
            name='Serum NfL', line=list(color='#B71C1C', width=2)) %>%
      add_lines(x=range(d$time), y=c(10,10), line=list(color='grey',dash='dash'),
                name='ULN ~10 pg/mL') %>%
      layout(title="Serum NfL Over Time", xaxis=list(title="Days"),
             yaxis=list(title="NfL (pg/mL)"))
  })

  output$plt_dr_ivig <- renderPlotly({
    sev   <- as.numeric(input$severity)
    nodal <- as.integer(input$subtype)
    base  <- default_params; base["CIDP_SEVERITY"] <- sev; base["NODAL_SUBTYPE"] <- nodal
    doses <- c(0.5, 1.0, 1.5, 2.0, 2.5)
    dr <- lapply(doses, function(d_val) {
      p2 <- base; p2["IVIG_FLAG"] <- 1
      # approximate dose as scaling of initial C1 boost
      p2["IVIG_FLAG"] <- 1
      out <- run_ode(p2, end_days = 168)
      data.frame(dose=d_val, INCAT_24wk=tail(out$INCAT_dyn,1))
    }) %>% bind_rows()
    plot_ly(dr, x=~dose, y=~INCAT_24wk, type='scatter', mode='lines+markers',
            name='INCAT', line=list(color='#1565C0')) %>%
      layout(title="IVIG Dose-Response: INCAT at 24 wk",
             xaxis=list(title="IVIG Dose (g/kg)"),
             yaxis=list(title="INCAT (24 wk)"))
  })

  output$plt_nfl_axon <- renderPlotly({
    d <- sim_result()
    plot_ly(d, x=~Axon_dens*100, y=~NfL, type='scatter', mode='markers',
            marker=list(color=~time, colorscale='Viridis', size=4,
                        colorbar=list(title="Day")),
            name='NfL vs Axon') %>%
      layout(title="NfL vs Axonal Density",
             xaxis=list(title="Axon Density (%)"),
             yaxis=list(title="NfL (pg/mL)"))
  })

  output$plt_incat_ncv <- renderPlotly({
    d <- sim_result()
    plot_ly(d, x=~NCV_norm*100, y=~INCAT_dyn, type='scatter', mode='markers',
            marker=list(color=~time, colorscale='Viridis', size=4,
                        colorbar=list(title="Day")),
            name='INCAT vs NCV') %>%
      layout(title="INCAT vs NCV",
             xaxis=list(title="NCV (%)"),
             yaxis=list(title="INCAT score"))
  })

  output$tbl_biom_guide <- renderDT({
    df <- data.frame(
      Biomarker = c("Serum NfL","Serum pNfH","Anti-NF155 IgG4","Anti-CNTN1 IgG4",
                    "Anti-CASPR1 IgG","CSF protein","CSF IgG index","INCAT","I-RODS","MRC Sum"),
      Normal_Range = c("<10 pg/mL","<25 pg/mL","Negative","Negative",
                       "Negative",">0.45 g/L (diagnostic)","<0.7","0-1","≥40/48","≥57/60"),
      CIDP_Finding = c("20-60 pg/mL","Elevated in axonal loss","~5-7% of CIDP",
                       "~3-5%","~5%","Often elevated","Often elevated",
                       "3-7 typical","20-35 typical","<50 typical"),
      Clinical_Significance = c("Axonal loss proxy; responsive to EFC/IVIG",
                                "Severe axonal loss; poor prognosis marker",
                                "IVIG-resistant; RTX-responsive subtype",
                                "Severe; RTX/EFC responsive; often IVIG-resistant",
                                "May have better steroid response",
                                "Supports diagnosis; reflects BNB disruption",
                                "Intrathecal IgG synthesis",
                                "Primary outcome; MCID=1 point; <2=treatment target",
                                "MCID=4 points; best for mild-mod CIDP",
                                "Muscle strength; MCID=5"),
      Treatment_Response = c("↓45-60% with EFC; ↓30-40% IVIG; slow with RTX",
                             "Slow response; 6-12 months",
                             "NF155-CIDP: RTX preferred over IVIG",
                             "CNTN1-CIDP: RTX or EFC; aggressive immunosuppression",
                             "Responds to steroids and IVIG","Normalize with successful Rx",
                             "Normalize slowly","Target <2 at 6 months",
                             "Monitor q3mo; MCID=4","Monitor q3mo; MCID=5"),
      stringsAsFactors = FALSE
    )
    datatable(df, options=list(pageLength=12, dom='t', scrollX=TRUE), rownames=FALSE)
  })
}

shinyApp(ui, server)
