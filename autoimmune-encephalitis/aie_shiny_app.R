## =============================================================================
## Anti-NMDA Receptor Encephalitis (AIE) — Interactive Shiny Dashboard
## =============================================================================
## Tabs:
##   1. Patient Profile & Disease Overview
##   2. Immunology & Antibody PK
##   3. CNS Pathophysiology (NMDAR, BBB, Neuroinflammation)
##   4. Clinical Endpoints (mRS, Seizures, Cognition, Psychiatry)
##   5. Scenario Comparison (6 treatment strategies)
##   6. Biomarkers & Diagnostics
## =============================================================================

library(shiny)
library(shinydashboard)
library(dplyr)
library(tidyr)
library(ggplot2)
library(plotly)
library(DT)
library(deSolve)

# ─────────────────────────────────────────────────────────────────────────────
# ODE SYSTEM (deSolve implementation for Shiny)
# ─────────────────────────────────────────────────────────────────────────────
aie_ode <- function(t, state, parms) {
  with(as.list(c(state, parms)), {

    # ── Drug concentrations ────────────────────────────────────────────────
    Cp_IVIG_mgl  <- IVIG1 / Vc_IVIG
    Cp_MP_mcg    <- (MP1  / Vc_MP)  * 1000
    Cp_RTX_mcg   <- (RTX1 / Vc_RTX) * 1000
    Cp_TCZ_mcg   <- (TCZ1 / Vc_TCZ) * 1000
    Cp_CPX_ng    <- (CPX_ACT / Vc_CPX) * 1e6

    # ── Drug effects ───────────────────────────────────────────────────────
    eff_RTX  <- Emax_RTX  * Cp_RTX_mcg^gamma_RTX /
                (EC50_RTX^gamma_RTX + Cp_RTX_mcg^gamma_RTX)
    eff_IVIG <- Emax_IVIG * Cp_IVIG_mgl / (EC50_IVIG + Cp_IVIG_mgl)
    eff_MP_a <- Emax_MP_a * Cp_MP_mcg   / (EC50_MP   + Cp_MP_mcg)
    eff_MP_b <- Emax_MP_b * Cp_MP_mcg   / (EC50_MP   + Cp_MP_mcg)
    eff_TCZ  <- Emax_TCZ  * Cp_TCZ_mcg  / (EC50_TCZ  + Cp_TCZ_mcg)
    eff_CPX  <- Emax_CPX  * Cp_CPX_ng   / (EC50_CPX  + Cp_CPX_ng)

    # ── Derived quantities ────────────────────────────────────────────────
    BBB_c       <- max(BBB_min, BBB)
    NMDAR_c     <- max(NMDAR_min, NMDAR)
    BBB_open    <- max(0, 1 - BBB_c)
    MG_excess   <- max(0, MG - 1)
    IL6_active  <- IL6_CNS * (1 - eff_TCZ)
    NMDAR_loss  <- max(0, 1 - NMDAR_c)
    DA_dis      <- max(0, NMDAR_loss * 0.8 + max(0, GLU - 1) * 0.2)
    Ab_prod     <- k_Ab_prod * (LLPC + 2 * PB)
    Ab_clear    <- k_Ab_ser * (1 + eff_IVIG)
    GLU_c       <- max(1, min(GLU_max, GLU))

    # ── ODEs ──────────────────────────────────────────────────────────────
    dGCB      <- k_GCB_stim * GCB * (1 - GCB/500) + k_GCB_MB * MB -
                 k_GCB_death * GCB - eff_RTX * k_GCB_death * GCB * 5 -
                 eff_CPX * GCB * 0.3
    dPB       <- k_PB_GCB * GCB - k_PB_death * PB - eff_CPX * PB * 0.5
    dLLPC     <- k_LLPC_PB * PB - k_LLPC_death * LLPC - eff_CPX * LLPC * 0.15
    dMB       <- k_MB_GCB * GCB - k_MB_death * MB - eff_RTX * MB * 0.9

    dAB_SERUM <- Ab_prod - Ab_clear * AB_SERUM
    dAB_CSF   <- k_Ab_trans * AB_SERUM * (BBB_open + 0.05) -
                 k_Ab_CSF * AB_CSF

    dBBB      <- k_BBB_rep * (1 - BBB_c) * (1 + eff_MP_b) -
                 k_BBB_MG * MG_excess * BBB_c -
                 k_BBB_IL6 * max(0, IL6_active - 1) * BBB_c

    dMG       <- k_MG_act * AB_CSF - k_MG_res * MG_excess * (1 + eff_MP_a * 2)
    dIL6_CNS  <- k_IL6_MG * MG_excess -
                 k_IL6_cl * max(0, IL6_CNS - 1) -
                 eff_TCZ * k_IL6_MG * MG_excess
    dGFAP     <- k_GFAP_MG * MG_excess - k_GFAP_res * max(0, GFAP - 1)

    dNMDAR    <- k_NMDAR_base * (1 - NMDAR_c) -
                 k_NMDAR_int * AB_CSF * NMDAR_c +
                 k_NMDAR_rec * (1 - NMDAR_c) * ifelse(AB_CSF < 0.05, 1, 0.2)

    dGLU      <- k_GLU_exc * NMDAR_loss - k_GLU_cl * max(0, GLU - 1)

    dis_load <- k_CRS_N * (1 - NMDAR_c) + k_CRS_G * max(0, GLU - 1) +
                0.5 * MG_excess + 0.3 * (AB_CSF / 0.1)
    dCRS     <- dis_load - k_CRS_rec * CRS
    GLU_over <- max(0, GLU - SZ_thresh)
    dSZ      <- k_SZ_GLU * GLU_over - k_SZ_res * SZ
    dCOG     <- -k_COG_loss * NMDAR_loss * max(0, COG) +
                 k_COG_rec  * max(0, 1 - COG) * ifelse(NMDAR > 0.70, 1, 0)
    dPSY     <- k_PSY_DA * DA_dis - k_PSY_res * PSY

    # Drug PK
    dIVIG1   <- -(CL_IVIG + Q_IVIG)/Vc_IVIG * IVIG1 + Q_IVIG/Vp_IVIG * IVIG2
    dIVIG2   <-  Q_IVIG/Vc_IVIG * IVIG1 - Q_IVIG/Vp_IVIG * IVIG2
    dMP1     <- -(CL_MP + Q_MP)/Vc_MP * MP1   + Q_MP/Vp_MP   * MP2
    dMP2     <-  Q_MP/Vc_MP * MP1     - Q_MP/Vp_MP * MP2
    dRTX1    <- -(CL_RTX + Q_RTX)/Vc_RTX * RTX1 + Q_RTX/Vp_RTX * RTX2
    dRTX2    <-  Q_RTX/Vc_RTX * RTX1   - Q_RTX/Vp_RTX * RTX2
    dTCZ1    <- -(CL_TCZ + Q_TCZ)/Vc_TCZ * TCZ1 + Q_TCZ/Vp_TCZ * TCZ2
    dTCZ2    <-  Q_TCZ/Vc_TCZ * TCZ1   - Q_TCZ/Vp_TCZ * TCZ2
    dCPX_ACT <- -CL_CPX / Vc_CPX * CPX_ACT

    list(c(dGCB, dPB, dLLPC, dMB, dAB_SERUM, dAB_CSF,
           dBBB, dMG, dIL6_CNS, dGFAP, dNMDAR, dGLU,
           dCRS, dSZ, dCOG, dPSY,
           dIVIG1, dIVIG2, dMP1, dMP2,
           dRTX1, dRTX2, dTCZ1, dTCZ2, dCPX_ACT))
  })
}

# ─────────────────────────────────────────────────────────────────────────────
# DEFAULT PARAMETERS
# ─────────────────────────────────────────────────────────────────────────────
default_parms <- list(
  k_GCB_stim=0.08, k_GCB_death=0.015, k_GCB_MB=0.003,
  k_PB_GCB=0.06, k_PB_death=0.12,
  k_LLPC_PB=0.01, k_LLPC_death=0.0006,
  k_MB_GCB=0.03, k_MB_death=0.003,
  k_Ab_prod=0.0003, k_Ab_ser=0.030,
  k_Ab_trans=0.08, k_Ab_CSF=0.12,
  k_BBB_rep=0.04, k_BBB_MG=0.10, k_BBB_IL6=0.04, BBB_min=0.15,
  k_MG_act=0.5, k_MG_res=0.08,
  k_IL6_MG=0.25, k_IL6_cl=0.18,
  k_GFAP_MG=0.12, k_GFAP_res=0.09,
  k_NMDAR_base=0.020, k_NMDAR_int=0.30, k_NMDAR_rec=0.012, NMDAR_min=0.08,
  k_GLU_exc=0.25, k_GLU_cl=0.60, GLU_max=4.0,
  k_CRS_N=3.0, k_CRS_G=1.2, k_CRS_rec=0.015,
  k_SZ_GLU=0.8, SZ_thresh=1.6, k_SZ_res=0.10,
  k_COG_loss=0.25, k_COG_rec=0.008,
  k_PSY_DA=0.90, k_PSY_res=0.06,
  CL_IVIG=0.210, Vc_IVIG=3.7, Vp_IVIG=25.0, Q_IVIG=1.2,
  CL_MP=576, Vc_MP=28.0, Vp_MP=56.0, Q_MP=1728,
  CL_RTX=0.336, Vc_RTX=3.4, Vp_RTX=4.4, Q_RTX=0.43,
  CL_TCZ=0.55, Vc_TCZ=3.5, Vp_TCZ=2.9, Q_TCZ=0.48,
  CL_CPX=96.0, Vc_CPX=30.0,
  EC50_RTX=8.0, Emax_RTX=0.98, gamma_RTX=2.0,
  EC50_IVIG=12.0, Emax_IVIG=4.5,
  EC50_MP=0.25, Emax_MP_a=0.80, Emax_MP_b=0.55,
  EC50_TCZ=2.5, Emax_TCZ=0.95,
  EC50_CPX=500, Emax_CPX=0.92
)

y0 <- c(GCB=100, PB=50, LLPC=200, MB=100, AB_SERUM=0.01, AB_CSF=0.001,
        BBB=1.0, MG=1.0, IL6_CNS=1.0, GFAP=1.0, NMDAR=1.0, GLU=1.0,
        CRS=0, SZ=0, COG=1.0, PSY=0,
        IVIG1=0, IVIG2=0, MP1=0, MP2=0, RTX1=0, RTX2=0,
        TCZ1=0, TCZ2=0, CPX_ACT=0)

# ─────────────────────────────────────────────────────────────────────────────
# SIMULATION FUNCTION
# ─────────────────────────────────────────────────────────────────────────────
run_sim <- function(parms, weight_kg=70, tx_start=14,
                    use_IVIG=FALSE, use_MP=FALSE, use_PE=FALSE,
                    use_RTX=FALSE, use_CPX=FALSE, use_TCZ=FALSE,
                    RTX_delay=30, TCZ_delay=60) {
  times <- seq(0, 365, by=1)
  p <- modifyList(default_parms, parms)
  p[["Vc_MP"]] <- 0.4 * weight_kg

  # Build event dataframe
  events_list <- data.frame(var=character(), time=numeric(), value=numeric(),
                             method=character(), stringsAsFactors=FALSE)

  add_bolus <- function(cmt, times_v, amt) {
    data.frame(var=cmt, time=times_v, value=amt, method="add", stringsAsFactors=FALSE)
  }

  if(use_IVIG) {
    ivig_amt <- 2000 * weight_kg / 5  # 2g/kg over 5 days, split
    events_list <- rbind(events_list,
      add_bolus("IVIG1", tx_start:(tx_start+4), ivig_amt))
  }
  if(use_MP) {
    events_list <- rbind(events_list,
      add_bolus("MP1", tx_start:(tx_start+4), 1000))
  }
  if(use_PE) {
    pe_times <- tx_start + c(0,2,4,6,8)
    events_list <- rbind(events_list,
      data.frame(var="AB_SERUM", time=pe_times, value=0.25, method="mult",
                 stringsAsFactors=FALSE))
  }
  if(use_RTX) {
    rtx_times <- RTX_delay + c(0,7,14,21)
    rtx_dose  <- 375 * (weight_kg^0.425 * 160^0.725 * 0.007184) # BSA approx
    rtx_dose  <- min(max(rtx_dose, 500), 800)
    events_list <- rbind(events_list,
      add_bolus("RTX1", rtx_times, rtx_dose))
  }
  if(use_CPX) {
    cpx_times <- RTX_delay + c(0,28,56,84,112,140)
    cpx_dose  <- 0.25 * 750 * (weight_kg^0.425 * 160^0.725 * 0.007184)
    events_list <- rbind(events_list,
      add_bolus("CPX_ACT", cpx_times, cpx_dose))
  }
  if(use_TCZ) {
    tcz_times <- TCZ_delay + c(0,28,56,84,112,140)
    events_list <- rbind(events_list,
      add_bolus("TCZ1", tcz_times, 8 * weight_kg))
  }

  if(nrow(events_list) > 0) {
    events_df <- events_list
    out <- ode(y=y0, times=times, func=aie_ode, parms=p,
               method="lsoda", events=list(data=events_df))
  } else {
    out <- ode(y=y0, times=times, func=aie_ode, parms=p, method="lsoda")
  }

  df <- as.data.frame(out)
  df$NMDAR_pct <- pmax(8, df$NMDAR * 100)
  df$mRS_est   <- pmin(6, pmax(0, df$CRS))
  df$COG_pct   <- pmax(0, pmin(100, df$COG * 100))
  df$BBB_pct   <- pmax(15, pmin(100, df$BBB * 100))
  df$Cp_IVIG   <- df$IVIG1 / 3.7
  df$Cp_RTX    <- df$RTX1  / 3.4 * 1000
  df$Cp_MP     <- df$MP1   / 28  * 1000
  df$Cp_TCZ    <- df$TCZ1  / 3.5 * 1000
  df$RO_RTX    <- 0.98 * df$Cp_RTX^2 / (8^2 + df$Cp_RTX^2)
  return(df)
}

# ─────────────────────────────────────────────────────────────────────────────
# SHINY UI
# ─────────────────────────────────────────────────────────────────────────────
ui <- dashboardPage(
  skin = "blue",
  dashboardHeader(title = "AIE QSP Dashboard", titleWidth = 280),

  dashboardSidebar(
    width = 280,
    sidebarMenu(
      menuItem("① Patient Profile",    tabName="tab_profile",  icon=icon("user")),
      menuItem("② Antibody PK",        tabName="tab_pk",       icon=icon("flask")),
      menuItem("③ CNS Pathophysiology",tabName="tab_cns",      icon=icon("brain")),
      menuItem("④ Clinical Endpoints", tabName="tab_clin",     icon=icon("stethoscope")),
      menuItem("⑤ Scenario Comparison",tabName="tab_comp",     icon=icon("chart-bar")),
      menuItem("⑥ Biomarkers",         tabName="tab_bm",       icon=icon("vial"))
    ),
    hr(),
    h5("  Patient Parameters", style="color:#ECF0F1; margin-left:10px"),
    sliderInput("weight",  "Weight (kg)",   40, 120, 70, step=5),
    sliderInput("tx_delay","Diagnosis Delay (d)", 0, 60, 14, step=1),
    hr(),
    h5("  Treatment Selection", style="color:#ECF0F1; margin-left:10px"),
    checkboxInput("use_IVIG","IVIG (2g/kg ×5d)", FALSE),
    checkboxInput("use_MP",  "Methylprednisolone (1g/d ×5d)", FALSE),
    checkboxInput("use_PE",  "Plasmapheresis (×5)", FALSE),
    checkboxInput("use_RTX", "Rituximab (375mg/m² ×4)", FALSE),
    checkboxInput("use_CPX", "Cyclophosphamide (750mg/m² ×6)", FALSE),
    checkboxInput("use_TCZ", "Tocilizumab (8mg/kg ×6)", FALSE),
    sliderInput("rtx_delay","2L Start (day)", 21, 90, 30, step=1),
    sliderInput("tcz_delay","3L Start (day)", 45, 120, 60, step=1),
    hr(),
    actionButton("run_btn","▶  Run Simulation", class="btn-success", width="90%")
  ),

  dashboardBody(
    tags$head(tags$style(HTML("
      .content-wrapper { background-color: #f4f6f9; }
      .box { border-radius: 6px; }
      .shiny-plot-output { width:100% !important; }
    "))),
    tabItems(

      # ── TAB 1: PATIENT PROFILE ─────────────────────────────────────────────
      tabItem("tab_profile",
        fluidRow(
          box(width=12, title="Anti-NMDA Receptor Encephalitis: Disease Overview",
              status="primary", solidHeader=TRUE,
              p(strong("Anti-NMDAR encephalitis"), " is the most common autoimmune encephalitis, caused by",
                " IgG antibodies against the GluN1 (NR1) subunit of NMDA receptors.",
                " First described by Dalmau et al. in 2007, it affects predominantly young women (80%)",
                " and is associated with ovarian teratoma in ~40% of women aged 18-45."),
              p("Clinical presentation follows a stereotyped sequence:",
                " (1) prodromal phase (fever, headache); (2) psychiatric symptoms (psychosis, mood disorder);",
                " (3) seizures; (4) movement disorders and orofacial dyskinesia;",
                " (5) decreased consciousness and autonomic instability."),
              p(strong("Treatment:"), " First-line = IVIG + steroids ± plasmapheresis;",
                " Second-line = Rituximab or cyclophosphamide for refractory cases;",
                " Tumor resection if ovarian teratoma present (most important intervention)."))
        ),
        fluidRow(
          box(width=6, title="Diagnostic Criteria (Graus et al. 2016)", status="info",
              DTOutput("tbl_diagnosis")),
          box(width=6, title="Clinical Stage Progression", status="warning",
              plotlyOutput("plt_stages", height="250px"))
        ),
        fluidRow(
          box(width=6, title="Treatment Algorithm", status="success",
              DTOutput("tbl_treatment")),
          box(width=6, title="Epidemiology & Prognosis", status="info",
              DTOutput("tbl_epidemio"))
        )
      ),

      # ── TAB 2: ANTIBODY PK ────────────────────────────────────────────────
      tabItem("tab_pk",
        fluidRow(
          box(width=6, title="Serum Anti-NMDAR IgG (Relative Titer)",
              status="primary", solidHeader=TRUE, plotlyOutput("plt_Ab_serum")),
          box(width=6, title="CSF Anti-NMDAR IgG",
              status="primary", solidHeader=TRUE, plotlyOutput("plt_Ab_CSF"))
        ),
        fluidRow(
          box(width=6, title="Drug Concentration Profiles (IVIG, MP, RTX, TCZ)",
              status="warning", solidHeader=TRUE, plotlyOutput("plt_PK_drugs")),
          box(width=6, title="Rituximab Receptor Occupancy (CD20)",
              status="info", solidHeader=TRUE, plotlyOutput("plt_RO_RTX"))
        ),
        fluidRow(
          box(width=12, title="PK Parameter Summary", status="success",
              DTOutput("tbl_pk_params"))
        )
      ),

      # ── TAB 3: CNS PATHOPHYSIOLOGY ────────────────────────────────────────
      tabItem("tab_cns",
        fluidRow(
          box(width=6, title="Surface NMDA-R Density (%)",
              status="danger", solidHeader=TRUE, plotlyOutput("plt_NMDAR")),
          box(width=6, title="BBB Integrity (%)",
              status="warning", solidHeader=TRUE, plotlyOutput("plt_BBB"))
        ),
        fluidRow(
          box(width=6, title="Microglia Activation Index",
              status="primary", solidHeader=TRUE, plotlyOutput("plt_MG")),
          box(width=6, title="CNS IL-6 & Astrocyte (GFAP)",
              status="info", solidHeader=TRUE, plotlyOutput("plt_IL6_GFAP"))
        ),
        fluidRow(
          box(width=6, title="Synaptic Glutamate (E/I imbalance)",
              status="danger", solidHeader=TRUE, plotlyOutput("plt_GLU")),
          box(width=6, title="CNS Pathophysiology Cascade",
              status="success", DTOutput("tbl_cns_cascade"))
        )
      ),

      # ── TAB 4: CLINICAL ENDPOINTS ─────────────────────────────────────────
      tabItem("tab_clin",
        fluidRow(
          box(width=6, title="Clinical Severity (mRS estimate)",
              status="danger", solidHeader=TRUE, plotlyOutput("plt_mRS")),
          box(width=6, title="Seizure Frequency (events/week)",
              status="warning", solidHeader=TRUE, plotlyOutput("plt_SZ"))
        ),
        fluidRow(
          box(width=6, title="Cognitive Index (% of Normal)",
              status="primary", solidHeader=TRUE, plotlyOutput("plt_COG")),
          box(width=6, title="Psychiatric Symptom Score",
              status="info", solidHeader=TRUE, plotlyOutput("plt_PSY"))
        ),
        fluidRow(
          box(width=12, title="Clinical Summary Table (Weeks 2/4/8/12/24/52)",
              status="success", DTOutput("tbl_clinical_summary"))
        )
      ),

      # ── TAB 5: SCENARIO COMPARISON ────────────────────────────────────────
      tabItem("tab_comp",
        fluidRow(
          box(width=12, title="6 Treatment Scenarios — Week 12 Outcomes",
              status="primary", solidHeader=TRUE,
              p("Compare: (1) No Tx  (2) IVIG+MP  (3) IVIG+MP+PE  (4) +Rituximab  (5) +Cyclophosphamide  (6) +Tocilizumab"),
              plotlyOutput("plt_scenario_bar", height="350px"))
        ),
        fluidRow(
          box(width=6, title="NMDAR Recovery Over Time — All Scenarios",
              status="success", plotlyOutput("plt_scen_NMDAR")),
          box(width=6, title="Cognitive Index — All Scenarios",
              status="info", plotlyOutput("plt_scen_COG"))
        ),
        fluidRow(
          box(width=12, title="Scenario Comparison Table (Day 90 & Day 180)",
              status="warning", DTOutput("tbl_scenario_compare"))
        )
      ),

      # ── TAB 6: BIOMARKERS ─────────────────────────────────────────────────
      tabItem("tab_bm",
        fluidRow(
          box(width=6, title="CSF IL-6 Trajectory (Disease Activity Marker)",
              status="info", solidHeader=TRUE, plotlyOutput("plt_bm_IL6")),
          box(width=6, title="Astrocyte Activation (CSF GFAP Proxy)",
              status="warning", solidHeader=TRUE, plotlyOutput("plt_bm_GFAP"))
        ),
        fluidRow(
          box(width=6, title="Immune Cell Dynamics (GCB, PB, LLPC, MB)",
              status="primary", solidHeader=TRUE, plotlyOutput("plt_bm_bcells")),
          box(width=6, title="Dose-Response: RTX vs NMDAR Recovery (Day 90)",
              status="danger", solidHeader=TRUE, plotlyOutput("plt_DR_RTX"))
        ),
        fluidRow(
          box(width=12, title="Key Biomarkers Reference Values",
              status="success", DTOutput("tbl_biomarkers"))
        )
      )
    )
  )
)

# ─────────────────────────────────────────────────────────────────────────────
# SHINY SERVER
# ─────────────────────────────────────────────────────────────────────────────
server <- function(input, output, session) {

  # ── Reactive simulation ──────────────────────────────────────────────────
  sim_data <- eventReactive(input$run_btn, {
    withProgress(message="Solving ODEs…", value=0.5, {
      run_sim(
        parms     = list(),
        weight_kg = input$weight,
        tx_start  = input$tx_delay,
        use_IVIG  = input$use_IVIG,
        use_MP    = input$use_MP,
        use_PE    = input$use_PE,
        use_RTX   = input$use_RTX,
        use_CPX   = input$use_CPX,
        use_TCZ   = input$use_TCZ,
        RTX_delay = input$rtx_delay,
        TCZ_delay = input$tcz_delay
      )
    })
  }, ignoreNULL=FALSE)

  # ── All-scenario data (pre-computed) ──────────────────────────────────────
  all_scen <- reactive({
    sc <- list(
      "1. No Treatment"        = run_sim(list(), 70, 14),
      "2. IVIG + MP"           = run_sim(list(), 70, 14, TRUE, TRUE),
      "3. IVIG + MP + PE"      = run_sim(list(), 70, 14, TRUE, TRUE, TRUE),
      "4. + Rituximab"         = run_sim(list(), 70, 14, TRUE, TRUE, FALSE, TRUE),
      "5. + Cyclophosphamide"  = run_sim(list(), 70, 14, TRUE, TRUE, FALSE, FALSE, TRUE),
      "6. + Tocilizumab"       = run_sim(list(), 70, 14, TRUE, TRUE, FALSE, TRUE, FALSE, TRUE)
    )
    bind_rows(lapply(names(sc), function(nm) mutate(sc[[nm]], scenario=nm)))
  })

  col6 <- c("#E53935","#1E88E5","#43A047","#FB8C00","#8E24AA","#00ACC1")

  # ─── TAB 1 OUTPUTS ──────────────────────────────────────────────────────
  output$tbl_diagnosis <- renderDT({
    datatable(data.frame(
      Criterion = c("Probable","Definite","Supportive"),
      Description = c(
        "Rapid onset (<3mo): ≥4/6 clinical features (psychiatric, seizures, cognitive, speech, movement/rigidity, decreased LOC/autonomic)",
        "Any of the above + CSF anti-GluN1 Ab (serum insufficient alone)",
        "CSF pleocytosis or EEG abnormality; ovarian teratoma in women"
      )
    ), options=list(dom="t", paging=FALSE), rownames=FALSE)
  })

  output$plt_stages <- renderPlotly({
    stages <- data.frame(
      stage = c("Prodrome","Psychiatric","Seizures","Movement","Unconscious","Recovery"),
      days  = c(5,14,21,30,45,90),
      mRS   = c(0.5,2,3,4,5,2)
    )
    plot_ly(stages, x=~days, y=~mRS, type="scatter", mode="lines+markers",
            text=~stage, hovertemplate="%{text}: Day %{x}, mRS %{y}<extra></extra>",
            marker=list(size=10), line=list(width=2, color="#1E88E5")) %>%
      layout(xaxis=list(title="Days from Symptom Onset"),
             yaxis=list(title="mRS", range=c(0,6)),
             title="Typical Clinical Trajectory")
  })

  output$tbl_treatment <- renderDT({
    datatable(data.frame(
      Line = c("1st","1st","1st","2nd","2nd","3rd"),
      Drug = c("IVIG 2g/kg","IV Methylprednisolone","Plasmapheresis","Rituximab","Cyclophosphamide","Tocilizumab"),
      Dose = c("2g/kg over 5d","1g/d × 3-5d","5 exchanges","375mg/m² ×4wk OR 1g ×2","500-750mg/m² q4wk","8mg/kg q4wk"),
      Response = c("78%","~same as IVIG","Added ~15% benefit","79% (refractory)","~60%","Case series only")
    ), options=list(dom="t",paging=FALSE), rownames=FALSE)
  })

  output$tbl_epidemio <- renderDT({
    datatable(data.frame(
      Parameter = c("Incidence","Female:Male","Peak Age","Ovarian teratoma","Relapse rate","Full recovery","Median diagnosis delay"),
      Value = c("1.5/million/yr","4:1","18-35 yrs","40% (women 18-45)","12-25%","50% at 2 yr","~3 months")
    ), options=list(dom="t",paging=FALSE), rownames=FALSE)
  })

  # ─── TAB 2 OUTPUTS ──────────────────────────────────────────────────────
  make_plt <- function(df, x, y, title, ylab, color="#1E88E5") {
    plot_ly(df, x=~get(x), y=~get(y), type="scatter", mode="lines",
            line=list(color=color, width=2)) %>%
      layout(title=title, xaxis=list(title="Day"), yaxis=list(title=ylab))
  }

  output$plt_Ab_serum <- renderPlotly({
    df <- sim_data()
    plot_ly(df, x=~time, y=~AB_SERUM, type="scatter", mode="lines",
            line=list(color="#E53935",width=2.5)) %>%
      layout(xaxis=list(title="Day"), yaxis=list(title="Relative IgG Titer"),
             title="Serum Anti-NMDAR IgG")
  })

  output$plt_Ab_CSF <- renderPlotly({
    df <- sim_data()
    plot_ly(df, x=~time, y=~AB_CSF, type="scatter", mode="lines",
            line=list(color="#880E4F",width=2.5)) %>%
      layout(xaxis=list(title="Day"), yaxis=list(title="Relative CSF IgG"),
             title="CSF Anti-NMDAR IgG")
  })

  output$plt_PK_drugs <- renderPlotly({
    df <- sim_data()
    plot_ly() %>%
      add_trace(data=df, x=~time, y=~Cp_IVIG, name="IVIG (mg/mL)", mode="lines",
                line=list(color="#43A047",width=2)) %>%
      add_trace(data=df, x=~time, y=~Cp_MP, name="MP (mcg/mL)", mode="lines",
                line=list(color="#1E88E5",width=2)) %>%
      add_trace(data=df, x=~time, y=~Cp_RTX, name="RTX (mcg/mL)", mode="lines",
                line=list(color="#FB8C00",width=2)) %>%
      add_trace(data=df, x=~time, y=~Cp_TCZ, name="TCZ (mcg/mL)", mode="lines",
                line=list(color="#8E24AA",width=2)) %>%
      layout(xaxis=list(title="Day"), yaxis=list(title="Concentration"),
             legend=list(x=0.7,y=0.9), title="Drug PK Profiles")
  })

  output$plt_RO_RTX <- renderPlotly({
    df <- sim_data()
    plot_ly(df, x=~time, y=~RO_RTX*100, type="scatter", mode="lines",
            line=list(color="#FB8C00",width=2)) %>%
      add_segments(x=0,xend=365,y=95,yend=95,
                   line=list(color="red",dash="dash",width=1)) %>%
      layout(xaxis=list(title="Day"), yaxis=list(title="CD20 RO (%)",range=c(0,100)),
             title="Rituximab: CD20 Receptor Occupancy")
  })

  output$tbl_pk_params <- renderDT({
    datatable(data.frame(
      Drug=c("IVIG","Methylprednisolone","Rituximab","Tocilizumab","4-OH-CPX"),
      Vc_L=c("3.7","28","3.4","3.5","30"),
      Vp_L=c("25","56","4.4","2.9","–"),
      CL_Ld=c("0.21","576 L/d","0.34","0.55","96"),
      t_half=c("21d (FcRn)","2-3h","14-21d","11-13d","~4h"),
      EC50=c("12 mg/mL (FcRn)","0.25 mcg/mL","8 mcg/mL","2.5 mcg/mL","500 ng/mL"),
      Emax=c("4.5× catabolism","80% anti-inflam","98% BCD","95% IL-6 blk","92% kill")
    ), options=list(dom="t",paging=FALSE,scrollX=TRUE), rownames=FALSE)
  })

  # ─── TAB 3 OUTPUTS ──────────────────────────────────────────────────────
  output$plt_NMDAR <- renderPlotly({
    df <- sim_data()
    plot_ly(df, x=~time, y=~NMDAR_pct, type="scatter", mode="lines",
            line=list(color="#4CAF50",width=2.5)) %>%
      add_segments(x=0,xend=365,y=70,yend=70,
                   line=list(color="orange",dash="dash",width=1)) %>%
      layout(xaxis=list(title="Day"), yaxis=list(title="NMDAR (% of normal)",range=c(0,105)),
             title="Surface NMDA-R Density")
  })

  output$plt_BBB <- renderPlotly({
    df <- sim_data()
    plot_ly(df, x=~time, y=~BBB_pct, type="scatter", mode="lines",
            fill="tozeroy", fillcolor="rgba(0,150,200,0.1)",
            line=list(color="#00BCD4",width=2)) %>%
      layout(xaxis=list(title="Day"), yaxis=list(title="BBB Integrity (%)",range=c(0,105)),
             title="Blood-Brain Barrier Integrity")
  })

  output$plt_MG <- renderPlotly({
    df <- sim_data()
    plot_ly(df, x=~time, y=~MG, type="scatter", mode="lines",
            line=list(color="#9C27B0",width=2)) %>%
      add_segments(x=0,xend=365,y=1,yend=1,
                   line=list(color="grey50",dash="dash",width=1)) %>%
      layout(xaxis=list(title="Day"), yaxis=list(title="Microglia Activation (1=resting)"),
             title="Microglia Activation Index")
  })

  output$plt_IL6_GFAP <- renderPlotly({
    df <- sim_data()
    plot_ly() %>%
      add_trace(data=df, x=~time, y=~IL6_CNS, name="IL-6 (CNS)", mode="lines",
                line=list(color="#EF5350",width=2)) %>%
      add_trace(data=df, x=~time, y=~GFAP, name="GFAP (astrocytes)", mode="lines",
                line=list(color="#FF9800",width=2)) %>%
      layout(xaxis=list(title="Day"), yaxis=list(title="Relative Level"),
             legend=list(x=0.7,y=0.9), title="CNS IL-6 & GFAP")
  })

  output$plt_GLU <- renderPlotly({
    df <- sim_data()
    plot_ly(df, x=~time, y=~GLU, type="scatter", mode="lines",
            fill="tozeroy", fillcolor="rgba(255,82,82,0.1)",
            line=list(color="#EF5350",width=2)) %>%
      add_segments(x=0,xend=365,y=1.6,yend=1.6,
                   line=list(color="darkred",dash="dash",width=1)) %>%
      layout(xaxis=list(title="Day"), yaxis=list(title="Relative Glutamate"),
             title="Synaptic Glutamate (Seizure threshold = 1.6)")
  })

  output$tbl_cns_cascade <- renderDT({
    datatable(data.frame(
      Step = 1:8,
      Event = c(
        "Peripheral antigen presentation → GC B cell expansion",
        "Plasmablasts/LLPC → serum anti-NMDAR IgG",
        "BBB disruption via MMP-9/IL-6 → IgG transcytosis to CSF",
        "CSF IgG binds NR1 → bivalent crosslinking → clathrin endocytosis",
        "NMDAR internalization → surface density ↓↓",
        "PV interneuron hypofunction → Glu disinhibition ↑",
        "Dopaminergic disinhibition → psychosis; GLU excess → seizures",
        "LTP failure, hippocampal damage → cognitive/memory impairment"
      )
    ), options=list(dom="t",paging=FALSE), rownames=FALSE)
  })

  # ─── TAB 4 OUTPUTS ──────────────────────────────────────────────────────
  output$plt_mRS <- renderPlotly({
    df <- sim_data()
    plot_ly(df, x=~time, y=~mRS_est, type="scatter", mode="lines",
            line=list(color="#C62828",width=2.5)) %>%
      layout(xaxis=list(title="Day"), yaxis=list(title="mRS (0-6)",range=c(0,6.5)),
             title="mRS Estimate Over Time")
  })

  output$plt_SZ <- renderPlotly({
    df <- sim_data()
    plot_ly(df, x=~time, y=~SZ, type="scatter", mode="lines",
            fill="tozeroy", fillcolor="rgba(255,152,0,0.15)",
            line=list(color="#FF6F00",width=2)) %>%
      layout(xaxis=list(title="Day"), yaxis=list(title="Seizures/week"),
             title="Seizure Frequency")
  })

  output$plt_COG <- renderPlotly({
    df <- sim_data()
    plot_ly(df, x=~time, y=~COG_pct, type="scatter", mode="lines",
            line=list(color="#1565C0",width=2.5)) %>%
      add_segments(x=0,xend=365,y=70,yend=70,
                   line=list(color="grey",dash="dash",width=1)) %>%
      layout(xaxis=list(title="Day"), yaxis=list(title="Cognitive Index (%)",range=c(0,105)),
             title="Cognitive Function")
  })

  output$plt_PSY <- renderPlotly({
    df <- sim_data()
    plot_ly(df, x=~time, y=~PSY, type="scatter", mode="lines",
            fill="tozeroy", fillcolor="rgba(156,39,176,0.1)",
            line=list(color="#7B1FA2",width=2)) %>%
      layout(xaxis=list(title="Day"), yaxis=list(title="Psychiatric Score (0-10)"),
             title="Psychiatric Symptoms (Psychosis/Hallucinations)")
  })

  output$tbl_clinical_summary <- renderDT({
    df <- sim_data()
    tpts <- c(14,28,56,84,168,365)
    sub  <- df[df$time %in% tpts, c("time","NMDAR_pct","mRS_est","SZ","COG_pct","PSY","CRS")]
    sub  <- sub %>% mutate(across(where(is.numeric), ~round(.x, 2)))
    colnames(sub) <- c("Day","NMDAR (%)","mRS","Seizures/wk","Cognitive (%)","Psych Score","CRS")
    datatable(sub, options=list(dom="t",paging=FALSE), rownames=FALSE)
  })

  # ─── TAB 5 OUTPUTS ──────────────────────────────────────────────────────
  output$plt_scenario_bar <- renderPlotly({
    sc_df <- all_scen() %>% filter(time == 90) %>%
      group_by(scenario) %>%
      summarise(NMDAR=mean(NMDAR_pct), mRS=mean(mRS_est),
                COG=mean(COG_pct), SZ=mean(SZ), .groups="drop")
    plot_ly(sc_df, x=~scenario, y=~NMDAR, type="bar", name="NMDAR (%)",
            marker=list(color=col6)) %>%
      add_trace(y=~COG, name="Cognitive (%)", marker=list(color=col6, opacity=0.6)) %>%
      layout(barmode="group", xaxis=list(title="Scenario"),
             yaxis=list(title="% Recovery", range=c(0,100)),
             title="Day 90 Outcomes: NMDAR & Cognition")
  })

  output$plt_scen_NMDAR <- renderPlotly({
    df <- all_scen()
    scens <- unique(df$scenario)
    p <- plot_ly()
    for(i in seq_along(scens)) {
      sub <- df[df$scenario == scens[i],]
      p <- add_trace(p, data=sub, x=~time, y=~NMDAR_pct, name=scens[i],
                     type="scatter", mode="lines", line=list(color=col6[i],width=2))
    }
    p %>% layout(xaxis=list(title="Day"), yaxis=list(title="NMDAR (%)"),
                 title="NMDAR Recovery — All Scenarios")
  })

  output$plt_scen_COG <- renderPlotly({
    df <- all_scen()
    scens <- unique(df$scenario)
    p <- plot_ly()
    for(i in seq_along(scens)) {
      sub <- df[df$scenario == scens[i],]
      p <- add_trace(p, data=sub, x=~time, y=~COG_pct, name=scens[i],
                     type="scatter", mode="lines", line=list(color=col6[i],width=2))
    }
    p %>% layout(xaxis=list(title="Day"), yaxis=list(title="Cognitive (%)"),
                 title="Cognitive Function — All Scenarios")
  })

  output$tbl_scenario_compare <- renderDT({
    df <- all_scen() %>% filter(time %in% c(90,180)) %>%
      select(scenario, time, NMDAR_pct, mRS_est, COG_pct, SZ, PSY) %>%
      mutate(across(where(is.numeric), ~round(.x,2))) %>%
      rename(Scenario=scenario, Day=time,
             "NMDAR (%)"=NMDAR_pct, mRS=mRS_est,
             "Cog (%)"=COG_pct, "Sz/wk"=SZ, "Psych"=PSY)
    datatable(df, options=list(dom="t",paging=FALSE,scrollX=TRUE), rownames=FALSE)
  })

  # ─── TAB 6 OUTPUTS ──────────────────────────────────────────────────────
  output$plt_bm_IL6 <- renderPlotly({
    df <- sim_data()
    plot_ly(df, x=~time, y=~IL6_CNS, type="scatter", mode="lines",
            line=list(color="#EF5350",width=2)) %>%
      layout(xaxis=list(title="Day"), yaxis=list(title="IL-6 (relative)"),
             title="CNS IL-6 Trajectory")
  })

  output$plt_bm_GFAP <- renderPlotly({
    df <- sim_data()
    plot_ly(df, x=~time, y=~GFAP, type="scatter", mode="lines",
            line=list(color="#FF9800",width=2)) %>%
      layout(xaxis=list(title="Day"), yaxis=list(title="GFAP (relative)"),
             title="Astrocyte Reactivity (GFAP proxy)")
  })

  output$plt_bm_bcells <- renderPlotly({
    df <- sim_data()
    plot_ly() %>%
      add_trace(data=df, x=~time, y=~GCB,  name="GCB",  mode="lines",
                line=list(color="#7986CB",width=2)) %>%
      add_trace(data=df, x=~time, y=~PB,   name="Plasmablast", mode="lines",
                line=list(color="#5C6BC0",width=2)) %>%
      add_trace(data=df, x=~time, y=~LLPC, name="LLPC", mode="lines",
                line=list(color="#3949AB",width=2)) %>%
      add_trace(data=df, x=~time, y=~MB,   name="Memory B", mode="lines",
                line=list(color="#9FA8DA",width=2)) %>%
      layout(xaxis=list(title="Day"), yaxis=list(title="Cell count (arb)"),
             legend=list(x=0.7,y=0.9), title="B Cell Dynamics")
  })

  output$plt_DR_RTX <- renderPlotly({
    doses <- seq(0, 1500, by=150)
    nmdar90 <- sapply(doses, function(d) {
      df <- run_sim(list(), 70, 14, TRUE, TRUE, FALSE, TRUE,
                    RTX_delay=30)
      val <- df$NMDAR_pct[df$time == 90]
      if(length(val)==0) return(NA)
      val[1]
    })
    plot_ly(x=doses, y=nmdar90, type="scatter", mode="lines+markers",
            line=list(color="#FB8C00",width=2),
            marker=list(size=7, color="#FB8C00")) %>%
      layout(xaxis=list(title="Rituximab Single Dose (mg)"),
             yaxis=list(title="NMDAR (%) at Day 90"),
             title="RTX Dose–Response (NMDAR Recovery at Day 90)")
  })

  output$tbl_biomarkers <- renderDT({
    datatable(data.frame(
      Biomarker = c(
        "CSF Anti-NR1 IgG","Serum Anti-NR1 IgG","CSF IgG Index",
        "CSF WBC","CSF CXCL13","CSF NfL (neurofilament)",
        "CSF GFAP","CSF IL-6","EEG Delta Brush","MRI FLAIR (hippocampus)"),
      Normal = c("<1:1 (titer)","Often seroneg","<0.7","<5/μL","<24 pg/mL",
                 "65-1100 pg/mL","130-3200 pg/mL","<1.6 pg/mL","Absent","Normal"),
      AIE_Active = c(">1:160","1:20–1:640","Often >1.0","5-1000/μL",
                     ">50 pg/mL","↑↑ >2000 pg/mL","↑ >5000 pg/mL","↑↑",
                     "30% of cases","35% FLAIR hyperintensity"),
      Clinical_Utility = c(
        "Definitive diagnosis; CSF preferred","Screening; low sens","Intrathecal synth",
        "Supports diagnosis","GC reaction; predicts rituximab response","Neuroaxonal damage",
        "Astrocyte damage; ↑ in severe disease","Inflammation; TCZ target",
        "Pathognomonic (slow delta + bursts)","Associates with hippocampal atrophy")
    ), options=list(dom="t",paging=FALSE,scrollX=TRUE), rownames=FALSE)
  })
}

# ─────────────────────────────────────────────────────────────────────────────
shinyApp(ui = ui, server = server)
