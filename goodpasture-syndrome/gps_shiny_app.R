###############################################################################
# Goodpasture Syndrome QSP — Interactive Shiny Dashboard
#
# Tabs (7):
#   1. Patient Profile      — disease overview, risk factors, diagnosis criteria
#   2. Drug PK              — concentration-time per drug (plex/CY/pred/RTX/ava)
#   3. Antibody & Complement — anti-GBM titer, C5a, B cell dynamics
#   4. Renal Endpoints      — GFR, creatinine, proteinuria, dialysis risk
#   5. Pulmonary Endpoints  — DLCO, lung damage, DAH, hemoptysis
#   6. Scenario Comparison  — all 6 treatment scenarios head-to-head
#   7. References           — clinical trials, PK/PD papers, outcome data
###############################################################################

library(shiny)
library(shinydashboard)
library(dplyr)
library(ggplot2)
library(tidyr)
library(plotly)
library(DT)

# ── Inline ODE model (pure R, no mrgsolve required for Shiny standalone) ─────
#    Euler integration for speed; replace with mrgsolve for higher accuracy

gps_ode <- function(t, state, parms) {
  with(as.list(c(state, parms)), {
    # ── Drug PK ──────────────────────────────────────────────────────────────
    # Cyclophosphamide IV bolus → active metabolite 4-OH-CY
    dCY_C  <- -kCY * CY_C
    dCY_M  <-  kCY * CY_C - kCY_M * CY_M

    # Prednisolone oral
    dPRED_C <- -kPRED * PRED_C

    # Rituximab IV (2-cpt approximated as 1-cpt here)
    dRTX_C  <- -kRTX * RTX_C

    # Avacopan oral BID (rapid, short t½)
    dAVA_C  <- -kAVA * AVA_C

    # ── Drug effects ─────────────────────────────────────────────────────────
    # Cyclophosphamide: suppresses plasma cells (via active metabolite)
    E_CY   <- Emax_CY  * CY_M  / (IC50_CY  + CY_M)
    # Prednisolone: broad immunosuppression (reduces anti-GBM production)
    E_PRED <- Emax_PRED* PRED_C/ (IC50_PRED + PRED_C)
    # Rituximab: depletes B cells (Emax ~95% depletion)
    E_RTX  <- Emax_RTX * RTX_C / (IC50_RTX  + RTX_C)
    # Avacopan: blocks C5a signaling
    E_AVA  <- Emax_AVA * AVA_C / (IC50_AVA  + AVA_C)

    # ── Anti-GBM antibody (IgG) ───────────────────────────────────────────
    # Produced by plasma cells; plasmapheresis removes ~60% per session
    PlexEffect <- PLEX_rate * AntiGBM  # fractional removal during active plex
    dAntiGBM <- kin_Ab * Plasma_cells * (1 - E_CY) * (1 - E_PRED)
                - kout_Ab * AntiGBM
                - PlexEffect

    # ── B cells ───────────────────────────────────────────────────────────
    # Depleted by rituximab; recover naturally (t½ recovery ~6mo)
    dB_cells <- kin_B * (1 - E_RTX) - kout_B * B_cells

    # ── Plasma cells ──────────────────────────────────────────────────────
    # Fed by B cell differentiation; suppressed by CY + pred
    dPlasma_cells <- k_Bplasma * B_cells * (1 - E_RTX) * (1 - E_CY) * (1 - E_PRED)
                    - kout_plasma * Plasma_cells

    # ── C5a complement ────────────────────────────────────────────────────
    # Activated by IgG deposition; blocked by avacopan
    dC5a <- kin_C5a * (AntiGBM / AntiGBM0)  - kout_C5a * C5a * (1 + E_AVA)

    # ── Neutrophil infiltration (kidney) ─────────────────────────────────
    # Recruited by C5a and FcγR ligation
    dNeut <- kin_Neut * C5a * (AntiGBM / AntiGBM0) - kout_Neut * Neut

    # ── GBM damage score (0=normal, 100=complete destruction) ────────────
    Repair_GBM <- k_GBM_repair * (1 + E_PRED)
    dGBM_damage <- k_GBM_inflam * Neut * (AntiGBM / AntiGBM0)
                   - Repair_GBM * GBM_damage

    # ── GFR (mL/min/1.73m²) ──────────────────────────────────────────────
    GFR_loss_rate <- k_GFR_loss * GBM_damage / 100
    dGFR_c <- -GFR_loss_rate * GFR_c * (1 - pmin(0.8, E_PRED + E_CY))

    # ── Lung damage (0=normal, 100=complete DAH) ──────────────────────────
    Repair_lung <- k_Lung_repair * (1 + E_PRED)
    dLung_damage <- k_Lung_inflam * Neut * (AntiGBM / AntiGBM0)
                    - Repair_lung * Lung_damage

    # ── DLCO (% predicted, normal=100%) ──────────────────────────────────
    dDLCO_c <- -k_DLCO_loss * Lung_damage / 100 * DLCO_c
               + k_DLCO_repair * (100 - DLCO_c) * (1 - Lung_damage/100)

    # ── Proteinuria (g/day) ───────────────────────────────────────────────
    dProteinuria_c <- k_prot_inflam * GBM_damage / 100
                    - k_prot_repair * Proteinuria_c * (1 - GBM_damage/100)

    # ── CRP (mg/L) ────────────────────────────────────────────────────────
    dCRP_c <- k_CRP_prod * (Neut / 100) - k_CRP_clear * CRP_c

    list(c(dCY_C, dCY_M, dPRED_C, dRTX_C, dAVA_C,
           dAntiGBM, dB_cells, dPlasma_cells, dC5a, dNeut,
           dGBM_damage, dGFR_c, dLung_damage, dDLCO_c,
           dProteinuria_c, dCRP_c))
  })
}

# Default parameters
default_parms <- list(
  # Drug PK rate constants (1/day)
  kCY   = 2.4,    # CY t½~7h → 0.693/7*24~2.4/day
  kCY_M = 3.6,    # 4-OH-CY t½~4.6h
  kPRED = 6.6,    # Prednisolone t½~2.5h
  kRTX  = 0.033,  # Rituximab t½~21days = 0.693/21
  kAVA  = 4.8,    # Avacopan t½~3.5h

  # Drug Emax and IC50
  Emax_CY=0.80, IC50_CY=2.0,       # CY metabolite (μg/mL)
  Emax_PRED=0.70, IC50_PRED=50.0,  # Prednisolone (ng/mL)
  Emax_RTX=0.95, IC50_RTX=0.5,     # Rituximab (μg/mL)
  Emax_AVA=0.85, IC50_AVA=80.0,    # Avacopan (ng/mL)

  # Disease dynamics
  AntiGBM0 = 200,   # Baseline anti-GBM titer
  kin_Ab = 4.0,     # AU/mL/day - antibody production rate (plasma cells)
  kout_Ab = 0.02,   # 1/day - natural clearance (IgG t½~21 days)
  PLEX_rate = 0.0,  # Plasmapheresis removal rate (set by scenario)

  kin_B = 5.0,      # B cell replenishment rate
  kout_B = 0.05,    # 1/day - natural B cell turnover
  k_Bplasma = 0.1,  # B→plasma cell differentiation rate
  kout_plasma = 0.03, # Plasma cell death rate

  kin_C5a = 0.5,    # C5a production from IgG deposition
  kout_C5a = 0.3,   # C5a clearance
  kin_Neut = 2.0,   # Neutrophil recruitment
  kout_Neut = 1.5,  # Neutrophil clearance

  k_GBM_inflam = 0.08,  # GBM damage rate
  k_GBM_repair = 0.02,  # GBM repair rate

  k_GFR_loss = 0.03,    # GFR decline per unit GBM damage
  k_Lung_inflam = 0.06, # Lung damage rate
  k_Lung_repair = 0.04, # Lung repair rate
  k_DLCO_loss = 0.04,   # DLCO decline rate
  k_DLCO_repair = 0.02, # DLCO recovery rate

  k_prot_inflam = 0.5,  # Proteinuria increase
  k_prot_repair = 0.15, # Proteinuria resolution rate
  k_CRP_prod = 2.0,     # CRP production from neutrophils
  k_CRP_clear = 0.15    # CRP clearance
)

# Initial state
initial_state_default <- c(
  CY_C=0, CY_M=0, PRED_C=0, RTX_C=0, AVA_C=0,
  AntiGBM=200, B_cells=100, Plasma_cells=100, C5a=8, Neut=80,
  GBM_damage=40, GFR_c=25, Lung_damage=35, DLCO_c=45,
  Proteinuria_c=3.5, CRP_c=45
)

# Simple Euler simulator (avoids mrgsolve dependency for Shiny)
simulate_gps <- function(parms, init_state, times, scenario,
                         anti_gbm0=200, gfr0=25, dlco0=45) {
  # Adjust initial state
  init_state["AntiGBM"] <- anti_gbm0
  init_state["GFR_c"]   <- gfr0
  init_state["DLCO_c"]  <- dlco0
  parms["AntiGBM0"]     <- anti_gbm0

  results <- matrix(0, nrow=length(times), ncol=length(init_state)+1)
  colnames(results) <- c("time", names(init_state))
  state <- init_state

  dt <- times[2] - times[1]

  for (i in seq_along(times)) {
    t <- times[i]
    results[i,] <- c(t, state)

    # Apply drug doses at specific timepoints (simplified as continuous infusion proxies)
    p <- parms
    # Scenario-specific drug forcing
    if (scenario == 1) {
      # Plex (14 sessions over 14 days), CY 2mg/kg/day, Pred 1mg/kg/day
      if (t <= 14) p["PLEX_rate"] <- 0.10
      if (t <= 90) {
        state["CY_C"]   <- pmax(state["CY_C"]   + 140/24 * dt, 0)  # ~2mg/kg day→mg/L proxy
        state["PRED_C"] <- pmax(state["PRED_C"] + 70/24 * dt, 0)   # ~1mg/kg
      }
    } else if (scenario == 2) {
      if (t <= 14) p["PLEX_rate"] <- 0.10
      # RTX: 4 weekly doses (wk0,1,2,3 = days 0,7,14,21)
      if (t %in% c(0, 7, 14, 21)) state["RTX_C"] <- state["RTX_C"] + 375*0.00001 # scaled
      if (t <= 90) state["PRED_C"] <- pmax(state["PRED_C"] + 30/24 * dt, 0)
    } else if (scenario == 3) {
      if (t <= 14) p["PLEX_rate"] <- 0.10
      # Avacopan 30mg BID: 60mg/day
      if (t <= 180) state["AVA_C"] <- pmax(state["AVA_C"] + 15 * dt, 0)  # scaled ng/mL proxy
      if (t <= 90) state["PRED_C"] <- pmax(state["PRED_C"] + 20/24 * dt, 0)
    } else if (scenario == 4) {
      # No plex, CY + pred
      if (t <= 90) {
        state["CY_C"]   <- pmax(state["CY_C"]   + 140/24 * dt, 0)
        state["PRED_C"] <- pmax(state["PRED_C"] + 70/24 * dt, 0)
      }
    } else if (scenario == 5) {
      # Pred only
      if (t <= 90) state["PRED_C"] <- pmax(state["PRED_C"] + 70/24 * dt, 0)
    }

    # Clamp plasmapheresis
    p_mod <- p
    if (!("PLEX_rate" %in% names(p))) p_mod[["PLEX_rate"]] <- 0

    # Euler step
    derivs <- tryCatch(
      gps_ode(t, state, p_mod)[[1]],
      error = function(e) rep(0, length(state))
    )

    state <- pmax(0, state + derivs * dt)
    # Clamp physiological limits
    state["GBM_damage"]   <- pmin(100, state["GBM_damage"])
    state["Lung_damage"]  <- pmin(100, state["Lung_damage"])
    state["DLCO_c"]       <- pmin(100, pmax(0, state["DLCO_c"]))
    state["GFR_c"]        <- pmax(0, state["GFR_c"])
    state["B_cells"]      <- pmax(0, pmin(200, state["B_cells"]))
    state["Plasma_cells"] <- pmax(0, pmin(200, state["Plasma_cells"]))
  }
  as.data.frame(results)
}

# Scenario names and colors
scenario_names  <- c("0"="No Treatment","1"="Plex+CY+Pred","2"="Plex+RTX+Pred",
                      "3"="Plex+Avacopan+Pred","4"="CY+Pred (no Plex)","5"="Pred only")
scenario_colors <- c(
  "No Treatment"       = "#999999",
  "Plex+CY+Pred"       = "#E41A1C",
  "Plex+RTX+Pred"      = "#377EB8",
  "Plex+Avacopan+Pred" = "#4DAF4A",
  "CY+Pred (no Plex)"  = "#984EA3",
  "Pred only"          = "#FF7F00"
)

times_default <- seq(0, 365, by=1)

# ── UI ────────────────────────────────────────────────────────────────────────
ui <- dashboardPage(
  skin = "red",

  dashboardHeader(
    title = "Goodpasture Syndrome QSP",
    titleWidth = 300
  ),

  dashboardSidebar(
    width = 270,
    sidebarMenu(
      menuItem("Patient Profile",      tabName="tab_patient",    icon=icon("user")),
      menuItem("Drug PK",              tabName="tab_pk",         icon=icon("pills")),
      menuItem("Antibody & Complement",tabName="tab_antibody",   icon=icon("dna")),
      menuItem("Renal Endpoints",      tabName="tab_renal",      icon=icon("kidneys")),
      menuItem("Pulmonary Endpoints",  tabName="tab_pulmonary",  icon=icon("lungs")),
      menuItem("Scenario Comparison",  tabName="tab_comparison", icon=icon("balance-scale")),
      menuItem("References",           tabName="tab_refs",       icon=icon("book-open"))
    ),
    hr(),
    h4("Patient Parameters", style="padding-left:15px;color:#aaa;"),
    sliderInput("anti_gbm_init", "Initial Anti-GBM Titer (AU/mL)", 50, 500, 200, 25),
    sliderInput("gfr_init",      "Initial GFR (mL/min/1.73m²)",    5,  60,  25, 5),
    sliderInput("dlco_init",     "Initial DLCO (% predicted)",      15, 90,  45, 5),
    sliderInput("bw_kg",         "Body Weight (kg)",                 40, 120, 70, 5),
    hr(),
    h4("Simulation", style="padding-left:15px;color:#aaa;"),
    sliderInput("sim_days", "Duration (days)", 56, 730, 365, 28),
    selectInput("scenario_pk", "Drug Scenario (PK/Ab/Renal/Lung tabs)",
                choices=list("No Treatment"="0","Plex+CY+Pred"="1","Plex+RTX+Pred"="2",
                             "Plex+Avacopan+Pred"="3","CY+Pred"="4","Pred only"="5"),
                selected="1"),
    actionButton("run_btn","Run Simulation", icon=icon("play"),
                 style="margin:10px 15px;width:240px;background:#c0392b;color:white;"),
    br(),
    checkboxGroupInput("scenarios_compare","Compare Scenarios:",
                       choices=list("No Treatment"="0","Plex+CY+Pred"="1",
                                    "Plex+RTX+Pred"="2","Plex+Avacopan+Pred"="3",
                                    "CY+Pred"="4","Pred only"="5"),
                       selected=c("0","1","2")),
    actionButton("compare_btn","Run Comparison", icon=icon("sync"),
                 style="margin:5px 15px;width:240px;background:#2980b9;color:white;"),
    br(),
    small(style="padding-left:15px;color:#888;","Goodpasture Syndrome QSP v1.0 | CCR 2026-06-19")
  ),

  dashboardBody(
    tabItems(

      # ── Tab 1: Patient Profile ─────────────────────────────────────────────
      tabItem(tabName="tab_patient",
        fluidRow(
          box(title="Disease Overview", width=6, status="danger", solidHeader=TRUE,
            h4("Goodpasture Syndrome (Anti-GBM Disease)"),
            p("A rare autoimmune pulmonary-renal syndrome (incidence ~1/million/year) caused by
               autoantibodies targeting the α3 chain of type IV collagen (COL4A3 NC1 domain)
               in the glomerular (GBM) and alveolar (ABM) basement membranes."),
            tags$ul(
              tags$li("Linear IgG deposition on GBM → rapidly progressive glomerulonephritis (RPGN)"),
              tags$li("Alveolar capillary BM involvement → diffuse alveolar hemorrhage (DAH)"),
              tags$li("~30% ANCA double-positive (usually MPO-ANCA)"),
              tags$li("Bimodal age distribution: 20–30y (male-dominated) and 60–70y (female-dominated)"),
              tags$li("Key triggers: smoking, hydrocarbons, viral infections, lithium exposure")
            ),
            hr(),
            h5("Current Patient Profile:"),
            tableOutput("patient_summary_tbl")
          ),
          box(title="Diagnosis & Classification", width=6, status="warning", solidHeader=TRUE,
            h5("Diagnostic Criteria:"),
            tableOutput("diagnosis_tbl"),
            hr(),
            h5("Baseline Disease Severity:"),
            plotOutput("severity_gauge", height="200px")
          )
        ),
        fluidRow(
          box(title="Treatment Algorithm", width=6, status="primary",
            h5("Standard Management:"),
            tags$ol(
              tags$li("Immediate plasmapheresis (daily ×14 sessions, 4L exchange)"),
              tags$li("Cyclophosphamide 2mg/kg/day PO or IV pulsed"),
              tags$li("Prednisolone 1mg/kg/day → taper"),
              tags$li("Anti-GBM titer monitoring every 2 weeks"),
              tags$li("Rituximab if CY contraindicated or relapse"),
              tags$li("Avacopan (C5aR-inhibitor) — investigational for anti-GBM")
            )
          ),
          box(title="Prognostic Factors", width=6, status="info",
            tableOutput("prognosis_tbl")
          )
        )
      ),

      # ── Tab 2: Drug PK ─────────────────────────────────────────────────────
      tabItem(tabName="tab_pk",
        fluidRow(
          box(title="Drug PK Parameters", width=4, status="primary", solidHeader=TRUE,
            selectInput("pk_drug_select","Show Drug PK:", choices=list(
              "Cyclophosphamide (CY)"="CY","Prednisolone (Pred)"="PRED",
              "Rituximab (RTX)"="RTX","Avacopan (Ava)"="AVA")),
            tableOutput("pk_param_detail")
          ),
          box(title="Concentration-Time Profile", width=8, status="primary", solidHeader=TRUE,
            plotlyOutput("pk_conc_plot", height="350px"),
            br(),
            plotlyOutput("drug_effect_plot_pk", height="250px")
          )
        ),
        fluidRow(
          box(title="Plasmapheresis Schedule & Effect", width=12, status="warning",
            plotlyOutput("plex_schedule_plot", height="250px")
          )
        )
      ),

      # ── Tab 3: Antibody & Complement ─────────────────────────────────────
      tabItem(tabName="tab_antibody",
        fluidRow(
          valueBoxOutput("vbox_antiGBM",  width=3),
          valueBoxOutput("vbox_Bcells",   width=3),
          valueBoxOutput("vbox_C5a",      width=3),
          valueBoxOutput("vbox_Neut",     width=3)
        ),
        fluidRow(
          box(title="Anti-GBM Antibody Titer", width=6, status="danger", solidHeader=TRUE,
            plotlyOutput("ab_titer_plot", height="300px"),
            p(style="font-size:0.85em;", "Target: Anti-GBM negative (<5 AU/mL) by 3 months")
          ),
          box(title="B Cells & Plasma Cells", width=6, status="warning", solidHeader=TRUE,
            plotlyOutput("bcell_plot", height="300px")
          )
        ),
        fluidRow(
          box(title="Complement C5a & Neutrophil Infiltration", width=6, status="info",
            plotlyOutput("c5a_neut_plot", height="280px")
          ),
          box(title="CRP over Time", width=6, status="info",
            plotlyOutput("crp_ab_plot", height="280px")
          )
        )
      ),

      # ── Tab 4: Renal Endpoints ────────────────────────────────────────────
      tabItem(tabName="tab_renal",
        fluidRow(
          valueBoxOutput("vbox_gfr",       width=3),
          valueBoxOutput("vbox_proteinuria",width=3),
          valueBoxOutput("vbox_gbmdamage", width=3),
          valueBoxOutput("vbox_dialysis",  width=3)
        ),
        fluidRow(
          box(title="GFR Trajectory", width=6, status="primary", solidHeader=TRUE,
            plotlyOutput("gfr_plot", height="300px"),
            p(style="font-size:0.85em;", "GFR 60+ = normal | 30-60 = stage 3 CKD | <15 = ESRD risk")
          ),
          box(title="GBM Damage Score", width=6, status="danger", solidHeader=TRUE,
            plotlyOutput("gbm_damage_plot", height="300px")
          )
        ),
        fluidRow(
          box(title="Proteinuria (g/day)", width=6, status="warning",
            plotlyOutput("proteinuria_plot", height="280px")
          ),
          box(title="Renal Outcomes Summary", width=6, status="info",
            DT::dataTableOutput("renal_summary_tbl")
          )
        )
      ),

      # ── Tab 5: Pulmonary Endpoints ────────────────────────────────────────
      tabItem(tabName="tab_pulmonary",
        fluidRow(
          valueBoxOutput("vbox_dlco",      width=3),
          valueBoxOutput("vbox_lungdamage",width=3),
          valueBoxOutput("vbox_resp",      width=3),
          valueBoxOutput("vbox_pulm_remission",width=3)
        ),
        fluidRow(
          box(title="DLCO (% Predicted)", width=6, status="primary", solidHeader=TRUE,
            plotlyOutput("dlco_plot", height="300px"),
            p(style="font-size:0.85em;", "DLCO >70% = normal | 40-70% = moderate impairment | <40% = severe")
          ),
          box(title="Alveolar/Lung Damage Score", width=6, status="danger", solidHeader=TRUE,
            plotlyOutput("lung_damage_plot", height="300px")
          )
        ),
        fluidRow(
          box(title="Pulmonary Remission Timeline", width=12, status="success",
            plotlyOutput("pulm_remission_timeline", height="250px")
          )
        )
      ),

      # ── Tab 6: Scenario Comparison ────────────────────────────────────────
      tabItem(tabName="tab_comparison",
        fluidRow(
          box(title="Anti-GBM Titer — All Scenarios", width=6, status="danger", solidHeader=TRUE,
            plotlyOutput("comp_antiGBM", height="300px")
          ),
          box(title="GFR — All Scenarios", width=6, status="primary", solidHeader=TRUE,
            plotlyOutput("comp_gfr", height="300px")
          )
        ),
        fluidRow(
          box(title="DLCO — All Scenarios", width=4, status="info", solidHeader=TRUE,
            plotlyOutput("comp_dlco", height="280px")
          ),
          box(title="GBM Damage — All Scenarios", width=4, status="warning",
            plotlyOutput("comp_gbmdamage", height="280px")
          ),
          box(title="B Cell Depletion", width=4, status="primary",
            plotlyOutput("comp_bcells", height="280px")
          )
        ),
        fluidRow(
          box(title="Outcomes at Key Timepoints", width=12, status="success",
            DT::dataTableOutput("comp_outcomes_tbl")
          )
        )
      ),

      # ── Tab 7: References ─────────────────────────────────────────────────
      tabItem(tabName="tab_refs",
        fluidRow(
          box(title="Clinical Trials & Cohorts", width=6, status="primary", solidHeader=TRUE,
            tableOutput("trials_tbl")
          ),
          box(title="PK/PD Modeling", width=6, status="info", solidHeader=TRUE,
            tableOutput("pkpd_refs_tbl")
          )
        ),
        fluidRow(
          box(title="Key Disease Biology References", width=12, status="warning",
            tableOutput("biology_refs_tbl")
          )
        )
      )
    )
  )
)

# ── SERVER ────────────────────────────────────────────────────────────────────
server <- function(input, output, session) {

  # Single simulation
  single_sim <- eventReactive(input$run_btn, {
    parms <- default_parms
    parms$AntiGBM0 <- input$anti_gbm_init
    init  <- initial_state_default
    init["AntiGBM"] <- input$anti_gbm_init
    init["GFR_c"]   <- input$gfr_init
    init["DLCO_c"]  <- input$dlco_init

    times <- seq(0, input$sim_days, by=1)
    sim <- simulate_gps(parms, init, times, as.integer(input$scenario_pk),
                        input$anti_gbm_init, input$gfr_init, input$dlco_init)
    sim$Scenario <- scenario_names[input$scenario_pk]
    sim$time_wk  <- sim$time / 7
    sim
  })

  # All scenario comparison
  all_sims <- eventReactive(c(input$run_btn, input$compare_btn), {
    parms <- default_parms
    parms$AntiGBM0 <- input$anti_gbm_init
    times <- seq(0, input$sim_days, by=1)
    scs <- input$scenarios_compare
    bind_rows(lapply(scs, function(sc) {
      init <- initial_state_default
      init["AntiGBM"] <- input$anti_gbm_init
      init["GFR_c"]   <- input$gfr_init
      init["DLCO_c"]  <- input$dlco_init
      sim <- simulate_gps(parms, init, times, as.integer(sc),
                          input$anti_gbm_init, input$gfr_init, input$dlco_init)
      sim$Scenario <- scenario_names[sc]
      sim$time_wk  <- sim$time / 7
      sim
    }))
  })

  # ── Tab 1 ─────────────────────────────────────────────────────────────────
  output$patient_summary_tbl <- renderTable({
    data.frame(
      Parameter = c("Initial Anti-GBM titer","Initial GFR","Initial DLCO","Body Weight"),
      Value     = c(paste(input$anti_gbm_init,"AU/mL"),paste(input$gfr_init,"mL/min/1.73m²"),
                    paste(input$dlco_init,"% predicted"),paste(input$bw_kg,"kg")),
      Status    = c(ifelse(input$anti_gbm_init>100,"High (+++)","Low (+)"),
                    ifelse(input$gfr_init<15,"ESRD risk",ifelse(input$gfr_init<30,"Severe","Moderate")),
                    ifelse(input$dlco_init<40,"Severe",ifelse(input$dlco_init<60,"Moderate","Mild")),
                    "")
    )
  })

  output$diagnosis_tbl <- renderTable({
    data.frame(
      Criterion = c("Anti-GBM antibody (serum ELISA)","Renal biopsy: linear IgG on GBM",
                    "Hemoptysis / DAH on BAL","Rapidly rising creatinine","ANCA (MPO) co-positivity"),
      Finding   = c("> 5–10 EU (positive)","Pathognomonic for anti-GBM disease",
                    "Hemosiderin-laden macrophages","RPGN: ≥50% rise in <3 months","~30% of cases")
    )
  })

  output$prognosis_tbl <- renderTable({
    data.frame(
      Factor = c("Dialysis at presentation","Anti-GBM titer","Crescent %","Treatment timing","ANCA co-positive"),
      `Renal prognosis` = c("Very poor (ESRD ~90%)","Higher titer = worse",">85% crescents = worst","Early treatment = better","Better response to Rx")
    )
  })

  output$severity_gauge <- renderPlot({
    df <- data.frame(
      Metric = c("Anti-GBM\n(÷5)","GFR decline\n(75-GFR)","DLCO deficit\n(100-DLCO)"),
      Value  = c(input$anti_gbm_init/5, 75-input$gfr_init, 100-input$dlco_init),
      Limit  = c(100, 70, 85)
    )
    ggplot(df, aes(Metric, Value, fill=Metric)) +
      geom_col(width=0.5) +
      geom_errorbar(aes(ymin=0, ymax=Limit), width=0, color="grey40") +
      coord_flip() +
      scale_fill_manual(values=c("#E41A1C","#984EA3","#377EB8")) +
      labs(title="Disease Burden (scaled)", x=NULL, y="Scaled severity") +
      theme_minimal(base_size=10) + theme(legend.position="none")
  })

  # ── Tab 2: PK ─────────────────────────────────────────────────────────────
  output$pk_param_detail <- renderTable({
    pk_df <- list(
      CY  = data.frame(Param=c("Route","Dose","Vd","CL","t½","Active metabolite"),
                       Value=c("IV/PO","2mg/kg/day","45L","6L/hr","~7hr","4-OH-CY (alkylating)")),
      PRED= data.frame(Param=c("Route","Dose","F","Vd","CL","t½"),
                       Value=c("PO","1mg/kg/day","82%","35L","8L/hr","~2.5hr")),
      RTX = data.frame(Param=c("Route","Dose","Vd","CL","t½","Mechanism"),
                       Value=c("IV","375mg/m² ×4","4.5L","0.23L/day","~21 days","anti-CD20 ADCC/CDC")),
      AVA = data.frame(Param=c("Route","Dose","F","Vd","CL","t½"),
                       Value=c("PO","30mg BID","69%","70L","15L/hr","~3.5hr"))
    )
    pk_df[[input$pk_drug_select]]
  })

  output$pk_conc_plot <- renderPlotly({
    sim <- single_sim()
    if (nrow(sim) == 0) return(plotly_empty())
    col_map <- c(CY="CY_M", PRED="PRED_C", RTX="RTX_C", AVA="AVA_C")
    col <- col_map[[input$pk_drug_select]]
    if (!col %in% names(sim)) return(plotly_empty())
    p <- ggplot(sim, aes_string("time_wk", col)) +
      geom_line(color="#E41A1C", size=1.2) +
      labs(title=paste(input$pk_drug_select,"Concentration"), x="Time (weeks)", y="Cp (model units)") +
      theme_bw(base_size=11)
    ggplotly(p)
  })

  output$drug_effect_plot_pk <- renderPlotly({
    sim <- single_sim()
    if (nrow(sim) == 0) return(plotly_empty())
    # Compute % suppression of anti-GBM vs no-treatment
    sup_frac <- pmax(0, (input$anti_gbm_init - sim$AntiGBM) / input$anti_gbm_init)
    df <- data.frame(time_wk=sim$time_wk, Suppression=sup_frac)
    p <- ggplot(df, aes(time_wk, Suppression)) +
      geom_line(color="steelblue", size=1.2) +
      scale_y_continuous(labels=scales::percent, limits=c(0,1)) +
      labs(title="Anti-GBM Suppression (%)", x="Time (weeks)", y="Suppression (%)") +
      theme_bw(base_size=11)
    ggplotly(p)
  })

  output$plex_schedule_plot <- renderPlotly({
    sc <- as.integer(input$scenario_pk)
    if (sc %in% c(0, 4, 5)) {
      p <- ggplot() + annotate("text",x=0,y=0,label="No plasmapheresis in this scenario",size=5) +
        theme_void()
    } else {
      days <- 1:14
      df <- data.frame(Day=days, Effect_pct=60)
      p <- ggplot(df, aes(Day, Effect_pct)) +
        geom_col(fill="#E41A1C", width=0.6) +
        labs(title="Plasmapheresis Sessions (Days 1–14; ~60% IgG removal per session)",
             x="Day", y="Anti-GBM removal (%)") +
        theme_bw(base_size=11)
    }
    ggplotly(p)
  })

  # ── Tab 3: Antibody ────────────────────────────────────────────────────────
  output$vbox_antiGBM <- renderValueBox({
    sim <- single_sim()
    val <- if(nrow(sim)>0) round(tail(sim$AntiGBM,1)) else "—"
    col <- if(is.numeric(val)&&val<10) "green" else "red"
    valueBox(paste(val,"AU/mL"), "Anti-GBM Titer", icon=icon("vial"), color=col)
  })
  output$vbox_Bcells <- renderValueBox({
    sim <- single_sim()
    val <- if(nrow(sim)>0) round(tail(sim$B_cells,1)) else "—"
    col <- if(is.numeric(val)&&val<30) "green" else "yellow"
    valueBox(paste(val,"%"), "B Cell Level", icon=icon("circle"), color=col)
  })
  output$vbox_C5a <- renderValueBox({
    sim <- single_sim()
    val <- if(nrow(sim)>0) round(tail(sim$C5a,2)) else "—"
    col <- if(is.numeric(val)&&val<2) "green" else "yellow"
    valueBox(paste(val,"ng/mL"), "C5a Complement", icon=icon("shield-alt"), color=col)
  })
  output$vbox_Neut <- renderValueBox({
    sim <- single_sim()
    val <- if(nrow(sim)>0) round(tail(sim$Neut,1)) else "—"
    col <- if(is.numeric(val)&&val<30) "green" else "red"
    valueBox(val, "Neutrophil Infiltration", icon=icon("bug"), color=col)
  })

  ab_plot <- function(sim, var, ylab, title, hline=NULL, color="#E41A1C") {
    if(nrow(sim)==0) return(plotly_empty())
    p <- ggplot(sim, aes_string("time_wk", var)) +
      geom_line(color=color, size=1.2)
    if(!is.null(hline)) p <- p + geom_hline(yintercept=hline$y, linetype=hline$lty, color=hline$col)
    p <- p + labs(title=title, x="Time (weeks)", y=ylab) + theme_bw(base_size=10)
    ggplotly(p)
  }

  output$ab_titer_plot <- renderPlotly({
    ab_plot(single_sim(), "AntiGBM", "Anti-GBM (AU/mL)", "Anti-GBM Antibody Titer",
            hline=list(y=10,lty="dashed",col="darkgreen"))
  })
  output$bcell_plot <- renderPlotly({
    sim <- single_sim()
    if(nrow(sim)==0) return(plotly_empty())
    df <- sim %>% select(time_wk, B_cells, Plasma_cells) %>%
      pivot_longer(-time_wk, names_to="Cell", values_to="Pct")
    p <- ggplot(df, aes(time_wk, Pct, color=Cell)) +
      geom_line(size=1.1) +
      scale_color_manual(values=c(B_cells="#377EB8",Plasma_cells="#E41A1C")) +
      geom_hline(yintercept=100, linetype="dashed", color="grey50") +
      labs(title="B Cells & Plasma Cells (% of normal)", x="Time (weeks)", y="% Normal") +
      theme_bw(base_size=10)
    ggplotly(p)
  })
  output$c5a_neut_plot <- renderPlotly({
    sim <- single_sim()
    if(nrow(sim)==0) return(plotly_empty())
    df <- sim %>% select(time_wk, C5a, Neut) %>%
      pivot_longer(-time_wk, names_to="Var", values_to="Value")
    p <- ggplot(df, aes(time_wk, Value, color=Var)) +
      geom_line(size=1.1) +
      scale_color_manual(values=c(C5a="#984EA3",Neut="#FF7F00")) +
      labs(title="C5a & Neutrophil Infiltration", x="Time (weeks)", y="Value") +
      theme_bw(base_size=10)
    ggplotly(p)
  })
  output$crp_ab_plot <- renderPlotly({
    ab_plot(single_sim(), "CRP_c", "CRP (mg/L)", "C-Reactive Protein",
            hline=list(y=5,lty="dashed",col="darkgreen"), color="#FF7F00")
  })

  # ── Tab 4: Renal ─────────────────────────────────────────────────────────
  output$vbox_gfr <- renderValueBox({
    sim <- single_sim()
    val <- if(nrow(sim)>0) round(tail(sim$GFR_c,1)) else "—"
    col <- if(is.numeric(val)&&val>30) "green" else if(is.numeric(val)&&val>15) "yellow" else "red"
    valueBox(paste(val,"mL/min"), "GFR", icon=icon("tint"), color=col)
  })
  output$vbox_proteinuria <- renderValueBox({
    sim <- single_sim()
    val <- if(nrow(sim)>0) round(tail(sim$Proteinuria_c,1),1) else "—"
    col <- if(is.numeric(val)&&val<0.3) "green" else if(is.numeric(val)&&val<1) "yellow" else "red"
    valueBox(paste(val,"g/day"), "Proteinuria", icon=icon("droplet"), color=col)
  })
  output$vbox_gbmdamage <- renderValueBox({
    sim <- single_sim()
    val <- if(nrow(sim)>0) round(tail(sim$GBM_damage,1)) else "—"
    col <- if(is.numeric(val)&&val<20) "green" else "red"
    valueBox(paste(val,"/100"), "GBM Damage", icon=icon("exclamation-triangle"), color=col)
  })
  output$vbox_dialysis <- renderValueBox({
    sim <- single_sim()
    gfr_end <- if(nrow(sim)>0) tail(sim$GFR_c,1) else 30
    lbl <- if(gfr_end<10) "High risk" else if(gfr_end<15) "Moderate" else "Low risk"
    col <- if(gfr_end<10) "red" else if(gfr_end<15) "yellow" else "green"
    valueBox(lbl, "Dialysis Risk", icon=icon("procedures"), color=col)
  })

  output$gfr_plot <- renderPlotly({
    sim <- single_sim()
    if(nrow(sim)==0) return(plotly_empty())
    p <- ggplot(sim, aes(time_wk, GFR_c)) +
      geom_line(color="#377EB8", size=1.2) +
      geom_hline(yintercept=15, linetype="dashed", color="red") +
      geom_hline(yintercept=30, linetype="dotted", color="orange") +
      geom_hline(yintercept=60, linetype="dotted", color="darkgreen") +
      labs(title="GFR Trajectory", x="Time (weeks)", y="GFR (mL/min/1.73m²)") +
      theme_bw(base_size=10)
    ggplotly(p)
  })
  output$gbm_damage_plot <- renderPlotly({
    ab_plot(single_sim(), "GBM_damage", "GBM Damage Score (0–100)", "GBM Damage",
            hline=list(y=20,lty="dashed",col="darkgreen"))
  })
  output$proteinuria_plot <- renderPlotly({
    ab_plot(single_sim(), "Proteinuria_c", "Proteinuria (g/day)", "Proteinuria",
            hline=list(y=0.3,lty="dashed",col="darkgreen"), color="#FF7F00")
  })
  output$renal_summary_tbl <- DT::renderDataTable({
    sim <- single_sim()
    if(nrow(sim)==0) return(data.frame())
    wks <- c(4,8,12,26,52)
    sim %>% filter(time %in% wks*7) %>%
      mutate(Week=time/7, `GFR (mL/min)`=round(GFR_c,1),
             `GBM Damage`=round(GBM_damage,1), `Proteinuria (g/day)`=round(Proteinuria_c,1),
             `Anti-GBM (AU/mL)`=round(AntiGBM,1)) %>%
      select(Week, `GFR (mL/min)`, `GBM Damage`, `Proteinuria (g/day)`, `Anti-GBM (AU/mL)`)
  }, options=list(pageLength=5))

  # ── Tab 5: Pulmonary ─────────────────────────────────────────────────────
  output$vbox_dlco <- renderValueBox({
    sim <- single_sim()
    val <- if(nrow(sim)>0) round(tail(sim$DLCO_c,1)) else "—"
    col <- if(is.numeric(val)&&val>70) "green" else if(is.numeric(val)&&val>40) "yellow" else "red"
    valueBox(paste(val,"%"), "DLCO (% predicted)", icon=icon("lungs"), color=col)
  })
  output$vbox_lungdamage <- renderValueBox({
    sim <- single_sim()
    val <- if(nrow(sim)>0) round(tail(sim$Lung_damage,1)) else "—"
    col <- if(is.numeric(val)&&val<15) "green" else "red"
    valueBox(paste(val,"/100"), "Lung Damage", icon=icon("wind"), color=col)
  })
  output$vbox_resp <- renderValueBox({
    sim <- single_sim()
    dlco_end <- if(nrow(sim)>0) tail(sim$DLCO_c,1) else 0
    lbl <- if(dlco_end>70) "Normal" else if(dlco_end>40) "Impaired" else "Severe"
    col <- if(dlco_end>70) "green" else if(dlco_end>40) "yellow" else "red"
    valueBox(lbl, "Respiratory Status", icon=icon("heartbeat"), color=col)
  })
  output$vbox_pulm_remission <- renderValueBox({
    sim <- single_sim()
    lung_end <- if(nrow(sim)>0) tail(sim$Lung_damage,1) else 100
    lbl <- if(lung_end<15) "Remission" else if(lung_end<30) "Partial" else "Active DAH"
    col <- if(lung_end<15) "green" else if(lung_end<30) "yellow" else "red"
    valueBox(lbl, "Pulmonary Remission", icon=icon("check-circle"), color=col)
  })

  output$dlco_plot <- renderPlotly({
    sim <- single_sim()
    if(nrow(sim)==0) return(plotly_empty())
    p <- ggplot(sim, aes(time_wk, DLCO_c)) +
      geom_line(color="#4DAF4A", size=1.2) +
      geom_hline(yintercept=70, linetype="dashed", color="darkgreen") +
      geom_hline(yintercept=40, linetype="dotted", color="orange") +
      labs(title="DLCO (% Predicted)", x="Time (weeks)", y="DLCO (%)") +
      theme_bw(base_size=10)
    ggplotly(p)
  })
  output$lung_damage_plot <- renderPlotly({
    ab_plot(single_sim(), "Lung_damage", "Lung Damage Score (0–100)", "Alveolar Damage",
            hline=list(y=15,lty="dashed",col="darkgreen"))
  })
  output$pulm_remission_timeline <- renderPlotly({
    sim <- single_sim()
    if(nrow(sim)==0) return(plotly_empty())
    sim$Pulm_remission <- as.integer(sim$Lung_damage < 15)
    p <- ggplot(sim, aes(time_wk, Pulm_remission)) +
      geom_area(fill="#4DAF4A", alpha=0.4) +
      geom_line(color="#4DAF4A", size=1.1) +
      scale_y_continuous(labels=c("Active DAH","Remission"), breaks=c(0,1)) +
      labs(title="Pulmonary Remission Status", x="Time (weeks)", y="Status") +
      theme_bw(base_size=10)
    ggplotly(p)
  })

  # ── Tab 6: Comparison ─────────────────────────────────────────────────────
  comp_line_plot <- function(all_sim, var, ylab, hlines=NULL, title="") {
    if(nrow(all_sim)==0) return(plotly_empty())
    p <- ggplot(all_sim, aes_string("time_wk", var, color="Scenario")) +
      geom_line(size=0.9) +
      scale_color_manual(values=scenario_colors) +
      labs(title=title, x="Time (weeks)", y=ylab, color=NULL) +
      theme_bw(base_size=9) + theme(legend.position="bottom") +
      guides(color=guide_legend(nrow=2))
    if(!is.null(hlines)) {
      for(h in hlines) p <- p + geom_hline(yintercept=h$y, linetype=h$lty, color=h$col)
    }
    ggplotly(p)
  }

  output$comp_antiGBM   <- renderPlotly({ comp_line_plot(all_sims(),"AntiGBM","Anti-GBM (AU/mL)",
    list(list(y=10,lty="dashed",col="green")), "Anti-GBM Titer") })
  output$comp_gfr       <- renderPlotly({ comp_line_plot(all_sims(),"GFR_c","GFR (mL/min)",
    list(list(y=15,lty="dashed",col="red")), "GFR Trajectory") })
  output$comp_dlco      <- renderPlotly({ comp_line_plot(all_sims(),"DLCO_c","DLCO (%)",
    list(list(y=70,lty="dashed",col="green")), "DLCO") })
  output$comp_gbmdamage <- renderPlotly({ comp_line_plot(all_sims(),"GBM_damage","GBM Score",NULL,"GBM Damage") })
  output$comp_bcells    <- renderPlotly({ comp_line_plot(all_sims(),"B_cells","B cells (%)",NULL,"B Cell Depletion") })

  output$comp_outcomes_tbl <- DT::renderDataTable({
    sim <- all_sims()
    if(nrow(sim)==0) return(data.frame())
    wks <- c(4,12,26,52)
    sim %>% filter(time %in% wks*7) %>%
      mutate(Week=paste0("Wk",time/7),
             `Anti-GBM`=round(AntiGBM,1), GFR=round(GFR_c,1),
             DLCO=round(DLCO_c,1), `GBM Damage`=round(GBM_damage,1),
             `Dialysis risk`=ifelse(GFR_c<10,"High",ifelse(GFR_c<15,"Moderate","Low"))) %>%
      select(Week, Scenario, `Anti-GBM`, GFR, DLCO, `GBM Damage`, `Dialysis risk`)
  }, options=list(pageLength=24))

  # ── Tab 7: References ─────────────────────────────────────────────────────
  output$trials_tbl <- renderTable({
    data.frame(
      Trial_Cohort = c("Levy et al. 2001","Alchi et al. 2015","Lazor et al. 2007","Segelmark et al. 2003","Griffiths et al. 1989"),
      n            = c(71, 88, 28, 53, 30),
      `1yr Renal Survival (%)` = c("60–80","55–75","N/A","72","67"),
      PMID         = c("11388816","25342255","17519713","12816369","2780498")
    )
  })
  output$pkpd_refs_tbl <- renderTable({
    data.frame(
      Drug   = c("Cyclophosphamide","Prednisolone","Rituximab","Avacopan"),
      PK_Ref = c("Uppugunduri 2013","Ternant 2019","Wang et al. 2008","Nester et al. 2020"),
      PMID   = c("23568534","30919462","18196377","32279213")
    )
  })
  output$biology_refs_tbl <- renderTable({
    data.frame(
      Topic  = c("Goodpasture antigen (COL4A3 NC1)","HLA association","Pathogenesis review",
                  "Complement in anti-GBM","Plasmapheresis mechanism","Linear IgG deposition"),
      Author = c("Hudson BG","Rees AJ","McAdoo & Pusey","Tang Z","Turner N","Lockwood CM"),
      Year   = c(2003,1978,2017,2021,1996,1976),
      PMID   = c("12815141","351459","28515156","33432762","10227016","56668")
    )
  })
}

# ── Launch ────────────────────────────────────────────────────────────────────
shinyApp(ui, server)
