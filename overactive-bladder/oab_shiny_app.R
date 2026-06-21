## ============================================================
##  Overactive Bladder (OAB) — Interactive Shiny Dashboard
##  6-Tab layout: Patient Profile · PK · Bladder PD ·
##                Clinical Endpoints · Scenario Comparison · Biomarkers
##  Requires: shiny, shinydashboard, plotly, dplyr, tidyr, deSolve
## ============================================================

library(shiny)
library(shinydashboard)
library(plotly)
library(dplyr)
library(tidyr)
library(DT)

# ────────────────────────────────────────────────────────────────
#  Pure R ODE solver (no mrgsolve dependency for Shiny deployment)
# ────────────────────────────────────────────────────────────────

oab_ode <- function(t, state, parms) {
  with(as.list(c(state, parms)), {
    # Plasma concentrations (ng/mL) via 1-CMT
    Cp_OXY  <- OXY_CENT  / Vc_OXY
    Cp_TOL  <- TOL_CENT  / Vc_TOL
    Cp_SOL  <- SOL_CENT  / Vc_SOL
    Cp_MIR  <- MIR_CENT  / Vc_MIR
    Cp_SOL2 <- SOL2_CENT / Vc_SOL2

    # Drug input rates (pulsed via events outside ODE; here = 0 between doses)
    dOXY_GUT   <- -KA_OXY  * OXY_GUT
    dOXY_CENT  <-  KA_OXY  * OXY_GUT * F_OXY  - (CL_OXY  / Vc_OXY)  * OXY_CENT
    dTOL_GUT   <- -KA_TOL  * TOL_GUT
    dTOL_CENT  <-  KA_TOL  * TOL_GUT * F_TOL  - (CL_TOL  / Vc_TOL)  * TOL_CENT
    dSOL_GUT   <- -KA_SOL  * SOL_GUT
    dSOL_CENT  <-  KA_SOL  * SOL_GUT * F_SOL  - (CL_SOL  / Vc_SOL)  * SOL_CENT
    dMIR_GUT   <- -KA_MIR  * MIR_GUT
    dMIR_CENT  <-  KA_MIR  * MIR_GUT * F_MIR  - (CL_MIR  / Vc_MIR)  * MIR_CENT
    dSOL2_GUT  <- -KA_SOL2 * SOL2_GUT
    dSOL2_CENT <-  KA_SOL2 * SOL2_GUT * F_SOL2 - (CL_SOL2 / Vc_SOL2) * SOL2_CENT

    # Receptor occupancies
    RO_OXY_eq  <- SCN_OXY  * Cp_OXY  / (Cp_OXY  + EC50_OXY_M3)
    RO_TOL_eq  <- SCN_TOL  * Cp_TOL  / (Cp_TOL  + EC50_TOL_M3)
    RO_SOL_eq  <- SCN_SOL  * Cp_SOL  / (Cp_SOL  + EC50_SOL_M3)
    RO_SOL2_eq <- SCN_COMB * Cp_SOL2 / (Cp_SOL2 + EC50_SOL2_M3)
    RO_M3_eq   <- 1 - (1-RO_OXY_eq)*(1-RO_TOL_eq)*(1-RO_SOL_eq)*(1-RO_SOL2_eq)
    if (RO_M3_eq > 0.99) RO_M3_eq <- 0.99
    RO_MIR_eq  <- (SCN_MIR + SCN_COMB * 0.5) * Cp_MIR / (Cp_MIR + EC50_MIR_B3)
    if (RO_MIR_eq > 0.99) RO_MIR_eq <- 0.99

    dRO_M3 <- 1.0 * (RO_M3_eq - RO_M3)
    dRO_B3 <- 0.8 * (RO_MIR_eq - RO_B3)

    # PD states
    inh_DO <- pmin(Emax_DO_M3 * RO_M3/(1+RO_M3) + Emax_DO_B3 * RO_B3/(1+RO_B3), 0.90)
    dDetAct  <- KOUT_DO  * BASE_DO  * (1 - inh_DO) - KOUT_DO  * DetAct

    stim_CAP <- pmin(Emax_CAP_M3 * RO_M3 + Emax_CAP_B3 * RO_B3, 0.50)
    dBladCap <- KOUT_CAP * BASE_CAP * (1 + stim_CAP) - KOUT_CAP * BladCap

    dVoidFreq <- KOUT_VOID * BASE_VOID * DetAct - KOUT_VOID * VoidFreq
    dUrgency  <- KOUT_URG  * BASE_URG  * DetAct - KOUT_URG  * Urgency
    dUUI      <- KOUT_UUI  * BASE_UUI  * (Urgency/BASE_URG) * DetAct - KOUT_UUI * UUI

    inh_NGF   <- Emax_NGF_M3 * RO_M3
    dNGF      <- KOUT_NGF * BASE_NGF * DetAct * (1 - inh_NGF) - KOUT_NGF * NGF
    dATP_bm   <- KOUT_ATP * BASE_ATP * DetAct - KOUT_ATP * ATP_bm

    dContScore <- (KOUT_DO * BASE_CONT / DetAct) - KOUT_DO * ContScore
    dNocturia  <- 0.025 * 2.5 * DetAct - 0.025 * Nocturia

    OABq_drv  <- (Urgency/BASE_URG + UUI/BASE_UUI + Nocturia/2.5) / 3
    dOABq     <- 0.04 * BASE_DO * 62 * OABq_drv - 0.04 * OABq

    list(c(dOXY_GUT, dOXY_CENT, dTOL_GUT, dTOL_CENT, dSOL_GUT, dSOL_CENT,
           dMIR_GUT, dMIR_CENT, dSOL2_GUT, dSOL2_CENT,
           dRO_M3, dRO_B3, dDetAct, dBladCap, dVoidFreq, dUrgency,
           dUUI, dNGF, dATP_bm, dContScore, dNocturia, dOABq))
  })
}

run_sim <- function(scenario, weight_kg, age_yr, severity, tend_days) {
  tend_h <- tend_days * 24
  dose_interval <- switch(scenario,
    "Oxybutynin IR 5mg TID"    = 8,
    "Tolterodine ER 4mg QD"    = 24,
    "Solifenacin 10mg QD"      = 24,
    "Mirabegron 50mg QD"       = 24,
    "Sofeli 5mg + Mirabe 25mg" = 24,
    24
  )

  # Scale parameters with patient covariates
  base_scale <- severity  # 1.0 = moderate OAB, 1.5 = severe

  parms <- c(
    KA_OXY=1.30, F_OXY=0.06, CL_OXY=85*(weight_kg/70)^0.75, Vc_OXY=193,
    KA_TOL=0.70, F_TOL=0.65, CL_TOL=46.2*(weight_kg/70)^0.75, Vc_TOL=310,
    KA_SOL=0.40, F_SOL=0.90, CL_SOL=9.8*(weight_kg/70)^0.75,  Vc_SOL=600,
    KA_MIR=0.50, F_MIR=0.32, CL_MIR=57.2*(weight_kg/70)^0.75, Vc_MIR=1670,
    KA_SOL2=0.40, F_SOL2=0.90, CL_SOL2=9.8*(weight_kg/70)^0.75, Vc_SOL2=600,
    EC50_OXY_M3=1.2, EC50_TOL_M3=4.8, EC50_SOL_M3=5.5, EC50_SOL2_M3=5.5,
    EC50_MIR_B3=22.0,
    BASE_DO=base_scale, BASE_CAP=250/base_scale, BASE_VOID=12*base_scale,
    BASE_URG=5.8*base_scale, BASE_UUI=3.5*base_scale, BASE_NGF=3.2*base_scale,
    BASE_ATP=2.5*base_scale, BASE_CONT=30/base_scale,
    KOUT_DO=0.05, KOUT_CAP=0.02, KOUT_VOID=0.03, KOUT_URG=0.04,
    KOUT_UUI=0.03, KOUT_NGF=0.008, KOUT_ATP=0.015,
    Emax_DO_M3=0.60, Emax_DO_B3=0.45, Emax_CAP_B3=0.35, Emax_CAP_M3=0.20,
    Emax_NGF_M3=0.30,
    SCN_OXY  = as.integer(scenario == "Oxybutynin IR 5mg TID"),
    SCN_TOL  = as.integer(scenario == "Tolterodine ER 4mg QD"),
    SCN_SOL  = as.integer(scenario == "Solifenacin 10mg QD"),
    SCN_MIR  = as.integer(scenario == "Mirabegron 50mg QD"),
    SCN_COMB = as.integer(scenario == "Sofeli 5mg + Mirabe 25mg")
  )

  state0 <- c(
    OXY_GUT=0, OXY_CENT=0, TOL_GUT=0, TOL_CENT=0,
    SOL_GUT=0, SOL_CENT=0, MIR_GUT=0, MIR_CENT=0,
    SOL2_GUT=0, SOL2_CENT=0,
    RO_M3=0, RO_B3=0,
    DetAct=parms["BASE_DO"], BladCap=parms["BASE_CAP"],
    VoidFreq=parms["BASE_VOID"], Urgency=parms["BASE_URG"],
    UUI=parms["BASE_UUI"], NGF=parms["BASE_NGF"], ATP_bm=parms["BASE_ATP"],
    ContScore=parms["BASE_CONT"], Nocturia=2.5*base_scale, OABq=62*base_scale
  )

  # Build dosing events
  dose_mg <- switch(scenario,
    "Oxybutynin IR 5mg TID"    = 5000,   # µg
    "Tolterodine ER 4mg QD"    = 4000,
    "Solifenacin 10mg QD"      = 10000,
    "Mirabegron 50mg QD"       = 50000,
    "Sofeli 5mg + Mirabe 25mg" = 5000,
    0
  )
  dose_cmt_sol <- c("SOL2_GUT", "MIR_GUT")

  # Use ode() with event data
  dose_times <- seq(0, tend_h, by = dose_interval)
  events_df  <- data.frame(
    var   = if (scenario == "Sofeli 5mg + Mirabe 25mg") {
              rep(c("SOL2_GUT","MIR_GUT"), length(dose_times))
            } else {
              rep(switch(scenario,
                "Oxybutynin IR 5mg TID" = "OXY_GUT",
                "Tolterodine ER 4mg QD" = "TOL_GUT",
                "Solifenacin 10mg QD"   = "SOL_GUT",
                "Mirabegron 50mg QD"    = "MIR_GUT",
                "SOL2_GUT"), length(dose_times))
            },
    time  = if (scenario == "Sofeli 5mg + Mirabe 25mg") {
              rep(dose_times, each = 2)
            } else { dose_times },
    value = if (scenario == "Sofeli 5mg + Mirabe 25mg") {
              rep(c(5000, 25000), length(dose_times))
            } else { rep(dose_mg, length(dose_times)) },
    method = "add"
  )

  times <- sort(unique(c(seq(0, tend_h, 4), events_df$time)))

  out <- deSolve::ode(
    y      = state0,
    times  = times,
    func   = oab_ode,
    parms  = parms,
    events = list(data = events_df),
    method = "lsoda"
  )
  df <- as.data.frame(out)
  df$time_day <- df$time / 24
  df$Cp_OXY  <- df$OXY_CENT  / parms["Vc_OXY"]
  df$Cp_TOL  <- df$TOL_CENT  / parms["Vc_TOL"]
  df$Cp_SOL  <- df$SOL_CENT  / parms["Vc_SOL"]
  df$Cp_MIR  <- df$MIR_CENT  / parms["Vc_MIR"]
  df$scenario <- scenario
  df
}

# ────────────────────────────────────────────────────────────────
#  UI
# ────────────────────────────────────────────────────────────────
ui <- dashboardPage(
  skin = "blue",
  dashboardHeader(title = "OAB QSP Dashboard"),

  dashboardSidebar(
    sidebarMenu(
      menuItem("① Patient Profile",    tabName = "profile",    icon = icon("user")),
      menuItem("② Drug PK",            tabName = "pk",         icon = icon("flask")),
      menuItem("③ Bladder PD",         tabName = "pd",         icon = icon("dna")),
      menuItem("④ Clinical Endpoints", tabName = "endpoints",  icon = icon("chart-line")),
      menuItem("⑤ Scenario Comparison",tabName = "comparison", icon = icon("balance-scale")),
      menuItem("⑥ Biomarkers",         tabName = "biomarkers", icon = icon("vial"))
    ),
    hr(),
    h5("Patient Parameters", style = "padding-left:10px; color:#aaa"),
    sliderInput("weight",    "Weight (kg)",    min=45, max=120, value=68, step=1),
    sliderInput("age",       "Age (years)",    min=18, max=85,  value=58, step=1),
    sliderInput("severity",  "OAB Severity",   min=0.7, max=2.0, value=1.0, step=0.1),
    hr(),
    h5("Treatment", style = "padding-left:10px; color:#aaa"),
    selectInput("scenario", "Drug / Scenario",
      choices = c("No Treatment",
                  "Oxybutynin IR 5mg TID",
                  "Tolterodine ER 4mg QD",
                  "Solifenacin 10mg QD",
                  "Mirabegron 50mg QD",
                  "Sofeli 5mg + Mirabe 25mg"),
      selected = "No Treatment"),
    sliderInput("tend", "Simulation Duration (days)", min=14, max=168, value=84, step=7),
    actionButton("run_btn", "Run Simulation", icon = icon("play"), class = "btn-success",
                 style = "margin:10px")
  ),

  dashboardBody(
    tags$head(tags$style(HTML(".box-header .box-title { font-size:14px; }"))),
    tabItems(

      # ── Tab 1: Patient Profile ──────────────────────────────────
      tabItem(tabName = "profile",
        fluidRow(
          valueBoxOutput("vbox_voids"),
          valueBoxOutput("vbox_uui"),
          valueBoxOutput("vbox_urgency"),
          valueBoxOutput("vbox_oabq")
        ),
        fluidRow(
          box(title = "OAB Diagnosis Criteria", width = 6, status = "primary",
            p("Overactive bladder (OAB) is defined by the ICS as:"),
            tags$ul(
              tags$li(strong("Urgency"), " — compelling desire to void that is difficult to defer"),
              tags$li(strong("Urgency Urinary Incontinence (UUI)"), " — involuntary leakage with urgency"),
              tags$li(strong("Frequency"), " — ≥8 voids/24h"),
              tags$li(strong("Nocturia"), " — ≥2 nocturnal voids")
            ),
            hr(),
            p(strong("OAB Types:")),
            tags$ul(
              tags$li(strong("OAB-wet"), " — urgency + UUI"),
              tags$li(strong("OAB-dry"), " — urgency without UUI"),
              tags$li(strong("Neurogenic OAB"), " — SCI, MS, PD, stroke-related")
            )
          ),
          box(title = "Drug Mechanism Reference", width = 6, status = "info",
            DTOutput("drug_mech_tbl")
          )
        ),
        fluidRow(
          box(title = "Key Clinical Trials", width = 12, status = "success",
            DTOutput("trial_tbl")
          )
        )
      ),

      # ── Tab 2: Drug PK ──────────────────────────────────────────
      tabItem(tabName = "pk",
        fluidRow(
          box(title = "Plasma Concentration–Time Profile", width = 12, status = "primary",
            plotlyOutput("pk_plot", height = "400px")
          )
        ),
        fluidRow(
          box(title = "Receptor Occupancy (Early Phase, 0–14 days)", width = 7, status = "info",
            plotlyOutput("ro_plot", height = "350px")
          ),
          box(title = "PK Parameter Summary", width = 5, status = "warning",
            DTOutput("pk_param_tbl")
          )
        )
      ),

      # ── Tab 3: Bladder PD ────────────────────────────────────────
      tabItem(tabName = "pd",
        fluidRow(
          box(title = "Detrusor Overactivity Index", width = 6, status = "danger",
            plotlyOutput("do_plot", height = "320px")
          ),
          box(title = "Bladder Cystometric Capacity (mL)", width = 6, status = "primary",
            plotlyOutput("cap_plot", height = "320px")
          )
        ),
        fluidRow(
          box(title = "M3/β3 Receptor Occupancy — Steady State", width = 6, status = "info",
            plotlyOutput("ro_steady_plot", height = "300px")
          ),
          box(title = "Smooth Muscle Contraction Cascade", width = 6,
            tags$table(class="table table-bordered table-condensed",
              tags$thead(tags$tr(tags$th("Step"), tags$th("Signal"), tags$th("Drug Effect"))),
              tags$tbody(
                tags$tr(tags$td("1"), tags$td("ACh → M3 → Gq"), tags$td("Antimuscarinics block M3")),
                tags$tr(tags$td("2"), tags$td("PLC → IP3 + DAG"), tags$td("↓IP3 → ↓Ca²⁺")),
                tags$tr(tags$td("3"), tags$td("[Ca²⁺]i → CaM"), tags$td("↓MLCK activation")),
                tags$tr(tags$td("4"), tags$td("MLCK → MLC-P"), tags$td("↓Contraction")),
                tags$tr(tags$td("5"), tags$td("β3-AR → cAMP → PKA"), tags$td("↑MLCP → Relaxation")),
                tags$tr(tags$td("6"), tags$td("RhoA/ROCK → MLCP↓"), tags$td("ROCK inhibitors target")  )
              )
            )
          )
        )
      ),

      # ── Tab 4: Clinical Endpoints ────────────────────────────────
      tabItem(tabName = "endpoints",
        fluidRow(
          box(title = "UUI Episodes / 24h", width = 6, status = "danger",
            plotlyOutput("uui_plot", height = "300px")
          ),
          box(title = "Urgency Episodes / 24h", width = 6, status = "warning",
            plotlyOutput("urg_plot", height = "300px")
          )
        ),
        fluidRow(
          box(title = "Voiding Frequency / 24h", width = 6, status = "info",
            plotlyOutput("void_plot", height = "300px")
          ),
          box(title = "OAB-q Symptom Bother Score", width = 6, status = "primary",
            plotlyOutput("oabq_plot", height = "300px")
          )
        ),
        fluidRow(
          box(title = "Weekly Summary Table", width = 12,
            DTOutput("endpoint_tbl")
          )
        )
      ),

      # ── Tab 5: Scenario Comparison ───────────────────────────────
      tabItem(tabName = "comparison",
        fluidRow(
          box(title = "UUI Reduction — All 6 Scenarios", width = 12, status = "danger",
            plotlyOutput("comp_uui_plot", height = "380px")
          )
        ),
        fluidRow(
          box(title = "Bladder Capacity — All 6 Scenarios", width = 6, status = "primary",
            plotlyOutput("comp_cap_plot", height = "320px")
          ),
          box(title = "OAB-q Score — All 6 Scenarios", width = 6, status = "warning",
            plotlyOutput("comp_oabq_plot", height = "320px")
          )
        ),
        fluidRow(
          box(title = "12-Week Summary: All Scenarios", width = 12,
            DTOutput("comp_tbl")
          )
        )
      ),

      # ── Tab 6: Biomarkers ────────────────────────────────────────
      tabItem(tabName = "biomarkers",
        fluidRow(
          box(title = "Urinary NGF (Normalized Ratio)", width = 6, status = "primary",
            plotlyOutput("ngf_plot", height = "300px")
          ),
          box(title = "Urinary ATP (Purinergic Biomarker)", width = 6, status = "info",
            plotlyOutput("atp_plot", height = "300px")
          )
        ),
        fluidRow(
          box(title = "Dose–Response: Solifenacin vs Urgency @ Week 12",
              width = 6, status = "warning",
            plotlyOutput("dr_sol_plot", height = "300px")
          ),
          box(title = "Biomarker Reference Values", width = 6,
            DTOutput("bm_ref_tbl")
          )
        ),
        fluidRow(
          box(title = "Nocturia Episodes", width = 6, status = "danger",
            plotlyOutput("noct_plot", height = "280px")
          ),
          box(title = "Continence Score (0–100)", width = 6, status = "success",
            plotlyOutput("cont_plot", height = "280px")
          )
        )
      )
    )
  )
)

# ────────────────────────────────────────────────────────────────
#  Server
# ────────────────────────────────────────────────────────────────
server <- function(input, output, session) {

  # ── Reactive simulation ──
  sim_data <- eventReactive(input$run_btn, {
    run_sim(input$scenario, input$weight, input$age, input$severity, input$tend)
  }, ignoreNULL = FALSE)

  # All 6 scenarios for comparison (run once at startup at default params)
  all_scn_data <- eventReactive(input$run_btn, {
    scn_list <- c("No Treatment", "Oxybutynin IR 5mg TID",
                  "Tolterodine ER 4mg QD", "Solifenacin 10mg QD",
                  "Mirabegron 50mg QD", "Sofeli 5mg + Mirabe 25mg")
    withProgress(message = "Running all 6 scenarios...", {
      res <- lapply(scn_list, function(s) {
        run_sim(s, input$weight, input$age, input$severity, input$tend)
      })
    })
    bind_rows(res)
  }, ignoreNULL = FALSE)

  pal_6 <- c("#7F7F7F","#E74C3C","#E67E22","#2980B9","#27AE60","#8E44AD")
  scn_names <- c("No Treatment", "Oxybutynin IR 5mg TID", "Tolterodine ER 4mg QD",
                 "Solifenacin 10mg QD", "Mirabegron 50mg QD", "Sofeli 5mg + Mirabe 25mg")

  # ── Value boxes ──
  output$vbox_voids    <- renderValueBox(
    valueBox(round(tail(sim_data()$VoidFreq, 1), 1), "Voids/24h", icon("toilet"), color="yellow"))
  output$vbox_uui      <- renderValueBox(
    valueBox(round(tail(sim_data()$UUI, 1), 1), "UUI/24h",   icon("tint"),    color="red"))
  output$vbox_urgency  <- renderValueBox(
    valueBox(round(tail(sim_data()$Urgency, 1), 1), "Urgency/24h", icon("bolt"), color="orange"))
  output$vbox_oabq     <- renderValueBox(
    valueBox(round(tail(sim_data()$OABq, 1), 1), "OAB-q Score", icon("star"), color="blue"))

  # ── Drug mechanism table ──
  output$drug_mech_tbl <- renderDT({
    datatable(data.frame(
      Drug = c("Oxybutynin IR","Tolterodine ER","Solifenacin","Darifenacin","Fesoterodine",
               "Mirabegron","Vibegron","OnaBotulinumtoxinA"),
      Class = c(rep("Antimuscarinic",5), rep("β3-Agonist",2), "Neurotoxin"),
      Selectivity = c("M1-M5 (non-sel)","M1-M5","M3>M1","M3>>>M1","M3>M2","β3>>β1,β2","β3>>β1,β2","SNAP25"),
      `t½` = c("2h (NDEO 8h)","12h (ER)","45-68h","13-19h","7h (5-HM)","50h","30.8h","Months"),
      `F (%)` = c(6,65,90,15,NA,32,60,NA)
    ), options = list(pageLength = 8, dom = 't'), rownames = FALSE)
  })

  # ── Clinical trials table ──
  output$trial_tbl <- renderDT({
    datatable(data.frame(
      Trial = c("OBJECT","ACET","STAR","SCORPIO","BESIDE","EMBARK","VIBRATO","Chapple 2013"),
      Drug = c("Oxy vs Tol","Tolterodine ER","Solifenacin","Mirabegron","Sol+Mir","OnaBoTox","Vibegron","Mirabegron"),
      N = c(378,1235,1033,1978,2125,548,1518,1306),
      Duration = c("12wk","12wk","12wk","12wk","12wk","24wk","12wk","12wk"),
      `UUI Change` = c("-4.2 vs -3.5","-1.4","−71% from BL","−47.8%","−3.23 vs −2.86","−2.65","-1.84","−47.8%"),
      Journal = c("BJC 2001","EurUrol 2005","BJCP 2005","EurUrol 2013","EurUrol 2016","NEJM 2020","JAMA IM 2019","EurUrol 2013")
    ), options = list(pageLength = 10, dom = 'ft'), rownames = FALSE)
  })

  # ── PK plots ──
  output$pk_plot <- renderPlotly({
    df <- sim_data()
    cp_col <- switch(input$scenario,
      "Oxybutynin IR 5mg TID"    = "Cp_OXY",
      "Tolterodine ER 4mg QD"    = "Cp_TOL",
      "Solifenacin 10mg QD"      = "Cp_SOL",
      "Mirabegron 50mg QD"       = "Cp_MIR",
      "Sofeli 5mg + Mirabe 25mg" = "Cp_SOL",
      NULL)
    if (is.null(cp_col) || input$scenario == "No Treatment") {
      return(plot_ly() %>% layout(title = "No drug selected"))
    }
    df_pk <- df %>% filter(time_day <= min(14, input$tend))
    plot_ly(df_pk, x = ~time_day, y = ~get(cp_col), type = "scatter", mode = "lines",
            name = cp_col, line = list(width = 2)) %>%
      layout(title = paste("Plasma Cp —", input$scenario),
             xaxis = list(title = "Day"), yaxis = list(title = "Cp (ng/mL)"))
  })

  output$ro_plot <- renderPlotly({
    df <- sim_data() %>% filter(time_day <= 14)
    plot_ly() %>%
      add_lines(data = df, x = ~time_day, y = ~RO_M3*100, name = "M3 RO (%)",
                line = list(color = "#2980B9", width = 2)) %>%
      add_lines(data = df, x = ~time_day, y = ~RO_B3*100, name = "β3 RO (%)",
                line = list(color = "#27AE60", width = 2, dash = "dash")) %>%
      layout(title = "Receptor Occupancy (0–14 days)",
             xaxis = list(title = "Day"), yaxis = list(title = "RO (%)"))
  })

  output$pk_param_tbl <- renderDT({
    datatable(data.frame(
      Drug       = c("Oxybutynin IR","Tolterodine ER","Solifenacin","Mirabegron"),
      `ka (/h)`  = c(1.30, 0.70, 0.40, 0.50),
      `F (%)`    = c(6, 65, 90, 32),
      `CL (L/h)` = c(85.0, 46.2, 9.8, 57.2),
      `Vc (L)`   = c(193, 310, 600, 1670),
      `t½ (h)`   = c(2, 12, 56, 50),
      `EC50 (ng/mL)` = c(1.2, 4.8, 5.5, "22 (β3)")
    ), options = list(dom = 't'), rownames = FALSE)
  })

  # ── PD plots ──
  output$do_plot <- renderPlotly({
    df <- sim_data()
    plot_ly(df, x = ~time_day, y = ~DetAct, type = "scatter", mode = "lines",
            name = "DO Index", line = list(color = "#E74C3C", width = 2)) %>%
      layout(title = "Detrusor Overactivity Index",
             xaxis = list(title = "Day"), yaxis = list(title = "DO (1 = severe OAB baseline)"))
  })

  output$cap_plot <- renderPlotly({
    df <- sim_data()
    plot_ly(df, x = ~time_day, y = ~BladCap, type = "scatter", mode = "lines",
            name = "Bladder Capacity", line = list(color = "#2980B9", width = 2)) %>%
      layout(title = "Cystometric Bladder Capacity",
             xaxis = list(title = "Day"), yaxis = list(title = "Capacity (mL)"))
  })

  output$ro_steady_plot <- renderPlotly({
    df <- sim_data() %>% filter(time_day >= input$tend * 0.8)
    ss_m3 <- mean(df$RO_M3) * 100
    ss_b3 <- mean(df$RO_B3) * 100
    plot_ly(x = c("M3 Occupancy", "β3 Occupancy"),
            y = c(ss_m3, ss_b3),
            type = "bar",
            marker = list(color = c("#2980B9", "#27AE60"))) %>%
      layout(title = "Steady-State Receptor Occupancy",
             xaxis = list(title = ""), yaxis = list(title = "RO (%)", range = c(0,100)))
  })

  # ── Clinical endpoint plots ──
  output$uui_plot <- renderPlotly({
    df <- sim_data()
    plot_ly(df, x = ~time_day, y = ~UUI, type = "scatter", mode = "lines",
            name = "UUI/24h", line = list(color = "#C0392B", width = 2)) %>%
      layout(xaxis = list(title = "Day"), yaxis = list(title = "Episodes/24h"))
  })

  output$urg_plot <- renderPlotly({
    df <- sim_data()
    plot_ly(df, x = ~time_day, y = ~Urgency, type = "scatter", mode = "lines",
            name = "Urgency/24h", line = list(color = "#E67E22", width = 2)) %>%
      layout(xaxis = list(title = "Day"), yaxis = list(title = "Episodes/24h"))
  })

  output$void_plot <- renderPlotly({
    df <- sim_data()
    plot_ly(df, x = ~time_day, y = ~VoidFreq, type = "scatter", mode = "lines",
            name = "Voids/24h", line = list(color = "#8E44AD", width = 2)) %>%
      layout(xaxis = list(title = "Day"), yaxis = list(title = "Voids/24h"))
  })

  output$oabq_plot <- renderPlotly({
    df <- sim_data()
    plot_ly(df, x = ~time_day, y = ~OABq, type = "scatter", mode = "lines",
            name = "OAB-q", line = list(color = "#2980B9", width = 2)) %>%
      layout(xaxis = list(title = "Day"), yaxis = list(title = "OAB-q Score"))
  })

  output$endpoint_tbl <- renderDT({
    df <- sim_data()
    wk_days <- c(4, 8, 12, 24) * 7
    tbl <- df %>%
      filter(round(time_day) %in% wk_days) %>%
      group_by(time_day) %>%
      summarise(
        Week     = round(time_day/7),
        `Voids/24h`  = round(mean(VoidFreq), 2),
        `UUI/24h`    = round(mean(UUI), 2),
        `Urgency/24h`= round(mean(Urgency), 2),
        `BladCap (mL)` = round(mean(BladCap), 1),
        `OABq Score`   = round(mean(OABq), 1),
        `M3 RO (%)`    = round(mean(RO_M3)*100, 1),
        `β3 RO (%)`    = round(mean(RO_B3)*100, 1),
        .groups = "drop"
      ) %>% select(-time_day)
    datatable(tbl, options = list(dom = 't', pageLength = 20), rownames = FALSE)
  })

  # ── Comparison plots ──
  output$comp_uui_plot <- renderPlotly({
    df <- all_scn_data()
    p  <- plot_ly()
    for (i in seq_along(scn_names)) {
      d <- df %>% filter(scenario == scn_names[i])
      p <- add_lines(p, x = d$time_day, y = d$UUI, name = scn_names[i],
                     line = list(color = pal_6[i], width = 2))
    }
    p %>% layout(title = "UUI Episodes / 24h — 6 Scenarios",
                 xaxis = list(title = "Day"), yaxis = list(title = "UUI/24h"))
  })

  output$comp_cap_plot <- renderPlotly({
    df <- all_scn_data()
    p  <- plot_ly()
    for (i in seq_along(scn_names)) {
      d <- df %>% filter(scenario == scn_names[i])
      p <- add_lines(p, x = d$time_day, y = d$BladCap, name = scn_names[i],
                     line = list(color = pal_6[i], width = 2))
    }
    p %>% layout(title = "Bladder Capacity (mL)",
                 xaxis = list(title = "Day"), yaxis = list(title = "mL"))
  })

  output$comp_oabq_plot <- renderPlotly({
    df <- all_scn_data()
    p  <- plot_ly()
    for (i in seq_along(scn_names)) {
      d <- df %>% filter(scenario == scn_names[i])
      p <- add_lines(p, x = d$time_day, y = d$OABq, name = scn_names[i],
                     line = list(color = pal_6[i], width = 2))
    }
    p %>% layout(title = "OAB-q Bother Score",
                 xaxis = list(title = "Day"), yaxis = list(title = "Score (higher=worse)"))
  })

  output$comp_tbl <- renderDT({
    df <- all_scn_data()
    tbl <- df %>%
      filter(round(time_day) == 84) %>%
      group_by(scenario) %>%
      summarise(
        `UUI/24h`     = round(mean(UUI), 2),
        `Urgency/24h` = round(mean(Urgency), 2),
        `Voids/24h`   = round(mean(VoidFreq), 2),
        `BladCap (mL)`= round(mean(BladCap), 1),
        `OABq Score`  = round(mean(OABq), 1),
        `M3 RO (%)`   = round(mean(RO_M3)*100, 1),
        `β3 RO (%)`   = round(mean(RO_B3)*100, 1),
        .groups = "drop"
      )
    datatable(tbl, options = list(dom = 't'), rownames = FALSE) %>%
      formatStyle('UUI/24h', background = styleColorBar(range(tbl$`UUI/24h`), 'lightblue'))
  })

  # ── Biomarker plots ──
  output$ngf_plot <- renderPlotly({
    df <- sim_data()
    plot_ly(df, x = ~time_day, y = ~NGF, type = "scatter", mode = "lines",
            name = "Urinary NGF", line = list(color = "#8E44AD", width = 2)) %>%
      add_lines(x = df$time_day, y = rep(1.0, nrow(df)), name = "Normal",
                line = list(color = "gray", dash = "dash")) %>%
      layout(title = "Urinary NGF (Normalized Ratio)",
             xaxis = list(title = "Day"), yaxis = list(title = "NGF ratio"))
  })

  output$atp_plot <- renderPlotly({
    df <- sim_data()
    plot_ly(df, x = ~time_day, y = ~ATP_bm, type = "scatter", mode = "lines",
            name = "Urinary ATP", line = list(color = "#27AE60", width = 2)) %>%
      layout(title = "Urinary ATP (Purinergic)",
             xaxis = list(title = "Day"), yaxis = list(title = "ATP ratio"))
  })

  output$dr_sol_plot <- renderPlotly({
    doses <- c(2.5, 5, 7.5, 10, 15, 20)
    urg_at_wk12 <- sapply(doses, function(d) {
      df_dr <- run_sim("Solifenacin 10mg QD", input$weight, input$age, input$severity, 84)
      tail(df_dr$Urgency, 1)
    })
    plot_ly(x = doses, y = urg_at_wk12, type = "scatter", mode = "lines+markers",
            name = "Urgency @ Wk12", marker = list(size = 8),
            line = list(color = "#2980B9")) %>%
      layout(title = "Solifenacin Dose–Response (Week 12)",
             xaxis = list(title = "Dose (mg QD)"),
             yaxis = list(title = "Urgency Episodes/24h"))
  })

  output$bm_ref_tbl <- renderDT({
    datatable(data.frame(
      Biomarker  = c("Urinary NGF", "Urinary BDNF", "Urinary ATP",
                     "Urinary PGE2", "Bladder Wall Thickness", "Uroflow Qmax"),
      Normal     = c("0.8–1.2 (norm)", "50–150 pg/mg Cr", "Low",
                     "1–3 ng/mg Cr", "3–5 mm", ">15 mL/s"),
      `OAB High` = c(">3x normal", ">300 pg/mg Cr", "3–5x", ">10 ng/mg Cr",
                     "5–10 mm", "<10 mL/s"),
      Reference  = c("Yoshida 2006 JUrol", "Liu 2010 JUrol", "Sun 2009 BJU",
                     "Kim 2011 Urology", "Oelke 2013 EurUrol", "Drake 2014")
    ), options = list(dom = 't', pageLength = 10), rownames = FALSE)
  })

  output$noct_plot <- renderPlotly({
    df <- sim_data()
    plot_ly(df, x = ~time_day, y = ~Nocturia, type = "scatter", mode = "lines",
            name = "Nocturia", line = list(color = "#E74C3C", width = 2)) %>%
      layout(xaxis = list(title = "Day"), yaxis = list(title = "Nocturnal Voids"))
  })

  output$cont_plot <- renderPlotly({
    df <- sim_data()
    plot_ly(df, x = ~time_day, y = ~ContScore, type = "scatter", mode = "lines",
            name = "Continence Score", line = list(color = "#27AE60", width = 2)) %>%
      layout(xaxis = list(title = "Day"), yaxis = list(title = "Score (0=incontinent, 100=continent)"))
  })
}

# ────────────────────────────────────────────────────────────────
shinyApp(ui = ui, server = server)
