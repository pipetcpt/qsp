## ============================================================
## Alzheimer's Disease QSP — Interactive Shiny Dashboard
## Tabs: Patient Profile · Drug PK · Amyloid/Tau Biomarkers ·
##       Neuroinflammation & Synapse · Cognitive Endpoints ·
##       Scenario Comparison
## ============================================================

library(shiny)
library(shinydashboard)
library(ggplot2)
library(dplyr)
library(tidyr)
library(scales)
library(DT)

## ── INLINE QSP SIMULATION (simplified for Shiny, no mrgsolve dep.) ──
## Uses Euler integration for interactive speed
simulate_ad <- function(
    duration_years = 3,
    dt_h = 24,             # daily steps
    APOE4   = 0,           # 0/1/2
    FAD     = 0,           # 0/1
    # Drug doses (0 = off)
    don_dose   = 0,        # mg/day donepezil  (0/5/10)
    mem_dose   = 0,        # mg/day memantine  (0/20)
    lec_dose   = 0,        # mg/kg biweekly lecanemab (0/10)
    age_yrs    = 72,
    # initial disease severity
    init_mmse  = 22
) {
  steps <- round(duration_years * 365 * 24 / dt_h)
  t_h   <- seq(0, by = dt_h, length.out = steps + 1)

  ## ── Parameters ──────────────────────────────────────────────
  ApoeScale <- 1 - 0.40 * APOE4
  FadScale  <- 1 + 0.80 * FAD
  AgeScale  <- 1 + 0.012 * (age_yrs - 70)

  # PK simple 1-cpt approximations
  KA_DON <- 0.0115; V_DON <- 594; CL_DON <- 13.3; Kp_DON <- 0.18
  V_MEM  <- 520;    CL_MEM <- 8.4; Kp_MEM <- 2.5
  V_LEC  <- 3.2;    CL_LEC <- 0.0167; kBBB_LEC <- 8.3e-5; kCNS_LEC <- 0.0042

  IC50_DON <- 6.7e-6; IC50_MEM <- 0.8; KD_LEC <- 5.4e-5

  kprod_Ab <- 0.042 * FadScale * AgeScale
  kdeg_Ab  <- 0.037 * ApoeScale
  kn_oligo <- 0.0018; ke_oligo <- 0.0045; ke_fibril <- 0.0012
  kdep_plaq <- 0.00055; kphago <- 0.0003
  kLRP1 <- 0.014

  kprod_tau <- 0.028 * AgeScale; kphos_tau <- 0.0022; kAb_tau <- 0.0055
  kdphos_tau <- 0.0019; kagg_tau <- 0.0009; kNFT <- 0.00035
  kdeg_tau <- 0.016; kspread <- 0.000045

  kact_micro <- 0.0065; kinact_micro <- 0.0045

  ksynth_ACh <- 0.15; kdeg_ACh_base <- 0.12; kACh_BFCN <- 0.008
  ACh_ss <- 1.0

  kloss_Ab <- 0.0012; kloss_tau <- 0.00085; kloss_inf <- 0.00055
  kregen <- 0.00018

  wSyn <- 0.50; wACh <- 0.25; wTau <- 0.15; wAb <- 0.10

  ## ── Initial States ─────────────────────────────────────────
  severity_shift <- (28 - init_mmse) / 18  # 0=mild, 1=severe
  st <- list(
    DON_GUT = 0, DON_CENT = 0, DON_CNS = 0,
    MEM_CENT = 0, MEM_CNS = 0,
    LEC_CENT = 0, LEC_PERI = 0, LEC_CNS = 0,
    AB_MONO  = 0.8 + 0.4 * severity_shift,
    AB_OLIGO = 0.5 + 0.3 * severity_shift,
    AB_PROTO = 0.3 + 0.2 * severity_shift,
    AB_PLAQUE= 2.5 + 2.0 * severity_shift,
    TAU_SOL  = 1.0,
    TAU_PHOS = 1.5 + 0.8 * severity_shift,
    TAU_AGG  = 0.8 + 0.5 * severity_shift,
    NEURO    = 0.4 + 0.2 * severity_shift,
    ACH      = max(0.3, 0.65 - 0.25 * severity_shift),
    SYN      = max(0.3, 0.72 - 0.25 * severity_shift),
    COG      = init_mmse
  )

  ## ── Storage ─────────────────────────────────────────────────
  n_out <- length(t_h)
  out <- data.frame(
    time_y    = t_h / 8760,
    MMSE      = numeric(n_out),
    ADAS_Cog  = numeric(n_out),
    CDR_SB    = numeric(n_out),
    CSF_Ab42  = numeric(n_out),
    CSF_pTau  = numeric(n_out),
    AmylPET   = numeric(n_out),
    TauPET    = numeric(n_out),
    NfL       = numeric(n_out),
    AB_PLAQUE = numeric(n_out),
    AB_OLIGO  = numeric(n_out),
    TAU_AGG   = numeric(n_out),
    NEURO     = numeric(n_out),
    SYN       = numeric(n_out),
    ACH       = numeric(n_out),
    Cp_DON    = numeric(n_out),
    Cb_DON    = numeric(n_out),
    Cp_MEM    = numeric(n_out),
    Cp_LEC    = numeric(n_out),
    AChE_inh  = numeric(n_out),
    NMDA_occ  = numeric(n_out)
  )

  ## Dosing schedule
  dose_don <- don_dose  # daily mg
  dose_mem <- mem_dose  # daily mg
  bw_kg    <- 70        # assume 70 kg
  dose_lec_mg <- lec_dose * bw_kg  # biweekly dose in mg

  record <- function(i, s) {
    Cp_don <- s$DON_CENT / V_DON
    inh_ache <- (Cp_don^1.2) / (IC50_DON^1.2 + Cp_don^1.2)
    Cp_mem <- s$MEM_CENT / V_MEM
    inh_nmda <- Cp_mem / (IC50_MEM + Cp_mem)
    inh_lec <- s$LEC_CNS / (KD_LEC + s$LEC_CNS)

    mmse <- s$COG
    adas <- max(0, 70 - 2.33 * mmse)
    cdr  <- max(0, 18 * (1 - mmse / 30))

    out[i, "MMSE"]     <<- mmse
    out[i, "ADAS_Cog"] <<- adas
    out[i, "CDR_SB"]   <<- cdr
    out[i, "CSF_Ab42"] <<- 1200 * exp(-0.35 * s$AB_PLAQUE)
    out[i, "CSF_pTau"] <<- 15 + 28 * s$TAU_PHOS
    out[i, "AmylPET"]  <<- 20 + 14 * s$AB_PLAQUE
    out[i, "TauPET"]   <<- 1 + 0.18 * s$TAU_AGG
    out[i, "NfL"]      <<- 10 + 22 * (1 - s$SYN)
    out[i, "AB_PLAQUE"]<<- s$AB_PLAQUE
    out[i, "AB_OLIGO"] <<- s$AB_OLIGO
    out[i, "TAU_AGG"]  <<- s$TAU_AGG
    out[i, "NEURO"]    <<- s$NEURO
    out[i, "SYN"]      <<- s$SYN
    out[i, "ACH"]      <<- s$ACH
    out[i, "Cp_DON"]   <<- Cp_don * 1000
    out[i, "Cb_DON"]   <<- s$DON_CNS * 1000
    out[i, "Cp_MEM"]   <<- Cp_mem * 1000
    out[i, "Cp_LEC"]   <<- s$LEC_CENT / V_LEC * 1000
    out[i, "AChE_inh"] <<- inh_ache * 100
    out[i, "NMDA_occ"] <<- inh_nmda * 100
  }

  record(1, st)
  lec_day_counter <- 0

  for (i in seq_len(steps)) {
    s <- st

    # Dosing events
    if (dose_don > 0) s$DON_GUT <- s$DON_GUT + dose_don / (24 / dt_h)
    if (dose_mem > 0) s$MEM_CENT <- s$MEM_CENT + dose_mem / (24 / dt_h)
    # Lecanemab biweekly
    lec_day_counter <- lec_day_counter + dt_h / 24
    if (dose_lec_mg > 0 && lec_day_counter >= 14) {
      s$LEC_CENT <- s$LEC_CENT + dose_lec_mg
      lec_day_counter <- 0
    }

    # Current PK-derived quantities
    Cp_don   <- s$DON_CENT / V_DON
    inh_ache <- (Cp_don^1.2) / (IC50_DON^1.2 + Cp_don^1.2)
    Cp_mem   <- s$MEM_CENT / V_MEM
    inh_nmda <- Cp_mem / (IC50_MEM + Cp_mem)
    inh_lec  <- s$LEC_CNS / (KD_LEC + s$LEC_CNS)
    kdeg_eff <- kdeg_ACh_base * (1 - inh_ache)

    # Amyloid
    dAB_MONO  <- kprod_Ab - kdeg_Ab * s$AB_MONO - kn_oligo * s$AB_MONO^2 - kLRP1 * s$AB_MONO
    dAB_OLIGO <- kn_oligo * s$AB_MONO^2 - ke_oligo * s$AB_OLIGO - 0.5 * kdeg_Ab * s$AB_OLIGO
    dAB_PROTO <- ke_oligo * s$AB_OLIGO - ke_fibril * s$AB_PROTO - inh_lec * ke_fibril * s$AB_PROTO
    kphago_eff <- kphago * (1 - 0.5 * s$NEURO)
    dAB_PLAQUE<- ke_fibril * (1 - inh_lec) * s$AB_PROTO - kphago_eff * s$AB_PLAQUE

    # Tau
    kphos_eff <- kphos_tau + kAb_tau * s$AB_OLIGO
    dTAU_SOL  <- kprod_tau - kphos_eff * s$TAU_SOL - kdeg_tau * s$TAU_SOL
    dTAU_PHOS <- kphos_eff * s$TAU_SOL - kdphos_tau * s$TAU_PHOS - kagg_tau * s$TAU_PHOS - 0.2 * kdeg_tau * s$TAU_PHOS
    dTAU_AGG  <- kagg_tau * s$TAU_PHOS + kspread * s$TAU_AGG * (1 - s$TAU_AGG / 20) - kNFT * s$TAU_AGG

    # Neuroinflam
    dNEURO <- kact_micro * s$AB_PLAQUE * (1 - s$NEURO) + 8e-4 * s$TAU_AGG * (1 - s$NEURO) - kinact_micro * s$NEURO

    # ACh
    BFCN_loss  <- kACh_BFCN * (s$AB_PLAQUE * 0.5 + s$TAU_AGG * 0.5)
    dACH <- ksynth_ACh * max(0, 1 - BFCN_loss) - kdeg_eff * s$ACH

    # Synapse
    kregen_eff <- kregen * (s$ACH / ACh_ss)
    dSYN <- kregen_eff * s$SYN * (1 - s$SYN) - kloss_Ab * s$AB_OLIGO * s$SYN - kloss_tau * s$TAU_PHOS * s$SYN - kloss_inf * s$NEURO * s$SYN + 4e-4 * inh_nmda * s$SYN

    # Cognition
    nSyn <- max(0, min(1, s$SYN))
    nACh <- max(0, min(2, s$ACH / ACh_ss))
    nTau <- max(0, min(1, s$TAU_AGG / 5))
    nAb  <- max(0, min(1, s$AB_PLAQUE / 10))
    MMSE_tgt <- 30 * (wSyn * nSyn + wACh * nACh * 0.5 - wTau * nTau - wAb * nAb)
    MMSE_tgt <- max(0, MMSE_tgt)
    dCOG <- 2.5e-4 * (MMSE_tgt - s$COG)

    # Euler update
    dt <- dt_h
    st$AB_MONO  <- max(0, s$AB_MONO  + dAB_MONO  * dt)
    st$AB_OLIGO <- max(0, s$AB_OLIGO + dAB_OLIGO * dt)
    st$AB_PROTO <- max(0, s$AB_PROTO + dAB_PROTO * dt)
    st$AB_PLAQUE<- max(0, s$AB_PLAQUE+ dAB_PLAQUE* dt)
    st$TAU_SOL  <- max(0, s$TAU_SOL  + dTAU_SOL  * dt)
    st$TAU_PHOS <- max(0, s$TAU_PHOS + dTAU_PHOS * dt)
    st$TAU_AGG  <- max(0, s$TAU_AGG  + dTAU_AGG  * dt)
    st$NEURO    <- max(0, min(1, s$NEURO + dNEURO * dt))
    st$ACH      <- max(0.01, s$ACH   + dACH      * dt)
    st$SYN      <- max(0.01, min(1, s$SYN + dSYN * dt))
    st$COG      <- max(0, min(30, s$COG   + dCOG  * dt))

    # PK updates (simple)
    st$DON_GUT   <- max(0, s$DON_GUT  - 0.0115 * s$DON_GUT  * dt)
    st$DON_CENT  <- max(0, s$DON_CENT + (0.0115 * s$DON_GUT - (CL_DON / V_DON) * s$DON_CENT) * dt)
    st$DON_CNS   <- max(0, s$DON_CNS  + 0.5 * (Kp_DON * s$DON_CENT / V_DON - s$DON_CNS) * dt)
    st$MEM_CENT  <- max(0, s$MEM_CENT - (CL_MEM / V_MEM) * s$MEM_CENT * dt)
    st$MEM_CNS   <- max(0, s$MEM_CNS  + 0.3 * (Kp_MEM * s$MEM_CENT / V_MEM - s$MEM_CNS) * dt)
    st$LEC_CENT  <- max(0, s$LEC_CENT - (CL_LEC / V_LEC + kBBB_LEC) * s$LEC_CENT * dt)
    st$LEC_CNS   <- max(0, s$LEC_CNS  + (kBBB_LEC * s$LEC_CENT - kCNS_LEC * s$LEC_CNS) * dt)

    record(i + 1, st)
  }
  out
}

## ═══════════════════════════════════════════════════════════════
## UI
## ═══════════════════════════════════════════════════════════════
ui <- dashboardPage(
  skin = "purple",
  dashboardHeader(title = "AD QSP Dashboard"),

  dashboardSidebar(
    sidebarMenu(
      menuItem("Patient Profile",     tabName = "patient",  icon = icon("user-md")),
      menuItem("Drug PK",             tabName = "pk",       icon = icon("pills")),
      menuItem("Amyloid & Tau",       tabName = "biomarker",icon = icon("vials")),
      menuItem("Neuro & Synapse",     tabName = "neuro",    icon = icon("brain")),
      menuItem("Cognitive Endpoints", tabName = "cognition",icon = icon("chart-line")),
      menuItem("Scenario Comparison", tabName = "scenario", icon = icon("sliders-h"))
    ),
    hr(),
    h5("  Patient & Disease", style = "padding-left:15px; color:#ddd;"),
    sliderInput("age",       "Age (years)",          60, 90, 72, step = 1),
    sliderInput("init_mmse", "Baseline MMSE",         5, 28, 22, step = 1),
    selectInput("apoe4",     "APOE ε4 Status",
                choices = list("Non-carrier (ε3/ε3)" = 0,
                               "Heterozygous (ε3/ε4)" = 1,
                               "Homozygous (ε4/ε4)" = 2)),
    checkboxInput("fad",     "FAD Mutation (PSEN1/APP)", FALSE),
    sliderInput("duration",  "Simulation Duration (years)", 1, 5, 3, step = 0.5),
    hr(),
    h5("  Drug Treatment", style = "padding-left:15px; color:#ddd;"),
    selectInput("donepezil", "Donepezil",
                choices = list("None" = 0, "5 mg/day" = 5, "10 mg/day" = 10)),
    checkboxInput("memantine","Memantine 20 mg/day", FALSE),
    selectInput("lecanemab", "Lecanemab (mg/kg biweekly)",
                choices = list("None" = 0, "2.5 mg/kg" = 2.5, "5 mg/kg" = 5, "10 mg/kg" = 10)),
    actionButton("simulate", "Run Simulation", icon = icon("play"),
                 style = "background:#8A2BE2; color:white; margin:10px 15px; width:180px;")
  ),

  dashboardBody(
    tags$head(tags$style(HTML("
      .skin-purple .sidebar { background-color: #2d1a42; }
      .box-header { background-color: #6a1b9a !important; color: white !important; }
      .small-box.bg-purple { background-color: #8e24aa !important; }
    "))),

    tabItems(

      ## ── TAB 1: PATIENT PROFILE ──────────────────────────────
      tabItem("patient",
        fluidRow(
          valueBoxOutput("vb_mmse",    width = 3),
          valueBoxOutput("vb_stage",   width = 3),
          valueBoxOutput("vb_apoe4",   width = 3),
          valueBoxOutput("vb_therapy", width = 3)
        ),
        fluidRow(
          box(title = "Disease Stage Characterization", width = 6, status = "primary",
            tableOutput("stage_table")
          ),
          box(title = "ATN Biomarker Framework", width = 6, status = "warning",
            plotOutput("atn_radar", height = 280)
          )
        ),
        fluidRow(
          box(title = "Model Parameters Summary", width = 12, status = "info",
            DTOutput("param_table")
          )
        )
      ),

      ## ── TAB 2: DRUG PK ──────────────────────────────────────
      tabItem("pk",
        fluidRow(
          box(title = "Donepezil PK — Plasma & CNS", width = 6, status = "primary",
              plotOutput("pk_don", height = 300)),
          box(title = "Memantine PK — Plasma (if active)", width = 6, status = "info",
              plotOutput("pk_mem", height = 300))
        ),
        fluidRow(
          box(title = "Lecanemab PK — Plasma & CNS", width = 6, status = "warning",
              plotOutput("pk_lec", height = 300)),
          box(title = "Pharmacodynamic Occupancy", width = 6, status = "success",
              plotOutput("pk_pd",  height = 300))
        )
      ),

      ## ── TAB 3: AMYLOID & TAU ────────────────────────────────
      tabItem("biomarker",
        fluidRow(
          box(title = "Amyloid PET (Centiloids)", width = 6, status = "danger",
              plotOutput("amyl_pet", height = 280)),
          box(title = "Tau PET (SUVr)", width = 6, status = "warning",
              plotOutput("tau_pet",  height = 280))
        ),
        fluidRow(
          box(title = "CSF Aβ42 (pg/mL)", width = 6, status = "primary",
              plotOutput("csf_ab", height = 280)),
          box(title = "CSF p-Tau181 (pg/mL)", width = 6, status = "warning",
              plotOutput("csf_tau", height = 280))
        ),
        fluidRow(
          box(title = "Plasma NfL (pg/mL) — Neurodegeneration", width = 6, status = "info",
              plotOutput("nfl", height = 280)),
          box(title = "Amyloid Cascade State Variables", width = 6, status = "danger",
              plotOutput("ab_pool", height = 280))
        )
      ),

      ## ── TAB 4: NEUROINFLAMMATION & SYNAPSE ──────────────────
      tabItem("neuro",
        fluidRow(
          box(title = "Neuroinflammation Index (0–1)", width = 6, status = "danger",
              plotOutput("neuro_inflam", height = 300)),
          box(title = "Synaptic Integrity (0–1, normalized)", width = 6, status = "success",
              plotOutput("synapse", height = 300))
        ),
        fluidRow(
          box(title = "Cholinergic Tone (ACh, normalized)", width = 6, status = "info",
              plotOutput("ach_plot", height = 300)),
          box(title = "Aggregated Tau — NFT Burden", width = 6, status = "warning",
              plotOutput("tau_agg_plot", height = 300))
        )
      ),

      ## ── TAB 5: COGNITIVE ENDPOINTS ──────────────────────────
      tabItem("cognition",
        fluidRow(
          box(title = "MMSE Trajectory (0–30, higher = better)", width = 6,
              status = "primary", solidHeader = TRUE,
              plotOutput("mmse_plot", height = 320)),
          box(title = "CDR-SB Trajectory (0–18, lower = better)", width = 6,
              status = "warning", solidHeader = TRUE,
              plotOutput("cdr_plot", height = 320))
        ),
        fluidRow(
          box(title = "ADAS-Cog13 Trajectory (0–70, lower = better)", width = 6,
              status = "danger",
              plotOutput("adas_plot", height = 280)),
          box(title = "6-Month Summary Statistics", width = 6, status = "info",
              DTOutput("summary_tbl"))
        )
      ),

      ## ── TAB 6: SCENARIO COMPARISON ──────────────────────────
      tabItem("scenario",
        fluidRow(
          box(title = "Scenario Settings", width = 3, status = "primary",
            checkboxGroupInput("scen_sel", "Select Scenarios:",
              choices = list(
                "No Treatment"          = "No Treatment",
                "Donepezil 10 mg"       = "Donepezil 10mg",
                "Memantine 20 mg"       = "Memantine 20mg",
                "Combo (Don + Mem)"     = "Don+Mem Combo",
                "Lecanemab 10 mg/kg"    = "Lecanemab 10mpk",
                "Leca + Donepezil"      = "Leca+Don Combo"
              ),
              selected = c("No Treatment","Donepezil 10mg","Lecanemab 10mpk")
            ),
            selectInput("scen_endpoint", "Primary Endpoint:",
              choices = list(
                "MMSE"         = "MMSE",
                "CDR-SB"       = "CDR_SB",
                "ADAS-Cog"     = "ADAS_Cog",
                "Amyloid PET"  = "AmylPET",
                "CSF Aβ42"     = "CSF_Ab42",
                "CSF p-Tau181" = "CSF_pTau",
                "NfL"          = "NfL"
              )),
            actionButton("run_scenarios", "Compare Scenarios",
                         icon = icon("play"),
                         style = "background:#8A2BE2; color:white; width:100%;")
          ),
          box(title = "Scenario Comparison Plot", width = 9, status = "primary",
              solidHeader = TRUE,
              plotOutput("scen_plot", height = 450))
        ),
        fluidRow(
          box(title = "Scenario Comparison Table (36-month endpoint)",
              width = 12, status = "info",
              DTOutput("scen_table"))
        )
      )
    )
  )
)

## ═══════════════════════════════════════════════════════════════
## SERVER
## ═══════════════════════════════════════════════════════════════
server <- function(input, output, session) {

  ## Reactive simulation result
  sim_data <- eventReactive(input$simulate, {
    withProgress(message = "Simulating...", {
      simulate_ad(
        duration_years = input$duration,
        APOE4  = as.numeric(input$apoe4),
        FAD    = if(input$fad) 1 else 0,
        don_dose = as.numeric(input$donepezil),
        mem_dose = if(input$memantine) 20 else 0,
        lec_dose = as.numeric(input$lecanemab),
        age_yrs  = input$age,
        init_mmse= input$init_mmse
      )
    })
  }, ignoreNULL = FALSE)

  ## ── TAB 1: VALUE BOXES ───────────────────────────────────────
  output$vb_mmse <- renderValueBox({
    df <- sim_data(); mmse <- round(df$MMSE[1])
    valueBox(mmse, "Baseline MMSE", icon = icon("brain"),
             color = ifelse(mmse >= 20, "green", ifelse(mmse >= 10, "yellow", "red")))
  })
  output$vb_stage <- renderValueBox({
    mmse <- input$init_mmse
    stage <- if (mmse >= 26) "Preclinical/MCI"
             else if (mmse >= 20) "Mild AD"
             else if (mmse >= 10) "Moderate AD"
             else "Severe AD"
    valueBox(stage, "Disease Stage", icon = icon("stethoscope"), color = "purple")
  })
  output$vb_apoe4 <- renderValueBox({
    apoe4_label <- c("ε3/ε3","ε3/ε4","ε4/ε4")[as.numeric(input$apoe4) + 1]
    valueBox(apoe4_label, "APOE Genotype", icon = icon("dna"),
             color = ifelse(as.numeric(input$apoe4) == 0, "green", "red"))
  })
  output$vb_therapy <- renderValueBox({
    drugs <- c(
      if(as.numeric(input$donepezil) > 0) "Donepezil",
      if(input$memantine) "Memantine",
      if(as.numeric(input$lecanemab) > 0) "Lecanemab"
    )
    lab <- if(length(drugs) == 0) "Untreated" else paste(drugs, collapse = " + ")
    valueBox(lab, "Current Therapy", icon = icon("pills"),
             color = if(length(drugs) == 0) "red" else "blue")
  })

  output$stage_table <- renderTable({
    data.frame(
      Stage = c("Preclinical AD","MCI due to AD","Mild AD","Moderate AD","Severe AD"),
      MMSE  = c("27–30","20–26","20–26","10–19","< 10"),
      CDR   = c("0","0.5","1","2","3"),
      Duration_Approx = c("~15 yr","~2–4 yr","~2–3 yr","~2–3 yr","~1–3 yr"),
      Key_Biomarker = c("Aβ+, tau−","Aβ+, tau+","NFTs entorhinal","Neocortical NFT","Widespread NFT")
    )
  })

  output$param_table <- renderDT({
    df <- data.frame(
      Parameter = c("kprod_Ab","kdeg_Ab","IC50_DON","IC50_MEM","KD_LEC",
                    "kphos_tau","kspread_tau","kact_micro","ksynth_ACh","kloss_syn"),
      Value = c(0.042,0.037,"6.7e-6",0.8,"5.4e-5",0.0022,"4.5e-5",0.0065,0.15,0.0012),
      Unit  = c("1/h","1/h","mg/L","mg/L","mg/L","1/h","1/h","1/h","1/h","1/h"),
      Reference = c("Bateman 2006","Mawuenyega 2010","Tiseo 1998","Parsons 1999",
                    "Swanson 2021","Iqbal 2016","Clavaguera 2009","Heneka 2015",
                    "Whitehouse 1982","Terry 1991")
    )
    datatable(df, options = list(pageLength = 10, dom = "t"), rownames = FALSE)
  })

  output$atn_radar <- renderPlot({
    df <- sim_data()
    end_row <- df[nrow(df), ]
    # Normalize to 0–10 scale for display
    atn <- data.frame(
      Component = c("Aβ (PET)", "p-Tau (CSF)", "NfL (N)", "Plaque Burden", "Syn Loss"),
      Score = c(
        min(10, end_row$AmylPET / 12),
        min(10, end_row$CSF_pTau / 10),
        min(10, end_row$NfL / 8),
        min(10, end_row$AB_PLAQUE * 1.5),
        min(10, (1 - end_row$SYN) * 10)
      )
    )
    ggplot(atn, aes(x = reorder(Component, Score), y = Score, fill = Score)) +
      geom_col(show.legend = FALSE) +
      scale_fill_gradient(low = "#90EE90", high = "#CC0000") +
      coord_flip() +
      labs(title = "ATN Biomarker Burden at End of Simulation",
           x = NULL, y = "Burden Score (0–10)") +
      theme_minimal(base_size = 11)
  })

  ## ── TAB 2: DRUG PK ───────────────────────────────────────────
  pk_theme <- function() theme_bw(base_size = 11) + theme(legend.position = "none")

  output$pk_don <- renderPlot({
    df <- sim_data()
    df %>% select(time_y, Cp_DON, Cb_DON) %>%
      pivot_longer(-time_y, names_to = "compartment") %>%
      ggplot(aes(x = time_y, y = value, color = compartment)) +
      geom_line(linewidth = 1.0) +
      scale_color_manual(values = c(Cp_DON = "#0044CC", Cb_DON = "#CC4400"),
                         labels = c(Cp_DON = "Plasma", Cb_DON = "CNS")) +
      labs(x = "Time (years)", y = "Concentration (ng/mL)",
           color = "Compartment") + pk_theme() + theme(legend.position = "top")
  })

  output$pk_mem <- renderPlot({
    df <- sim_data()
    ggplot(df, aes(x = time_y, y = Cp_MEM)) +
      geom_line(color = "#009900", linewidth = 1.1) +
      labs(x = "Time (years)", y = "Plasma Conc (ng/mL)",
           subtitle = if(input$memantine) "Active" else "Not Administered") +
      pk_theme()
  })

  output$pk_lec <- renderPlot({
    df <- sim_data()
    ggplot(df, aes(x = time_y, y = Cp_LEC)) +
      geom_line(color = "#990099", linewidth = 1.1) +
      labs(x = "Time (years)", y = "Plasma Conc (μg/mL)") + pk_theme()
  })

  output$pk_pd <- renderPlot({
    df <- sim_data()
    df %>% select(time_y, AChE_inh, NMDA_occ) %>%
      pivot_longer(-time_y) %>%
      ggplot(aes(x = time_y, y = value, color = name)) +
      geom_line(linewidth = 1.0) +
      geom_hline(yintercept = 70, linetype = "dashed", color = "#880000") +
      scale_color_manual(values = c(AChE_inh = "#0044CC", NMDA_occ = "#009900"),
                         labels = c(AChE_inh = "AChE Inhibition",
                                    NMDA_occ = "NMDA Occupancy")) +
      labs(x = "Time (years)", y = "Occupancy / Inhibition (%)", color = NULL) +
      theme_bw(base_size = 11) + theme(legend.position = "top")
  })

  ## ── TAB 3: AMYLOID & TAU ─────────────────────────────────────
  output$amyl_pet <- renderPlot({
    df <- sim_data()
    ggplot(df, aes(x = time_y, y = AmylPET)) +
      geom_line(color = "#CC0000", linewidth = 1.2) +
      geom_hline(yintercept = 24, linetype = "dashed", color = "gray40") +
      annotate("text", x = max(df$time_y) * 0.6, y = 26,
               label = "Positivity cutoff (24 CL)", size = 3) +
      labs(x = "Time (years)", y = "Centiloids") + theme_bw(base_size = 11)
  })

  output$tau_pet <- renderPlot({
    df <- sim_data()
    ggplot(df, aes(x = time_y, y = TauPET)) +
      geom_line(color = "#FF6600", linewidth = 1.2) +
      geom_hline(yintercept = 1.33, linetype = "dashed", color = "gray40") +
      labs(x = "Time (years)", y = "SUVr") + theme_bw(base_size = 11)
  })

  output$csf_ab <- renderPlot({
    df <- sim_data()
    ggplot(df, aes(x = time_y, y = CSF_Ab42)) +
      geom_line(color = "#0044AA", linewidth = 1.2) +
      geom_hline(yintercept = 1000, linetype = "dashed", color = "red") +
      annotate("text", x = max(df$time_y) * 0.6, y = 1050,
               label = "Cutoff ~1000 pg/mL", size = 3, color = "red") +
      labs(x = "Time (years)", y = "CSF Aβ42 (pg/mL)") + theme_bw(base_size = 11)
  })

  output$csf_tau <- renderPlot({
    df <- sim_data()
    ggplot(df, aes(x = time_y, y = CSF_pTau)) +
      geom_line(color = "#884400", linewidth = 1.2) +
      geom_hline(yintercept = 23, linetype = "dashed", color = "red") +
      labs(x = "Time (years)", y = "CSF p-Tau181 (pg/mL)") + theme_bw(base_size = 11)
  })

  output$nfl <- renderPlot({
    df <- sim_data()
    ggplot(df, aes(x = time_y, y = NfL)) +
      geom_line(color = "#555555", linewidth = 1.2) +
      labs(x = "Time (years)", y = "Plasma NfL (pg/mL)") + theme_bw(base_size = 11)
  })

  output$ab_pool <- renderPlot({
    df <- sim_data()
    df %>% select(time_y, AB_MONO, AB_OLIGO, AB_PROTO, AB_PLAQUE) %>%
      pivot_longer(-time_y, names_to = "pool") %>%
      ggplot(aes(x = time_y, y = value, color = pool)) +
      geom_line(linewidth = 0.9) +
      scale_color_manual(values = c(AB_MONO = "#FFAA55", AB_OLIGO = "#FF6600",
                                    AB_PROTO = "#CC2200", AB_PLAQUE = "#880000")) +
      labs(x = "Time (years)", y = "Amyloid Pool (normalized)", color = "Pool") +
      theme_bw(base_size = 11) + theme(legend.position = "top")
  })

  ## ── TAB 4: NEURO & SYNAPSE ───────────────────────────────────
  output$neuro_inflam <- renderPlot({
    df <- sim_data()
    ggplot(df, aes(x = time_y, y = NEURO)) +
      geom_line(color = "#CC0044", linewidth = 1.2) +
      geom_ribbon(aes(ymin = 0, ymax = NEURO), fill = "#CC0044", alpha = 0.2) +
      scale_y_continuous(limits = c(0, 1)) +
      labs(x = "Time (years)", y = "Neuroinflammation Index") + theme_bw(base_size = 11)
  })

  output$synapse <- renderPlot({
    df <- sim_data()
    ggplot(df, aes(x = time_y, y = SYN)) +
      geom_line(color = "#008800", linewidth = 1.2) +
      geom_ribbon(aes(ymin = 0, ymax = SYN), fill = "#008800", alpha = 0.15) +
      scale_y_continuous(limits = c(0, 1)) +
      labs(x = "Time (years)", y = "Synaptic Integrity (0–1)") + theme_bw(base_size = 11)
  })

  output$ach_plot <- renderPlot({
    df <- sim_data()
    ggplot(df, aes(x = time_y, y = ACH)) +
      geom_line(color = "#0044CC", linewidth = 1.2) +
      labs(x = "Time (years)", y = "ACh (normalized)") + theme_bw(base_size = 11)
  })

  output$tau_agg_plot <- renderPlot({
    df <- sim_data()
    ggplot(df, aes(x = time_y, y = TAU_AGG)) +
      geom_line(color = "#884400", linewidth = 1.2) +
      labs(x = "Time (years)", y = "Tau Aggregate (normalized)") + theme_bw(base_size = 11)
  })

  ## ── TAB 5: COGNITIVE ENDPOINTS ───────────────────────────────
  output$mmse_plot <- renderPlot({
    df <- sim_data()
    ggplot(df, aes(x = time_y, y = MMSE)) +
      geom_line(color = "#4400CC", linewidth = 1.3) +
      geom_hline(yintercept = c(10, 20, 26), linetype = "dashed",
                 color = c("red","orange","blue"), alpha = 0.7) +
      annotate("text", x = max(df$time_y), y = c(11, 21, 27),
               label = c("Severe","Moderate","MCI"), hjust = 1,
               color = c("red","orange","blue"), size = 3) +
      scale_y_continuous(limits = c(0, 30)) +
      labs(x = "Time (years)", y = "MMSE (0–30)") + theme_bw(base_size = 12)
  })

  output$cdr_plot <- renderPlot({
    df <- sim_data()
    ggplot(df, aes(x = time_y, y = CDR_SB)) +
      geom_line(color = "#CC8800", linewidth = 1.3) +
      geom_hline(yintercept = c(4.5, 9), linetype = "dashed",
                 color = c("orange","red"), alpha = 0.7) +
      labs(x = "Time (years)", y = "CDR-SB (0–18)") + theme_bw(base_size = 12)
  })

  output$adas_plot <- renderPlot({
    df <- sim_data()
    ggplot(df, aes(x = time_y, y = ADAS_Cog)) +
      geom_line(color = "#CC0000", linewidth = 1.2) +
      labs(x = "Time (years)", y = "ADAS-Cog (0–70)") + theme_bw(base_size = 11)
  })

  output$summary_tbl <- renderDT({
    df <- sim_data()
    checkpoints <- c(0, 0.5, 1, 1.5, 2, 2.5, 3)
    checkpoints <- checkpoints[checkpoints <= input$duration]
    tbl <- purrr::map_dfr(checkpoints, function(yr) {
      row <- df[which.min(abs(df$time_y - yr)), ]
      data.frame(
        `Time (yr)` = yr,
        MMSE        = round(row$MMSE, 1),
        `ADAS-Cog`  = round(row$ADAS_Cog, 1),
        `CDR-SB`    = round(row$CDR_SB, 1),
        `PET (CL)`  = round(row$AmylPET, 0),
        SYN         = round(row$SYN, 2),
        check.names = FALSE
      )
    })
    datatable(tbl, options = list(dom = "t", pageLength = 10), rownames = FALSE)
  })

  ## ── TAB 6: SCENARIO COMPARISON ───────────────────────────────
  scenario_data <- eventReactive(input$run_scenarios, {
    scen_map <- list(
      "No Treatment"      = list(don=0,  mem=FALSE, lec=0),
      "Donepezil 10mg"    = list(don=10, mem=FALSE, lec=0),
      "Memantine 20mg"    = list(don=0,  mem=TRUE,  lec=0),
      "Don+Mem Combo"     = list(don=10, mem=TRUE,  lec=0),
      "Lecanemab 10mpk"   = list(don=0,  mem=FALSE, lec=10),
      "Leca+Don Combo"    = list(don=10, mem=FALSE, lec=10)
    )
    sels <- input$scen_sel
    withProgress(message = "Running scenarios...", {
      purrr::map_dfr(sels, function(nm) {
        p <- scen_map[[nm]]
        d <- simulate_ad(
          duration_years = input$duration,
          APOE4 = as.numeric(input$apoe4),
          FAD   = if(input$fad) 1 else 0,
          don_dose  = p$don, mem_dose = if(p$mem) 20 else 0,
          lec_dose  = p$lec, age_yrs = input$age,
          init_mmse = input$init_mmse
        )
        d$scenario <- nm
        d
      })
    })
  })

  output$scen_plot <- renderPlot({
    df <- scenario_data()
    ep <- input$scen_endpoint
    ggplot(df, aes_string(x = "time_y", y = ep, color = "scenario")) +
      geom_line(linewidth = 1.2) +
      scale_color_brewer(palette = "Dark2") +
      labs(x = "Time (years)", y = ep, color = "Scenario",
           title = paste("Scenario Comparison —", ep)) +
      theme_bw(base_size = 13) +
      theme(legend.position = "bottom", legend.text = element_text(size = 9))
  })

  output$scen_table <- renderDT({
    df <- scenario_data()
    ep <- input$scen_endpoint
    tbl <- df %>%
      group_by(scenario) %>%
      summarise(
        Start  = round(.data[[ep]][1], 2),
        End    = round(.data[[ep]][n()], 2),
        Change = round(.data[[ep]][n()] - .data[[ep]][1], 2),
        `Pct Change` = round(100 * (.data[[ep]][n()] - .data[[ep]][1]) /
                               max(.001, .data[[ep]][1]), 1),
        .groups = "drop"
      )
    datatable(tbl, options = list(dom = "t"), rownames = FALSE)
  })
}

## ── RUN ────────────────────────────────────────────────────────
shinyApp(ui = ui, server = server)
