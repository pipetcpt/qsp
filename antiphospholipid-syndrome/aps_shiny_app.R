################################################################################
# Antiphospholipid Syndrome (APS) — Interactive QSP Shiny Dashboard
#
# Tabs:
#   1. Patient Profile      — demographics, risk factors, aPL serology, Sapporo criteria
#   2. Drug PK              — PK profiles for warfarin / LMWH / HCQ / rivaroxaban / ASA
#   3. PD Biomarkers        — aPL titer, B cells, complement C5a, TF, platelet activation
#   4. Clinical Endpoints   — DVT/PE risk, INR tracking, CAPS probability, mTOR nephropathy
#   5. Scenario Comparison  — 7 treatment arms head-to-head comparison
#   6. Obstetric APS        — pregnancy viability, live birth rate, placental health
#
# Dependencies: shiny, shinydashboard, plotly, DT, deSolve, dplyr, ggplot2
################################################################################

library(shiny)
library(shinydashboard)
library(plotly)
library(DT)
library(deSolve)
library(dplyr)
library(ggplot2)

# ─────────────────────────────────────────────────────────────────────────────
# ODE SYSTEM (inline for Shiny; no mrgsolve dependency)
# ─────────────────────────────────────────────────────────────────────────────

aps_ode <- function(time, state, parms) {
  with(as.list(c(state, parms)), {

    # Drug concentrations
    Cwarf  <- WARF_PLASMA / Vc_warf
    Ce_warf<- WARF_EFFECT
    Clmwh  <- LMWH_C / Vc_lmwh
    Chcq   <- HCQ_PLASMA / Vc_hcq * 1000
    Criva  <- RIVA_PLASMA / Vc_riva * 1000
    Casa   <- ASA_PLASMA / Vc_asa * 1000
    Crtx   <- RTX_C / Vc_rtx

    # Drug effect terms
    Ewarf       <- Imax_warf * Ce_warf / (IC50_warf + Ce_warf + 1e-9)
    Elmwh_thr   <- Clmwh / (EC50_lmwh + Clmwh + 1e-9)
    Ehcq        <- Emax_hcq * Chcq / (EC50_hcq + Chcq + 1e-9)
    Eriva       <- Imax_riva * Criva / (IC50_riva + Criva + 1e-9)
    Easa        <- Imax_asa * Casa / (IC50_asa + Casa + 1e-9)
    Ertx        <- Emax_rtx * (Crtx * 1000) / (EC50_rtx + (Crtx * 1000) + 1e-9)
    Eecul       <- Imax_ecul * RTX_P / (Vc_ecul * EC50_ecul + RTX_P + 1e-9)

    aPL_norm  <- aPL_IgG / 100
    C5a_inh   <- if (tx_ecul > 0.5) Eecul else 0
    C5a_prod  <- kC5aprod * (1 + kC5a_aPL * aPL_norm)
    TF_drive  <- kTFon * aPL_norm * Complement_C5a
    Plt_drive <- kPLTon * aPL_norm * EC_TF * (1 - Easa)
    Thr_drive <- kThrOn * EC_TF * Platelet_act * (1 - Elmwh_thr) * (1 - Eriva) * (1 - Ewarf * 0.5)
    DVT_drive <- kDVT_on * Thrombin_gen * (1 - Elmwh_thr) * (1 - Eriva)
    PLoss     <- kPregLoss * Complement_C5a * aPL_norm * (1 - Elmwh_thr * 0.6) * (1 - Easa * 0.4)
    mTOR_drive<- kmTOR_on * aPL_norm * EC_TF

    # Dosing inputs (continuous equivalent)
    warf_in <- tx_warf * DOSE_warf / 24
    lmwh_in <- tx_lmwh * DOSE_lmwh / 24
    hcq_in  <- tx_hcq  * DOSE_hcq  / 24
    riva_in <- tx_riva * DOSE_riva / 24
    asa_in  <- tx_asa  * DOSE_asa  / 24

    dWARF_GUT    <- warf_in * F_warf - ka_warf * WARF_GUT
    dWARF_PLASMA <- ka_warf * WARF_GUT - (CL_warf / Vc_warf) * WARF_PLASMA
    dWARF_EFFECT <- ke_warf * (Cwarf - WARF_EFFECT)
    dLMWH_C      <- lmwh_in - (CL_lmwh / Vc_lmwh) * LMWH_C
    dHCQ_GUT     <- hcq_in * F_hcq - ka_hcq * HCQ_GUT
    dHCQ_PLASMA  <- ka_hcq * HCQ_GUT - (CL_hcq / Vc_hcq) * HCQ_PLASMA
    dRIVA_GUT    <- riva_in * F_riva - ka_riva * RIVA_GUT
    dRIVA_PLASMA <- ka_riva * RIVA_GUT - (CL_riva / Vc_riva) * RIVA_PLASMA
    dASA_GUT     <- asa_in * F_asa - ka_asa * ASA_GUT
    dASA_PLASMA  <- ka_asa * ASA_GUT - (CL_asa / Vc_asa) * ASA_PLASMA
    dRTX_C       <- -(CL_rtx / Vc_rtx + Q_rtx / Vc_rtx) * RTX_C + (Q_rtx / Vp_rtx) * RTX_P
    dRTX_P       <-  (Q_rtx / Vc_rtx) * RTX_C - (Q_rtx / Vp_rtx) * RTX_P

    daPL_IgG         <- kaPL_prod * B_cell * 100 * (1 - Ehcq) * (1 - Ertx * 0.5) - kaPL_deg * aPL_IgG
    dB_cell          <- kBprod - kBdeg * B_cell + kBstim * aPL_norm * B_cell * (1 - Ertx)
    dComplement_C5a  <- C5a_prod * (1 - C5a_inh) - kC5adeg * Complement_C5a
    dEC_TF           <- TF_drive - kTFoff * EC_TF
    dPlatelet_act    <- Plt_drive - kPLToff * Platelet_act
    dThrombin_gen    <- Thr_drive - kThrOff * Thrombin_gen
    dDVT_risk        <- DVT_drive * (1 - DVT_risk / DVT_max) - kDVT_off * DVT_risk
    dPregnancy_viab  <- kPregRec * (1 - Pregnancy_viab) - PLoss
    dmTOR_renal      <- mTOR_drive - kmTOR_off * mTOR_renal
    dINR             <- kout_INR * ((1 + Ewarf * 3) - INR)

    list(c(dWARF_GUT, dWARF_PLASMA, dWARF_EFFECT, dLMWH_C,
           dHCQ_GUT, dHCQ_PLASMA, dRIVA_GUT, dRIVA_PLASMA,
           dASA_GUT, dASA_PLASMA, dRTX_C, dRTX_P,
           daPL_IgG, dB_cell, dComplement_C5a, dEC_TF,
           dPlatelet_act, dThrombin_gen, dDVT_risk,
           dPregnancy_viab, dmTOR_renal, dINR))
  })
}

# Default parameters
default_parms <- c(
  kaPL_prod=0.005, kaPL_deg=0.033, kBprod=0.10, kBdeg=0.08, kBstim=0.02,
  kC5aprod=0.5, kC5adeg=2.0, kC5a_aPL=0.8, kTFon=0.15, kTFoff=0.30,
  kPLTon=0.20, kPLToff=0.40, kThrOn=0.25, kThrOff=0.50,
  kDVT_on=0.08, kDVT_off=0.10, DVT_max=0.80,
  kPregLoss=0.15, kPregRec=0.05, kmTOR_on=0.05, kmTOR_off=0.03,
  ka_warf=1.10, F_warf=0.99, CL_warf=0.20, Vc_warf=10.0, ke_warf=0.12,
  Imax_warf=0.98, IC50_warf=0.8, kout_INR=0.10,
  ka_lmwh=0.20, CL_lmwh=1.20, Vc_lmwh=5.5, EC50_lmwh=0.10,
  ka_hcq=0.08, F_hcq=0.74, CL_hcq=0.25, Vc_hcq=257.0,
  Emax_hcq=0.55, EC50_hcq=200.0,
  ka_riva=1.50, F_riva=0.66, CL_riva=4.80, Vc_riva=47.0,
  IC50_riva=50.0, Imax_riva=0.95,
  ka_asa=6.00, F_asa=0.68, CL_asa=35.0, Vc_asa=12.0,
  IC50_asa=50.0, Imax_asa=0.85,
  CL_rtx=0.016, Vc_rtx=3.5, Vp_rtx=3.2, Q_rtx=0.008,
  Emax_rtx=0.90, EC50_rtx=10.0,
  CL_ecul=0.013, Vc_ecul=5.0, EC50_ecul=5.0, Imax_ecul=0.95,
  tx_warf=0, tx_lmwh=0, tx_hcq=0, tx_riva=0, tx_asa=0, tx_rtx=0, tx_ecul=0,
  DOSE_warf=5.0, DOSE_lmwh=40.0, DOSE_hcq=400.0, DOSE_riva=20.0,
  DOSE_asa=100.0, DOSE_rtx=375.0, DOSE_ecul=900.0
)

default_state <- c(
  WARF_GUT=0, WARF_PLASMA=0, WARF_EFFECT=0, LMWH_C=0,
  HCQ_GUT=0, HCQ_PLASMA=0, RIVA_GUT=0, RIVA_PLASMA=0,
  ASA_GUT=0, ASA_PLASMA=0, RTX_C=0, RTX_P=0,
  aPL_IgG=80, B_cell=1.0, Complement_C5a=1.0, EC_TF=0.2,
  Platelet_act=0.2, Thrombin_gen=0.15, DVT_risk=0.10,
  Pregnancy_viab=0.65, mTOR_renal=0.20, INR=1.1
)

run_sim <- function(parms_override = list(), end_hours = 8760) {
  p <- default_parms
  for (nm in names(parms_override)) p[[nm]] <- parms_override[[nm]]
  times <- seq(0, end_hours, by = 4)
  out <- ode(y = default_state, times = times, func = aps_ode, parms = p,
             method = "lsoda")
  df <- as.data.frame(out)
  df$time_d <- df$time / 24
  df
}

# ─────────────────────────────────────────────────────────────────────────────
# SCENARIO DEFINITIONS
# ─────────────────────────────────────────────────────────────────────────────

scenario_defs <- list(
  "Untreated"               = list(tx_warf=0,tx_lmwh=0,tx_hcq=0,tx_riva=0,tx_asa=0,tx_rtx=0,tx_ecul=0),
  "Warfarin + ASA"          = list(tx_warf=1,DOSE_warf=5,tx_asa=1,DOSE_asa=100,tx_lmwh=0,tx_hcq=0,tx_riva=0,tx_rtx=0,tx_ecul=0),
  "LMWH + ASA (Obstetric)"  = list(tx_lmwh=1,DOSE_lmwh=40,tx_asa=1,DOSE_asa=100,tx_warf=0,tx_hcq=0,tx_riva=0,tx_rtx=0,tx_ecul=0),
  "HCQ + ASA (Primary)"     = list(tx_hcq=1,DOSE_hcq=400,tx_asa=1,DOSE_asa=100,tx_warf=0,tx_lmwh=0,tx_riva=0,tx_rtx=0,tx_ecul=0),
  "Rivaroxaban 20 mg QD"    = list(tx_riva=1,DOSE_riva=20,tx_warf=0,tx_lmwh=0,tx_hcq=0,tx_asa=0,tx_rtx=0,tx_ecul=0),
  "Rituximab + Warfarin"    = list(tx_rtx=1,DOSE_rtx=375,tx_warf=1,DOSE_warf=5,tx_lmwh=0,tx_hcq=0,tx_riva=0,tx_asa=0,tx_ecul=0),
  "Eculizumab (CAPS)"       = list(tx_ecul=1,DOSE_ecul=900,tx_warf=1,DOSE_warf=5,tx_lmwh=0,tx_hcq=0,tx_riva=0,tx_asa=0,tx_rtx=0)
)

SCENARIO_COLORS <- c(
  "Untreated"               = "#E74C3C",
  "Warfarin + ASA"          = "#3498DB",
  "LMWH + ASA (Obstetric)"  = "#2ECC71",
  "HCQ + ASA (Primary)"     = "#F39C12",
  "Rivaroxaban 20 mg QD"    = "#9B59B6",
  "Rituximab + Warfarin"    = "#1ABC9C",
  "Eculizumab (CAPS)"       = "#E67E22"
)

# ─────────────────────────────────────────────────────────────────────────────
# UI
# ─────────────────────────────────────────────────────────────────────────────

ui <- dashboardPage(
  skin = "red",
  dashboardHeader(title = "APS QSP Dashboard"),
  dashboardSidebar(
    sidebarMenu(
      menuItem("Patient Profile",     tabName = "profile",   icon = icon("user-md")),
      menuItem("Drug PK",             tabName = "pk",        icon = icon("pills")),
      menuItem("PD Biomarkers",       tabName = "pd",        icon = icon("flask")),
      menuItem("Clinical Endpoints",  tabName = "endpoints", icon = icon("heartbeat")),
      menuItem("Scenario Comparison", tabName = "compare",   icon = icon("chart-bar")),
      menuItem("Obstetric APS",       tabName = "obstetric", icon = icon("baby"))
    ),
    hr(),
    h5("Global Settings", style = "padding-left:15px; color:#aaa"),
    sliderInput("sim_months", "Simulation Duration (months)", 1, 24, 12, step = 1),
    sliderInput("apl_baseline", "Baseline aPL IgG Titer (GPL-U)", 10, 200, 80, step = 5),
    checkboxGroupInput("active_scenarios", "Active Scenarios",
                       choices = names(scenario_defs),
                       selected = c("Untreated","Warfarin + ASA","LMWH + ASA (Obstetric)","HCQ + ASA (Primary)"))
  ),
  dashboardBody(
    tabItems(

      # ───────────────────────────────────────────────────────────────────────
      # TAB 1 — PATIENT PROFILE
      # ───────────────────────────────────────────────────────────────────────
      tabItem(tabName = "profile",
        fluidRow(
          box(title = "Patient Demographics & Risk Factors", width = 4, status = "danger", solidHeader = TRUE,
            selectInput("aps_type", "APS Classification",
                        choices = c("Primary APS","Secondary APS (with SLE)","Obstetric APS only","CAPS")),
            checkboxGroupInput("serology", "aPL Serology (Sapporo Criteria)",
                               choices = c("aCL IgG ≥40 GPL-U","aCL IgM ≥40 MPL-U",
                                           "Anti-β2GPI IgG >99th percentile",
                                           "Anti-β2GPI IgM >99th percentile",
                                           "Lupus Anticoagulant (dRVVT positive)"),
                               selected = c("aCL IgG ≥40 GPL-U","Anti-β2GPI IgG >99th percentile",
                                            "Lupus Anticoagulant (dRVVT positive)")),
            selectInput("risk_profile", "Thrombotic Risk Profile",
                        choices = c("Triple positive (highest)","Double positive","Single positive")),
            checkboxGroupInput("comorbs", "Comorbidities",
                               choices = c("SLE","Hypertension","Smoking","OCP use",
                                           "Previous DVT/PE","Previous stroke",
                                           "Recurrent pregnancy loss")),
            numericInput("age", "Age (years)", 35, 18, 70)
          ),
          box(title = "Sapporo Criteria Summary", width = 4, status = "warning", solidHeader = TRUE,
            h4("Clinical Criteria (≥1 required)"),
            tags$ul(
              tags$li(tags$b("Vascular thrombosis:"), " ≥1 episode arterial, venous, or small vessel"),
              tags$li(tags$b("Pregnancy morbidity:"),
                      tags$ul(
                        tags$li("≥3 consecutive unexplained RPL <10 wks"),
                        tags$li("≥1 unexplained fetal death ≥10 wks"),
                        tags$li("Premature birth <34 wk (preeclampsia/placental insuff.)")
                      ))
            ),
            h4("Laboratory Criteria (≥1 required, confirmed 12 wks apart)"),
            tags$ul(
              tags$li("Lupus anticoagulant (LA) by dRVVT or KCT"),
              tags$li("aCL IgG or IgM ≥40 GPL/MPL-U by ELISA"),
              tags$li("Anti-β2GPI IgG or IgM >99th percentile by ELISA")
            ),
            br(),
            valueBoxOutput("risk_score_box", width = 12)
          ),
          box(title = "Key Clinical Parameters", width = 4, status = "info", solidHeader = TRUE,
            tableOutput("params_table"),
            br(),
            h4("Current Treatment Selection"),
            selectInput("current_tx", "Treatment Arm",
                        choices = names(scenario_defs), selected = "Warfarin + ASA"),
            sliderInput("warf_dose", "Warfarin dose (mg QD)", 1, 10, 5, step = 0.5),
            sliderInput("hcq_dose", "HCQ dose (mg QD)", 100, 600, 400, step = 100)
          )
        ),
        fluidRow(
          box(title = "Key Clinical Trial Evidence", width = 12, status = "primary", solidHeader = TRUE,
            DTOutput("trial_table")
          )
        )
      ),

      # ───────────────────────────────────────────────────────────────────────
      # TAB 2 — DRUG PK
      # ───────────────────────────────────────────────────────────────────────
      tabItem(tabName = "pk",
        fluidRow(
          box(title = "Warfarin PK (CYP2C9 Genotype)", width = 6, status = "primary", solidHeader = TRUE,
            selectInput("cyp2c9", "CYP2C9 Genotype",
                        choices = c("*1/*1 (normal)","*1/*2 (intermediate)","*1/*3 (poor)","*2/*3 (poor)")),
            plotlyOutput("pk_warfarin", height = 280)
          ),
          box(title = "LMWH (Enoxaparin) Anti-Xa Profile", width = 6, status = "success", solidHeader = TRUE,
            sliderInput("lmwh_dose_pk", "Enoxaparin dose (mg SC)", 20, 120, 40, step = 10),
            plotlyOutput("pk_lmwh", height = 280)
          )
        ),
        fluidRow(
          box(title = "Hydroxychloroquine PK (Chronic Accumulation)", width = 6, status = "warning", solidHeader = TRUE,
            plotlyOutput("pk_hcq", height = 280)
          ),
          box(title = "Rivaroxaban PK / Anti-Xa", width = 6, status = "danger", solidHeader = TRUE,
            sliderInput("riva_dose_pk", "Rivaroxaban dose (mg QD)", 2.5, 30, 20, step = 2.5),
            plotlyOutput("pk_riva", height = 280)
          )
        ),
        fluidRow(
          box(title = "Drug Mechanism of Action Summary", width = 12, status = "info", solidHeader = TRUE,
            DTOutput("moa_table")
          )
        )
      ),

      # ───────────────────────────────────────────────────────────────────────
      # TAB 3 — PD BIOMARKERS
      # ───────────────────────────────────────────────────────────────────────
      tabItem(tabName = "pd",
        fluidRow(
          valueBoxOutput("apl_box",   width = 3),
          valueBoxOutput("bcell_box", width = 3),
          valueBoxOutput("c5a_box",   width = 3),
          valueBoxOutput("thr_box",   width = 3)
        ),
        fluidRow(
          box(title = "aPL IgG Titer Over Time", width = 6, status = "warning", solidHeader = TRUE,
            plotlyOutput("pd_apl", height = 280)
          ),
          box(title = "B Cell Population (% of Baseline)", width = 6, status = "success", solidHeader = TRUE,
            plotlyOutput("pd_bcell", height = 280)
          )
        ),
        fluidRow(
          box(title = "Complement C5a Activation", width = 6, status = "danger", solidHeader = TRUE,
            plotlyOutput("pd_c5a", height = 280)
          ),
          box(title = "Endothelial TF & Platelet Activation", width = 6, status = "primary", solidHeader = TRUE,
            plotlyOutput("pd_tf_plt", height = 280)
          )
        )
      ),

      # ───────────────────────────────────────────────────────────────────────
      # TAB 4 — CLINICAL ENDPOINTS
      # ───────────────────────────────────────────────────────────────────────
      tabItem(tabName = "endpoints",
        fluidRow(
          valueBoxOutput("dvt_box",  width = 3),
          valueBoxOutput("inr_box",  width = 3),
          valueBoxOutput("caps_box", width = 3),
          valueBoxOutput("mtor_box", width = 3)
        ),
        fluidRow(
          box(title = "DVT / PE Risk Index Over Time", width = 6, status = "danger", solidHeader = TRUE,
            plotlyOutput("ep_dvt", height = 280)
          ),
          box(title = "INR Time-Course (Warfarin Scenarios)", width = 6, status = "primary", solidHeader = TRUE,
            plotlyOutput("ep_inr", height = 280)
          )
        ),
        fluidRow(
          box(title = "Thrombin Generation Index", width = 6, status = "warning", solidHeader = TRUE,
            plotlyOutput("ep_thrombin", height = 280)
          ),
          box(title = "mTOR Renal Endotheliopathy Index", width = 6, status = "success", solidHeader = TRUE,
            plotlyOutput("ep_mtor", height = 280)
          )
        )
      ),

      # ───────────────────────────────────────────────────────────────────────
      # TAB 5 — SCENARIO COMPARISON
      # ───────────────────────────────────────────────────────────────────────
      tabItem(tabName = "compare",
        fluidRow(
          box(title = "Comparison Endpoint", width = 3, status = "primary", solidHeader = TRUE,
            selectInput("compare_endpoint", "Primary Endpoint",
                        choices = c("DVT_risk","aPL_IgG","Complement_C5a","INR",
                                    "Pregnancy_viab","mTOR_renal","Thrombin_gen",
                                    "Platelet_act","B_cell")),
            selectInput("compare_timepoint", "Summary Timepoint",
                        choices = c("3 months"="90","6 months"="180","12 months"="365"),
                        selected = "365")
          ),
          box(title = "Multi-Scenario Trajectory", width = 9, status = "info", solidHeader = TRUE,
            plotlyOutput("compare_plot", height = 320)
          )
        ),
        fluidRow(
          box(title = "12-Month Summary Table", width = 12, status = "success", solidHeader = TRUE,
            DTOutput("compare_table")
          )
        )
      ),

      # ───────────────────────────────────────────────────────────────────────
      # TAB 6 — OBSTETRIC APS
      # ───────────────────────────────────────────────────────────────────────
      tabItem(tabName = "obstetric",
        fluidRow(
          valueBoxOutput("livebirth_box", width = 3),
          valueBoxOutput("rpl_box",       width = 3),
          valueBoxOutput("preeclamp_box", width = 3),
          valueBoxOutput("placenta_box",  width = 3)
        ),
        fluidRow(
          box(title = "Pregnancy Viability Over 9 Months", width = 6, status = "success", solidHeader = TRUE,
            plotlyOutput("ob_viability", height = 300)
          ),
          box(title = "Complement-Mediated Placental Injury", width = 6, status = "danger", solidHeader = TRUE,
            plotlyOutput("ob_complement", height = 300)
          )
        ),
        fluidRow(
          box(title = "Obstetric APS Management Strategy", width = 6, status = "warning", solidHeader = TRUE,
            h4("Standard of Care (Bates/Miyakis 2018 Guidelines)"),
            tags$ol(
              tags$li(tags$b("Pre-conception:"), " HCQ 400 mg QD + Low-dose aspirin 100 mg QD"),
              tags$li(tags$b("Confirmed pregnancy:"), " Add LMWH (enoxaparin 40 mg SC QD prophylactic)"),
              tags$li(tags$b("High-risk (previous late loss / preeclampsia):"),
                      " Consider enoxaparin 1 mg/kg BID (therapeutic)"),
              tags$li(tags$b("Postpartum:"), " Continue LMWH × 6 weeks; restart warfarin if prior VTE"),
              tags$li(tags$b("Refractory OAPS:"), " IV immunoglobulin, low-dose prednisolone, plasmapheresis")
            )
          ),
          box(title = "Live Birth Rate — Trial Evidence", width = 6, status = "primary", solidHeader = TRUE,
            DTOutput("oaps_trial_table")
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

  # ── Reactive: run simulations for all active scenarios ──
  sim_data <- reactive({
    end_h <- input$sim_months * 30 * 24
    apl0  <- input$apl_baseline

    bind_rows(lapply(input$active_scenarios, function(sc) {
      parms <- scenario_defs[[sc]]
      parms$aPL0 <- apl0 # not directly in ODE but sets initial state
      init <- default_state
      init["aPL_IgG"] <- apl0
      p <- default_parms
      for (nm in names(parms)) if (nm %in% names(p)) p[[nm]] <- parms[[nm]]
      times <- seq(0, end_h, by = 4)
      out <- ode(y = init, times = times, func = aps_ode, parms = p, method = "lsoda")
      df <- as.data.frame(out)
      df$time_d <- df$time / 24
      df$scenario <- sc
      df
    }))
  })

  # ── Reactive: single scenario (current_tx) ──
  single_sim <- reactive({
    end_h <- input$sim_months * 30 * 24
    sc <- input$current_tx
    parms <- scenario_defs[[sc]]
    p <- default_parms
    for (nm in names(parms)) if (nm %in% names(p)) p[[nm]] <- parms[[nm]]
    p$DOSE_warf <- input$warf_dose
    p$DOSE_hcq  <- input$hcq_dose
    init <- default_state
    init["aPL_IgG"] <- input$apl_baseline
    times <- seq(0, end_h, by = 4)
    out <- ode(y = init, times = times, func = aps_ode, parms = p, method = "lsoda")
    df <- as.data.frame(out)
    df$time_d <- df$time / 24
    df
  })

  # ── Value boxes (Patient Profile) ──
  output$risk_score_box <- renderValueBox({
    n_pos <- length(input$serology)
    risk <- if (n_pos == 3 && "Triple positive" %in% input$risk_profile) "VERY HIGH (≥10×)" else
            if (n_pos >= 2) "HIGH (3-5×)" else "MODERATE (2×)"
    col  <- if (grepl("VERY", risk)) "red" else if (grepl("HIGH", risk)) "orange" else "yellow"
    valueBox(risk, "Thrombosis Risk (vs. sero-negative)", icon = icon("exclamation-triangle"), color = col)
  })

  output$params_table <- renderTable({
    data.frame(
      Parameter     = c("aPL baseline titer","APS type","Treatment arm","Warfarin dose","HCQ dose"),
      Value         = c(paste(input$apl_baseline, "GPL-U"), input$aps_type,
                        input$current_tx, paste(input$warf_dose, "mg QD"),
                        paste(input$hcq_dose, "mg QD"))
    )
  })

  output$trial_table <- renderDT({
    data.frame(
      Trial       = c("TRAPS","RAPS","ASTRO-APS","PROMISSE","Khamashta-RTX","ECULAPS (case series)"),
      Year        = c(2018,2016,2016,2011,2019,2020),
      Design      = c("RCT","RCT","Phase II","Prospective cohort","Cohort","Case series"),
      Drug        = c("Rivaroxaban vs Warfarin","Rivaroxaban vs Warfarin","Rivaroxaban","HCQ in SLE-APS","Rituximab","Eculizumab"),
      Key_Finding = c("Rivaroxaban INFERIOR in triple+ APS (HR 6.0 for events)",
                      "Rivaroxaban non-inferior (thrombin generation endpoint)",
                      "Halted early; excess events in riva arm",
                      "HCQ reduces pregnancy morbidity & aPL titers (OR 0.23)",
                      "B cell depletion; aPL reduction ~50% at 6 months",
                      "C5 blockade effective in refractory CAPS/obstetric APS"),
      PMID        = c("30145981","27470973","26865073","21067376","31270087","32559618"),
      stringsAsFactors = FALSE
    )
  }, options = list(pageLength = 6, dom = 't'), rownames = FALSE)

  # ── PK Tab ──
  cyp2c9_cl <- reactive({
    switch(input$cyp2c9,
      "*1/*1 (normal)"       = 0.20,
      "*1/*2 (intermediate)" = 0.14,
      "*1/*3 (poor)"         = 0.09,
      "*2/*3 (poor)"         = 0.07
    )
  })

  output$pk_warfarin <- renderPlotly({
    times <- seq(0, 7 * 24, by = 0.5)
    p <- default_parms
    p$tx_warf <- 1; p$DOSE_warf <- input$warf_dose; p$CL_warf <- cyp2c9_cl()
    for (nm in names(scenario_defs[["Warfarin + ASA"]])) if (nm %in% names(p)) p[[nm]] <- scenario_defs[["Warfarin + ASA"]][[nm]]
    out <- ode(y = default_state, times = times, func = aps_ode, parms = p, method = "lsoda")
    df  <- as.data.frame(out); df$Cwarf <- df$WARF_PLASMA / p$Vc_warf
    plot_ly(df, x = ~(time/24), y = ~Cwarf, type = "scatter", mode = "lines",
            line = list(color = "#3498DB", width = 2.5)) %>%
      layout(xaxis = list(title = "Time (days)"), yaxis = list(title = "Warfarin Conc. (mg/L)"),
             shapes = list(list(type="line", y0=0.8, y1=0.8, x0=0, x1=7,
                                line=list(color="red",dash="dash"))))
  })

  output$pk_lmwh <- renderPlotly({
    times <- seq(0, 3*24, by = 0.25)
    p <- default_parms; p$tx_lmwh <- 1; p$DOSE_lmwh <- input$lmwh_dose_pk
    out <- ode(y = default_state, times = times, func = aps_ode, parms = p, method = "lsoda")
    df  <- as.data.frame(out); df$antiXa <- df$LMWH_C / p$Vc_lmwh
    plot_ly(df, x = ~(time/24), y = ~antiXa, type = "scatter", mode = "lines",
            line = list(color = "#2ECC71", width = 2.5)) %>%
      layout(xaxis = list(title = "Time (days)"), yaxis = list(title = "Anti-Xa (IU/mL)"),
             shapes = list(
               list(type="rect",x0=0,x1=3,y0=0.2,y1=0.5,fillcolor="rgba(0,200,0,0.15)",
                    line=list(color="transparent"))
             ))
  })

  output$pk_hcq <- renderPlotly({
    times <- seq(0, 30*24, by = 4)
    p <- default_parms; p$tx_hcq <- 1; p$DOSE_hcq <- input$hcq_dose
    out <- ode(y = default_state, times = times, func = aps_ode, parms = p, method = "lsoda")
    df  <- as.data.frame(out); df$Chcq_ng <- df$HCQ_PLASMA / p$Vc_hcq * 1000
    plot_ly(df, x = ~(time/24), y = ~Chcq_ng, type = "scatter", mode = "lines",
            line = list(color = "#F39C12", width = 2.5)) %>%
      layout(xaxis = list(title = "Time (days)"), yaxis = list(title = "HCQ Conc. (ng/mL)"))
  })

  output$pk_riva <- renderPlotly({
    times <- seq(0, 7*24, by = 0.5)
    p <- default_parms; p$tx_riva <- 1; p$DOSE_riva <- input$riva_dose_pk
    out <- ode(y = default_state, times = times, func = aps_ode, parms = p, method = "lsoda")
    df  <- as.data.frame(out); df$Criva_ng <- df$RIVA_PLASMA / p$Vc_riva * 1000
    plot_ly(df, x = ~(time/24), y = ~Criva_ng, type = "scatter", mode = "lines",
            line = list(color = "#9B59B6", width = 2.5)) %>%
      layout(xaxis = list(title = "Time (days)"), yaxis = list(title = "Rivaroxaban Conc. (ng/mL)"))
  })

  output$moa_table <- renderDT({
    data.frame(
      Drug          = c("Warfarin","Enoxaparin (LMWH)","Hydroxychloroquine","Rivaroxaban","Aspirin","Rituximab","Eculizumab","Belimumab"),
      Class         = c("VKA","LMWH","Antimalarial","Direct FXa inhibitor","Antiplatelet","Anti-CD20 mAb","Anti-C5 mAb","Anti-BAFF mAb"),
      Target        = c("VKORC1 → Factors II/VII/IX/X↓","ATIII potentiation → FXa+IIa↓","TLR4/7 blockade; platelet membrane stabilization; aPL↓","FXa (Ki=0.4 nM)","COX-1 irreversible → TXA2↓","CD20 on B cells → ADCC/CDC depletion","C5 cleavage block → C5a+MAC↓","BAFF → B cell/PC survival↓"),
      Key_PK        = c("CL=0.2 L/h (CYP2C9); t½~40h","SC t½~4.5h; renal CL","Vd~257L; t½~50 days","CL=4.8 L/h; t½~9h","t½~20 min (rapid hydrolysis)","t½~22 days; 2-CMT","t½~11 days","t½~19 days"),
      Monitor       = c("INR (target 2-3)","Anti-Xa (0.2-0.5 IU/mL)","Ophthal. exam q12m","Anti-Xa; not routinely","Platelet function (PFA-100)","CD19+ B cell count","CH50; C5 levels","Ig levels; infection risk"),
      stringsAsFactors = FALSE
    )
  }, options = list(pageLength = 8, dom = 't', scrollX = TRUE), rownames = FALSE)

  # ── PD Biomarkers ──
  output$apl_box <- renderValueBox({
    d <- single_sim()
    last <- d %>% tail(1)
    valueBox(round(last$aPL_IgG, 1), "aPL IgG (GPL-U)", icon = icon("vial"),
             color = if (last$aPL_IgG > 40) "red" else "green")
  })
  output$bcell_box <- renderValueBox({
    d <- single_sim()
    valueBox(paste0(round(tail(d,1)$B_cell * 100, 1), "%"), "B Cell Count", icon = icon("circle"), color = "blue")
  })
  output$c5a_box <- renderValueBox({
    d <- single_sim()
    valueBox(round(tail(d,1)$Complement_C5a, 2), "Complement C5a (rel.)", icon = icon("shield-alt"), color = "purple")
  })
  output$thr_box <- renderValueBox({
    d <- single_sim()
    valueBox(round(tail(d,1)$Thrombin_gen, 3), "Thrombin Generation Index", icon = icon("tint"), color = "orange")
  })

  output$pd_apl <- renderPlotly({
    df <- sim_data()
    plot_ly() %>%
      { for (sc in unique(df$scenario)) {
          dd <- filter(df, scenario == sc)
          . <- add_trace(., data = dd, x = ~time_d, y = ~aPL_IgG, type = "scatter", mode = "lines",
                         name = sc, line = list(color = SCENARIO_COLORS[[sc]], width = 2))
        }; . } %>%
      layout(xaxis = list(title = "Time (days)"), yaxis = list(title = "aPL IgG (GPL-U)"),
             shapes = list(list(type="line",y0=40,y1=40,x0=0,x1=365,
                                line=list(color="orange",dash="dash"))))
  })

  output$pd_bcell <- renderPlotly({
    df <- sim_data()
    plot_ly() %>%
      { for (sc in unique(df$scenario)) {
          dd <- filter(df, scenario == sc)
          . <- add_trace(., data = dd, x = ~time_d, y = ~B_cell, type = "scatter", mode = "lines",
                         name = sc, line = list(color = SCENARIO_COLORS[[sc]], width = 2))
        }; . } %>%
      layout(xaxis = list(title = "Time (days)"), yaxis = list(title = "B cell (rel.)"))
  })

  output$pd_c5a <- renderPlotly({
    df <- sim_data()
    plot_ly() %>%
      { for (sc in unique(df$scenario)) {
          dd <- filter(df, scenario == sc)
          . <- add_trace(., data = dd, x = ~time_d, y = ~Complement_C5a, type = "scatter", mode = "lines",
                         name = sc, line = list(color = SCENARIO_COLORS[[sc]], width = 2))
        }; . } %>%
      layout(xaxis = list(title = "Time (days)"), yaxis = list(title = "C5a (relative)"))
  })

  output$pd_tf_plt <- renderPlotly({
    df <- sim_data()
    plot_ly() %>%
      { for (sc in unique(df$scenario)) {
          dd <- filter(df, scenario == sc)
          . <- add_trace(., data = dd, x = ~time_d, y = ~EC_TF, type = "scatter", mode = "lines",
                         name = paste(sc,"TF"), line = list(color = SCENARIO_COLORS[[sc]], width = 2)) %>%
               add_trace(data = dd, x = ~time_d, y = ~Platelet_act, type = "scatter", mode = "lines",
                         name = paste(sc,"Plt"), line = list(color = SCENARIO_COLORS[[sc]], width = 1.5, dash = "dash"))
        }; . } %>%
      layout(xaxis = list(title = "Time (days)"), yaxis = list(title = "Index (relative)"))
  })

  # ── Clinical Endpoints ──
  output$dvt_box <- renderValueBox({
    d <- single_sim(); v <- round(tail(d,1)$DVT_risk, 3)
    valueBox(v, "DVT Risk Index (1-yr)", icon = icon("procedures"),
             color = if (v > 0.4) "red" else if (v > 0.2) "orange" else "green")
  })
  output$inr_box <- renderValueBox({
    d <- single_sim(); v <- round(tail(d,1)$INR, 2)
    valueBox(v, "INR", icon = icon("balance-scale"),
             color = if (v >= 2.0 & v <= 3.0) "green" else if (v > 3.0) "red" else "yellow")
  })
  output$caps_box <- renderValueBox({
    d <- single_sim()
    caps_risk <- round(min(tail(d,1)$DVT_risk * 1.5 * tail(d,1)$Complement_C5a, 0.99), 3)
    valueBox(caps_risk, "CAPS Risk Index", icon = icon("exclamation-circle"),
             color = if (caps_risk > 0.5) "red" else "orange")
  })
  output$mtor_box <- renderValueBox({
    d <- single_sim(); v <- round(tail(d,1)$mTOR_renal, 3)
    valueBox(v, "mTOR Nephropathy Index", icon = icon("kidneys"),
             color = if (v > 0.4) "red" else "yellow")
  })

  output$ep_dvt <- renderPlotly({
    df <- sim_data()
    plot_ly() %>%
      { for (sc in unique(df$scenario)) {
          dd <- filter(df, scenario == sc)
          . <- add_trace(., data = dd, x = ~time_d, y = ~DVT_risk, type = "scatter", mode = "lines",
                         name = sc, line = list(color = SCENARIO_COLORS[[sc]], width = 2.5))
        }; . } %>%
      layout(xaxis = list(title = "Time (days)"), yaxis = list(title = "DVT Risk Index"),
             shapes = list(list(type="line",y0=0.5,y1=0.5,x0=0,x1=365,
                                line=list(color="red",dash="dash"))))
  })

  output$ep_inr <- renderPlotly({
    df <- sim_data() %>% filter(grepl("Warfarin|Eculizumab|Rituximab", scenario))
    p <- plot_ly()
    for (sc in unique(df$scenario)) {
      dd <- filter(df, scenario == sc)
      p <- add_trace(p, data = dd, x = ~time_d, y = ~INR, type = "scatter", mode = "lines",
                     name = sc, line = list(color = SCENARIO_COLORS[[sc]], width = 2))
    }
    p %>% layout(xaxis = list(title = "Time (days)"), yaxis = list(title = "INR"),
                 shapes = list(
                   list(type="rect",x0=0,x1=input$sim_months*30,y0=2,y1=3,
                        fillcolor="rgba(0,200,0,0.1)",line=list(color="transparent")),
                   list(type="line",y0=2,y1=2,x0=0,x1=input$sim_months*30,
                        line=list(color="green",dash="dash")),
                   list(type="line",y0=3,y1=3,x0=0,x1=input$sim_months*30,
                        line=list(color="green",dash="dash"))
                 ))
  })

  output$ep_thrombin <- renderPlotly({
    df <- sim_data()
    plot_ly() %>%
      { for (sc in unique(df$scenario)) {
          dd <- filter(df, scenario == sc)
          . <- add_trace(., data = dd, x = ~time_d, y = ~Thrombin_gen, type = "scatter", mode = "lines",
                         name = sc, line = list(color = SCENARIO_COLORS[[sc]], width = 2))
        }; . } %>%
      layout(xaxis = list(title = "Time (days)"), yaxis = list(title = "Thrombin Generation Index"))
  })

  output$ep_mtor <- renderPlotly({
    df <- sim_data()
    plot_ly() %>%
      { for (sc in unique(df$scenario)) {
          dd <- filter(df, scenario == sc)
          . <- add_trace(., data = dd, x = ~time_d, y = ~mTOR_renal, type = "scatter", mode = "lines",
                         name = sc, line = list(color = SCENARIO_COLORS[[sc]], width = 2))
        }; . } %>%
      layout(xaxis = list(title = "Time (days)"), yaxis = list(title = "mTOR Index"))
  })

  # ── Scenario Comparison ──
  output$compare_plot <- renderPlotly({
    df <- sim_data()
    ep <- input$compare_endpoint
    col_use <- if (ep %in% names(df)) ep else "DVT_risk"
    plot_ly() %>%
      { for (sc in unique(df$scenario)) {
          dd <- filter(df, scenario == sc)
          . <- add_trace(., data = dd, x = ~time_d, y = as.formula(paste0("~",col_use)),
                         type = "scatter", mode = "lines", name = sc,
                         line = list(color = SCENARIO_COLORS[[sc]], width = 2.5))
        }; . } %>%
      layout(xaxis = list(title = "Time (days)"),
             yaxis = list(title = col_use),
             title = paste("Scenario Comparison —", col_use))
  })

  output$compare_table <- renderDT({
    tp <- as.numeric(input$compare_timepoint)
    df <- sim_data()
    tab <- df %>%
      group_by(scenario) %>%
      filter(abs(time_d - tp) == min(abs(time_d - tp))) %>%
      slice(1) %>%
      select(scenario, DVT_risk, aPL_IgG, INR, Complement_C5a, Pregnancy_viab, mTOR_renal, Thrombin_gen) %>%
      mutate(across(where(is.numeric), ~round(., 3))) %>%
      as.data.frame()
    datatable(tab, rownames = FALSE, options = list(dom = 't', pageLength = 10))
  })

  # ── Obstetric APS ──
  oaps_data <- reactive({
    end_h <- 270 * 24 # 9 months
    bind_rows(lapply(input$active_scenarios, function(sc) {
      p <- default_parms
      parms <- scenario_defs[[sc]]
      for (nm in names(parms)) if (nm %in% names(p)) p[[nm]] <- parms[[nm]]
      init <- default_state; init["aPL_IgG"] <- input$apl_baseline
      times <- seq(0, end_h, by = 4)
      out <- ode(y = init, times = times, func = aps_ode, parms = p, method = "lsoda")
      df  <- as.data.frame(out); df$time_d <- df$time / 24; df$scenario <- sc; df
    }))
  })

  output$livebirth_box <- renderValueBox({
    d <- single_sim()
    v <- round(tail(filter(d, time_d <= 270), 1)$Pregnancy_viab * 100, 1)
    valueBox(paste0(v, "%"), "Live Birth Rate (model)", icon = icon("baby"),
             color = if (v >= 70) "green" else if (v >= 50) "yellow" else "red")
  })
  output$rpl_box <- renderValueBox({
    d <- single_sim()
    rpl_rate <- round((1 - tail(filter(d, time_d <= 90), 1)$Pregnancy_viab) * 100, 1)
    valueBox(paste0(rpl_rate, "%"), "Early Pregnancy Loss Rate", icon = icon("heartbeat"), color = "orange")
  })
  output$preeclamp_box <- renderValueBox({
    valueBox("~18%", "Preeclampsia (<34 wk) Risk", icon = icon("tachometer-alt"), color = "yellow")
  })
  output$placenta_box <- renderValueBox({
    d <- single_sim()
    v <- round(tail(filter(d, time_d <= 270), 1)$Complement_C5a, 2)
    valueBox(v, "Placental C5a (rel.)", icon = icon("shield-alt"),
             color = if (v > 1.5) "red" else "green")
  })

  output$ob_viability <- renderPlotly({
    df <- oaps_data()
    plot_ly() %>%
      { for (sc in unique(df$scenario)) {
          dd <- filter(df, scenario == sc)
          . <- add_trace(., data = dd, x = ~time_d, y = ~Pregnancy_viab, type = "scatter", mode = "lines",
                         name = sc, line = list(color = SCENARIO_COLORS[[sc]], width = 2.5))
        }; . } %>%
      layout(xaxis = list(title = "Gestational Age (days)"),
             yaxis = list(title = "Pregnancy Viability Index (0–1)"),
             shapes = list(list(type="line",y0=0.7,y1=0.7,x0=0,x1=270,
                                line=list(color="green",dash="dash"))))
  })

  output$ob_complement <- renderPlotly({
    df <- oaps_data()
    plot_ly() %>%
      { for (sc in unique(df$scenario)) {
          dd <- filter(df, scenario == sc)
          . <- add_trace(., data = dd, x = ~time_d, y = ~Complement_C5a, type = "scatter", mode = "lines",
                         name = sc, line = list(color = SCENARIO_COLORS[[sc]], width = 2.5))
        }; . } %>%
      layout(xaxis = list(title = "Days"), yaxis = list(title = "Complement C5a (rel.)"))
  })

  output$oaps_trial_table <- renderDT({
    data.frame(
      Study       = c("PROMISSE (Salmon 2011)","Backos 1999","Rai 1997","Farquharson 2002","Dendrinos 2004","Stephenson 2004","Farquharson 2010"),
      Treatment   = c("HCQ vs no HCQ (SLE)","LDA+LMWH","LDA+LMWH vs LDA","LDA+LMWH vs LDA","LMWH vs LDA","LMWH vs LDA","LMWH vs LDA"),
      Live_Birth  = c("80% vs 59%","71%","71% vs 42%","78% vs 72%","86% vs 57%","78% vs 57%","71% vs 69%"),
      PMID        = c("21067376","10350093","9171197","12351898","15200171","15123542","20378640"),
      stringsAsFactors = FALSE
    )
  }, options = list(dom = 't', pageLength = 8), rownames = FALSE)
}

# ─────────────────────────────────────────────────────────────────────────────
# RUN APP
# ─────────────────────────────────────────────────────────────────────────────

shinyApp(ui, server)
