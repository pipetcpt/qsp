# ============================================================================
# PAD QSP Shiny Dashboard
# Peripheral Arterial Disease — Interactive PK/PD Simulator
#
# Tabs:
#   1. Patient Profile & Disease Overview
#   2. Drug PK & Target Engagement
#   3. PD Biomarkers (Platelet / Coagulation / Lipids / Inflammation)
#   4. Clinical Endpoints (ABI · Walking Distance · CLI Risk)
#   5. Scenario Comparison (7 treatment regimens head-to-head)
#   6. Risk Stratification & Biomarker Panel
# ============================================================================

library(shiny)
library(shinydashboard)
library(deSolve)
library(ggplot2)
library(dplyr)
library(tidyr)
library(DT)
library(plotly)
library(scales)

# ─────────────────────────────────────────────────────────────────────────────
# ODE SYSTEM  (simplified, self-contained; no mrgsolve dependency for Shiny)
# ─────────────────────────────────────────────────────────────────────────────
pad_ode <- function(time, state, parms) {
  with(as.list(c(state, parms)), {

    # ── Clopidogrel PK ──────────────────────────────────────────────────────
    dA_clopi  <- -ka_clopi * A_clopi
    dC_clopi  <- ka_clopi * A_clopi / Vc_clopi -
                  (CL_clopi / Vc_clopi) * C_clopi - kact * C_clopi
    dC_am     <- kact * C_clopi * Vc_clopi / Vc_am -
                  (CL_am / Vc_am) * C_am

    # ── Aspirin PK ──────────────────────────────────────────────────────────
    dA_asp    <- -ka_asp * A_asp
    dC_asp    <- ka_asp * A_asp / Vc_asp - (CL_asp / Vc_asp) * C_asp

    # ── Ticagrelor PK ───────────────────────────────────────────────────────
    dA_tica   <- -ka_tica * A_tica
    dC_tica   <- ka_tica * A_tica / Vc_tica - (CL_tica / Vc_tica) * C_tica

    # ── Rivaroxaban PK ──────────────────────────────────────────────────────
    dA_riva   <- -ka_riva * A_riva
    dC_riva   <- ka_riva * A_riva / Vc_riva - (CL_riva / Vc_riva) * C_riva

    # ── Cilostazol PK ───────────────────────────────────────────────────────
    dA_cilo   <- -ka_cilo * A_cilo
    dC_cilo   <- ka_cilo * A_cilo / Vc_cilo - (CL_cilo / Vc_cilo) * C_cilo

    # ── Atorvastatin PK ─────────────────────────────────────────────────────
    dA_atst   <- -ka_atst * A_atst
    dC_atst   <- ka_atst * A_atst / Vc_atst - (CL_atst / Vc_atst) * C_atst

    # ── Drug Effects ────────────────────────────────────────────────────────
    Inh_AM    <- Emax_P2Y12 * C_am   / (EC50_P2Y12 + C_am   + 1e-12)
    Inh_tica  <- Emax_tica  * C_tica / (EC50_tica  + C_tica + 1e-12)
    Inh_P2Y12 <- max(Inh_AM, Inh_tica)
    Inh_COX1  <- Emax_COX1  * C_asp  / (EC50_COX1  + C_asp  + 1e-12)
    Inh_FXa   <- Emax_FXa   * C_riva / (EC50_FXa   + C_riva + 1e-12)
    Eff_PDE3  <- Emax_PDE3  * C_cilo / (EC50_PDE3  + C_cilo + 1e-12)
    Eff_HMG   <- Emax_HMG   * C_atst / (EC50_HMG   + C_atst + 1e-12)

    # ── PD: Platelet aggregation ─────────────────────────────────────────────
    Plt_inhib  <- 1 - (1 - Inh_P2Y12) * (1 - 0.65*Inh_COX1) * (1 - 0.35*Eff_PDE3)
    Plt_target <- Plt_agg0 * (1 - Plt_inhib)
    dPlt_agg   <- kplt_rec * (Plt_target - Plt_agg)

    # ── PD: Thrombin index ───────────────────────────────────────────────────
    Thrombin_target <- Thrombin0 * (1 - 0.85*Inh_FXa) * (1 - 0.20*Inh_COX1)
    dThrombin       <- kthrombin * (Thrombin_target - Thrombin)

    # ── PD: LDL-C ────────────────────────────────────────────────────────────
    LDL_cl_total <- k_ldl_cl * (1 + 3.0*Eff_HMG) * LDL_C
    dLDL_C       <- k_ldl_prod * LDL_base - LDL_cl_total

    # ── PD: Plaque volume ────────────────────────────────────────────────────
    plaque_driver <- (LDL_C / 130) * (1 - 0.25*Eff_HMG) * (1 - 0.10*Inh_P2Y12)
    dPlaqueVol    <- k_plaque * PlaqueVol * plaque_driver * 100

    # ── PD: ABI ──────────────────────────────────────────────────────────────
    ABI_prog  <- k_abi_prog * (PlaqueVol / 40) * ABI
    ABI_ben   <- 0.00003 * Eff_PDE3 * max(0, 0.90 - ABI)
    dABI      <- -ABI_prog + ABI_ben

    # ── PD: Collateral vessels ───────────────────────────────────────────────
    isch_signal <- max(0, 1 - ABI / 0.9)
    dCollat     <- k_coll * isch_signal * (100 - Collat) * (1 + 0.30*Eff_PDE3)

    # ── PD: Walking distance ─────────────────────────────────────────────────
    WalkD_pot  <- WalkD0 * (1 + Emax_walk*Eff_PDE3) *
                  (1 + 0.30*(Collat - 20)/80) * (ABI / ABI0)
    dWalkDist  <- k_walk * (WalkD_pot - WalkDist)

    # ── PD: Endothelial function ─────────────────────────────────────────────
    EF_target  <- EF0 * (1 + 0.50*Eff_HMG) * (0.70 + 0.30*ABI/ABI0)
    dEF_idx    <- k_EF * (EF_target - EF_idx)

    # ── PD: hs-CRP ───────────────────────────────────────────────────────────
    CRP_driver <- (LDL_C / 130) * (PlaqueVol / 40) * (1 - 0.35*Eff_HMG)
    dhsCRP     <- k_CRP_prod * CRP_driver * CRP0 - k_CRP_cl * hsCRP

    list(c(dA_clopi, dC_clopi, dC_am,
           dA_asp, dC_asp,
           dA_tica, dC_tica,
           dA_riva, dC_riva,
           dA_cilo, dC_cilo,
           dA_atst, dC_atst,
           dPlt_agg, dThrombin, dLDL_C,
           dABI, dWalkDist, dCollat, dEF_idx, dPlaqueVol, dhsCRP))
  })
}

# Default parameters
default_parms <- list(
  # Clopidogrel
  ka_clopi=0.8, CL_clopi=1400, Vc_clopi=80, kact=0.055,
  CL_am=8, Vc_am=5, EC50_P2Y12=0.006, Emax_P2Y12=0.88,
  # Aspirin
  ka_asp=6.0, CL_asp=55, Vc_asp=12, EC50_COX1=0.10, Emax_COX1=0.99,
  # Ticagrelor
  ka_tica=1.4, CL_tica=22, Vc_tica=88, EC50_tica=0.12, Emax_tica=0.90,
  # Rivaroxaban
  ka_riva=1.5, CL_riva=4.8, Vc_riva=47, EC50_FXa=0.05, Emax_FXa=0.95,
  # Cilostazol
  ka_cilo=0.7, CL_cilo=12, Vc_cilo=115, EC50_PDE3=0.35, Emax_PDE3=0.78,
  Emax_walk=0.82,
  # Atorvastatin
  ka_atst=1.2, CL_atst=28, Vc_atst=340, EC50_HMG=0.002, Emax_HMG=0.55,
  # PD
  kplt_rec=0.005, kthrombin=0.05, Thrombin0=100, Plt_agg0=80,
  k_plaque=0.000020, k_ldl_prod=0.0024, LDL_base=130, k_ldl_cl=0.0000185,
  k_abi_prog=0.00000167, ABI0=0.70, k_walk=0.0002, WalkD0=120,
  k_coll=0.000030, CollD0=20, k_EF=0.0005, EF0=45,
  k_CRP_prod=0.10, k_CRP_cl=0.0289, CRP0=3.5, PlaqueV0=40
)

state_names <- c(
  "A_clopi","C_clopi","C_am","A_asp","C_asp",
  "A_tica","C_tica","A_riva","C_riva",
  "A_cilo","C_cilo","A_atst","C_atst",
  "Plt_agg","Thrombin","LDL_C","ABI","WalkDist","Collat","EF_idx","PlaqueVol","hsCRP"
)

# Build dosing events as data.frame of {time_h, cmt, amount}
build_doses <- function(clopi_dose, asp_dose, tica_dose, riva_dose, cilo_dose,
                        atst_dose, sim_days = 365) {
  times_qd  <- seq(0, (sim_days - 1) * 24, by = 24)
  times_bid <- seq(0, sim_days * 24 - 12,   by = 12)
  doses <- data.frame(time = numeric(), cmt = integer(), amt = numeric())
  if (clopi_dose > 0) doses <- rbind(doses, data.frame(time = times_qd,  cmt = 1L, amt = clopi_dose))
  if (asp_dose   > 0) doses <- rbind(doses, data.frame(time = times_qd,  cmt = 4L, amt = asp_dose))
  if (tica_dose  > 0) doses <- rbind(doses, data.frame(time = times_bid, cmt = 6L, amt = tica_dose))
  if (riva_dose  > 0) doses <- rbind(doses, data.frame(time = times_bid, cmt = 8L, amt = riva_dose))
  if (cilo_dose  > 0) doses <- rbind(doses, data.frame(time = times_bid, cmt = 10L,amt = cilo_dose))
  if (atst_dose  > 0) doses <- rbind(doses, data.frame(time = times_qd,  cmt = 12L,amt = atst_dose))
  doses[order(doses$time), ]
}

# ODE solver with dosing events
run_ode <- function(parms, doses, sim_days = 365) {
  y0 <- c(
    A_clopi=0, C_clopi=0, C_am=0, A_asp=0, C_asp=0,
    A_tica=0, C_tica=0, A_riva=0, C_riva=0,
    A_cilo=0, C_cilo=0, A_atst=0, C_atst=0,
    Plt_agg=parms$Plt_agg0, Thrombin=parms$Thrombin0,
    LDL_C=parms$LDL_base, ABI=parms$ABI0, WalkDist=parms$WalkD0,
    Collat=parms$CollD0, EF_idx=parms$EF0,
    PlaqueVol=parms$PlaqueV0, hsCRP=parms$CRP0
  )
  all_times <- sort(unique(c(
    seq(0, sim_days * 24, by = 6),
    if (nrow(doses) > 0) doses$time else NULL
  )))

  results <- data.frame()
  current_y <- y0
  dose_times <- if (nrow(doses) > 0) doses$time else c()

  for (i in seq_along(all_times)) {
    t_now <- all_times[i]
    if (t_now > 0) {
      t_prev <- all_times[i - 1]
      if (t_prev < t_now) {
        sol <- ode(y = current_y, times = c(t_prev, t_now),
                   func = pad_ode, parms = parms, method = "lsoda")
        current_y <- sol[nrow(sol), -1]
        names(current_y) <- state_names
      }
    }
    if (t_now %in% dose_times) {
      d_rows <- doses[doses$time == t_now, ]
      for (j in seq_len(nrow(d_rows))) {
        current_y[d_rows$cmt[j]] <- current_y[d_rows$cmt[j]] + d_rows$amt[j]
      }
    }
    results <- rbind(results, c(time = t_now, current_y))
  }
  colnames(results) <- c("time", state_names)
  as.data.frame(results) %>%
    mutate(
      time_days  = time / 24,
      ABI_val    = pmax(0.1, pmin(1.5, ABI)),
      AM_ngmL    = C_am   * 1000,
      Riva_ngmL  = C_riva * 1000,
      Cilo_ngmL  = C_cilo * 1000,
      Atst_ngmL  = C_atst * 1000,
      Inh_P2Y12  = pmax(0, pmin(100, Emax_P2Y12 * C_am / (parms$EC50_P2Y12 + C_am + 1e-12) * 100)),
      Inh_FXa    = pmax(0, pmin(100, parms$Emax_FXa * C_riva / (parms$EC50_FXa + C_riva + 1e-12) * 100)),
      Eff_PDE3   = pmax(0, pmin(100, parms$Emax_PDE3 * C_cilo / (parms$EC50_PDE3 + C_cilo + 1e-12) * 100)),
      LDL_red_pct= pmax(0, pmin(100, parms$Emax_HMG * C_atst / (parms$EC50_HMG + C_atst + 1e-12) * 100)),
      MACE_risk  = pmin(100, (Plt_agg/80) * (Thrombin/100) * pmax(0, 1 - ABI) * 60 +
                          (LDL_C/130) * 20 + (hsCRP/3.5) * 20),
      MALE_risk  = pmin(100, (1 - ABI_val) * 80 + (PlaqueVol/100) * 20),
      Rutherford = case_when(
        ABI > 0.9 ~ "0: Normal",
        ABI > 0.7 ~ "I: Mild",
        ABI > 0.5 & WalkDist > 200 ~ "IIa: Moderate claudication",
        ABI > 0.5 ~ "IIb: Severe claudication",
        ABI > 0.4 ~ "III: Ischemic rest pain",
        ABI > 0.2 ~ "IV: Ulceration/Gangrene",
        TRUE ~ "V: Major tissue loss"
      )
    )
}

# ─────────────────────────────────────────────────────────────────────────────
# UI
# ─────────────────────────────────────────────────────────────────────────────
ui <- dashboardPage(
  skin = "blue",
  dashboardHeader(
    title = tags$span(
      style = "font-size:14px; font-weight:bold;",
      "PAD QSP Simulator"
    ),
    titleWidth = 240
  ),

  dashboardSidebar(
    width = 240,
    sidebarMenu(
      menuItem("Patient Profile",     tabName = "patient",   icon = icon("user-circle")),
      menuItem("Drug PK",             tabName = "pk",        icon = icon("pills")),
      menuItem("PD Biomarkers",       tabName = "pd",        icon = icon("heartbeat")),
      menuItem("Clinical Endpoints",  tabName = "endpoints", icon = icon("chart-line")),
      menuItem("Scenario Comparison", tabName = "scenario",  icon = icon("balance-scale")),
      menuItem("Risk Stratification", tabName = "risk",      icon = icon("exclamation-triangle"))
    ),
    hr(),
    tags$div(style = "padding:10px; font-size:11px; color:#AAA;",
      "PAD QSP Model v1.0",
      tags$br(), "Calibrated to CAPRIE, COMPASS,",
      tags$br(), "CASTLE, EUCLID trials"
    )
  ),

  dashboardBody(
    tags$head(tags$style(HTML("
      .content-wrapper { background-color: #F5F7FA; }
      .box { border-radius: 8px; }
      .small-box { border-radius: 8px; }
    "))),

    tabItems(

      # ══════════════════════════════════════════════════════════════════════
      # TAB 1: Patient Profile
      # ══════════════════════════════════════════════════════════════════════
      tabItem(tabName = "patient",
        fluidRow(
          box(title = "Patient Demographics & Baseline", width = 4, status = "primary",
            numericInput("age",      "Age (years):", 65, 40, 90, 5),
            selectInput("sex",       "Sex:",        c("Male", "Female")),
            numericInput("BMI",      "BMI (kg/m²):", 27, 18, 45, 1),
            selectInput("smoking",   "Smoking:",    c("Current", "Former", "Never")),
            selectInput("dm",        "Diabetes:",   c("No", "Yes — well controlled", "Yes — poorly controlled")),
            selectInput("htn",       "Hypertension:", c("Controlled", "Uncontrolled", "No")),
            numericInput("abi_base", "Baseline ABI:", 0.70, 0.10, 0.90, 0.05),
            numericInput("ldl_base", "Baseline LDL-C (mg/dL):", 130, 50, 250, 5)
          ),
          box(title = "Disease Staging", width = 4, status = "warning",
            selectInput("rutherford",
              "Rutherford Classification:",
              choices = c(
                "I — Asymptomatic (ABI 0.7–0.9)",
                "IIa — Mild claudication (ABI 0.5–0.7)",
                "IIb — Moderate-to-severe claudication",
                "III — Ischemic rest pain (ABI < 0.4)",
                "IV — Ischemic ulceration / gangrene"
              ),
              selected = "IIb — Moderate-to-severe claudication"
            ),
            numericInput("walk_base", "Baseline max walking (m):", 120, 10, 1000, 10),
            selectInput("prior_revasc", "Prior Revascularization:",
                        c("None", "Angioplasty", "Bypass surgery")),
            numericInput("crp_base",    "Baseline hs-CRP (mg/L):", 3.5, 0.1, 30, 0.5),
            hr(),
            h5("Cardiovascular History"),
            checkboxInput("hx_mi",     "Prior MI", FALSE),
            checkboxInput("hx_stroke", "Prior Stroke", FALSE),
            checkboxInput("hx_cad",    "Coronary Artery Disease", FALSE)
          ),
          box(title = "Treatment Selection", width = 4, status = "success",
            h5(icon("prescription"), " Antiplatelet Agent"),
            radioButtons("antiplatelet",
              label = NULL,
              choices = c(
                "None",
                "Aspirin 100 mg QD",
                "Clopidogrel 75 mg QD",
                "Ticagrelor 90 mg BID",
                "Aspirin + Clopidogrel (DAPT)"
              ),
              selected = "Aspirin 100 mg QD"
            ),
            hr(),
            h5(icon("syringe"), " Anticoagulant"),
            checkboxInput("use_riva",  "Rivaroxaban 2.5 mg BID (COMPASS)", FALSE),
            hr(),
            h5(icon("walking"), " Vasodilator"),
            checkboxInput("use_cilo",  "Cilostazol 100 mg BID", FALSE),
            hr(),
            h5(icon("flask"), " Statin"),
            checkboxInput("use_atst",  "Atorvastatin 40 mg QD", FALSE),
            hr(),
            numericInput("sim_months", "Simulation (months):", 12, 1, 24, 1)
          )
        ),
        fluidRow(
          valueBoxOutput("vbox_abi"),
          valueBoxOutput("vbox_walk"),
          valueBoxOutput("vbox_risk")
        ),
        fluidRow(
          box(title = "Disease Mechanism Overview", width = 12, status = "info",
            fluidRow(
              column(4,
                tags$h5("Pathophysiology"),
                tags$ul(
                  tags$li("Atherosclerosis causes progressive luminal stenosis of lower limb arteries"),
                  tags$li("Plaque rupture or erosion triggers acute thrombosis → acute limb ischemia"),
                  tags$li("Chronic ischemia causes skeletal muscle atrophy and mitochondrial dysfunction"),
                  tags$li("CLI (critical limb ischemia) → rest pain → tissue loss → amputation"),
                  tags$li("Systemic atherosclerosis confers high MACE risk (MI, stroke, CV death)")
                )
              ),
              column(4,
                tags$h5("Key Drug Targets"),
                tags$table(class = "table table-condensed",
                  tags$thead(tags$tr(tags$th("Drug"), tags$th("Target"), tags$th("Benefit"))),
                  tags$tbody(
                    tags$tr(tags$td("Clopidogrel"), tags$td("P2Y12 (ADP)"), tags$td("MALE↓, MACE↓")),
                    tags$tr(tags$td("Aspirin"),     tags$td("COX-1 (TXA2)"), tags$td("MACE↓")),
                    tags$tr(tags$td("Rivaroxaban"), tags$td("FXa (direct)"), tags$td("MACE↓, MALE↓")),
                    tags$tr(tags$td("Cilostazol"),  tags$td("PDE3 (cAMP↑)"), tags$td("Walking↑")),
                    tags$tr(tags$td("Atorvastatin"),tags$td("HMG-CoA"), tags$td("LDL↓, plaque↓"))
                  )
                )
              ),
              column(4,
                tags$h5("Key Clinical Trials"),
                tags$table(class = "table table-condensed",
                  tags$thead(tags$tr(tags$th("Trial"), tags$th("n"), tags$th("Key Finding"))),
                  tags$tbody(
                    tags$tr(tags$td("CAPRIE (1996)"), tags$td("19,185"), tags$td("Clopi superior to ASA: RRR 8.7%")),
                    tags$tr(tags$td("EUCLID (2016)"), tags$td("13,885"), tags$td("Tica = Clopi for MACE in PAD")),
                    tags$tr(tags$td("COMPASS (2018)"), tags$td("7,470"), tags$td("Riva+ASA: MACE↓ 28%, MALE↓ 46%")),
                    tags$tr(tags$td("CASTLE (2008)"), tags$td("1,439"), tags$td("Cilostazol +40% walking distance")),
                    tags$tr(tags$td("CHARISMA (2006)"),tags$td("15,603"),tags$td("DAPT ↑bleeding vs mono in symptomatic PAD"))
                  )
                )
              )
            )
          )
        )
      ),

      # ══════════════════════════════════════════════════════════════════════
      # TAB 2: Drug PK
      # ══════════════════════════════════════════════════════════════════════
      tabItem(tabName = "pk",
        fluidRow(
          box(title = "PK Parameter Adjustments", width = 3, status = "primary",
            selectInput("pk_display_hours", "Display window:",
                        c("48 h (first dose)", "168 h (1 week)", "720 h (1 month)"),
                        selected = "48 h (first dose)"),
            hr(),
            h5("CYP2C19 Status (Clopidogrel)"),
            selectInput("cyp2c19",
              label = NULL,
              choices = c(
                "Normal metabolizer (kact = 0.055/h)",
                "Poor metabolizer   (kact = 0.020/h)",
                "Rapid metabolizer  (kact = 0.080/h)"
              )
            ),
            hr(),
            h5("Renal Function"),
            selectInput("renal_fn",
              label = NULL,
              choices = c(
                "Normal (eGFR ≥ 60)",
                "Moderate CKD (eGFR 30–59)",
                "Severe CKD (eGFR < 30 — riva dose ↓)"
              )
            )
          ),
          box(title = "Plasma Concentration Profiles", width = 9, status = "info",
            plotlyOutput("pk_plot", height = "480px")
          )
        ),
        fluidRow(
          box(title = "Target Engagement Over Time", width = 6, status = "success",
            plotlyOutput("target_engagement_plot", height = "300px")
          ),
          box(title = "Drug Mechanism Summary", width = 6, status = "warning",
            DT::dataTableOutput("moa_table")
          )
        )
      ),

      # ══════════════════════════════════════════════════════════════════════
      # TAB 3: PD Biomarkers
      # ══════════════════════════════════════════════════════════════════════
      tabItem(tabName = "pd",
        fluidRow(
          box(title = "Platelet Aggregation (ADP-induced)", width = 6, status = "danger",
            plotlyOutput("plt_plot", height = "280px")
          ),
          box(title = "Thrombin Generation Index", width = 6, status = "warning",
            plotlyOutput("thrombin_plot", height = "280px")
          )
        ),
        fluidRow(
          box(title = "LDL-C Trajectory", width = 6, status = "info",
            plotlyOutput("ldl_plot", height = "280px")
          ),
          box(title = "hs-CRP (Inflammation)", width = 6, status = "primary",
            plotlyOutput("crp_plot", height = "280px")
          )
        ),
        fluidRow(
          box(title = "Endothelial Function (FMD %)", width = 6, status = "success",
            plotlyOutput("ef_plot", height = "280px")
          ),
          box(title = "Collateral Vessel Index", width = 6, status = "info",
            plotlyOutput("collat_plot", height = "280px")
          )
        )
      ),

      # ══════════════════════════════════════════════════════════════════════
      # TAB 4: Clinical Endpoints
      # ══════════════════════════════════════════════════════════════════════
      tabItem(tabName = "endpoints",
        fluidRow(
          box(title = "Ankle-Brachial Index (ABI) Trajectory", width = 6, status = "danger",
            plotlyOutput("abi_plot", height = "300px")
          ),
          box(title = "Max Walking Distance", width = 6, status = "success",
            plotlyOutput("walk_plot", height = "300px")
          )
        ),
        fluidRow(
          box(title = "MACE Risk Index", width = 6, status = "warning",
            plotlyOutput("mace_plot", height = "300px")
          ),
          box(title = "MALE Risk Index", width = 6, status = "danger",
            plotlyOutput("male_plot", height = "300px")
          )
        ),
        fluidRow(
          box(title = "Rutherford Classification Over Time", width = 12, status = "info",
            plotlyOutput("rutherford_plot", height = "200px")
          )
        )
      ),

      # ══════════════════════════════════════════════════════════════════════
      # TAB 5: Scenario Comparison
      # ══════════════════════════════════════════════════════════════════════
      tabItem(tabName = "scenario",
        fluidRow(
          box(title = "Scenario Definitions", width = 12, status = "info",
            tags$p("All scenarios simulated with identical patient baseline. Duration: 12 months."),
            fluidRow(
              column(6,
                checkboxGroupInput("scen_select",
                  "Select scenarios to compare:",
                  choices = c(
                    "1. No treatment",
                    "2. Aspirin 100 mg QD",
                    "3. Clopidogrel 75 mg QD",
                    "4. DAPT (Clopi + ASA)",
                    "5. COMPASS (Riva 2.5 mg BID + ASA)",
                    "6. Cilostazol 100 mg BID + ASA",
                    "7. Optimal (Clopi+ASA+Riva+Statin)"
                  ),
                  selected = c("2. Aspirin 100 mg QD", "3. Clopidogrel 75 mg QD",
                               "5. COMPASS (Riva 2.5 mg BID + ASA)",
                               "6. Cilostazol 100 mg BID + ASA")
                )
              ),
              column(6,
                selectInput("scen_endpoint",
                  "Primary endpoint to compare:",
                  choices = c("ABI", "Walking distance (m)", "Platelet aggregation (%)",
                              "MACE risk index", "MALE risk index", "LDL-C (mg/dL)", "hs-CRP (mg/L)"),
                  selected = "MACE risk index"
                ),
                actionButton("run_comparison", "Run All Scenarios", class = "btn-primary btn-lg",
                             icon = icon("play"))
              )
            )
          )
        ),
        fluidRow(
          box(title = "Head-to-Head Comparison", width = 8, status = "primary",
            plotlyOutput("comparison_plot", height = "400px")
          ),
          box(title = "Summary at 12 Months", width = 4, status = "success",
            DT::dataTableOutput("summary_table")
          )
        )
      ),

      # ══════════════════════════════════════════════════════════════════════
      # TAB 6: Risk Stratification
      # ══════════════════════════════════════════════════════════════════════
      tabItem(tabName = "risk",
        fluidRow(
          valueBoxOutput("vbox_rutherford"),
          valueBoxOutput("vbox_mace_risk"),
          valueBoxOutput("vbox_male_risk")
        ),
        fluidRow(
          box(title = "Biomarker Dashboard at Simulation End", width = 6, status = "primary",
            plotlyOutput("radar_plot", height = "350px")
          ),
          box(title = "Risk Trend Over Simulation", width = 6, status = "warning",
            plotlyOutput("risk_trend_plot", height = "350px")
          )
        ),
        fluidRow(
          box(title = "Plaque Volume Progression", width = 6, status = "danger",
            plotlyOutput("plaque_plot", height = "280px")
          ),
          box(title = "ABI by Fontaine Stage", width = 6, status = "info",
            DT::dataTableOutput("fontaine_table")
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

  # Reactive: parse treatment selection
  dose_inputs <- reactive({
    clopi <- if (grepl("Clopidogrel|DAPT", input$antiplatelet)) 75 else 0
    asp   <- if (grepl("Aspirin|DAPT", input$antiplatelet)) 100 else 0
    tica  <- if (grepl("Ticagrelor", input$antiplatelet)) 90 else 0
    riva  <- if (isTRUE(input$use_riva)) 2.5 else 0
    cilo  <- if (isTRUE(input$use_cilo)) 100 else 0
    atst  <- if (isTRUE(input$use_atst)) 40  else 0
    list(clopi = clopi, asp = asp, tica = tica,
         riva = riva, cilo = cilo, atst = atst)
  })

  # Reactive: adjust parameters from UI
  parms_reactive <- reactive({
    p <- default_parms
    p$ABI0     <- input$abi_base
    p$WalkD0   <- input$walk_base
    p$LDL_base <- input$ldl_base
    p$CRP0     <- input$crp_base
    p$PlaqueV0 <- 40
    # CYP2C19
    p$kact <- switch(input$cyp2c19,
      "Normal metabolizer (kact = 0.055/h)" = 0.055,
      "Poor metabolizer   (kact = 0.020/h)" = 0.020,
      "Rapid metabolizer  (kact = 0.080/h)" = 0.080
    )
    p
  })

  sim_days_rv <- reactive({ input$sim_months * 30 })

  # Main simulation
  sim_result <- reactive({
    withProgress(message = "Running ODE simulation...", value = 0, {
      d  <- dose_inputs()
      p  <- parms_reactive()
      sd <- sim_days_rv()
      incProgress(0.3)
      doses <- build_doses(d$clopi, d$asp, d$tica, d$riva, d$cilo, d$atst, sd)
      incProgress(0.3)
      result <- run_ode(p, doses, sd)
      incProgress(0.4)
      result
    })
  })

  # ── Value Boxes ──────────────────────────────────────────────────────────
  output$vbox_abi <- renderValueBox({
    df <- sim_result()
    final_abi <- tail(df$ABI_val, 1)
    color <- if (final_abi > 0.9) "green" else if (final_abi > 0.5) "yellow" else "red"
    valueBox(round(final_abi, 3), "ABI at Simulation End", icon = icon("heartbeat"), color = color)
  })

  output$vbox_walk <- renderValueBox({
    df  <- sim_result()
    val <- round(tail(df$WalkDist, 1), 0)
    base_val <- head(df$WalkDist, 1)
    pct_change <- round((val - base_val) / base_val * 100, 1)
    valueBox(paste0(val, " m (", ifelse(pct_change >= 0, "+", ""), pct_change, "%)"),
             "Max Walking Distance", icon = icon("walking"), color = "green")
  })

  output$vbox_risk <- renderValueBox({
    df  <- sim_result()
    val <- round(tail(df$MACE_risk, 1), 1)
    color <- if (val < 20) "green" else if (val < 50) "yellow" else "red"
    valueBox(paste0(val, "/100"), "MACE Risk Index", icon = icon("exclamation-triangle"), color = color)
  })

  # ── PK Plot ──────────────────────────────────────────────────────────────
  output$pk_plot <- renderPlotly({
    df  <- sim_result()
    win <- switch(input$pk_display_hours,
      "48 h (first dose)"  = 48,
      "168 h (1 week)"     = 168,
      "720 h (1 month)"    = 720
    )
    pk_df <- df %>%
      filter(time <= win) %>%
      select(time, AM_ngmL, Riva_ngmL, Cilo_ngmL, Atst_ngmL) %>%
      pivot_longer(-time, names_to = "Drug", values_to = "Conc") %>%
      mutate(Drug = recode(Drug,
        AM_ngmL   = "Clopidogrel AM",
        Riva_ngmL = "Rivaroxaban",
        Cilo_ngmL = "Cilostazol",
        Atst_ngmL = "Atorvastatin"
      ))
    p <- ggplot(pk_df, aes(x = time, y = Conc, color = Drug)) +
      geom_line(size = 0.8) +
      facet_wrap(~Drug, scales = "free_y") +
      labs(x = "Time (h)", y = "Plasma Concentration (ng/mL)",
           title = "Drug PK Profiles") +
      theme_bw() + theme(legend.position = "none")
    ggplotly(p)
  })

  output$target_engagement_plot <- renderPlotly({
    df <- sim_result()
    te_df <- df %>%
      filter(time_days <= sim_days_rv()) %>%
      select(time_days, Inh_P2Y12, Inh_FXa, Eff_PDE3, LDL_red_pct) %>%
      pivot_longer(-time_days, names_to = "Effect", values_to = "Pct") %>%
      mutate(Effect = recode(Effect,
        Inh_P2Y12  = "P2Y12 Inhibition",
        Inh_FXa    = "FXa Inhibition",
        Eff_PDE3   = "PDE3 Inhibition",
        LDL_red_pct= "LDL Reduction"
      ))
    p <- ggplot(te_df, aes(x = time_days, y = Pct, color = Effect)) +
      geom_line(size = 0.8) +
      scale_y_continuous(limits = c(0, 105), labels = percent_format(scale = 1)) +
      labs(x = "Time (days)", y = "Effect (%)", title = "Target Engagement") +
      theme_bw()
    ggplotly(p)
  })

  output$moa_table <- DT::renderDataTable({
    data.frame(
      Drug        = c("Clopidogrel", "Aspirin", "Ticagrelor", "Rivaroxaban", "Cilostazol", "Atorvastatin"),
      Mechanism   = c("CYP2C19 prodrug → AM → P2Y12 block",
                      "Irreversible COX-1 → TXA2↓",
                      "Reversible P2Y12 antagonist (direct)",
                      "Direct FXa inhibitor",
                      "PDE3 inhibition → cAMP↑",
                      "HMG-CoA reductase inhibitor"),
      PD_Effect   = c("Plt aggregation↓ ~50%", "TXB2↓ >90%", "Plt agg↓ ~60%",
                      "Thrombin gen↓ ~80%", "Walking distance↑ ~40%", "LDL-C↓ ~40–55%"),
      Trial       = c("CAPRIE", "CAPRIE", "EUCLID", "COMPASS", "CASTLE", "REACH Registry")
    )
  }, options = list(pageLength = 6, dom = "t"), rownames = FALSE)

  # ── PD Plots ─────────────────────────────────────────────────────────────
  make_pd_plot <- function(df, yvar, ylabel, color, ref_line = NULL, ref_label = "") {
    p <- df %>%
      ggplot(aes_string(x = "time_days", y = yvar)) +
      geom_line(color = color, size = 0.9)
    if (!is.null(ref_line)) {
      p <- p + geom_hline(yintercept = ref_line, linetype = "dashed", color = "gray40") +
               annotate("text", x = min(df$time_days) + 3, y = ref_line * 1.02,
                        label = ref_label, size = 3)
    }
    p + labs(x = "Time (days)", y = ylabel) + theme_bw()
  }

  output$plt_plot <- renderPlotly({
    df <- sim_result()
    ggplotly(make_pd_plot(df, "Plt_agg", "Platelet Aggregation (%)", "steelblue",
                          80, "Baseline (80%)"))
  })
  output$thrombin_plot <- renderPlotly({
    df <- sim_result()
    ggplotly(make_pd_plot(df, "Thrombin", "Thrombin Generation Index (nM·min)",
                          "darkorange", 100, "Baseline (100)"))
  })
  output$ldl_plot <- renderPlotly({
    df <- sim_result()
    ggplotly(make_pd_plot(df, "LDL_C", "LDL-C (mg/dL)", "firebrick",
                          70, "ESC target < 70 mg/dL"))
  })
  output$crp_plot <- renderPlotly({
    df <- sim_result()
    ggplotly(make_pd_plot(df, "hsCRP", "hs-CRP (mg/L)", "purple",
                          1.0, "Low-risk < 1.0 mg/L"))
  })
  output$ef_plot <- renderPlotly({
    df <- sim_result()
    ggplotly(make_pd_plot(df, "EF_idx", "Endothelial FMD (%)", "green4"))
  })
  output$collat_plot <- renderPlotly({
    df <- sim_result()
    ggplotly(make_pd_plot(df, "Collat", "Collateral Vessel Index (0-100)", "teal"))
  })

  # ── Clinical Endpoint Plots ───────────────────────────────────────────────
  output$abi_plot <- renderPlotly({
    df <- sim_result()
    p <- ggplot(df, aes(x = time_days, y = ABI_val)) +
      geom_line(color = "firebrick", size = 1) +
      geom_hline(yintercept = 0.9, linetype = "dashed", color = "gray50") +
      geom_hline(yintercept = 0.4, linetype = "dashed", color = "red") +
      annotate("text", x = 5, y = 0.92, label = "Normal (0.9)", size = 3) +
      annotate("text", x = 5, y = 0.38, label = "CLI threshold (0.4)", color = "red", size = 3) +
      labs(x = "Time (days)", y = "ABI", title = "ABI Trajectory") +
      theme_bw()
    ggplotly(p)
  })
  output$walk_plot <- renderPlotly({
    df <- sim_result()
    p <- ggplot(df, aes(x = time_days, y = WalkDist)) +
      geom_line(color = "darkgreen", size = 1) +
      geom_hline(yintercept = input$walk_base, linetype = "dashed", color = "gray50") +
      labs(x = "Time (days)", y = "Max Walking Distance (m)",
           title = "Walking Distance") + theme_bw()
    ggplotly(p)
  })
  output$mace_plot <- renderPlotly({
    df <- sim_result()
    p <- ggplot(df, aes(x = time_days, y = MACE_risk)) +
      geom_line(color = "darkorange", size = 1) +
      scale_y_continuous(limits = c(0, 100)) +
      labs(x = "Time (days)", y = "MACE Risk Index (0-100)") + theme_bw()
    ggplotly(p)
  })
  output$male_plot <- renderPlotly({
    df <- sim_result()
    p <- ggplot(df, aes(x = time_days, y = MALE_risk)) +
      geom_line(color = "red3", size = 1) +
      scale_y_continuous(limits = c(0, 100)) +
      labs(x = "Time (days)", y = "MALE Risk Index (0-100)") + theme_bw()
    ggplotly(p)
  })
  output$rutherford_plot <- renderPlotly({
    df <- sim_result()
    p <- ggplot(df, aes(x = time_days, y = ABI_val, fill = Rutherford)) +
      geom_area(alpha = 0.6) +
      scale_fill_brewer(palette = "RdYlGn") +
      labs(x = "Time (days)", y = "ABI", fill = "Rutherford Class",
           title = "Rutherford Classification by ABI") + theme_bw()
    ggplotly(p)
  })

  # ── Scenario Comparison ───────────────────────────────────────────────────
  scenario_defs <- list(
    "1. No treatment"                       = list(0, 0, 0, 0, 0, 0),
    "2. Aspirin 100 mg QD"                  = list(0, 100, 0, 0, 0, 0),
    "3. Clopidogrel 75 mg QD"               = list(75, 0, 0, 0, 0, 0),
    "4. DAPT (Clopi + ASA)"                 = list(75, 100, 0, 0, 0, 0),
    "5. COMPASS (Riva 2.5 mg BID + ASA)"    = list(0, 100, 0, 2.5, 0, 0),
    "6. Cilostazol 100 mg BID + ASA"        = list(0, 100, 0, 0, 100, 0),
    "7. Optimal (Clopi+ASA+Riva+Statin)"    = list(75, 100, 0, 2.5, 0, 40)
  )

  scenario_results <- eventReactive(input$run_comparison, {
    req(input$scen_select)
    p  <- parms_reactive()
    sd <- 365
    withProgress(message = "Running all scenarios...", value = 0, {
      imap_dfr(scenario_defs[input$scen_select], function(d, nm) {
        doses  <- build_doses(d[[1]], d[[2]], d[[3]], d[[4]], d[[5]], d[[6]], sd)
        result <- run_ode(p, doses, sd)
        result$scenario <- nm
        incProgress(1 / length(input$scen_select))
        result
      })
    })
  })

  output$comparison_plot <- renderPlotly({
    df <- scenario_results()
    yvar <- switch(input$scen_endpoint,
      "ABI"                     = "ABI_val",
      "Walking distance (m)"    = "WalkDist",
      "Platelet aggregation (%)"= "Plt_agg",
      "MACE risk index"         = "MACE_risk",
      "MALE risk index"         = "MALE_risk",
      "LDL-C (mg/dL)"           = "LDL_C",
      "hs-CRP (mg/L)"           = "hsCRP"
    )
    p <- ggplot(df, aes_string(x = "time_days", y = yvar, color = "scenario")) +
      geom_line(size = 0.9) +
      labs(x = "Time (days)", y = input$scen_endpoint, color = "Scenario",
           title = paste0("Scenario Comparison: ", input$scen_endpoint)) +
      theme_bw() + theme(legend.position = "bottom")
    ggplotly(p) %>% layout(legend = list(orientation = "h", y = -0.2))
  })

  output$summary_table <- DT::renderDataTable({
    df <- scenario_results()
    df %>%
      filter(abs(time - 8760) < 7) %>%
      group_by(scenario) %>%
      slice_tail(n = 1) %>%
      ungroup() %>%
      transmute(
        Scenario   = scenario,
        ABI        = round(ABI_val, 3),
        Walk_m     = round(WalkDist, 0),
        MACE       = round(MACE_risk, 1),
        MALE       = round(MALE_risk, 1),
        LDL        = round(LDL_C, 0)
      )
  }, options = list(pageLength = 8, dom = "t"), rownames = FALSE)

  # ── Risk Stratification ───────────────────────────────────────────────────
  output$vbox_rutherford <- renderValueBox({
    df  <- sim_result()
    val <- tail(df$Rutherford, 1)
    col <- if (grepl("0|I\\b", val)) "green" else if (grepl("II", val)) "yellow" else "red"
    valueBox(val, "Rutherford Class", icon = icon("user-md"), color = col)
  })
  output$vbox_mace_risk <- renderValueBox({
    df  <- sim_result()
    val <- round(tail(df$MACE_risk, 1), 1)
    col <- if (val < 20) "green" else if (val < 50) "yellow" else "red"
    valueBox(paste0(val, "/100"), "MACE Risk Index", icon = icon("heart"), color = col)
  })
  output$vbox_male_risk <- renderValueBox({
    df  <- sim_result()
    val <- round(tail(df$MALE_risk, 1), 1)
    col <- if (val < 20) "green" else if (val < 50) "yellow" else "red"
    valueBox(paste0(val, "/100"), "MALE Risk Index", icon = icon("walking"), color = col)
  })

  output$radar_plot <- renderPlotly({
    df   <- sim_result()
    last <- tail(df, 1)
    cats <- c("ABI (×10)", "Walk/100m", "Plt Agg/80", "LDL/130", "hs-CRP/3.5",
              "Collateral", "EF index")
    vals <- c(
      last$ABI_val * 10,
      last$WalkDist / 100,
      last$Plt_agg / 80,
      last$LDL_C / 130,
      last$hsCRP / 3.5,
      last$Collat / 100,
      last$EF_idx / 100
    )
    plot_ly(type = "scatterpolar", r = c(vals, vals[1]),
            theta = c(cats, cats[1]), fill = "toself",
            fillcolor = "rgba(31, 119, 180, 0.3)",
            line = list(color = "rgb(31, 119, 180)")) %>%
      layout(polar = list(radialaxis = list(range = c(0, 2))),
             title = "Biomarker Radar (relative to baseline)")
  })

  output$risk_trend_plot <- renderPlotly({
    df <- sim_result()
    p  <- ggplot(df, aes(x = time_days)) +
      geom_line(aes(y = MACE_risk, color = "MACE Risk"), size = 0.9) +
      geom_line(aes(y = MALE_risk, color = "MALE Risk"), size = 0.9) +
      scale_color_manual(values = c("MACE Risk" = "darkorange", "MALE Risk" = "red3")) +
      scale_y_continuous(limits = c(0, 100)) +
      labs(x = "Time (days)", y = "Risk Index (0-100)", color = "") + theme_bw()
    ggplotly(p)
  })

  output$plaque_plot <- renderPlotly({
    df <- sim_result()
    p  <- ggplot(df, aes(x = time_days, y = PlaqueVol)) +
      geom_line(color = "saddlebrown", size = 1) +
      labs(x = "Time (days)", y = "Plaque Burden Index (0-100)",
           title = "Atherosclerotic Plaque Progression") + theme_bw()
    ggplotly(p)
  })

  output$fontaine_table <- DT::renderDataTable({
    data.frame(
      Stage  = c("I", "IIa", "IIb", "III", "IV"),
      ABI    = c("> 0.9", "0.7–0.9", "0.5–0.7", "0.3–0.5", "< 0.3"),
      Symptoms = c("Asymptomatic", "Mild claudication (> 200m)",
                   "Mod-severe claudication (< 200m)", "Ischemic rest pain",
                   "Tissue loss / gangrene"),
      MACErate = c("3%/yr", "4%/yr", "5%/yr", "8%/yr", "15%/yr"),
      Amputation = c("< 1%/yr", "< 1%/yr", "2%/yr", "10%/yr", "30%/yr")
    )
  }, options = list(pageLength = 5, dom = "t"), rownames = FALSE)
}

# ─────────────────────────────────────────────────────────────────────────────
# RUN
# ─────────────────────────────────────────────────────────────────────────────
shinyApp(ui, server)
