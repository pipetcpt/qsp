## ============================================================
## ITP QSP Shiny App — Interactive Dashboard
## Immune Thrombocytopenic Purpura
## ============================================================
##
## Tabs:
##   1. Patient Profile & Disease Overview
##   2. Pharmacokinetics (drug concentration-time profiles)
##   3. Platelet & Immune Biomarkers (PD core)
##   4. Clinical Endpoints (response rates, bleeding risk)
##   5. Treatment Scenario Comparison (6 scenarios)
##   6. Mechanistic Panel (Ab, Treg/Th17, TPO, MK, Macrophage)
## ============================================================

library(shiny)
library(shinydashboard)
library(dplyr)
library(ggplot2)
library(tidyr)
library(plotly)
library(DT)

## ── Embedded ODE solver (deSolve, no mrgsolve dependency needed for app) ──
library(deSolve)

# ─── ODE function ────────────────────────────────────────────────────────────
itp_ode <- function(t, state, parms) {
  with(as.list(c(state, parms)), {

    # Drug concentrations
    PRED_conc  <- PRED_c  / V_pred
    IVIG_conc  <- IVIG_c  / V_IVIG
    ROMI_conc  <- ROMI_c  / V_romi
    RTX_conc   <- RTX_c   / V1_rtx
    R788_conc  <- R788_c  / V_R788
    EFGAR_conc <- EFGAR_c / V_efgar
    ELTP_conc  <- ELTP_c  / V_eltp

    # PD effects
    E_ster_Ab   <- Emax_ster_Ab   * PRED_conc / (EC50_ster_Ab   + PRED_conc + 1e-9)
    E_ster_Treg <- Emax_ster_Treg * PRED_conc / (EC50_ster_Treg + PRED_conc + 1e-9)
    E_ster_Mac  <- Emax_ster_Mac  * PRED_conc / (EC50_ster_Mac  + PRED_conc + 1e-9)
    E_ivig_Mac  <- Emax_ivig_Mac  * IVIG_conc / (EC50_ivig_Mac  + IVIG_conc + 1e-9)
    E_ivig_FcRn <- Emax_ivig_FcRn * IVIG_conc / (EC50_ivig_FcRn + IVIG_conc + 1e-9)
    E_efgar_FcRn<- Emax_FcRn      * EFGAR_conc/ (EC50_FcRn      + EFGAR_conc+ 1e-9)
    E_romi      <- Emax_romi * ROMI_conc^n_romi /
                   (EC50_romi^n_romi + ROMI_conc^n_romi + 1e-9)
    E_eltp      <- Emax_eltp * ELTP_conc^n_eltp /
                   (EC50_eltp^n_eltp + ELTP_conc^n_eltp + 1e-9)
    E_TPORA     <- max(E_romi, E_eltp)
    E_rtx_Bc    <- Emax_rtx * RTX_conc / (EC50_rtx + RTX_conc + 1e-9)
    E_syk       <- R788_conc^n_syk / (IC50_syk^n_syk + R788_conc^n_syk + 1e-9)

    FcRn_factor <- 1.0 + E_ivig_FcRn + E_efgar_FcRn
    f_TPO       <- TPO^TPO_n / (TPO_EC50^TPO_n + TPO^TPO_n + 1e-9)
    Mac_eff     <- Mac * (1.0 - E_ster_Mac) * (1.0 - E_ivig_Mac) * (1.0 - E_syk)
    Mac_eff     <- max(Mac_eff, 0)
    kdes        <- kdes_spleen * Ab * Mac_eff
    Ab_BM       <- f_Ab_MK * Ab
    MK_AB_inh   <- Ab_BM / (km_MK + Ab_BM + 1e-9)
    Th17_stim   <- stim_Th17Bc * Th17
    Treg_inh_Bc <- inh_Treg * Treg

    # Dosing inputs (handled via event list → state resets)
    dPLT    <- kMK_rel * MK * PLT0 - kplt_kel * PLT - kdes * PLT -
               kplt_sp_in * PLT + kplt_sp_out * PLT_SP - kdes_liver * Ab * PLT
    dPLT_SP <- kplt_sp_in * PLT - kplt_sp_out * PLT_SP - kdes * PLT_SP
    dTPO    <- kTPO_prod * (PLT0 / (PLT + 0.1)) - kTPO_clear * TPO - kTPO_plt * PLT * TPO
    dMKP    <- kMKP_prod * (1.0 + E_TPORA) * f_TPO - kMKP_mat * MKP
    dMK     <- kMKP_mat * MKP - kMK_die * MK - kMK_rel * MK * (1.0 - MK_AB_inh)
    dAb     <- kAb_prod * Bc * (1.0 - E_ster_Ab) - kAb_clear * FcRn_factor * Ab
    dBc     <- ITP_sev * kBc_stim * (1.0 + Th17_stim - Treg_inh_Bc) -
               kBc_die * Bc * (1.0 + E_rtx_Bc) +
               k_Bc_rec * (Bc0 - Bc) * (1.0 - E_rtx_Bc)
    dTreg   <- kTreg_prod * (1.0 + E_ster_Treg) * (1.0 - inh_Treg * Th17 / (Th170 + 1e-9)) -
               kTreg_die * Treg
    dTh17   <- kTh17_prod * (1.0 - inh_Treg * Treg) - kTh17_die * Th17
    dMac    <- kMac_act * Ab * (PLT + PLT_SP) / (PLT0 + 66.0) - kMac_die * Mac

    # Drug PK
    dPRED_c   <- -CL_pred * PRED_c / V_pred
    dIVIG_c   <- -kIVIG_el * IVIG_c
    dROMI_sc  <- -ka_romi * ROMI_sc
    dROMI_c   <-  ka_romi * ROMI_sc - CL_romi * ROMI_c / V_romi
    dRTX_c    <- -(CL_rtx + Q_rtx) * RTX_c / V1_rtx + Q_rtx * RTX_p / V2_rtx
    dRTX_p    <-  Q_rtx * RTX_c / V1_rtx - Q_rtx * RTX_p / V2_rtx
    dFOSTA_gut<- -ka_fosta * FOSTA_gut
    dR788_c   <-  ka_fosta * FOSTA_gut * 0.8 - CL_R788 * R788_c / V_R788
    dEFGAR_c  <- -kEl_efgar * EFGAR_c
    dELTP_c   <- -CL_eltp * ELTP_c / V_eltp

    list(c(dPLT, dPLT_SP, dTPO, dMKP, dMK, dAb, dBc, dTreg, dTh17, dMac,
           dPRED_c, dIVIG_c, dROMI_sc, dROMI_c, dRTX_c, dRTX_p,
           dFOSTA_gut, dR788_c, dEFGAR_c, dELTP_c))
  })
}

# ─── Default parameters ──────────────────────────────────────────────────────
default_parms <- c(
  PLT0=200, kplt_kel=0.116, kplt_sp_in=0.5, kplt_sp_out=1.0,
  kdes_spleen=3.0, kdes_liver=0.8, f_Ab_MK=0.15, km_MK=0.5,
  MKP0=1.0, MK0=1.0, kMKP_prod=0.2, kMKP_mat=0.2,
  kMK_rel=0.116, kMK_die=0.1,
  TPO_EC50=0.5, TPO_n=1.5, kTPO_prod=0.2, kTPO_clear=0.2, kTPO_plt=0.05, TPO0=1.0,
  Ab0=0.2, kAb_prod=0.04, kAb_clear=0.033,
  Bc0=1.0, kBc_stim=0.02, kBc_die=0.02,
  Treg0=1.0, kTreg_prod=0.1, kTreg_die=0.1,
  Th170=1.0, kTh17_prod=0.1, kTh17_die=0.1,
  Mac0=1.0, kMac_act=0.4, kMac_die=0.4,
  inh_Treg=0.3, stim_Th17Bc=0.15,
  ITP_sev=4.0,
  ka_pred=2.0, CL_pred=480, V_pred=40,
  kIVIG_el=0.033, V_IVIG=5,
  ka_romi=0.7, CL_romi=0.163, V_romi=1.4,
  Emax_romi=3.0, EC50_romi=2.0, n_romi=1.5,
  ka_eltp=19.2, CL_eltp=9.12, V_eltp=40,
  Emax_eltp=2.5, EC50_eltp=600, n_eltp=1.3,
  CL_rtx=0.008, V1_rtx=3.1, V2_rtx=4.1, Q_rtx=0.009,
  Emax_rtx=0.95, EC50_rtx=5.0, k_Bc_rec=0.007,
  ka_fosta=12.0, ka_R788=2.4, CL_R788=360, V_R788=250,
  IC50_syk=0.3, n_syk=1.2,
  kEl_efgar=0.462, V_efgar=3.5, Emax_FcRn=0.80, EC50_FcRn=30, n_FcRn=1.0,
  Emax_ster_Ab=0.75, EC50_ster_Ab=1.0,
  Emax_ster_Treg=0.6, EC50_ster_Treg=0.8,
  Emax_ster_Mac=0.7, EC50_ster_Mac=1.5,
  Emax_ivig_Mac=0.85, EC50_ivig_Mac=8.0,
  Emax_ivig_FcRn=0.6, EC50_ivig_FcRn=6.0
)

default_state <- c(
  PLT=200, PLT_SP=66, TPO=1.0, MKP=1.0, MK=1.0, Ab=0.2,
  Bc=1.0, Treg=1.0, Th17=1.0, Mac=1.0,
  PRED_c=0, IVIG_c=0, ROMI_sc=0, ROMI_c=0,
  RTX_c=0, RTX_p=0, FOSTA_gut=0, R788_c=0, EFGAR_c=0, ELTP_c=0
)

# ─── ODE solver wrapper ──────────────────────────────────────────────────────
run_sim <- function(parms, state, sim_days, drug_events) {
  times  <- seq(0, sim_days, by = 0.5)
  tryCatch({
    out <- ode(y = state, times = times, func = itp_ode, parms = parms,
               method = "lsoda", events = list(data = drug_events))
    as.data.frame(out)
  }, error = function(e) {
    data.frame(time = times, PLT = rep(NA, length(times)))
  })
}

build_events <- function(sim_days, wt_kg,
                         pred_dose, pred_days,
                         ivig_g_kg,
                         romi_mcg_kg,
                         eltp_mg,
                         rtx_mg, rtx_start,
                         fosta_mg,
                         efgar_mg, efgar_start) {
  ev_list <- list()

  if (pred_dose > 0) {
    d <- seq(0, pred_days, by = 1)
    dose_v <- pmax(pred_dose * 0.5^(pmax(d - 14, 0) / 7), 5)
    ev_list[[1]] <- data.frame(var="PRED_c", time=d, value=dose_v, method="add")
  }
  if (ivig_g_kg > 0) {
    total_g <- ivig_g_kg * wt_kg
    ev_list[[2]] <- data.frame(var="IVIG_c", time=c(0,1),
                               value=total_g/2, method="add")
  }
  if (romi_mcg_kg > 0) {
    t_romi <- seq(0, sim_days, by=7)
    ev_list[[3]] <- data.frame(var="ROMI_sc",
                               time=t_romi,
                               value=romi_mcg_kg*wt_kg*1000,
                               method="add")
  }
  if (eltp_mg > 0) {
    t_eltp <- seq(0, sim_days, by=1)
    ev_list[[4]] <- data.frame(var="ELTP_c",
                               time=t_eltp,
                               value=eltp_mg*1e6*0.52,
                               method="add")
  }
  if (rtx_mg > 0) {
    t_rtx <- rtx_start + c(0,7,14,21)
    t_rtx <- t_rtx[t_rtx <= sim_days]
    ev_list[[5]] <- data.frame(var="RTX_c", time=t_rtx,
                               value=rtx_mg, method="add")
  }
  if (fosta_mg > 0) {
    t_fosta_am <- seq(0, sim_days, by=1)
    t_fosta_pm <- seq(0.5, sim_days, by=1)
    ev_list[[6]] <- data.frame(var="FOSTA_gut",
                               time=c(t_fosta_am, t_fosta_pm),
                               value=fosta_mg/2, method="add")
  }
  if (efgar_mg > 0) {
    t_efgar <- c(seq(efgar_start, efgar_start+21, by=7),
                 seq(efgar_start+42, sim_days, by=21))
    t_efgar <- t_efgar[t_efgar <= sim_days]
    ev_list[[7]] <- data.frame(var="EFGAR_c", time=t_efgar,
                               value=efgar_mg, method="add")
  }

  df <- do.call(rbind, ev_list)
  if (is.null(df) || nrow(df) == 0) {
    return(data.frame(var=character(), time=numeric(), value=numeric(), method=character()))
  }
  df[order(df$time), ]
}

# ── Scenario definitions ──────────────────────────────────────────────────────
scenario_defs <- list(
  "Untreated ITP" = list(
    ITP_sev=4.0, pred_dose=0, pred_days=0, ivig_g_kg=0,
    romi_mcg_kg=0, eltp_mg=0, rtx_mg=0, rtx_start=0,
    fosta_mg=0, efgar_mg=0, efgar_start=0),
  "Prednisone 1 mg/kg/d" = list(
    ITP_sev=4.0, pred_dose=70, pred_days=56, ivig_g_kg=0,
    romi_mcg_kg=0, eltp_mg=0, rtx_mg=0, rtx_start=0,
    fosta_mg=0, efgar_mg=0, efgar_start=0),
  "IVIG rescue" = list(
    ITP_sev=4.0, pred_dose=0, pred_days=0, ivig_g_kg=1.0,
    romi_mcg_kg=0, eltp_mg=0, rtx_mg=0, rtx_start=0,
    fosta_mg=0, efgar_mg=0, efgar_start=0),
  "Romiplostim 3 μg/kg" = list(
    ITP_sev=3.5, pred_dose=0, pred_days=0, ivig_g_kg=0,
    romi_mcg_kg=3, eltp_mg=0, rtx_mg=0, rtx_start=0,
    fosta_mg=0, efgar_mg=0, efgar_start=0),
  "Rituximab 375 mg/m²" = list(
    ITP_sev=4.0, pred_dose=70, pred_days=28, ivig_g_kg=0,
    romi_mcg_kg=0, eltp_mg=0, rtx_mg=700, rtx_start=28,
    fosta_mg=0, efgar_mg=0, efgar_start=0),
  "Fostamatinib 150 mg BID" = list(
    ITP_sev=3.5, pred_dose=0, pred_days=0, ivig_g_kg=0,
    romi_mcg_kg=0, eltp_mg=0, rtx_mg=0, rtx_start=0,
    fosta_mg=300, efgar_mg=0, efgar_start=0),
  "Efgartigimod 700 mg q1w" = list(
    ITP_sev=4.0, pred_dose=0, pred_days=0, ivig_g_kg=0,
    romi_mcg_kg=0, eltp_mg=0, rtx_mg=0, rtx_start=0,
    fosta_mg=0, efgar_mg=700, efgar_start=0),
  "Eltrombopag 50 mg/d" = list(
    ITP_sev=3.5, pred_dose=0, pred_days=0, ivig_g_kg=0,
    romi_mcg_kg=0, eltp_mg=50, rtx_mg=0, rtx_start=0,
    fosta_mg=0, efgar_mg=0, efgar_start=0)
)

## ── UI ───────────────────────────────────────────────────────────────────────
ui <- dashboardPage(
  skin = "purple",
  dashboardHeader(title = "ITP QSP Dashboard"),
  dashboardSidebar(
    sidebarMenu(
      menuItem("Patient Profile",   tabName = "tab_profile",   icon = icon("user")),
      menuItem("Pharmacokinetics",  tabName = "tab_pk",        icon = icon("flask")),
      menuItem("PD Biomarkers",     tabName = "tab_pd",        icon = icon("chart-line")),
      menuItem("Clinical Endpoints",tabName = "tab_endpoints", icon = icon("heartbeat")),
      menuItem("Scenario Comparison",tabName="tab_scenario",   icon = icon("sliders-h")),
      menuItem("Mechanistic Panel", tabName = "tab_mech",      icon = icon("network-wired"))
    ),
    hr(),
    h5("  Global Settings", style = "color:white; margin-left:10px;"),
    sliderInput("sim_days", "Simulation Duration (days)", 30, 365, 180, step = 10),
    sliderInput("ITP_sev",  "ITP Severity (1–6)", 1.0, 6.0, 4.0, step = 0.5),
    sliderInput("wt_kg",    "Patient weight (kg)", 40, 120, 70, step = 5)
  ),
  dashboardBody(
    tags$head(tags$style(HTML(
      ".main-header .logo { font-size: 16px; }
       .info-box { min-height: 70px; }
       .info-box-icon { height: 70px; line-height: 70px; }"))),
    tabItems(

      ## ── TAB 1: Patient Profile ─────────────────────────────────────────────
      tabItem("tab_profile",
        fluidRow(
          box(title="Disease Overview", status="purple", solidHeader=TRUE, width=7,
            p("Immune Thrombocytopenic Purpura (ITP) is an acquired autoimmune disorder
              characterised by isolated thrombocytopenia (platelet count < 100 × 10⁹/L) in
              the absence of other causes. It results from anti-platelet antibody-mediated
              destruction predominantly in the spleen (FcγR phagocytosis) and from impaired
              platelet production due to megakaryocyte suppression."),
            p(strong("Phases: "), "Newly diagnosed (<3 mo), Persistent (3–12 mo),
              Chronic (>12 mo). Severe ITP: PLT < 20 × 10⁹/L."),
            p(strong("First-line:"), " Corticosteroids (prednisone 1 mg/kg/d); IVIG for urgent rise."),
            p(strong("Second-line:"), " TPO-RAs (romiplostim, eltrombopag, avatrombopag),
              Rituximab (anti-CD20), Fostamatinib (SYK inhibitor), Splenectomy."),
            p(strong("Emerging:"), " FcRn inhibitors (efgartigimod, rozanolixizumab),
              Daratumumab (anti-CD38 BMPC depletion).")
          ),
          box(title="QSP Model Structure", status="purple", solidHeader=TRUE, width=5,
            tags$ul(
              tags$li(strong("20-compartment ODE"), ": platelet kinetics, megakaryopoiesis,
                immune biology, multi-drug PK"),
              tags$li(strong("Biological pathways:"), " TPO → c-Mpl → MK → PLT;
                B cells → Anti-GPIIb/IIIa IgG → FcγR-dependent splenic destruction"),
              tags$li(strong("Immune regulation:"), " Treg (↓) / Th17 (↑) imbalance,
                macrophage activation, complement co-opsonization"),
              tags$li(strong("6 drug classes:"), " CS, IVIG, TPO-RA, RTX, SYKi, FcRni"),
              tags$li(strong("Parameters calibrated to:"),
                " RAISE (Lancet 2008/2011), FIT 1+2 (AJH 2018), ADVANCE IV (LH 2022)")
            )
          )
        ),
        fluidRow(
          valueBoxOutput("vbox_plt",  width=3),
          valueBoxOutput("vbox_ab",   width=3),
          valueBoxOutput("vbox_tpo",  width=3),
          valueBoxOutput("vbox_treg", width=3)
        ),
        fluidRow(
          box(title="Patient Parameters", status="info", solidHeader=TRUE, width=4,
            selectInput("patient_type", "Disease Severity",
                        choices = c("Mild (PLT 50–100)","Moderate (PLT 20–50)",
                                    "Severe (PLT < 20)","Normal control"),
                        selected = "Severe (PLT < 20)"),
            numericInput("age_yr",   "Age (years)",   45, 18, 90),
            selectInput("sex",       "Sex",     c("Female","Male")),
            checkboxInput("prev_steroids", "Prior steroid exposure", TRUE),
            checkboxInput("prev_rtx",      "Prior rituximab", FALSE)
          ),
          box(title="Diagnosis Summary", status="info", solidHeader=TRUE, width=8,
            DTOutput("diag_table")
          )
        )
      ),

      ## ── TAB 2: Pharmacokinetics ────────────────────────────────────────────
      tabItem("tab_pk",
        fluidRow(
          box(title="Drug PK Settings", status="navy", solidHeader=TRUE, width=3,
            selectInput("pk_drug", "Select Drug",
                        choices = c("Romiplostim","Eltrombopag","Rituximab",
                                    "Fostamatinib / R788","Efgartigimod","Prednisone","IVIG"),
                        selected = "Romiplostim"),
            sliderInput("pk_dose_romi",  "Romiplostim (μg/kg/week)", 1, 10, 3, 0.5),
            sliderInput("pk_dose_eltp",  "Eltrombopag (mg/day)", 25, 75, 50, 25),
            sliderInput("pk_dose_rtx",   "Rituximab dose (mg)", 200, 1000, 700, 100),
            sliderInput("pk_dose_pred",  "Prednisone (mg/day)", 10, 100, 70, 10),
            sliderInput("pk_dose_efgar", "Efgartigimod (mg)", 400, 1200, 700, 100),
            sliderInput("pk_dose_fosta", "Fostamatinib total mg/day", 100, 300, 300, 100)
          ),
          box(title="Concentration-Time Profile", status="navy", solidHeader=TRUE, width=9,
            plotlyOutput("pk_plot", height="450px")
          )
        ),
        fluidRow(
          box(title="PK Parameter Summary", status="info", width=12,
            DTOutput("pk_table"))
        )
      ),

      ## ── TAB 3: PD Biomarkers ──────────────────────────────────────────────
      tabItem("tab_pd",
        fluidRow(
          box(title="Treatment Inputs", status="olive", solidHeader=TRUE, width=3,
            selectInput("pd_drug", "Treatment",
                        choices = names(scenario_defs),
                        selected = "Prednisone 1 mg/kg/d"),
            hr(),
            h5("Custom overrides"),
            sliderInput("pd_pred",  "Prednisone (mg/day)", 0, 100, 0, 10),
            sliderInput("pd_romi",  "Romiplostim (μg/kg/wk)", 0, 10, 0, 0.5),
            sliderInput("pd_eltp",  "Eltrombopag (mg/day)", 0, 75, 0, 25),
            checkboxInput("pd_ivig", "IVIG 1 g/kg day 0–1", FALSE),
            checkboxInput("pd_rtx",  "Rituximab ×4 weekly", FALSE),
            checkboxInput("pd_fosta","Fostamatinib 150 mg BID", FALSE),
            checkboxInput("pd_efgar","Efgartigimod 700 mg q1w", FALSE),
            actionButton("pd_run", "Run Simulation", class="btn-success btn-block")
          ),
          box(title="Core PD: Platelet Count + Key Biomarkers", status="olive",
              solidHeader=TRUE, width=9,
              plotlyOutput("pd_plot_main", height="500px"))
        )
      ),

      ## ── TAB 4: Clinical Endpoints ──────────────────────────────────────────
      tabItem("tab_endpoints",
        fluidRow(
          box(title="Clinical Response Thresholds", status="red", solidHeader=TRUE, width=12,
            plotlyOutput("endpoint_plot", height="400px"))
        ),
        fluidRow(
          box(title="Response Rate Summary (Week 12 / Week 24)", status="red",
              solidHeader=TRUE, width=7,
              DTOutput("resp_table")),
          box(title="Bleeding Risk Assessment", status="red",
              solidHeader=TRUE, width=5,
              plotlyOutput("bleed_gauge", height="300px"))
        )
      ),

      ## ── TAB 5: Scenario Comparison ────────────────────────────────────────
      tabItem("tab_scenario",
        fluidRow(
          box(title="Select Scenarios to Compare", status="maroon",
              solidHeader=TRUE, width=3,
            checkboxGroupInput("scen_select", "Scenarios",
                               choices = names(scenario_defs),
                               selected = c("Untreated ITP","Prednisone 1 mg/kg/d",
                                            "Romiplostim 3 μg/kg",
                                            "Fostamatinib 150 mg BID",
                                            "Efgartigimod 700 mg q1w",
                                            "Rituximab 375 mg/m²")),
            hr(),
            actionButton("scen_run", "Run All Scenarios", class="btn-warning btn-block")
          ),
          box(title="Platelet Count — All Scenarios", status="maroon",
              solidHeader=TRUE, width=9,
              plotlyOutput("scen_plt_plot", height="500px"))
        ),
        fluidRow(
          box(title="Anti-platelet Ab", status="light", width=4,
              plotlyOutput("scen_ab_plot", height="300px")),
          box(title="Treg:Th17 Ratio", status="light", width=4,
              plotlyOutput("scen_treg_plot", height="300px")),
          box(title="Macrophage Activation", status="light", width=4,
              plotlyOutput("scen_mac_plot", height="300px"))
        )
      ),

      ## ── TAB 6: Mechanistic Panel ───────────────────────────────────────────
      tabItem("tab_mech",
        fluidRow(
          box(title="Select Scenario", status="teal", solidHeader=TRUE, width=3,
            selectInput("mech_scen", "Scenario",
                        choices = names(scenario_defs),
                        selected = "Romiplostim 3 μg/kg"),
            actionButton("mech_run", "Simulate", class="btn-info btn-block"),
            hr(),
            p("Panel shows all mechanistic variables from the 20-CMT ODE system.")
          ),
          box(title="Platelet + TPO + Megakaryocytes", status="teal",
              solidHeader=TRUE, width=9,
              plotlyOutput("mech_plt_tpo", height="350px"))
        ),
        fluidRow(
          box(title="Immune: Ab / Bc / Treg / Th17 / Mac",
              status="teal", solidHeader=TRUE, width=6,
              plotlyOutput("mech_immune", height="350px")),
          box(title="Megakaryopoiesis (MKP → MK) + BM Inhibition",
              status="teal", solidHeader=TRUE, width=6,
              plotlyOutput("mech_mk", height="350px"))
        )
      )

    ) # end tabItems
  )  # end dashboardBody
)   # end dashboardPage

## ── Server ───────────────────────────────────────────────────────────────────
server <- function(input, output, session) {

  wt_kg_r <- reactive(input$wt_kg)

  # ── Helper: run a named scenario ─────────────────────────────────────────
  sim_scenario <- function(sname, sim_days, ITP_sev_override = NULL) {
    def <- scenario_defs[[sname]]
    p   <- default_parms
    p["ITP_sev"] <- if (!is.null(ITP_sev_override)) ITP_sev_override else def$ITP_sev
    wt  <- wt_kg_r()
    evts <- build_events(sim_days, wt,
                         def$pred_dose, def$pred_days,
                         def$ivig_g_kg, def$romi_mcg_kg,
                         def$eltp_mg, def$rtx_mg, def$rtx_start,
                         def$fosta_mg, def$efgar_mg, def$efgar_start)
    df <- run_sim(p, default_state, sim_days, evts)
    df$scenario <- sname
    df$Treg_Th17 <- df$Treg / (df$Th17 + 1e-4)
    df
  }

  # ── Value boxes ─────────────────────────────────────────────────────────
  output$vbox_plt  <- renderValueBox(
    valueBox("< 20 ×10⁹/L", "Platelet Count (Severe ITP)", color="red", icon=icon("tint")))
  output$vbox_ab   <- renderValueBox(
    valueBox("↑ 3–5×", "Anti-platelet IgG", color="orange", icon=icon("bacteria")))
  output$vbox_tpo  <- renderValueBox(
    valueBox("↑ 2–4×", "Thrombopoietin", color="blue", icon=icon("arrow-up")))
  output$vbox_treg <- renderValueBox(
    valueBox("↓ 50%", "Treg:Th17 Ratio", color="purple", icon=icon("balance-scale")))

  # ── Diagnosis table ──────────────────────────────────────────────────────
  output$diag_table <- renderDT({
    data.frame(
      Parameter   = c("Platelet count","Bleeding time","PT/APTT","Bone marrow",
                      "Anti-platelet Ab (MAIPA)","Anti-GPIIb/IIIa IgG",
                      "Treg (CD4+FoxP3+)","Th17 (CD4+IL-17+)","BAFF","TPO"),
      FindingNormal= c("150–400","Normal","Normal","Normal megakaryocytes",
                       "Negative","Undetectable","~5–7% CD4","~1–2% CD4",
                       "Normal","Normal"),
      FindingITP   = c("< 30 (severe)","Prolonged","Normal","MK ↑ / abnormal maturation",
                       "Positive (60–80%)","↑ 3–5× baseline","↓ 50–70%","↑ 2–3×",
                       "Elevated","↑ (PLT feedback)"),
      Significance = c("Primary endpoint","Bleeding risk","Rule out coagulopathy",
                       "Megakaryocyte dysfunction","Diagnostic","Most specific Ab",
                       "Immune dysregulation","Drives B cell activation",
                       "B cell survival","Compensatory ↑"),
      stringsAsFactors = FALSE
    )
  }, options=list(pageLength=10, dom='t'), rownames=FALSE)

  # ── PK simulation ─────────────────────────────────────────────────────────
  pk_sim <- reactive({
    wt <- wt_kg_r()
    sd <- input$sim_days
    p  <- default_parms
    drug <- input$pk_drug

    evts <- switch(drug,
      "Romiplostim"       = build_events(sd, wt, 0,0,0, input$pk_dose_romi, 0,0,0,0,0,0),
      "Eltrombopag"       = build_events(sd, wt, 0,0,0, 0, input$pk_dose_eltp, 0,0,0,0,0),
      "Rituximab"         = build_events(sd, wt, 0,0,0, 0, 0, input$pk_dose_rtx, 0,0,0,0),
      "Fostamatinib / R788"=build_events(sd, wt, 0,0,0, 0, 0,0,0, input$pk_dose_fosta,0,0),
      "Efgartigimod"      = build_events(sd, wt, 0,0,0, 0, 0,0,0,0, input$pk_dose_efgar, 0),
      "Prednisone"        = build_events(sd, wt, input$pk_dose_pred, 28, 0,0,0,0,0,0,0,0),
      "IVIG"              = build_events(sd, wt, 0,0, 1, 0,0,0,0,0,0,0),
      build_events(sd, wt, 0,0,0,0,0,0,0,0,0,0)
    )
    df <- run_sim(p, default_state, sd, evts)
    df
  })

  output$pk_plot <- renderPlotly({
    df   <- pk_sim()
    drug <- input$pk_drug
    yvar <- switch(drug,
      "Romiplostim"        = "ROMI_c",
      "Eltrombopag"        = "ELTP_c",
      "Rituximab"          = "RTX_c",
      "Fostamatinib / R788"= "R788_c",
      "Efgartigimod"       = "EFGAR_c",
      "Prednisone"         = "PRED_c",
      "IVIG"               = "IVIG_c",
      "ROMI_c"
    )
    p <- ggplot(df, aes_string("time", yvar)) +
      geom_line(color="#3C5488", linewidth=1.2) +
      labs(title=paste("PK:", drug), x="Time (days)",
           y=paste(drug, "Amount / Conc.")) +
      theme_bw(base_size=12)
    ggplotly(p)
  })

  output$pk_table <- renderDT({
    data.frame(
      Drug=c("Prednisone","IVIG","Romiplostim","Eltrombopag","Rituximab",
             "Fostamatinib (R788)","Efgartigimod"),
      Route=c("PO","IV","SC","PO","IV","PO","IV"),
      Dose=c("1 mg/kg/d × 4w","1 g/kg × 1–2d","1–10 μg/kg/wk","25–75 mg/d",
             "375 mg/m² × 4w","100–150 mg BID","10 mg/kg q1w × 4"),
      `t½`=c("2–3h","21d","~1–4d","~16–35h","22d (terminal)","~13–14h R788","~34h"),
      F=c("90%","100% (IV)","~66%","~52%","100% (IV)","~70% as R788","100% (IV)"),
      `Mechanism (PD)`=c("B cell suppress, FcR↓, Treg↑","FcRn sat, FcγR block, Treg↑",
                          "c-Mpl agonist → MK↑","c-Mpl agonist (non-peptide) → MK↑",
                          "CD20 → B cell depletion","SYK inhibitor → phagocytosis↓",
                          "FcRn blockade → IgG catabolism↑"),
      stringsAsFactors=FALSE
    )
  }, options=list(pageLength=10, dom='t'), rownames=FALSE)

  # ── PD simulation ─────────────────────────────────────────────────────────
  pd_data <- eventReactive(input$pd_run, {
    sname <- input$pd_scen %||% input$pd_drug
    def   <- scenario_defs[[sname]]
    p     <- default_parms
    p["ITP_sev"] <- input$ITP_sev
    wt    <- wt_kg_r()
    sd    <- input$sim_days
    evts  <- build_events(sd, wt,
      if (input$pd_pred > 0) input$pd_pred else def$pred_dose,
      if (input$pd_pred > 0) 56 else def$pred_days,
      if (input$pd_ivig) 1.0 else def$ivig_g_kg,
      input$pd_romi, input$pd_eltp,
      if (input$pd_rtx)   700 else def$rtx_mg, def$rtx_start,
      if (input$pd_fosta) 300 else def$fosta_mg,
      if (input$pd_efgar) 700 else def$efgar_mg, def$efgar_start)
    df <- run_sim(p, default_state, sd, evts)
    df$Treg_Th17 <- df$Treg / (df$Th17 + 1e-4)
    df
  }, ignoreNULL=FALSE)

  output$pd_plot_main <- renderPlotly({
    df <- pd_data()
    if (is.null(df) || !"PLT" %in% names(df)) return(NULL)
    df_long <- df %>%
      select(time, PLT, Ab, Treg_Th17, TPO, MK) %>%
      pivot_longer(-time, names_to="var", values_to="val")
    p <- ggplot(df_long, aes(time, val, color=var)) +
      geom_line(linewidth=1.0) +
      facet_wrap(~var, scales="free_y", nrow=2) +
      scale_color_brewer(palette="Set1") +
      labs(title="Core PD Biomarkers", x="Time (days)", y="Value", color=NULL) +
      theme_bw(base_size=11) + theme(legend.position="none")
    ggplotly(p)
  })

  # ── Scenario comparison ───────────────────────────────────────────────────
  scen_data <- eventReactive(input$scen_run, {
    sels <- input$scen_select
    sd   <- input$sim_days
    isev <- input$ITP_sev
    lapply(sels, function(s) sim_scenario(s, sd, isev)) %>% bind_rows()
  }, ignoreNULL=FALSE)

  cols8 <- c("#E64B35","#4DBBD5","#00A087","#3C5488",
             "#F39B7F","#8491B4","#91D1C2","#DC0000")

  output$scen_plt_plot <- renderPlotly({
    df <- scen_data()
    if (is.null(df) || !"PLT" %in% names(df)) return(NULL)
    p <- ggplot(df, aes(time, PLT, color=scenario)) +
      geom_line(linewidth=1.1) +
      geom_hline(yintercept=c(20,30,100), linetype="dashed",
                 color=c("#ee2222","#ffaa00","#009900")) +
      scale_color_manual(values=cols8) +
      labs(title="Platelet Count by Scenario", x="Days", y="PLT (×10⁹/L)", color=NULL) +
      theme_bw(base_size=11)
    ggplotly(p)
  })

  output$scen_ab_plot <- renderPlotly({
    df <- scen_data()
    if (is.null(df)) return(NULL)
    p <- ggplot(df, aes(time, Ab, color=scenario)) +
      geom_line(linewidth=1.0) +
      scale_color_manual(values=cols8) +
      labs(title="Anti-platelet Ab", x="Days", y="Ab (AU)", color=NULL) +
      theme_bw(base_size=10) + theme(legend.position="none")
    ggplotly(p)
  })

  output$scen_treg_plot <- renderPlotly({
    df <- scen_data()
    if (is.null(df)) return(NULL)
    p <- ggplot(df, aes(time, Treg_Th17, color=scenario)) +
      geom_line(linewidth=1.0) +
      geom_hline(yintercept=1, linetype="dashed", color="grey50") +
      scale_color_manual(values=cols8) +
      labs(title="Treg:Th17 Ratio", x="Days", y="Treg/Th17", color=NULL) +
      theme_bw(base_size=10) + theme(legend.position="none")
    ggplotly(p)
  })

  output$scen_mac_plot <- renderPlotly({
    df <- scen_data()
    if (is.null(df)) return(NULL)
    p <- ggplot(df, aes(time, Mac, color=scenario)) +
      geom_line(linewidth=1.0) +
      scale_color_manual(values=cols8) +
      labs(title="Macrophage Activation", x="Days", y="Mac (AU)", color=NULL) +
      theme_bw(base_size=10) + theme(legend.position="none")
    ggplotly(p)
  })

  # ── Endpoints ─────────────────────────────────────────────────────────────
  output$endpoint_plot <- renderPlotly({
    sels <- isolate(input$scen_select) %||% c("Untreated ITP","Prednisone 1 mg/kg/d",
                                               "Romiplostim 3 μg/kg")
    df <- lapply(sels, function(s)
      sim_scenario(s, input$sim_days, input$ITP_sev)) %>% bind_rows()
    p <- ggplot(df, aes(time, PLT, color=scenario)) +
      geom_line(linewidth=1.1) +
      geom_ribbon(aes(ymin=0, ymax=20), fill="#ee222220", color=NA) +
      geom_hline(yintercept=c(20,30,100), linetype="dashed",
                 color=c("#ee2222","#ffaa00","#009900")) +
      scale_color_manual(values=cols8) +
      labs(title="PLT Response vs. Threshold", x="Days", y="PLT ×10⁹/L", color=NULL) +
      theme_bw(base_size=11)
    ggplotly(p)
  })

  output$resp_table <- renderDT({
    sels <- c("Untreated ITP","Prednisone 1 mg/kg/d","Romiplostim 3 μg/kg",
              "Fostamatinib 150 mg BID","Rituximab 375 mg/m²","Efgartigimod 700 mg q1w")
    df <- lapply(sels, function(s)
      sim_scenario(s, input$sim_days, input$ITP_sev)) %>% bind_rows()
    df %>%
      filter(time %in% c(84, 168)) %>%
      mutate(Week = ifelse(time <= 90, "Wk 12", "Wk 24"),
             PLT_r = round(PLT, 0),
             CR = PLT >= 100, R = PLT >= 30 & PLT < 100,
             NR = PLT < 30) %>%
      group_by(Scenario=scenario, Week) %>%
      summarise(PLT=round(mean(PLT_r),0),
                CR=any(CR), R=any(R), NR=any(NR), .groups="drop") %>%
      datatable(options=list(pageLength=12, dom='t'), rownames=FALSE)
  })

  output$bleed_gauge <- renderPlotly({
    df <- sim_scenario("Untreated ITP", input$sim_days, input$ITP_sev)
    plt_last <- tail(df$PLT, 1)
    risk_pct  <- max(0, min(100, (30 - plt_last) / 30 * 100))
    plot_ly(type="indicator", mode="gauge+number",
            value=round(risk_pct,1),
            title=list(text="Estimated Bleeding Risk (%)"),
            gauge=list(axis=list(range=list(0,100)),
                       bar=list(color="darkred"),
                       steps=list(
                         list(range=c(0,30),   color="green"),
                         list(range=c(30,60),  color="orange"),
                         list(range=c(60,100), color="red")
                       )))
  })

  # ── Mechanistic panel ────────────────────────────────────────────────────
  mech_data <- eventReactive(input$mech_run, {
    sim_scenario(input$mech_scen, input$sim_days, input$ITP_sev)
  }, ignoreNULL=FALSE)

  output$mech_plt_tpo <- renderPlotly({
    df <- mech_data()
    df_long <- df %>%
      select(time, PLT, TPO, MK, MKP) %>%
      pivot_longer(-time, names_to="var", values_to="val")
    p <- ggplot(df_long, aes(time, val, color=var)) +
      geom_line(linewidth=1.1) +
      facet_wrap(~var, scales="free_y") +
      scale_color_manual(values=c("#E64B35","#4DBBD5","#00A087","#8491B4")) +
      labs(title="Platelet / TPO / MK / MKP", x="Days", y="", color=NULL) +
      theme_bw(base_size=11) + theme(legend.position="none")
    ggplotly(p)
  })

  output$mech_immune <- renderPlotly({
    df <- mech_data()
    df$Treg_Th17 <- df$Treg / (df$Th17 + 1e-4)
    df_long <- df %>%
      select(time, Ab, Bc, Treg, Th17, Mac) %>%
      pivot_longer(-time, names_to="var", values_to="val")
    p <- ggplot(df_long, aes(time, val, color=var)) +
      geom_line(linewidth=1.0) +
      facet_wrap(~var, scales="free_y") +
      scale_color_brewer(palette="Set1") +
      labs(title="Immune Variables", x="Days", y="", color=NULL) +
      theme_bw(base_size=10) + theme(legend.position="none")
    ggplotly(p)
  })

  output$mech_mk <- renderPlotly({
    df <- mech_data()
    p <- ggplot(df) +
      geom_line(aes(time, MKP, color="MKP"),  linewidth=1.0) +
      geom_line(aes(time, MK,  color="MK"),   linewidth=1.0) +
      geom_line(aes(time, Ab*10, color="Ab×10"), linewidth=1.0, linetype="dashed") +
      scale_color_manual(values=c(MKP="#4DBBD5",MK="#E64B35","Ab×10"="#F39B7F")) +
      labs(title="Megakaryopoiesis vs. Ab Inhibition",
           x="Days", y="Normalized units", color=NULL) +
      theme_bw(base_size=11)
    ggplotly(p)
  })

}

`%||%` <- function(x, y) if (!is.null(x) && length(x) > 0) x else y

shinyApp(ui, server)
