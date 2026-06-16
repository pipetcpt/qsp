## =============================================================================
## PAH QSP Shiny Application
## Interactive simulation of Pulmonary Arterial Hypertension treatment
##
## Tabs:
##   1. Disease Overview  — pathophysiology diagram and key concepts
##   2. PK Profiles       — drug concentration vs time for each drug class
##   3. PD Response       — cGMP, cAMP, ET-1 dynamics
##   4. Hemodynamics      — PVR, mPAP, CO, RV coupling over time
##   5. Clinical Outcomes — 6MWD, NT-proBNP, WHO FC over time
##   6. Treatment Compare — multi-arm comparison plots
##   7. Sensitivity Anal. — one-at-a-time parameter sensitivity
##   8. References        — key literature
## =============================================================================

library(shiny)
library(shinydashboard)
library(mrgsolve)
library(dplyr)
library(ggplot2)
library(tidyr)
library(DT)
library(plotly)

## ─────────────────────────────────────────────────────────────────────────────
## Source the mrgsolve model (assumes pah_mrgsolve.R is in same directory)
## The pah_code string and mod object should be defined there.
## For standalone Shiny use, we re-embed the model code here.
## ─────────────────────────────────────────────────────────────────────────────
pah_model_code <- '
$PROB PAH QSP Model (Shiny embedded)

$PARAM @annotated
  KA_SIL:0.90:SIL absorption (1/h)
  CL_SIL:40.0:SIL clearance (L/h)
  V2_SIL:100.0:SIL central vol (L)
  Q_SIL:18.0:SIL Q (L/h)
  V3_SIL:200.0:SIL periph vol (L)
  F1_SIL:0.40:SIL bioavailability
  EC50_SIL:30.0:SIL PDE5i EC50 (ng/mL)
  EMAX_SIL:0.95:SIL Emax

  KA_BOS:0.60:BOS absorption (1/h)
  CL_BOS:18.0:BOS clearance (L/h)
  V2_BOS:42.0:BOS central vol (L)
  Q_BOS:8.0:BOS Q (L/h)
  V3_BOS:50.0:BOS periph vol (L)
  F1_BOS:0.50:BOS bioavailability
  EC50_BOS:150.0:BOS ERA EC50 (ng/mL)
  EMAX_BOS:0.90:BOS Emax

  KA_TRE:0.80:TRE absorption (1/h)
  CL_TRE:38.0:TRE clearance (L/h)
  V2_TRE:14.0:TRE vol (L)
  F1_TRE:1.00:TRE bioavailability
  EC50_TRE:0.50:TRE IP EC50 (ng/mL)
  EMAX_TRE:0.80:TRE Emax

  KSY_ET1:0.60:ET-1 synthesis (pmol/L/h)
  KDG_ET1:0.12:ET-1 clearance (1/h)
  STIM_ET1:1.50:ET-1 disease multiplier
  ETB_FRAC:0.30:ETB clearance fraction
  ET1_REF:5.0:Reference ET-1 (pmol/L)

  KPROD_CGMP:0.50:cGMP production (nmol/L/h)
  KDEG_CGMP:0.40:cGMP degradation (1/h)
  CGMP_REF:1.25:Reference cGMP (nmol/L)
  NO_EFF:0.60:NO-driven sGC activity

  KPROD_CAMP:0.80:cAMP production (nmol/L/h)
  KDEG_CAMP:0.60:cAMP degradation (1/h)
  CAMP_REF:1.33:Reference cAMP (nmol/L)
  PGI2_EFF:0.50:Endogenous PGI2 effect

  VT_BASE:0.70:Vascular tone baseline
  KET1_VT:0.15:ET-1 vasoconstriction coeff
  KCGMP_VT:0.20:cGMP vasodilation coeff
  KCAMP_VT:0.15:cAMP vasodilation coeff

  VRI_0:0.40:Initial VRI
  KGROW_VRI:0.002:VRI growth rate (1/h)
  KREG_VRI:0.001:VRI regression rate (1/h)
  VRI_MAX:0.95:Maximum VRI

  PVR_BASE:800.0:Baseline PVR (dyn.s/cm5)
  PAWP_0:8.0:PAWP (mmHg)
  CO_0:4.0:Baseline CO (L/min)
  ALPHA_CO:0.30:CO sensitivity to coupling

  Ees_0:0.70:RV Ees (mmHg/mL)
  Ea_BASE:0.75:RV Ea (mmHg/mL)
  Ees_HYP:0.40:Max Ees from hypertrophy
  TAU_HYP:720.0:RV hypertrophy time constant (h)
  COUP_THRESH:0.80:Decompensation threshold

  BNP_BASE:1000.0:Baseline NT-proBNP (pg/mL)
  BNP_SENS:2.5:NT-proBNP RVSP sensitivity
  RVSP_REF:60.0:Reference RVSP (mmHg)

  WALK_BASE:350.0:Baseline 6MWD (m)
  WALK_MAX:550.0:Max 6MWD (m)
  WALK_SENS:1.8:6MWD PVR sensitivity

$CMT @annotated
  DEPOT_SIL:SIL GI depot
  CENT_SIL:SIL central
  PERI_SIL:SIL peripheral
  DEPOT_BOS:BOS GI depot
  CENT_BOS:BOS central
  PERI_BOS:BOS peripheral
  SC_TRE:TRE SC depot
  CENT_TRE:TRE central
  ET1_PD:ET-1
  CGMP_PD:cGMP
  CAMP_PD:cAMP
  VRI_PD:Vascular Remodeling Index
  RV_HYP:RV hypertrophy index

$MAIN
  double C_SIL = CENT_SIL / V2_SIL * 1000.0;
  double C_BOS = CENT_BOS / V2_BOS * 1000.0;
  double C_TRE = CENT_TRE / V2_TRE * 1000.0;

  double EFF_SIL = EMAX_SIL * pow(C_SIL,1.5)/(pow(EC50_SIL,1.5)+pow(C_SIL,1.5));
  double EFF_BOS = EMAX_BOS * C_BOS/(EC50_BOS+C_BOS);
  double EFF_TRE = EMAX_TRE * C_TRE/(EC50_TRE+C_TRE);

  double ET1_ratio  = ET1_PD  / ET1_REF;
  double CGMP_ratio = CGMP_PD / CGMP_REF;
  double CAMP_ratio = CAMP_PD / CAMP_REF;
  double ETA_block  = 1.0 - EFF_BOS;
  double VT = VT_BASE
    + KET1_VT  * (ET1_ratio  - 1.0) * ETA_block
    - KCGMP_VT * (CGMP_ratio - 1.0)
    - KCAMP_VT * (CAMP_ratio - 1.0);
  if(VT < 0.05) VT = 0.05;
  if(VT > 1.0)  VT = 1.0;

  double PVR   = PVR_BASE * (VT / VT_BASE) * (1.0 + 0.5 * VRI_PD);
  double Ees   = Ees_0 * (1.0 + Ees_HYP * RV_HYP);
  double Ea    = Ea_BASE * (PVR / PVR_BASE);
  double COUP  = (Ea > 0) ? Ees/Ea : 1.0;
  if(COUP > 3.0) COUP = 3.0;

  double CO_adj = CO_0;
  if(COUP < 1.0) CO_adj = CO_0 * (0.4 + 0.6 * COUP);
  else           CO_adj = CO_0 * (1.0 + ALPHA_CO * (COUP - 1.0));
  if(CO_adj < 1.0) CO_adj = 1.0;

  double mPAP  = PAWP_0 + (PVR / 80.0) * CO_adj;
  double RVSP  = mPAP + 5.0;
  double NT_proBNP = BNP_BASE * exp(BNP_SENS * (RVSP - RVSP_REF)/RVSP_REF);
  if(NT_proBNP < 50) NT_proBNP = 50;
  double WALK6 = WALK_MAX - (WALK_MAX - WALK_BASE) * pow(PVR/PVR_BASE, WALK_SENS) /
                 (1.0 + ALPHA_CO * ((COUP > 1.0) ? COUP-1.0 : 0));
  if(WALK6 < 50)  WALK6 = 50;
  if(WALK6 > 650) WALK6 = 650;

  double WHO_FC;
  if(mPAP < 25 && WALK6 > 500)      WHO_FC = 1;
  else if(mPAP < 35 && WALK6 > 400) WHO_FC = 2;
  else if(mPAP < 50 && WALK6 > 250) WHO_FC = 3;
  else                                WHO_FC = 4;

  ET1_PD_0  = (KSY_ET1 * STIM_ET1) / (KDG_ET1 * (1.0 - ETB_FRAC * 0.5));
  CGMP_PD_0 = KPROD_CGMP * NO_EFF / KDEG_CGMP;
  CAMP_PD_0 = KPROD_CAMP / KDEG_CAMP;
  VRI_PD_0  = VRI_0;
  RV_HYP_0  = 0.30;

$ODE
  double C_SIL_o = CENT_SIL / V2_SIL * 1000.0;
  double C_BOS_o = CENT_BOS / V2_BOS * 1000.0;
  double C_TRE_o = CENT_TRE / V2_TRE * 1000.0;
  double E_SIL = EMAX_SIL*pow(C_SIL_o,1.5)/(pow(EC50_SIL,1.5)+pow(C_SIL_o,1.5));
  double E_BOS = EMAX_BOS*C_BOS_o/(EC50_BOS+C_BOS_o);
  double E_TRE = EMAX_TRE*C_TRE_o/(EC50_TRE+C_TRE_o);

  dxdt_DEPOT_SIL = -KA_SIL * DEPOT_SIL;
  dxdt_CENT_SIL  =  KA_SIL*DEPOT_SIL*F1_SIL - (CL_SIL+Q_SIL)/V2_SIL*CENT_SIL + Q_SIL/V3_SIL*PERI_SIL;
  dxdt_PERI_SIL  =  Q_SIL/V2_SIL*CENT_SIL - Q_SIL/V3_SIL*PERI_SIL;
  dxdt_DEPOT_BOS = -KA_BOS*DEPOT_BOS;
  dxdt_CENT_BOS  =  KA_BOS*DEPOT_BOS*F1_BOS - (CL_BOS+Q_BOS)/V2_BOS*CENT_BOS + Q_BOS/V3_BOS*PERI_BOS;
  dxdt_PERI_BOS  =  Q_BOS/V2_BOS*CENT_BOS - Q_BOS/V3_BOS*PERI_BOS;
  dxdt_SC_TRE    = -KA_TRE*SC_TRE;
  dxdt_CENT_TRE  =  KA_TRE*SC_TRE*F1_TRE - CL_TRE/V2_TRE*CENT_TRE;

  double ETB_clr = ETB_FRAC*KDG_ET1*(1.0-E_BOS*0.5);
  dxdt_ET1_PD = KSY_ET1*STIM_ET1 - (KDG_ET1*(1.0-ETB_FRAC)+ETB_clr)*ET1_PD;

  dxdt_CGMP_PD = KPROD_CGMP*NO_EFF - KDEG_CGMP*(1.0-E_SIL)*CGMP_PD;

  double IP_stim = PGI2_EFF + E_TRE;
  dxdt_CAMP_PD = KPROD_CAMP*(IP_stim/PGI2_EFF) - KDEG_CAMP*CAMP_PD;

  double ET1_r = ET1_PD/ET1_REF, CG_r = CGMP_PD/CGMP_REF, CA_r = CAMP_PD/CAMP_REF;
  double PROL = fmax(0.0, 0.5*(ET1_r-1.0)-0.3*(CG_r-1.0)-0.2*(CA_r-1.0));
  double REGR = E_BOS*0.4 + E_SIL*0.3 + E_TRE*0.3;
  dxdt_VRI_PD = KGROW_VRI*PROL*(VRI_MAX-VRI_PD) - KREG_VRI*REGR*VRI_PD;

  double PVR_o = PVR_BASE*(((VT_BASE + KET1_VT*(ET1_r-1.0)*(1.0-E_BOS) - KCGMP_VT*(CG_r-1.0) - KCAMP_VT*(CA_r-1.0))/VT_BASE))*(1.0+0.5*VRI_PD);
  double Ea_o  = Ea_BASE*(PVR_o/PVR_BASE);
  double Ees_o = Ees_0*(1.0+Ees_HYP*RV_HYP);
  double COUP_o = (Ea_o > 0) ? Ees_o/Ea_o : 1.0;
  if(COUP_o > 3.0) COUP_o = 3.0;
  double HS = fmax(0.0, 1.0-COUP_o);
  double HR = fmax(0.0, COUP_o-1.0)*0.5;
  dxdt_RV_HYP = (1.0/TAU_HYP)*(HS*(1.0-RV_HYP) - HR*RV_HYP);

$CAPTURE
  PVR CGMP_PD CAMP_PD ET1_PD VRI_PD RV_HYP
  mPAP RVSP CO_adj COUP NT_proBNP WALK6 WHO_FC
  E_SIL E_BOS E_TRE C_SIL_o C_BOS_o C_TRE_o
'

## ─────────────────────────────────────────────────────────────────────────────
## Compile model once at startup
## ─────────────────────────────────────────────────────────────────────────────
mod_global <- mcode("PAH_Shiny", pah_model_code, soloc = tempdir())

## ─────────────────────────────────────────────────────────────────────────────
## Shared theme
## ─────────────────────────────────────────────────────────────────────────────
theme_qsp <- theme_bw(base_size = 13) +
  theme(
    legend.position  = "bottom",
    legend.title     = element_blank(),
    panel.grid.minor = element_blank(),
    strip.background = element_rect(fill = "#2C3E50"),
    strip.text       = element_text(colour = "white", face = "bold")
  )

## ─────────────────────────────────────────────────────────────────────────────
## Helper: run one simulation arm
## ─────────────────────────────────────────────────────────────────────────────
run_sim <- function(mod, n_weeks,
                    dose_sil = 0, dose_bos = 0, dose_tre = 0,
                    extra_params = list()) {
  hours <- n_weeks * 7 * 24
  ev_list <- list()
  if (dose_sil > 0)
    ev_list[["sil"]] <- ev(amt = dose_sil, ii = 8,
                           addl = ceiling(hours / 8),
                           cmt = "DEPOT_SIL", evid = 1)
  if (dose_bos > 0)
    ev_list[["bos"]] <- ev(amt = dose_bos, ii = 12,
                           addl = ceiling(hours / 12),
                           cmt = "DEPOT_BOS", evid = 1)
  if (dose_tre > 0)
    ev_list[["tre"]] <- ev(amt = dose_tre * 24, ii = 24,
                           addl = ceiling(hours / 24),
                           cmt = "SC_TRE", evid = 1)
  ev_final <- if (length(ev_list) > 0) do.call(c, ev_list) else ev(amt = 0, cmt = "DEPOT_SIL", evid = 2)

  p_list <- c(list(STIM_ET1 = 1.5), extra_params)
  do.call(param, c(list(mod), p_list)) %>%
    mrgsim(ev = ev_final, end = hours, delta = 24) %>%
    as_tibble() %>%
    mutate(Week = time / (7 * 24))
}

## ─────────────────────────────────────────────────────────────────────────────
## UI
## ─────────────────────────────────────────────────────────────────────────────
ui <- dashboardPage(
  skin = "blue",

  dashboardHeader(
    title = span(icon("lungs"), "PAH QSP Simulator"),
    titleWidth = 300
  ),

  dashboardSidebar(
    width = 280,
    sidebarMenu(
      id = "tabs",
      menuItem("Disease Overview",   tabName = "overview",   icon = icon("info-circle")),
      menuItem("PK Profiles",        tabName = "pk",         icon = icon("chart-line")),
      menuItem("PD Markers",         tabName = "pd",         icon = icon("flask")),
      menuItem("Hemodynamics",       tabName = "hemo",       icon = icon("heartbeat")),
      menuItem("Clinical Outcomes",  tabName = "clinical",   icon = icon("user-md")),
      menuItem("Treatment Compare",  tabName = "compare",    icon = icon("balance-scale")),
      menuItem("Sensitivity Anal.",  tabName = "sensitivity",icon = icon("sliders-h")),
      menuItem("References",         tabName = "refs",       icon = icon("book"))
    ),

    hr(),
    ## ── Global Controls ───────────────────────────────────────────
    h4("Treatment Settings", style = "color:white; padding-left:15px"),

    ## ERA
    checkboxInput("use_era",   "ERA (Bosentan)", value = FALSE),
    conditionalPanel(
      "input.use_era",
      sliderInput("dose_bos", "Bosentan dose (mg BID)",
                  min = 62.5, max = 250, value = 125, step = 62.5)
    ),

    ## PDE5i
    checkboxInput("use_pde5",  "PDE5i (Sildenafil)", value = FALSE),
    conditionalPanel(
      "input.use_pde5",
      sliderInput("dose_sil", "Sildenafil dose (mg TID)",
                  min = 5, max = 80, value = 20, step = 5)
    ),

    ## Prostacyclin
    checkboxInput("use_tre",   "SC Treprostinil", value = FALSE),
    conditionalPanel(
      "input.use_tre",
      sliderInput("dose_tre", "Treprostinil rate (ng/kg/min)",
                  min = 1.25, max = 40, value = 10, step = 1.25)
    ),

    hr(),
    ## ── Simulation Duration ─────────────────────────────────────
    sliderInput("n_weeks", "Simulation duration (weeks)",
                min = 4, max = 104, value = 52, step = 4),

    ## ── Disease Severity ────────────────────────────────────────
    selectInput("severity", "PAH Severity",
                choices = c("Mild (WHO FC II)"    = "mild",
                            "Moderate (WHO FC III)" = "moderate",
                            "Severe (WHO FC IV)"  = "severe"),
                selected = "moderate"),

    br(),
    actionButton("run_sim", "Run Simulation",
                 icon = icon("play"),
                 class = "btn-success btn-block")
  ),

  dashboardBody(
    tags$head(tags$style(HTML("
      .content-wrapper { background-color: #F4F6F9; }
      .box { border-top-color: #2980B9; }
      .info-box-icon { background-color: #2980B9; }
      .kpi-value { font-size: 28px; font-weight: bold; }
      .kpi-unit  { font-size: 14px; color: #7f8c8d; }
    "))),

    tabItems(

      ## ════════════════════════════════════════════════════════════
      ## TAB 1: Disease Overview
      ## ════════════════════════════════════════════════════════════
      tabItem(
        tabName = "overview",
        fluidRow(
          box(width = 12, title = "Pulmonary Arterial Hypertension — Disease Overview",
              status = "primary", solidHeader = TRUE,
              h4("Definition"),
              p("PAH is a severe, progressive vasculopathy characterised by pulmonary vascular
                remodelling, elevated pulmonary vascular resistance (PVR > 240 dyne·s/cm⁵),
                and right ventricular (RV) failure. Haemodynamically defined as mean PAP
                ≥ 25 mmHg with PAWP ≤ 15 mmHg (pre-capillary PH)."),
              h4("Three Core Pathological Pathways"),
              tags$ul(
                tags$li(strong("ET-1 (Endothelin-1): "), "Overproduced in PAH endothelium;
                  activates ETₐ receptors on PA smooth muscle cells (PASMC) → RhoA/ROCK activation
                  → vasoconstriction + PASMC proliferation. Targeted by ERAs."),
                tags$li(strong("NO/cGMP: "), "eNOS activity and NO bioavailability are
                  reduced in PAH (↑ROS → BH4 oxidation → eNOS uncoupling). Downstream sGC/cGMP/PKG
                  axis mediates vasodilation and antiproliferative effects.
                  Targeted by PDE5 inhibitors and sGC stimulators."),
                tags$li(strong("PGI2/cAMP: "), "PGIS is inactivated by ROS; prostacyclin
                  deficiency shifts the balance towards thromboxane A2 (TXA2) → vasoconstriction.
                  Replaced pharmacologically by exogenous prostacyclin or IP receptor agonists.")
              ),
              h4("Key Pathophysiological Processes"),
              p("Vascular remodelling (intimal thickening, medial hypertrophy, plexiform lesions),
                perivascular inflammation (macrophages, T cells, mast cells), in-situ thrombosis,
                oxidative/nitrosative stress, BMPR2 pathway loss-of-function (heritable PAH)."),
              h4("Clinical Course"),
              p("Progressive: 5-year survival ~57% in the modern treatment era. WHO FC IV median
                survival < 6 months without treatment. RV failure is the final common pathway.")
          )
        ),
        fluidRow(
          valueBox(width = 3, value = "≥ 25 mmHg", subtitle = "mPAP threshold for PH",
                   icon = icon("arrow-up"), color = "red"),
          valueBox(width = 3, value = "800 dyne·s/cm⁵", subtitle = "Typical PAH PVR (normal ~80)",
                   icon = icon("filter"), color = "orange"),
          valueBox(width = 3, value = "350 m", subtitle = "Typical 6MWD (normal ~550 m)",
                   icon = icon("walking"), color = "yellow"),
          valueBox(width = 3, value = "1000 pg/mL", subtitle = "Typical NT-proBNP elevation",
                   icon = icon("tint"), color = "purple")
        ),
        fluidRow(
          box(width = 6, title = "Approved Drug Classes", status = "info",
              tableOutput("drug_table")),
          box(width = 6, title = "Model Architecture Summary", status = "info",
              tableOutput("model_arch_table"))
        )
      ),

      ## ════════════════════════════════════════════════════════════
      ## TAB 2: PK Profiles
      ## ════════════════════════════════════════════════════════════
      tabItem(
        tabName = "pk",
        fluidRow(
          box(width = 12, title = "Pharmacokinetic Profiles",
              status = "primary", solidHeader = TRUE,
              tabBox(
                width = 12,
                tabPanel("Sildenafil (PDE5i)",
                         plotlyOutput("pk_sil_plot", height = "420px")),
                tabPanel("Bosentan (ERA)",
                         plotlyOutput("pk_bos_plot", height = "420px")),
                tabPanel("Treprostinil (IP)",
                         plotlyOutput("pk_tre_plot", height = "420px"))
              )
          )
        ),
        fluidRow(
          box(width = 12, title = "PK Parameter Summary",
              status = "info",
              DTOutput("pk_param_table"))
        )
      ),

      ## ════════════════════════════════════════════════════════════
      ## TAB 3: PD Markers
      ## ════════════════════════════════════════════════════════════
      tabItem(
        tabName = "pd",
        fluidRow(
          box(width = 6, title = "cGMP (NO/PDE5 Pathway)",
              status = "success", solidHeader = TRUE,
              plotlyOutput("cgmp_plot", height = "350px")),
          box(width = 6, title = "cAMP (Prostacyclin/IP Pathway)",
              status = "warning", solidHeader = TRUE,
              plotlyOutput("camp_plot", height = "350px"))
        ),
        fluidRow(
          box(width = 6, title = "Plasma ET-1",
              status = "danger", solidHeader = TRUE,
              plotlyOutput("et1_plot", height = "350px")),
          box(width = 6, title = "Vascular Remodeling Index (VRI)",
              status = "primary", solidHeader = TRUE,
              plotlyOutput("vri_plot", height = "350px"))
        )
      ),

      ## ════════════════════════════════════════════════════════════
      ## TAB 4: Hemodynamics
      ## ════════════════════════════════════════════════════════════
      tabItem(
        tabName = "hemo",
        fluidRow(
          valueBoxOutput("pvr_box",  width = 3),
          valueBoxOutput("mpap_box", width = 3),
          valueBoxOutput("co_box",   width = 3),
          valueBoxOutput("coup_box", width = 3)
        ),
        fluidRow(
          box(width = 6, title = "Pulmonary Vascular Resistance (PVR)",
              status = "primary", solidHeader = TRUE,
              plotlyOutput("pvr_plot", height = "350px")),
          box(width = 6, title = "Mean Pulmonary Arterial Pressure (mPAP)",
              status = "danger", solidHeader = TRUE,
              plotlyOutput("mpap_plot", height = "350px"))
        ),
        fluidRow(
          box(width = 6, title = "Cardiac Output (CO)",
              status = "success", solidHeader = TRUE,
              plotlyOutput("co_plot", height = "350px")),
          box(width = 6, title = "RV-PA Coupling (Ees/Ea Ratio)",
              status = "warning", solidHeader = TRUE,
              plotlyOutput("coup_plot", height = "350px"))
        )
      ),

      ## ════════════════════════════════════════════════════════════
      ## TAB 5: Clinical Outcomes
      ## ════════════════════════════════════════════════════════════
      tabItem(
        tabName = "clinical",
        fluidRow(
          valueBoxOutput("bnp_box",  width = 3),
          valueBoxOutput("mwd_box",  width = 3),
          valueBoxOutput("fc_box",   width = 3),
          valueBoxOutput("rvh_box",  width = 3)
        ),
        fluidRow(
          box(width = 6, title = "6-Minute Walk Distance (6MWD)",
              status = "success", solidHeader = TRUE,
              plotlyOutput("walk_plot", height = "350px")),
          box(width = 6, title = "NT-proBNP Biomarker",
              status = "danger", solidHeader = TRUE,
              plotlyOutput("bnp_plot", height = "350px"))
        ),
        fluidRow(
          box(width = 6, title = "WHO Functional Class",
              status = "warning", solidHeader = TRUE,
              plotlyOutput("fc_plot", height = "350px")),
          box(width = 6, title = "RV Hypertrophy Index",
              status = "primary", solidHeader = TRUE,
              plotlyOutput("rvh_plot", height = "350px"))
        )
      ),

      ## ════════════════════════════════════════════════════════════
      ## TAB 6: Treatment Comparison
      ## ════════════════════════════════════════════════════════════
      tabItem(
        tabName = "compare",
        fluidRow(
          box(width = 12, title = "Multi-Arm Comparison",
              status = "primary", solidHeader = TRUE,
              p("Compares: No treatment vs ERA monotherapy vs PDE5i monotherapy
                vs ERA+PDE5i combination vs Triple therapy."),
              fluidRow(
                column(6, plotlyOutput("comp_pvr",  height = "320px")),
                column(6, plotlyOutput("comp_6mwd", height = "320px"))
              ),
              fluidRow(
                column(6, plotlyOutput("comp_bnp",  height = "320px")),
                column(6, plotlyOutput("comp_coup", height = "320px"))
              )
          )
        ),
        fluidRow(
          box(width = 12, title = "Week 12 & Week 52 Summary Table",
              status = "info",
              DTOutput("summary_table"))
        )
      ),

      ## ════════════════════════════════════════════════════════════
      ## TAB 7: Sensitivity Analysis
      ## ════════════════════════════════════════════════════════════
      tabItem(
        tabName = "sensitivity",
        fluidRow(
          box(width = 4, title = "Parameter Selection",
              status = "info", solidHeader = TRUE,
              selectInput("sa_param", "Select parameter to vary:",
                          choices = c(
                            "Sildenafil EC50 (EC50_SIL)" = "EC50_SIL",
                            "Bosentan EC50 (EC50_BOS)"   = "EC50_BOS",
                            "PVR baseline (PVR_BASE)"    = "PVR_BASE",
                            "ET-1 synthesis stim (STIM_ET1)" = "STIM_ET1",
                            "VRI growth rate (KGROW_VRI)" = "KGROW_VRI",
                            "RV Ees baseline (Ees_0)"    = "Ees_0"
                          ), selected = "EC50_SIL"),
              sliderInput("sa_pct", "Range (% of baseline)",
                          min = 10, max = 500, value = 200, step = 10),
              numericInput("sa_n", "Number of levels", value = 5, min = 3, max = 10),
              selectInput("sa_outcome", "Outcome variable:",
                          choices = c("PVR", "mPAP", "WALK6", "NT_proBNP",
                                      "CGMP_PD", "CAMP_PD", "COUP"),
                          selected = "WALK6"),
              actionButton("run_sa", "Run SA", icon = icon("chart-bar"),
                           class = "btn-primary btn-block")
          ),
          box(width = 8, title = "Sensitivity Analysis Results",
              status = "success", solidHeader = TRUE,
              plotlyOutput("sa_plot", height = "460px"))
        ),
        fluidRow(
          box(width = 12, title = "Tornado Plot (Week 52 endpoint)",
              status = "warning",
              plotlyOutput("tornado_plot", height = "360px"))
        )
      ),

      ## ════════════════════════════════════════════════════════════
      ## TAB 8: References
      ## ════════════════════════════════════════════════════════════
      tabItem(
        tabName = "refs",
        fluidRow(
          box(width = 12, title = "Key Literature References",
              status = "primary", solidHeader = TRUE,
              p("See full references.md file in the project directory."),
              DTOutput("refs_table"))
        )
      )
    ) # end tabItems
  ) # end dashboardBody
) # end dashboardPage

## ─────────────────────────────────────────────────────────────────────────────
## SERVER
## ─────────────────────────────────────────────────────────────────────────────
server <- function(input, output, session) {

  ## ── Severity presets ──────────────────────────────────────────
  severity_params <- reactive({
    switch(input$severity,
           "mild"     = list(PVR_BASE = 400, STIM_ET1 = 1.2, VRI_0 = 0.20,
                             WALK_BASE = 450, BNP_BASE = 400),
           "moderate" = list(PVR_BASE = 800, STIM_ET1 = 1.5, VRI_0 = 0.40,
                             WALK_BASE = 350, BNP_BASE = 1000),
           "severe"   = list(PVR_BASE = 1400, STIM_ET1 = 2.0, VRI_0 = 0.65,
                             WALK_BASE = 200, BNP_BASE = 3000)
    )
  })

  ## ── Reactive: run simulation when button pressed ──────────────
  sim_data <- eventReactive(input$run_sim, {
    showNotification("Running simulation...", type = "message", duration = 3)
    sp <- severity_params()
    dose_sil_mg <- if (input$use_pde5) input$dose_sil else 0
    dose_bos_mg <- if (input$use_era)  input$dose_bos else 0
    dose_tre_ng <- if (input$use_tre)  input$dose_tre else 0
    # Convert treprostinil ng/kg/min → mcg/h (assume 70 kg patient)
    dose_tre_mcg_h <- dose_tre_ng * 70 * 60 / 1000

    run_sim(mod_global,
            n_weeks   = input$n_weeks,
            dose_sil  = dose_sil_mg,
            dose_bos  = dose_bos_mg,
            dose_tre  = dose_tre_mcg_h,
            extra_params = sp)
  })

  ## ── Helper: make plotly from ggplot ───────────────────────────
  make_ply <- function(p) ggplotly(p) %>% layout(legend = list(orientation = "h"))

  ## ── KPI boxes (current endpoint at final week) ────────────────
  final_vals <- reactive({
    df <- sim_data(); tail(df, 1)
  })

  output$pvr_box  <- renderValueBox(
    valueBox(round(final_vals()$PVR, 0),
             subtitle = "PVR (dyne·s/cm⁵)", icon = icon("filter"),
             color = if (final_vals()$PVR > 600) "red" else "yellow"))

  output$mpap_box <- renderValueBox(
    valueBox(round(final_vals()$mPAP, 1),
             subtitle = "mPAP (mmHg)", icon = icon("heartbeat"),
             color = if (final_vals()$mPAP > 40) "red" else "yellow"))

  output$co_box   <- renderValueBox(
    valueBox(round(final_vals()$CO_adj, 2),
             subtitle = "CO (L/min)", icon = icon("tint"),
             color = if (final_vals()$CO_adj < 3.5) "red" else "green"))

  output$coup_box <- renderValueBox(
    valueBox(round(final_vals()$COUP, 2),
             subtitle = "Ees/Ea Coupling", icon = icon("compress"),
             color = if (final_vals()$COUP < 0.8) "red" else "green"))

  output$bnp_box  <- renderValueBox(
    valueBox(round(final_vals()$NT_proBNP, 0),
             subtitle = "NT-proBNP (pg/mL)", icon = icon("vial"),
             color = if (final_vals()$NT_proBNP > 1000) "red" else "yellow"))

  output$mwd_box  <- renderValueBox(
    valueBox(round(final_vals()$WALK6, 0),
             subtitle = "6MWD (m)", icon = icon("walking"),
             color = if (final_vals()$WALK6 < 300) "red" else "green"))

  output$fc_box   <- renderValueBox(
    valueBox(paste0("FC ", round(final_vals()$WHO_FC, 0)),
             subtitle = "WHO Functional Class", icon = icon("user-md"),
             color = if (final_vals()$WHO_FC >= 3) "red" else "green"))

  output$rvh_box  <- renderValueBox(
    valueBox(round(final_vals()$RV_HYP, 2),
             subtitle = "RV Hypertrophy Index", icon = icon("heart"),
             color = if (final_vals()$RV_HYP > 0.5) "orange" else "green"))

  ## ── PK plots ──────────────────────────────────────────────────
  pk_short <- function(cmt_depot, conc_col, dose_amt, dose_ii,
                       title, ylab, color) {
    e <- ev(amt = dose_amt, ii = dose_ii, addl = 10, cmt = cmt_depot)
    df <- mod_global %>% mrgsim(ev = e, end = 5 * dose_ii, delta = 0.25) %>%
      as_tibble() %>% mutate(h = time)
    p <- ggplot(df, aes(x = h, y = .data[[conc_col]])) +
      geom_line(size = 1.2, colour = color) +
      labs(title = title, x = "Time (h)", y = ylab) + theme_qsp
    make_ply(p)
  }

  output$pk_sil_plot <- renderPlotly(
    pk_short("DEPOT_SIL", "C_SIL_o", input$dose_sil, 8,
             "Sildenafil PK — First 5 doses", "Concentration (ng/mL)", "#1ABC9C"))

  output$pk_bos_plot <- renderPlotly(
    pk_short("DEPOT_BOS", "C_BOS_o", input$dose_bos, 12,
             "Bosentan PK — First 5 doses", "Concentration (ng/mL)", "#2980B9"))

  output$pk_tre_plot <- renderPlotly({
    dose_tre_mg <- input$dose_tre * 70 * 60 / 1000
    e <- ev(amt = dose_tre_mg * 24, ii = 24, addl = 5, cmt = "SC_TRE")
    df <- mod_global %>% mrgsim(ev = e, end = 6 * 24, delta = 1) %>%
      as_tibble() %>% mutate(h = time)
    p <- ggplot(df, aes(x = h, y = C_TRE_o)) +
      geom_line(size = 1.2, colour = "#8E44AD") +
      labs(title = "Treprostinil PK — SC infusion",
           x = "Time (h)", y = "Concentration (ng/mL)") + theme_qsp
    make_ply(p)
  })

  ## ── PK parameter table ─────────────────────────────────────────
  output$pk_param_table <- renderDT({
    data.frame(
      Drug = c("Sildenafil","Tadalafil","Bosentan","Treprostinil"),
      `t½ (h)` = c(4, 35, 5, 4),
      `Bioavailability` = c("40%","80%","50%","100%"),
      `Route` = c("PO TID","PO QD","PO BID","SC/IV/inh"),
      `Target` = c("PDE5","PDE5","ETₐ/ETB","IP receptor"),
      `EC50` = c("30 ng/mL","5 nM","150 ng/mL","0.5 ng/mL")
    )
  }, options = list(dom = "t"), rownames = FALSE)

  ## ── PD plots ──────────────────────────────────────────────────
  pd_line_plot <- function(col, title, ylab, color, ref_val = NULL) {
    df <- sim_data()
    p  <- ggplot(df, aes(x = Week, y = .data[[col]])) +
      geom_line(size = 1.2, colour = color)
    if (!is.null(ref_val))
      p <- p + geom_hline(yintercept = ref_val, linetype = "dashed", colour = "grey50")
    p <- p + labs(title = title, x = "Week", y = ylab) + theme_qsp
    make_ply(p)
  }

  output$cgmp_plot <- renderPlotly(
    pd_line_plot("CGMP_PD", "cGMP Level", "cGMP (nmol/L)", "#1ABC9C", 1.25))
  output$camp_plot <- renderPlotly(
    pd_line_plot("CAMP_PD", "cAMP Level", "cAMP (nmol/L)", "#8E44AD", 1.33))
  output$et1_plot  <- renderPlotly(
    pd_line_plot("ET1_PD", "Plasma ET-1", "ET-1 (pmol/L)", "#E74C3C", 5.0))
  output$vri_plot  <- renderPlotly(
    pd_line_plot("VRI_PD", "Vascular Remodeling Index", "VRI (0–1)", "#E67E22"))

  ## ── Hemodynamics plots ────────────────────────────────────────
  output$pvr_plot  <- renderPlotly(
    pd_line_plot("PVR",  "PVR", "PVR (dyne·s/cm⁵)", "#2980B9", 80))
  output$mpap_plot <- renderPlotly(
    pd_line_plot("mPAP", "mPAP", "mPAP (mmHg)", "#E74C3C", 20))
  output$co_plot   <- renderPlotly(
    pd_line_plot("CO_adj", "Cardiac Output", "CO (L/min)", "#27AE60", 5.5))
  output$coup_plot <- renderPlotly(
    pd_line_plot("COUP", "RV-PA Coupling", "Ees/Ea", "#9B59B6", 1.0))

  ## ── Clinical outcomes plots ───────────────────────────────────
  output$walk_plot <- renderPlotly(
    pd_line_plot("WALK6", "6-Minute Walk Distance", "6MWD (m)", "#27AE60", 550))
  output$bnp_plot  <- renderPlotly({
    df <- sim_data()
    p  <- ggplot(df, aes(x = Week, y = NT_proBNP)) +
      geom_line(size = 1.2, colour = "#E74C3C") +
      scale_y_log10() +
      geom_hline(yintercept = 300, linetype = "dashed", colour = "grey50") +
      labs(title = "NT-proBNP", x = "Week", y = "NT-proBNP (pg/mL, log10)") +
      theme_qsp
    make_ply(p)
  })
  output$fc_plot   <- renderPlotly(
    pd_line_plot("WHO_FC", "WHO FC (1–4)", "Functional Class", "#F39C12"))
  output$rvh_plot  <- renderPlotly(
    pd_line_plot("RV_HYP", "RV Hypertrophy Index", "Index (0–1)", "#9B59B6"))

  ## ── Multi-arm comparison ──────────────────────────────────────
  multi_sim <- eventReactive(input$run_sim, {
    sp <- severity_params()
    nw <- input$n_weeks
    arms <- list(
      "No treatment"        = list(0, 0, 0),
      "ERA mono"            = list(0, 125, 0),
      "PDE5i mono"          = list(20, 0, 0),
      "ERA + PDE5i"         = list(20, 125, 0),
      "Triple therapy"      = list(20, 125, 5 * 70 * 60 / 1000)
    )
    purrr::map2_dfr(names(arms), arms, function(nm, doses) {
      run_sim(mod_global,
              n_weeks = nw,
              dose_sil = doses[[1]],
              dose_bos = doses[[2]],
              dose_tre = doses[[3]],
              extra_params = sp) %>%
        mutate(Group = nm)
    }) %>%
      mutate(Group = factor(Group, levels = names(arms)))
  })

  color_pal <- c(
    "No treatment"   = "#C0392B",
    "ERA mono"       = "#2980B9",
    "PDE5i mono"     = "#1ABC9C",
    "ERA + PDE5i"    = "#8E44AD",
    "Triple therapy" = "#27AE60"
  )

  multi_plot <- function(col, ylab, ref = NULL) {
    df <- multi_sim()
    p  <- ggplot(df, aes(x = Week, y = .data[[col]], colour = Group)) +
      geom_line(size = 1) +
      scale_colour_manual(values = color_pal) +
      labs(x = "Week", y = ylab) + theme_qsp
    if (!is.null(ref))
      p <- p + geom_hline(yintercept = ref, linetype = "dashed", colour = "grey50")
    make_ply(p)
  }

  output$comp_pvr  <- renderPlotly(multi_plot("PVR",  "PVR (dyne·s/cm⁵)", 80))
  output$comp_6mwd <- renderPlotly(multi_plot("WALK6","6MWD (m)", 550))
  output$comp_bnp  <- renderPlotly({
    df <- multi_sim()
    p  <- ggplot(df, aes(x = Week, y = NT_proBNP, colour = Group)) +
      geom_line(size = 1) + scale_y_log10() +
      scale_colour_manual(values = color_pal) +
      geom_hline(yintercept = 300, linetype = "dashed", colour = "grey50") +
      labs(x = "Week", y = "NT-proBNP (log10, pg/mL)") + theme_qsp
    make_ply(p)
  })
  output$comp_coup <- renderPlotly(multi_plot("COUP", "Ees/Ea Ratio", 1.0))

  output$summary_table <- renderDT({
    df <- multi_sim()
    make_tab <- function(df, wk) {
      df %>%
        filter(abs(Week - wk) < 0.6) %>%
        group_by(Group) %>% slice(1) %>%
        select(Group, PVR, mPAP, WALK6, NT_proBNP, COUP, WHO_FC) %>%
        mutate(Week = wk, across(where(is.numeric), ~round(.x, 1)))
    }
    bind_rows(make_tab(df, 12), make_tab(df, 52)) %>%
      arrange(Week, Group)
  }, options = list(pageLength = 10), rownames = FALSE)

  ## ── Sensitivity analysis ──────────────────────────────────────
  sa_result <- eventReactive(input$run_sa, {
    sp   <- severity_params()
    base <- param(mod_global)[[input$sa_param]]
    mults <- seq(0.1, input$sa_pct / 100, length.out = input$sa_n)
    vals  <- base * mults

    purrr::map_dfr(vals, function(v) {
      extra <- c(sp, setNames(list(v), input$sa_param))
      run_sim(mod_global,
              n_weeks  = input$n_weeks,
              dose_sil = if (input$use_pde5) input$dose_sil else 0,
              dose_bos = if (input$use_era) input$dose_bos else 0,
              dose_tre = if (input$use_tre) input$dose_tre * 70 * 60 / 1000 else 0,
              extra_params = extra) %>%
        mutate(ParamVal = v, ParamLabel = paste0(input$sa_param, " = ", round(v, 2)))
    })
  })

  output$sa_plot <- renderPlotly({
    df <- sa_result()
    p  <- ggplot(df, aes(x = Week, y = .data[[input$sa_outcome]],
                         colour = factor(round(ParamVal, 2)))) +
      geom_line(size = 1) +
      scale_colour_viridis_d(option = "plasma", name = input$sa_param) +
      labs(title = paste("SA:", input$sa_param, "→", input$sa_outcome),
           x = "Week", y = input$sa_outcome) + theme_qsp
    make_ply(p)
  })

  output$tornado_plot <- renderPlotly({
    sp <- severity_params()
    params_to_vary <- c("EC50_SIL","EC50_BOS","PVR_BASE","STIM_ET1","Ees_0","KGROW_VRI")
    base_run <- run_sim(mod_global, n_weeks = input$n_weeks,
                        dose_sil = if (input$use_pde5) input$dose_sil else 0,
                        dose_bos = if (input$use_era) input$dose_bos else 0,
                        dose_tre = if (input$use_tre) input$dose_tre * 70 * 60 / 1000 else 0,
                        extra_params = sp)
    base_val <- tail(base_run$WALK6, 1)

    tornado_df <- purrr::map_dfr(params_to_vary, function(p_name) {
      base_p <- param(mod_global)[[p_name]]
      run_lo <- run_sim(mod_global, n_weeks = input$n_weeks,
                        dose_sil = if (input$use_pde5) input$dose_sil else 0,
                        dose_bos = if (input$use_era) input$dose_bos else 0,
                        dose_tre = if (input$use_tre) input$dose_tre * 70 * 60 / 1000 else 0,
                        extra_params = c(sp, setNames(list(base_p * 0.5), p_name)))
      run_hi <- run_sim(mod_global, n_weeks = input$n_weeks,
                        dose_sil = if (input$use_pde5) input$dose_sil else 0,
                        dose_bos = if (input$use_era) input$dose_bos else 0,
                        dose_tre = if (input$use_tre) input$dose_tre * 70 * 60 / 1000 else 0,
                        extra_params = c(sp, setNames(list(base_p * 2.0), p_name)))
      data.frame(
        Param = p_name,
        Lo = tail(run_lo$WALK6, 1) - base_val,
        Hi = tail(run_hi$WALK6, 1) - base_val
      )
    }) %>%
      mutate(Range = abs(Hi - Lo)) %>%
      arrange(Range) %>%
      mutate(Param = factor(Param, levels = Param))

    p <- ggplot(tornado_df) +
      geom_segment(aes(x = Lo, xend = Hi, y = Param, yend = Param),
                   size = 8, colour = "#2980B9", alpha = 0.7) +
      geom_vline(xintercept = 0, colour = "red", linetype = "dashed") +
      labs(title = "Tornado Plot — Impact on 6MWD at Final Week",
           x = "Change in 6MWD from baseline (m)", y = "Parameter") +
      theme_qsp
    make_ply(p)
  })

  ## ── Overview tables ───────────────────────────────────────────
  output$drug_table <- renderTable({
    data.frame(
      `Drug Class` = c("PDE5 Inhibitors","PDE5 Inhibitors",
                       "sGC Stimulator",
                       "ERAs","ERAs","ERAs",
                       "Prostacyclin","Prostacyclin","Prostacyclin",
                       "IP Agonist"),
      Drug = c("Sildenafil","Tadalafil","Riociguat",
               "Bosentan","Ambrisentan","Macitentan",
               "Epoprostenol","Treprostinil","Iloprost","Selexipag"),
      Route = c("PO","PO","PO","PO","PO","PO","IV","SC/IV/inh","Inh","PO"),
      Target = c("PDE5","PDE5","sGC","ETₐ/ETB","ETₐ","ETₐ/ETB",
                 "IP","IP","IP","IP"),
      stringsAsFactors = FALSE
    )
  })

  output$model_arch_table <- renderTable({
    data.frame(
      Module = c("PK (Sildenafil)","PK (Bosentan)","PK (Treprostinil)",
                 "ET-1 dynamics","cGMP pathway","cAMP pathway",
                 "Vascular tone","Vascular remodeling",
                 "Hemodynamics","RV mechanics","Clinical outputs"),
      `Compartments` = c(3,3,2,1,1,1,0,1,0,1,0),
      `ODEs` = c(3,3,2,1,1,1,"algebraic",1,"algebraic",1,"algebraic"),
      stringsAsFactors = FALSE
    )
  })

  ## ── References table ─────────────────────────────────────────
  output$refs_table <- renderDT({
    data.frame(
      `#` = 1:12,
      Topic = c("Disease overview","ET-1 pathway","ET-1 pathway","ERA trial",
                "ERA trial","NO/sGC","sGC stimulator","PDE5i trial",
                "Prostacyclin","IP agonist trial","BMPR2","Combo therapy"),
      Reference = c(
        "Humbert M et al. Eur Respir J 2019;53:1801887",
        "Channick RN et al. Lancet 2001;358:1119",
        "Pulido T et al. NEJM 2013;369:809 (SERAPHIN)",
        "Galiè N et al. Circulation 2008;117:3010 (ARIES)",
        "Ghofrani HA et al. NEJM 2013;369:330 (PATENT-1)",
        "Stasch JP et al. J Clin Invest 2006;116:2552",
        "Galiè N et al. NEJM 2005;353:2148 (SUPER-1)",
        "Galiè N et al. Circulation 2009;119:2894 (PHIRST)",
        "Barst RJ et al. NEJM 1996;334:296",
        "Sitbon O et al. NEJM 2015;373:2522 (GRIPHON)",
        "Lane KB et al. Nat Genet 2000;26:81",
        "Galiè N et al. NEJM 2015;373:834 (AMBITION)"
      ),
      stringsAsFactors = FALSE
    )
  }, options = list(pageLength = 12), rownames = FALSE)
}

## ─────────────────────────────────────────────────────────────────────────────
## Launch
## ─────────────────────────────────────────────────────────────────────────────
shinyApp(ui = ui, server = server)
