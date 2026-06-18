## =============================================================================
## Metabolic Syndrome QSP — Interactive Shiny Dashboard
## 6 Tabs: Patient Profile · PK Profiles · PD Biomarkers ·
##         Clinical Endpoints · Scenario Comparison · Sensitivity Analysis
##
## Required: shiny, shinydashboard, mrgsolve, dplyr, ggplot2, plotly,
##           DT, tidyr, purrr, patchwork
##
## Author: CCR — Claude Code Routine (2026-06-18)
## =============================================================================

library(shiny)
library(shinydashboard)
library(mrgsolve)
library(dplyr)
library(ggplot2)
library(plotly)
library(DT)
library(tidyr)
library(purrr)

## ── Model code (same as mrgsolve model file) ──────────────────────────────────
ms_model_code <- '
$PROB Metabolic Syndrome QSP — Interactive Shiny Dashboard Model

$PARAM
BW=90 VAT0=5.0 SAT0=20.0 G_b=100 I_b=25
Sg=0.021 Si=4.5e-4 k_HGP=0.18 k_gut=0.07
GFR=120 TmG0=375 Kg1=0.5 Kg2=0.02 tau_B=90
Beta0=1.0 k_BetaDeath=0.0003 Gluc0=80 kGluc=0.08
GLP1_e0=5.0 kGLP1=0.12
VLDL0=120 LDL0=140 HDL0=38 TG0=220
kVLDL=0.05 kLDL=0.03 kHDL=0.015 kCETP=0.02
TNFa0=15.0 IL6_0=3.0 IL1b0=2.0 CRP0=3.0
kTNF=0.06 kIL6=0.04 kCRP=0.10
AngII0=0.2 MAP0=100 SVR0=1200 CO0=5.0 kRAAS=0.03
Lep0=25 Adip0=5.0 kLep=0.8 kAdip=0.3
AMPK0=1.0 kAMPK=0.15
MET_Fg=0.55 MET_ka=0.6 MET_CL=28 MET_Vc=380 MET_Q=60 MET_Vt=2000
GLP_ka=0.004 GLP_CL=0.056 GLP_Vc=8.0
SGi_ka=1.5 SGi_CL=12.0 SGI_Vc=73.0 SGi_F=0.86
ST_ka=0.4 ST_CL=25 ST_Vc=134 ST_F=0.20
ARB_ka=1.0 ARB_CL=75 ARB_Vc=34 ARB_F=0.33
MET_EC50=1.5 MET_Emax=0.30
GLP_EC50=0.8 GLP_Emax=2.0 GLP_WT50=2.0 GLP_WTmax=0.12
SGi_EC50=15 SGi_Emax=0.90
ST_EC50=0.05 ST_Emax=0.70
ARB_EC50=0.3 ARB_Emax=0.85
ADIPON0=5.0 VAT0_ref=5.0 Adip0_ref=5.0 Lep0_ref=25.0

$CMT GGUT GPLAS BETA IPLAS GLUCPLAS GLP1E
     VLDLC LDLC HDLC TRIGLY
     VAT SAT LEP ADIPON TNFA IL6C IL1BC CRPC
     ANGII MAPC AMPKC
     MET_GUT MET_CEN MET_PER
     GLP_SC GLP_CEN SGI_GUT SGI_CEN STA_GUT STA_CEN ARB_GUT ARB_CEN

$INIT GGUT=0 GPLAS=100 BETA=1.0 IPLAS=25 GLUCPLAS=80 GLP1E=5.0
      VLDLC=120 LDLC=140 HDLC=38 TRIGLY=220
      VAT=5.0 SAT=20.0 LEP=25 ADIPON=5.0
      TNFA=15 IL6C=3.0 IL1BC=2.0 CRPC=3.0
      ANGII=0.2 MAPC=100 AMPKC=1.0
      MET_GUT=0 MET_CEN=0 MET_PER=0
      GLP_SC=0 GLP_CEN=0 SGI_GUT=0 SGI_CEN=0 STA_GUT=0 STA_CEN=0 ARB_GUT=0 ARB_CEN=0

$ODE
double Cp_MET  = MET_CEN / MET_Vc;
double Cp_GLP  = GLP_CEN / GLP_Vc;
double Cp_SGI  = SGI_CEN / SGI_Vc;
double Cp_STA  = STA_CEN / ST_Vc;
double Cp_ARB  = ARB_CEN / ARB_Vc;
double E_MET_HGP  = MET_Emax * Cp_MET  / (MET_EC50 + Cp_MET);
double E_GLP_INS  = GLP_Emax * Cp_GLP  / (GLP_EC50 + Cp_GLP);
double E_GLP_WT   = GLP_WTmax* Cp_GLP  / (GLP_WT50 + Cp_GLP);
double E_SGI_TMG  = SGi_Emax * Cp_SGI  / (SGi_EC50 + Cp_SGI);
double E_STA_CHOL = ST_Emax  * Cp_STA  / (ST_EC50  + Cp_STA);
double E_ARB_BP   = ARB_Emax * Cp_ARB  / (ARB_EC50 + Cp_ARB);
double IR_FFA  = 1.0 + 0.15 * (VAT / VAT0 - 1);
double IR_TNF  = 1.0 + 0.08 * (TNFA / TNFa0 - 1);
double IR_AngII = 1.0 + 0.04 * (ANGII / AngII0 - 1);
double IR_total = IR_FFA * IR_TNF * IR_AngII;
double f_Ins_HGP = 1.0 / (1.0 + 0.03 * IPLAS);
double HGP = k_HGP * BW * f_Ins_HGP * IR_total * (1.0 - E_MET_HGP)
           * (1.0 + 0.2 * (GLUCPLAS / Gluc0 - 1));
double Rd  = (Sg + Si * IPLAS / IR_total) * GPLAS;
double TmG_eff = TmG0 * (1.0 - E_SGI_TMG);
double UGE  = (GFR * GPLAS / 100.0 > TmG_eff) ? (GFR * GPLAS / 100.0 - TmG_eff) : 0.0;
double Ra_gut = k_gut * GGUT;
dxdt_GGUT  = -k_gut * GGUT;
dxdt_GPLAS = (HGP + Ra_gut - Rd - UGE) / BW;
double GSIS_first  = Kg1 * (GPLAS - G_b);
double GSIS_second = Kg2 * GPLAS * BETA;
double Ins_sec = (GSIS_first + GSIS_second > 0) ? (GSIS_first + GSIS_second) * (1.0 + E_GLP_INS) : 0.0;
double k_BetaGrowth = 0.0002 * (GPLAS / G_b);
double k_BetaLoss   = k_BetaDeath * IL1BC / IL1b0;
dxdt_BETA  = (k_BetaGrowth - k_BetaLoss) * BETA;
dxdt_IPLAS = Ins_sec - 0.05 * IPLAS;
double Gluc_target = Gluc0 / (1.0 + 0.015 * IPLAS);
dxdt_GLUCPLAS = kGluc * (Gluc_target - GLUCPLAS);
dxdt_GLP1E = kGLP1 * (GPLAS / G_b) - 0.15 * GLP1E;
double VLDL_prod = kVLDL * (VLDLC + IR_total * 15) * (ADIPON0 / ADIPON) * (1.0 - E_MET_HGP * 0.3);
dxdt_VLDLC = VLDL_prod - 0.06 * VLDLC;
double LDL_elim = (0.04 + E_STA_CHOL * 0.08) * LDLC;
dxdt_LDLC = kLDL * TRIGLY - LDL_elim;
double HDL_prod = 0.015 * ADIPON;
dxdt_HDLC = HDL_prod - (kHDL * HDLC + kCETP * TRIGLY / TG0);
double TG_prod = 0.04 * VLDLC;
dxdt_TRIGLY = TG_prod - 0.025 * TRIGLY * (1.0 + 0.1 * ADIPON / Adip0) * (1.0 + 0.1 * E_GLP_WT);
double VAT_accum = 0.0005 * (IPLAS / I_b - 1) * VAT;
double VAT_loss  = E_GLP_WT * 0.005 * VAT + 0.002 * VAT;
dxdt_VAT = VAT_accum - VAT_loss;
dxdt_SAT = -0.0003 * (1.0 - E_GLP_WT) * SAT;
dxdt_LEP = 0.05 * (kLep * (VAT + SAT) / (VAT0 + SAT0) * Lep0 - LEP);
dxdt_ADIPON = 0.03 * (Adip0 * (VAT0 / VAT) * (1.0 + 0.5 * E_ARB_BP) - ADIPON);
double TNF_prod = kTNF * (VAT / VAT0) * (LEP / Lep0);
dxdt_TNFA = TNF_prod - 0.05 * TNFA - 0.02 * (ADIPON / Adip0) * TNFA;
dxdt_IL6C = kIL6 * (VAT / VAT0) * (TNFA / TNFa0) - 0.08 * IL6C;
dxdt_IL1BC = 0.015 * (VAT / VAT0) * (TNFA / TNFa0) - 0.06 * IL1BC;
dxdt_CRPC = kCRP * IL6C / IL6_0 - 0.04 * CRPC;
double AngII_prod = kRAAS * (VAT / VAT0) * (TNFA / TNFa0);
dxdt_ANGII = AngII_prod - 0.15 * ANGII * (1.0 - E_ARB_BP * 0.7);
double MAP_target = MAP0 + 15*(ANGII/AngII0-1) + 8*(IL6C/IL6_0-1) - 10*E_ARB_BP - 4*E_SGI_TMG;
dxdt_MAPC = 0.01 * (MAP_target - MAPC);
double AMPK_target = AMPK0 * (1.0 + 2.0*E_MET_HGP) * (1.0 + 0.5*ADIPON/Adip0) / IR_total;
dxdt_AMPKC = 0.05 * (AMPK_target - AMPKC);
dxdt_MET_GUT = -MET_ka * MET_GUT;
dxdt_MET_CEN =  MET_ka*MET_Fg*MET_GUT - (MET_CL/MET_Vc)*MET_CEN - (MET_Q/MET_Vc)*MET_CEN + (MET_Q/MET_Vt)*MET_PER;
dxdt_MET_PER =  (MET_Q/MET_Vc)*MET_CEN - (MET_Q/MET_Vt)*MET_PER;
dxdt_GLP_SC  = -GLP_ka * GLP_SC;
dxdt_GLP_CEN =  GLP_ka*GLP_SC - (GLP_CL/GLP_Vc)*GLP_CEN;
dxdt_SGI_GUT = -SGi_ka * SGI_GUT;
dxdt_SGI_CEN =  SGi_ka*SGi_F*SGI_GUT - (SGi_CL/SGI_Vc)*SGI_CEN;
dxdt_STA_GUT = -ST_ka * STA_GUT;
dxdt_STA_CEN =  ST_ka*ST_F*STA_GUT - (ST_CL/ST_Vc)*STA_CEN;
dxdt_ARB_GUT = -ARB_ka * ARB_GUT;
dxdt_ARB_CEN =  ARB_ka*ARB_F*ARB_GUT - (ARB_CL/ARB_Vc)*ARB_CEN;

$TABLE
double HbA1c      = 5.0 + GPLAS / 30.0;
double HOMA_IR    = GPLAS * IPLAS / 405.0;
double FPG        = GPLAS;
double TG_HDL     = TRIGLY / HDLC;
double SBP        = MAPC + 40;
double DBP        = MAPC - 10;
double MetS_Zscore = (GPLAS-100)/20 + (MAPC-93)/13 + (TRIGLY-150)/50 - (HDLC-50)/12 + (VAT-4.5)/1.0;
double Cp_MET_out = MET_CEN / MET_Vc;
double Cp_GLP_out = GLP_CEN / GLP_Vc;
double Cp_SGI_out = SGI_CEN / SGI_Vc;
double Cp_STA_out = STA_CEN / ST_Vc;
double Cp_ARB_out = ARB_CEN / ARB_Vc;

$CAPTURE GPLAS IPLAS GLUCPLAS GLP1E
         VLDLC LDLC HDLC TRIGLY
         VAT SAT LEP ADIPON TNFA IL6C IL1BC CRPC ANGII MAPC AMPKC BETA
         HbA1c HOMA_IR FPG TG_HDL SBP DBP MetS_Zscore
         Cp_MET_out Cp_GLP_out Cp_SGI_out Cp_STA_out Cp_ARB_out
'

ms_mod <- mcode("MetabolicSyndrome_Shiny", ms_model_code)

## ── Helper: run_scenario ─────────────────────────────────────────────────────
run_scenario <- function(scenario_inputs, sim_weeks, bw, vat0, sat0) {
  mod <- ms_mod %>%
    param(BW = bw, VAT0 = vat0, SAT0 = sat0) %>%
    init(VAT = vat0, SAT = sat0,
         LEP = 0.8 * (vat0 + sat0) / 25 * 25)

  ev_list <- list()

  if (scenario_inputs$use_met) {
    ev_list[["met"]] <- ev(ID = 1,
                           amt = scenario_inputs$met_dose,
                           cmt = "MET_GUT",
                           ii = 12,
                           addl = 2 * 7 * sim_weeks - 1)
  }
  if (scenario_inputs$use_glp) {
    ev_list[["glp"]] <- ev(ID = 1,
                           amt = scenario_inputs$glp_dose,
                           cmt = "GLP_SC",
                           ii = 168,
                           addl = sim_weeks - 1)
  }
  if (scenario_inputs$use_sgi) {
    ev_list[["sgi"]] <- ev(ID = 1,
                           amt = scenario_inputs$sgi_dose,
                           cmt = "SGI_GUT",
                           ii = 24,
                           addl = 7 * sim_weeks - 1)
  }
  if (scenario_inputs$use_sta) {
    ev_list[["sta"]] <- ev(ID = 1,
                           amt = scenario_inputs$sta_dose,
                           cmt = "STA_GUT",
                           ii = 24,
                           addl = 7 * sim_weeks - 1)
  }
  if (scenario_inputs$use_arb) {
    ev_list[["arb"]] <- ev(ID = 1,
                           amt = scenario_inputs$arb_dose,
                           cmt = "ARB_GUT",
                           ii = 24,
                           addl = 7 * sim_weeks - 1)
  }

  ev_combined <- if (length(ev_list) > 0) {
    do.call(c, ev_list)
  } else {
    ev(ID = 1, time = 0, amt = 0, cmt = 1)
  }

  out <- mod %>%
    ev(ev_combined) %>%
    mrgsim(end = sim_weeks * 168, delta = 24) %>%
    as_tibble() %>%
    mutate(Week = time / 168)

  return(out)
}

## ─────────────────────────────────────────────────────────────────────────────
## UI
## ─────────────────────────────────────────────────────────────────────────────
ui <- dashboardPage(
  skin = "purple",

  dashboardHeader(
    title = "MetS QSP Dashboard",
    titleWidth = 280
  ),

  dashboardSidebar(
    width = 280,
    sidebarMenu(
      menuItem("Patient Profile",       tabName = "patient",   icon = icon("user-md")),
      menuItem("PK Profiles",           tabName = "pk",        icon = icon("pills")),
      menuItem("PD Biomarkers",         tabName = "pd",        icon = icon("heartbeat")),
      menuItem("Clinical Endpoints",    tabName = "endpoints", icon = icon("chart-line")),
      menuItem("Scenario Comparison",   tabName = "scenarios", icon = icon("balance-scale")),
      menuItem("Sensitivity Analysis",  tabName = "sensitivity", icon = icon("sliders-h"))
    ),
    hr(),
    h5("Patient Parameters", style = "color:white; padding-left:15px;"),
    sliderInput("bw",   "Body Weight (kg)",     min = 60, max = 140, value = 90),
    sliderInput("vat0", "Visceral Fat (kg)",    min = 2,  max = 12,  value = 5, step = 0.5),
    sliderInput("sat0", "Subcutaneous Fat (kg)",min = 10, max = 35,  value = 20),
    sliderInput("sim_weeks", "Simulation (weeks)", min = 12, max = 104, value = 52, step = 4),
    hr(),
    h5("Treatment Selection", style = "color:white; padding-left:15px;"),
    checkboxInput("use_met", "Metformin",       value = FALSE),
    conditionalPanel(condition = "input.use_met",
      numericInput("met_dose", "Dose (mg BID)", value = 1000, min = 250, max = 2000, step = 250)),
    checkboxInput("use_glp", "GLP-1 RA (weekly SC)", value = FALSE),
    conditionalPanel(condition = "input.use_glp",
      numericInput("glp_dose", "Dose (mg)", value = 1.0, min = 0.25, max = 2.4, step = 0.25)),
    checkboxInput("use_sgi", "SGLT2 Inhibitor (QD)", value = FALSE),
    conditionalPanel(condition = "input.use_sgi",
      numericInput("sgi_dose", "Dose (mg)", value = 10, min = 5, max = 25, step = 5)),
    checkboxInput("use_sta", "Statin (QD)",     value = FALSE),
    conditionalPanel(condition = "input.use_sta",
      numericInput("sta_dose", "Dose (mg)", value = 10, min = 5, max = 40, step = 5)),
    checkboxInput("use_arb", "ARB / ACEi (QD)", value = FALSE),
    conditionalPanel(condition = "input.use_arb",
      numericInput("arb_dose", "Dose (mg)", value = 50, min = 25, max = 100, step = 25)),
    actionButton("run_sim", "Run Simulation", class = "btn-success btn-block",
                 style = "margin:10px; width:250px;")
  ),

  dashboardBody(
    tags$head(tags$style(HTML("
      .content-wrapper { background-color: #f4f4f4; }
      .box { border-top: 3px solid #605ca8; }
    "))),

    tabItems(

      ## ── Tab 1: Patient Profile ───────────────────────────────────────────
      tabItem(tabName = "patient",
        fluidRow(
          valueBoxOutput("vbox_bmi",    width = 3),
          valueBoxOutput("vbox_mets",   width = 3),
          valueBoxOutput("vbox_homa",   width = 3),
          valueBoxOutput("vbox_hba1c",  width = 3)
        ),
        fluidRow(
          box(title = "Metabolic Syndrome NCEP-ATP III Criteria",
              width = 6, status = "primary",
              tableOutput("ncep_table")),
          box(title = "Baseline Biomarker Profile",
              width = 6, status = "info",
              plotlyOutput("radar_plot", height = "380px"))
        ),
        fluidRow(
          box(title = "Disease Pathophysiology Summary",
              width = 12, status = "warning",
              HTML('<p><b>Metabolic Syndrome</b> is defined by ≥3 of 5 NCEP-ATP III criteria:
                   central obesity, hypertriglyceridemia, low HDL-C, hypertension, and
                   impaired fasting glucose. The underlying pathophysiology centers on
                   <b>visceral adipose tissue (VAT)</b> expansion, leading to:</p>
                   <ul>
                     <li><b>Insulin resistance</b> via FFA flux and adipokine imbalance (↑Resistin, ↓Adiponectin)</li>
                     <li><b>Chronic low-grade inflammation</b> via M1 macrophage infiltration and TNF-α, IL-6 secretion</li>
                     <li><b>Dyslipidemia</b> via ↑VLDL secretion, CETP-mediated HDL reduction, ↑small dense LDL</li>
                     <li><b>Hypertension</b> via RAAS activation, endothelial dysfunction (↓NO), and ↑sympathetic tone</li>
                     <li><b>β-cell exhaustion</b> via glucotoxicity and IL-1β-mediated apoptosis → T2DM risk</li>
                   </ul>
                   <p>This QSP model integrates 22 ODEs covering glucose homeostasis, lipid metabolism,
                   adipose biology, inflammation, RAAS, and 5 drug mechanisms.</p>'))
        )
      ),

      ## ── Tab 2: PK Profiles ───────────────────────────────────────────────
      tabItem(tabName = "pk",
        fluidRow(
          box(title = "Drug Plasma Concentration — Metformin",
              width = 6, status = "primary",
              plotlyOutput("pk_met", height = "300px")),
          box(title = "Drug Plasma Concentration — GLP-1 RA",
              width = 6, status = "success",
              plotlyOutput("pk_glp", height = "300px"))
        ),
        fluidRow(
          box(title = "Drug Plasma Concentration — SGLT2 Inhibitor",
              width = 4, status = "info",
              plotlyOutput("pk_sgi", height = "300px")),
          box(title = "Drug Plasma Concentration — Statin",
              width = 4, status = "warning",
              plotlyOutput("pk_sta", height = "300px")),
          box(title = "Drug Plasma Concentration — ARB/ACEi",
              width = 4, status = "danger",
              plotlyOutput("pk_arb", height = "300px"))
        ),
        fluidRow(
          box(title = "PK Parameter Summary", width = 12, status = "primary",
              DTOutput("pk_table"))
        )
      ),

      ## ── Tab 3: PD Biomarkers ─────────────────────────────────────────────
      tabItem(tabName = "pd",
        fluidRow(
          box(title = "Glucose Dynamics",
              width = 6, status = "primary",
              plotlyOutput("pd_glucose", height = "300px")),
          box(title = "Insulin & β-cell Function",
              width = 6, status = "success",
              plotlyOutput("pd_insulin", height = "300px"))
        ),
        fluidRow(
          box(title = "Adipokines (Leptin & Adiponectin)",
              width = 6, status = "warning",
              plotlyOutput("pd_adipokines", height = "300px")),
          box(title = "Inflammatory Cytokines",
              width = 6, status = "danger",
              plotlyOutput("pd_inflam", height = "300px"))
        ),
        fluidRow(
          box(title = "RAAS & Blood Pressure",
              width = 6, status = "info",
              plotlyOutput("pd_bp", height = "300px")),
          box(title = "AMPK & Visceral Fat",
              width = 6, status = "primary",
              plotlyOutput("pd_ampk_vat", height = "300px"))
        )
      ),

      ## ── Tab 4: Clinical Endpoints ─────────────────────────────────────────
      tabItem(tabName = "endpoints",
        fluidRow(
          box(title = "HbA1c Trajectory",
              width = 6, status = "primary",
              plotlyOutput("ep_hba1c", height = "300px")),
          box(title = "Lipid Panel Over Time",
              width = 6, status = "warning",
              plotlyOutput("ep_lipids", height = "300px"))
        ),
        fluidRow(
          box(title = "Blood Pressure Response",
              width = 6, status = "danger",
              plotlyOutput("ep_bp", height = "300px")),
          box(title = "MetS Z-score (Composite Severity)",
              width = 6, status = "success",
              plotlyOutput("ep_metsz", height = "300px"))
        ),
        fluidRow(
          box(title = "Key Endpoint Summary at End of Simulation",
              width = 12, status = "primary",
              DTOutput("ep_summary_table"))
        )
      ),

      ## ── Tab 5: Scenario Comparison ────────────────────────────────────────
      tabItem(tabName = "scenarios",
        fluidRow(
          box(title = "5-Scenario Comparison: HbA1c",
              width = 6, status = "primary",
              plotlyOutput("sc_hba1c", height = "300px")),
          box(title = "5-Scenario Comparison: LDL-C",
              width = 6, status = "warning",
              plotlyOutput("sc_ldl", height = "300px"))
        ),
        fluidRow(
          box(title = "5-Scenario Comparison: SBP",
              width = 6, status = "danger",
              plotlyOutput("sc_sbp", height = "300px")),
          box(title = "5-Scenario Comparison: VAT",
              width = 6, status = "success",
              plotlyOutput("sc_vat", height = "300px"))
        ),
        fluidRow(
          box(title = "Week-52 Scenario Summary Table",
              width = 12, status = "primary",
              DTOutput("sc_summary_table"))
        )
      ),

      ## ── Tab 6: Sensitivity Analysis ───────────────────────────────────────
      tabItem(tabName = "sensitivity",
        fluidRow(
          box(title = "BMI vs HbA1c at End of Simulation",
              width = 6, status = "primary",
              plotlyOutput("sa_bmi_hba1c", height = "350px")),
          box(title = "BMI vs HOMA-IR at End of Simulation",
              width = 6, status = "warning",
              plotlyOutput("sa_bmi_homa", height = "350px"))
        ),
        fluidRow(
          box(title = "Drug Dose-Response at Week 52",
              width = 12, status = "info",
              selectInput("sa_drug", "Select Drug",
                          choices = c("Metformin", "GLP-1 RA", "SGLT2i", "Statin", "ARB"),
                          selected = "Metformin"),
              selectInput("sa_endpoint", "Select Endpoint",
                          choices = c("HbA1c", "LDLC", "SBP", "VAT", "MetS_Zscore"),
                          selected = "HbA1c"),
              plotlyOutput("sa_dose_response", height = "350px"))
        )
      )
    )
  )
)

## ─────────────────────────────────────────────────────────────────────────────
## SERVER
## ─────────────────────────────────────────────────────────────────────────────
server <- function(input, output, session) {

  ## ── Reactive simulation result ────────────────────────────────────────────
  sim_result <- eventReactive(input$run_sim, {
    withProgress(message = "Running QSP simulation...", value = 0, {
      incProgress(0.3, detail = "Compiling events...")
      inputs <- list(
        use_met   = input$use_met,
        met_dose  = if (input$use_met) input$met_dose else 0,
        use_glp   = input$use_glp,
        glp_dose  = if (input$use_glp) input$glp_dose else 0,
        use_sgi   = input$use_sgi,
        sgi_dose  = if (input$use_sgi) input$sgi_dose else 0,
        use_sta   = input$use_sta,
        sta_dose  = if (input$use_sta) input$sta_dose else 0,
        use_arb   = input$use_arb,
        arb_dose  = if (input$use_arb) input$arb_dose else 0
      )
      incProgress(0.5, detail = "Solving ODEs...")
      out <- run_scenario(inputs, input$sim_weeks, input$bw, input$vat0, input$sat0)
      incProgress(0.2, detail = "Done.")
      out
    })
  }, ignoreNULL = FALSE)

  ## Initial run on app load
  observe({ isolate({ input$run_sim }) })

  ## ── Tab 1: ValueBoxes ─────────────────────────────────────────────────────
  output$vbox_bmi <- renderValueBox({
    bmi <- input$bw / (1.70^2)
    valueBox(round(bmi, 1), "BMI (kg/m²)", icon = icon("weight"),
             color = if (bmi >= 30) "red" else "yellow")
  })
  output$vbox_mets <- renderValueBox({
    valueBox(paste(input$vat0, "kg"), "Visceral Fat (VAT)",
             icon = icon("apple-alt"), color = "orange")
  })
  output$vbox_homa <- renderValueBox({
    req(sim_result())
    homa <- tail(sim_result()$HOMA_IR, 1)
    valueBox(round(homa, 1), "HOMA-IR", icon = icon("tachometer-alt"),
             color = if (homa > 2.5) "red" else "green")
  })
  output$vbox_hba1c <- renderValueBox({
    req(sim_result())
    hba <- tail(sim_result()$HbA1c, 1)
    valueBox(paste0(round(hba, 1), "%"), "HbA1c", icon = icon("vial"),
             color = if (hba > 6.5) "red" else if (hba > 5.7) "yellow" else "green")
  })

  ## NCEP criteria table
  output$ncep_table <- renderTable({
    req(sim_result())
    df <- tail(sim_result(), 1)
    tribble(
      ~Criterion,              ~Value,           ~Threshold,   ~Status,
      "Waist Circumference",   paste0(round(input$vat0 * 12 + 85), " cm"), "> 102 cm (M)", if(input$vat0 > 1.5) "⚠ ABNORMAL" else "✓ NORMAL",
      "Triglycerides",         paste0(round(df$TRIGLY), " mg/dL"), "≥ 150 mg/dL", if(df$TRIGLY >= 150) "⚠ ABNORMAL" else "✓ NORMAL",
      "HDL-C",                 paste0(round(df$HDLC), " mg/dL"), "< 40 mg/dL (M)", if(df$HDLC < 40) "⚠ ABNORMAL" else "✓ NORMAL",
      "Blood Pressure",        paste0(round(df$SBP), "/", round(df$DBP), " mmHg"), "≥ 130/85 mmHg", if(df$SBP >= 130) "⚠ ABNORMAL" else "✓ NORMAL",
      "Fasting Glucose",       paste0(round(df$FPG), " mg/dL"), "≥ 100 mg/dL", if(df$FPG >= 100) "⚠ ABNORMAL" else "✓ NORMAL"
    )
  }, striped = TRUE, bordered = TRUE)

  ## Radar plot of baseline vs end-of-sim biomarkers
  output$radar_plot <- renderPlotly({
    req(sim_result())
    df_end <- tail(sim_result(), 1)
    cats <- c("HbA1c", "LDL-C", "HDL-C*", "TG", "SBP", "VAT")
    baseline <- c(6.0, 140, 50, 220, 130, input$vat0)
    current  <- c(round(df_end$HbA1c, 1),
                  round(df_end$LDLC),
                  round(df_end$HDLC),
                  round(df_end$TRIGLY),
                  round(df_end$SBP),
                  round(df_end$VAT, 1))
    plot_ly(type = "scatterpolar",
            r = c(baseline, baseline[1]),
            theta = c(cats, cats[1]),
            fill = "toself", name = "Baseline",
            line = list(color = "#E74C3C")) %>%
      add_trace(r = c(current, current[1]),
                theta = c(cats, cats[1]),
                fill = "toself", name = "End of Sim",
                line = list(color = "#2ECC71")) %>%
      layout(polar = list(radialaxis = list(visible = TRUE)),
             showlegend = TRUE,
             title = "Biomarker Radar Chart")
  })

  ## ── Tab 2: PK ────────────────────────────────────────────────────────────
  make_pk_plot <- function(cp_col, drug_name, color = "#3498DB") {
    req(sim_result())
    d <- sim_result()
    if (max(d[[cp_col]]) < 1e-6) {
      return(plotly_empty() %>% layout(title = paste(drug_name, "— Not administered")))
    }
    plot_ly(d, x = ~Week, y = ~.data[[cp_col]], type = "scatter", mode = "lines",
            line = list(color = color, width = 2)) %>%
      layout(title = drug_name, xaxis = list(title = "Week"),
             yaxis = list(title = "Concentration"))
  }
  output$pk_met <- renderPlotly({ make_pk_plot("Cp_MET_out", "Metformin (μg/mL)", "#E67E22") })
  output$pk_glp <- renderPlotly({ make_pk_plot("Cp_GLP_out", "GLP-1 RA (ng/mL)", "#2ECC71") })
  output$pk_sgi <- renderPlotly({ make_pk_plot("Cp_SGI_out", "SGLT2i (ng/mL)",   "#3498DB") })
  output$pk_sta <- renderPlotly({ make_pk_plot("Cp_STA_out", "Statin (μg/mL)",   "#F39C12") })
  output$pk_arb <- renderPlotly({ make_pk_plot("Cp_ARB_out", "ARB (μg/mL)",      "#9B59B6") })

  output$pk_table <- renderDT({
    tribble(
      ~Drug, ~Route, ~F_pct, ~CL_Lh, ~Vc_L, ~t_half_h,
      "Metformin",     "PO", "55%", 28,    380,  9.4,
      "Semaglutide",   "SC", "~95%", 0.056, 8.0, 168,
      "Empagliflozin", "PO", "86%", 12,    73,   12.4,
      "Rosuvastatin",  "PO", "20%", 25,    134,  19,
      "Losartan",      "PO", "33%", 75,    34,   6.2
    ) %>%
      datatable(options = list(dom = "t"), rownames = FALSE) %>%
      formatStyle("Drug", fontWeight = "bold")
  })

  ## ── Tab 3: PD Biomarkers ─────────────────────────────────────────────────
  output$pd_glucose <- renderPlotly({
    req(sim_result())
    d <- sim_result()
    plot_ly(d, x = ~Week) %>%
      add_lines(y = ~GPLAS, name = "Plasma Glucose (mg/dL)", line = list(color = "#E74C3C")) %>%
      add_lines(y = ~GLP1E, name = "Endogenous GLP-1 (pM)", yaxis = "y2",
                line = list(color = "#2ECC71", dash = "dash")) %>%
      layout(title = "Glucose & Endogenous GLP-1",
             yaxis  = list(title = "Glucose (mg/dL)"),
             yaxis2 = list(title = "GLP-1 (pM)", overlaying = "y", side = "right"))
  })

  output$pd_insulin <- renderPlotly({
    req(sim_result())
    d <- sim_result()
    plot_ly(d, x = ~Week) %>%
      add_lines(y = ~IPLAS, name = "Insulin (μU/mL)", line = list(color = "#3498DB")) %>%
      add_lines(y = ~BETA * 100, name = "β-cell Mass (%)", yaxis = "y2",
                line = list(color = "#E67E22", dash = "dash")) %>%
      layout(yaxis = list(title = "Insulin (μU/mL)"),
             yaxis2 = list(title = "β-cell Mass (%)", overlaying = "y", side = "right"))
  })

  output$pd_adipokines <- renderPlotly({
    req(sim_result())
    d <- sim_result()
    plot_ly(d, x = ~Week) %>%
      add_lines(y = ~LEP,    name = "Leptin (ng/mL)",      line = list(color = "#E74C3C")) %>%
      add_lines(y = ~ADIPON, name = "Adiponectin (μg/mL)", line = list(color = "#2ECC71")) %>%
      layout(title = "Adipokines", yaxis = list(title = "Concentration"))
  })

  output$pd_inflam <- renderPlotly({
    req(sim_result())
    d <- sim_result()
    plot_ly(d, x = ~Week) %>%
      add_lines(y = ~TNFA,  name = "TNF-α (pg/mL)", line = list(color = "#E74C3C")) %>%
      add_lines(y = ~IL6C,  name = "IL-6 (pg/mL)",  line = list(color = "#F39C12")) %>%
      add_lines(y = ~CRPC,  name = "hsCRP (mg/L)",  line = list(color = "#9B59B6")) %>%
      layout(title = "Inflammatory Markers", yaxis = list(title = "Concentration"))
  })

  output$pd_bp <- renderPlotly({
    req(sim_result())
    d <- sim_result()
    plot_ly(d, x = ~Week) %>%
      add_lines(y = ~SBP,  name = "SBP (mmHg)",        line = list(color = "#E74C3C")) %>%
      add_lines(y = ~DBP,  name = "DBP (mmHg)",        line = list(color = "#3498DB")) %>%
      add_lines(y = ~ANGII*100, name = "AngII ×100",   line = list(color = "#9B59B6", dash = "dash")) %>%
      layout(title = "Blood Pressure & RAAS", yaxis = list(title = "mmHg / scaled"))
  })

  output$pd_ampk_vat <- renderPlotly({
    req(sim_result())
    d <- sim_result()
    plot_ly(d, x = ~Week) %>%
      add_lines(y = ~AMPKC, name = "AMPK Activity (rel)", line = list(color = "#27AE60")) %>%
      add_lines(y = ~VAT,   name = "VAT (kg)",           yaxis = "y2",
                line = list(color = "#E74C3C", dash = "dash")) %>%
      layout(yaxis  = list(title = "AMPK (relative)"),
             yaxis2 = list(title = "VAT (kg)", overlaying = "y", side = "right"))
  })

  ## ── Tab 4: Clinical Endpoints ─────────────────────────────────────────────
  output$ep_hba1c <- renderPlotly({
    req(sim_result())
    d <- sim_result()
    plot_ly(d, x = ~Week, y = ~HbA1c, type = "scatter", mode = "lines",
            line = list(color = "#E74C3C", width = 2)) %>%
      add_segments(x = 0, xend = max(d$Week), y = 7, yend = 7,
                   line = list(dash = "dot", color = "grey"), name = "Target 7%") %>%
      layout(title = "HbA1c (%)", yaxis = list(title = "HbA1c (%)"))
  })

  output$ep_lipids <- renderPlotly({
    req(sim_result())
    d <- sim_result()
    plot_ly(d, x = ~Week) %>%
      add_lines(y = ~LDLC,   name = "LDL-C",  line = list(color = "#E74C3C")) %>%
      add_lines(y = ~HDLC,   name = "HDL-C",  line = list(color = "#2ECC71")) %>%
      add_lines(y = ~TRIGLY, name = "TG",      line = list(color = "#F39C12")) %>%
      layout(title = "Lipid Panel", yaxis = list(title = "mg/dL"))
  })

  output$ep_bp <- renderPlotly({
    req(sim_result())
    d <- sim_result()
    plot_ly(d, x = ~Week) %>%
      add_lines(y = ~SBP, name = "SBP", line = list(color = "#E74C3C")) %>%
      add_lines(y = ~DBP, name = "DBP", line = list(color = "#3498DB")) %>%
      layout(title = "BP Response", yaxis = list(title = "mmHg"))
  })

  output$ep_metsz <- renderPlotly({
    req(sim_result())
    d <- sim_result()
    plot_ly(d, x = ~Week, y = ~MetS_Zscore, type = "scatter", mode = "lines",
            line = list(color = "#9B59B6", width = 2)) %>%
      layout(title = "MetS Z-score", yaxis = list(title = "Z-score"))
  })

  output$ep_summary_table <- renderDT({
    req(sim_result())
    df <- tail(sim_result(), 1)
    tribble(
      ~Endpoint, ~Value, ~Target, ~Status,
      "HbA1c (%)",          round(df$HbA1c, 1),   "< 7.0",   if(df$HbA1c < 7) "✓" else "✗",
      "FPG (mg/dL)",        round(df$FPG),          "< 100",   if(df$FPG < 100) "✓" else "✗",
      "LDL-C (mg/dL)",      round(df$LDLC),         "< 100",   if(df$LDLC < 100) "✓" else "✗",
      "HDL-C (mg/dL)",      round(df$HDLC),         "> 40",    if(df$HDLC > 40) "✓" else "✗",
      "TG (mg/dL)",         round(df$TRIGLY),        "< 150",   if(df$TRIGLY < 150) "✓" else "✗",
      "SBP (mmHg)",         round(df$SBP),           "< 130",   if(df$SBP < 130) "✓" else "✗",
      "DBP (mmHg)",         round(df$DBP),           "< 85",    if(df$DBP < 85) "✓" else "✗",
      "HOMA-IR",            round(df$HOMA_IR, 1),   "< 2.5",   if(df$HOMA_IR < 2.5) "✓" else "✗",
      "VAT (kg)",           round(df$VAT, 2),        "< 4.0",   if(df$VAT < 4) "✓" else "✗",
      "MetS Z-score",       round(df$MetS_Zscore,1), "< 0",     if(df$MetS_Zscore < 0) "✓" else "✗"
    ) %>%
      datatable(options = list(dom = "t"), rownames = FALSE) %>%
      formatStyle("Status", color = styleEqual(c("✓","✗"), c("green","red")))
  })

  ## ── Tab 5: Scenario Comparison ────────────────────────────────────────────
  scenarios_data <- reactive({
    withProgress(message = "Running all 5 scenarios...", {
      wks <- input$sim_weeks
      bw  <- input$bw
      v0  <- input$vat0
      s0  <- input$sat0

      sc1 <- run_scenario(list(use_met=F,use_glp=F,use_sgi=F,use_sta=F,use_arb=F,
                               met_dose=0,glp_dose=0,sgi_dose=0,sta_dose=0,arb_dose=0),
                          wks, bw, v0, s0) %>% mutate(Scenario = "No Treatment")
      incProgress(0.2)
      sc2 <- run_scenario(list(use_met=T,use_glp=F,use_sgi=F,use_sta=F,use_arb=F,
                               met_dose=1000,glp_dose=0,sgi_dose=0,sta_dose=0,arb_dose=0),
                          wks, bw, v0, s0) %>% mutate(Scenario = "Metformin")
      incProgress(0.2)
      sc3 <- run_scenario(list(use_met=T,use_glp=T,use_sgi=F,use_sta=F,use_arb=F,
                               met_dose=1000,glp_dose=1.0,sgi_dose=0,sta_dose=0,arb_dose=0),
                          wks, bw, v0, s0) %>% mutate(Scenario = "GLP-1 RA + Met")
      incProgress(0.2)
      sc4 <- run_scenario(list(use_met=T,use_glp=F,use_sgi=T,use_sta=F,use_arb=F,
                               met_dose=1000,glp_dose=0,sgi_dose=10,sta_dose=0,arb_dose=0),
                          wks, bw, v0, s0) %>% mutate(Scenario = "SGLT2i + Met")
      incProgress(0.2)
      sc5 <- run_scenario(list(use_met=T,use_glp=T,use_sgi=F,use_sta=T,use_arb=T,
                               met_dose=1000,glp_dose=1.0,sgi_dose=0,sta_dose=10,arb_dose=50),
                          wks, bw, v0, s0) %>% mutate(Scenario = "Quadruple Therapy")
      incProgress(0.2)
      bind_rows(sc1, sc2, sc3, sc4, sc5)
    })
  })

  sc_cols <- c("No Treatment"="red","Metformin"="orange",
               "GLP-1 RA + Met"="green","SGLT2i + Met"="steelblue","Quadruple Therapy"="purple")

  output$sc_hba1c <- renderPlotly({
    d <- scenarios_data()
    plot_ly(d, x = ~Week, y = ~HbA1c, color = ~Scenario, colors = sc_cols,
            type = "scatter", mode = "lines") %>%
      layout(title = "HbA1c Comparison", yaxis = list(title = "HbA1c (%)"))
  })
  output$sc_ldl <- renderPlotly({
    d <- scenarios_data()
    plot_ly(d, x = ~Week, y = ~LDLC, color = ~Scenario, colors = sc_cols,
            type = "scatter", mode = "lines") %>%
      layout(title = "LDL-C Comparison", yaxis = list(title = "mg/dL"))
  })
  output$sc_sbp <- renderPlotly({
    d <- scenarios_data()
    plot_ly(d, x = ~Week, y = ~SBP, color = ~Scenario, colors = sc_cols,
            type = "scatter", mode = "lines") %>%
      layout(title = "SBP Comparison", yaxis = list(title = "mmHg"))
  })
  output$sc_vat <- renderPlotly({
    d <- scenarios_data()
    plot_ly(d, x = ~Week, y = ~VAT, color = ~Scenario, colors = sc_cols,
            type = "scatter", mode = "lines") %>%
      layout(title = "VAT Comparison", yaxis = list(title = "kg"))
  })

  output$sc_summary_table <- renderDT({
    req(scenarios_data())
    d <- scenarios_data()
    wk_max <- max(d$Week)
    d %>%
      filter(abs(Week - wk_max) < 0.5) %>%
      group_by(Scenario) %>%
      summarise(
        `HbA1c (%)` = round(mean(HbA1c), 1),
        `LDL-C (mg/dL)` = round(mean(LDLC)),
        `HDL-C (mg/dL)` = round(mean(HDLC)),
        `TG (mg/dL)` = round(mean(TRIGLY)),
        `SBP (mmHg)` = round(mean(SBP)),
        `VAT (kg)` = round(mean(VAT), 2),
        `HOMA-IR` = round(mean(HOMA_IR), 1),
        `MetS Z-score` = round(mean(MetS_Zscore), 2),
        .groups = "drop"
      ) %>%
      datatable(options = list(dom = "t"), rownames = FALSE)
  })

  ## ── Tab 6: Sensitivity Analysis ───────────────────────────────────────────
  sens_bmi <- reactive({
    bmi_vals <- seq(25, 45, by = 2.5)
    map_dfr(bmi_vals, function(bmi_i) {
      ht <- 1.70
      bw_i <- bmi_i * ht^2
      vat_i <- 2.0 + 0.08 * (bmi_i - 25)
      run_scenario(list(use_met=T,use_glp=F,use_sgi=F,use_sta=F,use_arb=F,
                        met_dose=1000,glp_dose=0,sgi_dose=0,sta_dose=0,arb_dose=0),
                   input$sim_weeks, bw_i, vat_i, 18) %>%
        filter(time == max(time)) %>%
        transmute(BMI = bmi_i, HbA1c, HOMA_IR, LDLC, SBP, VAT)
    })
  })

  output$sa_bmi_hba1c <- renderPlotly({
    d <- sens_bmi()
    plot_ly(d, x = ~BMI, y = ~HbA1c, type = "scatter", mode = "lines+markers",
            line = list(color = "#E74C3C")) %>%
      layout(title = "BMI vs HbA1c (Metformin)",
             xaxis = list(title = "BMI (kg/m²)"),
             yaxis = list(title = "HbA1c (%)"))
  })
  output$sa_bmi_homa <- renderPlotly({
    d <- sens_bmi()
    plot_ly(d, x = ~BMI, y = ~HOMA_IR, type = "scatter", mode = "lines+markers",
            line = list(color = "#F39C12")) %>%
      layout(title = "BMI vs HOMA-IR (Metformin)",
             xaxis = list(title = "BMI (kg/m²)"),
             yaxis = list(title = "HOMA-IR"))
  })

  output$sa_dose_response <- renderPlotly({
    drug <- input$sa_drug
    endpoint <- input$sa_endpoint
    doses <- switch(drug,
      "Metformin"  = seq(250, 2000, by = 250),
      "GLP-1 RA"   = seq(0.25, 2.0, by = 0.25),
      "SGLT2i"     = c(5, 10, 15, 20, 25),
      "Statin"      = c(5, 10, 20, 40),
      "ARB"         = c(25, 50, 100)
    )
    cmt_name <- switch(drug,
      "Metformin"  = "MET_GUT",
      "GLP-1 RA"   = "GLP_SC",
      "SGLT2i"     = "SGI_GUT",
      "Statin"     = "STA_GUT",
      "ARB"        = "ARB_GUT"
    )
    ii_val <- switch(drug,
      "Metformin"  = 12,
      "GLP-1 RA"   = 168,
      "SGLT2i"     = 24,
      "Statin"     = 24,
      "ARB"        = 24
    )
    addl_val <- switch(drug,
      "Metformin"  = 2 * 7 * input$sim_weeks - 1,
      "GLP-1 RA"   = input$sim_weeks - 1,
      "SGLT2i"     = 7 * input$sim_weeks - 1,
      "Statin"     = 7 * input$sim_weeks - 1,
      "ARB"        = 7 * input$sim_weeks - 1
    )
    dr_results <- map_dfr(doses, function(d) {
      ev_d <- ev(ID = 1, amt = d, cmt = cmt_name, ii = ii_val, addl = addl_val)
      ms_mod %>%
        param(BW = input$bw, VAT0 = input$vat0, SAT0 = input$sat0) %>%
        init(VAT = input$vat0, SAT = input$sat0) %>%
        ev(ev_d) %>%
        mrgsim(end = input$sim_weeks * 168, delta = 168) %>%
        as_tibble() %>%
        filter(time == max(time)) %>%
        transmute(Dose = d, Value = .data[[endpoint]])
    })
    plot_ly(dr_results, x = ~Dose, y = ~Value, type = "scatter", mode = "lines+markers",
            line = list(color = "#3498DB")) %>%
      layout(title = paste(drug, "Dose-Response:", endpoint),
             xaxis = list(title = paste("Dose —", drug)),
             yaxis = list(title = endpoint))
  })
}

## ─────────────────────────────────────────────────────────────────────────────
## Run App
## ─────────────────────────────────────────────────────────────────────────────
shinyApp(ui = ui, server = server)
