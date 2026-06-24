##############################################################################
## VTE QSP Interactive Dashboard — Shiny App
## Venous Thromboembolism: Coagulation, PK/PD, Clinical Simulation
##
## Tabs:
##   1. Patient Profile & Risk Assessment
##   2. Drug Pharmacokinetics
##   3. Coagulation PD & Biomarkers
##   4. Thrombus Dynamics
##   5. Treatment Scenario Comparison
##   6. Biomarker Dashboard (INR/Anti-Xa/D-dimer)
##############################################################################

library(shiny)
library(shinydashboard)
library(mrgsolve)
library(dplyr)
library(ggplot2)
library(tidyr)
library(plotly)
library(DT)

## ── Embedded minimal model (for Shiny) ────────────────────────────────────
vte_code <- '
$PARAM
KA_RIV=1.2 F_RIV=0.93 CL_RIV=4.8 V1_RIV=33 Q_RIV=3.2 V2_RIV=20.5
KA_APIX=0.78 F_APIX=0.50 CL_APIX=3.3 V1_APIX=23
KA_DABI=0.35 F_DABI=0.065 CL_DABI=8.5 V1_DABI=80
KA_WARF=0.9 CL_WARF=0.20 V1_WARF=9.5 IC50_WARF=0.65 HILL_W=1.0
KA_ENOX=0.23 CL_ENOX=0.82 V1_ENOX=4.5
KG_FXa=0.08 KF_FXa=0.12 FXa_base=1.0
EMAX_RIV=0.97 EC50_RIV=12.0 HILL_RIV=1.3
EMAX_APIX=0.97 EC50_APIX=5.0 HILL_APX=1.2
EMAX_ENOX=0.85 EC50_ENOX=0.35
KG_FIIa=0.15 KF_FIIa=0.20 FIIa_base=1.0
EMAX_DABI=0.95 EC50_DABI=35 HILL_DABI=1.0
KF_FBR=0.05 KL_FBR=0.03
KG_CLOT=0.04 KD_CLOT=0.008 CLOT_init=80
KP_FORM=0.06 KP_DECAY=0.15 PAI1_eff=0.7 PLASMIN_base=0.5
K_DDIMER=0.10 K_DDIMER_CL=0.025
KSYN_FII=0.012 KDEG_FII=0.012 FII_init=100
KSYN_FVII=0.12 KDEG_FVII=0.12 FVII_init=100
KSYN_FX=0.017 KDEG_FX=0.017 FX_init=100
KIN_VK=0.14 KOUT_VK=0.18 VK0_ox=1.0
eGFR_pat=90 eGFR_ref=90 BWT=70 CLOT_SIZE_INIT=80

$CMT RIV_GUT RIV_CENT RIV_PERIPH APIX_CENT DABI_CENT WARF_CENT ENOX_CENT
     FXa_ACT FIIa_ACT FIBRIN_FORM CLOT_SIZE PLASMIN_ACT DDIMER_CONC
     VK_OX VK_RED FVII_POOL FX_POOL FII_POOL

$INIT RIV_GUT=0 RIV_CENT=0 RIV_PERIPH=0 APIX_CENT=0 DABI_CENT=0
      WARF_CENT=0 ENOX_CENT=0
      FXa_ACT=1 FIIa_ACT=1 FIBRIN_FORM=0 CLOT_SIZE=80 PLASMIN_ACT=0.5
      DDIMER_CONC=2.5 VK_OX=1 VK_RED=0.778
      FVII_POOL=100 FX_POOL=100 FII_POOL=100

$MAIN
double VK0_red = VK0_ox * KIN_VK / KOUT_VK;
double RF = eGFR_pat/eGFR_ref;
double Cp_RIV  = RIV_CENT / V1_RIV;
double Cp_APIX = APIX_CENT / V1_APIX;
double Cp_DABI = DABI_CENT / V1_DABI;
double Cp_WARF = WARF_CENT / V1_WARF;
double Cp_ENOX = ENOX_CENT / V1_ENOX;

double INH_FXa_RIV  = EMAX_RIV * pow(Cp_RIV,HILL_RIV)/(pow(EC50_RIV,HILL_RIV)+pow(Cp_RIV,HILL_RIV));
double INH_FXa_APIX = EMAX_APIX*pow(Cp_APIX,HILL_APX)/(pow(EC50_APIX,HILL_APX)+pow(Cp_APIX,HILL_APX));
double INH_FXa_ENOX = EMAX_ENOX*Cp_ENOX/(EC50_ENOX+Cp_ENOX);
double INH_FXa_TOT  = 1.0-(1.0-INH_FXa_RIV)*(1.0-INH_FXa_APIX)*(1.0-INH_FXa_ENOX);
double INH_FIIa_DABI= EMAX_DABI*pow(Cp_DABI,HILL_DABI)/(pow(EC50_DABI,HILL_DABI)+pow(Cp_DABI,HILL_DABI));
double INH_FIIa_TOT = 1.0-(1.0-INH_FIIa_DABI)*(1.0-INH_FXa_TOT*0.7);

double WARF_INH = pow(Cp_WARF,HILL_W)/(pow(IC50_WARF,HILL_W)+pow(Cp_WARF,HILL_W));
double VK_ratio = VK_RED/VK0_red;
double PT_pct = (FVII_POOL/100.0)*0.5+(FX_POOL/100.0)*0.3+(FII_POOL/100.0)*0.2;
double INR = 12.0/(12.0*(PT_pct+0.001));
double aPTT = 32.0*(1.0+2.5*INH_FIIa_DABI);
double ANTI_XA = Cp_ENOX;

$ODE
dxdt_RIV_GUT   = -KA_RIV*RIV_GUT;
dxdt_RIV_CENT  =  KA_RIV*RIV_GUT-(CL_RIV/V1_RIV)*RIV_CENT-(Q_RIV/V1_RIV)*RIV_CENT+(Q_RIV/V2_RIV)*RIV_PERIPH;
dxdt_RIV_PERIPH=  (Q_RIV/V1_RIV)*RIV_CENT-(Q_RIV/V2_RIV)*RIV_PERIPH;
dxdt_APIX_CENT = -KA_APIX*APIX_CENT-(CL_APIX/V1_APIX)*APIX_CENT;
dxdt_DABI_CENT = -KA_DABI*DABI_CENT-(CL_DABI/V1_DABI)/RF*DABI_CENT;
dxdt_WARF_CENT = -KA_WARF*WARF_CENT-(CL_WARF/V1_WARF)*WARF_CENT;
dxdt_ENOX_CENT = -KA_ENOX*ENOX_CENT-(CL_ENOX/V1_ENOX)/RF*ENOX_CENT;
double FXa_gen = KG_FXa*FXa_base*(FX_POOL/100.0);
dxdt_FXa_ACT   = FXa_gen*(1.0-INH_FXa_TOT)-KF_FXa*FXa_ACT;
double FIIa_gen = KG_FIIa*FXa_ACT*(FII_POOL/100.0);
dxdt_FIIa_ACT  = FIIa_gen*(1.0-INH_FIIa_TOT)-KF_FIIa*FIIa_ACT;
dxdt_FIBRIN_FORM= KF_FBR*FIIa_ACT*FIIa_base-KL_FBR*PLASMIN_ACT*FIBRIN_FORM;
dxdt_CLOT_SIZE  = KG_CLOT*FIBRIN_FORM-KD_CLOT*PLASMIN_ACT*CLOT_SIZE;
dxdt_PLASMIN_ACT= KP_FORM*(1.0-PAI1_eff*0.3)-KP_DECAY*PLASMIN_ACT;
dxdt_DDIMER_CONC= K_DDIMER*PLASMIN_ACT*FIBRIN_FORM-K_DDIMER_CL*DDIMER_CONC;
double VK0r = VK0_ox*KIN_VK/KOUT_VK;
dxdt_VK_OX  = KOUT_VK*VK_RED-KIN_VK*VK_OX;
dxdt_VK_RED = KIN_VK*VK_OX*(1.0-WARF_INH)-KOUT_VK*VK_RED;
double STIM_VK = VK_RED/VK0r;
dxdt_FVII_POOL = KSYN_FVII*STIM_VK*FVII_init-KDEG_FVII*FVII_POOL;
dxdt_FX_POOL   = KSYN_FX*STIM_VK*FX_init-KDEG_FX*FX_POOL;
dxdt_FII_POOL  = KSYN_FII*STIM_VK*FII_init-KDEG_FII*FII_POOL;

$TABLE
capture Cp_RIV=Cp_RIV; capture Cp_APIX=Cp_APIX; capture Cp_DABI=Cp_DABI;
capture Cp_WARF=Cp_WARF; capture ANTI_XA=ANTI_XA;
capture INR=INR; capture aPTT_out=aPTT;
capture INH_FXa=INH_FXa_TOT*100; capture INH_FIIa=INH_FIIa_TOT*100;
capture TG_ETP=FIIa_ACT*100;
capture DDIMER=DDIMER_CONC; capture CLOT_PCT=CLOT_SIZE/80*100;
capture FVII_pct=FVII_POOL; capture FX_pct=FX_POOL; capture FII_pct=FII_POOL;
$CAPTURE Cp_RIV Cp_APIX Cp_DABI Cp_WARF ANTI_XA INR aPTT_out
         INH_FXa INH_FIIa TG_ETP DDIMER CLOT_PCT
         FVII_pct FX_pct FII_pct
'

mod_shiny <- mcode("vte_shiny", vte_code, quiet = TRUE)

## ── UI ─────────────────────────────────────────────────────────────────────

ui <- dashboardPage(
  skin = "blue",
  dashboardHeader(
    title = "VTE QSP Dashboard",
    titleWidth = 280
  ),

  dashboardSidebar(
    width = 280,
    sidebarMenu(
      menuItem("① Patient Profile",        tabName = "tab_patient",  icon = icon("user-md")),
      menuItem("② Drug Pharmacokinetics",  tabName = "tab_pk",       icon = icon("pills")),
      menuItem("③ Coagulation PD",         tabName = "tab_pd",       icon = icon("tint")),
      menuItem("④ Thrombus Dynamics",      tabName = "tab_clot",     icon = icon("circle")),
      menuItem("⑤ Scenario Comparison",    tabName = "tab_compare",  icon = icon("chart-bar")),
      menuItem("⑥ Biomarker Dashboard",    tabName = "tab_bm",       icon = icon("microscope"))
    ),
    hr(),
    h5("  Global Settings", style = "color:#aaa; padding-left:15px;"),
    sliderInput("sim_days", "Simulation Duration (days)",
                min = 7, max = 365, value = 90, step = 7),
    sliderInput("bwt", "Body Weight (kg)", min = 40, max = 150, value = 70),
    sliderInput("egfr", "eGFR (mL/min/1.73m²)", min = 15, max = 130, value = 90),
    selectInput("drug_choice", "Primary Drug",
                choices = c("Rivaroxaban DVT" = "riv_dvt",
                             "Rivaroxaban PE"  = "riv_pe",
                             "Apixaban PE"     = "apix_pe",
                             "Dabigatran DVT"  = "dabi_dvt",
                             "Warfarin+Bridge" = "warf",
                             "Enoxaparin Proph"= "enox_proph"),
                selected = "riv_dvt")
  ),

  dashboardBody(
    tags$head(tags$style(HTML("
      .content-wrapper { background-color: #f4f6f9; }
      .box-header h3 { font-size: 14px; }
      .info-box-content .info-box-number { font-size: 22px; }
    "))),

    tabItems(

      ## ── TAB 1: Patient Profile & Risk ───────────────────────────────────
      tabItem("tab_patient",
        fluidRow(
          box(title = "Patient Demographics", status = "primary",
              solidHeader = TRUE, width = 4,
              sliderInput("age_p", "Age (years)", 18, 90, 60),
              selectInput("sex_p", "Sex", c("Male", "Female")),
              sliderInput("bwt_p", "Weight (kg)", 40, 150, 70),
              sliderInput("ht_p", "Height (cm)", 140, 200, 170),
              numericInput("scr_p", "Serum Creatinine (mg/dL)", 1.0, min = 0.3, max = 10)
          ),
          box(title = "VTE Risk Factors (Virchow's Triad)", status = "warning",
              solidHeader = TRUE, width = 4,
              checkboxGroupInput("risk_stasis",
                "Venous Stasis:",
                choices = c("Immobility >3 days" = "imm",
                             "Major surgery"      = "surg",
                             "Long-haul flight"   = "flight",
                             "Hospitalization"    = "hosp"),
                selected = "surg"),
              checkboxGroupInput("risk_endo",
                "Endothelial Injury:",
                choices = c("Prior DVT/PE" = "prior",
                             "Trauma"      = "trauma",
                             "Central line" = "cvc"),
                selected = NULL),
              checkboxGroupInput("risk_hyper",
                "Hypercoagulability:",
                choices = c("Active malignancy" = "malig",
                             "Factor V Leiden"  = "fvl",
                             "Oral contraceptives" = "ocp",
                             "Antiphospholipid Ab" = "aps",
                             "Protein C/S defic." = "pcs",
                             "Pregnancy"           = "preg"),
                selected = NULL)
          ),
          box(title = "VTE Risk Score & Clinical Decision", status = "danger",
              solidHeader = TRUE, width = 4,
              h4("Wells Score (DVT)"),
              verbatimTextOutput("wells_dvt"),
              hr(),
              h4("Simplified Wells (PE)"),
              verbatimTextOutput("wells_pe"),
              hr(),
              h4("Anticoagulation Recommendation"),
              verbatimTextOutput("anticoag_rec"),
              hr(),
              h4("Renal Function Summary"),
              verbatimTextOutput("renal_summary")
          )
        ),
        fluidRow(
          box(title = "Diagnosis: VTE Type", status = "info",
              solidHeader = TRUE, width = 6,
              selectInput("vte_type", "VTE Presentation",
                          choices = c("Proximal DVT (Femoral/Popliteal)" = "prox_dvt",
                                       "Distal DVT (Calf)"               = "dist_dvt",
                                       "Low-Risk PE (PESI I/II)"          = "low_pe",
                                       "Submassive PE (RV Dysfunction)"   = "submass_pe",
                                       "Massive PE (Hemodynamic Unstable)" = "massive_pe"),
                          selected = "prox_dvt"),
              selectInput("vte_provoked", "Provoked / Unprovoked",
                          choices = c("Provoked (transient risk factor)" = "provoked",
                                       "Unprovoked"                      = "unprovoked",
                                       "Cancer-associated"               = "cancer")),
              h4("Duration of Anticoagulation Recommendation:"),
              verbatimTextOutput("duration_rec")
          ),
          box(title = "Contraindications & Precautions", status = "danger",
              solidHeader = TRUE, width = 6,
              checkboxGroupInput("contraind",
                "Contraindications:",
                choices = c("Active bleeding / High ICH risk" = "bleed",
                             "Severe renal impairment (eGFR<15)" = "renal_sev",
                             "Mechanical heart valve"           = "valve",
                             "Pregnancy (warfarin contraindicated)" = "preg_warf",
                             "Antiphospholipid Syndrome (use warfarin)" = "aps_warf")),
              h4("Recommended Anticoagulant:"),
              verbatimTextOutput("drug_rec")
          )
        )
      ),

      ## ── TAB 2: Drug Pharmacokinetics ────────────────────────────────────
      tabItem("tab_pk",
        fluidRow(
          box(title = "Dosing Parameters", status = "primary",
              solidHeader = TRUE, width = 3,
              h5("Rivaroxaban"),
              numericInput("riv_dose", "Dose (mg)", 15, min = 2.5, max = 30, step = 2.5),
              selectInput("riv_reg", "Regimen",
                          c("QD"="qd","BID"="bid","DVT Phase (15BID→20QD)"="dvt_phase")),
              hr(),
              h5("Apixaban"),
              numericInput("apix_dose", "Dose (mg)", 5, min = 2.5, max = 10, step = 2.5),
              selectInput("apix_reg", "Regimen",
                          c("BID"="bid","QD"="qd","PE Phase (10BID→5BID)"="pe_phase")),
              hr(),
              h5("Enoxaparin"),
              numericInput("enox_dose_mg", "Dose (mg)", 40, min = 20, max = 160, step = 10),
              selectInput("enox_reg", "Regimen", c("QD SC"="qd","BID SC"="bid"))
          ),
          box(title = "Rivaroxaban PK Profile", status = "primary",
              solidHeader = TRUE, width = 9,
              plotlyOutput("pk_riv_plot", height = "300px"),
              fluidRow(
                infoBoxOutput("ib_riv_cmax", width = 4),
                infoBoxOutput("ib_riv_trough", width = 4),
                infoBoxOutput("ib_riv_auc", width = 4)
              )
          )
        ),
        fluidRow(
          box(title = "Apixaban PK Profile", status = "success",
              solidHeader = TRUE, width = 6,
              plotlyOutput("pk_apix_plot", height = "250px")),
          box(title = "Enoxaparin Anti-Xa Profile", status = "warning",
              solidHeader = TRUE, width = 6,
              plotlyOutput("pk_enox_plot", height = "250px"),
              p("Therapeutic anti-Xa target (treatment): 0.5–1.0 IU/mL",
                style = "font-size:11px; color:gray;"))
        ),
        fluidRow(
          box(title = "Dabigatran PK — Normal vs. Renal Impairment",
              status = "danger", solidHeader = TRUE, width = 6,
              numericInput("dabi_dose", "Dabigatran Dose (mg)", 150, min = 75, max = 150),
              plotlyOutput("pk_dabi_renalcomp", height = "250px")),
          box(title = "Warfarin PK + VK Cycle",
              status = "info", solidHeader = TRUE, width = 6,
              numericInput("warf_dose", "Warfarin Dose (mg/day)", 5, min = 1, max = 15),
              plotlyOutput("pk_warf_plot", height = "250px"))
        )
      ),

      ## ── TAB 3: Coagulation PD & Biomarkers ──────────────────────────────
      tabItem("tab_pd",
        fluidRow(
          box(title = "FXa & Thrombin Inhibition (% Inhibition)", status = "danger",
              solidHeader = TRUE, width = 6,
              plotlyOutput("pd_inh_plot", height = "300px"),
              p("Peak inhibition reflects coagulation cascade suppression",
                style = "font-size:11px;")),
          box(title = "Thrombin Generation (ETP Proxy)", status = "warning",
              solidHeader = TRUE, width = 6,
              plotlyOutput("pd_tg_plot", height = "300px"))
        ),
        fluidRow(
          box(title = "INR & Factor Levels (Warfarin)", status = "primary",
              solidHeader = TRUE, width = 6,
              plotlyOutput("pd_inr_plot", height = "300px"),
              p("Dotted lines: INR therapeutic range 2.0–3.0",
                style = "font-size:11px;")),
          box(title = "DOAC Concentration–Effect Relationship",
              status = "success", solidHeader = TRUE, width = 6,
              selectInput("doac_pd_drug", "Drug", c("Rivaroxaban","Apixaban","Dabigatran")),
              plotlyOutput("pd_emax_plot", height = "300px"))
        ),
        fluidRow(
          box(title = "PK-PD Parameters Summary",
              status = "info", solidHeader = TRUE, width = 12,
              DTOutput("pd_params_table"))
        )
      ),

      ## ── TAB 4: Thrombus Dynamics ─────────────────────────────────────────
      tabItem("tab_clot",
        fluidRow(
          box(title = "Thrombus Dynamics Controls", status = "primary",
              solidHeader = TRUE, width = 3,
              sliderInput("clot_init_size", "Initial Clot Size (%)", 10, 100, 80),
              sliderInput("fibrinolysis_str", "Endogenous Fibrinolysis", 0, 1, 0.5,
                          step = 0.1,
                          post = " (0=low, 1=high)"),
              sliderInput("pai1_level", "PAI-1 Level (fibrinolysis inhib.)",
                          0, 1, 0.7, step = 0.1),
              checkboxInput("add_thrombolytic", "Add Thrombolytic (Alteplase)", FALSE),
              conditionalPanel("input.add_thrombolytic",
                numericInput("alteplase_dose", "Alteplase Dose (mg)", 100, min = 0.6, max = 100),
                numericInput("alteplase_time", "Admin Time (hours after dx)", 0)
              )
          ),
          box(title = "Thrombus Resolution Over Time", status = "danger",
              solidHeader = TRUE, width = 9,
              plotlyOutput("clot_res_plot", height = "350px"),
              fluidRow(
                infoBoxOutput("ib_d7_clot", width = 4),
                infoBoxOutput("ib_d30_clot", width = 4),
                infoBoxOutput("ib_d90_clot", width = 4)
              )
          )
        ),
        fluidRow(
          box(title = "Fibrin Formation & Plasmin Dynamics",
              status = "warning", solidHeader = TRUE, width = 6,
              plotlyOutput("fibrin_plasmin_plot", height = "280px")),
          box(title = "D-dimer Trajectory (VTE Treatment Biomarker)",
              status = "success", solidHeader = TRUE, width = 6,
              plotlyOutput("ddimer_traj_plot", height = "280px"),
              p("D-dimer <500 ng/mL: consider VTE excluded (with low pre-test prob.)",
                style = "font-size:11px; color:gray;"))
        )
      ),

      ## ── TAB 5: Scenario Comparison ──────────────────────────────────────
      tabItem("tab_compare",
        fluidRow(
          box(title = "Select Scenarios to Compare", status = "primary",
              solidHeader = TRUE, width = 3,
              checkboxGroupInput("scenarios",
                "Treatment Scenarios:",
                choices = c("Rivaroxaban DVT (15→20 mg)"  = "riv_dvt",
                             "Apixaban PE (10→5 mg BID)"  = "apix_pe",
                             "Warfarin + LMWH Bridge"     = "warf",
                             "Dabigatran 150 BID"         = "dabi",
                             "Enoxaparin 1 mg/kg BID"     = "enox",
                             "Extended: Riv 10 QD"        = "riv_ext"),
                selected = c("riv_dvt","apix_pe","warf")),
              hr(),
              selectInput("compare_endpoint",
                          "Comparison Endpoint:",
                          c("Clot Resolution (%)"  = "CLOT_PCT",
                            "D-dimer (ng/mL)"      = "DDIMER",
                            "FXa Inhibition (%)"   = "INH_FXa",
                            "FIIa Inhibition (%)"  = "INH_FIIa",
                            "INR"                  = "INR",
                            "Thrombin Generation (ETP)" = "TG_ETP"))
          ),
          box(title = "Comparative Efficacy Plot", status = "primary",
              solidHeader = TRUE, width = 9,
              plotlyOutput("compare_plot", height = "400px"))
        ),
        fluidRow(
          box(title = "Efficacy & Safety Summary Table",
              status = "success", solidHeader = TRUE, width = 12,
              DTOutput("summary_table"))
        )
      ),

      ## ── TAB 6: Biomarker Dashboard ───────────────────────────────────────
      tabItem("tab_bm",
        fluidRow(
          infoBoxOutput("ib_inr_ss",   width = 3),
          infoBoxOutput("ib_antixa_pk",width = 3),
          infoBoxOutput("ib_ddimer_d7",width = 3),
          infoBoxOutput("ib_clot_d30", width = 3)
        ),
        fluidRow(
          box(title = "Lab Monitoring Dashboard",
              status = "primary", solidHeader = TRUE, width = 6,
              plotlyOutput("bm_labs_plot", height = "350px"),
              p("Therapeutic ranges: INR 2-3 (warfarin); Anti-Xa 0.5-1.0 IU/mL (LMWH)",
                style = "font-size:11px; color:gray;")),
          box(title = "Coagulation Cascade Activity Over Treatment",
              status = "danger", solidHeader = TRUE, width = 6,
              plotlyOutput("bm_coag_plot", height = "350px"))
        ),
        fluidRow(
          box(title = "D-dimer Normalization & Risk Stratification",
              status = "warning", solidHeader = TRUE, width = 6,
              plotlyOutput("bm_ddimer_plot", height = "280px"),
              p("Elevated D-dimer 3 months post-treatment → ↑ recurrence risk",
                style = "font-size:11px; color:gray;")),
          box(title = "Biomarker Interpretation Guide",
              status = "info", solidHeader = TRUE, width = 6,
              DTOutput("bm_interp_table"))
        )
      )
    )
  )
)

## ── SERVER ─────────────────────────────────────────────────────────────────

server <- function(input, output, session) {

  ## Reactive: run primary simulation
  sim_primary <- reactive({
    dur_h <- input$sim_days * 24
    bwt   <- input$bwt
    egfr  <- input$egfr

    drug <- input$drug_choice
    dose_ev <- switch(drug,
      riv_dvt  = c(
        ev(amt = 15 * 0.93, cmt = "RIV_GUT", ii = 12, addl = 41, time = 0),
        ev(amt = 20 * 0.93, cmt = "RIV_GUT", ii = 24,
           addl = max(0, input$sim_days - 22), time = 21 * 24)
      ),
      riv_pe   = ev(amt = 20 * 0.93, cmt = "RIV_GUT", ii = 24,
                    addl = input$sim_days - 1),
      apix_pe  = c(
        ev(amt = 10 * 0.50, cmt = "APIX_CENT", ii = 12, addl = 13, time = 0),
        ev(amt = 5  * 0.50, cmt = "APIX_CENT", ii = 12,
           addl = max(0, (input$sim_days - 7) * 2 - 1), time = 7 * 24)
      ),
      dabi_dvt = ev(amt = 150 * 0.065, cmt = "DABI_CENT", ii = 12,
                    addl = input$sim_days * 2 - 1),
      warf     = c(
        ev(amt = 5 * 1.0, cmt = "WARF_CENT", ii = 24, addl = input$sim_days - 1),
        ev(amt = 70 * 0.92, cmt = "ENOX_CENT", ii = 12, addl = 19, time = 0)
      ),
      enox_proph = ev(amt = 40 * 0.92, cmt = "ENOX_CENT", ii = 24,
                      addl = min(13, input$sim_days - 1))
    )

    mod_shiny %>%
      param(BWT = bwt, eGFR_pat = egfr,
            CLOT_init = 80) %>%
      mrgsim(events = dose_ev, end = dur_h, delta = 2) %>%
      as.data.frame() %>%
      mutate(time_d = time / 24)
  })

  ## ── Tab 2: PK Plots ──────────────────────────────────────────────────────
  output$pk_riv_plot <- renderPlotly({
    sim <- sim_primary()
    p <- sim %>%
      filter(time_d <= min(10, input$sim_days)) %>%
      ggplot(aes(x = time_d, y = Cp_RIV)) +
      geom_line(color = "#e53935", linewidth = 1) +
      geom_hline(yintercept = c(20, 200), linetype = "dashed",
                 color = c("blue", "red"), alpha = 0.7) +
      labs(x = "Time (days)", y = "Rivaroxaban Cp (ng/mL)",
           title = "Rivaroxaban Plasma Concentration") +
      theme_bw(base_size = 11)
    ggplotly(p)
  })

  output$pk_apix_plot <- renderPlotly({
    dur_h <- min(10, input$sim_days) * 24
    dose_apix <- c(
      ev(amt = 10 * 0.50, cmt = "APIX_CENT", ii = 12, addl = 13, time = 0),
      ev(amt = 5  * 0.50, cmt = "APIX_CENT", ii = 12, addl = 79, time = 7 * 24)
    )
    sim_a <- mod_shiny %>%
      param(BWT = input$bwt, eGFR_pat = input$egfr) %>%
      mrgsim(events = dose_apix, end = 240, delta = 0.5) %>%
      as.data.frame() %>%
      mutate(time_d = time / 24)
    p <- sim_a %>%
      ggplot(aes(x = time_d, y = Cp_APIX)) +
      geom_line(color = "#1e88e5", linewidth = 1) +
      labs(x = "Time (days)", y = "Apixaban Cp (ng/mL)",
           title = "Apixaban: 10 BID (d1-7) → 5 BID") +
      theme_bw(base_size = 11)
    ggplotly(p)
  })

  output$pk_enox_plot <- renderPlotly({
    dose_en <- ev(amt = input$enox_dose_mg * 0.92, cmt = "ENOX_CENT",
                  ii = ifelse(input$enox_reg == "bid", 12, 24),
                  addl = 13, time = 0)
    sim_en <- mod_shiny %>%
      param(BWT = input$bwt, eGFR_pat = input$egfr) %>%
      mrgsim(events = dose_en, end = 72, delta = 0.25) %>%
      as.data.frame() %>%
      mutate(time_d = time / 24)
    p <- sim_en %>%
      ggplot(aes(x = time_d, y = ANTI_XA)) +
      geom_line(color = "#43a047", linewidth = 1) +
      geom_hline(yintercept = c(0.5, 1.0), linetype = "dashed",
                 color = "red", alpha = 0.7) +
      annotate("rect", xmin = -Inf, xmax = Inf, ymin = 0.5, ymax = 1.0,
               fill = "green", alpha = 0.1) +
      labs(x = "Time (days)", y = "Anti-Xa (IU/mL)",
           title = paste("Enoxaparin", input$enox_dose_mg, "mg", input$enox_reg)) +
      theme_bw(base_size = 11)
    ggplotly(p)
  })

  output$pk_dabi_renalcomp <- renderPlotly({
    dose_d <- ev(amt = input$dabi_dose * 0.065, cmt = "DABI_CENT",
                 ii = 12, addl = 9 * 2, time = 0)
    sim_norm <- mod_shiny %>% param(eGFR_pat = 90) %>%
      mrgsim(events = dose_d, end = 120, delta = 0.5) %>%
      as.data.frame() %>% mutate(gfr = "Normal (90)", time_d = time/24)
    sim_ckd3 <- mod_shiny %>% param(eGFR_pat = 30) %>%
      mrgsim(events = dose_d, end = 120, delta = 0.5) %>%
      as.data.frame() %>% mutate(gfr = "CKD3 (30)", time_d = time/24)
    df_r <- bind_rows(sim_norm, sim_ckd3)
    p <- df_r %>%
      ggplot(aes(x = time_d, y = Cp_DABI, color = gfr)) +
      geom_line(linewidth = 1) +
      scale_color_manual(values = c("Normal (90)" = "steelblue",
                                     "CKD3 (30)"   = "tomato")) +
      labs(x = "Time (days)", y = "Dabigatran Cp (ng/mL)", color = "eGFR",
           title = "Dabigatran: Renal Sensitivity") +
      theme_bw(base_size = 11)
    ggplotly(p)
  })

  output$pk_warf_plot <- renderPlotly({
    dose_w <- ev(amt = input$warf_dose, cmt = "WARF_CENT",
                 ii = 24, addl = input$sim_days - 1)
    sim_w <- mod_shiny %>%
      param(BWT = input$bwt) %>%
      mrgsim(events = dose_w, end = input$sim_days * 24, delta = 2) %>%
      as.data.frame() %>%
      mutate(time_d = time / 24)
    p <- sim_w %>%
      ggplot(aes(x = time_d, y = Cp_WARF)) +
      geom_line(color = "#8e24aa", linewidth = 1) +
      labs(x = "Time (days)", y = "Warfarin Cp (mg/L)",
           title = paste("Warfarin", input$warf_dose, "mg QD")) +
      theme_bw(base_size = 11)
    ggplotly(p)
  })

  ## Info boxes
  output$ib_riv_cmax <- renderInfoBox({
    sim <- sim_primary()
    infoBox("Rivaroxaban Cmax",
            paste(round(max(sim$Cp_RIV[sim$time_d <= 1]), 1), "ng/mL"),
            icon = icon("arrow-up"), color = "red")
  })
  output$ib_riv_trough <- renderInfoBox({
    sim <- sim_primary()
    infoBox("Rivaroxaban Trough",
            paste(round(min(sim$Cp_RIV[sim$time_d >= 2 & sim$time_d <= 3]), 1), "ng/mL"),
            icon = icon("arrow-down"), color = "blue")
  })
  output$ib_riv_auc <- renderInfoBox({
    sim <- sim_primary()
    auc_24h <- sum(diff(sim$time[sim$time_d <= 1]) *
                   (head(sim$Cp_RIV[sim$time_d <= 1], -1) +
                    tail(sim$Cp_RIV[sim$time_d <= 1], -1)) / 2)
    infoBox("Rivaroxaban AUC₀₋₂₄h",
            paste(round(auc_24h, 0), "ng·h/mL"),
            icon = icon("chart-area"), color = "purple")
  })

  ## ── Tab 3: PD Plots ──────────────────────────────────────────────────────
  output$pd_inh_plot <- renderPlotly({
    sim <- sim_primary()
    p <- sim %>%
      select(time_d, INH_FXa, INH_FIIa) %>%
      pivot_longer(c(INH_FXa, INH_FIIa)) %>%
      ggplot(aes(x = time_d, y = value, color = name)) +
      geom_line(linewidth = 1) +
      scale_color_manual(values = c(INH_FXa = "#e53935", INH_FIIa = "#8e24aa"),
                         labels = c("FXa Inhibition (%)", "FIIa Inhibition (%)")) +
      labs(x = "Time (days)", y = "Inhibition (%)", color = "Target",
           title = "Coagulation Factor Inhibition Over Time") +
      theme_bw(base_size = 11)
    ggplotly(p)
  })

  output$pd_tg_plot <- renderPlotly({
    sim <- sim_primary()
    p <- sim %>%
      ggplot(aes(x = time_d, y = TG_ETP)) +
      geom_line(color = "#ff6f00", linewidth = 1) +
      geom_hline(yintercept = 100, linetype = "dashed", color = "gray50") +
      labs(x = "Time (days)", y = "ETP Proxy (% of Baseline)",
           title = "Thrombin Generation (ETP proxy)") +
      theme_bw(base_size = 11)
    ggplotly(p)
  })

  output$pd_inr_plot <- renderPlotly({
    dose_w <- ev(amt = 5, cmt = "WARF_CENT", ii = 24, addl = input$sim_days - 1)
    sim_w <- mod_shiny %>%
      mrgsim(events = dose_w, end = input$sim_days * 24, delta = 2) %>%
      as.data.frame() %>%
      mutate(time_d = time / 24) %>%
      select(time_d, INR, FVII_pct, FX_pct, FII_pct) %>%
      pivot_longer(-time_d)
    p <- sim_w %>%
      ggplot(aes(x = time_d, y = value, color = name)) +
      geom_line(linewidth = 1) +
      geom_hline(yintercept = c(2, 3), linetype = "dashed",
                 color = c("red", "red"), alpha = 0.6) +
      scale_color_manual(values = c(INR = "#8e24aa", FVII_pct = "#e53935",
                                     FX_pct = "#1e88e5", FII_pct = "#43a047"),
                         labels = c("FII (%)", "FVII (%)", "FX (%)", "INR")) +
      labs(x = "Time (days)", y = "Value (% or INR)",
           title = "Warfarin: Factor Depletion & INR", color = "Measure") +
      theme_bw(base_size = 11)
    ggplotly(p)
  })

  output$pd_emax_plot <- renderPlotly({
    conc_range <- switch(input$doac_pd_drug,
      Rivaroxaban = seq(0, 300, by = 1),
      Apixaban    = seq(0, 150, by = 0.5),
      Dabigatran  = seq(0, 500, by = 2)
    )
    params <- list(
      Rivaroxaban = list(EMAX = 0.97, EC50 = 12,  HILL = 1.3),
      Apixaban    = list(EMAX = 0.97, EC50 = 5,   HILL = 1.2),
      Dabigatran  = list(EMAX = 0.95, EC50 = 35,  HILL = 1.0)
    )
    p_sel <- params[[input$doac_pd_drug]]
    effect <- p_sel$EMAX * conc_range^p_sel$HILL /
              (p_sel$EC50^p_sel$HILL + conc_range^p_sel$HILL)
    df_e <- data.frame(Cp = conc_range, Effect = effect * 100)
    p <- df_e %>%
      ggplot(aes(x = Cp, y = Effect)) +
      geom_line(color = "darkblue", linewidth = 1) +
      geom_vline(xintercept = p_sel$EC50, linetype = "dotted", color = "red") +
      annotate("text", x = p_sel$EC50 * 1.1, y = 10,
               label = paste0("EC₅₀ = ", p_sel$EC50, " ng/mL"),
               color = "red", size = 3.5) +
      labs(x = paste(input$doac_pd_drug, "Concentration (ng/mL)"),
           y = "% Inhibition",
           title = paste(input$doac_pd_drug, "Emax Concentration-Effect Curve")) +
      theme_bw(base_size = 11)
    ggplotly(p)
  })

  output$pd_params_table <- renderDT({
    tibble(
      Drug       = c("Rivaroxaban","Apixaban","Edoxaban","Dabigatran","Enoxaparin","Warfarin"),
      Target     = c("FXa","FXa","FXa","FIIa","FXa>FIIa via AT","VitK-epoxide reductase (VKORC1)"),
      `Emax (%)`  = c(97, 97, 97, 95, 85, 97),
      `EC50 (ng/mL or IU/mL)` = c("12 ng/mL","5 ng/mL","11 ng/mL","35 ng/mL","0.35 IU/mL","0.65 mg/L"),
      `t½ (h)`   = c("5-9","12","10-14","12-17","4-7","36-42"),
      `F (%)`    = c("80-100","50","62","3-7","92 (SC)","~100"),
      `Renal (%)` = c("33","27","50","80","65","<5")
    ) %>%
      datatable(rownames = FALSE, options = list(dom = "t", pageLength = 10))
  })

  ## ── Tab 4: Thrombus ───────────────────────────────────────────────────────
  sim_clot <- reactive({
    sim <- sim_primary()
    clot_init_use <- input$clot_init_size
    sim %>%
      mutate(
        CLOT_adj = pmax(0, CLOT_PCT - (100 - clot_init_use)),
        FIBRIN_disp = FIBRIN_FORM * (1 + input$fibrinolysis_str * 0.5),
        PLASMIN_disp = PLASMIN_ACT * (1 + input$fibrinolysis_str * 0.8)
      )
  })

  output$clot_res_plot <- renderPlotly({
    sim <- sim_clot()
    p <- sim %>%
      ggplot(aes(x = time_d, y = CLOT_PCT)) +
      geom_line(color = "#e53935", linewidth = 1.5) +
      geom_hline(yintercept = 50, linetype = "dashed", color = "gray50") +
      annotate("text", x = input$sim_days * 0.5, y = 52,
               label = "50% Clot Reduction", size = 3.5) +
      labs(x = "Time (days)", y = "Residual Thrombus (%)",
           title = "Thrombus Resolution During Anticoagulation") +
      ylim(0, 100) + theme_bw(base_size = 11)
    ggplotly(p)
  })

  output$ib_d7_clot  <- renderInfoBox({
    sim <- sim_clot()
    v <- sim$CLOT_PCT[which.min(abs(sim$time_d - 7))]
    infoBox("Day 7 Residual Clot", paste0(round(v, 1), "%"),
            icon = icon("circle"), color = "red")
  })
  output$ib_d30_clot <- renderInfoBox({
    sim <- sim_clot()
    v <- sim$CLOT_PCT[which.min(abs(sim$time_d - 30))]
    infoBox("Day 30 Residual Clot", paste0(round(v, 1), "%"),
            icon = icon("circle-half-stroke"), color = "yellow")
  })
  output$ib_d90_clot <- renderInfoBox({
    sim <- sim_clot()
    v <- sim$CLOT_PCT[which.min(abs(sim$time_d - min(90, input$sim_days)))]
    infoBox("Day 90 Residual Clot", paste0(round(v, 1), "%"),
            icon = icon("check-circle"), color = "green")
  })

  output$fibrin_plasmin_plot <- renderPlotly({
    sim <- sim_clot()
    p <- sim %>%
      select(time_d, FIBRIN_FORM, PLASMIN_ACT) %>%
      pivot_longer(-time_d) %>%
      ggplot(aes(x = time_d, y = value, color = name)) +
      geom_line(linewidth = 1) +
      scale_color_manual(values = c(FIBRIN_FORM = "#e53935", PLASMIN_ACT = "#1e88e5"),
                         labels = c("Fibrin Formation", "Plasmin Activity")) +
      labs(x = "Time (days)", y = "Relative Activity",
           title = "Fibrin vs. Plasmin Balance", color = "Pathway") +
      theme_bw(base_size = 11)
    ggplotly(p)
  })

  output$ddimer_traj_plot <- renderPlotly({
    sim <- sim_primary()
    p <- sim %>%
      ggplot(aes(x = time_d, y = DDIMER)) +
      geom_line(color = "#43a047", linewidth = 1) +
      geom_hline(yintercept = 0.5, linetype = "dashed", color = "darkgreen") +
      annotate("text", x = input$sim_days * 0.5, y = 0.55,
               label = "Normal limit (500 ng/mL)", color = "darkgreen", size = 3.5) +
      labs(x = "Time (days)", y = "D-dimer (μg/mL)",
           title = "D-dimer Normalization During VTE Treatment") +
      theme_bw(base_size = 11)
    ggplotly(p)
  })

  ## ── Tab 5: Scenario Comparison ───────────────────────────────────────────
  all_scenario_sims <- reactive({
    scen_list <- input$scenarios
    if (length(scen_list) == 0) scen_list <- "riv_dvt"
    dur_h <- input$sim_days * 24

    dose_map <- list(
      riv_dvt  = c(ev(amt=15*0.93,cmt="RIV_GUT",ii=12,addl=41,time=0),
                   ev(amt=20*0.93,cmt="RIV_GUT",ii=24,addl=max(0,input$sim_days-22),time=21*24)),
      apix_pe  = c(ev(amt=10*0.50,cmt="APIX_CENT",ii=12,addl=13,time=0),
                   ev(amt=5*0.50,cmt="APIX_CENT",ii=12,addl=max(0,(input$sim_days-7)*2-1),time=7*24)),
      warf     = c(ev(amt=5,cmt="WARF_CENT",ii=24,addl=input$sim_days-1,time=0),
                   ev(amt=70*0.92,cmt="ENOX_CENT",ii=12,addl=19,time=0)),
      dabi     = ev(amt=150*0.065,cmt="DABI_CENT",ii=12,addl=input$sim_days*2-1),
      enox     = ev(amt=70*0.92,cmt="ENOX_CENT",ii=12,addl=input$sim_days*2-1),
      riv_ext  = ev(amt=10*0.93,cmt="RIV_GUT",ii=24,addl=input$sim_days-1)
    )
    name_map <- c(riv_dvt="Rivaroxaban DVT", apix_pe="Apixaban PE",
                   warf="Warfarin+Bridge", dabi="Dabigatran", enox="Enoxaparin",
                   riv_ext="Rivaroxaban 10 QD")

    purrr::map_dfr(scen_list, function(sc) {
      mod_shiny %>%
        param(BWT = input$bwt, eGFR_pat = input$egfr, CLOT_init = 80) %>%
        mrgsim(events = dose_map[[sc]], end = dur_h, delta = 2) %>%
        as.data.frame() %>%
        mutate(time_d = time / 24, scenario = name_map[sc])
    })
  })

  output$compare_plot <- renderPlotly({
    df <- all_scenario_sims()
    ep <- input$compare_endpoint
    p <- df %>%
      ggplot(aes_string(x = "time_d", y = ep, color = "scenario")) +
      geom_line(linewidth = 1) +
      scale_color_brewer(palette = "Set1") +
      labs(x = "Time (days)", y = ep, color = "Scenario",
           title = paste("Treatment Comparison:", ep)) +
      theme_bw(base_size = 11) +
      theme(legend.position = "bottom")
    ggplotly(p)
  })

  output$summary_table <- renderDT({
    df <- all_scenario_sims()
    df %>%
      group_by(scenario) %>%
      summarise(
        `D7 Clot (%)`  = round(CLOT_PCT[which.min(abs(time_d-7))],1),
        `D30 Clot (%)` = round(CLOT_PCT[which.min(abs(time_d-30))],1),
        `Peak FXa Inh (%)` = round(max(INH_FXa), 1),
        `Peak FIIa Inh (%)` = round(max(INH_FIIa), 1),
        `D7 D-dimer (ng/mL)` = round(DDIMER[which.min(abs(time_d-7))], 2),
        `D30 D-dimer` = round(DDIMER[which.min(abs(time_d-30))], 2),
        .groups = "drop"
      ) %>%
      datatable(rownames = FALSE, options = list(dom = "t"))
  })

  ## ── Tab 6: Biomarker Dashboard ────────────────────────────────────────────
  output$ib_inr_ss <- renderInfoBox({
    sim <- sim_primary()
    inr_ss <- tail(sim$INR, 50) %>% mean() %>% round(2)
    color <- if (inr_ss >= 2 && inr_ss <= 3) "green" else if (inr_ss < 2) "yellow" else "red"
    infoBox("Steady-State INR", inr_ss,
            icon = icon("balance-scale"), color = color,
            subtitle = "Target: 2.0–3.0")
  })
  output$ib_antixa_pk <- renderInfoBox({
    sim <- sim_primary()
    peak_xa <- round(max(sim$ANTI_XA), 2)
    color <- if (peak_xa >= 0.5 && peak_xa <= 1.0) "green" else "yellow"
    infoBox("Peak Anti-Xa", paste(peak_xa, "IU/mL"),
            icon = icon("shield"), color = color,
            subtitle = "Therapeutic: 0.5–1.0 IU/mL")
  })
  output$ib_ddimer_d7 <- renderInfoBox({
    sim <- sim_primary()
    d7_dd <- sim$DDIMER[which.min(abs(sim$time_d - 7))] %>% round(2)
    infoBox("D-dimer Day 7", paste(d7_dd, "μg/mL"),
            icon = icon("vial"), color = if (d7_dd > 0.5) "red" else "green",
            subtitle = "Normal <0.5 μg/mL")
  })
  output$ib_clot_d30 <- renderInfoBox({
    sim <- sim_primary()
    c30 <- sim$CLOT_PCT[which.min(abs(sim$time_d - 30))] %>% round(1)
    infoBox("Day 30 Residual Clot", paste0(c30, "%"),
            icon = icon("circle-notch"), color = if (c30 < 50) "green" else "yellow",
            subtitle = "Goal: <50% at 30 days")
  })

  output$bm_labs_plot <- renderPlotly({
    sim <- sim_primary()
    p <- sim %>%
      select(time_d, INR, ANTI_XA, DDIMER) %>%
      pivot_longer(-time_d) %>%
      ggplot(aes(x = time_d, y = value, color = name)) +
      geom_line(linewidth = 1) +
      facet_wrap(~name, scales = "free_y",
                 labeller = labeller(name = c(INR="INR", ANTI_XA="Anti-Xa (IU/mL)",
                                              DDIMER="D-dimer (μg/mL)"))) +
      scale_color_manual(values = c(INR="#8e24aa", ANTI_XA="#43a047", DDIMER="#e53935")) +
      labs(x = "Time (days)", y = "Value", color = "Lab Test",
           title = "Laboratory Monitoring Over Treatment Course") +
      theme_bw(base_size = 11) + theme(legend.position = "none")
    ggplotly(p)
  })

  output$bm_coag_plot <- renderPlotly({
    sim <- sim_primary()
    p <- sim %>%
      select(time_d, INH_FXa, INH_FIIa, TG_ETP) %>%
      pivot_longer(-time_d) %>%
      ggplot(aes(x = time_d, y = value, color = name)) +
      geom_line(linewidth = 1) +
      scale_color_manual(values = c(INH_FXa="#e53935", INH_FIIa="#8e24aa", TG_ETP="#1e88e5"),
                         labels = c("FXa Inhib (%)","FIIa Inhib (%)","Thrombin Gen (ETP %)")) +
      labs(x = "Time (days)", y = "% (vs. Baseline)",
           color = "Parameter",
           title = "Coagulation Cascade Suppression") +
      theme_bw(base_size = 11)
    ggplotly(p)
  })

  output$bm_ddimer_plot <- renderPlotly({
    sim <- sim_primary()
    p <- sim %>%
      ggplot(aes(x = time_d, y = DDIMER)) +
      geom_ribbon(aes(ymin = 0, ymax = DDIMER), fill = "#e53935", alpha = 0.3) +
      geom_line(color = "#e53935", linewidth = 1) +
      geom_hline(yintercept = 0.5, linetype = "dashed", color = "darkgreen", linewidth = 0.8) +
      annotate("text", x = input$sim_days * 0.5, y = 0.6,
               label = "Negative threshold (500 ng/mL)", color = "darkgreen", size = 3) +
      labs(x = "Time (days)", y = "D-dimer (μg/mL FEU)",
           title = "D-dimer Normalization Trajectory") +
      theme_bw(base_size = 11)
    ggplotly(p)
  })

  output$bm_interp_table <- renderDT({
    tibble(
      Biomarker    = c("D-dimer","INR","Anti-Xa (LMWH)","aPTT","F1+2 fragment","TAT complex"),
      `Normal Range`=c("<500 ng/mL","0.8–1.2","N/A","25–37s","<1.1 nM","<4 μg/L"),
      `Therapeutic Target`=c(">500 at VTE; <500=exclude","2.0–3.0 (warfarin)","0.5–1.0 IU/mL","1.5–2.5×(UFH)","↓ with anticoag","↓ with anticoag"),
      Interpretation=c("↑ in VTE; ↓ with treatment","Warfarin efficacy monitor","LMWH efficacy","UFH efficacy","Thrombin generation","Active thrombin")
    ) %>%
      datatable(rownames = FALSE, options = list(dom = "t", pageLength = 10))
  })

  ## ── Tab 1: Risk Assessment ────────────────────────────────────────────────
  output$wells_dvt <- renderText({
    pts <- 0
    if ("surg" %in% input$risk_stasis || "imm" %in% input$risk_stasis) pts <- pts + 1
    if ("cvc" %in% input$risk_endo) pts <- pts + 1
    if ("prior" %in% input$risk_endo) pts <- pts + 1
    cat_str <- if (pts >= 2) "High Probability (≥2 pts)" else
               if (pts == 1) "Moderate (1 pt)" else "Low (0 pts)"
    paste0("Score: ", pts, " — ", cat_str, "\n",
           "Recommendation: ", if (pts >= 2) "Compression ultrasound + anticoag"
                               else "D-dimer first; if ↑ → imaging")
  })

  output$wells_pe <- renderText({
    pts <- 0
    if ("prior" %in% input$risk_endo) pts <- pts + 1.5
    if ("surg" %in% input$risk_stasis) pts <- pts + 1.5
    cat_str <- if (pts > 4) "High (>4 pts)" else if (pts > 0) "Intermediate (1-4 pts)" else "Low (0 pts)"
    paste0("Score: ", pts, " — ", cat_str, "\n",
           if (pts > 4) "CTPA directly" else "D-dimer → if ↑ → CTPA")
  })

  output$anticoag_rec <- renderText({
    renal_ok <- input$egfr >= 30
    "Standard DOAC preferred (DOAC preferred over warfarin per ESC 2019)\nRivaroxaban or Apixaban — no need for bridging"
  })

  output$renal_summary <- renderText({
    egfr <- input$egfr
    stage <- if (egfr >= 90) "G1 (Normal)" else if (egfr >= 60) "G2 (Mildly ↓)" else
             if (egfr >= 30) "G3 (Moderately ↓)" else if (egfr >= 15) "G4 (Severely ↓)" else "G5 (Failure)"
    note <- if (egfr < 30) "Avoid dabigatran; use warfarin or reduced-dose LMWH" else
            if (egfr < 50) "Reduce dabigatran dose; use anti-Xa DOACs with caution" else "Standard dosing acceptable"
    paste0("eGFR: ", egfr, " mL/min/1.73m²\nCKD Stage: ", stage, "\n", note)
  })

  output$duration_rec <- renderText({
    vte    <- input$vte_type
    prov   <- input$vte_provoked
    dur <- if (prov == "provoked") "3 months (minimum)" else
           if (prov == "cancer")  "Indefinite (cancer-associated → LMWH or rivaroxaban)" else
                                   "Minimum 3-6 months; consider indefinite if unprovoked"
    paste0(dur, "\n",
           "Reassess bleeding risk at 3 months\n",
           if (prov == "unprovoked") "(D-dimer-guided extension: HERDOO2/Vienna Prediction Score)" else "")
  })

  output$drug_rec <- renderText({
    contr <- input$contraind
    if ("valve" %in% contr) {
      "Mechanical valve: Warfarin ONLY (DOACs contraindicated)"
    } else if ("preg_warf" %in% contr) {
      "Pregnancy: LMWH (enoxaparin 1 mg/kg BID SC)\nWarfarin contraindicated in 1st trimester"
    } else if ("aps_warf" %in% contr) {
      "Triple-positive APS: Warfarin preferred (INR 2-3)\n(DOACs ↑ recurrence in triple-positive APS)"
    } else if ("renal_sev" %in% contr) {
      "Severe CKD (eGFR<15): Warfarin preferred\nUnfractionated heparin (IV) for acute"
    } else {
      "First-line: Rivaroxaban (15mg BID×21d→20mg QD) OR\nApixaban (10mg BID×7d→5mg BID)\nConvenient: No bridging needed"
    }
  })
}

shinyApp(ui = ui, server = server)
