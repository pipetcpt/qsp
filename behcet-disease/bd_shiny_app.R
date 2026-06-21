## ============================================================
## Behcet's Disease (BD) — QSP Interactive Shiny App
## ============================================================
## Tabs:
##   1. Patient Profile & Disease Overview
##   2. Drug Pharmacokinetics (PK)
##   3. Cytokine & Immune Cell Dynamics
##   4. Clinical Endpoints & Organ Activity
##   5. Treatment Scenario Comparison
##   6. Biomarker & HLA-B51 Stratification
##   7. Dose-Response Optimization
## ============================================================

library(shiny)
library(shinydashboard)
library(mrgsolve)
library(dplyr)
library(tidyr)
library(ggplot2)
library(plotly)
library(DT)
library(patchwork)

## ============================================================
## Embedded mrgsolve Model
## ============================================================
bd_code <- '
$PARAM @annotated
ka_col   : 1.2   : Colchicine ka (1/h)
Foral_col: 0.44  : Colchicine F
CL_col   : 17.0  : Colchicine CL (L/h)
V1_col   : 28.0  : Colchicine V1 (L)
V2_col   : 5500  : Colchicine V2 (L)
Q_col    : 30.0  : Colchicine Q (L/h)
CL_pred  : 8.4   : Prednisolone CL (L/h)
V1_pred  : 22.0  : Prednisolone V1 (L)
V2_pred  : 50.0  : Prednisolone V2 (L)
Q_pred   : 10.0  : Prednisolone Q (L/h)
CL_ada   : 0.012 : Adalimumab CL (L/h)
V1_ada   : 2.8   : Adalimumab V1 (L)
V2_ada   : 3.4   : Adalimumab V2 (L)
Q_ada    : 0.003 : Adalimumab Q (L/h)
ka_apr   : 0.58  : Apremilast ka (1/h)
CL_apr   : 10.0  : Apremilast CL (L/h)
V1_apr   : 87.0  : Apremilast V1 (L)
CL_can   : 0.007 : Canakinumab CL (L/h)
V1_can   : 3.0   : Canakinumab V1 (L)
V2_can   : 3.2   : Canakinumab V2 (L)
Q_can    : 0.0015: Canakinumab Q (L/h)
kNEU_in  : 0.05  : NEU synthesis (1/h)
kNEU_out : 0.05  : NEU decay (1/h)
kTH1_in  : 0.03  : TH1 synthesis (1/h)
kTH1_out : 0.03  : TH1 decay (1/h)
kTH17_in : 0.03  : TH17 synthesis (1/h)
kTH17_out: 0.03  : TH17 decay (1/h)
kTREG_in : 0.02  : TREG synthesis (1/h)
kTREG_out: 0.02  : TREG decay (1/h)
kTNFA_syn: 0.08  : TNFa synthesis (1/h)
kTNFA_deg: 0.08  : TNFa decay (1/h)
kIL1B_syn: 0.06  : IL1B synthesis (1/h)
kIL1B_deg: 0.06  : IL1B decay (1/h)
kIL6_syn : 0.07  : IL6 synthesis (1/h)
kIL6_deg : 0.07  : IL6 decay (1/h)
kIL17_syn: 0.05  : IL17A synthesis (1/h)
kIL17_deg: 0.05  : IL17A decay (1/h)
kEA_on   : 0.04  : EA onset (1/h)
kEA_off  : 0.04  : EA offset (1/h)
kOUL_on  : 0.02  : OUL onset (1/h)
kOUL_off : 0.01  : OUL heal (1/h)
kOCI_on  : 0.015 : OCI onset (1/h)
kOCI_off : 0.008 : OCI heal (1/h)
a_TNFA_TH1 : 0.3  : TNFa→TH1 amplify
a_TNFA_NEU : 0.4  : TNFa→NEU prime
a_IL17_NEU : 0.3  : IL17→NEU recruit
a_TH17_IL17: 0.5  : TH17→IL17A
a_TH1_TNFA : 0.4  : TH1→TNFa
a_NEU_TNFA : 0.3  : NEU→TNFa
a_NEU_IL1B : 0.3  : NEU→IL1B
a_IL1B_IL6 : 0.3  : IL1B→IL6
a_TNFA_IL6 : 0.2  : TNFa→IL6
a_IL6_TH17 : 0.2  : IL6→TH17
a_TREG_inh : 0.3  : TREG suppress TH17
IC50_col_NEU : 0.8  : Colch IC50 NEU (ng/mL)
IC50_col_IL1 : 1.5  : Colch IC50 IL1B (ng/mL)
EC50_pred_TNF: 50.0 : Pred EC50 TNFa (ng/mL)
EC50_pred_IL6: 30.0 : Pred EC50 IL6 (ng/mL)
Emax_pred    : 0.80 : Pred Emax
EC50_ada_TNF : 1500 : Ada EC50 TNFa (ng/mL)
Emax_ada     : 0.90 : Ada Emax
EC50_apr_TNF : 200  : Aprem EC50 TNFa (ng/mL)
EC50_apr_IL17: 300  : Aprem EC50 IL17 (ng/mL)
Emax_apr     : 0.65 : Aprem Emax
EC50_can_IL1 : 800  : Canaki EC50 IL1B (ng/mL)
Emax_can     : 0.92 : Canaki Emax
HLAB51_factor: 1.4  : HLA-B51 disease severity
$CMT @annotated
AGUT_COL : Colchicine gut (mg)
ACOL     : Colchicine central (mg)
ACOL_T   : Colchicine tissue (mg)
APRED    : Prednisolone central (mg)
APRED_T  : Prednisolone tissue (mg)
AADA     : Adalimumab central (mg)
AADA_P   : Adalimumab peripheral (mg)
AAPR     : Apremilast central (mg)
ACAN     : Canakinumab central (mg)
ACAN_P   : Canakinumab peripheral (mg)
NEU      : Neutrophil activation
TH1      : Th1 cells
TH17     : Th17 cells
TREG     : Treg cells
TNFA     : TNF-alpha
IL1B     : IL-1beta
IL6C     : IL-6
IL17A    : IL-17A
EA       : Endothelial activation
OUL      : Oral ulcer index
OCI      : Ocular inflammation index
BDCAF    : BDCAF composite score
$MAIN
double Cp_col  = ACOL  / V1_col  * 1000;
double Cp_pred = APRED / V1_pred * 1000;
double Cp_ada  = AADA  / V1_ada  * 1000;
double Cp_apr  = AAPR  / V1_apr  * 1000;
double Cp_can  = ACAN  / V1_can  * 1000;
double Ecol_NEU = Cp_col / (IC50_col_NEU + Cp_col);
double Ecol_IL1 = Cp_col / (IC50_col_IL1 + Cp_col);
double Epred_TNF = Emax_pred * (Cp_pred / (EC50_pred_TNF + Cp_pred));
double Epred_IL6 = Emax_pred * (Cp_pred / (EC50_pred_IL6 + Cp_pred));
double Eada_TNF  = Emax_ada  * (Cp_ada  / (EC50_ada_TNF  + Cp_ada));
double Eapr_TNF  = Emax_apr  * (Cp_apr  / (EC50_apr_TNF  + Cp_apr));
double Eapr_IL17 = Emax_apr  * (Cp_apr  / (EC50_apr_IL17 + Cp_apr));
double Ecan_IL1  = Emax_can  * (Cp_can  / (EC50_can_IL1  + Cp_can));
double DIS = HLAB51_factor;
double TNFA_syn = kTNFA_syn * DIS * (1.0 + a_TH1_TNFA*(TH1-1.0) + a_NEU_TNFA*(NEU-1.0));
double IL1B_syn = kIL1B_syn * DIS * (1.0 + a_NEU_IL1B*(NEU-1.0));
double IL6_syn  = kIL6_syn  * DIS * (1.0 + a_IL1B_IL6*(IL1B-1.0) + a_TNFA_IL6*(TNFA-1.0));
double IL17_syn = kIL17_syn * DIS * (1.0 + a_TH17_IL17*(TH17-1.0));
double NEU_syn  = kNEU_in   * DIS * (1.0 + a_TNFA_NEU*(TNFA-1.0) + a_IL17_NEU*(IL17A-1.0));
double TH1_syn  = kTH1_in   * DIS * (1.0 + a_TNFA_TH1*(TNFA-1.0));
double TH17_syn = kTH17_in  * DIS * (1.0 + a_IL6_TH17*(IL6C-1.0)) / (1.0 + a_TREG_inh*TREG);
double TREG_syn = kTREG_in  / (1.0 + 0.2*(TNFA-1.0));
double EA_on_rate  = kEA_on  * DIS * (TNFA + IL1B + IL17A) / 3.0;
double OUL_on_rate = kOUL_on * DIS * (IL17A*0.4 + TNFA*0.4 + NEU*0.2);
double OCI_on_rate = kOCI_on * DIS * EA;
$ODE
dxdt_AGUT_COL = -ka_col * AGUT_COL;
dxdt_ACOL     =  ka_col * AGUT_COL * Foral_col - (CL_col+Q_col)/V1_col*ACOL + Q_col/V2_col*ACOL_T;
dxdt_ACOL_T   =  Q_col/V1_col*ACOL - Q_col/V2_col*ACOL_T;
dxdt_APRED    = -(CL_pred+Q_pred)/V1_pred*APRED + Q_pred/V2_pred*APRED_T;
dxdt_APRED_T  =  Q_pred/V1_pred*APRED - Q_pred/V2_pred*APRED_T;
dxdt_AADA     = -(CL_ada+Q_ada)/V1_ada*AADA + Q_ada/V2_ada*AADA_P;
dxdt_AADA_P   =  Q_ada/V1_ada*AADA - Q_ada/V2_ada*AADA_P;
dxdt_AAPR     = -CL_apr/V1_apr*AAPR;
dxdt_ACAN     = -(CL_can+Q_can)/V1_can*ACAN + Q_can/V2_can*ACAN_P;
dxdt_ACAN_P   =  Q_can/V1_can*ACAN - Q_can/V2_can*ACAN_P;
dxdt_NEU  = NEU_syn  - kNEU_out  * NEU  * (1.0 + Ecol_NEU);
dxdt_TH1  = TH1_syn  - kTH1_out  * TH1  * (1.0 + Epred_TNF*0.5);
dxdt_TH17 = TH17_syn - kTH17_out * TH17 * (1.0 + Epred_TNF*0.3 + Eapr_IL17*0.4);
dxdt_TREG = TREG_syn - kTREG_out * TREG;
dxdt_TNFA = TNFA_syn * (1.0-Eada_TNF) * (1.0-Epred_TNF) * (1.0-Eapr_TNF*0.6) - kTNFA_deg*TNFA;
dxdt_IL1B = IL1B_syn * (1.0-Ecol_IL1) * (1.0-Ecan_IL1)  * (1.0-Epred_TNF*0.4) - kIL1B_deg*IL1B;
dxdt_IL6C = IL6_syn  * (1.0-Epred_IL6) - kIL6_deg*IL6C;
dxdt_IL17A = IL17_syn * (1.0-Eapr_IL17) * (1.0-Epred_TNF*0.3) - kIL17_deg*IL17A;
dxdt_EA   = EA_on_rate  * (1.0-Eada_TNF*0.5-Epred_TNF*0.3) - kEA_off*EA;
dxdt_OUL  = OUL_on_rate * (1.0-Epred_TNF*0.4-Eada_TNF*0.4-Eapr_TNF*0.3) - kOUL_off*OUL;
dxdt_OCI  = OCI_on_rate * (1.0-Eada_TNF*0.6-Epred_TNF*0.5) - kOCI_off*OCI;
dxdt_BDCAF = 0.25*OUL + 0.25*OCI + 0.25*EA + 0.25*(TNFA+IL1B+IL6C+IL17A)/4.0 - 0.1*BDCAF;
$TABLE
double Cp_COL  = ACOL  / V1_col  * 1000;
double Cp_PRED = APRED / V1_pred * 1000;
double Cp_ADA  = AADA  / V1_ada  * 1000;
double Cp_APR  = AAPR  / V1_apr  * 1000;
double Cp_CAN  = ACAN  / V1_can  * 1000;
$CAPTURE @annotated
Cp_COL  : Colchicine Cp (ng/mL)
Cp_PRED : Prednisolone Cp (ng/mL)
Cp_ADA  : Adalimumab Cp (ng/mL)
Cp_APR  : Apremilast Cp (ng/mL)
Cp_CAN  : Canakinumab Cp (ng/mL)
NEU : Neutrophil activation
TH1 : Th1 cells
TH17 : Th17 cells
TREG : Treg cells
TNFA : TNF-alpha
IL1B : IL-1beta
IL6C : IL-6
IL17A : IL-17A
EA : Endothelial activation
OUL : Oral ulcer index
OCI : Ocular inflammation
BDCAF : BDCAF score
'

bd_mod <- mcode("BehcetDisease_Shiny", bd_code, quiet = TRUE)

## ============================================================
## Helper Functions
## ============================================================
make_init <- function(disease_severity = "Moderate") {
  sev <- switch(disease_severity,
    "Mild"   = list(NEU=1.5, TH1=1.5, TH17=1.8, TREG=0.8, TNFA=2.0, IL1B=1.5, IL6C=1.5, IL17A=2.0, EA=1.5, OUL=2.0, OCI=1.0, BDCAF=4.0),
    "Moderate" = list(NEU=2.5, TH1=2.0, TH17=2.8, TREG=0.5, TNFA=3.0, IL1B=2.5, IL6C=2.0, IL17A=2.8, EA=2.0, OUL=3.0, OCI=2.0, BDCAF=8.0),
    "Severe"   = list(NEU=4.0, TH1=3.5, TH17=4.5, TREG=0.3, TNFA=5.0, IL1B=4.0, IL6C=3.5, IL17A=4.5, EA=3.5, OUL=5.0, OCI=4.0, BDCAF=12.0)
  )
  base <- list(AGUT_COL=0, ACOL=0, ACOL_T=0, APRED=0, APRED_T=0,
               AADA=0, AADA_P=0, AAPR=0, ACAN=0, ACAN_P=0)
  c(base, sev)
}

run_bd_sim <- function(mod, init_vals, dose_events, hlab51_pos, sim_days = 90) {
  sim_end <- sim_days * 24
  mod2 <- mod %>%
    init(init_vals) %>%
    param(HLAB51_factor = ifelse(hlab51_pos, 1.4, 1.0))

  tryCatch({
    if (is.null(dose_events)) {
      sim <- mrgsim(mod2, end = sim_end, delta = 4)
    } else {
      sim <- mrgsim(mod2, events = dose_events, end = sim_end, delta = 4)
    }
    as.data.frame(sim) %>% mutate(time_days = time / 24)
  }, error = function(e) {
    data.frame(time = 0, time_days = 0, BDCAF = NA, OUL = NA, OCI = NA,
               TNFA = NA, IL1B = NA, IL6C = NA, IL17A = NA, NEU = NA,
               TH1 = NA, TH17 = NA, TREG = NA, EA = NA,
               Cp_COL = NA, Cp_PRED = NA, Cp_ADA = NA, Cp_APR = NA, Cp_CAN = NA)
  })
}

build_events <- function(drug, dose, freq_days, duration_days) {
  cmt_map <- c(colchicine = 1, prednisolone = 4, adalimumab = 6, apremilast = 8, canakinumab = 9)
  cmt_id <- cmt_map[drug]
  if (is.na(cmt_id)) return(NULL)
  times <- seq(0, (duration_days - 1) * 24, by = freq_days * 24)
  ev(data.frame(time = times, cmt = cmt_id, amt = dose, evid = 1))
}

## ============================================================
## UI
## ============================================================
ui <- dashboardPage(
  skin = "red",
  dashboardHeader(title = "Behcet's Disease QSP"),
  dashboardSidebar(
    sidebarMenu(
      menuItem("Patient Profile",         tabName = "tab_profile",   icon = icon("user-md")),
      menuItem("Drug PK",                 tabName = "tab_pk",        icon = icon("pills")),
      menuItem("Cytokine & Immune Dynamics", tabName = "tab_cytokine", icon = icon("microscope")),
      menuItem("Clinical Endpoints",      tabName = "tab_endpoints", icon = icon("stethoscope")),
      menuItem("Scenario Comparison",     tabName = "tab_compare",   icon = icon("chart-line")),
      menuItem("Biomarker & HLA-B51",     tabName = "tab_biomarker", icon = icon("dna")),
      menuItem("Dose Optimization",       tabName = "tab_dose",      icon = icon("sliders-h"))
    )
  ),
  dashboardBody(
    tags$head(tags$style(HTML("
      .skin-red .main-header .navbar { background-color: #c1121f; }
      .skin-red .main-header .logo  { background-color: #780000; }
      .content-wrapper { background-color: #f5f5f5; }
      .box { border-radius: 8px; }
    "))),
    tabItems(

      ## ======= TAB 1: Patient Profile =======
      tabItem(tabName = "tab_profile",
        fluidRow(
          box(title = "Patient Demographics & Disease Characteristics", width = 4, status = "danger",
            selectInput("severity",  "Disease Severity",
                        choices = c("Mild", "Moderate", "Severe"), selected = "Moderate"),
            checkboxInput("hlab51", "HLA-B51 Positive", value = TRUE),
            selectInput("sex",       "Sex",     choices = c("Male", "Female")),
            numericInput("age",      "Age (years)", value = 32, min = 15, max = 70),
            numericInput("weight",   "Weight (kg)",  value = 65, min = 40, max = 120),
            checkboxGroupInput("organs", "Organs Affected",
                               choices = c("Oral mucosa", "Genital", "Ocular", "Skin",
                                           "Vascular", "Neurological", "GI"),
                               selected = c("Oral mucosa", "Ocular", "Skin"))
          ),
          box(title = "Disease Description", width = 8, status = "danger",
            HTML('
              <h4 style="color:#c1121f;">Behcet\'s Disease (BD)</h4>
              <p>Behcet\'s Disease is a systemic, multifocal inflammatory vasculitis of unknown etiology,
              primarily affecting individuals from the "Silk Road" geographic belt (Turkey, Middle East, East Asia).
              Characterized by recurrent oral aphthous ulcers (hallmark), genital ulcers, ocular inflammation
              (uveitis/retinal vasculitis), skin lesions, and vascular/neurological manifestations.</p>
              <hr/>
              <table class="table table-bordered table-sm">
                <tr><th>Feature</th><th>Detail</th></tr>
                <tr><td>Prevalence</td><td>~20-370/100,000 (Turkey), ~5/100,000 (Japan)</td></tr>
                <tr><td>Peak Onset</td><td>20-40 years (both sexes)</td></tr>
                <tr><td>Genetic Risk</td><td>HLA-B51 (OR 5-6x), IL-10, STAT4 polymorphisms</td></tr>
                <tr><td>Key Pathomechanism</td><td>Neutrophil hyperactivation + Th1/Th17 dysregulation</td></tr>
                <tr><td>Diagnosis</td><td>ISG criteria (oral ulcers + 2 of: genital, ocular, skin, pathergy)</td></tr>
                <tr><td>Mortality</td><td>Low overall; vascular/neurological BD carries significant morbidity</td></tr>
                <tr><td>QoL Impact</td><td>Ocular BD → blindness 25% untreated; frequent relapses</td></tr>
              </table>
              <h5 style="color:#c1121f; margin-top:12px;">Immunopathology at a Glance</h5>
              <p>
                <b>Innate:</b> Neutrophil hyperactivation (HLA-B51 dependent), NLRP3 inflammasome activation (→IL-1β),
                NETs formation, monocyte/macrophage overactivation.<br/>
                <b>Adaptive:</b> Th1 (IFN-γ, TNF-α) and Th17 (IL-17A) predominance; Treg dysfunction; CD8+ CTL damage.<br/>
                <b>Vascular:</b> Endothelial cell activation, adhesion molecule upregulation, pro-coagulant shift →
                venous/arterial thrombosis, aneurysm.
              </p>
            ')
          )
        ),
        fluidRow(
          box(title = "Behcet's Disease Activity (BDCAF) Scale", width = 12, status = "warning",
            HTML('<table class="table table-condensed table-hover" style="font-size:11px;">
              <tr style="background:#c1121f; color:white;">
                <th>BDCAF Component</th><th>Score Range</th><th>Assessment</th>
              </tr>
              <tr><td>Oral ulcers</td><td>0-3</td><td>None / Mild / Moderate / Severe</td></tr>
              <tr><td>Genital ulcers</td><td>0-2</td><td>None / Present / Severe</td></tr>
              <tr><td>Ocular involvement</td><td>0-3</td><td>None / Anterior UV / Post UV / Retinal vasculitis</td></tr>
              <tr><td>Skin lesions</td><td>0-3</td><td>None / Pseudofolliculitis / EN / Both</td></tr>
              <tr><td>Neurological</td><td>0-2</td><td>None / Peripheral / CNS</td></tr>
              <tr><td>Arthritis</td><td>0-1</td><td>None / Present</td></tr>
              <tr><td>Vascular</td><td>0-2</td><td>None / Venous / Arterial</td></tr>
              <tr style="background:#fdf3cd;"><td><b>Total BDCAF</b></td><td><b>0-12</b></td>
              <td><b>Remission &lt;3 | Active 3-6 | Severe ≥7</b></td></tr>
            </table>')
          )
        )
      ),

      ## ======= TAB 2: Drug PK =======
      tabItem(tabName = "tab_pk",
        fluidRow(
          box(title = "Drug Selection & Dosing", width = 4, status = "primary",
            selectInput("pk_drug",  "Drug",
                        choices = c("Colchicine", "Prednisolone", "Adalimumab (anti-TNF)",
                                    "Apremilast (PDE4i)", "Canakinumab (anti-IL-1B)"),
                        selected = "Adalimumab (anti-TNF)"),
            numericInput("pk_dose",  "Dose (mg)", value = 40,  min = 1,   max = 1000),
            numericInput("pk_freq",  "Frequency (days)", value = 14, min = 0.5, max = 90),
            numericInput("pk_days",  "Duration (days)",  value = 90, min = 7,   max = 365),
            selectInput("pk_sev",   "Disease Severity",
                        choices = c("Mild", "Moderate", "Severe"), selected = "Moderate"),
            checkboxInput("pk_hlab51", "HLA-B51 Positive", value = TRUE),
            actionButton("run_pk", "Run PK Simulation", class = "btn-primary btn-block")
          ),
          box(title = "PK Reference Values", width = 8, status = "primary",
            HTML('<table class="table table-sm table-hover">
              <tr style="background:#003049; color:white;">
                <th>Drug</th><th>Route</th><th>F (%)</th><th>t1/2</th>
                <th>CL (L/h)</th><th>Vd (L)</th><th>Cmax Target</th>
              </tr>
              <tr><td>Colchicine</td><td>Oral</td><td>44</td><td>27 h</td>
                  <td>17.0</td><td>588</td><td>2-4 ng/mL</td></tr>
              <tr><td>Prednisolone</td><td>Oral</td><td>80</td><td>3.5 h</td>
                  <td>8.4</td><td>22-70</td><td>100-1000 ng/mL</td></tr>
              <tr><td>Adalimumab</td><td>SC Q2W</td><td>64</td><td>~2 wk</td>
                  <td>0.012</td><td>6.1</td><td>5000-8000 ng/mL</td></tr>
              <tr><td>Apremilast</td><td>Oral BID</td><td>73</td><td>9 h</td>
                  <td>10.0</td><td>87</td><td>500-800 ng/mL</td></tr>
              <tr><td>Canakinumab</td><td>SC Q8W</td><td>66</td><td>26 d</td>
                  <td>0.007</td><td>6.2</td><td>1000-3000 ng/mL</td></tr>
            </table>')
          )
        ),
        fluidRow(
          box(title = "Concentration-Time Profile",        width = 8, status = "info",
            plotlyOutput("pk_plot", height = "380px")),
          box(title = "PK Summary Statistics",             width = 4, status = "info",
            DTOutput("pk_summary_tbl"))
        )
      ),

      ## ======= TAB 3: Cytokine & Immune Dynamics =======
      tabItem(tabName = "tab_cytokine",
        fluidRow(
          box(title = "Treatment Parameters", width = 3, status = "warning",
            selectInput("cyt_drug", "Drug",
                        choices = c("None (untreated)",
                                    "Colchicine 0.5 mg BID",
                                    "Prednisolone 40 mg/day",
                                    "Adalimumab 40 mg Q2W",
                                    "Apremilast 30 mg BID",
                                    "Canakinumab 150 mg Q8W"),
                        selected = "Adalimumab 40 mg Q2W"),
            selectInput("cyt_sev", "Baseline Disease Severity",
                        choices = c("Mild", "Moderate", "Severe"), selected = "Moderate"),
            checkboxInput("cyt_hlab51", "HLA-B51 Positive", value = TRUE),
            numericInput("cyt_days", "Simulation Duration (days)", value = 90, min = 14, max = 180),
            actionButton("run_cyt", "Run Simulation", class = "btn-warning btn-block")
          ),
          box(title = "Cytokine Dynamics", width = 9, status = "warning",
            plotlyOutput("cyt_cytokine_plot", height = "400px"))
        ),
        fluidRow(
          box(title = "Immune Cell Dynamics",              width = 6, status = "success",
            plotlyOutput("cyt_cells_plot",     height = "320px")),
          box(title = "Th17/Treg Balance (Key for BD)",   width = 6, status = "success",
            plotlyOutput("cyt_th17treg_plot",  height = "320px"))
        )
      ),

      ## ======= TAB 4: Clinical Endpoints =======
      tabItem(tabName = "tab_endpoints",
        fluidRow(
          box(title = "Treatment Settings", width = 3, status = "danger",
            selectInput("ep_drug", "Primary Treatment",
                        choices = c("None (untreated)",
                                    "Colchicine 0.5 mg BID",
                                    "Prednisolone 40 mg/day",
                                    "Adalimumab 40 mg Q2W",
                                    "Apremilast 30 mg BID",
                                    "Canakinumab 150 mg Q8W",
                                    "Colchicine + Prednisolone",
                                    "Adalimumab + Apremilast"),
                        selected = "Adalimumab 40 mg Q2W"),
            selectInput("ep_sev",  "Disease Severity",
                        choices = c("Mild", "Moderate", "Severe"), selected = "Moderate"),
            checkboxInput("ep_hlab51", "HLA-B51 Positive", value = TRUE),
            numericInput("ep_days", "Duration (days)", value = 90, min = 14, max = 365),
            actionButton("run_ep", "Simulate", class = "btn-danger btn-block")
          ),
          box(title = "BDCAF Disease Activity Score", width = 9, status = "danger",
            plotlyOutput("ep_bdcaf_plot", height = "400px"))
        ),
        fluidRow(
          box(title = "Oral Ulcer Activity",              width = 4, status = "warning",
            plotlyOutput("ep_oral_plot",    height = "280px")),
          box(title = "Ocular Inflammation Score",        width = 4, status = "warning",
            plotlyOutput("ep_ocular_plot",  height = "280px")),
          box(title = "Endothelial Activation (Vascular Risk)", width = 4, status = "warning",
            plotlyOutput("ep_ea_plot",      height = "280px"))
        ),
        fluidRow(
          box(title = "Day-90 Endpoint Summary", width = 12, status = "primary",
            DTOutput("ep_summary_tbl"))
        )
      ),

      ## ======= TAB 5: Scenario Comparison =======
      tabItem(tabName = "tab_compare",
        fluidRow(
          box(title = "Comparison Settings", width = 3, status = "info",
            selectInput("cmp_sev",  "Disease Severity",
                        choices = c("Mild", "Moderate", "Severe"), selected = "Moderate"),
            checkboxInput("cmp_hlab51", "HLA-B51 Positive", value = TRUE),
            numericInput("cmp_days", "Simulation Duration (days)", value = 90, min = 14, max = 180),
            checkboxGroupInput("cmp_scenarios", "Treatment Scenarios to Compare",
                               choices = c(
                                 "Untreated",
                                 "Colchicine",
                                 "Prednisolone",
                                 "Adalimumab (anti-TNF)",
                                 "Apremilast (PDE4i)",
                                 "Canakinumab (anti-IL-1B)",
                                 "Colchicine + Prednisolone",
                                 "Adalimumab + Apremilast (Refractory)"
                               ),
                               selected = c("Untreated", "Colchicine", "Adalimumab (anti-TNF)",
                                            "Apremilast (PDE4i)", "Canakinumab (anti-IL-1B)")),
            actionButton("run_cmp", "Compare All", class = "btn-info btn-block")
          ),
          box(title = "BDCAF Comparison — All Scenarios", width = 9, status = "info",
            plotlyOutput("cmp_bdcaf_plot", height = "400px"))
        ),
        fluidRow(
          box(title = "Oral Ulcer Comparison",            width = 6, status = "primary",
            plotlyOutput("cmp_oral_plot",  height = "300px")),
          box(title = "Ocular Inflammation Comparison",   width = 6, status = "primary",
            plotlyOutput("cmp_ocular_plot", height = "300px"))
        ),
        fluidRow(
          box(title = "Day-90 Endpoint Summary Table", width = 12, status = "warning",
            DTOutput("cmp_summary_tbl"))
        )
      ),

      ## ======= TAB 6: Biomarker & HLA-B51 =======
      tabItem(tabName = "tab_biomarker",
        fluidRow(
          box(title = "Biomarker Stratification Settings", width = 3, status = "success",
            selectInput("bio_drug", "Treatment",
                        choices = c("Adalimumab (anti-TNF)", "Canakinumab (anti-IL-1B)",
                                    "Apremilast (PDE4i)", "Colchicine"),
                        selected = "Adalimumab (anti-TNF)"),
            numericInput("bio_days", "Duration (days)", value = 90, min = 30, max = 180),
            selectInput("bio_sev",   "Disease Severity",
                        choices = c("Mild", "Moderate", "Severe"), selected = "Moderate"),
            hr(),
            h5("HLA-B51 Impact Analysis"),
            p("Compare HLA-B51 positive vs negative patient response to the same treatment."),
            actionButton("run_bio", "Analyze Biomarkers", class = "btn-success btn-block")
          ),
          box(title = "HLA-B51 Stratified BDCAF Response",  width = 9, status = "success",
            plotlyOutput("bio_hlab51_plot", height = "350px"))
        ),
        fluidRow(
          box(title = "Serum Cytokine Biomarkers (TNF-α, IL-1β, IL-6, IL-17A)", width = 6, status = "primary",
            plotlyOutput("bio_cytokine_bio_plot", height = "330px")),
          box(title = "Time to BDCAF < 3 (Remission) by HLA-B51 Status", width = 6, status = "warning",
            plotlyOutput("bio_remission_plot", height = "330px"))
        ),
        fluidRow(
          box(title = "Biomarker Clinical Reference Ranges", width = 12, status = "danger",
            HTML('<table class="table table-bordered table-sm">
              <tr style="background:#1b4332; color:white;">
                <th>Biomarker</th><th>Normal</th><th>Mild BD</th>
                <th>Active BD</th><th>Severe BD</th><th>Interpretation</th>
              </tr>
              <tr><td>TNF-α (pg/mL)</td><td>&lt;8</td><td>10-25</td><td>30-80</td><td>&gt;100</td>
                  <td>Anti-TNF target; correlates with BDCAF</td></tr>
              <tr><td>IL-1β (pg/mL)</td><td>&lt;5</td><td>8-15</td><td>20-60</td><td>&gt;80</td>
                  <td>Canakinumab/anakinra target</td></tr>
              <tr><td>IL-6 (pg/mL)</td><td>&lt;7</td><td>10-20</td><td>25-70</td><td>&gt;100</td>
                  <td>Drives APR (CRP, SAA)</td></tr>
              <tr><td>IL-17A (pg/mL)</td><td>&lt;20</td><td>25-50</td><td>60-150</td><td>&gt;200</td>
                  <td>Apremilast/secukinumab target</td></tr>
              <tr><td>CRP (mg/L)</td><td>&lt;3</td><td>5-15</td><td>20-80</td><td>&gt;100</td>
                  <td>General inflammation marker</td></tr>
              <tr><td>Neutrophil:Lymphocyte Ratio</td><td>1-3</td><td>3-5</td><td>5-10</td><td>&gt;10</td>
                  <td>Neutrophil hyperactivation proxy</td></tr>
            </table>')
          )
        )
      ),

      ## ======= TAB 7: Dose Optimization =======
      tabItem(tabName = "tab_dose",
        fluidRow(
          box(title = "Dose-Response Analysis Settings", width = 3, status = "purple",
            selectInput("dr_drug", "Drug for Dose-Response",
                        choices = c("Adalimumab (anti-TNF)", "Canakinumab (anti-IL-1B)",
                                    "Apremilast (PDE4i)", "Prednisolone"),
                        selected = "Adalimumab (anti-TNF)"),
            sliderInput("dr_dose_min", "Min Dose (mg)", min = 1, max = 50, value = 10),
            sliderInput("dr_dose_max", "Max Dose (mg)", min = 20, max = 500, value = 80),
            numericInput("dr_n_doses", "Number of Doses to Test", value = 8, min = 4, max = 20),
            selectInput("dr_sev",  "Disease Severity",
                        choices = c("Mild", "Moderate", "Severe"), selected = "Moderate"),
            checkboxInput("dr_hlab51", "HLA-B51 Positive", value = TRUE),
            numericInput("dr_days",  "Duration (days)", value = 90, min = 30, max = 180),
            actionButton("run_dr", "Run Dose-Response", class = "btn-block",
                         style = "background-color:#6a0572; color:white;")
          ),
          box(title = "Dose-Response Curve: BDCAF at Day 90", width = 9, status = "info",
            plotlyOutput("dr_curve_plot", height = "400px"))
        ),
        fluidRow(
          box(title = "Dose-Response: Oral Ulcer Reduction", width = 4, status = "warning",
            plotlyOutput("dr_oral_plot",    height = "280px")),
          box(title = "Dose-Response: Ocular Inflammation",  width = 4, status = "warning",
            plotlyOutput("dr_ocular_plot",  height = "280px")),
          box(title = "Optimal Dose Summary Table",          width = 4, status = "success",
            DTOutput("dr_summary_tbl"))
        )
      )
    )  # end tabItems
  )
)

## ============================================================
## SERVER
## ============================================================
server <- function(input, output, session) {

  ## ---- DRUG PRESET HELPER ----
  drug_events <- function(drug_name, sim_days) {
    switch(drug_name,
      "Colchicine 0.5 mg BID" = , "Colchicine" =
        ev(data.frame(time = seq(0, (sim_days-1)*24, by = 12), cmt = 1, amt = 0.5, evid = 1)),
      "Prednisolone 40 mg/day" = , "Prednisolone" =
        ev(data.frame(time = seq(0, (sim_days-1)*24, by = 24), cmt = 4, amt = 40, evid = 1, rate = -2)),
      "Adalimumab (anti-TNF)" = , "Adalimumab 40 mg Q2W" =
        ev(data.frame(time = seq(0, sim_days*24, by = 336), cmt = 6, amt = 40, evid = 1, rate = -2)),
      "Apremilast (PDE4i)" = , "Apremilast 30 mg BID" =
        ev(data.frame(time = seq(0, (sim_days-1)*24, by = 12), cmt = 8, amt = 30, evid = 1)),
      "Canakinumab (anti-IL-1B)" = , "Canakinumab 150 mg Q8W" =
        ev(data.frame(time = seq(0, sim_days*24, by = 1344), cmt = 9, amt = 150, evid = 1, rate = -2)),
      "Colchicine + Prednisolone" = {
        e1 <- data.frame(time = seq(0, (sim_days-1)*24, by = 12), cmt = 1, amt = 0.5, evid = 1, rate = 0)
        e2 <- data.frame(time = seq(0, (sim_days-1)*24, by = 24), cmt = 4, amt = 40,  evid = 1, rate = -2)
        ev(dplyr::bind_rows(e1, e2) %>% dplyr::arrange(time))
      },
      "Adalimumab + Apremilast (Refractory)" = {
        e1 <- data.frame(time = seq(0, sim_days*24, by = 336), cmt = 6, amt = 40, evid = 1, rate = -2)
        e2 <- data.frame(time = seq(0, (sim_days-1)*24, by = 12), cmt = 8, amt = 30, evid = 1, rate = 0)
        ev(dplyr::bind_rows(e1, e2) %>% dplyr::arrange(time))
      },
      NULL
    )
  }

  scen_colors <- c(
    "Untreated"                        = "#e63946",
    "Colchicine"                       = "#f4a261",
    "Prednisolone"                     = "#2a9d8f",
    "Adalimumab (anti-TNF)"            = "#457b9d",
    "Apremilast (PDE4i)"               = "#52b788",
    "Canakinumab (anti-IL-1B)"         = "#9d4edd",
    "Colchicine + Prednisolone"        = "#f48c06",
    "Adalimumab + Apremilast (Refractory)" = "#1d3557"
  )

  ## ---- TAB 2: PK ----
  pk_sim <- eventReactive(input$run_pk, {
    drug <- input$pk_drug
    d_name <- switch(drug,
      "Colchicine"           = "colchicine",
      "Prednisolone"         = "prednisolone",
      "Adalimumab (anti-TNF)"        = "adalimumab",
      "Apremilast (PDE4i)"           = "apremilast",
      "Canakinumab (anti-IL-1B)"     = "canakinumab"
    )
    ev_drug <- build_events(d_name, input$pk_dose, input$pk_freq, input$pk_days)
    init_v  <- make_init(input$pk_sev)
    run_bd_sim(bd_mod, init_v, ev_drug, input$pk_hlab51, input$pk_days)
  })

  output$pk_plot <- renderPlotly({
    req(pk_sim())
    df <- pk_sim()
    cp_col <- switch(input$pk_drug,
      "Colchicine"               = "Cp_COL",
      "Prednisolone"             = "Cp_PRED",
      "Adalimumab (anti-TNF)"    = "Cp_ADA",
      "Apremilast (PDE4i)"       = "Cp_APR",
      "Canakinumab (anti-IL-1B)" = "Cp_CAN"
    )
    df$cp_val <- df[[cp_col]]
    p <- ggplot(df, aes(x = time_days, y = cp_val)) +
      geom_line(color = "#457b9d", linewidth = 1.2) +
      labs(title = paste0(input$pk_drug, " PK Profile"),
           x = "Time (days)", y = "Plasma Concentration (ng/mL)") +
      theme_minimal(base_size = 12)
    ggplotly(p)
  })

  output$pk_summary_tbl <- renderDT({
    req(pk_sim())
    df <- pk_sim()
    cp_col <- switch(input$pk_drug,
      "Colchicine"               = "Cp_COL",
      "Prednisolone"             = "Cp_PRED",
      "Adalimumab (anti-TNF)"    = "Cp_ADA",
      "Apremilast (PDE4i)"       = "Cp_APR",
      "Canakinumab (anti-IL-1B)" = "Cp_CAN"
    )
    df$cp_val <- df[[cp_col]]
    tbl <- data.frame(
      Metric = c("Cmax (ng/mL)", "Cmin (trough, ng/mL)", "Cavg (ng/mL)",
                 "AUC (ng·h/mL)", "t1/2 (approx, h)"),
      Value  = round(c(max(df$cp_val), min(df$cp_val[df$time_days > 1]),
                       mean(df$cp_val), sum(df$cp_val) * 4,
                       log(2) / (tail(df$cp_val[df$cp_val > 0], 1) / max(df$cp_val)) * 4), 2)
    )
    datatable(tbl, options = list(dom = "t", pageLength = 10), rownames = FALSE)
  })

  ## ---- TAB 3: Cytokine Dynamics ----
  cyt_sim <- eventReactive(input$run_cyt, {
    ev_d <- if (input$cyt_drug == "None (untreated)") NULL
            else drug_events(input$cyt_drug, input$cyt_days)
    init_v <- make_init(input$cyt_sev)
    run_bd_sim(bd_mod, init_v, ev_d, input$cyt_hlab51, input$cyt_days)
  })

  output$cyt_cytokine_plot <- renderPlotly({
    req(cyt_sim())
    df <- cyt_sim() %>% select(time_days, TNFA, IL1B, IL6C, IL17A) %>%
      pivot_longer(-time_days, names_to = "Cytokine", values_to = "Level")
    p <- ggplot(df, aes(x = time_days, y = Level, color = Cytokine)) +
      geom_line(linewidth = 1) +
      scale_color_manual(values = c(TNFA="#e63946", IL1B="#f48c06", IL6C="#457b9d", IL17A="#9d4edd")) +
      geom_hline(yintercept = 1, linetype = "dashed", color = "gray50", alpha = 0.7) +
      labs(title = "Cytokine Dynamics", x = "Time (days)", y = "Normalized Level") +
      theme_minimal(base_size = 12)
    ggplotly(p)
  })

  output$cyt_cells_plot <- renderPlotly({
    req(cyt_sim())
    df <- cyt_sim() %>% select(time_days, NEU, TH1, TH17, TREG) %>%
      pivot_longer(-time_days, names_to = "Cell", values_to = "Activity")
    p <- ggplot(df, aes(x = time_days, y = Activity, color = Cell)) +
      geom_line(linewidth = 1) +
      scale_color_manual(values = c(NEU="#52b788", TH1="#c77dff", TH17="#e63946", TREG="#2a9d8f")) +
      geom_hline(yintercept = 1, linetype = "dashed", color = "gray50", alpha = 0.7) +
      labs(title = "Immune Cell Dynamics", x = "Time (days)", y = "Normalized Activity") +
      theme_minimal(base_size = 12)
    ggplotly(p)
  })

  output$cyt_th17treg_plot <- renderPlotly({
    req(cyt_sim())
    df <- cyt_sim()
    p <- ggplot(df, aes(x = time_days)) +
      geom_line(aes(y = TH17, color = "Th17"), linewidth = 1.2) +
      geom_line(aes(y = TREG, color = "Treg"),  linewidth = 1.2) +
      geom_ribbon(aes(ymin = TREG, ymax = TH17, fill = "Th17 > Treg imbalance"), alpha = 0.2) +
      scale_color_manual(values = c(Th17 = "#e63946", Treg = "#2a9d8f")) +
      scale_fill_manual(values = c("Th17 > Treg imbalance" = "#e63946")) +
      labs(title = "Th17/Treg Balance", x = "Time (days)", y = "Normalized Level", color = "", fill = "") +
      theme_minimal(base_size = 12)
    ggplotly(p)
  })

  ## ---- TAB 4: Clinical Endpoints ----
  ep_sim <- eventReactive(input$run_ep, {
    ev_d <- if (grepl("untreated", input$ep_drug, ignore.case = TRUE)) NULL
            else drug_events(input$ep_drug, input$ep_days)
    init_v <- make_init(input$ep_sev)
    run_bd_sim(bd_mod, init_v, ev_d, input$ep_hlab51, input$ep_days)
  })

  make_ep_plot <- function(df, yvar, ytitle, yintercept = NULL) {
    req(!is.null(df))
    p <- ggplot(df, aes_string(x = "time_days", y = yvar)) +
      geom_line(color = "#c1121f", linewidth = 1.2) +
      labs(x = "Time (days)", y = ytitle) +
      theme_minimal(base_size = 12)
    if (!is.null(yintercept))
      p <- p + geom_hline(yintercept = yintercept, linetype = "dashed", color = "gray50")
    ggplotly(p)
  }

  output$ep_bdcaf_plot  <- renderPlotly({ req(ep_sim()); make_ep_plot(ep_sim(), "BDCAF",  "BDCAF Score", 3) })
  output$ep_oral_plot   <- renderPlotly({ req(ep_sim()); make_ep_plot(ep_sim(), "OUL",    "Oral Ulcer Index") })
  output$ep_ocular_plot <- renderPlotly({ req(ep_sim()); make_ep_plot(ep_sim(), "OCI",    "Ocular Score") })
  output$ep_ea_plot     <- renderPlotly({ req(ep_sim()); make_ep_plot(ep_sim(), "EA",     "Endothelial Activation") })

  output$ep_summary_tbl <- renderDT({
    req(ep_sim())
    df <- ep_sim() %>% filter(time_days >= max(time_days) - 2) %>%
      summarise(BDCAF = round(mean(BDCAF), 2), OralUlcer = round(mean(OUL), 2),
                OcularScore = round(mean(OCI), 2), EndothelialAct = round(mean(EA), 2),
                TNFA = round(mean(TNFA), 2), IL1B = round(mean(IL1B), 2),
                IL17A = round(mean(IL17A), 2),
                Remission = ifelse(mean(BDCAF) < 3, "YES", "No"))
    datatable(df, options = list(dom = "t"), rownames = FALSE)
  })

  ## ---- TAB 5: Scenario Comparison ----
  cmp_sim <- eventReactive(input$run_cmp, {
    scens <- input$cmp_scenarios
    req(length(scens) >= 1)
    init_v <- make_init(input$cmp_sev)
    purrr::map_dfr(scens, function(s) {
      ev_d <- if (s == "Untreated") NULL else drug_events(s, input$cmp_days)
      run_bd_sim(bd_mod, init_v, ev_d, input$cmp_hlab51, input$cmp_days) %>%
        mutate(scenario = s)
    })
  })

  make_cmp_plot <- function(df, yvar, ytitle, yintercept = NULL) {
    p <- ggplot(df, aes_string(x = "time_days", y = yvar, color = "scenario")) +
      geom_line(linewidth = 1) +
      scale_color_manual(values = scen_colors, name = "Treatment") +
      labs(x = "Time (days)", y = ytitle) +
      theme_minimal(base_size = 11)
    if (!is.null(yintercept))
      p <- p + geom_hline(yintercept = yintercept, linetype = "dashed", color = "gray50")
    ggplotly(p)
  }

  output$cmp_bdcaf_plot  <- renderPlotly({ req(cmp_sim()); make_cmp_plot(cmp_sim(), "BDCAF", "BDCAF Score", 3) })
  output$cmp_oral_plot   <- renderPlotly({ req(cmp_sim()); make_cmp_plot(cmp_sim(), "OUL",   "Oral Ulcer Index") })
  output$cmp_ocular_plot <- renderPlotly({ req(cmp_sim()); make_cmp_plot(cmp_sim(), "OCI",   "Ocular Score") })

  output$cmp_summary_tbl <- renderDT({
    req(cmp_sim())
    cmp_sim() %>% filter(time_days >= max(time_days) - 2) %>%
      group_by(scenario) %>%
      summarise(BDCAF_d90     = round(mean(BDCAF), 2),
                OralUlcer_d90 = round(mean(OUL), 2),
                Ocular_d90    = round(mean(OCI), 2),
                TNFA_d90      = round(mean(TNFA), 2),
                IL17A_d90     = round(mean(IL17A), 2),
                Remission     = ifelse(mean(BDCAF) < 3, "YES", "No"),
                .groups = "drop") %>%
      arrange(BDCAF_d90) %>%
      datatable(options = list(pageLength = 10), rownames = FALSE) %>%
      formatStyle("Remission", backgroundColor = styleEqual("YES", "#d4edda"))
  })

  ## ---- TAB 6: Biomarker & HLA-B51 ----
  bio_sim <- eventReactive(input$run_bio, {
    ev_d <- drug_events(input$bio_drug, input$bio_days)
    init_v <- make_init(input$bio_sev)
    pos <- run_bd_sim(bd_mod, init_v, ev_d, TRUE,  input$bio_days) %>% mutate(hlab51 = "HLA-B51 Positive")
    neg <- run_bd_sim(bd_mod, init_v, ev_d, FALSE, input$bio_days) %>% mutate(hlab51 = "HLA-B51 Negative")
    bind_rows(pos, neg)
  })

  output$bio_hlab51_plot <- renderPlotly({
    req(bio_sim())
    p <- ggplot(bio_sim(), aes(x = time_days, y = BDCAF, color = hlab51)) +
      geom_line(linewidth = 1.3) +
      scale_color_manual(values = c("HLA-B51 Positive" = "#e63946", "HLA-B51 Negative" = "#457b9d")) +
      geom_hline(yintercept = 3, linetype = "dashed", color = "gray50") +
      labs(title = paste("BDCAF Response:", input$bio_drug),
           x = "Time (days)", y = "BDCAF Score", color = "HLA-B51") +
      theme_minimal(base_size = 12)
    ggplotly(p)
  })

  output$bio_cytokine_bio_plot <- renderPlotly({
    req(bio_sim())
    df <- bio_sim() %>% filter(hlab51 == "HLA-B51 Positive") %>%
      select(time_days, TNFA, IL1B, IL6C, IL17A) %>%
      pivot_longer(-time_days, names_to = "Cytokine", values_to = "Level")
    p <- ggplot(df, aes(x = time_days, y = Level, color = Cytokine)) +
      geom_line(linewidth = 1) +
      scale_color_manual(values = c(TNFA="#e63946", IL1B="#f48c06", IL6C="#457b9d", IL17A="#9d4edd")) +
      geom_hline(yintercept = 1, linetype = "dashed", color = "gray50") +
      labs(title = "Cytokine Biomarkers (HLA-B51+)", x = "Time (days)", y = "Normalized Level") +
      theme_minimal(base_size = 11)
    ggplotly(p)
  })

  output$bio_remission_plot <- renderPlotly({
    req(bio_sim())
    df <- bio_sim() %>%
      group_by(hlab51) %>%
      filter(BDCAF < 3) %>%
      slice(1) %>%
      ungroup() %>%
      select(hlab51, time_days) %>%
      rename(Days_to_Remission = time_days)
    if (nrow(df) == 0) {
      p <- ggplot() + labs(title = "No remission achieved in simulation period") + theme_minimal()
    } else {
      p <- ggplot(df, aes(x = hlab51, y = Days_to_Remission, fill = hlab51)) +
        geom_bar(stat = "identity", width = 0.5) +
        scale_fill_manual(values = c("HLA-B51 Positive" = "#e63946", "HLA-B51 Negative" = "#457b9d")) +
        labs(title = "Days to BDCAF < 3 (Remission)", x = "", y = "Days", fill = "") +
        theme_minimal(base_size = 12) +
        theme(legend.position = "none")
    }
    ggplotly(p)
  })

  ## ---- TAB 7: Dose Optimization ----
  dr_sim <- eventReactive(input$run_dr, {
    doses <- seq(input$dr_dose_min, input$dr_dose_max, length.out = input$dr_n_doses)
    d_name <- switch(input$dr_drug,
      "Adalimumab (anti-TNF)"    = "adalimumab",
      "Canakinumab (anti-IL-1B)" = "canakinumab",
      "Apremilast (PDE4i)"       = "apremilast",
      "Prednisolone"             = "prednisolone"
    )
    init_v <- make_init(input$dr_sev)
    purrr::map_dfr(doses, function(d) {
      freq <- switch(d_name,
        adalimumab  = 14, canakinumab = 56, apremilast = 0.5, prednisolone = 1, 1
      )
      ev_d <- build_events(d_name, d, freq, input$dr_days)
      sim <- run_bd_sim(bd_mod, init_v, ev_d, input$dr_hlab51, input$dr_days)
      sim %>% filter(time_days >= max(time_days) - 2) %>%
        summarise(dose_mg = d,
                  BDCAF_d90 = mean(BDCAF), OUL_d90 = mean(OUL), OCI_d90 = mean(OCI),
                  TNFA_d90 = mean(TNFA), IL17A_d90 = mean(IL17A),
                  Remission = mean(BDCAF) < 3)
    })
  })

  output$dr_curve_plot <- renderPlotly({
    req(dr_sim())
    p <- ggplot(dr_sim(), aes(x = dose_mg, y = BDCAF_d90)) +
      geom_line(color = "#6a0572", linewidth = 1.3) +
      geom_point(aes(color = Remission), size = 3) +
      scale_color_manual(values = c("TRUE" = "#2dc653", "FALSE" = "#e63946"),
                         labels = c("No", "Yes"), name = "Remission") +
      geom_hline(yintercept = 3, linetype = "dashed", color = "gray50") +
      labs(title = paste0(input$dr_drug, " Dose-Response (Day 90 BDCAF)"),
           x = "Dose (mg)", y = "BDCAF Score at Day 90") +
      theme_minimal(base_size = 12)
    ggplotly(p)
  })

  output$dr_oral_plot <- renderPlotly({
    req(dr_sim())
    p <- ggplot(dr_sim(), aes(x = dose_mg, y = OUL_d90)) +
      geom_line(color = "#f48c06", linewidth = 1.3) + geom_point(color = "#f48c06") +
      labs(title = "Oral Ulcer Score", x = "Dose (mg)", y = "Oral Ulcer Index") +
      theme_minimal(base_size = 11)
    ggplotly(p)
  })

  output$dr_ocular_plot <- renderPlotly({
    req(dr_sim())
    p <- ggplot(dr_sim(), aes(x = dose_mg, y = OCI_d90)) +
      geom_line(color = "#457b9d", linewidth = 1.3) + geom_point(color = "#457b9d") +
      labs(title = "Ocular Inflammation", x = "Dose (mg)", y = "Ocular Score") +
      theme_minimal(base_size = 11)
    ggplotly(p)
  })

  output$dr_summary_tbl <- renderDT({
    req(dr_sim())
    dr_sim() %>%
      mutate(across(where(is.numeric), ~ round(.x, 3)),
             Remission = ifelse(Remission, "YES", "No")) %>%
      arrange(dose_mg) %>%
      datatable(options = list(pageLength = 10, dom = "t"), rownames = FALSE) %>%
      formatStyle("Remission", backgroundColor = styleEqual("YES", "#d4edda"))
  })
}

## ============================================================
## Launch
## ============================================================
shinyApp(ui, server)
