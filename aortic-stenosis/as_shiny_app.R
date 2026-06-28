## ============================================================
## Calcific Aortic Valve Stenosis (CAVD/AS) — QSP Shiny App
## Interactive QSP Dashboard
## 6 Tabs: Patient Profile · PK · PD Key Metrics · Clinical Endpoints
##         · Scenario Comparison · Biomarker Panel
## ============================================================

library(shiny)
library(shinydashboard)
library(mrgsolve)
library(tidyverse)
library(ggplot2)
library(plotly)
library(DT)
library(patchwork)
library(scales)
library(viridis)

## ============================================================
## INLINE MODEL CODE (self-contained for Shiny)
## ============================================================
code_as_shiny <- '
$PARAM
AGE=70 BMI=27 BSA=1.80 HR=70 SV0=75 SVR0=1200 MAP0=93 LVEDP0=8
k_calc=0.003 CS_max=3000 hill_cs=2.0 AVA0=2.5 AVA_min=0.3 CS_50=500
LVMI0=85 LV_hyp_max=200 LV_hyp_k=0.15 k_col_form=0.06 k_col_deg=0.02
collagen_max=0.35 LVEF0=65 LVEF_min=20 k_LVEF_loss=0.03
AngII0=1.0 BNP0=50 k_BNP=0.25 LVEDP_0=8.0
LDL0=3.5 LPA0=60 RANKL0=1.0 IL6_0=2.5 MGP0=0.30
STATIN_F=0.20 STATIN_ka=1.5 STATIN_CL=50 STATIN_V1=134 STATIN_V2=100 STATIN_Q=20
EC50_statin_LDL=20 Emax_statin_LDL=0.55 EC50_statin_IL6=40 Emax_statin_IL6=0.25
PCSK9i_F=0.72 PCSK9i_ka=0.012 PCSK9i_CL=0.30 PCSK9i_V=3.5
EC50_pcsk9i_LDL=5 Emax_pcsk9i_LDL=0.70 EC50_pcsk9i_LPA=15 Emax_pcsk9i_LPA=0.30
DENO_F=0.62 DENO_ka=0.004 DENO_CL=0.18 DENO_V=3.0 DENO_ke_target=0.02
EC50_deno_RANKL=2.0 Emax_deno_RANKL=0.80 Emax_deno_calc=0.30
VK2_F=0.85 VK2_ka=0.5 VK2_CL=2.5 VK2_V=35
EC50_vk2_MGP=0.3 Emax_vk2_MGP=0.65 MGP_effect_calc=0.20
ACEi_F=0.56 ACEi_ka=0.8 ACEi_CL=3.2 ACEi_V=57
EC50_acei_AngII=0.2 Emax_acei_AngII=0.70 Emax_acei_fib=0.40

$CMT STATIN_GUT STATIN_CENTRAL STATIN_PERIPH
     PCSK9I_DEPOT PCSK9I_CENTRAL
     DENO_DEPOT DENO_CENTRAL
     VK2_DEPOT VK2_CENTRAL
     ACEI_DEPOT ACEI_CENTRAL
     CS LDL_C LPA RANKL IL6 MGP_carbox AngII LVMI COLLAGEN LVEF NTproBNP

$MAIN
double Cp_statin = STATIN_CENTRAL / STATIN_V1;
double Cp_pcsk9i = PCSK9I_CENTRAL / PCSK9i_V;
double Cp_deno   = DENO_CENTRAL / DENO_V;
double Cp_vk2    = VK2_CENTRAL / VK2_V;
double Cp_acei   = ACEI_CENTRAL / ACEi_V;

double E_statin_LDL = (Emax_statin_LDL * Cp_statin)/(EC50_statin_LDL + Cp_statin);
double E_statin_IL6 = (Emax_statin_IL6 * Cp_statin)/(EC50_statin_IL6 + Cp_statin);
double E_pcsk9i_LDL = (Emax_pcsk9i_LDL * Cp_pcsk9i)/(EC50_pcsk9i_LDL + Cp_pcsk9i);
double E_pcsk9i_LPA = (Emax_pcsk9i_LPA * Cp_pcsk9i)/(EC50_pcsk9i_LPA + Cp_pcsk9i);
double E_total_LDL  = fmin(E_statin_LDL + E_pcsk9i_LDL, 0.85);
double E_deno_RANKL = (Emax_deno_RANKL * Cp_deno)/(EC50_deno_RANKL + Cp_deno);
double E_deno_calc  = E_deno_RANKL * Emax_deno_calc;
double E_vk2_MGP    = (Emax_vk2_MGP * Cp_vk2)/(EC50_vk2_MGP + Cp_vk2);
double E_vk2_calc   = E_vk2_MGP * MGP_effect_calc;
double E_acei_AngII = (Emax_acei_AngII * Cp_acei)/(EC50_acei_AngII + Cp_acei);
double E_acei_fib   = E_acei_AngII * Emax_acei_fib;

double AVA = AVA_min + (AVA0 - AVA_min) / (1.0 + pow(CS / CS_50, hill_cs));
double CO  = HR * (SV0 * (LVEF / LVEF0)) / 1000.0;
double Vmax = CO / (fmax(AVA, 0.3) * 60.0 * 0.785);
double MeanPG = 2.4 * Vmax * Vmax;
double afterload = fmax(MeanPG, 0.0) + MAP0;
double wall_stress = afterload / (LVEF / 100.0 + 0.1);

$ODE
dxdt_STATIN_GUT     = -STATIN_ka * STATIN_GUT;
dxdt_STATIN_CENTRAL = STATIN_F * STATIN_ka * STATIN_GUT
                      - (STATIN_CL/STATIN_V1)*STATIN_CENTRAL
                      - (STATIN_Q/STATIN_V1)*STATIN_CENTRAL
                      + (STATIN_Q/STATIN_V2)*STATIN_PERIPH;
dxdt_STATIN_PERIPH  = (STATIN_Q/STATIN_V1)*STATIN_CENTRAL
                      - (STATIN_Q/STATIN_V2)*STATIN_PERIPH;
dxdt_PCSK9I_DEPOT   = -PCSK9i_ka * PCSK9I_DEPOT;
dxdt_PCSK9I_CENTRAL = PCSK9i_F * PCSK9i_ka * PCSK9I_DEPOT
                      - (PCSK9i_CL/PCSK9i_V)*PCSK9I_CENTRAL;
dxdt_DENO_DEPOT   = -DENO_ka * DENO_DEPOT;
dxdt_DENO_CENTRAL = DENO_F * DENO_ka * DENO_DEPOT
                    - (DENO_CL/DENO_V)*DENO_CENTRAL
                    - DENO_ke_target * RANKL * DENO_CENTRAL;
dxdt_VK2_DEPOT   = -VK2_ka * VK2_DEPOT;
dxdt_VK2_CENTRAL = VK2_F * VK2_ka * VK2_DEPOT - (VK2_CL/VK2_V)*VK2_CENTRAL;
dxdt_ACEI_DEPOT   = -ACEi_ka * ACEI_DEPOT;
dxdt_ACEI_CENTRAL = ACEi_F * ACEi_ka * ACEI_DEPOT - (ACEi_CL/ACEi_V)*ACEI_CENTRAL;

double k_LDL = 0.1;
dxdt_LDL_C = k_LDL * (LDL0 * (1.0 - E_total_LDL) - LDL_C);

dxdt_LPA = 0.08 * (LPA0 * (1.0 - E_pcsk9i_LPA) - LPA);

double RANKL_prod = RANKL0 * (1.0 + 0.05 * CS / CS_50);
dxdt_RANKL = RANKL_prod - 0.5*RANKL - E_deno_RANKL*RANKL;

dxdt_IL6 = IL6_0*(1.0+0.3*(CS/CS_50))*(1.0-E_statin_IL6) - 0.8*IL6;

double MGP_target = MGP0 + E_vk2_MGP * (1.0 - MGP0);
dxdt_MGP_carbox = 0.3 * (MGP_target - MGP_carbox);

double calc_driver  = (LDL_C/LDL0)*0.3 + (RANKL/RANKL0)*0.4 + (LPA/LPA0)*0.3;
double calc_inhibit = E_deno_calc*0.4 + E_vk2_calc*0.3 + MGP_carbox*0.15 + E_statin_LDL*0.05;
double k_cs = k_calc * calc_driver * (1.0 - fmin(calc_inhibit, 0.7));
dxdt_CS = k_cs * fmax(CS, 100.0) * (1.0 - CS/CS_max);

dxdt_AngII = AngII0*(1.0+0.2*(afterload/130.0)) - (0.5+E_acei_AngII)*AngII;

double LVMI_target = fmin(LVMI0*(1.0+0.8*fmax(afterload-120.0,0.0)/100.0), LV_hyp_max);
dxdt_LVMI = LV_hyp_k * (LVMI_target - LVMI);

double col_form = k_col_form*(AngII/AngII0)*(wall_stress/100.0)*(1.0-E_acei_fib);
dxdt_COLLAGEN = col_form - k_col_deg*COLLAGEN;
if (COLLAGEN > collagen_max) dxdt_COLLAGEN = -k_col_deg*COLLAGEN;
if (COLLAGEN < 0) dxdt_COLLAGEN = 0;

dxdt_LVEF = -k_LVEF_loss*(COLLAGEN/collagen_max)*(fmax(afterload-130.0,0.0)/50.0)*LVEF;
if (LVEF < LVEF_min) dxdt_LVEF = 0;
if (LVEF > 80) dxdt_LVEF = 0;

double lvedp_proxy = LVEDP0*(1.0+5.0*COLLAGEN)*(LVEF0/fmax(LVEF,25.0));
dxdt_NTproBNP = k_BNP*(lvedp_proxy/LVEDP0)*BNP0 - 0.5*NTproBNP;

$TABLE
double AVA_o  = AVA_min + (AVA0-AVA_min)/(1.0+pow(CS/CS_50, hill_cs));
double CO_o   = HR*(SV0*(LVEF/LVEF0))/1000.0;
double PG_o   = 2.4 * pow(CO_o/(fmax(AVA_o,0.3)*60.0*0.785), 2);
double NYHA = 1.0;
if (NTproBNP > 125 && LVEF < 60) NYHA = 2.0;
if (NTproBNP > 600 && LVEF < 50) NYHA = 3.0;
if (NTproBNP > 1800 && LVEF < 40) NYHA = 4.0;
double AS_grade = (AVA_o >= 1.5) ? 1.0 : ((AVA_o >= 1.0) ? 2.0 : 3.0);
double lvedp_o = LVEDP0*(1.0+5.0*COLLAGEN)*(LVEF0/fmax(LVEF,25.0));
capture Cp_statin_ngmL = STATIN_CENTRAL/STATIN_V1;
capture Cp_pcsk9i_ugmL = PCSK9I_CENTRAL/PCSK9i_V;
capture Cp_deno_ugmL   = DENO_CENTRAL/DENO_V;
capture AVA_cm2     = AVA_o;
capture MeanPG_mmHg = PG_o;
capture CO_Lmin     = CO_o;
capture LVEF_pct    = LVEF;
capture LVMI_gm2    = LVMI;
capture Collagen_f  = COLLAGEN;
capture NTproBNP_pg = NTproBNP;
capture LVEDP_mmHg  = lvedp_o;
capture LDL_mmolL   = LDL_C;
capture LPA_mgdL    = LPA;
capture IL6_pgmL    = IL6;
capture RANKL_norm  = RANKL;
capture AngII_norm  = AngII;
capture MGP_f       = MGP_carbox;
capture CS_au       = CS;
capture NYHA_class  = NYHA;
capture AS_grade    = AS_grade;

$INIT
STATIN_GUT=0 STATIN_CENTRAL=0 STATIN_PERIPH=0
PCSK9I_DEPOT=0 PCSK9I_CENTRAL=0
DENO_DEPOT=0 DENO_CENTRAL=0
VK2_DEPOT=0 VK2_CENTRAL=0
ACEI_DEPOT=0 ACEI_CENTRAL=0
CS=100 LDL_C=3.5 LPA=60 RANKL=1.0 IL6=2.5 MGP_carbox=0.30
AngII=1.0 LVMI=85 COLLAGEN=0.05 LVEF=65 NTproBNP=50
'

## ============================================================
## COMPILE MODEL (once, outside UI/Server)
## ============================================================
mod_shiny <- mread_cache("as_shiny", tempdir(), code_as_shiny,
                         quiet = TRUE)

## ============================================================
## HELPER: Build event table
## ============================================================
build_events <- function(use_statin, statin_dose_mg,
                          use_pcsk9i, pcsk9i_dose_mg, pcsk9i_freq_days,
                          use_deno, deno_dose_mg, deno_freq_days,
                          use_vk2, vk2_dose_ug,
                          use_acei, acei_dose_mg,
                          sim_years = 10) {
  sim_hours <- sim_years * 8760
  ev_list <- list()

  if (use_statin && statin_dose_mg > 0) {
    ev_list[["statin"]] <- ev(
      time = seq(0, sim_hours - 1, by = 24),
      cmt  = "STATIN_GUT",
      amt  = statin_dose_mg * 1e6  # mg → ng
    )
  }
  if (use_pcsk9i && pcsk9i_dose_mg > 0) {
    ev_list[["pcsk9i"]] <- ev(
      time = seq(0, sim_hours - 1, by = 24 * pcsk9i_freq_days),
      cmt  = "PCSK9I_DEPOT",
      amt  = pcsk9i_dose_mg * 1000  # mg → μg
    )
  }
  if (use_deno && deno_dose_mg > 0) {
    ev_list[["deno"]] <- ev(
      time = seq(0, sim_hours - 1, by = 24 * deno_freq_days),
      cmt  = "DENO_DEPOT",
      amt  = deno_dose_mg * 1000  # mg → μg
    )
  }
  if (use_vk2 && vk2_dose_ug > 0) {
    ev_list[["vk2"]] <- ev(
      time = seq(0, sim_hours - 1, by = 24),
      cmt  = "VK2_DEPOT",
      amt  = vk2_dose_ug  # μg
    )
  }
  if (use_acei && acei_dose_mg > 0) {
    ev_list[["acei"]] <- ev(
      time = seq(0, sim_hours - 1, by = 24),
      cmt  = "ACEI_DEPOT",
      amt  = acei_dose_mg * 1e6  # mg → ng
    )
  }

  if (length(ev_list) == 0) {
    return(ev(time = 0, cmt = "STATIN_GUT", amt = 0))
  }
  Reduce(`+`, ev_list)
}

## ============================================================
## UI
## ============================================================
ui <- dashboardPage(
  skin = "blue",
  dashboardHeader(
    title = span(icon("heartbeat"), "AS / CAVD QSP Dashboard"),
    titleWidth = 320
  ),

  dashboardSidebar(
    width = 300,
    sidebarMenu(
      id = "tabs",
      menuItem("Patient Profile",     tabName = "tab_patient",  icon = icon("user-md")),
      menuItem("Drug PK",             tabName = "tab_pk",       icon = icon("flask")),
      menuItem("PD Key Metrics",      tabName = "tab_pd",       icon = icon("chart-line")),
      menuItem("Clinical Endpoints",  tabName = "tab_clinical", icon = icon("stethoscope")),
      menuItem("Scenario Comparison", tabName = "tab_scenarios",icon = icon("balance-scale")),
      menuItem("Biomarker Panel",     tabName = "tab_biomarker",icon = icon("dna"))
    ),
    hr(),
    h4("Patient Parameters", style = "padding-left:15px; color:#ccc;"),
    sliderInput("p_age",  "Age (years):", 40, 90, 70, step = 1),
    sliderInput("p_ldl",  "LDL-C (mmol/L):", 1.5, 7.0, 3.5, step = 0.1),
    sliderInput("p_lpa",  "Lp(a) (mg/dL):", 5, 250, 60, step = 5),
    sliderInput("p_lvef", "Baseline LVEF (%):", 30, 75, 65, step = 1),
    sliderInput("p_ava",  "Baseline AVA (cm²):", 1.0, 2.8, 2.5, step = 0.1),
    sliderInput("p_cs",   "Baseline Ca Score (AU):", 50, 500, 100, step = 25),
    sliderInput("p_years","Simulation (years):", 3, 15, 10, step = 1),
    hr(),
    h4("Drug Treatment", style = "padding-left:15px; color:#ccc;"),
    checkboxInput("use_statin",  "Statin (Rosuvastatin)", TRUE),
    conditionalPanel("input.use_statin",
      sliderInput("statin_dose", "Dose (mg/d):", 5, 40, 20, step = 5)
    ),
    checkboxInput("use_pcsk9i", "PCSK9i (Evolocumab)", FALSE),
    conditionalPanel("input.use_pcsk9i",
      sliderInput("pcsk9i_dose", "Dose (mg):", 70, 420, 140, step = 70),
      sliderInput("pcsk9i_freq", "Every (days):", 7, 28, 14, step = 7)
    ),
    checkboxInput("use_deno",   "Denosumab (anti-RANKL)*", FALSE),
    conditionalPanel("input.use_deno",
      sliderInput("deno_dose",  "Dose (mg):", 30, 120, 60, step = 30),
      sliderInput("deno_freq",  "Every (days):", 90, 365, 182, step = 30)
    ),
    checkboxInput("use_vk2",    "Vitamin K2 (MK-7)", FALSE),
    conditionalPanel("input.use_vk2",
      sliderInput("vk2_dose",   "Dose (μg/d):", 90, 360, 180, step = 90)
    ),
    checkboxInput("use_acei",   "ACEi/ARB (Ramipril)", FALSE),
    conditionalPanel("input.use_acei",
      sliderInput("acei_dose",  "Dose (mg/d):", 1.25, 10, 5, step = 1.25)
    ),
    hr(),
    em("* Denosumab for CAVD is investigational / hypothesis", style = "padding-left:10px; font-size:10px; color:#aaa;"),
    actionButton("run_sim", "Run Simulation", class = "btn-success btn-block",
                 style = "margin:10px;"),
    actionButton("reset_params", "Reset Defaults", class = "btn-default btn-block",
                 style = "margin:10px;")
  ),

  dashboardBody(
    tags$head(tags$style(HTML("
      .box-header {font-weight: bold;}
      .severity-box {text-align:center; padding:15px; border-radius:8px; margin:5px;}
      .severity-normal {background:#E8F5E9; border:2px solid #4CAF50;}
      .severity-mild   {background:#FFF8E1; border:2px solid #FF9800;}
      .severity-moderate{background:#FFF3E0; border:2px solid #FF5722;}
      .severity-severe  {background:#FFEBEE; border:2px solid #F44336;}
    "))),

    tabItems(
      ## ----------------------------------------------------------
      ## TAB 1: Patient Profile
      ## ----------------------------------------------------------
      tabItem(tabName = "tab_patient",
        fluidRow(
          valueBoxOutput("vbox_ava",   width = 3),
          valueBoxOutput("vbox_pg",    width = 3),
          valueBoxOutput("vbox_lvef",  width = 3),
          valueBoxOutput("vbox_nyha",  width = 3)
        ),
        fluidRow(
          box(title = "Current AS Severity Assessment", width = 12, status = "primary",
            fluidRow(
              column(3, uiOutput("severity_box_ava")),
              column(3, uiOutput("severity_box_pg")),
              column(3, uiOutput("severity_box_lvmi")),
              column(3, uiOutput("severity_box_bnp"))
            )
          )
        ),
        fluidRow(
          box(title = "Simulated Patient Timeline — AVA & Gradient", width = 8, status = "info",
            plotlyOutput("plot_patient_overview", height = "380px")
          ),
          box(title = "Patient Risk Summary", width = 4, status = "warning",
            h4("Key Risk Factors:"),
            tableOutput("tbl_risk_profile"),
            br(),
            h4("10-Year Prognosis (without AVR):"),
            uiOutput("prognosis_text")
          )
        )
      ),

      ## ----------------------------------------------------------
      ## TAB 2: Drug PK
      ## ----------------------------------------------------------
      tabItem(tabName = "tab_pk",
        fluidRow(
          box(title = "Plasma Concentration — Statin (Rosuvastatin)", width = 6, status = "primary",
            plotlyOutput("plot_pk_statin", height = "320px"),
            p("2-compartment oral PK. Steady state achieved ~5 days.", style = "font-size:11px; color:#666;")
          ),
          box(title = "Plasma Concentration — PCSK9i (Evolocumab SC)", width = 6, status = "info",
            plotlyOutput("plot_pk_pcsk9i", height = "320px"),
            p("1-compartment SC. Tmax ~3 days. Biweekly dosing maintains >10 μg/mL trough.",
              style = "font-size:11px; color:#666;")
          )
        ),
        fluidRow(
          box(title = "Plasma Concentration — Denosumab SC", width = 6, status = "success",
            plotlyOutput("plot_pk_deno", height = "320px"),
            p("1-compartment SC with target-mediated disposition. Tmax ~10 days. Q6M dosing.",
              style = "font-size:11px; color:#666;")
          ),
          box(title = "PK Summary Table (Steady State)", width = 6, status = "warning",
            DTOutput("tbl_pk_summary"),
            br(),
            p("Note: Times shown in first 90 days for PK visualization.",
              style = "font-size:11px; color:#666;")
          )
        )
      ),

      ## ----------------------------------------------------------
      ## TAB 3: PD Key Metrics
      ## ----------------------------------------------------------
      tabItem(tabName = "tab_pd",
        fluidRow(
          box(title = "LDL Cholesterol Response", width = 6, status = "primary",
            plotlyOutput("plot_pd_ldl", height = "300px")
          ),
          box(title = "Lipoprotein(a) Response", width = 6, status = "warning",
            plotlyOutput("plot_pd_lpa", height = "300px")
          )
        ),
        fluidRow(
          box(title = "RANKL Activity & Calcification Driver", width = 6, status = "danger",
            plotlyOutput("plot_pd_rankl", height = "300px")
          ),
          box(title = "MGP Carboxylation & IL-6", width = 6, status = "success",
            plotlyOutput("plot_pd_mgp_il6", height = "300px")
          )
        ),
        fluidRow(
          box(title = "PD Effect Summary at Year 10", width = 12, status = "info",
            DTOutput("tbl_pd_effects")
          )
        )
      ),

      ## ----------------------------------------------------------
      ## TAB 4: Clinical Endpoints
      ## ----------------------------------------------------------
      tabItem(tabName = "tab_clinical",
        fluidRow(
          box(title = "Aortic Valve Area (AVA) — Primary Endpoint", width = 6, status = "danger",
            plotlyOutput("plot_clin_ava", height = "320px")
          ),
          box(title = "Mean Transvalvular Gradient (ΔP)", width = 6, status = "warning",
            plotlyOutput("plot_clin_pg", height = "320px")
          )
        ),
        fluidRow(
          box(title = "LV Ejection Fraction & LVMI", width = 6, status = "info",
            plotlyOutput("plot_clin_lv", height = "320px")
          ),
          box(title = "NYHA Functional Class & NT-proBNP", width = 6, status = "success",
            plotlyOutput("plot_clin_bnp", height = "320px")
          )
        ),
        fluidRow(
          box(title = "Projected Time to Severe AS Onset", width = 12, status = "primary",
            fluidRow(
              column(6, uiOutput("time_to_severe")),
              column(6, uiOutput("time_to_avr"))
            )
          )
        )
      ),

      ## ----------------------------------------------------------
      ## TAB 5: Scenario Comparison
      ## ----------------------------------------------------------
      tabItem(tabName = "tab_scenarios",
        fluidRow(
          box(title = "Scenario Comparison — Preset Treatment Strategies",
              width = 12, status = "primary",
            p("Compare 6 pre-defined treatment scenarios. Uses current patient parameters."),
            fluidRow(
              column(4, plotlyOutput("plot_scen_ava", height = "280px")),
              column(4, plotlyOutput("plot_scen_cs",  height = "280px")),
              column(4, plotlyOutput("plot_scen_lvef",height = "280px"))
            ),
            fluidRow(
              column(4, plotlyOutput("plot_scen_pg",   height = "280px")),
              column(4, plotlyOutput("plot_scen_lvmi", height = "280px")),
              column(4, plotlyOutput("plot_scen_bnp",  height = "280px"))
            )
          )
        ),
        fluidRow(
          box(title = "Scenario Comparison Table — Year 10 Outcomes",
              width = 12, status = "info",
            DTOutput("tbl_scenario_comparison")
          )
        )
      ),

      ## ----------------------------------------------------------
      ## TAB 6: Biomarker Panel
      ## ----------------------------------------------------------
      tabItem(tabName = "tab_biomarker",
        fluidRow(
          box(title = "Valve Calcium Score Trajectory", width = 6, status = "danger",
            plotlyOutput("plot_bio_cs", height = "300px"),
            p("CT Agatston Units. Annual increase ~200-300 AU untreated.",
              style = "font-size:11px; color:#666;")
          ),
          box(title = "Fibrosis Index (Interstitial Collagen Fraction)", width = 6, status = "warning",
            plotlyOutput("plot_bio_collagen", height = "300px"),
            p("Surrogate for myocardial fibrosis. Correlates with T1 mapping / ECV on CMR.",
              style = "font-size:11px; color:#666;")
          )
        ),
        fluidRow(
          box(title = "Cardiac Biomarker Heatmap (Year 10)", width = 6, status = "info",
            plotlyOutput("plot_bio_heatmap", height = "320px")
          ),
          box(title = "Biomarker Reference Table", width = 6, status = "success",
            DTOutput("tbl_biomarker_ref"),
            br(),
            p("Reference ranges from published guidelines (AHA/ACC 2021 AS Guidelines).",
              style = "font-size:11px; color:#666;")
          )
        ),
        fluidRow(
          box(title = "AngII & LV Remodeling Axis", width = 12, status = "primary",
            plotlyOutput("plot_bio_angii_lv", height = "280px")
          )
        )
      )
    )  # end tabItems
  )  # end dashboardBody
)  # end dashboardPage

## ============================================================
## SERVER
## ============================================================
server <- function(input, output, session) {

  ## ------ Reactive simulation --------------------------------
  sim_data <- eventReactive(input$run_sim, {
    req(input$p_age, input$p_ldl, input$p_lpa)

    # Custom patient parameters
    mod_patient <- param(mod_shiny,
                         AGE   = input$p_age,
                         LDL0  = input$p_ldl,
                         LPA0  = input$p_lpa,
                         LVEF0 = input$p_lvef,
                         AVA0  = input$p_ava)

    init_vals <- c(CS = input$p_cs,
                   LDL_C = input$p_ldl,
                   LPA   = input$p_lpa,
                   RANKL = 1.0, IL6 = 2.5,
                   MGP_carbox = 0.30, AngII = 1.0,
                   LVMI = 85, COLLAGEN = 0.05,
                   LVEF = input$p_lvef,
                   NTproBNP = 50)

    mod_patient <- init(mod_patient, .list = init_vals)

    # Build dosing events
    ev_tbl <- build_events(
      use_statin       = input$use_statin,
      statin_dose_mg   = if (input$use_statin) input$statin_dose else 0,
      use_pcsk9i       = input$use_pcsk9i,
      pcsk9i_dose_mg   = if (input$use_pcsk9i) input$pcsk9i_dose else 0,
      pcsk9i_freq_days = if (input$use_pcsk9i) input$pcsk9i_freq else 14,
      use_deno         = input$use_deno,
      deno_dose_mg     = if (input$use_deno) input$deno_dose else 0,
      deno_freq_days   = if (input$use_deno) input$deno_freq else 182,
      use_vk2          = input$use_vk2,
      vk2_dose_ug      = if (input$use_vk2) input$vk2_dose else 0,
      use_acei         = input$use_acei,
      acei_dose_mg     = if (input$use_acei) input$acei_dose else 0,
      sim_years        = input$p_years
    )

    sim_hrs <- input$p_years * 8760
    result <- mrgsim(mod_patient, events = ev_tbl,
                     end = sim_hrs, delta = 24) %>%
      as_tibble() %>%
      mutate(year = time / 8760,
             AS_label = case_when(
               AVA_cm2 >= 1.5 ~ "Mild/Normal",
               AVA_cm2 >= 1.0 ~ "Moderate",
               TRUE           ~ "Severe"
             ))
    result
  }, ignoreNULL = FALSE)

  ## ------ Scenario comparison data ---------------------------
  scenario_data <- eventReactive(input$run_sim, {
    req(input$p_age)

    mod_patient <- param(mod_shiny,
                         AGE = input$p_age, LDL0 = input$p_ldl,
                         LPA0 = input$p_lpa, LVEF0 = input$p_lvef,
                         AVA0 = input$p_ava)
    init_vals <- c(CS = input$p_cs, LDL_C = input$p_ldl, LPA = input$p_lpa,
                   RANKL = 1.0, IL6 = 2.5, MGP_carbox = 0.30,
                   AngII = 1.0, LVMI = 85, COLLAGEN = 0.05,
                   LVEF = input$p_lvef, NTproBNP = 50)
    mod_patient <- init(mod_patient, .list = init_vals)
    sim_hrs <- input$p_years * 8760

    scenarios <- list(
      "1. No Treatment"        = build_events(FALSE,0,FALSE,0,14,FALSE,0,182,FALSE,0,FALSE,0,input$p_years),
      "2. Statin Only"         = build_events(TRUE,20,FALSE,0,14,FALSE,0,182,FALSE,0,FALSE,0,input$p_years),
      "3. Statin + PCSK9i"     = build_events(TRUE,20,TRUE,140,14,FALSE,0,182,FALSE,0,FALSE,0,input$p_years),
      "4. Statin + Denosumab*" = build_events(TRUE,20,FALSE,0,14,TRUE,60,182,FALSE,0,FALSE,0,input$p_years),
      "5. Statin + VK2"        = build_events(TRUE,20,FALSE,0,14,FALSE,0,182,TRUE,180,FALSE,0,input$p_years),
      "6. Max Medical Therapy" = build_events(TRUE,20,TRUE,140,14,TRUE,60,182,TRUE,180,TRUE,5,input$p_years)
    )

    bind_rows(lapply(names(scenarios), function(nm) {
      mrgsim(mod_patient, events = scenarios[[nm]], end = sim_hrs, delta = 24) %>%
        as_tibble() %>%
        mutate(scenario = nm, year = time / 8760)
    }))
  }, ignoreNULL = FALSE)

  scen_colors <- c(
    "1. No Treatment"        = "#E53935",
    "2. Statin Only"         = "#FB8C00",
    "3. Statin + PCSK9i"     = "#8E24AA",
    "4. Statin + Denosumab*" = "#1E88E5",
    "5. Statin + VK2"        = "#43A047",
    "6. Max Medical Therapy" = "#00ACC1"
  )

  ## ------ Helper: last row of simulation ---------------------
  last_vals <- reactive({
    tail(sim_data(), 1)
  })

  ## ------ Value boxes ----------------------------------------
  output$vbox_ava <- renderValueBox({
    lv <- last_vals()
    ava <- round(lv$AVA_cm2, 2)
    color <- if (ava >= 1.5) "green" else if (ava >= 1.0) "yellow" else "red"
    label <- if (ava >= 1.5) "Mild/Normal" else if (ava >= 1.0) "Moderate" else "Severe"
    valueBox(paste0(ava, " cm²"), paste("AVA —", label), icon = icon("heart"),
             color = color)
  })
  output$vbox_pg <- renderValueBox({
    lv <- last_vals()
    pg <- round(lv$MeanPG_mmHg, 1)
    color <- if (pg < 25) "green" else if (pg < 40) "yellow" else "red"
    valueBox(paste0(pg, " mmHg"), "Mean Gradient", icon = icon("tachometer-alt"),
             color = color)
  })
  output$vbox_lvef <- renderValueBox({
    lv <- last_vals()
    ef <- round(lv$LVEF_pct, 1)
    color <- if (ef >= 50) "green" else if (ef >= 40) "yellow" else "red"
    valueBox(paste0(ef, "%"), "LVEF", icon = icon("heartbeat"), color = color)
  })
  output$vbox_nyha <- renderValueBox({
    lv <- last_vals()
    nyha <- round(lv$NYHA_class, 0)
    color <- if (nyha <= 1) "green" else if (nyha <= 2) "yellow" else "red"
    label <- c("Class I", "Class II", "Class III", "Class IV")[min(nyha, 4)]
    valueBox(label, "NYHA Function", icon = icon("walking"), color = color)
  })

  ## ------ Patient overview plot ------------------------------
  output$plot_patient_overview <- renderPlotly({
    d <- sim_data()
    p <- ggplot(d, aes(x = year)) +
      geom_line(aes(y = AVA_cm2, color = "AVA (cm²)"), linewidth = 1.3) +
      geom_hline(yintercept = c(1.5, 1.0), linetype = "dashed",
                 color = c("#FF9800", "#F44336")) +
      scale_color_manual(values = c("AVA (cm²)" = "#1E88E5")) +
      labs(x = "Year", y = "AVA (cm²)", color = "") +
      theme_bw(base_size = 12)
    ggplotly(p)
  })

  ## ------ Risk profile table ---------------------------------
  output$tbl_risk_profile <- renderTable({
    tibble(
      "Parameter"  = c("Age", "LDL-C", "Lp(a)", "LVEF", "Baseline AVA"),
      "Value"      = c(input$p_age, input$p_ldl, input$p_lpa, input$p_lvef, input$p_ava),
      "Unit"       = c("years", "mmol/L", "mg/dL", "%", "cm²"),
      "Risk Level" = c(
        if (input$p_age > 75) "High" else if (input$p_age > 65) "Moderate" else "Low",
        if (input$p_ldl > 4.5) "High" else if (input$p_ldl > 3.0) "Moderate" else "Low",
        if (input$p_lpa > 100) "High" else if (input$p_lpa > 50) "Moderate" else "Low",
        if (input$p_lvef < 50) "Reduced" else if (input$p_lvef < 60) "Mildly Reduced" else "Normal",
        if (input$p_ava < 1.0) "Severe" else if (input$p_ava < 1.5) "Moderate" else "Mild/Normal"
      )
    )
  }, striped = TRUE, hover = TRUE, bordered = TRUE, width = "100%")

  ## ------ PK Plots -------------------------------------------
  output$plot_pk_statin <- renderPlotly({
    d <- sim_data() %>% filter(year <= 0.25)
    p <- ggplot(d, aes(x = time/24, y = Cp_statin_ngmL)) +
      geom_line(color = "#E53935", linewidth = 1.2) +
      labs(x = "Days", y = "Rosuvastatin (ng/mL)", title = "") +
      theme_bw(base_size = 11)
    ggplotly(p)
  })

  output$plot_pk_pcsk9i <- renderPlotly({
    d <- sim_data() %>% filter(year <= 0.25)
    p <- ggplot(d, aes(x = time/24, y = Cp_pcsk9i_ugmL)) +
      geom_line(color = "#8E24AA", linewidth = 1.2) +
      labs(x = "Days", y = "Evolocumab (μg/mL)", title = "") +
      theme_bw(base_size = 11)
    ggplotly(p)
  })

  output$plot_pk_deno <- renderPlotly({
    d <- sim_data() %>% filter(year <= 1.0)
    p <- ggplot(d, aes(x = time/24, y = Cp_deno_ugmL)) +
      geom_line(color = "#1E88E5", linewidth = 1.2) +
      labs(x = "Days", y = "Denosumab (μg/mL)", title = "") +
      theme_bw(base_size = 11)
    ggplotly(p)
  })

  output$tbl_pk_summary <- renderDT({
    d <- sim_data() %>% filter(year >= (input$p_years - 0.1))
    lv <- tail(d, 1)
    tibble(
      Drug = c("Rosuvastatin", "Evolocumab", "Denosumab"),
      `Cmax (trough)` = c(
        sprintf("%.1f ng/mL", lv$Cp_statin_ngmL),
        sprintf("%.2f μg/mL", lv$Cp_pcsk9i_ugmL),
        sprintf("%.2f μg/mL", lv$Cp_deno_ugmL)
      ),
      `t½` = c("~19 hr", "~11 days", "~26 days"),
      `Route` = c("Oral QD", "SC Q2W", "SC Q6M")
    )
  }, options = list(dom = 't'), rownames = FALSE)

  ## ------ PD Plots -------------------------------------------
  output$plot_pd_ldl <- renderPlotly({
    d <- sim_data()
    p <- ggplot(d, aes(x = year, y = LDL_mmolL)) +
      geom_line(color = "#E53935", linewidth = 1.3) +
      geom_hline(yintercept = 1.8, linetype = "dashed", color = "#1E88E5") +
      annotate("text", x = max(d$year)*0.05, y = 1.95,
               label = "ESC Target", size = 3, color = "#1E88E5") +
      labs(x = "Year", y = "LDL-C (mmol/L)") +
      theme_bw(base_size = 11)
    ggplotly(p)
  })

  output$plot_pd_lpa <- renderPlotly({
    d <- sim_data()
    p <- ggplot(d, aes(x = year, y = LPA_mgdL)) +
      geom_line(color = "#FB8C00", linewidth = 1.3) +
      geom_hline(yintercept = 50, linetype = "dashed", color = "#E53935") +
      labs(x = "Year", y = "Lp(a) (mg/dL)") +
      theme_bw(base_size = 11)
    ggplotly(p)
  })

  output$plot_pd_rankl <- renderPlotly({
    d <- sim_data()
    p <- ggplot(d, aes(x = year, y = RANKL_norm)) +
      geom_line(color = "#8E24AA", linewidth = 1.3) +
      labs(x = "Year", y = "RANKL Activity (normalized)") +
      theme_bw(base_size = 11)
    ggplotly(p)
  })

  output$plot_pd_mgp_il6 <- renderPlotly({
    d <- sim_data()
    p <- ggplot(d) +
      geom_line(aes(x = year, y = MGP_f, color = "MGP Carboxylation (fraction)"),
                linewidth = 1.2) +
      geom_line(aes(x = year, y = IL6_pgmL / 10, color = "IL-6/10 (pg/mL)"),
                linewidth = 1.2) +
      scale_color_manual(values = c("MGP Carboxylation (fraction)" = "#43A047",
                                    "IL-6/10 (pg/mL)" = "#E53935")) +
      labs(x = "Year", y = "Value", color = "") +
      theme_bw(base_size = 11)
    ggplotly(p)
  })

  ## ------ Clinical Endpoints ---------------------------------
  output$plot_clin_ava <- renderPlotly({
    d <- sim_data()
    p <- ggplot(d, aes(x = year, y = AVA_cm2)) +
      geom_line(color = "#1E88E5", linewidth = 1.5) +
      geom_hline(yintercept = c(1.5, 1.0), linetype = "dashed",
                 color = c("#FF9800", "#F44336")) +
      labs(x = "Year", y = "AVA (cm²)") +
      theme_bw(base_size = 12)
    ggplotly(p)
  })

  output$plot_clin_pg <- renderPlotly({
    d <- sim_data()
    p <- ggplot(d, aes(x = year, y = MeanPG_mmHg)) +
      geom_line(color = "#E53935", linewidth = 1.5) +
      geom_hline(yintercept = c(25, 40), linetype = "dashed",
                 color = c("#FF9800", "#B71C1C")) +
      labs(x = "Year", y = "Mean Gradient (mmHg)") +
      theme_bw(base_size = 12)
    ggplotly(p)
  })

  output$plot_clin_lv <- renderPlotly({
    d <- sim_data()
    p <- ggplot(d) +
      geom_line(aes(x = year, y = LVEF_pct, color = "LVEF (%)"), linewidth = 1.3) +
      geom_line(aes(x = year, y = LVMI_gm2 / 2, color = "LVMI/2 (g/m²)"),
                linewidth = 1.3) +
      scale_color_manual(values = c("LVEF (%)" = "#4CAF50", "LVMI/2 (g/m²)" = "#FF5722")) +
      labs(x = "Year", y = "Value", color = "") +
      theme_bw(base_size = 12)
    ggplotly(p)
  })

  output$plot_clin_bnp <- renderPlotly({
    d <- sim_data()
    p <- ggplot(d) +
      geom_line(aes(x = year, y = NTproBNP_pg, color = "NT-proBNP (pg/mL)"),
                linewidth = 1.3) +
      scale_y_log10() +
      geom_hline(yintercept = c(125, 900), linetype = "dashed",
                 color = c("#FF9800", "#F44336")) +
      scale_color_manual(values = c("NT-proBNP (pg/mL)" = "#9C27B0")) +
      labs(x = "Year", y = "NT-proBNP (log scale)") +
      theme_bw(base_size = 12)
    ggplotly(p)
  })

  output$time_to_severe <- renderUI({
    d <- sim_data()
    time_severe <- d %>% filter(AVA_cm2 < 1.0) %>% slice(1) %>% pull(year)
    if (length(time_severe) == 0 || is.na(time_severe)) {
      HTML('<div class="severity-box severity-normal"><h4>Severe AS Not Reached</h4>
            <p>Within simulation window</p></div>')
    } else {
      HTML(sprintf('<div class="severity-box severity-severe">
        <h4>Time to Severe AS: %.1f years</h4>
        <p>AVA < 1.0 cm² — AVR indicated</p></div>', time_severe))
    }
  })

  output$time_to_avr <- renderUI({
    d <- sim_data()
    avr_time <- d %>% filter(NYHA_class >= 3 | (AVA_cm2 < 1.0 & NTproBNP_pg > 600)) %>%
      slice(1) %>% pull(year)
    if (length(avr_time) == 0 || is.na(avr_time)) {
      HTML('<div class="severity-box severity-normal"><h4>AVR Threshold Not Reached</h4></div>')
    } else {
      HTML(sprintf('<div class="severity-box severity-moderate">
        <h4>AVR Indication: ~%.1f years</h4>
        <p>Based on NYHA III+ or severe AS + elevated BNP</p></div>', avr_time))
    }
  })

  ## ------ Scenario Comparison Plots -------------------------
  make_scen_plot <- function(y_var, y_lab, y_lines = NULL, y_labels = NULL,
                              log_y = FALSE) {
    d <- scenario_data()
    p <- ggplot(d, aes_string(x = "year", y = y_var, color = "scenario")) +
      geom_line(linewidth = 1.0) +
      scale_color_manual(values = scen_colors) +
      labs(x = "Year", y = y_lab, color = "") +
      theme_bw(base_size = 10) +
      theme(legend.position = "none")
    if (!is.null(y_lines)) {
      for (i in seq_along(y_lines)) {
        p <- p + geom_hline(yintercept = y_lines[i], linetype = "dashed", alpha = 0.5)
      }
    }
    if (log_y) p <- p + scale_y_log10()
    ggplotly(p)
  }

  output$plot_scen_ava  <- renderPlotly(make_scen_plot("AVA_cm2",     "AVA (cm²)",     c(1.5,1.0)))
  output$plot_scen_cs   <- renderPlotly(make_scen_plot("CS_au",       "Calcium Score", c(300,800), log_y=FALSE))
  output$plot_scen_lvef <- renderPlotly(make_scen_plot("LVEF_pct",    "LVEF (%)",      c(50,40)))
  output$plot_scen_pg   <- renderPlotly(make_scen_plot("MeanPG_mmHg", "Mean PG (mmHg)",c(25,40)))
  output$plot_scen_lvmi <- renderPlotly(make_scen_plot("LVMI_gm2",    "LVMI (g/m²)",   c(95,115)))
  output$plot_scen_bnp  <- renderPlotly(make_scen_plot("NTproBNP_pg", "NT-proBNP",     log_y=TRUE))

  output$tbl_scenario_comparison <- renderDT({
    d <- scenario_data() %>%
      group_by(scenario) %>%
      slice_tail(n = 1) %>%
      ungroup() %>%
      transmute(
        Scenario   = scenario,
        `AVA (cm²)` = round(AVA_cm2, 2),
        `PG (mmHg)` = round(MeanPG_mmHg, 1),
        `Ca Score`  = round(CS_au, 0),
        `LVEF (%)`  = round(LVEF_pct, 1),
        `LVMI (g/m²)` = round(LVMI_gm2, 1),
        `NT-proBNP` = round(NTproBNP_pg, 0),
        `NYHA`      = round(NYHA_class, 0),
        `LDL (mmol/L)` = round(LDL_mmolL, 2)
      )
    datatable(d, options = list(pageLength = 10, dom = 't'),
              rownames = FALSE) %>%
      formatStyle("AVA (cm²)",
                  backgroundColor = styleInterval(c(1.0, 1.5), c("#FFCDD2","#FFF8E1","#E8F5E9")))
  })

  ## ------ Biomarker Panel ------------------------------------
  output$plot_bio_cs <- renderPlotly({
    d <- sim_data()
    p <- ggplot(d, aes(x = year, y = CS_au)) +
      geom_line(color = "#E53935", linewidth = 1.4) +
      geom_hline(yintercept = c(100, 300, 800, 2000),
                 linetype = "dashed", color = c("#9E9E9E","#FF9800","#F44336","#880E4F"),
                 alpha = 0.7) +
      scale_y_continuous(labels = comma) +
      labs(x = "Year", y = "Calcium Score (AU)") +
      theme_bw(base_size = 11)
    ggplotly(p)
  })

  output$plot_bio_collagen <- renderPlotly({
    d <- sim_data()
    p <- ggplot(d, aes(x = year, y = Collagen_f * 100)) +
      geom_line(color = "#795548", linewidth = 1.4) +
      geom_hline(yintercept = c(10, 20), linetype = "dashed",
                 color = c("#FF9800", "#F44336"), alpha = 0.7) +
      labs(x = "Year", y = "Collagen Fraction (%)") +
      theme_bw(base_size = 11)
    ggplotly(p)
  })

  output$plot_bio_heatmap <- renderPlotly({
    d <- scenario_data() %>%
      filter(year %in% round(seq(1, input$p_years, length.out = 5))) %>%
      group_by(scenario, year) %>%
      summarize(
        AVA   = mean(AVA_cm2),
        CS    = mean(CS_au),
        LVEF  = mean(LVEF_pct),
        BNP   = mean(NTproBNP_pg),
        LVMI  = mean(LVMI_gm2),
        .groups = "drop"
      ) %>%
      pivot_longer(c(AVA, CS, LVEF, BNP, LVMI), names_to = "Biomarker")

    p <- ggplot(d, aes(x = scenario, y = Biomarker, fill = value)) +
      geom_tile(color = "white") +
      facet_wrap(~paste("Year", round(year)), ncol = 5) +
      scale_fill_viridis_c(option = "plasma") +
      labs(x = "", y = "", fill = "Value") +
      theme_bw(base_size = 9) +
      theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 7))
    ggplotly(p)
  })

  output$tbl_biomarker_ref <- renderDT({
    tibble(
      Biomarker     = c("AVA", "Mean PG", "LVEF", "LVMI (M/F)", "NT-proBNP", "Ca Score"),
      Normal        = c(">2.0 cm²", "<20 mmHg", ">50%", "<95/<75 g/m²", "<125 pg/mL", "0-100 AU"),
      Mild          = c("1.5-2.0", "20-25", "50-55%", "95-115/<75-95", "125-400", "100-300"),
      Moderate      = c("1.0-1.5", "25-40", "40-50%", "115-145/95-115", "400-900", "300-800"),
      Severe        = c("<1.0", ">40", "<40%", ">145/>115", ">900", ">800"),
      Reference     = c("AHA/ACC 2021","AHA/ACC 2021","Lang 2015",
                        "Devereux 2004","ESC HF 2021","Baumgartner 2017")
    )
  }, options = list(dom = 't'), rownames = FALSE)

  output$plot_bio_angii_lv <- renderPlotly({
    d <- sim_data()
    p <- ggplot(d) +
      geom_line(aes(x = year, y = AngII_norm, color = "AngII (norm)"), linewidth = 1.2) +
      geom_line(aes(x = year, y = LVMI_gm2/100, color = "LVMI/100 (g/m²)"), linewidth = 1.2) +
      scale_color_manual(values = c("AngII (norm)" = "#E53935", "LVMI/100 (g/m²)" = "#1E88E5")) +
      labs(x = "Year", y = "Value (normalized)", color = "") +
      theme_bw(base_size = 12)
    ggplotly(p)
  })

  ## ------ Severity boxes -------------------------------------
  output$severity_box_ava <- renderUI({
    lv <- last_vals()
    ava <- round(lv$AVA_cm2, 2)
    cls <- if (ava >= 1.5) "severity-normal" else if (ava >= 1.0) "severity-mild" else "severity-severe"
    HTML(sprintf('<div class="severity-box %s"><b>AVA</b><br><h3>%.2f cm²</h3></div>', cls, ava))
  })
  output$severity_box_pg <- renderUI({
    lv <- last_vals()
    pg <- round(lv$MeanPG_mmHg, 1)
    cls <- if (pg < 25) "severity-normal" else if (pg < 40) "severity-mild" else "severity-severe"
    HTML(sprintf('<div class="severity-box %s"><b>Mean PG</b><br><h3>%.1f mmHg</h3></div>', cls, pg))
  })
  output$severity_box_lvmi <- renderUI({
    lv <- last_vals()
    lvmi <- round(lv$LVMI_gm2, 0)
    cls <- if (lvmi < 95) "severity-normal" else if (lvmi < 115) "severity-mild" else "severity-severe"
    HTML(sprintf('<div class="severity-box %s"><b>LVMI</b><br><h3>%d g/m²</h3></div>', cls, lvmi))
  })
  output$severity_box_bnp <- renderUI({
    lv <- last_vals()
    bnp <- round(lv$NTproBNP_pg, 0)
    cls <- if (bnp < 125) "severity-normal" else if (bnp < 900) "severity-mild" else "severity-severe"
    HTML(sprintf('<div class="severity-box %s"><b>NT-proBNP</b><br><h3>%d pg/mL</h3></div>', cls, bnp))
  })

  output$prognosis_text <- renderUI({
    d <- sim_data()
    time_severe <- d %>% filter(AVA_cm2 < 1.0) %>% slice(1) %>% pull(year)
    time_death  <- d %>% filter(NYHA_class >= 4) %>% slice(1) %>% pull(year)

    msgs <- character(0)
    if (length(time_severe) > 0 && !is.na(time_severe))
      msgs <- c(msgs, sprintf("⚠️ Severe AS in %.1f years", time_severe))
    else
      msgs <- c(msgs, "✅ Severe AS not reached in simulation window")

    if (length(time_death) > 0 && !is.na(time_death))
      msgs <- c(msgs, sprintf("🚨 NYHA IV in %.1f years", time_death))

    HTML(paste(msgs, collapse = "<br>"))
  })

  output$tbl_pd_effects <- renderDT({
    lv <- last_vals()
    tibble(
      Parameter        = c("LDL-C", "Lp(a)", "RANKL", "IL-6", "MGP Carboxylation",
                            "AngII", "Ca Score", "LVEF", "NT-proBNP"),
      `Current Value`  = c(round(lv$LDL_mmolL, 2),
                            round(lv$LPA_mgdL, 1),
                            round(lv$RANKL_norm, 2),
                            round(lv$IL6_pgmL, 1),
                            round(lv$MGP_f, 2),
                            round(lv$AngII_norm, 2),
                            round(lv$CS_au, 0),
                            round(lv$LVEF_pct, 1),
                            round(lv$NTproBNP_pg, 0)),
      Unit             = c("mmol/L", "mg/dL", "norm", "pg/mL", "fraction",
                            "norm", "AU", "%", "pg/mL"),
      `Normal Range`   = c("<1.8", "<50", "1.0", "<3.0", ">0.7",
                            "1.0", "<100", ">50", "<125")
    )
  }, options = list(dom = 't'), rownames = FALSE)

  # Auto-run on startup
  observeEvent(TRUE, {
    shinyjs::click("run_sim")
  }, once = TRUE)
}

## ============================================================
## RUN APP
## ============================================================
shinyApp(ui = ui, server = server)
