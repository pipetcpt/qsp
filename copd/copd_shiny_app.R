## =============================================================================
## COPD QSP Interactive Shiny Dashboard
## =============================================================================
## Tabs:
##   1. Patient Profile & GOLD Classification
##   2. Drug Pharmacokinetics (LAMA / LABA / ICS / PDE4i)
##   3. Lung Function (FEV1, FVC, GOLD Spirometry)
##   4. Inflammatory Biomarkers (IL-8, NE, CRP, Eosinophils)
##   5. Clinical Endpoints (CAT, mMRC, Exacerbations, BODE)
##   6. Scenario Comparison (All 6 treatment arms)
##   7. Exacerbation Risk Calculator
## =============================================================================

library(shiny)
library(shinydashboard)
library(mrgsolve)
library(dplyr)
library(ggplot2)
library(plotly)
library(DT)
library(tidyr)
library(scales)

# ── Embedded mrgsolve model code ─────────────────────────────────────────────
copd_code <- '
$PARAM @annotated
ka_LAMA:0.5:LAMA lung absorption (1/h)
F_lung_LAMA:0.20:LAMA lung fraction
F_sys_LAMA:0.03:LAMA systemic fraction
CL_LAMA:3.92:LAMA clearance (L/h)
Vc_LAMA:36:LAMA central volume (L)
Q_LAMA:0.5:LAMA inter-CL (L/h)
Vp_LAMA:8:LAMA peripheral (L)
EC50_LAMA:0.10:LAMA EC50 (ng/mL lung)
Emax_LAMA:0.18:LAMA max FEV1 gain
ka_LABA:1.2:LABA lung absorption (1/h)
F_lung_LABA:0.20:LABA lung fraction
F_sys_LABA:0.05:LABA systemic fraction
CL_LABA:23:LABA clearance (L/h)
Vc_LABA:45:LABA central (L)
Q_LABA:2:LABA inter-CL (L/h)
Vp_LABA:12:LABA peripheral (L)
EC50_LABA:0.08:LABA EC50 (ng/mL lung)
Emax_LABA:0.15:LABA max FEV1 gain
ka_ICS:0.60:ICS lung absorption (1/h)
F_lung_ICS:0.22:ICS lung fraction
F_sys_ICS:0.10:ICS systemic fraction
CL_ICS:78:ICS clearance (L/h)
Vc_ICS:180:ICS central (L)
Q_ICS:8:ICS inter-CL (L/h)
Vp_ICS:60:ICS peripheral (L)
EC50_ICS:0.15:ICS EC50 (ng/mL lung)
Emax_ICS:0.30:ICS max anti-inflam
Hill_ICS:1.5:ICS Hill coefficient
Eos_thresh:300:Eos threshold (cells/uL)
ka_PDE4i:0.70:PDE4i absorption (1/h)
F_PDE4i:0.80:PDE4i bioavailability
CL_PDE4i:9.1:PDE4i clearance (L/h)
Vc_PDE4i:210:PDE4i central (L)
Emax_PDE4i:0.20:PDE4i max anti-inflam
EC50_PDE4i:2.5:PDE4i EC50 (ng/mL)
ksyn_IL8:15:IL-8 synthesis (pg/mL/h)
kout_IL8:0.12:IL-8 elimination (1/h)
IL8_0:125:Baseline IL-8 (pg/mL)
ksyn_NE:4:NE synthesis (ug/mL/h)
kout_NE:0.08:NE elimination (1/h)
NE_0:50:Baseline NE (ug/mL)
ksyn_CRP:0.10:CRP synthesis (mg/L/h)
kout_CRP:0.020:CRP elimination (1/h)
CRP_0:5:Baseline CRP (mg/L)
ksyn_Eos:2:Eos production (cells/uL/h)
kout_Eos:0.003:Eos removal (1/h)
Eos_0:200:Baseline Eos (cells/uL)
FEV1_0:55:Baseline FEV1 pct pred
k_FEV1_dec:0.000055:FEV1 decline rate (/h)
k_AE_FEV1:0.015:FEV1 loss per AE
k_inflam_FEV1:0.0002:Inflam FEV1 penalty
Emph_0:20:Baseline emphysema (%)
k_emph_prog:0.000020:Emphysema progression
lambda_AE_0:1.8:Baseline AE rate (/yr)
k_AE:0.000205:AE accumulation
k_post_AE:0.50:Post-AE rate increase
PVR_0:280:Baseline PVR (dyn.s.cm-5)
mPAP_0:26:Baseline mPAP (mmHg)
k_PVR_prog:0.000015:PVR progression
CO_0:4.5:Cardiac output (L/min)
PAWP_0:8:Wedge pressure (mmHg)
DOSE_LAMA:0:LAMA active (0/1)
DOSE_LABA:0:LABA active (0/1)
DOSE_ICS:0:ICS active (0/1)
DOSE_PDE4i:0:PDE4i active (0/1)

$CMT LAMA_LUNG LAMA_C LAMA_P LABA_LUNG LABA_C LABA_P ICS_LUNG ICS_C ICS_P PDE4i_C
     IL8 NE_sput CRP Eos FEV1 Emph PVR AE_cum AE_rate_ann

$INIT LAMA_LUNG=0,LAMA_C=0,LAMA_P=0,LABA_LUNG=0,LABA_C=0,LABA_P=0,
      ICS_LUNG=0,ICS_C=0,ICS_P=0,PDE4i_C=0,
      IL8=125,NE_sput=50,CRP=5,Eos=200,FEV1=55,Emph=20,PVR=280,
      AE_cum=0,AE_rate_ann=1.8

$ODE
double LAMA_conc_lung=LAMA_LUNG/(F_lung_LAMA*200.0);
dxdt_LAMA_LUNG=-ka_LAMA*LAMA_LUNG;
dxdt_LAMA_C=ka_LAMA*LAMA_LUNG*F_sys_LAMA-(CL_LAMA+Q_LAMA)*(LAMA_C/Vc_LAMA)+Q_LAMA*(LAMA_P/Vp_LAMA);
dxdt_LAMA_P=Q_LAMA*(LAMA_C/Vc_LAMA)-Q_LAMA*(LAMA_P/Vp_LAMA);
double Cp_LAMA=LAMA_C/Vc_LAMA;
double LABA_conc_lung=LABA_LUNG/(F_lung_LABA*200.0);
dxdt_LABA_LUNG=-ka_LABA*LABA_LUNG;
dxdt_LABA_C=ka_LABA*LABA_LUNG*F_sys_LABA-(CL_LABA+Q_LABA)*(LABA_C/Vc_LABA)+Q_LABA*(LABA_P/Vp_LABA);
dxdt_LABA_P=Q_LABA*(LABA_C/Vc_LABA)-Q_LABA*(LABA_P/Vp_LABA);
double Cp_LABA=LABA_C/Vc_LABA;
double ICS_conc_lung=ICS_LUNG/(F_lung_ICS*200.0);
dxdt_ICS_LUNG=-ka_ICS*ICS_LUNG;
dxdt_ICS_C=ka_ICS*ICS_LUNG*F_sys_ICS-(CL_ICS+Q_ICS)*(ICS_C/Vc_ICS)+Q_ICS*(ICS_P/Vp_ICS);
dxdt_ICS_P=Q_ICS*(ICS_C/Vc_ICS)-Q_ICS*(ICS_P/Vp_ICS);
double Cp_ICS=ICS_C/Vc_ICS;
double Cp_PDE4i=PDE4i_C/Vc_PDE4i;
dxdt_PDE4i_C=-CL_PDE4i*Cp_PDE4i;
double E_LAMA=DOSE_LAMA*Emax_LAMA*LAMA_conc_lung/(EC50_LAMA+LAMA_conc_lung);
double E_LABA=DOSE_LABA*Emax_LABA*LABA_conc_lung/(EC50_LABA+LABA_conc_lung);
double Eos_factor=(Eos>Eos_thresh)?1.0:(Eos/Eos_thresh);
double E_ICS_lung=DOSE_ICS*Emax_ICS*pow(ICS_conc_lung,Hill_ICS)/(pow(EC50_ICS,Hill_ICS)+pow(ICS_conc_lung,Hill_ICS));
double E_ICS=E_ICS_lung*Eos_factor;
double E_PDE4i=DOSE_PDE4i*Emax_PDE4i*Cp_PDE4i/(EC50_PDE4i+Cp_PDE4i);
double IL8_stim_factor=(FEV1<50.0)?1.5:1.0;
dxdt_IL8=ksyn_IL8*IL8_stim_factor*(1.0-E_ICS*0.5-E_PDE4i*0.4)-kout_IL8*IL8;
dxdt_NE_sput=ksyn_NE*(IL8/IL8_0)*(1.0-E_ICS*0.3)-kout_NE*NE_sput;
dxdt_CRP=ksyn_CRP*(IL8/IL8_0)*(1.0-E_ICS*0.4)-kout_CRP*CRP;
double Eos_stim=(DOSE_ICS==1.0)?0.3:1.0;
dxdt_Eos=ksyn_Eos*Eos_stim-kout_Eos*Eos;
double AE_rate_h=(lambda_AE_0*(1.0-E_ICS*0.35-E_PDE4i*0.15))/(24.0*365.0);
double E_BD=E_LAMA+E_LABA-0.3*E_LAMA*E_LABA;
double FEV1_target=FEV1_0*(1.0+E_BD);
double inflam_penalty=k_inflam_FEV1*(IL8/IL8_0-1.0)*FEV1;
dxdt_FEV1=0.05*(FEV1_target-FEV1)-k_FEV1_dec*FEV1-inflam_penalty-k_AE_FEV1*AE_rate_h*FEV1;
dxdt_Emph=k_emph_prog*(NE_sput/NE_0)*(1.0-E_ICS*0.05)*(100.0-Emph);
double FEV1_frac=(FEV1>0.1?FEV1:0.1)/FEV1_0;
dxdt_PVR=k_PVR_prog*(1.0/FEV1_frac-1.0)*PVR_0-k_PVR_prog*0.1*(PVR-PVR_0);
dxdt_AE_cum=AE_rate_h;
dxdt_AE_rate_ann=0.001*(lambda_AE_0*(1.0-E_ICS*0.35-E_PDE4i*0.15)*(1.0+k_post_AE*AE_cum/100.0)-AE_rate_ann);

$TABLE
double mPAP=CO_0*PVR/80.0+PAWP_0;
double GOLD_stage=(FEV1>=80)?1:(FEV1>=50)?2:(FEV1>=30)?3:4;
double CAT_approx=10.0+(55.0-FEV1)*0.35+(CRP/5.0)*1.2;
if(CAT_approx>40) CAT_approx=40.0;
if(CAT_approx<0)  CAT_approx=0.0;
double SpO2=98.0-0.12*(PVR/280.0-1.0)*5.0-0.15*(100.0-FEV1)/50.0;
if(SpO2>99.0) SpO2=99.0;
if(SpO2<85.0) SpO2=85.0;
double E_BD_cap=E_LAMA+E_LABA-0.3*E_LAMA*E_LABA;

$CAPTURE Cp_LAMA Cp_LABA Cp_ICS Cp_PDE4i LAMA_conc_lung LABA_conc_lung ICS_conc_lung
         IL8 NE_sput CRP Eos FEV1 Emph PVR AE_cum AE_rate_ann
         mPAP GOLD_stage CAT_approx SpO2 E_LAMA E_LABA E_ICS E_PDE4i E_BD_cap
'

mod_base <- mcode("COPD_SHINY", copd_code, quiet = TRUE)

# ── Helper functions ─────────────────────────────────────────────────────────
build_events <- function(use_lama, use_laba, use_ics, use_pde4i) {
  evs <- NULL
  if (use_lama) evs <- c(evs, ev(cmt="LAMA_LUNG", amt=3600, ii=24, addl=364*2))
  if (use_laba) evs <- c(evs, ev(cmt="LABA_LUNG", amt=10000, ii=12, addl=728*2))
  if (use_ics)  evs <- c(evs, ev(cmt="ICS_LUNG",  amt=88000, ii=12, addl=728*2))
  if (use_pde4i) evs <- c(evs, ev(cmt="PDE4i_C",  amt=400000, ii=24, addl=364*2))
  evs
}

run_sim <- function(p_list, evs, duration_days = 365) {
  mod_tmp <- param(mod_base, p_list)
  tg <- seq(0, duration_days * 24, by = 6)
  if (is.null(evs)) {
    out <- mrgsim(mod_tmp, tgrid = tg, delta = 6)
  } else {
    out <- mrgsim(mod_tmp, events = evs, tgrid = tg, delta = 6)
  }
  as.data.frame(out) %>% mutate(time_day = time / 24, time_wk = time / 168)
}

theme_shiny <- theme_bw(base_size = 12) +
  theme(
    panel.grid.minor = element_blank(),
    legend.position  = "bottom"
  )

gold_color <- function(stage) {
  c("1"="#1B5E20","2"="#F57F17","3"="#E65100","4"="#B71C1C")[as.character(stage)]
}

# ── UI ───────────────────────────────────────────────────────────────────────
ui <- dashboardPage(
  skin = "blue",

  dashboardHeader(
    title = span(icon("lungs"), "COPD QSP Dashboard"),
    titleWidth = 300
  ),

  dashboardSidebar(
    width = 280,
    sidebarMenu(
      menuItem("Patient Profile",         tabName = "tab_patient",  icon = icon("user-md")),
      menuItem("Pharmacokinetics",        tabName = "tab_pk",       icon = icon("flask")),
      menuItem("Lung Function",           tabName = "tab_lung",     icon = icon("wind")),
      menuItem("Inflammatory Biomarkers", tabName = "tab_inflam",   icon = icon("vial")),
      menuItem("Clinical Endpoints",      tabName = "tab_clinical", icon = icon("chart-line")),
      menuItem("Scenario Comparison",     tabName = "tab_compare",  icon = icon("th")),
      menuItem("Exacerbation Risk",       tabName = "tab_ae",       icon = icon("exclamation-triangle"))
    ),
    hr(),

    h5("Patient Parameters", style="color:#90CAF9; padding-left:15px;"),
    sliderInput("fev1_0",    "Baseline FEV1 (% pred):",  min=20, max=80,  value=55, step=1),
    sliderInput("eos_0",     "Blood Eos (cells/µL):",    min=50, max=800, value=200, step=25),
    sliderInput("ae_rate_0", "Baseline AE rate (/yr):",  min=0.5, max=5, value=1.8, step=0.1),
    sliderInput("emph_0",    "Baseline Emphysema (%):",  min=5, max=60,  value=20, step=5),
    sliderInput("duration",  "Simulation duration (yr):", min=0.5, max=2, value=1, step=0.5),

    hr(),
    h5("Treatment Selection", style="color:#90CAF9; padding-left:15px;"),
    checkboxInput("use_lama",  "LAMA (Tiotropium 18µg/d)",      value = TRUE),
    checkboxInput("use_laba",  "LABA (Salmeterol 50µg bid)",     value = FALSE),
    checkboxInput("use_ics",   "ICS (Budesonide 400µg bid)",     value = FALSE),
    checkboxInput("use_pde4i", "PDE4i (Roflumilast 500µg/d)",   value = FALSE),
    br(),
    actionButton("run_btn", "Run Simulation", icon = icon("play-circle"),
                 style = "color:#fff; background-color:#1565C0; border-color:#0D47A1; width:100%; font-weight:bold;")
  ),

  dashboardBody(
    tags$head(tags$style(HTML("
      .content-wrapper, .right-side { background-color: #F5F7FA; }
      .box.box-primary { border-top-color: #1565C0; }
      .value-box .inner { padding: 10px 15px; }
    "))),

    tabItems(

      # ── TAB 1: Patient Profile ─────────────────────────────────────────────
      tabItem(tabName = "tab_patient",
        fluidRow(
          valueBoxOutput("vbox_gold",   width = 3),
          valueBoxOutput("vbox_fev1",   width = 3),
          valueBoxOutput("vbox_eos",    width = 3),
          valueBoxOutput("vbox_cat",    width = 3)
        ),
        fluidRow(
          box(title = "GOLD ABE Group Classification",
              status = "primary", solidHeader = TRUE, width = 6,
              plotlyOutput("plot_gold_gauge", height = "320px")),
          box(title = "Patient Phenotype Summary",
              status = "primary", solidHeader = TRUE, width = 6,
              tableOutput("tbl_profile"))
        ),
        fluidRow(
          box(title = "GOLD 2023 ABCD → ABE Group Framework",
              status = "info", solidHeader = TRUE, width = 12,
              HTML('<table class="table table-bordered table-sm" style="font-size:0.9em;">
                <thead class="thead-dark">
                  <tr><th>GOLD ABE Group</th><th>mMRC Dyspnea</th><th>CAT Score</th>
                      <th>Exacerbations/yr</th><th>Treatment Initiation</th></tr>
                </thead><tbody>
                  <tr><td><b>Group A</b></td><td>0-1</td><td>&lt;10</td><td>0-1 (no hosp.)</td>
                      <td>Bronchodilator monotherapy</td></tr>
                  <tr><td><b>Group B</b></td><td>≥2</td><td>≥10</td><td>0-1 (no hosp.)</td>
                      <td>LAMA/LABA dual therapy</td></tr>
                  <tr><td><b>Group E</b></td><td>Any</td><td>Any</td><td>≥2 or ≥1 hosp.</td>
                      <td>LAMA/LABA ± ICS (if Eos≥300); consider Roflumilast</td></tr>
                </tbody></table>'))
        )
      ),

      # ── TAB 2: Pharmacokinetics ────────────────────────────────────────────
      tabItem(tabName = "tab_pk",
        fluidRow(
          box(title = "Drug Lung Concentration (ng/mL) — First 48 Hours",
              status = "primary", solidHeader = TRUE, width = 6,
              plotlyOutput("plot_pk_lung", height = "380px")),
          box(title = "Drug Systemic Plasma Concentration — First 48 Hours",
              status = "primary", solidHeader = TRUE, width = 6,
              plotlyOutput("plot_pk_sys", height = "380px"))
        ),
        fluidRow(
          box(title = "Steady-State Lung Concentrations & PD Effect Estimates",
              status = "info", solidHeader = TRUE, width = 12,
              tableOutput("tbl_pk_ss"))
        ),
        fluidRow(
          box(title = "PK Parameters Summary",
              status = "warning", solidHeader = TRUE, width = 12,
              HTML('<table class="table table-sm table-bordered">
              <thead class="thead-dark"><tr><th>Drug</th><th>F_lung</th><th>F_sys</th>
                <th>CL (L/h)</th><th>Vc (L)</th><th>t½ (h)</th><th>Route</th><th>Dosing</th></tr></thead>
              <tbody>
                <tr><td>Tiotropium (LAMA)</td><td>20%</td><td>3%</td><td>3.92</td>
                    <td>36</td><td>~120h (5d)</td><td>DPI (HandiHaler)</td><td>18 µg qd</td></tr>
                <tr><td>Salmeterol (LABA)</td><td>20%</td><td>5%</td><td>23</td>
                    <td>45</td><td>~5.5h</td><td>MDI/DPI</td><td>50 µg bid</td></tr>
                <tr><td>Budesonide (ICS)</td><td>22%</td><td>10%</td><td>78</td>
                    <td>180</td><td>~2.8h</td><td>MDI/DPI</td><td>400 µg bid</td></tr>
                <tr><td>Roflumilast (PDE4i)</td><td>—</td><td>80%</td><td>9.1</td>
                    <td>210</td><td>~17h</td><td>Oral</td><td>500 µg qd</td></tr>
              </tbody></table>'))
        )
      ),

      # ── TAB 3: Lung Function ───────────────────────────────────────────────
      tabItem(tabName = "tab_lung",
        fluidRow(
          box(title = "FEV1 (% predicted) Over Time",
              status = "primary", solidHeader = TRUE, width = 8,
              plotlyOutput("plot_fev1", height = "380px")),
          box(title = "GOLD Stage Distribution",
              status = "primary", solidHeader = TRUE, width = 4,
              plotlyOutput("plot_gold_bar", height = "380px"))
        ),
        fluidRow(
          box(title = "Emphysema Index (% LAA) — Structural Damage",
              status = "warning", solidHeader = TRUE, width = 6,
              plotlyOutput("plot_emph", height = "320px")),
          box(title = "Pulmonary Vascular Resistance (PVR)",
              status = "warning", solidHeader = TRUE, width = 6,
              plotlyOutput("plot_pvr", height = "320px"))
        )
      ),

      # ── TAB 4: Inflammatory Biomarkers ─────────────────────────────────────
      tabItem(tabName = "tab_inflam",
        fluidRow(
          box(title = "Sputum IL-8 (pg/mL)",
              status = "danger", solidHeader = TRUE, width = 6,
              plotlyOutput("plot_il8", height = "320px")),
          box(title = "Sputum Neutrophil Elastase (µg/mL)",
              status = "danger", solidHeader = TRUE, width = 6,
              plotlyOutput("plot_ne", height = "320px"))
        ),
        fluidRow(
          box(title = "Serum CRP (mg/L)",
              status = "warning", solidHeader = TRUE, width = 6,
              plotlyOutput("plot_crp", height = "320px")),
          box(title = "Blood Eosinophils (cells/µL)",
              status = "info", solidHeader = TRUE, width = 6,
              plotlyOutput("plot_eos", height = "320px"))
        ),
        fluidRow(
          box(title = "ICS Eosinophil-Guided Treatment Rationale",
              status = "info", solidHeader = TRUE, width = 12,
              HTML('<ul>
                <li><b>Eos ≥300 cells/µL:</b> Full ICS benefit — ≥35% exacerbation reduction (IMPACT: RR 0.65)</li>
                <li><b>Eos 100-299 cells/µL:</b> Partial ICS benefit; consider based on exacerbation history</li>
                <li><b>Eos &lt;100 cells/µL:</b> ICS benefit unlikely; may increase pneumonia risk</li>
                <li>Blood Eos is a pharmacodynamic biomarker (predicts ICS response, not disease severity)</li>
              </ul>'))
        )
      ),

      # ── TAB 5: Clinical Endpoints ──────────────────────────────────────────
      tabItem(tabName = "tab_clinical",
        fluidRow(
          box(title = "CAT Score Approximation",
              status = "primary", solidHeader = TRUE, width = 6,
              plotlyOutput("plot_cat", height = "320px")),
          box(title = "Cumulative Exacerbations",
              status = "danger", solidHeader = TRUE, width = 6,
              plotlyOutput("plot_ae", height = "320px"))
        ),
        fluidRow(
          box(title = "Annualized Exacerbation Rate",
              status = "warning", solidHeader = TRUE, width = 6,
              plotlyOutput("plot_ae_rate", height = "320px")),
          box(title = "SpO2 Estimate",
              status = "info", solidHeader = TRUE, width = 6,
              plotlyOutput("plot_spo2", height = "320px"))
        ),
        fluidRow(
          box(title = "mPAP (Mean Pulmonary Artery Pressure)",
              status = "danger", solidHeader = TRUE, width = 12,
              plotlyOutput("plot_mpap", height = "280px"))
        )
      ),

      # ── TAB 6: Scenario Comparison ─────────────────────────────────────────
      tabItem(tabName = "tab_compare",
        fluidRow(
          box(title = "FEV1 — All 6 Treatment Scenarios",
              status = "primary", solidHeader = TRUE, width = 12,
              plotlyOutput("plot_compare_fev1", height = "380px"))
        ),
        fluidRow(
          box(title = "IL-8 — All Scenarios",
              status = "danger", solidHeader = TRUE, width = 6,
              plotlyOutput("plot_compare_il8", height = "320px")),
          box(title = "Cumulative AEs — All Scenarios",
              status = "warning", solidHeader = TRUE, width = 6,
              plotlyOutput("plot_compare_ae", height = "320px"))
        ),
        fluidRow(
          box(title = "End-of-Year Summary Table",
              status = "info", solidHeader = TRUE, width = 12,
              DTOutput("tbl_compare"))
        )
      ),

      # ── TAB 7: Exacerbation Risk ───────────────────────────────────────────
      tabItem(tabName = "tab_ae",
        fluidRow(
          box(title = "AECOPD Risk Factors & Model Inputs",
              status = "danger", solidHeader = TRUE, width = 4,
              sliderInput("ae_fev1", "FEV1 % predicted:", min=20, max=80, value=45, step=1),
              sliderInput("ae_eos", "Blood Eos (cells/µL):", min=50, max=800, value=300, step=50),
              sliderInput("ae_history", "Prior AE/yr (history):", min=0, max=6, value=2, step=0.5),
              selectInput("ae_triple", "Treatment:",
                          choices = c("None","LAMA","LABA+LAMA","ICS+LABA",
                                      "Triple (LAMA+LABA+ICS)","Triple+Roflumilast"),
                          selected = "LAMA"),
              actionButton("ae_calc", "Calculate Risk", icon=icon("calculator"),
                           style="background-color:#E53935; color:white; width:100%;")
          ),
          box(title = "Estimated Annual AECOPD Rate",
              status = "danger", solidHeader = TRUE, width = 4,
              plotlyOutput("plot_ae_risk", height = "350px")),
          box(title = "AECOPD Phenotype Guide",
              status = "warning", solidHeader = TRUE, width = 4,
              HTML('<h5>Bacterial vs Viral Triggers</h5>
              <table class="table table-sm">
              <tr><th>Pathogen</th><th>Frequency</th><th>Antibiotic?</th></tr>
              <tr><td>Rhinovirus (RV)</td><td>40-50%</td><td>No (viral)</td></tr>
              <tr><td>H. influenzae</td><td>20-30%</td><td>Amox-clav / Doxyc.</td></tr>
              <tr><td>S. pneumoniae</td><td>10-15%</td><td>Amoxicillin</td></tr>
              <tr><td>M. catarrhalis</td><td>5-10%</td><td>Amox-clav</td></tr>
              <tr><td>P. aeruginosa</td><td>5-10% (severe)</td><td>Cipro/Piperacillin</td></tr>
              </table>
              <hr/>
              <b>Biomarker-guided antibiotic:</b><br/>
              PCT &lt;0.1 ng/mL → antibiotics unlikely needed<br/>
              PCT &gt;0.25 ng/mL → likely bacterial → treat<br/>
              Eos &lt;100 → non-eosinophilic → reduced systemic CS benefit'))
        ),
        fluidRow(
          box(title = "Exacerbation Risk Reduction by Treatment (vs Placebo, GOLD Group E)",
              status = "info", solidHeader = TRUE, width = 12,
              plotlyOutput("plot_rr_waterfall", height = "320px"))
        )
      )
    )
  )
)

# ── SERVER ───────────────────────────────────────────────────────────────────
server <- function(input, output, session) {

  # Reactive: single treatment arm simulation
  sim_data <- eventReactive(input$run_btn, {
    req(input$fev1_0, input$eos_0)
    p_list <- list(
      FEV1_0    = input$fev1_0,
      Eos_0     = input$eos_0,
      lambda_AE_0 = input$ae_rate_0,
      Emph_0    = input$emph_0,
      DOSE_LAMA  = as.numeric(input$use_lama),
      DOSE_LABA  = as.numeric(input$use_laba),
      DOSE_ICS   = as.numeric(input$use_ics),
      DOSE_PDE4i = as.numeric(input$use_pde4i)
    )
    evs <- build_events(input$use_lama, input$use_laba,
                        input$use_ics, input$use_pde4i)
    run_sim(p_list, evs, duration_days = input$duration * 365)
  }, ignoreNULL = FALSE)

  # Reactive: all 6 scenarios (for comparison)
  all_scenarios <- eventReactive(input$run_btn, {
    base_params <- list(
      FEV1_0     = input$fev1_0,
      Eos_0      = input$eos_0,
      lambda_AE_0 = input$ae_rate_0,
      Emph_0     = input$emph_0
    )
    scen_list <- list(
      list(name="1. Placebo",       lama=0, laba=0, ics=0, pde4i=0),
      list(name="2. LAMA",          lama=1, laba=0, ics=0, pde4i=0),
      list(name="3. LABA+LAMA",     lama=1, laba=1, ics=0, pde4i=0),
      list(name="4. ICS+LABA",      lama=0, laba=1, ics=1, pde4i=0),
      list(name="5. Triple",        lama=1, laba=1, ics=1, pde4i=0),
      list(name="6. Triple+PDE4i",  lama=1, laba=1, ics=1, pde4i=1)
    )
    bind_rows(lapply(scen_list, function(s) {
      p <- c(base_params, list(
        DOSE_LAMA=s$lama, DOSE_LABA=s$laba,
        DOSE_ICS=s$ics, DOSE_PDE4i=s$pde4i))
      evs <- build_events(s$lama==1, s$laba==1, s$ics==1, s$pde4i==1)
      run_sim(p, evs, input$duration * 365) %>%
        mutate(scenario = s$name)
    }))
  }, ignoreNULL = FALSE)

  # ── Value Boxes ─────────────────────────────────────────────────────────────
  output$vbox_gold <- renderValueBox({
    d <- sim_data()
    gold <- round(tail(d$GOLD_stage, 1))
    colors <- c("1"="green","2"="yellow","3"="orange","4"="red")
    valueBox(paste0("GOLD ", gold), "Spirometric Stage",
             icon = icon("lungs"), color = colors[as.character(gold)])
  })
  output$vbox_fev1 <- renderValueBox({
    d <- sim_data()
    valueBox(paste0(round(tail(d$FEV1, 1), 1), "%"),
             "FEV1 (% predicted)", icon = icon("chart-bar"), color = "blue")
  })
  output$vbox_eos <- renderValueBox({
    d <- sim_data()
    eos_val <- round(tail(d$Eos, 1))
    col <- ifelse(eos_val >= 300, "green", ifelse(eos_val >= 100, "yellow", "red"))
    valueBox(paste0(eos_val, " cells/µL"),
             "Blood Eosinophils", icon = icon("tint"), color = col)
  })
  output$vbox_cat <- renderValueBox({
    d <- sim_data()
    cat_val <- round(tail(d$CAT_approx, 1), 1)
    col <- ifelse(cat_val < 10, "green", ifelse(cat_val < 20, "yellow",
                  ifelse(cat_val < 30, "orange", "red")))
    valueBox(cat_val, "CAT Score", icon = icon("clipboard-list"), color = col)
  })

  # ── Tab 1: Profile gauge ────────────────────────────────────────────────────
  output$plot_gold_gauge <- renderPlotly({
    d <- sim_data()
    fev_end <- tail(d$FEV1, 1)
    plot_ly(
      type = "indicator", mode = "gauge+number",
      value = round(fev_end, 1),
      title = list(text = "FEV1 % predicted", font = list(size = 16)),
      gauge = list(
        axis = list(range = list(0, 100)),
        bar  = list(color = "#1565C0"),
        steps = list(
          list(range=c(0,30),  color="#FFCDD2"),
          list(range=c(30,50), color="#FFE0B2"),
          list(range=c(50,80), color="#FFF9C4"),
          list(range=c(80,100),color="#C8E6C9")
        ),
        threshold = list(line=list(color="#B71C1C", width=3), thickness=0.75,
                         value = round(fev_end, 1))
      )
    ) %>% layout(margin = list(t=50))
  })

  output$tbl_profile <- renderTable({
    d <- sim_data()
    last <- tail(d, 1)
    data.frame(
      Parameter = c("FEV1 (% pred)", "GOLD Stage", "CAT Score",
                    "Blood Eos (cells/µL)", "CRP (mg/L)", "SpO2 (%)",
                    "AE cum (1 yr)", "AE rate (/yr)",
                    "Emphysema (%)", "PVR (dyn)","mPAP (mmHg)"),
      Value = c(
        round(last$FEV1,1), round(last$GOLD_stage,0), round(last$CAT_approx,1),
        round(last$Eos,0), round(last$CRP,2), round(last$SpO2,1),
        round(last$AE_cum,2), round(last$AE_rate_ann,2),
        round(last$Emph,1), round(last$PVR,0), round(last$mPAP,1)
      )
    )
  }, striped=TRUE, hover=TRUE)

  # ── Tab 2: PK plots ─────────────────────────────────────────────────────────
  output$plot_pk_lung <- renderPlotly({
    d <- sim_data() %>% filter(time_day <= 2)
    p <- plot_ly(d, x=~time_day) %>%
      add_lines(y=~LAMA_conc_lung, name="LAMA lung", line=list(color="#1565C0")) %>%
      add_lines(y=~LABA_conc_lung, name="LABA lung", line=list(color="#2E7D32")) %>%
      add_lines(y=~ICS_conc_lung,  name="ICS lung",  line=list(color="#880E4F")) %>%
      layout(xaxis=list(title="Day"), yaxis=list(title="Lung concentration (ng/mL)"),
             legend=list(orientation="h"))
    p
  })
  output$plot_pk_sys <- renderPlotly({
    d <- sim_data() %>% filter(time_day <= 2)
    p <- plot_ly(d, x=~time_day) %>%
      add_lines(y=~Cp_LAMA,   name="Tiotropium", line=list(color="#1565C0")) %>%
      add_lines(y=~Cp_LABA,   name="Salmeterol", line=list(color="#2E7D32")) %>%
      add_lines(y=~Cp_ICS,    name="Budesonide", line=list(color="#880E4F")) %>%
      add_lines(y=~Cp_PDE4i,  name="Roflumilast",line=list(color="#6A1B9A")) %>%
      layout(xaxis=list(title="Day"), yaxis=list(title="Plasma Cp (ng/mL)"),
             legend=list(orientation="h"))
    p
  })
  output$tbl_pk_ss <- renderTable({
    d <- sim_data() %>% filter(time_day >= (input$duration*365 - 2))
    data.frame(
      Drug = c("Tiotropium (LAMA)","Salmeterol (LABA)","Budesonide (ICS)","Roflumilast (PDE4i)"),
      `Lung Css (ng/mL)` = c(mean(d$LAMA_conc_lung), mean(d$LABA_conc_lung),
                               mean(d$ICS_conc_lung),  NA),
      `Plasma Css (ng/mL)` = c(mean(d$Cp_LAMA), mean(d$Cp_LABA),
                                mean(d$Cp_ICS),  mean(d$Cp_PDE4i)),
      `E_drug (%)` = c(mean(d$E_LAMA)*100, mean(d$E_LABA)*100,
                       mean(d$E_ICS)*100,  mean(d$E_PDE4i)*100),
      check.names = FALSE
    )
  }, digits=3, striped=TRUE)

  # ── Tab 3: Lung Function ────────────────────────────────────────────────────
  output$plot_fev1 <- renderPlotly({
    d <- sim_data()
    plot_ly(d, x=~time_day, y=~FEV1, type="scatter", mode="lines",
            line=list(color="#1565C0", width=2.5)) %>%
      add_lines(y=rep(80, nrow(d)), name="GOLD 1 (80%)", line=list(dash="dash",color="#1B5E20")) %>%
      add_lines(y=rep(50, nrow(d)), name="GOLD 2 (50%)", line=list(dash="dash",color="#F57F17")) %>%
      add_lines(y=rep(30, nrow(d)), name="GOLD 3 (30%)", line=list(dash="dash",color="#E65100")) %>%
      layout(xaxis=list(title="Day"), yaxis=list(title="FEV1 (% predicted)", range=c(10,90)),
             legend=list(orientation="h"))
  })
  output$plot_gold_bar <- renderPlotly({
    d <- sim_data()
    d$gold_f <- factor(round(d$GOLD_stage),
                       levels=c(1,2,3,4),
                       labels=c("GOLD 1","GOLD 2","GOLD 3","GOLD 4"))
    cnt <- as.data.frame(table(d$gold_f))
    plot_ly(cnt, x=~Var1, y=~Freq, type="bar",
            marker=list(color=c("#1B5E20","#F9A825","#E65100","#B71C1C"))) %>%
      layout(xaxis=list(title="GOLD Stage"), yaxis=list(title="Time points"),
             showlegend=FALSE)
  })
  output$plot_emph <- renderPlotly({
    d <- sim_data()
    plot_ly(d, x=~time_day, y=~Emph, type="scatter", mode="lines",
            line=list(color="#FF6F00", width=2)) %>%
      layout(xaxis=list(title="Day"), yaxis=list(title="Emphysema index (% LAA)"))
  })
  output$plot_pvr <- renderPlotly({
    d <- sim_data()
    plot_ly(d, x=~time_day, y=~PVR, type="scatter", mode="lines",
            line=list(color="#880E4F", width=2)) %>%
      add_lines(y=rep(240, nrow(d)), name="PH threshold", line=list(dash="dash",color="#B71C1C")) %>%
      layout(xaxis=list(title="Day"), yaxis=list(title="PVR (dyn·s·cm⁻⁵)"))
  })

  # ── Tab 4: Biomarkers ───────────────────────────────────────────────────────
  make_biomarker_plot <- function(var, title, ylab, refline=NULL, col="#E53935") {
    d <- sim_data()
    p <- plot_ly(d, x=~time_day, y=~get(var), type="scatter", mode="lines",
                 line=list(color=col, width=2), name=ylab)
    if (!is.null(refline)) {
      p <- p %>% add_lines(y=rep(refline, nrow(d)), name="Baseline",
                           line=list(dash="dash", color="#616161"))
    }
    p %>% layout(xaxis=list(title="Day"), yaxis=list(title=ylab), showlegend=FALSE)
  }
  output$plot_il8 <- renderPlotly(make_biomarker_plot("IL8","IL-8","IL-8 (pg/mL)",125,"#E53935"))
  output$plot_ne  <- renderPlotly(make_biomarker_plot("NE_sput","NE","NE (µg/mL)",50,"#C62828"))
  output$plot_crp <- renderPlotly(make_biomarker_plot("CRP","CRP","CRP (mg/L)",5,"#F57C00"))
  output$plot_eos <- renderPlotly({
    d <- sim_data()
    plot_ly(d, x=~time_day, y=~Eos, type="scatter", mode="lines",
            line=list(color="#1565C0",width=2)) %>%
      add_lines(y=rep(300, nrow(d)), name="ICS benefit (300)",
                line=list(dash="dash",color="#1B5E20")) %>%
      add_lines(y=rep(100, nrow(d)), name="ICS harm threshold (100)",
                line=list(dash="dash",color="#B71C1C")) %>%
      layout(xaxis=list(title="Day"), yaxis=list(title="Blood Eos (cells/µL)"))
  })

  # ── Tab 5: Clinical ─────────────────────────────────────────────────────────
  output$plot_cat <- renderPlotly({
    d <- sim_data()
    plot_ly(d, x=~time_day, y=~CAT_approx, type="scatter", mode="lines",
            line=list(color="#1565C0",width=2)) %>%
      layout(xaxis=list(title="Day"), yaxis=list(title="CAT Score",range=c(0,40)))
  })
  output$plot_ae <- renderPlotly({
    d <- sim_data()
    plot_ly(d, x=~time_day, y=~AE_cum, type="scatter", mode="lines",
            line=list(color="#E53935",width=2)) %>%
      layout(xaxis=list(title="Day"), yaxis=list(title="Cumulative Exacerbations"))
  })
  output$plot_ae_rate <- renderPlotly({
    d <- sim_data()
    plot_ly(d, x=~time_day, y=~AE_rate_ann, type="scatter", mode="lines",
            line=list(color="#FF6F00",width=2)) %>%
      layout(xaxis=list(title="Day"), yaxis=list(title="Annualized AE Rate (/yr)"))
  })
  output$plot_spo2 <- renderPlotly({
    d <- sim_data()
    plot_ly(d, x=~time_day, y=~SpO2, type="scatter", mode="lines",
            line=list(color="#0277BD",width=2)) %>%
      add_lines(y=rep(88, nrow(d)), name="LTOT threshold (88%)",
                line=list(dash="dash",color="#B71C1C")) %>%
      layout(xaxis=list(title="Day"), yaxis=list(title="SpO2 (%)", range=c(80,100)))
  })
  output$plot_mpap <- renderPlotly({
    d <- sim_data()
    plot_ly(d, x=~time_day, y=~mPAP, type="scatter", mode="lines",
            line=list(color="#880E4F",width=2)) %>%
      add_lines(y=rep(25, nrow(d)), name="PH threshold (25 mmHg)",
                line=list(dash="dash",color="#E53935")) %>%
      layout(xaxis=list(title="Day"), yaxis=list(title="mPAP (mmHg)"))
  })

  # ── Tab 6: Scenario Comparison ─────────────────────────────────────────────
  scen_pal <- c(
    "1. Placebo"="grey40","2. LAMA"="#1565C0","3. LABA+LAMA"="#2E7D32",
    "4. ICS+LABA"="#AD1457","5. Triple"="#E65100","6. Triple+PDE4i"="#6A1B9A"
  )
  output$plot_compare_fev1 <- renderPlotly({
    d <- all_scenarios()
    p <- plot_ly()
    for (s in unique(d$scenario)) {
      ds <- filter(d, scenario == s)
      p <- add_lines(p, x=ds$time_day, y=ds$FEV1, name=s,
                     line=list(color=scen_pal[s], width=2))
    }
    layout(p, xaxis=list(title="Day"), yaxis=list(title="FEV1 (% predicted)"),
           legend=list(orientation="h"))
  })
  output$plot_compare_il8 <- renderPlotly({
    d <- all_scenarios()
    p <- plot_ly()
    for (s in unique(d$scenario)) {
      ds <- filter(d, scenario == s)
      p <- add_lines(p, x=ds$time_day, y=ds$IL8, name=s,
                     line=list(color=scen_pal[s], width=2))
    }
    layout(p, xaxis=list(title="Day"), yaxis=list(title="IL-8 (pg/mL)"),
           legend=list(orientation="h", font=list(size=9)))
  })
  output$plot_compare_ae <- renderPlotly({
    d <- all_scenarios()
    p <- plot_ly()
    for (s in unique(d$scenario)) {
      ds <- filter(d, scenario == s)
      p <- add_lines(p, x=ds$time_day, y=ds$AE_cum, name=s,
                     line=list(color=scen_pal[s], width=2))
    }
    layout(p, xaxis=list(title="Day"), yaxis=list(title="Cumulative AEs"),
           legend=list(orientation="h", font=list(size=9)))
  })
  output$tbl_compare <- renderDT({
    d <- all_scenarios()
    tbl <- d %>%
      group_by(scenario) %>%
      summarise(
        `FEV1 end (%)` = round(last(FEV1), 1),
        `GOLD Stage`   = round(last(GOLD_stage), 0),
        `IL-8 (pg/mL)` = round(last(IL8), 1),
        `CRP (mg/L)`   = round(last(CRP), 2),
        `Eos (cells/µL)` = round(last(Eos), 0),
        `Cum AEs`      = round(last(AE_cum), 2),
        `AE rate/yr`   = round(last(AE_rate_ann), 2),
        `Emphysema %`  = round(last(Emph), 1),
        `PVR (dyn)`    = round(last(PVR), 0),
        `CAT`          = round(last(CAT_approx), 1),
        `SpO2 (%)`     = round(last(SpO2), 1),
        .groups = "drop"
      )
    datatable(tbl, options=list(scrollX=TRUE, pageLength=6), rownames=FALSE) %>%
      formatStyle("GOLD Stage",
                  backgroundColor = styleEqual(c(1,2,3,4),
                                               c("#C8E6C9","#FFF9C4","#FFE0B2","#FFCDD2")))
  })

  # ── Tab 7: AE Risk ──────────────────────────────────────────────────────────
  ae_risk_data <- eventReactive(input$ae_calc, {
    base_rate <- input$ae_history * (1 + (60 - input$ae_fev1) / 60)
    eos_mult  <- ifelse(input$ae_eos >= 300, 1.0,
                        ifelse(input$ae_eos >= 100, 1.2, 1.5))
    base_rate <- base_rate * eos_mult

    reductions <- c(
      "None"=1.0, "LAMA"=0.85, "LABA+LAMA"=0.78,
      "ICS+LABA"=0.75, "Triple (LAMA+LABA+ICS)"=0.63,
      "Triple+Roflumilast"=0.52
    )
    sel <- input$ae_triple
    estimated_rate <- base_rate * reductions[sel]
    list(base_rate = base_rate, estimated_rate = estimated_rate, reductions = reductions)
  })

  output$plot_ae_risk <- renderPlotly({
    req(ae_risk_data())
    rd <- ae_risk_data()
    plot_ly(
      type="indicator", mode="gauge+number+delta",
      value = round(rd$estimated_rate, 2),
      delta = list(reference = round(rd$base_rate, 2), decreasing=list(color="green")),
      title = list(text="Estimated AE Rate (/yr)"),
      gauge = list(
        axis = list(range=list(0, 8)),
        bar  = list(color="#E53935"),
        steps = list(
          list(range=c(0,1), color="#C8E6C9"),
          list(range=c(1,2), color="#FFF9C4"),
          list(range=c(2,4), color="#FFE0B2"),
          list(range=c(4,8), color="#FFCDD2")
        )
      )
    )
  })

  output$plot_rr_waterfall <- renderPlotly({
    rr_data <- data.frame(
      Treatment = c("LAMA\n(UPLIFT)",
                    "LABA+LAMA\n(FLAME)",
                    "ICS+LABA\n(TORCH)",
                    "Triple\n(IMPACT)",
                    "Triple+PDE4i\n(ETHOS+EINSTEIN)"),
      RRR_pct = c(14, 22, 25, 37, 48),
      Trial = c("UPLIFT 2008","FLAME 2016","TORCH 2007","IMPACT 2018","Combined")
    )
    plot_ly(rr_data, x=~Treatment, y=~RRR_pct, type="bar",
            marker=list(color=c("#1565C0","#2E7D32","#AD1457","#E65100","#6A1B9A")),
            text=~paste0(RRR_pct,"% RRR (", Trial,")"),
            textposition="outside") %>%
      layout(
        xaxis=list(title="Treatment Strategy"),
        yaxis=list(title="Relative Risk Reduction (%) vs Placebo", range=c(0,60)),
        showlegend=FALSE,
        title="AECOPD Relative Risk Reduction by Treatment Strategy"
      )
  })
}

shinyApp(ui, server)
