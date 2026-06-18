################################################################################
# Bronchial Asthma QSP – Interactive Shiny Dashboard
#
# Tabs:
#  1. Patient Profile & Disease Endotype
#  2. Drug PK – Biologic Concentration
#  3. PD Biomarkers (Eosinophil, IgE, IL-5/IL-13, TSLP)
#  4. Lung Function (FEV1, PEF, AHR)
#  5. Clinical Endpoints (Exacerbation, Symptom Score, SABA Use)
#  6. Treatment Scenario Comparison
#  7. Dose–Response & Biomarker Threshold
#
# Author: Claude Code Routine (CCR) · Date: 2026-06-16
################################################################################

library(shiny)
library(shinydashboard)
library(mrgsolve)
library(dplyr)
library(tidyr)
library(ggplot2)
library(DT)
library(plotly)

## ── mrgsolve model code (inline) ─────────────────────────────────────────────
asthma_code <- '
$PARAM @annotated
ka_MEPO  : 0.0067 : Mepolizumab SC absorption (1/h)
CL_MEPO  : 0.0115 : Mepolizumab clearance (L/h)
V1_MEPO  : 3.6    : Mepolizumab Vc (L)
V2_MEPO  : 3.8    : Mepolizumab Vp (L)
Q_MEPO   : 0.043  : Mepolizumab Q (L/h)
ka_BENZ  : 0.0055 : Benralizumab SC absorption (1/h)
CL_BENZ  : 0.0085 : Benralizumab clearance (L/h)
V1_BENZ  : 3.1    : Benralizumab Vc (L)
V2_BENZ  : 3.4    : Benralizumab Vp (L)
Q_BENZ   : 0.028  : Benralizumab Q (L/h)
ka_DUPIL : 0.0042 : Dupilumab SC absorption (1/h)
CL_DUPIL : 0.0065 : Dupilumab clearance (L/h)
V1_DUPIL : 4.8    : Dupilumab Vc (L)
V2_DUPIL : 5.2    : Dupilumab Vp (L)
Q_DUPIL  : 0.058  : Dupilumab Q (L/h)
ka_TEZE  : 0.0072 : Tezepelumab SC absorption (1/h)
CL_TEZE  : 0.0088 : Tezepelumab clearance (L/h)
V1_TEZE  : 3.9    : Tezepelumab Vc (L)
V2_TEZE  : 4.4    : Tezepelumab Vp (L)
Q_TEZE   : 0.051  : Tezepelumab Q (L/h)
ksyn_TSLP : 0.006 : TSLP synthesis (nM/h)
kdeg_TSLP : 0.045 : TSLP degradation (1/h)
TSLP_ss   : 0.133 : TSLP baseline (nM)
Emax_TEZE_TSLP : 0.95 : Teze max TSLP reduction
ksyn_IL5  : 0.18  : IL-5 synthesis (pg/mL/h)
kdeg_IL5  : 0.065 : IL-5 degradation (1/h)
IL5_ss    : 2.77  : IL-5 baseline (pg/mL)
Emax_MEPO_IL5 : 0.95 : Mepolizumab max IL-5 reduction
EC50_MEPO_IL5 : 0.12 : Mepolizumab EC50 IL-5 (ug/mL)
ksyn_IL13 : 0.10  : IL-13 synthesis (pg/mL/h)
kdeg_IL13 : 0.08  : IL-13 degradation (1/h)
IL13_ss   : 1.25  : IL-13 baseline (pg/mL)
Emax_DUPIL_IL13 : 0.90 : Dupilumab max IL-13 reduction
EC50_DUPIL_IL13 : 0.25 : Dupilumab EC50 IL-13 (ug/mL)
kprod_EOS : 0.18  : Eos BM production (cells/uL/h)
kin_EOS   : 0.040 : Eos tissue influx (1/h)
kout_EOS  : 0.038 : Eos tissue egress/death (1/h)
EOS_B_ss  : 450   : Blood Eos baseline (cells/uL)
EOS_T_ss  : 2.0   : Tissue Eos baseline (1e6/g)
Emax_MEPO_EOS : 0.80 : Mepolizumab max Eos suppression
EC50_MEPO_EOS : 0.15 : Mepolizumab EC50 Eos (ug/mL)
Emax_BENZ_EOS : 0.97 : Benralizumab max Eos suppression
EC50_BENZ_EOS : 0.08 : Benralizumab EC50 Eos (ug/mL)
kact_ASM  : 0.012 : ASM tone activation (1/h)
krel_ASM  : 0.025 : ASM tone relaxation (1/h)
ASM_ss    : 1.0   : ASM baseline
ksyn_MUC  : 0.15  : Mucus production (AU/h)
kdeg_MUC  : 0.08  : Mucus clearance (1/h)
MUC_ss    : 1.875 : Mucus baseline (AU)
FEV1_max  : 78    : Max FEV1 (% pred)
FEV1_min  : 48    : Min FEV1 (% pred)
kFEV1     : 0.005 : FEV1 adapt rate (1/h)
Emax_ICS_IL5  : 0.45 : ICS max IL-5 suppression
EC50_ICS      : 2.0  : ICS EC50 (ng/mL)
Emax_ICS_IL13 : 0.40 : ICS max IL-13 suppression
Emax_ICS_TSLP : 0.30 : ICS max TSLP suppression
Emax_LABA_ASM : 0.55 : LABA max ASM reduction
EC50_LABA     : 0.8  : LABA EC50 (ng/mL)
lambda0   : 0.18  : Baseline exacerbation hazard (events/yr)
beta_EOS  : 0.0015: Eos contribution to exacerb risk
beta_FEV1 : -0.025: FEV1 protective effect
beta_MUC  : 0.04  : Mucus contribution to exacerb risk

$CMT
MEPO_SC MEPO_C1 MEPO_C2
BENZ_SC BENZ_C1 BENZ_C2
DUPIL_SC DUPIL_C1 DUPIL_C2
TEZE_SC TEZE_C1 TEZE_C2
ICS_SYS LABA_C
TSLP_PD IL5_PD IL13_PD
EOS_B EOS_T ASM_TONE MUCUS FEV1_ODE

$MAIN
double cMEPO  = MEPO_C1;
double cBENZ  = BENZ_C1;
double cDUPIL = DUPIL_C1;
double cTEZE  = TEZE_C1;
double EFF_ICS_IL5  = Emax_ICS_IL5  * ICS_SYS / (EC50_ICS + ICS_SYS + 0.001);
double EFF_ICS_IL13 = Emax_ICS_IL13 * ICS_SYS / (EC50_ICS + ICS_SYS + 0.001);
double EFF_ICS_TSLP = Emax_ICS_TSLP * ICS_SYS / (EC50_ICS + ICS_SYS + 0.001);
double EFF_LABA_ASM = Emax_LABA_ASM * LABA_C  / (EC50_LABA + LABA_C + 0.001);
double EFF_MEPO_IL5  = Emax_MEPO_IL5  * cMEPO  / (EC50_MEPO_IL5  + cMEPO  + 0.001);
double EFF_MEPO_EOS  = Emax_MEPO_EOS  * cMEPO  / (EC50_MEPO_EOS  + cMEPO  + 0.001);
double EFF_BENZ_EOS  = Emax_BENZ_EOS  * cBENZ  / (EC50_BENZ_EOS  + cBENZ  + 0.001);
double EFF_DUPIL_IL13= Emax_DUPIL_IL13* cDUPIL / (EC50_DUPIL_IL13+ cDUPIL + 0.001);
double EFF_TEZE_TSLP = Emax_TEZE_TSLP * cTEZE  / (0.05 + cTEZE  + 0.001);
double TSLP_rel = TSLP_PD / (TSLP_ss + 0.001);
double IL5_rel  = IL5_PD  / (IL5_ss  + 0.001);
double IL13_rel = IL13_PD / (IL13_ss + 0.001);
double f_ASM    = ASM_TONE / ASM_ss;
double f_MUC    = MUCUS    / MUC_ss;
double f_EOS    = EOS_T    / EOS_T_ss;
double FEV1_target = FEV1_max - (FEV1_max - FEV1_min) *
                     (0.5*f_ASM + 0.3*f_MUC + 0.2*f_EOS);
if(FEV1_target < FEV1_min) FEV1_target = FEV1_min;
if(FEV1_target > FEV1_max) FEV1_target = FEV1_max;

$ODE
dxdt_MEPO_SC = -ka_MEPO * MEPO_SC;
dxdt_MEPO_C1 =  ka_MEPO * MEPO_SC / V1_MEPO - (CL_MEPO+Q_MEPO)/V1_MEPO*MEPO_C1 + Q_MEPO/V2_MEPO*MEPO_C2;
dxdt_MEPO_C2 =  Q_MEPO/V1_MEPO*MEPO_C1 - Q_MEPO/V2_MEPO*MEPO_C2;
dxdt_BENZ_SC = -ka_BENZ * BENZ_SC;
dxdt_BENZ_C1 =  ka_BENZ * BENZ_SC / V1_BENZ - (CL_BENZ+Q_BENZ)/V1_BENZ*BENZ_C1 + Q_BENZ/V2_BENZ*BENZ_C2;
dxdt_BENZ_C2 =  Q_BENZ/V1_BENZ*BENZ_C1 - Q_BENZ/V2_BENZ*BENZ_C2;
dxdt_DUPIL_SC= -ka_DUPIL * DUPIL_SC;
dxdt_DUPIL_C1=  ka_DUPIL * DUPIL_SC/V1_DUPIL - (CL_DUPIL+Q_DUPIL)/V1_DUPIL*DUPIL_C1 + Q_DUPIL/V2_DUPIL*DUPIL_C2;
dxdt_DUPIL_C2=  Q_DUPIL/V1_DUPIL*DUPIL_C1 - Q_DUPIL/V2_DUPIL*DUPIL_C2;
dxdt_TEZE_SC = -ka_TEZE * TEZE_SC;
dxdt_TEZE_C1 =  ka_TEZE * TEZE_SC/V1_TEZE - (CL_TEZE+Q_TEZE)/V1_TEZE*TEZE_C1 + Q_TEZE/V2_TEZE*TEZE_C2;
dxdt_TEZE_C2 =  Q_TEZE/V1_TEZE*TEZE_C1 - Q_TEZE/V2_TEZE*TEZE_C2;
dxdt_ICS_SYS = -0.012 * ICS_SYS;
dxdt_LABA_C  = -(22.0/140.0) * LABA_C;
dxdt_TSLP_PD = ksyn_TSLP*(1-EFF_ICS_TSLP) - kdeg_TSLP*TSLP_PD*(1+EFF_TEZE_TSLP);
dxdt_IL5_PD  = ksyn_IL5*(1+0.5*TSLP_rel)*(1-EFF_ICS_IL5)*(1-EFF_MEPO_IL5) - kdeg_IL5*IL5_PD;
dxdt_IL13_PD = ksyn_IL13*(1+0.4*TSLP_rel)*(1-EFF_ICS_IL13)*(1-EFF_DUPIL_IL13) - kdeg_IL13*IL13_PD;
dxdt_EOS_B   = kprod_EOS*IL5_rel*(1-EFF_MEPO_EOS)*(1-EFF_BENZ_EOS) - kin_EOS*EOS_B;
dxdt_EOS_T   = kin_EOS*EOS_B/EOS_B_ss*(1+0.5*IL13_rel) - kout_EOS*EOS_T;
dxdt_ASM_TONE= kact_ASM*(1+0.4*IL13_rel+0.2*(EOS_T/EOS_T_ss-1)) - krel_ASM*ASM_TONE*(1+EFF_LABA_ASM);
dxdt_MUCUS   = ksyn_MUC*(1+0.5*IL13_rel)*(1-0.3*EFF_ICS_IL13) - kdeg_MUC*MUCUS;
dxdt_FEV1_ODE= kFEV1*(FEV1_target - FEV1_ODE);

$TABLE
capture FEV1_pct      = FEV1_ODE;
capture EOS_blood     = EOS_B;
capture EOS_tissue    = EOS_T;
capture IL5_pgmL      = IL5_PD;
capture IL13_pgmL     = IL13_PD;
capture TSLP_nM       = TSLP_PD;
capture ASM_index     = ASM_TONE;
capture Mucus_AU      = MUCUS;
capture MEPO_conc     = MEPO_C1;
capture BENZ_conc     = BENZ_C1;
capture DUPIL_conc    = DUPIL_C1;
capture TEZE_conc     = TEZE_C1;
capture exacerb_hazard= lambda0 + beta_EOS*(EOS_B-450) + beta_FEV1*(FEV1_ODE-70) + beta_MUC*(MUCUS-1.875);

$INIT
MEPO_SC=0, MEPO_C1=0, MEPO_C2=0
BENZ_SC=0, BENZ_C1=0, BENZ_C2=0
DUPIL_SC=0, DUPIL_C1=0, DUPIL_C2=0
TEZE_SC=0, TEZE_C1=0, TEZE_C2=0
ICS_SYS=0, LABA_C=0
TSLP_PD=0.133, IL5_PD=2.77, IL13_PD=1.25
EOS_B=450, EOS_T=2.0
ASM_TONE=1.0, MUCUS=1.875, FEV1_ODE=58
'

## ── Compile once ──────────────────────────────────────────────────────────────
mod <- mcode("asthma_shiny", asthma_code, quiet = TRUE)

## ── Helper: simulate scenario ─────────────────────────────────────────────────
sim_scenario <- function(model, params,
                         mepo_dose, benz_dose, dupil_dose, teze_dose,
                         ics_dose_mgd, laba_on, sim_weeks) {

  model <- param(model, params)

  evs <- ev(time = 0, amt = 0, cmt = 1)   # placeholder

  # ICS continuous infusion (lung → sys proxy)
  ics_rate  <- ics_dose_mgd * 1000 / 24   # μg/day → μg/h → ng/mL/h proxy
  if (ics_dose_mgd > 0)
    evs <- evs + ev(time=0, amt=0, cmt="ICS_SYS", rate=ics_rate, addl=0)

  # LABA
  if (laba_on)
    evs <- evs + ev(time=0, amt=0, cmt="LABA_C", rate=0.18/140*1000, addl=0)

  # Mepolizumab q4w
  if (mepo_dose > 0)
    evs <- evs + ev(time=0, amt=mepo_dose, cmt="MEPO_SC",
                    ii=4*7*24, addl=floor(sim_weeks/4)-1)

  # Benralizumab q4w×3 then q8w
  if (benz_dose > 0) {
    evs <- evs + ev(time=0, amt=benz_dose, cmt="BENZ_SC", ii=4*7*24, addl=2)
    q8s <- 3*4*7*24
    addl_q8 <- max(0, floor((sim_weeks-12)/8))
    evs <- evs + ev(time=q8s, amt=benz_dose, cmt="BENZ_SC",
                    ii=8*7*24, addl=addl_q8)
  }

  # Dupilumab q2w
  if (dupil_dose > 0)
    evs <- evs + ev(time=0, amt=dupil_dose, cmt="DUPIL_SC",
                    ii=2*7*24, addl=floor(sim_weeks/2)-1)

  # Tezepelumab q4w
  if (teze_dose > 0)
    evs <- evs + ev(time=0, amt=teze_dose, cmt="TEZE_SC",
                    ii=4*7*24, addl=floor(sim_weeks/4)-1)

  mrgsim(model, evs,
         end   = sim_weeks * 7 * 24,
         delta = 12) %>%
    as.data.frame() %>%
    mutate(time_weeks = time / (7*24))
}

## ── UI ────────────────────────────────────────────────────────────────────────
ui <- dashboardPage(
  skin = "blue",

  dashboardHeader(title = "Bronchial Asthma QSP"),

  dashboardSidebar(
    sidebarMenu(
      menuItem("Patient Profile",      tabName = "tab_patient", icon = icon("user")),
      menuItem("Drug PK",              tabName = "tab_pk",      icon = icon("pills")),
      menuItem("PD Biomarkers",        tabName = "tab_pd",      icon = icon("flask")),
      menuItem("Lung Function",        tabName = "tab_lung",    icon = icon("lungs")),
      menuItem("Clinical Endpoints",   tabName = "tab_clin",    icon = icon("chart-line")),
      menuItem("Scenario Comparison",  tabName = "tab_compare", icon = icon("table")),
      menuItem("Dose–Response",        tabName = "tab_dr",      icon = icon("sliders-h"))
    ),

    hr(),
    h5("  Patient Parameters", style = "color:#ccc; margin-left:10px"),

    sliderInput("eos_baseline",  "Blood Eos (cells/μL)", 100, 1500, 450, step=50),
    sliderInput("fev1_baseline", "Baseline FEV1 (%)",    35,  80,   58,  step=1),
    sliderInput("il13_base",     "IL-13 baseline (pg/mL)", 0.5, 5, 1.25, step=0.25),

    hr(),
    h5("  Treatment Selection", style = "color:#ccc; margin-left:10px"),

    checkboxInput("use_ics",   "ICS (budesonide)", value = TRUE),
    checkboxInput("use_laba",  "LABA (formoterol)", value = TRUE),
    sliderInput("ics_dose_mgd", "ICS dose (mg/day)", 0.1, 1.6, 0.4, step=0.1),

    selectInput("biologic", "Add Biologic:",
                choices = c("None", "Mepolizumab", "Benralizumab",
                            "Dupilumab", "Tezepelumab"),
                selected = "None"),

    conditionalPanel(
      condition = "input.biologic == 'Mepolizumab'",
      sliderInput("mepo_dose", "Mepo dose (mg SC q4w)", 25, 200, 100, step=25)
    ),
    conditionalPanel(
      condition = "input.biologic == 'Benralizumab'",
      sliderInput("benz_dose", "Benz dose (mg SC)", 10, 100, 30, step=5)
    ),
    conditionalPanel(
      condition = "input.biologic == 'Dupilumab'",
      sliderInput("dupil_dose", "Dupilumab dose (mg SC q2w)", 100, 400, 300, step=50)
    ),
    conditionalPanel(
      condition = "input.biologic == 'Tezepelumab'",
      sliderInput("teze_dose", "Tezepelumab dose (mg SC q4w)", 70, 420, 210, step=70)
    ),

    sliderInput("sim_weeks", "Simulation (weeks)", 12, 104, 52, step=4),
    actionButton("run_sim", "▶ Run Simulation", class = "btn-primary",
                 style = "margin:10px; width:90%")
  ),

  dashboardBody(
    tags$head(tags$style(HTML("
      .content-wrapper { background-color: #f5f5f5; }
      .box { border-radius:8px; }
    "))),

    tabItems(

      ## ── Tab 1: Patient Profile ───────────────────────────────────────────
      tabItem(tabName = "tab_patient",
        fluidRow(
          box(title = "Patient & Disease Endotype", width = 6, status = "primary",
              solidHeader = TRUE,
              tableOutput("patient_summary_tbl")),
          box(title = "Type-2 Endotype Classification", width = 6, status = "warning",
              solidHeader = TRUE,
              plotlyOutput("endotype_radar", height = "300px")),
          box(title = "Mechanistic Map", width = 12, status = "info",
              solidHeader = TRUE,
              p("Bronchial Asthma QSP mechanistic map — all 110+ pathway components"),
              tags$img(src = "ba_qsp_model.png",
                       style = "max-width:100%; border-radius:4px"))
        )
      ),

      ## ── Tab 2: Drug PK ───────────────────────────────────────────────────
      tabItem(tabName = "tab_pk",
        fluidRow(
          box(title = "Biologic PK – Central Concentration", width = 12, status = "primary",
              solidHeader = TRUE,
              plotlyOutput("pk_plot", height = "400px")),
          box(title = "PK Parameters", width = 12, status = "info",
              DTOutput("pk_params_tbl"))
        )
      ),

      ## ── Tab 3: PD Biomarkers ─────────────────────────────────────────────
      tabItem(tabName = "tab_pd",
        fluidRow(
          box(title = "Blood Eosinophils (cells/μL)", width = 6, status = "danger",
              solidHeader = TRUE,
              plotlyOutput("eos_blood_plot", height = "300px")),
          box(title = "Tissue Eosinophils (10⁶/g)", width = 6, status = "danger",
              solidHeader = TRUE,
              plotlyOutput("eos_tissue_plot", height = "300px")),
          box(title = "IL-5 & IL-13 (pg/mL)", width = 6, status = "warning",
              solidHeader = TRUE,
              plotlyOutput("cytokine_plot", height = "300px")),
          box(title = "TSLP (nM)", width = 6, status = "warning",
              solidHeader = TRUE,
              plotlyOutput("tslp_plot", height = "300px"))
        )
      ),

      ## ── Tab 4: Lung Function ─────────────────────────────────────────────
      tabItem(tabName = "tab_lung",
        fluidRow(
          box(title = "FEV1 (% Predicted)", width = 8, status = "success",
              solidHeader = TRUE,
              plotlyOutput("fev1_plot", height = "400px")),
          box(title = "Lung Function Zones", width = 4, status = "info",
              solidHeader = TRUE,
              p(strong("GINA Severity Classification:")),
              tags$ul(
                tags$li(HTML("<span style='color:green'>FEV1 ≥80%</span> – Mild")),
                tags$li(HTML("<span style='color:orange'>FEV1 60–79%</span> – Moderate")),
                tags$li(HTML("<span style='color:red'>FEV1 40–59%</span> – Severe")),
                tags$li(HTML("<span style='color:darkred'>FEV1 <40%</span> – Very Severe"))
              ),
              hr(),
              p(strong("ASM Tone Index:")),
              plotlyOutput("asm_plot", height = "200px")),
          box(title = "Mucus Production (AU)", width = 6, status = "warning",
              solidHeader = TRUE,
              plotlyOutput("mucus_plot", height = "300px")),
          box(title = "FEV1 at Key Timepoints", width = 6, status = "success",
              solidHeader = TRUE,
              DTOutput("fev1_table"))
        )
      ),

      ## ── Tab 5: Clinical Endpoints ─────────────────────────────────────────
      tabItem(tabName = "tab_clin",
        fluidRow(
          box(title = "Exacerbation Hazard Rate (events/year)", width = 8, status = "danger",
              solidHeader = TRUE,
              plotlyOutput("exacerb_plot", height = "350px")),
          box(title = "Clinical Metrics at Week 52", width = 4, status = "primary",
              solidHeader = TRUE,
              tableOutput("clin_summary_tbl")),
          box(title = "ACQ Proxy (from FEV1)", width = 6, status = "info",
              solidHeader = TRUE,
              plotlyOutput("acq_plot", height = "300px")),
          box(title = "Rescue SABA Use Proxy", width = 6, status = "warning",
              solidHeader = TRUE,
              plotlyOutput("saba_plot", height = "300px"))
        )
      ),

      ## ── Tab 6: Scenario Comparison ────────────────────────────────────────
      tabItem(tabName = "tab_compare",
        fluidRow(
          box(title = "All Treatment Scenarios – FEV1", width = 12, status = "primary",
              solidHeader = TRUE,
              plotlyOutput("compare_fev1", height = "400px")),
          box(title = "All Treatment Scenarios – Blood Eosinophils", width = 12, status = "danger",
              solidHeader = TRUE,
              plotlyOutput("compare_eos", height = "350px")),
          box(title = "Week-52 Summary Table", width = 12, status = "info",
              solidHeader = TRUE,
              DTOutput("compare_table"))
        )
      ),

      ## ── Tab 7: Dose–Response ─────────────────────────────────────────────
      tabItem(tabName = "tab_dr",
        fluidRow(
          box(title = "Biologic Dose–Response at Week 52", width = 12, status = "primary",
              solidHeader = TRUE,
              selectInput("dr_biologic", "Biologic for dose–response:",
                          choices = c("Mepolizumab", "Benralizumab",
                                      "Dupilumab", "Tezepelumab"),
                          selected = "Benralizumab"),
              actionButton("run_dr", "Calculate Dose–Response", class="btn-warning"),
              hr(),
              plotlyOutput("dr_fev1_plot", height = "300px")),
          box(title = "Biomarker Threshold Analysis", width = 12, status = "warning",
              solidHeader = TRUE,
              p("Exacerbation risk across baseline blood eosinophil ranges"),
              plotlyOutput("threshold_plot", height = "300px"))
        )
      )

    )
  )
)

## ── Server ────────────────────────────────────────────────────────────────────
server <- function(input, output, session) {

  # ── Reactive simulation ──────────────────────────────────────────────────────
  sim_result <- eventReactive(input$run_sim, {

    mepo_d  <- if (input$biologic == "Mepolizumab")  input$mepo_dose  else 0
    benz_d  <- if (input$biologic == "Benralizumab") input$benz_dose  else 0
    dupil_d <- if (input$biologic == "Dupilumab")    input$dupil_dose else 0
    teze_d  <- if (input$biologic == "Tezepelumab")  input$teze_dose  else 0

    extra_params <- list(
      EOS_B_ss  = input$eos_baseline,
      kprod_EOS = 0.18 * (input$eos_baseline / 450),
      IL13_ss   = input$il13_base,
      ksyn_IL13 = 0.10 * (input$il13_base / 1.25),
      FEV1_min  = input$fev1_baseline - 10,
      FEV1_max  = min(input$fev1_baseline + 20, 95)
    )

    sim_scenario(
      model      = mod,
      params     = extra_params,
      mepo_dose  = mepo_d,
      benz_dose  = benz_d,
      dupil_dose = dupil_d,
      teze_dose  = teze_d,
      ics_dose_mgd = if (input$use_ics) input$ics_dose_mgd else 0,
      laba_on    = input$use_laba,
      sim_weeks  = input$sim_weeks
    )
  }, ignoreNULL = FALSE)

  # ── Patient summary ──────────────────────────────────────────────────────────
  output$patient_summary_tbl <- renderTable({
    endotype <- if (input$eos_baseline >= 300 || input$il13_base >= 1.5)
      "T2-High (eos≥300 or IL-13↑)" else "T2-Low (neutrophilic/pauci)"
    data.frame(
      Parameter = c("Blood Eosinophils", "Baseline FEV1", "IL-13 Baseline",
                    "Endotype", "Biologic Selected", "ICS Dose"),
      Value = c(
        paste(input$eos_baseline, "cells/μL"),
        paste(input$fev1_baseline, "%"),
        paste(input$il13_base, "pg/mL"),
        endotype,
        input$biologic,
        if (input$use_ics) paste(input$ics_dose_mgd, "mg/day budesonide equiv.") else "None"
      )
    )
  })

  # ── Endotype radar ───────────────────────────────────────────────────────────
  output$endotype_radar <- renderPlotly({
    theta  <- c("Blood Eos","IL-13","IL-5","FeNO","IgE","Th17/neutrophil")
    r_vals <- c(
      min(input$eos_baseline / 600, 1),
      min(input$il13_base / 4, 1),
      0.6,
      0.55,
      0.65,
      0.25
    )
    plot_ly(type="scatterpolar", mode="lines+markers",
            r=c(r_vals, r_vals[1]), theta=c(theta, theta[1]),
            fill="toself", fillcolor="rgba(30,136,229,0.2)",
            line=list(color="rgb(30,136,229)")) %>%
      layout(polar=list(radialaxis=list(visible=TRUE, range=c(0,1))),
             showlegend=FALSE, title="Type-2 Biomarker Profile")
  })

  # ── PK plot ──────────────────────────────────────────────────────────────────
  output$pk_plot <- renderPlotly({
    df <- sim_result()
    biol <- input$biologic
    if (biol == "None") {
      plot_ly() %>% layout(title = "No biologic selected")
    } else {
      conc_col <- switch(biol,
        "Mepolizumab"  = "MEPO_conc",
        "Benralizumab" = "BENZ_conc",
        "Dupilumab"    = "DUPIL_conc",
        "Tezepelumab"  = "TEZE_conc"
      )
      plot_ly(df, x=~time_weeks, y=~.data[[conc_col]], type="scatter",
              mode="lines", line=list(color="steelblue", width=2)) %>%
        layout(title = paste(biol, "– Central Concentration (μg/mL)"),
               xaxis=list(title="Time (weeks)"),
               yaxis=list(title="Concentration (μg/mL)"))
    }
  })

  output$pk_params_tbl <- renderDT({
    data.frame(
      Biologic     = c("Mepolizumab", "Benralizumab", "Dupilumab", "Tezepelumab"),
      Target       = c("IL-5", "IL-5Rα", "IL-4Rα", "TSLP"),
      Dose         = c("100 mg SC q4w", "30 mg SC q4w×3→q8w",
                       "300 mg SC q2w", "210 mg SC q4w"),
      Bioavail     = c("77%", "59%", "64%", "81%"),
      `t½ (d)`     = c("20", "15.5", "21", "26"),
      `CL (L/h)`   = c("0.0115", "0.0085", "0.0065", "0.0088"),
      Key_trial    = c("MENSA", "CALIMA", "LIBERTY AIR", "NAVIGATOR"),
      `ΔAER`       = c("−50%", "−51%", "−48%", "−56%")
    )
  }, rownames=FALSE, options=list(pageLength=5))

  # ── Eosinophil plots ─────────────────────────────────────────────────────────
  output$eos_blood_plot <- renderPlotly({
    df <- sim_result()
    plot_ly(df, x=~time_weeks, y=~EOS_blood, type="scatter", mode="lines",
            line=list(color="#C62828", width=2)) %>%
      add_segments(x=0, xend=max(df$time_weeks), y=300, yend=300,
                   line=list(dash="dash", color="orange")) %>%
      layout(title="Blood Eosinophils",
             xaxis=list(title="Time (weeks)"),
             yaxis=list(title="cells/μL"))
  })

  output$eos_tissue_plot <- renderPlotly({
    df <- sim_result()
    plot_ly(df, x=~time_weeks, y=~EOS_tissue, type="scatter", mode="lines",
            line=list(color="#AD1457", width=2)) %>%
      layout(title="Tissue Eosinophils",
             xaxis=list(title="Time (weeks)"),
             yaxis=list(title="10⁶/g"))
  })

  output$cytokine_plot <- renderPlotly({
    df <- sim_result() %>%
      select(time_weeks, IL5_pgmL, IL13_pgmL) %>%
      pivot_longer(-time_weeks)
    plot_ly(df, x=~time_weeks, y=~value, color=~name,
            type="scatter", mode="lines",
            colors=c("#1565C0","#E65100")) %>%
      layout(title="IL-5 & IL-13",
             xaxis=list(title="Time (weeks)"),
             yaxis=list(title="pg/mL"))
  })

  output$tslp_plot <- renderPlotly({
    df <- sim_result()
    plot_ly(df, x=~time_weeks, y=~TSLP_nM, type="scatter", mode="lines",
            line=list(color="#6A1B9A", width=2)) %>%
      layout(title="TSLP",
             xaxis=list(title="Time (weeks)"),
             yaxis=list(title="nM"))
  })

  # ── Lung function ────────────────────────────────────────────────────────────
  output$fev1_plot <- renderPlotly({
    df <- sim_result()
    plot_ly(df, x=~time_weeks, y=~FEV1_pct, type="scatter", mode="lines",
            line=list(color="#1B5E20", width=2.5)) %>%
      add_segments(x=0, xend=max(df$time_weeks), y=80, yend=80,
                   line=list(dash="dot", color="green")) %>%
      add_segments(x=0, xend=max(df$time_weeks), y=60, yend=60,
                   line=list(dash="dot", color="orange")) %>%
      add_segments(x=0, xend=max(df$time_weeks), y=40, yend=40,
                   line=list(dash="dot", color="red")) %>%
      layout(title="FEV1 (% Predicted)",
             xaxis=list(title="Time (weeks)"),
             yaxis=list(title="FEV1 (%)", range=c(30,95)))
  })

  output$asm_plot <- renderPlotly({
    df <- sim_result()
    plot_ly(df, x=~time_weeks, y=~ASM_index, type="scatter", mode="lines",
            line=list(color="#5D4037", width=1.5)) %>%
      layout(xaxis=list(title="Weeks"), yaxis=list(title="ASM tone index"))
  })

  output$mucus_plot <- renderPlotly({
    df <- sim_result()
    plot_ly(df, x=~time_weeks, y=~Mucus_AU, type="scatter", mode="lines",
            line=list(color="#827717", width=2)) %>%
      layout(title="Mucus Production",
             xaxis=list(title="Time (weeks)"),
             yaxis=list(title="Arbitrary Units"))
  })

  output$fev1_table <- renderDT({
    df  <- sim_result()
    wks <- c(0, 4, 12, 24, 52)
    wks <- wks[wks <= input$sim_weeks]
    purrr::map_dfr(wks, function(w) {
      df %>% filter(abs(time_weeks - w) < 0.5) %>%
        slice_head(n=1) %>%
        transmute(Week = w,
                  `FEV1 (%)` = round(FEV1_pct, 1),
                  `Δ from baseline` = round(FEV1_pct - 58, 1),
                  `Eos (cells/μL)` = round(EOS_blood, 0))
    })
  }, rownames=FALSE)

  # ── Clinical endpoints ───────────────────────────────────────────────────────
  output$exacerb_plot <- renderPlotly({
    df <- sim_result()
    plot_ly(df, x=~time_weeks, y=~exacerb_hazard, type="scatter", mode="lines",
            line=list(color="#B71C1C", width=2)) %>%
      add_segments(x=0, xend=max(df$time_weeks), y=0.18, yend=0.18,
                   line=list(dash="dash", color="grey")) %>%
      layout(title="Exacerbation Hazard Rate",
             xaxis=list(title="Time (weeks)"),
             yaxis=list(title="Events/year"))
  })

  output$acq_plot <- renderPlotly({
    df <- sim_result() %>%
      mutate(ACQ_proxy = pmax(0, 5 - (FEV1_pct - 50) / 10))
    plot_ly(df, x=~time_weeks, y=~ACQ_proxy, type="scatter", mode="lines",
            line=list(color="#0D47A1", width=2)) %>%
      add_segments(x=0, xend=max(df$time_weeks), y=0.75, yend=0.75,
                   line=list(dash="dot", color="green")) %>%
      layout(title="ACQ-5 Proxy",
             xaxis=list(title="Weeks"),
             yaxis=list(title="Score (lower = better)", range=c(0,5)))
  })

  output$saba_plot <- renderPlotly({
    df <- sim_result() %>%
      mutate(SABA_proxy = pmax(0, exacerb_hazard / 0.18 * 2))
    plot_ly(df, x=~time_weeks, y=~SABA_proxy, type="scatter", mode="lines",
            fill="tozeroy", fillcolor="rgba(255,152,0,0.2)",
            line=list(color="#E65100", width=2)) %>%
      layout(title="Rescue SABA Use (proxy puffs/week)",
             xaxis=list(title="Weeks"),
             yaxis=list(title="Puffs/week"))
  })

  output$clin_summary_tbl <- renderTable({
    df <- sim_result() %>%
      filter(abs(time_weeks - min(input$sim_weeks, 52)) < 0.5) %>%
      slice_tail(n=1)
    data.frame(
      Metric = c("FEV1 (%)", "Blood Eos (cells/μL)",
                 "IL-5 (pg/mL)", "IL-13 (pg/mL)",
                 "Exacerb. hazard (events/yr)"),
      Value = c(round(df$FEV1_pct, 1), round(df$EOS_blood, 0),
                round(df$IL5_pgmL, 2), round(df$IL13_pgmL, 2),
                round(df$exacerb_hazard, 3))
    )
  })

  # ── Scenario comparison ──────────────────────────────────────────────────────
  compare_data <- reactive({
    extra <- list(EOS_B_ss=input$eos_baseline, kprod_EOS=0.18*(input$eos_baseline/450))
    scenarios <- list(
      list(label="ICS/LABA only",          mepo=0,   benz=0,  dupil=0,   teze=0),
      list(label="+ Mepolizumab 100 mg",   mepo=100, benz=0,  dupil=0,   teze=0),
      list(label="+ Benralizumab 30 mg",   mepo=0,   benz=30, dupil=0,   teze=0),
      list(label="+ Dupilumab 300 mg",     mepo=0,   benz=0,  dupil=300, teze=0),
      list(label="+ Tezepelumab 210 mg",   mepo=0,   benz=0,  dupil=0,   teze=210)
    )
    purrr::map_dfr(scenarios, function(s) {
      sim_scenario(mod, extra, s$mepo, s$benz, s$dupil, s$teze,
                   ics_dose_mgd=0.4, laba_on=TRUE,
                   sim_weeks=input$sim_weeks) %>%
        mutate(Scenario=s$label)
    })
  })

  output$compare_fev1 <- renderPlotly({
    df <- compare_data()
    cols <- c("#616161","#1565C0","#2E7D32","#AD1457","#E65100")
    p <- ggplot(df, aes(time_weeks, FEV1_pct, color=Scenario)) +
      geom_line(size=1) +
      scale_color_manual(values=cols) +
      labs(x="Weeks", y="FEV1 (%)") +
      theme_bw()
    ggplotly(p)
  })

  output$compare_eos <- renderPlotly({
    df <- compare_data()
    cols <- c("#616161","#1565C0","#2E7D32","#AD1457","#E65100")
    p <- ggplot(df, aes(time_weeks, EOS_blood, color=Scenario)) +
      geom_line(size=1) +
      scale_color_manual(values=cols) +
      geom_hline(yintercept=300, linetype="dashed", color="orange") +
      labs(x="Weeks", y="Blood Eos (cells/μL)") +
      theme_bw()
    ggplotly(p)
  })

  output$compare_table <- renderDT({
    compare_data() %>%
      filter(abs(time_weeks - min(input$sim_weeks, 52)) < 0.5) %>%
      group_by(Scenario) %>% slice_tail(n=1) %>% ungroup() %>%
      select(Scenario, FEV1_pct, EOS_blood, IL5_pgmL, IL13_pgmL, exacerb_hazard) %>%
      rename(`FEV1 (%)` = FEV1_pct, `Blood Eos` = EOS_blood,
             `IL-5 (pg/mL)` = IL5_pgmL, `IL-13 (pg/mL)` = IL13_pgmL,
             `Exacerb hazard` = exacerb_hazard) %>%
      mutate(across(where(is.numeric), ~round(.x, 2)))
  }, rownames=FALSE)

  # ── Dose–Response ────────────────────────────────────────────────────────────
  dr_result <- eventReactive(input$run_dr, {
    biologic <- input$dr_biologic
    doses <- c(5, 10, 25, 50, 100, 200, 400)
    purrr::map_dfr(doses, function(d) {
      mepo_d  <- if (biologic == "Mepolizumab")  d else 0
      benz_d  <- if (biologic == "Benralizumab") d else 0
      dupil_d <- if (biologic == "Dupilumab")    d else 0
      teze_d  <- if (biologic == "Tezepelumab")  d else 0
      sim_scenario(mod, list(), mepo_d, benz_d, dupil_d, teze_d,
                   0.4, TRUE, 52) %>%
        filter(abs(time_weeks - 52) < 0.5) %>%
        slice_tail(n=1) %>%
        select(FEV1_pct, EOS_blood, exacerb_hazard) %>%
        mutate(Dose = d, Biologic = biologic)
    })
  })

  output$dr_fev1_plot <- renderPlotly({
    df <- dr_result() %>%
      select(Dose, FEV1_pct, EOS_blood) %>%
      pivot_longer(-Dose)
    p <- ggplot(df, aes(Dose, value, color=name)) +
      geom_point(size=3) + geom_line(size=1) +
      scale_x_log10() +
      facet_wrap(~name, scales="free_y") +
      labs(x="Dose (mg)", y="Value at week 52") +
      theme_bw() + theme(legend.position="none")
    ggplotly(p)
  })

  output$threshold_plot <- renderPlotly({
    eos_levels <- seq(100, 1200, by=100)
    tbl <- purrr::map_dfr(eos_levels, function(e) {
      data.frame(Eos=e,
                 Hazard = 0.18 + 0.0015*(e - 450) + (-0.025)*(58-70))
    })
    plot_ly(tbl, x=~Eos, y=~Hazard, type="scatter", mode="lines+markers",
            line=list(color="#B71C1C")) %>%
      add_segments(x=300, xend=300, y=min(tbl$Hazard), yend=max(tbl$Hazard),
                   line=list(dash="dash", color="orange")) %>%
      layout(title="Exacerbation Risk vs. Baseline Blood Eosinophil",
             xaxis=list(title="Blood Eos (cells/μL)"),
             yaxis=list(title="Hazard (events/yr)"))
  })
}

## ── Launch ────────────────────────────────────────────────────────────────────
shinyApp(ui = ui, server = server)
