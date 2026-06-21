################################################################################
## Epilepsy QSP Shiny Dashboard
## 뇌전증 정량적 시스템 약리학 인터랙티브 앱
## Version 1.0 | 2026-06-21
##
## Tabs:
##   1. Patient Profile     — 환자 설정 & 발작 분류
##   2. PK Profiles         — AED 혈장 농도-시간 곡선
##   3. PD Biomarkers       — GABA · 글루타메이트 · Nav 차단 · SV2A
##   4. Clinical Endpoints  — 발작 빈도 · 반응률 · 발작 무 rate
##   5. Scenario Comparison — 단독요법 vs 병용요법 비교
##   6. Drug Resistance & Risk — P-gp 발현 · SUDEP 위험 · TSC/mTOR
################################################################################

library(shiny)
library(mrgsolve)
library(dplyr)
library(ggplot2)
library(tidyr)
library(DT)
library(shinydashboard)
library(plotly)
library(scales)

# ─── Embed QSP model code ─────────────────────────────────────────────────────
epi_code <- '
$PARAM
ka_VPA=1.20, Vc_VPA=9.10, Vp_VPA=12.0, CL_VPA=0.47, Q_VPA=0.23, F_VPA=0.90, fu_VPA=0.10,
ka_LEV=1.50, Vc_LEV=42.0, CL_LEV=3.80, F_LEV=1.00,
ka_CBZ=0.50, Vc_CBZ=51.0, CL_CBZ0=3.00, CL_CBZ_max=6.50, EC50_auto=3.00,
F_CBZ=0.75, fm_CBZ=0.40, CL_CBZE=2.00, Vc_CBZE=20.0,
ka_LTG=0.80, Vc_LTG=77.0, CL_LTG0=1.50, F_LTG=0.98,
DDI_VPA=0.0, DDI_CBZ=0.0,
kin_GABA=0.100, kout_GABA=0.100, GABA0=1.000,
IC50_VPA_GABA=50.0, Imax_VPA_GABA=0.55, n_VPA_GABA=1.20,
kin_SYNAP=0.100, kout_SYNAP=0.100, SYNAP0=1.000,
IC50_LEV_SYNAP=5.0, Imax_LEV_SYNAP=0.45,
Kd_SV2A=4.0, kon_SV2A=0.50, koff_SV2A=2.00,
IC50_CBZ_Nav=3.5, IC50_LTG_Nav=2.0, Emax_Nav=0.90, n_Nav=1.50,
STHRES0=1.000, kthres_rec=0.020,
alpha_GABA=0.350, alpha_Nav=0.500, alpha_SV2A=0.250,
SeizBasal=8.00, k_seiz=12.0,
PGP0=1.000, k_PGP_ind=0.003, k_PGP_deg=0.050,
pgp_CBZ=0.55, pgp_VPA=0.40, pgp_LTG=0.35, pgp_LEV=0.15,
PHT_DOSE=0.0, BZD_BOLUS=0.0,
mTOR_activ=1.0, Emax_mTOR=0.50, EC50_mTOR=0.5, ever_dose=0.0

$CMT AGUT ACENT APER BGUT BCENT CGUT CCENT CMETA DGUT DCENT
     GABA SYNAP SV2A_OCC NAV_BLOCK STHRES PGP

$GLOBAL
#define C_VPA   (ACENT/Vc_VPA)
#define C_LEV   (BCENT/Vc_LEV)
#define C_CBZ   (CCENT/Vc_CBZ)
#define C_CBZE  (CMETA/Vc_CBZE)
#define C_LTG   (DCENT/Vc_LTG)
#define fmaxd(a,b) ((a)>(b)?(a):(b))

$MAIN
double CL_LTG = CL_LTG0;
if(DDI_VPA>0.5) CL_LTG = CL_LTG0*0.50;
if(DDI_CBZ>0.5) CL_LTG = CL_LTG*2.00;
double CL_CBZ = CL_CBZ0 + (CL_CBZ_max-CL_CBZ0)*C_CBZ/(EC50_auto+C_CBZ);
double pgp_factor = fmaxd(1.0, PGP);
double CNS_VPA = (fu_VPA*C_VPA)/(1.0+pgp_VPA*(pgp_factor-1.0));
double CNS_LEV = C_LEV/(1.0+pgp_LEV*(pgp_factor-1.0));
double CNS_CBZ = C_CBZ/(1.0+pgp_CBZ*(pgp_factor-1.0));
double CNS_LTG = C_LTG/(1.0+pgp_LTG*(pgp_factor-1.0));
double CNS_CBZE= C_CBZE/(1.0+pgp_CBZ*(pgp_factor-1.0));
double VPA_pow = pow(CNS_VPA, n_VPA_GABA);
double IC50_pow= pow(IC50_VPA_GABA, n_VPA_GABA);
double Imax_GABA_eff = Imax_VPA_GABA*VPA_pow/(IC50_pow+VPA_pow);
double Imax_SYNAP_eff = Imax_LEV_SYNAP*CNS_LEV/(IC50_LEV_SYNAP+CNS_LEV);
double Nav_CBZ_eff = Emax_Nav*pow(CNS_CBZ+0.5*CNS_CBZE,n_Nav)/(pow(IC50_CBZ_Nav,n_Nav)+pow(CNS_CBZ+0.5*CNS_CBZE,n_Nav));
double Nav_LTG_eff = Emax_Nav*pow(CNS_LTG,n_Nav)/(pow(IC50_LTG_Nav,n_Nav)+pow(CNS_LTG,n_Nav));
double Nav_combined = 1.0-(1.0-Nav_CBZ_eff)*(1.0-Nav_LTG_eff);
double BZD_GABA_boost = (BZD_BOLUS>0.5)?0.60:0.0;
double mTOR_effect = Emax_mTOR*ever_dose/(EC50_mTOR+ever_dose);
double mTOR_thresh_adj = (mTOR_activ>1.0)?mTOR_effect*(mTOR_activ-1.0):0.0;
double SeizFreq = SeizBasal*exp(-k_seiz*(STHRES-STHRES0));
if(SeizFreq<0.0) SeizFreq=0.0;

$ODE
dxdt_AGUT  = -ka_VPA*AGUT;
dxdt_ACENT =  ka_VPA*AGUT - (CL_VPA/Vc_VPA)*ACENT - (Q_VPA/Vc_VPA)*ACENT + (Q_VPA/Vp_VPA)*APER;
dxdt_APER  =  (Q_VPA/Vc_VPA)*ACENT - (Q_VPA/Vp_VPA)*APER;
dxdt_BGUT  = -ka_LEV*BGUT;
dxdt_BCENT =  ka_LEV*BGUT - (CL_LEV/Vc_LEV)*BCENT;
dxdt_CGUT  = -ka_CBZ*CGUT;
dxdt_CCENT =  ka_CBZ*CGUT - (CL_CBZ/Vc_CBZ)*CCENT;
dxdt_CMETA =  fm_CBZ*(CL_CBZ/Vc_CBZ)*CCENT - (CL_CBZE/Vc_CBZE)*CMETA;
dxdt_DGUT  = -ka_LTG*DGUT;
dxdt_DCENT =  ka_LTG*DGUT - (CL_LTG/Vc_LTG)*DCENT;
dxdt_GABA  =  kin_GABA*(1.0+Imax_GABA_eff+BZD_GABA_boost) - kout_GABA*GABA;
dxdt_SYNAP =  kin_SYNAP*(1.0-Imax_SYNAP_eff) - kout_SYNAP*SYNAP;
dxdt_SV2A_OCC = kon_SV2A*CNS_LEV*(1.0-SV2A_OCC) - koff_SV2A*SV2A_OCC;
dxdt_NAV_BLOCK = 5.0*(Nav_combined-NAV_BLOCK);
double STHRES_target = STHRES0+alpha_GABA*(GABA-GABA0)+alpha_Nav*NAV_BLOCK+alpha_SV2A*SV2A_OCC+mTOR_thresh_adj;
dxdt_STHRES = kthres_rec*(STHRES_target-STHRES);
double PGP_drive = k_PGP_ind*SeizFreq;
dxdt_PGP = PGP_drive - k_PGP_deg*(PGP-PGP0);

$TABLE
capture C_VPA_mcg=C_VPA; capture C_LEV_mcg=C_LEV; capture C_CBZ_mcg=C_CBZ;
capture C_CBZE_mcg=C_CBZE; capture C_LTG_mcg=C_LTG;
capture CNS_VPA_f=CNS_VPA; capture CNS_CBZ_f=CNS_CBZ;
capture GABA_norm=GABA; capture SYNAP_norm=SYNAP;
capture SV2A_frac=SV2A_OCC; capture Nav_frac=NAV_BLOCK;
capture Thresh=STHRES; capture PGP_exp=PGP;
capture SeizFreq_obs=SeizFreq;
capture Responder=(SeizFreq<SeizBasal*0.50)?1.0:0.0;
capture SeizFree=(SeizFreq<0.1)?1.0:0.0;
'

mod <- mcode("epi_app", epi_code, quiet=TRUE)

run_sim <- function(dose_VPA, dose_LEV, dose_CBZ, dose_LTG,
                    ddi_VPA, ddi_CBZ, pgp_init, mTOR_val, ever_val,
                    seiz_basal, duration_days=180) {
  e <- ev()
  if (dose_VPA > 0) e <- e + ev(cmt="AGUT", amt=dose_VPA/2, ii=12, addl=duration_days*2-1)
  if (dose_LEV > 0) e <- e + ev(cmt="BGUT", amt=dose_LEV/2, ii=12, addl=duration_days*2-1)
  if (dose_CBZ > 0) e <- e + ev(cmt="CGUT", amt=dose_CBZ/2, ii=12, addl=duration_days*2-1)
  if (dose_LTG > 0) e <- e + ev(cmt="DGUT", amt=dose_LTG/2, ii=12, addl=duration_days*2-1)
  ini <- init(mod, GABA=1.0, SYNAP=1.0, STHRES=1.0,
              PGP=pgp_init, SV2A_OCC=0, NAV_BLOCK=0)
  mrgsim(ini, events=e, end=duration_days*24, delta=2,
         param(mod, DDI_VPA=ddi_VPA, DDI_CBZ=ddi_CBZ,
               mTOR_activ=mTOR_val, ever_dose=ever_val,
               SeizBasal=seiz_basal)) %>%
    as_tibble() %>%
    mutate(day = time / 24)
}

# ─── UI ───────────────────────────────────────────────────────────────────────
ui <- dashboardPage(
  skin = "purple",
  dashboardHeader(title = "Epilepsy QSP Model | 뇌전증 정량적 시스템 약리학"),
  dashboardSidebar(
    sidebarMenu(
      menuItem("Patient Profile",     tabName="tab_patient",   icon=icon("user-circle")),
      menuItem("PK Profiles",         tabName="tab_pk",        icon=icon("chart-line")),
      menuItem("PD Biomarkers",       tabName="tab_pd",        icon=icon("brain")),
      menuItem("Clinical Endpoints",  tabName="tab_clinical",  icon=icon("heartbeat")),
      menuItem("Scenario Comparison", tabName="tab_scenario",  icon=icon("balance-scale")),
      menuItem("Resistance & Risk",   tabName="tab_resistance",icon=icon("shield-alt"))
    )
  ),
  dashboardBody(
    tags$head(tags$style(HTML("
      .box { border-top: 3px solid #8B5CF6 !important; }
      .info-box { background: linear-gradient(135deg, #667eea 0%, #764ba2 100%); }
    "))),
    tabItems(

      # ── TAB 1: PATIENT PROFILE ──────────────────────────────────────────────
      tabItem(tabName="tab_patient",
        fluidRow(
          box(title="Patient Demographics & Epilepsy Classification",
              status="primary", solidHeader=TRUE, width=4,
              selectInput("seizure_type", "Epilepsy Type / Seizure Classification",
                choices=c("Focal Onset Impaired Awareness (TLE)"="focal_tia",
                          "Focal to Bilateral Tonic-Clonic"="focal_bilat",
                          "Generalized Tonic-Clonic (GTCS)"="gen_tc",
                          "Absence Epilepsy (CAE/JAE)"="absence",
                          "Juvenile Myoclonic Epilepsy (JME)"="jme",
                          "Dravet Syndrome (SCN1A)"="dravet",
                          "Lennox-Gastaut Syndrome"="lgs",
                          "TSC-related Epilepsy"="tsc"),
                selected="focal_tia"),
              numericInput("age",    "Age (years)",      value=35, min=2,  max=80),
              numericInput("weight", "Weight (kg)",      value=70, min=20, max=150),
              selectInput("sex", "Sex",
                choices=c("Male"="M","Female"="F"), selected="M"),
              sliderInput("seiz_history", "Pre-treatment seizure frequency (per month)",
                min=1, max=30, value=8, step=1)
          ),
          box(title="Genetic & Biomarker Profile",
              status="warning", solidHeader=TRUE, width=4,
              selectInput("genetics", "Genetic Testing Result",
                choices=c("Not tested / Unknown"="none",
                          "SCN1A pathogenic variant (Dravet)"="scn1a",
                          "KCNQ2 pathogenic variant (BFNE)"="kcnq2",
                          "TSC1/TSC2 mutation (TSC)"="tsc",
                          "DEPDC5 variant (focal)"="depdc5",
                          "GABRG2 variant (GEFS+)"="gabrg2",
                          "ABCB1 C3435T (P-gp high)"="abcb1_tt",
                          "Structural (MRI positive)"="structural"),
                selected="none"),
              selectInput("eeg_type", "EEG Pattern",
                choices=c("Focal temporal IED"="focal_temp",
                          "Generalized 3Hz spike-wave (absence)"="gen_3hz",
                          "Multifocal / diffuse"="multifocal",
                          "Normal EEG"="normal"),
                selected="focal_temp"),
              selectInput("mri_result", "MRI Finding",
                choices=c("Normal"="normal",
                          "Mesial temporal sclerosis (MTS)"="mts",
                          "Focal cortical dysplasia (FCD)"="fcd",
                          "Tuberous sclerosis (TSC)"="tsc_mri",
                          "Post-stroke lesion"="stroke",
                          "Cortical tumor"="tumor"),
                selected="normal"),
              checkboxInput("prev_aed_failed", "Previous AED trial failed (potential DRE)", FALSE),
              conditionalPanel("input.prev_aed_failed",
                numericInput("n_aed_failed", "Number of AEDs failed", value=2, min=1, max=10))
          ),
          box(title="Comorbidities & Risk Factors",
              status="danger", solidHeader=TRUE, width=4,
              checkboxGroupInput("comorbidities", "Active Comorbidities",
                choices=c("Depression"="depression",
                          "Anxiety disorder"="anxiety",
                          "Cognitive impairment"="cognitive",
                          "Osteoporosis (AED-induced)"="osteoporosis",
                          "Pregnancy (women of childbearing age)"="pregnancy",
                          "Renal impairment (CKD)"="ckd",
                          "Hepatic impairment"="liver")),
              sliderInput("sudep_risk_base", "Background SUDEP risk factor",
                min=1, max=10, value=3, step=1,
                post=" × baseline"),
              selectInput("aed_history", "Previous AED Experience",
                choices=c("Drug-naïve"="naive",
                          "1 prior AED failed"="1aed",
                          "2+ prior AEDs failed (DRE)"="dre",
                          "Post-surgical"="postsurg")),
              verbatimTextOutput("patient_summary")
          )
        ),
        fluidRow(
          valueBoxOutput("vbox_seiz_type"),
          valueBoxOutput("vbox_dre_risk"),
          valueBoxOutput("vbox_sudep")
        )
      ),

      # ── TAB 2: PK PROFILES ──────────────────────────────────────────────────
      tabItem(tabName="tab_pk",
        fluidRow(
          box(title="AED Dosing Configuration",
              status="primary", solidHeader=TRUE, width=3,
              h4("Valproate (VPA)"),
              sliderInput("dose_VPA", "Daily dose (mg/day)", 0, 3000, 1000, step=100),
              h4("Levetiracetam (LEV)"),
              sliderInput("dose_LEV", "Daily dose (mg/day)", 0, 4000, 0, step=250),
              h4("Carbamazepine (CBZ)"),
              sliderInput("dose_CBZ", "Daily dose (mg/day)", 0, 1800, 0, step=100),
              h4("Lamotrigine (LTG)"),
              sliderInput("dose_LTG", "Daily dose (mg/day)", 0, 600, 0, step=25),
              hr(),
              checkboxInput("ddi_VPA_flag", "VPA co-administered (DDI: LTG t½ doubles)", FALSE),
              checkboxInput("ddi_CBZ_flag", "CBZ co-administered (DDI: LTG CL doubles)", FALSE),
              sliderInput("sim_days", "Simulation duration (days)", 30, 365, 180, step=30),
              actionButton("run_btn", "Run Simulation", icon=icon("play"),
                           style="background-color:#8B5CF6; color:white; width:100%")
          ),
          box(title="Plasma Concentration-Time Profiles",
              status="info", solidHeader=TRUE, width=9,
              tabsetPanel(
                tabPanel("VPA",  plotlyOutput("pk_VPA",  height="320px")),
                tabPanel("LEV",  plotlyOutput("pk_LEV",  height="320px")),
                tabPanel("CBZ",  plotlyOutput("pk_CBZ",  height="320px")),
                tabPanel("LTG",  plotlyOutput("pk_LTG",  height="320px")),
                tabPanel("All",  plotlyOutput("pk_all",  height="320px")),
                tabPanel("SS Detail (last 48h)", plotlyOutput("pk_ss", height="320px"))
              )
          )
        )
      ),

      # ── TAB 3: PD BIOMARKERS ────────────────────────────────────────────────
      tabItem(tabName="tab_pd",
        fluidRow(
          box(title="Brain GABA Level (VPA — GABA-T Inhibition)",
              status="success", solidHeader=TRUE, width=6,
              plotlyOutput("pd_GABA", height="280px")),
          box(title="Synaptic Glutamate (LEV — SV2A/Vesicle Release)",
              status="warning", solidHeader=TRUE, width=6,
              plotlyOutput("pd_SYNAP", height="280px"))
        ),
        fluidRow(
          box(title="SV2A Occupancy — Levetiracetam",
              status="info", solidHeader=TRUE, width=6,
              plotlyOutput("pd_SV2A", height="280px")),
          box(title="Sodium Channel Blockade — CBZ + LTG + PHT",
              status="primary", solidHeader=TRUE, width=6,
              plotlyOutput("pd_Nav", height="280px"))
        ),
        fluidRow(
          box(title="PD Parameter Settings",
              status="default", solidHeader=TRUE, width=4,
              sliderInput("IC50_VPA_s", "VPA IC50 for GABA-T (mcg/mL)", 10, 150, 50, step=5),
              sliderInput("IC50_CBZ_s", "CBZ IC50 for Nav (mcg/mL)",     1,  10,  3.5, step=0.5),
              sliderInput("IC50_LTG_s", "LTG IC50 for Nav (mcg/mL)",     0.5, 8,  2.0, step=0.5),
              sliderInput("IC50_LEV_s", "LEV IC50 for SV2A (mcg/mL)",    1,  20,  4.0, step=1)
          ),
          box(title="PD Biomarker Summary (Steady-State)",
              status="default", solidHeader=TRUE, width=8,
              DTOutput("pd_summary_tbl"))
        )
      ),

      # ── TAB 4: CLINICAL ENDPOINTS ────────────────────────────────────────────
      tabItem(tabName="tab_clinical",
        fluidRow(
          valueBoxOutput("vbox_seiz_freq"),
          valueBoxOutput("vbox_responder"),
          valueBoxOutput("vbox_seizfree")
        ),
        fluidRow(
          box(title="Seizure Frequency Over Time (episodes/28 days)",
              status="danger", solidHeader=TRUE, width=8,
              plotlyOutput("clinical_seiz_freq", height="350px")),
          box(title="Seizure Threshold (Normalized)",
              status="warning", solidHeader=TRUE, width=4,
              plotlyOutput("clinical_threshold", height="350px"))
        ),
        fluidRow(
          box(title="AED Therapeutic Windows",
              status="info", solidHeader=TRUE, width=6,
              DTOutput("therapeutic_range_tbl")),
          box(title="Response Probability Over Time",
              status="success", solidHeader=TRUE, width=6,
              plotlyOutput("response_prob_plot", height="280px"))
        )
      ),

      # ── TAB 5: SCENARIO COMPARISON ──────────────────────────────────────────
      tabItem(tabName="tab_scenario",
        fluidRow(
          box(title="Select Scenarios for Comparison",
              status="primary", solidHeader=TRUE, width=3,
              checkboxGroupInput("scenarios_sel", "Include scenarios",
                choices=c("Untreated (baseline)"="s_untreated",
                          "VPA 1,000 mg/day"="s_vpa",
                          "LEV 3,000 mg/day"="s_lev",
                          "CBZ 600 mg/day"="s_cbz",
                          "LTG 200 mg/day"="s_ltg",
                          "VPA 500 + LTG 100 mg/day (DDI)"="s_vpa_ltg",
                          "CBZ 600 + LTG 400 mg/day (DDI)"="s_cbz_ltg",
                          "DRE: VPA (P-gp 3×)"="s_dre"),
                selected=c("s_untreated","s_vpa","s_lev","s_cbz")),
              hr(),
              radioButtons("comp_metric", "Primary comparison metric",
                choices=c("Seizure Frequency"="SeizFreq_obs",
                          "Seizure Threshold"="Thresh",
                          "Brain GABA"="GABA_norm",
                          "Nav Blockade"="Nav_frac",
                          "P-gp Expression"="PGP_exp"),
                selected="SeizFreq_obs"),
              actionButton("run_compare", "Compare Scenarios",
                           icon=icon("chart-bar"),
                           style="background-color:#7C3AED; color:white; width:100%")
          ),
          box(title="Scenario Comparison Plot",
              status="info", solidHeader=TRUE, width=9,
              plotlyOutput("scenario_compare_plot", height="450px"))
        ),
        fluidRow(
          box(title="Comparative Efficacy Summary (Days 150–180)",
              status="success", solidHeader=TRUE, width=12,
              DTOutput("scenario_summary_tbl"))
        )
      ),

      # ── TAB 6: RESISTANCE & RISK ─────────────────────────────────────────────
      tabItem(tabName="tab_resistance",
        fluidRow(
          box(title="Drug Resistance Settings",
              status="danger", solidHeader=TRUE, width=3,
              h4("P-glycoprotein (ABCB1/MDR1)"),
              sliderInput("pgp_init", "Initial P-gp expression level",
                min=1.0, max=5.0, value=1.0, step=0.5),
              helpText("1.0 = normal; 2.0 = moderate overexpression; 3+ = DRE-level"),
              sliderInput("pgp_ind_rate", "P-gp induction rate (seizure-driven)",
                min=0.001, max=0.010, value=0.003, step=0.001),
              hr(),
              h4("mTOR Pathway (TSC/FCD)"),
              sliderInput("mTOR_val", "mTOR pathway activity (1=normal, 2=TSC-level)",
                min=1.0, max=3.0, value=1.0, step=0.1),
              sliderInput("ever_dose", "Everolimus relative dose (0=none, 1=full)",
                min=0.0, max=1.0, value=0.0, step=0.1),
              hr(),
              h4("SUDEP Risk Calculator"),
              checkboxInput("gtcs_nocturnal", "Nocturnal GTCS present", FALSE),
              checkboxInput("poor_seizure_control", "Poor seizure control (>1/month)", TRUE),
              checkboxInput("non_compliant", "Non-adherent to AED", FALSE),
              selectInput("aed_combo_risk", "AED regimen",
                choices=c("Monotherapy"="mono","Polytherapy"="poly","None"="none"),
                selected="mono")
          ),
          box(title="P-glycoprotein Expression Trajectory",
              status="danger", solidHeader=TRUE, width=5,
              plotlyOutput("pgp_plot", height="320px")),
          box(title="BBB AED Exposure Ratio (CNS/Plasma)",
              status="warning", solidHeader=TRUE, width=4,
              plotlyOutput("bbb_plot", height="320px"))
        ),
        fluidRow(
          box(title="SUDEP Risk Score & Assessment",
              status="danger", solidHeader=TRUE, width=4,
              uiOutput("sudep_risk_ui")),
          box(title="mTOR Pathway: TSC/FCD Seizure Control",
              status="purple", solidHeader=TRUE, width=8,
              plotlyOutput("mtor_plot", height="300px"))
        ),
        fluidRow(
          box(title="AED Safety & Monitoring Summary",
              status="info", solidHeader=TRUE, width=12,
              DTOutput("safety_tbl"))
        )
      )

    )  # end tabItems
  )    # end dashboardBody
)

# ─── SERVER ───────────────────────────────────────────────────────────────────
server <- function(input, output, session) {

  # Reactive: run base simulation when button pressed
  sim_data <- eventReactive(input$run_btn, {
    withProgress(message="Running ODE simulation...", value=0.5, {
      run_sim(
        dose_VPA  = input$dose_VPA,
        dose_LEV  = input$dose_LEV,
        dose_CBZ  = input$dose_CBZ,
        dose_LTG  = input$dose_LTG,
        ddi_VPA   = as.numeric(input$ddi_VPA_flag),
        ddi_CBZ   = as.numeric(input$ddi_CBZ_flag),
        pgp_init  = input$pgp_init,
        mTOR_val  = input$mTOR_val,
        ever_val  = input$ever_dose,
        seiz_basal = input$seiz_history,
        duration_days = input$sim_days
      )
    })
  }, ignoreNULL=FALSE)

  # ── Tab 1: Value boxes ──────────────────────────────────────────────────────
  output$vbox_seiz_type <- renderValueBox({
    label <- switch(input$seizure_type,
      focal_tia="Focal TLE", focal_bilat="Focal→Bilateral",
      gen_tc="Generalized GTCS", absence="Absence",
      jme="JME", dravet="Dravet", lgs="LGS", tsc="TSC", "Unknown")
    valueBox(label, "Epilepsy Type", icon=icon("brain"), color="purple")
  })
  output$vbox_dre_risk <- renderValueBox({
    risk <- if (input$prev_aed_failed && input$n_aed_failed >= 2) "HIGH" else
            if (input$prev_aed_failed) "MODERATE" else "LOW"
    valueBox(risk, "DRE Risk", icon=icon("exclamation-triangle"),
             color=if(risk=="HIGH") "red" else if(risk=="MODERATE") "yellow" else "green")
  })
  output$vbox_sudep <- renderValueBox({
    base_risk <- input$sudep_risk_base
    seiz_mult <- if (input$seiz_history >= 12) 10 else if (input$seiz_history >= 4) 4 else 2
    sudep_rate <- round(base_risk * seiz_mult / 10, 1)
    valueBox(paste0(sudep_rate, "×"), "Relative SUDEP Risk",
             icon=icon("heartbeat"), color="red")
  })
  output$patient_summary <- renderText({
    paste0(
      "Patient: ", input$age, "y ", input$sex, ", ", input$weight, "kg\n",
      "Epilepsy: ", input$seizure_type, "\n",
      "Genetics: ", input$genetics, "\n",
      "EEG: ", input$eeg_type, "\n",
      "MRI: ", input$mri_result, "\n",
      "Pre-Rx seizures: ", input$seiz_history, "/month\n",
      "Comorbidities: ", paste(input$comorbidities, collapse=", "), "\n",
      "AED history: ", input$aed_history
    )
  })

  # ── Tab 2: PK plots ─────────────────────────────────────────────────────────
  output$pk_VPA <- renderPlotly({
    df <- sim_data()
    p <- ggplot(df, aes(x=day, y=C_VPA_mcg)) +
      geom_line(color="#8B5CF6") +
      geom_hline(yintercept=c(50,100), linetype="dashed", color=c("green","red")) +
      annotate("text", x=max(df$day)*0.9, y=55, label="Min TW (50)", size=3, color="darkgreen") +
      annotate("text", x=max(df$day)*0.9, y=105, label="Max TW (100)", size=3, color="red") +
      labs(title="VPA Plasma Concentration", x="Day", y="Concentration (mcg/mL)") +
      theme_bw(base_size=11)
    ggplotly(p, tooltip=c("x","y"))
  })
  output$pk_LEV <- renderPlotly({
    df <- sim_data()
    p <- ggplot(df, aes(x=day, y=C_LEV_mcg)) +
      geom_line(color="#2563EB") +
      geom_hline(yintercept=c(12,46), linetype="dashed", color=c("green","red")) +
      labs(title="LEV Plasma Concentration", x="Day", y="Concentration (mcg/mL)") +
      theme_bw(base_size=11)
    ggplotly(p)
  })
  output$pk_CBZ <- renderPlotly({
    df <- sim_data()
    p <- ggplot(df %>% pivot_longer(c(C_CBZ_mcg, C_CBZE_mcg),
                names_to="species", values_to="conc")) +
      geom_line(aes(x=day, y=conc, color=species)) +
      geom_hline(yintercept=c(4,12), linetype="dashed", color=c("green","red")) +
      scale_color_manual(values=c("#D97706","#B45309"), labels=c("CBZ","CBZ-E")) +
      labs(title="CBZ + CBZ-Epoxide Plasma", x="Day", y="Concentration (mcg/mL)") +
      theme_bw(base_size=11)
    ggplotly(p)
  })
  output$pk_LTG <- renderPlotly({
    df <- sim_data()
    label <- if (input$ddi_VPA_flag) " (DDI: VPA↑ LTG by 2×)" else
             if (input$ddi_CBZ_flag) " (DDI: CBZ↓ LTG by 0.5×)" else ""
    p <- ggplot(df, aes(x=day, y=C_LTG_mcg)) +
      geom_line(color="#059669") +
      geom_hline(yintercept=c(3,14), linetype="dashed", color=c("green","red")) +
      labs(title=paste0("LTG Plasma Concentration", label), x="Day", y="Concentration (mcg/mL)") +
      theme_bw(base_size=11)
    ggplotly(p)
  })
  output$pk_all <- renderPlotly({
    df <- sim_data() %>%
      select(day, C_VPA_mcg, C_LEV_mcg, C_CBZ_mcg, C_LTG_mcg) %>%
      pivot_longer(-day, names_to="AED", values_to="conc") %>%
      mutate(AED = gsub("_mcg","",AED) %>% gsub("C_","",.) )
    p <- ggplot(df, aes(x=day, y=conc, color=AED)) +
      geom_line() +
      scale_color_brewer(palette="Set1") +
      labs(title="All AED Plasma Concentrations", x="Day", y="mcg/mL") +
      theme_bw(base_size=11)
    ggplotly(p)
  })
  output$pk_ss <- renderPlotly({
    df <- sim_data()
    ss_day <- max(df$day, na.rm=TRUE) - 2
    df_ss  <- df %>% filter(day >= ss_day)
    p <- ggplot(df_ss %>%
      select(day, C_VPA_mcg, C_LEV_mcg, C_CBZ_mcg, C_LTG_mcg) %>%
      pivot_longer(-day, names_to="AED", values_to="conc") %>%
      mutate(AED=gsub("C_|_mcg","",AED)),
      aes(x=day, y=conc, color=AED)) +
      geom_line(linewidth=1.2) +
      scale_color_brewer(palette="Set1") +
      labs(title="Steady-State AED Concentrations (Last 48h)", x="Day", y="mcg/mL") +
      theme_bw(base_size=11)
    ggplotly(p)
  })

  # ── Tab 3: PD Biomarkers ─────────────────────────────────────────────────────
  output$pd_GABA <- renderPlotly({
    df <- sim_data()
    p <- ggplot(df, aes(x=day, y=GABA_norm)) +
      geom_line(color="#16A34A", linewidth=1.1) +
      geom_hline(yintercept=1.0, linetype="dashed", color="gray50") +
      annotate("text", x=max(df$day)*0.5, y=1.03, label="Baseline GABA", size=3) +
      labs(title="Brain GABA Level (VPA → GABA-T Inhibition)",
           x="Day", y="GABA (normalized)") +
      theme_bw(base_size=11)
    ggplotly(p)
  })
  output$pd_SYNAP <- renderPlotly({
    df <- sim_data()
    p <- ggplot(df, aes(x=day, y=SYNAP_norm)) +
      geom_line(color="#D97706", linewidth=1.1) +
      geom_hline(yintercept=1.0, linetype="dashed", color="gray50") +
      labs(title="Synaptic Glutamate (LEV → SV2A → Release↓)",
           x="Day", y="Synaptic Glu (normalized)") +
      theme_bw(base_size=11)
    ggplotly(p)
  })
  output$pd_SV2A <- renderPlotly({
    df <- sim_data()
    p <- ggplot(df, aes(x=day, y=SV2A_frac*100)) +
      geom_line(color="#2563EB", linewidth=1.1) +
      labs(title="SV2A Occupancy (LEV Binding Fraction)",
           x="Day", y="SV2A Occupancy (%)") +
      scale_y_continuous(limits=c(0,100)) +
      theme_bw(base_size=11)
    ggplotly(p)
  })
  output$pd_Nav <- renderPlotly({
    df <- sim_data()
    p <- ggplot(df, aes(x=day, y=Nav_frac*100)) +
      geom_line(color="#7C3AED", linewidth=1.1) +
      labs(title="Sodium Channel Blockade (CBZ + LTG Combined)",
           x="Day", y="Nav Channel Block (%)") +
      scale_y_continuous(limits=c(0,100)) +
      theme_bw(base_size=11)
    ggplotly(p)
  })
  output$pd_summary_tbl <- renderDT({
    df <- sim_data()
    ss <- df %>% filter(day >= max(day)*0.85) %>%
      summarise(
        `VPA C (mcg/mL)`   = round(mean(C_VPA_mcg), 1),
        `LEV C (mcg/mL)`   = round(mean(C_LEV_mcg), 1),
        `CBZ C (mcg/mL)`   = round(mean(C_CBZ_mcg), 1),
        `LTG C (mcg/mL)`   = round(mean(C_LTG_mcg), 1),
        `Brain GABA (norm)` = round(mean(GABA_norm), 3),
        `Synaptic Glu (norm)` = round(mean(SYNAP_norm), 3),
        `SV2A Occ (%)`     = round(mean(SV2A_frac)*100, 1),
        `Nav Block (%)`    = round(mean(Nav_frac)*100, 1),
        `Threshold`        = round(mean(Thresh), 3),
        `P-gp (norm)`      = round(mean(PGP_exp), 3)
      )
    datatable(ss, options=list(dom="t", scrollX=TRUE), rownames=FALSE)
  })

  # ── Tab 4: Clinical Endpoints ─────────────────────────────────────────────
  output$vbox_seiz_freq <- renderValueBox({
    df <- sim_data()
    ss_freq <- df %>% filter(day >= max(day)*0.85) %>%
      summarise(freq=mean(SeizFreq_obs)) %>% pull(freq)
    valueBox(round(ss_freq, 1), "SS Seizure Freq (per 28d)",
             icon=icon("bolt"), color=if(ss_freq<1)"green" else if(ss_freq<4)"yellow" else "red")
  })
  output$vbox_responder <- renderValueBox({
    df <- sim_data()
    resp <- df %>% filter(day >= max(day)*0.85) %>%
      summarise(r=mean(Responder)*100) %>% pull(r)
    valueBox(paste0(round(resp,0), "%"), "Responder Rate (≥50% reduction)",
             icon=icon("check-circle"), color=if(resp>50)"green" else "yellow")
  })
  output$vbox_seizfree <- renderValueBox({
    df <- sim_data()
    sf  <- df %>% filter(day >= max(day)*0.85) %>%
      summarise(sf=mean(SeizFree)*100) %>% pull(sf)
    valueBox(paste0(round(sf,0), "%"), "Seizure Freedom Rate",
             icon=icon("star"), color=if(sf>70)"green" else if(sf>30)"yellow" else "red")
  })
  output$clinical_seiz_freq <- renderPlotly({
    df <- sim_data()
    p  <- ggplot(df, aes(x=day, y=SeizFreq_obs)) +
      geom_line(color="#DC2626", linewidth=1.2) +
      geom_hline(yintercept=input$seiz_history*0.5, linetype="dashed", color="green3") +
      annotate("text", x=max(df$day)*0.7, y=input$seiz_history*0.5+0.3,
               label="50% reduction goal", size=3, color="darkgreen") +
      geom_hline(yintercept=0.1, linetype="dotted", color="blue3") +
      annotate("text", x=max(df$day)*0.7, y=0.4, label="Seizure freedom threshold", size=3, color="blue") +
      labs(title="Seizure Frequency Over Time",
           x="Day", y="Seizure frequency (episodes/28 days)") +
      theme_bw(base_size=11)
    ggplotly(p)
  })
  output$clinical_threshold <- renderPlotly({
    df <- sim_data()
    p  <- ggplot(df, aes(x=day, y=Thresh)) +
      geom_line(color="#7C3AED", linewidth=1.1) +
      geom_hline(yintercept=1.0, linetype="dashed", color="gray50") +
      labs(title="Seizure Threshold (STHRES)", x="Day", y="Threshold (normalized)") +
      theme_bw(base_size=11)
    ggplotly(p)
  })
  output$therapeutic_range_tbl <- renderDT({
    tbl <- data.frame(
      AED=c("VPA","LEV","CBZ","LTG","PHT","GBP"),
      `Min TW (mcg/mL)`=c(50,12,4,3,10,2),
      `Max TW (mcg/mL)`=c(100,46,12,14,20,20),
      `t1/2 (h)`=c(9:17,7,8:12,25:36,22,5:7),
      `Primary mechanism`=c("GABA-T inhib","SV2A","Nav slow inact",
                             "Nav state block","Nav inact","Cav α2δ-1"),
      stringsAsFactors=FALSE, check.names=FALSE
    )
    datatable(tbl, options=list(dom="t"), rownames=FALSE)
  })
  output$response_prob_plot <- renderPlotly({
    df <- sim_data()
    daily_resp <- df %>% group_by(day=round(day)) %>%
      summarise(resp_pct=mean(Responder)*100, .groups="drop")
    p <- ggplot(daily_resp, aes(x=day, y=resp_pct)) +
      geom_line(color="#059669", linewidth=1.1) +
      scale_y_continuous(limits=c(0,100), labels=percent_format(scale=1)) +
      labs(title="Responder Probability Over Time",
           x="Day", y="% Time as Responder") +
      theme_bw(base_size=11)
    ggplotly(p)
  })

  # ── Tab 5: Scenario Comparison ────────────────────────────────────────────
  compare_data <- eventReactive(input$run_compare, {
    withProgress(message="Running all scenarios...", value=0.2, {
      scenarios <- list()
      dur <- 180
      base <- input$seiz_history
      if ("s_untreated" %in% input$scenarios_sel) {
        scenarios[["Untreated"]] <- run_sim(0,0,0,0,0,0,1,1,0,base,dur)
      }
      if ("s_vpa" %in% input$scenarios_sel) {
        scenarios[["VPA 1,000"]] <- run_sim(1000,0,0,0,0,0,1,1,0,base,dur)
      }
      if ("s_lev" %in% input$scenarios_sel) {
        scenarios[["LEV 3,000"]] <- run_sim(0,3000,0,0,0,0,1,1,0,base,dur)
      }
      if ("s_cbz" %in% input$scenarios_sel) {
        scenarios[["CBZ 600"]]   <- run_sim(0,0,600,0,0,0,1,1,0,base,dur)
      }
      if ("s_ltg" %in% input$scenarios_sel) {
        scenarios[["LTG 200"]]   <- run_sim(0,0,0,200,0,0,1,1,0,base,dur)
      }
      if ("s_vpa_ltg" %in% input$scenarios_sel) {
        scenarios[["VPA+LTG(DDI)"]] <- run_sim(500,0,0,100,1,0,1,1,0,base,dur)
      }
      if ("s_cbz_ltg" %in% input$scenarios_sel) {
        scenarios[["CBZ+LTG(DDI)"]] <- run_sim(0,0,600,400,0,1,1,1,0,base,dur)
      }
      if ("s_dre" %in% input$scenarios_sel) {
        scenarios[["DRE(P-gp 3×)"]] <- run_sim(1000,0,0,0,0,0,3,1,0,base,dur)
      }
      bind_rows(lapply(names(scenarios), function(nm) {
        scenarios[[nm]] %>% mutate(scenario=nm)
      }))
    })
  })
  output$scenario_compare_plot <- renderPlotly({
    df <- compare_data()
    metric <- input$comp_metric
    p <- ggplot(df, aes_string(x="day", y=metric, color="scenario")) +
      geom_line(linewidth=0.9) +
      scale_color_brewer(palette="Dark2") +
      labs(title=paste("Scenario Comparison —", metric),
           x="Day", y=metric, color="Scenario") +
      theme_bw(base_size=11) +
      theme(legend.position="bottom") +
      guides(color=guide_legend(nrow=2))
    ggplotly(p) %>% layout(legend=list(orientation="h"))
  })
  output$scenario_summary_tbl <- renderDT({
    df <- compare_data()
    ss <- df %>% filter(day >= 150) %>%
      group_by(Scenario=scenario) %>%
      summarise(
        `Seizure Freq`      = round(mean(SeizFreq_obs), 2),
        `Responder (%)`     = round(mean(Responder)*100, 1),
        `Seizure Free (%)`  = round(mean(SeizFree)*100, 1),
        `STHRES`            = round(mean(Thresh), 3),
        `Brain GABA`        = round(mean(GABA_norm), 3),
        `Nav Block (%)`     = round(mean(Nav_frac)*100, 1),
        `SV2A Occ (%)`      = round(mean(SV2A_frac)*100, 1),
        `P-gp`              = round(mean(PGP_exp), 2),
        .groups="drop"
      )
    datatable(ss, options=list(dom="Bt", scrollX=TRUE, buttons=c("csv","excel")),
              extensions="Buttons", rownames=FALSE)
  })

  # ── Tab 6: Drug Resistance & SUDEP ──────────────────────────────────────────
  resist_data <- reactive({
    run_sim(input$dose_VPA, input$dose_LEV, input$dose_CBZ, input$dose_LTG,
            as.numeric(input$ddi_VPA_flag), as.numeric(input$ddi_CBZ_flag),
            input$pgp_init, input$mTOR_val, input$ever_dose,
            input$seiz_history, 180)
  })
  output$pgp_plot <- renderPlotly({
    df <- resist_data()
    p <- ggplot(df, aes(x=day, y=PGP_exp)) +
      geom_line(color="#DC2626", linewidth=1.1) +
      geom_hline(yintercept=c(2,3), linetype="dashed", color=c("orange","red")) +
      annotate("text",x=150,y=2.15,label="Moderate overexpression",size=3,color="orange") +
      annotate("text",x=150,y=3.15,label="DRE level (3×)",size=3,color="red") +
      labs(title="P-glycoprotein Expression (ABCB1/MDR1)",
           x="Day", y="P-gp (normalized)") +
      theme_bw(base_size=11)
    ggplotly(p)
  })
  output$bbb_plot <- renderPlotly({
    df <- resist_data() %>%
      mutate(
        BBB_VPA  = CNS_VPA_f / pmax(C_VPA_mcg * 0.1, 1e-6),
        BBB_CBZ  = CNS_CBZ_f / pmax(C_CBZ_mcg, 1e-6)
      ) %>%
      select(day, BBB_VPA, BBB_CBZ) %>%
      pivot_longer(-day, names_to="AED", values_to="ratio")
    p <- ggplot(df, aes(x=day, y=ratio, color=AED)) +
      geom_line(linewidth=1.1) +
      scale_color_manual(values=c("#8B5CF6","#0891B2")) +
      labs(title="CNS/Plasma Ratio Over Time (BBB Penetration)",
           x="Day", y="CNS/Plasma exposure ratio") +
      theme_bw(base_size=11)
    ggplotly(p)
  })
  output$sudep_risk_ui <- renderUI({
    base_rate <- 1 / 1000  # ~1 per 1,000 person-years in general epilepsy
    multiplier <- 1
    if (input$gtcs_nocturnal)     multiplier <- multiplier * 3.0
    if (input$poor_seizure_control) multiplier <- multiplier * 4.0
    if (input$non_compliant)       multiplier <- multiplier * 2.0
    if (input$aed_combo_risk == "none") multiplier <- multiplier * 5.0
    if (input$aed_combo_risk == "mono") multiplier <- multiplier * 1.5
    risk_pct <- round(base_rate * multiplier * 100, 3)
    color <- if (multiplier > 10) "red" else if (multiplier > 4) "orange" else "green"
    tagList(
      h4("SUDEP Risk Estimate"),
      p("Risk factors identified:", style=paste0("color:", color)),
      tags$ul(
        if (input$gtcs_nocturnal) tags$li("Nocturnal GTCS (×3.0)"),
        if (input$poor_seizure_control) tags$li("Poor seizure control (×4.0)"),
        if (input$non_compliant) tags$li("Non-adherence (×2.0)"),
        if (input$aed_combo_risk == "none") tags$li("No AED treatment (×5.0)")
      ),
      tags$p(paste0("Estimated annual SUDEP risk: ~", risk_pct, "% per year"),
             style=paste0("font-weight:bold; color:", color)),
      tags$p("Key prevention: seizure freedom, especially nocturnal GTCS control",
             style="font-style:italic; font-size:0.9em")
    )
  })
  output$mtor_plot <- renderPlotly({
    dur <- 180; base <- input$seiz_history
    d1 <- run_sim(input$dose_VPA,0,0,0,0,0,input$pgp_init,
                  input$mTOR_val, 0.0, base, dur) %>%
      mutate(scenario="Without Everolimus")
    d2 <- run_sim(input$dose_VPA,0,0,0,0,0,input$pgp_init,
                  input$mTOR_val, input$ever_dose, base, dur) %>%
      mutate(scenario=paste0("Everolimus (dose=", input$ever_dose, ")"))
    df <- bind_rows(d1, d2)
    p <- ggplot(df, aes(x=day, y=SeizFreq_obs, color=scenario)) +
      geom_line(linewidth=1.2) +
      scale_color_manual(values=c("darkorange","steelblue")) +
      labs(title=paste0("TSC/mTOR: mTOR activity=", input$mTOR_val),
           x="Day", y="Seizure frequency (per 28 days)", color="") +
      theme_bw(base_size=11) + theme(legend.position="bottom")
    ggplotly(p) %>% layout(legend=list(orientation="h"))
  })
  output$safety_tbl <- renderDT({
    tbl <- data.frame(
      AED       = c("VPA","LEV","CBZ","LTG","PHT","GBP","PER","LCS"),
      `Hepatotoxicity`   = c("++","−","−","−","+","−","−","−"),
      `Teratogenicity`   = c("+++","±","++","+","++","−","−","−"),
      `Cognitive AE`     = c("++","±","++","±","++","−","±","−"),
      `Hyponatremia`     = c("−","−","+++","−","±","±","−","±"),
      `Bone density↓`    = c("+","−","++","−","++","−","−","−"),
      `DDI (inducer)`    = c("−","−","+++","−","+++","−","−","−"),
      `DDI (inhibitor)`  = c("UGT inhibit","−","−","−","CYP2C9","−","−","−"),
      `Monitoring`       = c("LFT,NH3","Renal Cr","CBC,Na,LFT","Rash,LFT","CBC,LFT","Renal","Psych","ECG"),
      stringsAsFactors=FALSE, check.names=FALSE
    )
    datatable(tbl, options=list(dom="t", scrollX=TRUE), rownames=FALSE) %>%
      formatStyle("Teratogenicity",
        backgroundColor=styleEqual(c("+++","++","+","±","−"),
          c("#FFAAAA","#FFCCAA","#FFEEDD","#FFFFF0","white")))
  })
}

shinyApp(ui=ui, server=server)
