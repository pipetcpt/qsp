################################################################################
# Schizophrenia QSP Interactive Dashboard
# 6 tabs: ① Patient Profile · ② PK Profiles · ③ D2/5-HT2A Occupancy
#         ④ Clinical Endpoints (PANSS) · ⑤ Scenario Comparison · ⑥ Biomarkers
################################################################################

library(shiny)
library(shinydashboard)
library(mrgsolve)
library(dplyr)
library(ggplot2)
library(tidyr)
library(plotly)
library(DT)

# ─────────────────────────────────────────────────────────────────────────────
# Inline mrgsolve model (shared with sch_mrgsolve_model.R)
# ─────────────────────────────────────────────────────────────────────────────
sch_model_code <- '
$PARAM
ka_HAL=0.80 F_HAL=0.65 CL_HAL=15.0 Vc_HAL=20.0 Vp_HAL=250.0 Qp_HAL=30.0 Kp_HAL=12.0
ka_RIS=1.00 F_RIS=0.74 CL_RIS=25.0 Vc_RIS=30.0 Vp_RIS=100.0 Qp_RIS=20.0 Kp_RIS=7.0
CL_PALI=7.5 Vc_PALI=50.0 fm_RIS=0.75
ka_CLZ=0.60 F_CLZ=0.55 CL_CLZ=30.0 Vc_CLZ=50.0 Kp_CLZ=6.0
ka_ARI=0.30 F_ARI=0.87 CL_ARI=3.6 Vc_ARI=245.0 CL_dARI=1.5 Vc_dARI=245.0 fm_ARI=0.40
D2tot=100.0 Kd_HAL=1.0 Kd_RIS=3.0 Kd_PALI=3.0 Kd_CLZ=160.0 Kd_ARI=0.34 Kd_dARI=1.0
HT2Atot=100.0 Kd_RIS_HT=0.16 Kd_CLZ_HT=5.3 Kd_ARI_HT=3.4 Kd_HAL_HT=53.0
DA_MESOLIM_0=1.0 DA_MESOCORT_0=1.0 DA_NIGROSTR_0=1.0 kout_DA=0.5
SCZ_amp=1.8 SCZ_sup=0.6 SCZ_nigrostr=1.0
PANSS_pos_0=35.0 PANSS_neg_0=28.0 PANSS_gen_0=50.0 kout_PANSS=0.03
Emax_pos=0.70 Emax_neg_SGA=0.40 Emax_neg_FGA=0.20 EC50_D2=70.0 EC50_HT2A=60.0
PRL_base=12.0 kout_PRL=0.5 Emax_PRL=3.0
EPS_thresh=80.0 EPS_slope=0.05
IL6_0=4.0 BDNF_0=22.0 OxidStress_0=1.5
PV_0=1.0 NMDAhypo_sev=0.4 KYNA_0=1.5

$CMT
GUT_HAL CENT_HAL PERI_HAL
GUT_RIS CENT_RIS PERI_RIS CENT_PALI
CENT_CLZ
GUT_ARI CENT_ARI CENT_dARI
DA_MESOLIM DA_MESOCORT DA_NIGROSTR
PRL_CMPT PV_ACT
PANSS_POS PANSS_NEG PANSS_GEN
BDNF_CMPT IL6_CMPT EPS_RISK

$INIT
GUT_HAL=0 CENT_HAL=0 PERI_HAL=0
GUT_RIS=0 CENT_RIS=0 PERI_RIS=0 CENT_PALI=0
CENT_CLZ=0
GUT_ARI=0 CENT_ARI=0 CENT_dARI=0
DA_MESOLIM=1.0 DA_MESOCORT=1.0 DA_NIGROSTR=1.0
PRL_CMPT=12.0 PV_ACT=1.0
PANSS_POS=35.0 PANSS_NEG=28.0 PANSS_GEN=50.0
BDNF_CMPT=22.0 IL6_CMPT=4.0 EPS_RISK=0.0

$ODE
dxdt_GUT_HAL  = -ka_HAL * GUT_HAL;
double k10_HAL=CL_HAL/Vc_HAL, k12_HAL=Qp_HAL/Vc_HAL, k21_HAL=Qp_HAL/Vp_HAL;
dxdt_CENT_HAL = ka_HAL*GUT_HAL - k10_HAL*CENT_HAL - k12_HAL*CENT_HAL + k21_HAL*PERI_HAL;
dxdt_PERI_HAL = k12_HAL*CENT_HAL - k21_HAL*PERI_HAL;

dxdt_GUT_RIS  = -ka_RIS * GUT_RIS;
double k10_RIS=CL_RIS/Vc_RIS, k12_RIS=Qp_RIS/Vc_RIS, k21_RIS=Qp_RIS/Vp_RIS;
double k_PALI=CL_PALI/Vc_PALI;
dxdt_CENT_RIS  = ka_RIS*GUT_RIS - k10_RIS*CENT_RIS - k12_RIS*CENT_RIS
               + k21_RIS*PERI_RIS - fm_RIS*(CL_RIS/Vc_RIS)*CENT_RIS;
dxdt_PERI_RIS  = k12_RIS*CENT_RIS - k21_RIS*PERI_RIS;
dxdt_CENT_PALI = fm_RIS*(CL_RIS/Vc_RIS)*CENT_RIS*(Vc_RIS/Vc_PALI) - k_PALI*CENT_PALI;

double k10_CLZ=CL_CLZ/Vc_CLZ;
dxdt_CENT_CLZ  = -k10_CLZ*CENT_CLZ;

dxdt_GUT_ARI  = -ka_ARI * GUT_ARI;
double k10_ARI=CL_ARI/Vc_ARI, k10_dARI=CL_dARI/Vc_dARI;
dxdt_CENT_ARI  = ka_ARI*GUT_ARI - k10_ARI*CENT_ARI - fm_ARI*k10_ARI*CENT_ARI;
dxdt_CENT_dARI = fm_ARI*k10_ARI*CENT_ARI*(Vc_ARI/Vc_dARI) - k10_dARI*CENT_dARI;

double Cp_HAL=CENT_HAL/Vc_HAL, Cb_HAL=Kp_HAL*Cp_HAL;
double Cp_RIS=CENT_RIS/Vc_RIS, Cb_RIS=Kp_RIS*Cp_RIS;
double Cp_PALI=CENT_PALI/Vc_PALI, Cb_PALI=Kp_RIS*Cp_PALI;
double Cp_CLZ=CENT_CLZ/Vc_CLZ, Cb_CLZ=Kp_CLZ*Cp_CLZ;
double Cp_ARI=CENT_ARI/Vc_ARI, Cb_ARI=15.0*Cp_ARI;
double Cp_dARI=CENT_dARI/Vc_dARI, Cb_dARI=15.0*Cp_dARI;

double D2occ_num=Cb_HAL/Kd_HAL+Cb_RIS/Kd_RIS+Cb_PALI/Kd_PALI+Cb_CLZ/Kd_CLZ+Cb_ARI/Kd_ARI+Cb_dARI/Kd_dARI;
double D2_occ_frac=D2occ_num/(1.0+D2occ_num);
double D2_occ_pct=100.0*D2_occ_frac;
double HT2A_num=Cb_RIS/Kd_RIS_HT+Cb_CLZ/Kd_CLZ_HT+Cb_ARI/Kd_ARI_HT+Cb_HAL/Kd_HAL_HT;
double HT2A_occ_frac=HT2A_num/(1.0+HT2A_num);
double HT2A_occ_pct=100.0*HT2A_occ_frac;

double DA_mesolim_scz=DA_MESOLIM_0*SCZ_amp;
double D2_eff_meso=D2_occ_frac*0.8+HT2A_occ_frac*0.1;
double kin_MESO=kout_DA*DA_mesolim_scz*(1.0-D2_eff_meso);
dxdt_DA_MESOLIM = kin_MESO-kout_DA*DA_MESOLIM;

double DA_mesocort_scz=DA_MESOCORT_0*SCZ_sup;
double HT2A_eff_meso=HT2A_occ_frac*0.6;
double kin_CORT=kout_DA*(DA_mesocort_scz+HT2A_eff_meso*(DA_MESOCORT_0-DA_mesocort_scz));
dxdt_DA_MESOCORT=kin_CORT-kout_DA*DA_MESOCORT;

double kin_NIGRO=kout_DA*DA_NIGROSTR_0;
double HT2A_nigro_rel=HT2A_occ_frac*0.4;
dxdt_DA_NIGROSTR=kin_NIGRO*(1.0-D2_occ_frac+HT2A_nigro_rel)-kout_DA*DA_NIGROSTR;

double PRL_Emax_eff=Emax_PRL*D2_occ_frac;
double kin_PRL=kout_PRL*PRL_base*(1.0+PRL_Emax_eff);
double ari_pn=(Cb_ARI/Kd_ARI+Cb_dARI/Kd_dARI)/(1.0+Cb_ARI/Kd_ARI+Cb_dARI/Kd_dARI)*0.5;
dxdt_PRL_CMPT=kin_PRL*(1.0-ari_pn)-kout_PRL*PRL_CMPT;

double PV_scz_level=PV_0*(1.0-NMDAhypo_sev);
double HT1A_restore=0.0;
if(Cb_ARI>0){HT1A_restore=0.15*(Cb_ARI/(Cb_ARI+5.1));}
dxdt_PV_ACT=0.1*(PV_scz_level+HT1A_restore*PV_0-PV_ACT);

double E_D2_pos=Emax_pos*pow(D2_occ_pct,2)/(pow(EC50_D2,2)+pow(D2_occ_pct,2));
double target_pos=PANSS_pos_0*(1.0-E_D2_pos);
dxdt_PANSS_POS=kout_PANSS*(target_pos-PANSS_POS);

double E_HT2A_neg=Emax_neg_SGA*HT2A_occ_pct/(EC50_HT2A+HT2A_occ_pct);
double E_D2_neg=Emax_neg_FGA*D2_occ_pct/(EC50_D2+D2_occ_pct);
double target_neg=PANSS_neg_0*(1.0-E_HT2A_neg-E_D2_neg);
dxdt_PANSS_NEG=kout_PANSS*(target_neg-PANSS_NEG);

double E_gen=0.5*E_D2_pos+0.5*E_HT2A_neg;
double target_gen=PANSS_gen_0*(1.0-E_gen*0.6);
dxdt_PANSS_GEN=kout_PANSS*(target_gen-PANSS_GEN);

double BDNF_restore=HT2A_occ_frac*0.3*(28.0-BDNF_0);
dxdt_BDNF_CMPT=0.05*(BDNF_0+BDNF_restore-BDNF_CMPT);
dxdt_IL6_CMPT=0.1*(IL6_0*(1.0-D2_occ_frac*0.2)-IL6_CMPT);

double EPS_excess=(D2_occ_pct>EPS_thresh)?EPS_slope*(D2_occ_pct-EPS_thresh):0.0;
double EPS_SGA_benefit=HT2A_occ_frac*0.5*EPS_excess;
dxdt_EPS_RISK=0.5*(EPS_excess-EPS_SGA_benefit-EPS_RISK);

$TABLE
double CP_HAL=CENT_HAL/Vc_HAL;
double CP_RIS=CENT_RIS/Vc_RIS;
double CP_PALI=CENT_PALI/Vc_PALI;
double CP_CLZ=CENT_CLZ/Vc_CLZ;
double CP_ARI=CENT_ARI/Vc_ARI;
double CP_dARI=CENT_dARI/Vc_dARI;
double D2_OCC_PCT = 100.0*(CENT_HAL/Vc_HAL*Kp_HAL/Kd_HAL+CENT_RIS/Vc_RIS*Kp_RIS/Kd_RIS+CENT_PALI/Vc_PALI*Kp_RIS/Kd_PALI+CENT_CLZ/Vc_CLZ*Kp_CLZ/Kd_CLZ+CENT_ARI/Vc_ARI*15.0/Kd_ARI+CENT_dARI/Vc_dARI*15.0/Kd_dARI)/(1.0+CENT_HAL/Vc_HAL*Kp_HAL/Kd_HAL+CENT_RIS/Vc_RIS*Kp_RIS/Kd_RIS+CENT_PALI/Vc_PALI*Kp_RIS/Kd_PALI+CENT_CLZ/Vc_CLZ*Kp_CLZ/Kd_CLZ+CENT_ARI/Vc_ARI*15.0/Kd_ARI+CENT_dARI/Vc_dARI*15.0/Kd_dARI);
double HT2A_OCC_PCT = 100.0*(CENT_RIS/Vc_RIS*Kp_RIS/Kd_RIS_HT+CENT_CLZ/Vc_CLZ*Kp_CLZ/Kd_CLZ_HT+CENT_ARI/Vc_ARI*15.0/Kd_ARI_HT+CENT_HAL/Vc_HAL*Kp_HAL/Kd_HAL_HT)/(1.0+CENT_RIS/Vc_RIS*Kp_RIS/Kd_RIS_HT+CENT_CLZ/Vc_CLZ*Kp_CLZ/Kd_CLZ_HT+CENT_ARI/Vc_ARI*15.0/Kd_ARI_HT+CENT_HAL/Vc_HAL*Kp_HAL/Kd_HAL_HT);
double PANSS_TOTAL = PANSS_POS+PANSS_NEG+PANSS_GEN;

$CAPTURE
CP_HAL CP_RIS CP_PALI CP_CLZ CP_ARI CP_dARI
D2_OCC_PCT HT2A_OCC_PCT
DA_MESOLIM DA_MESOCORT DA_NIGROSTR
PRL_CMPT PV_ACT
PANSS_POS PANSS_NEG PANSS_GEN PANSS_TOTAL
BDNF_CMPT IL6_CMPT EPS_RISK
'

mod <- mread_cache("sch_shiny", inline = sch_model_code)

# ─────────────────────────────────────────────────────────────────────────────
# SIMULATION HELPER
# ─────────────────────────────────────────────────────────────────────────────
simulate_ap <- function(drug, dose, duration_days, freq_h = 24,
                        panss_pos_0 = 35, panss_neg_0 = 28, panss_gen_0 = 50,
                        scz_severity = "moderate") {
  scz_amp  <- switch(scz_severity, mild = 1.4, moderate = 1.8, severe = 2.2)
  scz_sup  <- switch(scz_severity, mild = 0.8, moderate = 0.6, severe = 0.4)
  nmda_sev <- switch(scz_severity, mild = 0.25, moderate = 0.40, severe = 0.55)

  cmt_map <- c(haloperidol = "GUT_HAL", risperidone = "GUT_RIS",
               clozapine = "CENT_CLZ", aripiprazole = "GUT_ARI")
  cmt <- cmt_map[[drug]]
  dose_times <- seq(0, (duration_days * 24 - freq_h), by = freq_h)
  ev_dose <- ev(time = dose_times, amt = dose, cmt = cmt)

  mod %>%
    param(PANSS_pos_0 = panss_pos_0, PANSS_neg_0 = panss_neg_0,
          PANSS_gen_0 = panss_gen_0,
          SCZ_amp = scz_amp, SCZ_sup = scz_sup, NMDAhypo_sev = nmda_sev) %>%
    mrgsim(events = ev_dose, end = duration_days * 24, delta = 2) %>%
    as_tibble() %>%
    mutate(day = time / 24)
}

# Color palette
pal <- c("#e74c3c","#3498db","#2ecc71","#9b59b6","#f39c12","#1abc9c","#e67e22")

# ─────────────────────────────────────────────────────────────────────────────
# UI
# ─────────────────────────────────────────────────────────────────────────────
ui <- dashboardPage(
  skin = "purple",
  dashboardHeader(title = "Schizophrenia QSP"),
  dashboardSidebar(
    sidebarMenu(
      menuItem("① Patient Profile",    tabName = "profile",   icon = icon("user")),
      menuItem("② PK Profiles",        tabName = "pk",        icon = icon("flask")),
      menuItem("③ D2/5-HT2A Occupancy",tabName = "occupancy", icon = icon("chart-line")),
      menuItem("④ PANSS Endpoints",    tabName = "panss",     icon = icon("brain")),
      menuItem("⑤ Scenario Comparison",tabName = "scenarios", icon = icon("balance-scale")),
      menuItem("⑥ Biomarkers",         tabName = "biomarkers",icon = icon("vial"))
    )
  ),
  dashboardBody(
    tabItems(

      # ── TAB 1: PATIENT PROFILE ─────────────────────────────────────────────
      tabItem("profile",
        fluidRow(
          box(title = "Patient Demographics", status = "primary", solidHeader = TRUE,
              width = 4,
              numericInput("age",    "Age (years)",     value = 28, min = 18, max = 80),
              selectInput("sex",    "Sex",             choices = c("Male","Female")),
              numericInput("weight","Body Weight (kg)", value = 70, min = 40, max = 150),
              selectInput("scz_severity","SCZ Severity",
                          choices = c("Mild","Moderate","Severe"), selected = "Moderate")
          ),
          box(title = "Disease Status (Baseline PANSS)", status = "warning", solidHeader = TRUE,
              width = 4,
              sliderInput("panss_pos_0","PANSS Positive (7-49)", 7, 49, 35),
              sliderInput("panss_neg_0","PANSS Negative (7-49)", 7, 49, 28),
              sliderInput("panss_gen_0","PANSS General (16-112)",16, 112, 50),
              verbatimTextOutput("panss_baseline_total")
          ),
          box(title = "Treatment Selection", status = "success", solidHeader = TRUE,
              width = 4,
              selectInput("drug","Antipsychotic Drug",
                          choices = c("Haloperidol (FGA)" = "haloperidol",
                                      "Risperidone (SGA)" = "risperidone",
                                      "Clozapine (TRS-SGA)" = "clozapine",
                                      "Aripiprazole (Partial D2)" = "aripiprazole")),
              numericInput("dose","Dose (mg/day)", value = 10, min = 1, max = 600),
              numericInput("duration","Treatment Duration (days)", value = 90, min = 7, max = 365),
              selectInput("freq","Dosing Frequency",
                          choices = c("Once daily (QD)" = "24",
                                      "Twice daily (BID)" = "12")),
              actionButton("run_sim","▶ Run Simulation", class = "btn-success btn-lg",
                           width = "100%")
          )
        ),
        fluidRow(
          box(title = "Disease Mechanism Summary", status = "info", solidHeader = TRUE,
              width = 12,
              div(style = "font-size:13px;",
                  h4("Schizophrenia Pathophysiology"),
                  tags$ul(
                    tags$li(strong("Dopamine Hypothesis:"),
                            " Mesolimbic DA hyperactivity → positive symptoms;
                              Mesocortical DA hypoactivity → negative + cognitive symptoms"),
                    tags$li(strong("Glutamate/NMDA Hypofunction:"),
                            " PV interneuron dysfunction → cortical disinhibition;
                              KYNA elevation antagonizes NMDA receptors"),
                    tags$li(strong("GABAergic Deficit:"),
                            " GAD67 ↓25-50%; PV+ fast-spiking interneurons impaired →
                              loss of gamma oscillation → working memory deficit"),
                    tags$li(strong("Neuroinflammation:"),
                            " IL-6, IL-1β, TNF-α elevated; C4A complement → synaptic pruning;
                              BDNF ↓ in SCZ"),
                    tags$li(strong("Treatment Strategy:"),
                            " FGA: D2 blockade (65-80% occupancy optimal);
                              SGA: D2 + 5-HT2A block → improved negative Sx + ↓ EPS;
                              ARI: Partial D2 agonist → dopamine stabilization")
                  )
              )
          )
        )
      ),

      # ── TAB 2: PK PROFILES ────────────────────────────────────────────────
      tabItem("pk",
        fluidRow(
          box(title = "Plasma Concentration–Time Profile", status = "primary",
              solidHeader = TRUE, width = 8,
              plotlyOutput("pk_plot", height = 450)
          ),
          box(title = "PK Parameters", status = "info", solidHeader = TRUE, width = 4,
              tableOutput("pk_params_table"),
              hr(),
              h5("Brain Concentration"),
              plotlyOutput("brain_conc_plot", height = 200)
          )
        ),
        fluidRow(
          box(title = "Multi-Dose Steady-State (Days 1–14)", status = "success",
              solidHeader = TRUE, width = 12,
              plotlyOutput("pk_ss_plot", height = 350)
          )
        )
      ),

      # ── TAB 3: RECEPTOR OCCUPANCY ─────────────────────────────────────────
      tabItem("occupancy",
        fluidRow(
          box(title = "D2 Receptor Occupancy (%) vs. Time", status = "primary",
              solidHeader = TRUE, width = 6,
              plotlyOutput("d2_occ_plot", height = 400),
              hr(),
              p(em("Therapeutic window: 65–80% D2 occupancy (Kapur & Seeman 2000)"))
          ),
          box(title = "5-HT2A Receptor Occupancy (%) vs. Time", status = "warning",
              solidHeader = TRUE, width = 6,
              plotlyOutput("ht2a_occ_plot", height = 400),
              hr(),
              p(em("SGA targets >80% 5-HT2A occupancy for improved negative Sx"))
          )
        ),
        fluidRow(
          box(title = "D2 vs 5-HT2A Occupancy — Drug Fingerprint", status = "info",
              solidHeader = TRUE, width = 12,
              plotlyOutput("occ_scatter", height = 350)
          )
        )
      ),

      # ── TAB 4: PANSS CLINICAL ENDPOINTS ───────────────────────────────────
      tabItem("panss",
        fluidRow(
          box(title = "PANSS Total Score", status = "danger", solidHeader = TRUE,
              width = 6,
              plotlyOutput("panss_total_plot", height = 350),
              verbatimTextOutput("panss_response_text")
          ),
          box(title = "PANSS Subscale Scores", status = "warning", solidHeader = TRUE,
              width = 6,
              plotlyOutput("panss_subscale_plot", height = 350)
          )
        ),
        fluidRow(
          box(title = "Dopamine Pathway Dynamics", status = "primary",
              solidHeader = TRUE, width = 6,
              plotlyOutput("da_pathway_plot", height = 350)
          ),
          box(title = "Clinical Response Summary", status = "success",
              solidHeader = TRUE, width = 6,
              plotlyOutput("response_gauge", height = 200),
              tableOutput("response_table")
          )
        )
      ),

      # ── TAB 5: SCENARIO COMPARISON ────────────────────────────────────────
      tabItem("scenarios",
        fluidRow(
          box(title = "Compare Treatment Scenarios", status = "primary",
              solidHeader = TRUE, width = 4,
              checkboxGroupInput("scen_drugs","Drugs to Compare",
                choices = c("Haloperidol 5mg"  = "hal5",
                            "Haloperidol 10mg" = "hal10",
                            "Risperidone 2mg"  = "ris2",
                            "Risperidone 4mg"  = "ris4",
                            "Clozapine 300mg"  = "clz300",
                            "Aripiprazole 15mg"= "ari15",
                            "Untreated"        = "none"),
                selected = c("hal10","ris4","clz300","ari15","none")),
              numericInput("scen_duration","Duration (days)", value = 90, min = 14, max = 365),
              actionButton("run_scen","▶ Compare All", class = "btn-primary", width = "100%")
          ),
          box(title = "PANSS Total Comparison", status = "danger",
              solidHeader = TRUE, width = 8,
              plotlyOutput("scen_panss_plot", height = 380)
          )
        ),
        fluidRow(
          box(title = "D2 Occupancy Comparison", status = "info",
              solidHeader = TRUE, width = 6,
              plotlyOutput("scen_d2_plot", height = 350)
          ),
          box(title = "EPS Risk Comparison", status = "warning",
              solidHeader = TRUE, width = 6,
              plotlyOutput("scen_eps_plot", height = 350)
          )
        ),
        fluidRow(
          box(title = "Summary Table (Day 90)", status = "success",
              solidHeader = TRUE, width = 12,
              DTOutput("scen_summary_table")
          )
        )
      ),

      # ── TAB 6: BIOMARKERS ────────────────────────────────────────────────
      tabItem("biomarkers",
        fluidRow(
          box(title = "Prolactin (ng/mL)", status = "warning", solidHeader = TRUE,
              width = 6,
              plotlyOutput("prl_plot", height = 350),
              p(em("ULN: ~25 ng/mL women / ~20 ng/mL men. RIS/HAL → hyperprolactinemia"))
          ),
          box(title = "EPS Risk Index", status = "danger", solidHeader = TRUE,
              width = 6,
              plotlyOutput("eps_plot", height = 350),
              p(em("EPS: D2>80% nigrostriatal; SGA 5-HT2A block mitigates risk"))
          )
        ),
        fluidRow(
          box(title = "BDNF Level (ng/mL)", status = "info", solidHeader = TRUE,
              width = 6,
              plotlyOutput("bdnf_plot", height = 350),
              p(em("Baseline BDNF ↓ in SCZ (~22 ng/mL vs ~28 ng/mL healthy). SGA partially restores."))
          ),
          box(title = "IL-6 (pg/mL) — Neuroinflammation", status = "primary",
              solidHeader = TRUE, width = 6,
              plotlyOutput("il6_plot", height = 350),
              p(em("IL-6 elevated in SCZ; antipsychotics have modest anti-inflammatory effect."))
          )
        ),
        fluidRow(
          box(title = "PV Interneuron Activity (normalized)", status = "success",
              solidHeader = TRUE, width = 12,
              plotlyOutput("pv_plot", height = 300),
              p(em("PV interneuron deficit (↓ 25-40%) underlies gamma oscillation loss and WM deficits in SCZ.
                    Aripiprazole 5-HT1A partial agonism provides modest PV restoration."))
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

  output$panss_baseline_total <- renderText({
    total <- input$panss_pos_0 + input$panss_neg_0 + input$panss_gen_0
    severity <- if (total < 60) "Mild" else if (total < 90) "Moderate" else "Severe"
    paste0("PANSS Total: ", total, " (", severity, ")")
  })

  # Reactive simulation result
  sim_result <- eventReactive(input$run_sim, {
    req(input$drug, input$dose, input$duration)
    withProgress(message = "Running simulation...", {
      simulate_ap(
        drug = input$drug,
        dose = input$dose,
        duration_days = input$duration,
        freq_h = as.numeric(input$freq),
        panss_pos_0 = input$panss_pos_0,
        panss_neg_0 = input$panss_neg_0,
        panss_gen_0 = input$panss_gen_0,
        scz_severity = tolower(input$scz_severity)
      )
    })
  }, ignoreNULL = FALSE)

  # ── TAB 2: PK ─────────────────────────────────────────────────────────────
  pk_col_map <- c(haloperidol = "CP_HAL", risperidone = "CP_RIS",
                  clozapine = "CP_CLZ", aripiprazole = "CP_ARI")

  output$pk_plot <- renderPlotly({
    df <- sim_result()
    req(df)
    col <- pk_col_map[[input$drug]]
    p <- ggplot(df, aes(day, .data[[col]])) +
      geom_line(color = "#3498db", size = 1.2) +
      labs(title = paste("Plasma Concentration:", tools::toTitleCase(input$drug)),
           x = "Day", y = "Concentration (ng/mL)") +
      theme_minimal()
    ggplotly(p)
  })

  output$brain_conc_plot <- renderPlotly({
    df <- sim_result()
    req(df)
    col <- pk_col_map[[input$drug]]
    kp <- switch(input$drug, haloperidol=12, risperidone=7, clozapine=6, aripiprazole=15)
    df$brain <- df[[col]] * kp
    p <- ggplot(df, aes(day, brain)) +
      geom_line(color = "#9b59b6", size = 1.2) +
      labs(x = "Day", y = "Brain Conc (ng/mL)") +
      theme_minimal()
    ggplotly(p)
  })

  output$pk_ss_plot <- renderPlotly({
    df <- sim_result() %>% filter(day <= 14)
    col <- pk_col_map[[input$drug]]
    p <- ggplot(df, aes(day, .data[[col]])) +
      geom_line(color = "#e74c3c", size = 1.2) +
      labs(title = "Approach to Steady-State (First 14 Days)",
           x = "Day", y = "Plasma Concentration (ng/mL)") +
      theme_minimal()
    ggplotly(p)
  })

  output$pk_params_table <- renderTable({
    data.frame(
      Parameter = c("Bioavailability (F)","CL (L/h)","Vc (L)","Half-life (h)","Brain Kp"),
      Value = switch(input$drug,
        haloperidol  = c("65%","15","20","18-24","12"),
        risperidone  = c("74%","25","30","3/21*","7"),
        clozapine    = c("55%","30","50","12","6"),
        aripiprazole = c("87%","3.6","245","75","~15")
      )
    )
  })

  # ── TAB 3: RECEPTOR OCCUPANCY ─────────────────────────────────────────────
  output$d2_occ_plot <- renderPlotly({
    df <- sim_result()
    req(df)
    p <- ggplot(df, aes(day, D2_OCC_PCT)) +
      geom_line(color = "#e74c3c", size = 1.2) +
      geom_hline(yintercept = 65, linetype = "dashed", color = "green4") +
      geom_hline(yintercept = 80, linetype = "dashed", color = "red3") +
      annotate("rect", xmin=-Inf, xmax=Inf, ymin=65, ymax=80,
               alpha=0.1, fill="green") +
      labs(title = "D2 Receptor Occupancy (%)",
           x = "Day", y = "D2 Occupancy (%)") +
      ylim(0, 100) + theme_minimal()
    ggplotly(p)
  })

  output$ht2a_occ_plot <- renderPlotly({
    df <- sim_result()
    req(df)
    p <- ggplot(df, aes(day, HT2A_OCC_PCT)) +
      geom_line(color = "#9b59b6", size = 1.2) +
      geom_hline(yintercept = 80, linetype = "dashed", color = "purple") +
      labs(title = "5-HT2A Receptor Occupancy (%)",
           x = "Day", y = "5-HT2A Occupancy (%)") +
      ylim(0, 100) + theme_minimal()
    ggplotly(p)
  })

  output$occ_scatter <- renderPlotly({
    df <- sim_result() %>%
      filter(day %in% c(1,7,14,30,60,90)) %>%
      group_by(day) %>% slice(1)
    p <- ggplot(df, aes(D2_OCC_PCT, HT2A_OCC_PCT, color = factor(day), label = paste("Day",day))) +
      geom_point(size = 4) +
      geom_path(arrow = arrow(length = unit(0.2,"cm"))) +
      geom_vline(xintercept = c(65,80), linetype = "dashed", color = c("green4","red")) +
      geom_hline(yintercept = 80, linetype = "dashed", color = "purple") +
      labs(title = "D2 vs 5-HT2A Occupancy Over Time",
           x = "D2 Occupancy (%)", y = "5-HT2A Occupancy (%)",
           color = "Day") +
      theme_minimal()
    ggplotly(p)
  })

  # ── TAB 4: PANSS ──────────────────────────────────────────────────────────
  output$panss_total_plot <- renderPlotly({
    df <- sim_result()
    req(df)
    p <- ggplot(df, aes(day, PANSS_TOTAL)) +
      geom_line(color = "#e74c3c", size = 1.3) +
      geom_hline(yintercept = input$panss_pos_0 + input$panss_neg_0 + input$panss_gen_0,
                 linetype = "dashed", color = "gray") +
      labs(title = "PANSS Total Score Over Time",
           x = "Day", y = "PANSS Total") +
      theme_minimal()
    ggplotly(p)
  })

  output$panss_subscale_plot <- renderPlotly({
    df <- sim_result()
    req(df)
    df_long <- df %>%
      select(day, PANSS_POS, PANSS_NEG, PANSS_GEN) %>%
      pivot_longer(-day, names_to = "Subscale", values_to = "Score")
    p <- ggplot(df_long, aes(day, Score, color = Subscale)) +
      geom_line(size = 1.1) +
      scale_color_manual(values = c("#e74c3c","#3498db","#2ecc71"),
                         labels = c("General","Negative","Positive")) +
      labs(title = "PANSS Subscale Scores", x = "Day", y = "Score") +
      theme_minimal()
    ggplotly(p)
  })

  output$panss_response_text <- renderText({
    df <- sim_result()
    req(df)
    baseline <- input$panss_pos_0 + input$panss_neg_0 + input$panss_gen_0
    final_row <- df %>% filter(day == max(day)) %>% slice(1)
    final_total <- final_row$PANSS_TOTAL
    reduction_pct <- (baseline - final_total) / baseline * 100
    resp <- if (reduction_pct >= 20) "Responder (≥20% reduction)"
            else if (reduction_pct >= 10) "Partial Responder"
            else "Non-Responder"
    paste0("PANSS Reduction: ", round(reduction_pct, 1), "% — ", resp)
  })

  output$da_pathway_plot <- renderPlotly({
    df <- sim_result()
    req(df)
    df_long <- df %>%
      select(day, DA_MESOLIM, DA_MESOCORT, DA_NIGROSTR) %>%
      pivot_longer(-day, names_to = "Pathway", values_to = "Activity")
    p <- ggplot(df_long, aes(day, Activity, color = Pathway)) +
      geom_line(size = 1.1) +
      geom_hline(yintercept = 1.0, linetype = "dashed", color = "gray40") +
      scale_color_manual(values = c("#e74c3c","#3498db","#2ecc71"),
                         labels = c("Mesolimbic","Mesocortical","Nigrostriatal")) +
      labs(title = "Dopamine Pathway Activity", x = "Day",
           y = "DA Activity (normalized)") +
      theme_minimal()
    ggplotly(p)
  })

  output$response_table <- renderTable({
    df <- sim_result()
    req(df)
    days_pts <- c(7, 14, 30, 90, 180)
    days_pts <- days_pts[days_pts <= max(df$day)]
    baseline <- input$panss_pos_0 + input$panss_neg_0 + input$panss_gen_0
    df %>%
      filter(day %in% sapply(days_pts, function(d) which.min(abs(df$day - d))[1])) %>%
      group_by(day) %>% slice(1) %>%
      mutate(
        Day = round(day),
        `PANSS Total` = round(PANSS_TOTAL, 1),
        `% Reduction` = round((baseline - PANSS_TOTAL)/baseline*100, 1),
        `D2 Occ (%)` = round(D2_OCC_PCT, 1)
      ) %>%
      select(Day, `PANSS Total`, `% Reduction`, `D2 Occ (%)`) %>%
      ungroup()
  })

  # ── TAB 5: SCENARIO COMPARISON ────────────────────────────────────────────
  scen_results <- eventReactive(input$run_scen, {
    drug_map <- list(
      hal5   = list(drug="haloperidol",  dose=5),
      hal10  = list(drug="haloperidol",  dose=10),
      ris2   = list(drug="risperidone",  dose=2),
      ris4   = list(drug="risperidone",  dose=4),
      clz300 = list(drug="clozapine",    dose=300),
      ari15  = list(drug="aripiprazole", dose=15)
    )
    dur <- input$scen_duration
    baseline_pos <- input$panss_pos_0
    baseline_neg <- input$panss_neg_0
    baseline_gen <- input$panss_gen_0
    scz_sev <- tolower(input$scz_severity)

    results <- list()
    withProgress(message = "Running scenario comparison...", {
      for (key in input$scen_drugs) {
        if (key == "none") {
          ev_none <- ev(time = 0, amt = 0, cmt = "GUT_HAL")
          df_none <- mod %>%
            param(PANSS_pos_0 = baseline_pos, PANSS_neg_0 = baseline_neg,
                  PANSS_gen_0 = baseline_gen) %>%
            mrgsim(events = ev_none, end = dur * 24, delta = 2) %>%
            as_tibble() %>%
            mutate(day = time/24, scenario = "Untreated")
          results[["none"]] <- df_none
        } else if (key %in% names(drug_map)) {
          dm <- drug_map[[key]]
          df_s <- simulate_ap(dm$drug, dm$dose, dur, 24,
                              baseline_pos, baseline_neg, baseline_gen, scz_sev)
          label_map <- c(
            hal5="HAL 5mg", hal10="HAL 10mg", ris2="RIS 2mg",
            ris4="RIS 4mg", clz300="CLZ 300mg", ari15="ARI 15mg"
          )
          df_s$scenario <- label_map[[key]]
          results[[key]] <- df_s
        }
      }
    })
    bind_rows(results)
  })

  output$scen_panss_plot <- renderPlotly({
    df <- scen_results()
    req(df)
    p <- ggplot(df, aes(day, PANSS_TOTAL, color = scenario)) +
      geom_line(size = 1.1) +
      scale_color_manual(values = pal) +
      labs(title = "PANSS Total Score Comparison",
           x = "Day", y = "PANSS Total", color = "Scenario") +
      theme_minimal()
    ggplotly(p)
  })

  output$scen_d2_plot <- renderPlotly({
    df <- scen_results()
    req(df)
    p <- ggplot(df, aes(day, D2_OCC_PCT, color = scenario)) +
      geom_line(size = 1.1) +
      geom_hline(yintercept = c(65,80), linetype="dashed",
                 color=c("green4","red"), alpha=0.8) +
      scale_color_manual(values = pal) +
      labs(title = "D2 Occupancy Comparison",
           x = "Day", y = "D2 Occupancy (%)", color = "Scenario") +
      ylim(0, 100) + theme_minimal()
    ggplotly(p)
  })

  output$scen_eps_plot <- renderPlotly({
    df <- scen_results()
    req(df)
    p <- ggplot(df, aes(day, EPS_RISK, color = scenario)) +
      geom_line(size = 1.1) +
      scale_color_manual(values = pal) +
      labs(title = "EPS Risk Index Comparison",
           x = "Day", y = "EPS Risk Index", color = "Scenario") +
      theme_minimal()
    ggplotly(p)
  })

  output$scen_summary_table <- renderDT({
    df <- scen_results()
    req(df)
    max_day <- max(df$day)
    tbl_day <- min(90, max_day)
    baseline <- input$panss_pos_0 + input$panss_neg_0 + input$panss_gen_0
    df %>%
      filter(abs(day - tbl_day) < 1) %>%
      group_by(scenario) %>%
      slice(1) %>%
      mutate(
        `PANSS Total` = round(PANSS_TOTAL, 1),
        `PANSS % Reduction` = round((baseline-PANSS_TOTAL)/baseline*100, 1),
        `D2 Occ (%)` = round(D2_OCC_PCT, 1),
        `HT2A Occ (%)` = round(HT2A_OCC_PCT, 1),
        `EPS Risk` = round(EPS_RISK, 3),
        `Prolactin` = round(PRL_CMPT, 1)
      ) %>%
      select(scenario, `PANSS Total`, `PANSS % Reduction`,
             `D2 Occ (%)`, `HT2A Occ (%)`, `EPS Risk`, `Prolactin`) %>%
      ungroup() %>%
      rename(Scenario = scenario) %>%
      as.data.frame()
  }, options = list(pageLength = 10))

  # ── TAB 6: BIOMARKERS ─────────────────────────────────────────────────────
  output$prl_plot <- renderPlotly({
    df <- sim_result(); req(df)
    p <- ggplot(df, aes(day, PRL_CMPT)) +
      geom_line(color="#f39c12", size=1.2) +
      geom_hline(yintercept=25, linetype="dashed", color="red") +
      labs(title="Prolactin Over Time", x="Day", y="Prolactin (ng/mL)") +
      theme_minimal()
    ggplotly(p)
  })

  output$eps_plot <- renderPlotly({
    df <- sim_result(); req(df)
    p <- ggplot(df, aes(day, EPS_RISK)) +
      geom_line(color="#e74c3c", size=1.2) +
      labs(title="EPS Risk Index", x="Day", y="EPS Risk") +
      theme_minimal()
    ggplotly(p)
  })

  output$bdnf_plot <- renderPlotly({
    df <- sim_result(); req(df)
    p <- ggplot(df, aes(day, BDNF_CMPT)) +
      geom_line(color="#3498db", size=1.2) +
      geom_hline(yintercept=28, linetype="dashed", color="green4") +
      annotate("text", x=max(df$day)*0.8, y=29,
               label="Healthy ~28 ng/mL", size=3, color="green4") +
      labs(title="BDNF Level", x="Day", y="BDNF (ng/mL)") +
      theme_minimal()
    ggplotly(p)
  })

  output$il6_plot <- renderPlotly({
    df <- sim_result(); req(df)
    p <- ggplot(df, aes(day, IL6_CMPT)) +
      geom_line(color="#9b59b6", size=1.2) +
      geom_hline(yintercept=2, linetype="dashed", color="green4") +
      annotate("text", x=max(df$day)*0.8, y=2.2,
               label="Healthy ~2 pg/mL", size=3, color="green4") +
      labs(title="IL-6 (Neuroinflammation Marker)",
           x="Day", y="IL-6 (pg/mL)") +
      theme_minimal()
    ggplotly(p)
  })

  output$pv_plot <- renderPlotly({
    df <- sim_result(); req(df)
    p <- ggplot(df, aes(day, PV_ACT)) +
      geom_line(color="#2ecc71", size=1.2) +
      geom_hline(yintercept=1.0, linetype="dashed", color="gray40") +
      annotate("text", x=5, y=1.02, label="Normal = 1.0",
               size=3, color="gray40", hjust=0) +
      labs(title="PV Interneuron Activity (normalized)",
           x="Day", y="PV Activity") +
      ylim(0, 1.2) +
      theme_minimal()
    ggplotly(p)
  })
}

shinyApp(ui, server)
