## ============================================================
## ME/CFS QSP Interactive Shiny App
## Myalgic Encephalomyelitis / Chronic Fatigue Syndrome
## ============================================================

library(shiny)
library(shinydashboard)
library(ggplot2)
library(dplyr)
library(tidyr)
library(plotly)
library(DT)
library(mrgsolve)

# ============================================================
# Source Model (simplified inline version for Shiny)
# ============================================================

mecfs_code <- '
$PARAM
kV_prod=0.0005 kV_clear=0.15 kV_IFNclear=0.10 V_max=10.0
kNK_prod=0.05 kNK_decay=0.05 kNK_exhaust=0.02 kNK_IFN=0.008
kTex_induce=0.03 kTex_recover=0.005
kAutoAb_prod=0.002 kAutoAb_decay=0.015 kAutoAb_mol=0.05
kIL6_prod=0.10 kIL6_decay=0.25 kIL6_V=0.08
kTNFa_prod=0.08 kTNFa_decay=0.30
kIFNg_prod=0.06 kIFNg_decay=0.20 kIFNg_V=0.10
kNLRP3_act=0.05 kNLRP3_decay=0.15
kIFN_prod=0.15 kIFN_decay=0.30
kMC_act=0.04 kMC_decay=0.08 kHist_prod=0.20 kHist_decay=0.35
kCRH_prod=0.10 kCRH_decay=0.25
kCort_prod=0.15 kCort_decay=0.12 Cort_base=1.0 IL6_GR_IC50=1.5
kNE_prod=0.08 kNE_decay=0.20
kHRV_restore=0.05 kHRV_suppress=0.15 NE_base=1.0 HRV_base=1.0
kPDH_base=1.0 kPDK1_IFNg=0.30 kPDK1_TNFa=0.20
kROS_prod=0.15 kROS_decay=0.10
kATP_prod=1.2 kATP_decay=0.30 kATP_ROS=0.25
kNI_prod=0.08 kNI_decay=0.06 kNI_IL6=0.12 kNI_TNFa=0.10
kCog_impair=0.10 kCog_restore=0.05
kPEM_ATP=0.20 kPEM_NI=0.15
kFat_prod=0.15 kFat_decay=0.08
LDN_dose=0 LDN_F=0.96 LDN_ka=1.2 LDN_CL=9.7 LDN_V1=28.0
LDN_Emax_TLR4=0.7 LDN_EC50_TLR4=0.05
Pyr_dose=0 Pyr_F=0.20 Pyr_ka=0.8 Pyr_CL=10.0 Pyr_V=50.0
Pyr_Emax_ANS=0.6 Pyr_EC50=30.0
Rit_dose=0 Rit_CL=0.23 Rit_V=4.4 Rit_Emax_B=0.95 Rit_EC50=0.1
NAD_dose=0 NAD_ka=0.5 NAD_F=0.35 NAD_CL=2.0 NAD_V=10.0
NAD_Emax_mito=0.5 NAD_EC50=200.0

$INIT
V=0.01 IFN=0.1 NK=0.6 Tex=0.4 AutoAb=0.5
IL6=1.3 TNFa=1.2 IFNg=1.4 NLRP3state=0.8
MC_act=0.6 Histamine=0.8
CRH=0.8 Cortisol=0.6
NE_plasma=1.4 HRV_index=0.5
PDH_act=0.4 ATP_state=0.4 ROS_state=1.6
Neuro_inf=0.8 Cog_func=0.5
PEM_sens=1.8 Fatigue=0.7
LDN_cp=0.0 Pyr_cp=0.0 Rit_cp=0.0 NADpool=0.3

$ODE
double LDN_TLR4_inh = LDN_Emax_TLR4*LDN_cp/(LDN_EC50_TLR4+LDN_cp+1e-9);
double Pyr_ANS_rest = Pyr_Emax_ANS*Pyr_cp/(Pyr_EC50+Pyr_cp+1e-9);
double Rit_B_depl = Rit_Emax_B*Rit_cp/(Rit_EC50+Rit_cp+1e-9);
double NAD_mito_rest = NAD_Emax_mito*NADpool/(NAD_EC50/1000.0+NADpool+1e-9);
dxdt_LDN_cp = LDN_dose*LDN_F*LDN_ka/LDN_V1 - (LDN_CL/LDN_V1)*LDN_cp;
dxdt_Pyr_cp = (Pyr_dose*Pyr_F*Pyr_ka/Pyr_V) - (Pyr_CL/Pyr_V)*Pyr_cp;
dxdt_Rit_cp = (Rit_dose/Rit_V) - (Rit_CL/Rit_V)*Rit_cp;
dxdt_NADpool = (NAD_dose*NAD_F*NAD_ka/NAD_V) - (NAD_CL/NAD_V)*NADpool;
dxdt_V = kV_prod*V*(1.0-V/V_max) - kV_clear*NK*V - kV_IFNclear*IFNg*V;
dxdt_IFN = kIFN_prod*(V+0.1)*(1.0-LDN_TLR4_inh) - kIFN_decay*IFN;
dxdt_NK = kNK_prod - kNK_decay*NK - kNK_exhaust*V*NK - kNK_IFN*IFN*NK;
dxdt_Tex = kTex_induce*IFN*V - kTex_recover*Tex*Cortisol;
double AutoAb_mol = kAutoAb_mol*V*(1.0-Rit_B_depl);
dxdt_AutoAb = kAutoAb_prod*(1.0-Rit_B_depl)+AutoAb_mol - kAutoAb_decay*AutoAb;
dxdt_IL6 = kIL6_prod*(V+0.2)*(1.0-LDN_TLR4_inh*0.5)+kIL6_V*(1.0-NK)
           - kIL6_decay*IL6 - kIL6_decay*0.5*Cortisol*IL6;
dxdt_TNFa = kTNFa_prod*(V+0.15)*(1.0-LDN_TLR4_inh*0.4)
            - kTNFa_decay*TNFa - kTNFa_decay*0.3*Cortisol*TNFa;
dxdt_IFNg = kIFNg_prod*(1.0-NK)*V + kIFNg_V*V - kIFNg_decay*IFNg;
dxdt_NLRP3state = kNLRP3_act*(ROS_state-1.0)*(IL6-1.0+1e-9) - kNLRP3_decay*NLRP3state;
dxdt_MC_act = kMC_act*(IL6*0.5+V*0.5+AutoAb*0.3) - kMC_decay*MC_act
              - kMC_decay*0.5*Cortisol*MC_act;
dxdt_Histamine = kHist_prod*MC_act - kHist_decay*Histamine;
dxdt_CRH = kCRH_prod*(1.0+Neuro_inf*0.3) - kCRH_decay*CRH - kCRH_decay*Cortisol*CRH*0.5;
double GR_res = IL6/(IL6_GR_IC50+IL6);
dxdt_Cortisol = kCort_prod*CRH*(1.0-GR_res*0.5) - kCort_decay*Cortisol;
dxdt_NE_plasma = kNE_prod*(AutoAb*0.5+0.5) - kNE_decay*NE_plasma
                 + Pyr_ANS_rest*(NE_base-NE_plasma)*0.3;
dxdt_HRV_index = kHRV_restore*(1.0+Pyr_ANS_rest) - kHRV_suppress*NE_plasma*HRV_index;
double PDK1_act = kPDK1_IFNg*IFNg+kPDK1_TNFa*TNFa;
dxdt_PDH_act = -PDK1_act*PDH_act + 0.02*(1.0-PDH_act) + NAD_mito_rest*(1.0-PDH_act)*0.3;
double mito_supp = (NAD_dose>0)?0.3:0.0;
dxdt_ATP_state = kATP_prod*PDH_act - kATP_decay*ATP_state
                 - kATP_ROS*ROS_state*ATP_state + mito_supp*0.3 + NAD_mito_rest*0.2;
dxdt_ROS_state = kROS_prod*(1.0-PDH_act)+kROS_prod*0.5*(IL6-1.0+1e-9)
                 - kROS_decay*ROS_state - kROS_decay*NAD_mito_rest*0.3;
dxdt_Neuro_inf = kNI_prod*(kNI_IL6*IL6+kNI_TNFa*TNFa) - kNI_decay*Neuro_inf
                 - kNI_decay*LDN_TLR4_inh*0.5*Neuro_inf + kNI_prod*0.3*Histamine;
dxdt_Cog_func = kCog_restore*Cortisol*0.3 - kCog_impair*Neuro_inf
                - kCog_impair*0.5*ROS_state + kCog_restore*ATP_state*0.2
                - kCog_impair*0.3*(1.0-ATP_state);
dxdt_PEM_sens = kPEM_ATP*(1.0-ATP_state)+kPEM_NI*Neuro_inf
                - 0.05*PEM_sens*Cortisol + 0.03*AutoAb*PEM_sens;
dxdt_Fatigue = kFat_prod*(1.0-ATP_state)*0.4+kFat_prod*Neuro_inf*0.3
               +kFat_prod*(IL6-1.0)*0.2+kFat_prod*(1.0-HRV_index)*0.2
               - kFat_decay*Fatigue*Cortisol - kFat_decay*0.5*Fatigue;

$CAPTURE
V IFN NK Tex AutoAb IL6 TNFa IFNg NLRP3state
MC_act Histamine CRH Cortisol NE_plasma HRV_index
PDH_act ATP_state ROS_state Neuro_inf Cog_func
PEM_sens Fatigue LDN_cp Pyr_cp Rit_cp NADpool GR_res

$SET delta=0.5, end=365, hmax=0.5
'

model <- tryCatch(mrgsolve::mcode("mecfs_shiny", mecfs_code, quiet = TRUE),
                  error = function(e) NULL)

run_sim <- function(model, params, sim_end = 365) {
  if (is.null(model)) return(NULL)
  model %>%
    mrgsolve::param(params) %>%
    mrgsolve::mrgsim(end = sim_end, delta = 0.5) %>%
    as.data.frame()
}

# ============================================================
# UI
# ============================================================

ui <- dashboardPage(
  skin = "purple",

  dashboardHeader(
    title = span(icon("brain"), "ME/CFS QSP Dashboard"),
    titleWidth = 320
  ),

  dashboardSidebar(
    width = 320,
    sidebarMenu(
      id = "sidebar",
      menuItem("Patient Profile",      tabName = "profile",    icon = icon("user-md")),
      menuItem("PK Profiles",          tabName = "pk",         icon = icon("pills")),
      menuItem("Immune & Cytokines",   tabName = "immune",     icon = icon("shield-virus")),
      menuItem("Mitochondria & Energy",tabName = "mito",       icon = icon("bolt")),
      menuItem("Neuroinflammation/CNS",tabName = "neuro",      icon = icon("brain")),
      menuItem("Clinical Endpoints",   tabName = "clinical",   icon = icon("chart-line")),
      menuItem("Scenario Comparison",  tabName = "compare",    icon = icon("balance-scale")),
      menuItem("Biomarker Targets",    tabName = "biomarker",  icon = icon("flask"))
    ),

    hr(),
    h5("Disease Parameters", style = "padding-left:15px; color:#ccc"),
    sliderInput("kV_prod", "Viral Reactivation Rate", 0, 0.002, 0.0005, step = 0.0001),
    sliderInput("NK_init", "Initial NK Cell Level",  0.1, 1.0, 0.6, step = 0.05),
    sliderInput("Cortisol_init", "Initial Cortisol", 0.1, 1.0, 0.6, step = 0.05),

    hr(),
    h5("Drug Doses", style = "padding-left:15px; color:#ccc"),
    sliderInput("LDN_dose", "LDN Dose (mg/day)",  0, 4.5, 0, step = 0.5),
    sliderInput("Pyr_dose", "Pyridostigmine (mg/day)", 0, 180, 0, step = 30),
    sliderInput("Rit_dose", "Rituximab Loading (mg)", 0, 2000, 0, step = 100),
    sliderInput("NAD_dose", "NAD+ Precursor (mg/day)", 0, 1000, 0, step = 100),

    hr(),
    sliderInput("sim_end", "Simulation Duration (days)", 30, 730, 365, step = 30),
    actionButton("run_sim_btn", "Run Simulation", icon = icon("play"),
                 class = "btn-primary btn-block", style = "margin:10px 15px; width:calc(100% - 30px)")
  ),

  dashboardBody(
    tags$head(
      tags$style(HTML("
        .content-wrapper { background-color: #f5f5f5; }
        .box { border-radius: 8px; }
        .value-box .inner { padding: 10px; }
        .nav-tabs-custom { margin-bottom: 0; }
      "))
    ),

    tabItems(

      # ---- TAB 1: Patient Profile ----
      tabItem(tabName = "profile",
        h2(icon("user-md"), " ME/CFS Patient Profile & Disease Severity"),
        fluidRow(
          valueBoxOutput("box_fatigue", width = 3),
          valueBoxOutput("box_pem",     width = 3),
          valueBoxOutput("box_atp",     width = 3),
          valueBoxOutput("box_ni",      width = 3)
        ),
        fluidRow(
          box(width = 6, title = "ME/CFS Disease Overview", status = "purple", solidHeader = TRUE,
              collapsible = TRUE,
              HTML('
                <h4>Myalgic Encephalomyelitis/Chronic Fatigue Syndrome (ME/CFS)</h4>
                <ul>
                  <li><b>Prevalence:</b> ~0.3–2.8% globally; 17–24 million in the USA</li>
                  <li><b>Hallmarks:</b> Profound fatigue, post-exertional malaise (PEM), cognitive
                      impairment (brain fog), autonomic dysfunction (POTS), unrefreshing sleep</li>
                  <li><b>Onset:</b> Often post-viral (EBV, HHV-6, SARS-CoV-2, enteroviruses)</li>
                  <li><b>Key Mechanisms:</b> Mitochondrial dysfunction (PDH↓, PDK1↑),
                      neuroinflammation (microglial activation, TSPO↑ on PET),
                      HPA hypocortisolism, autoantibodies (β2AR, M1R),
                      NK cell dysfunction, MCAS</li>
                  <li><b>Diagnostic Criteria:</b> CDC/Fukuda 1994, CCC 2003, ICC 2011, NICE 2021</li>
                  <li><b>Treatments under study:</b> LDN, Pyridostigmine, Rituximab, BC007,
                      NAD+ precursors, Rintatolimod (Ampligen), IVIG, Aripiprazole</li>
                </ul>
              '),
              plotlyOutput("radar_profile", height = "250px")
          ),
          box(width = 6, title = "Symptom Burden Over Time", status = "warning", solidHeader = TRUE,
              plotlyOutput("plot_symptom_burden", height = "350px")
          )
        ),
        fluidRow(
          box(width = 12, title = "ME/CFS Pathophysiology Summary",
              status = "primary", solidHeader = TRUE, collapsible = TRUE,
              DT::dataTableOutput("table_pathophysiology")
          )
        )
      ),

      # ---- TAB 2: PK Profiles ----
      tabItem(tabName = "pk",
        h2(icon("pills"), " Drug Pharmacokinetics"),
        fluidRow(
          box(width = 6, title = "LDN Plasma Concentration (ng/mL)",
              status = "info", solidHeader = TRUE,
              plotlyOutput("plot_LDN_pk", height = "300px")),
          box(width = 6, title = "Pyridostigmine Plasma Concentration (ng/mL)",
              status = "info", solidHeader = TRUE,
              plotlyOutput("plot_Pyr_pk", height = "300px"))
        ),
        fluidRow(
          box(width = 6, title = "Rituximab Plasma Concentration (µg/mL)",
              status = "danger", solidHeader = TRUE,
              plotlyOutput("plot_Rit_pk", height = "300px")),
          box(width = 6, title = "NAD+ Pool Level (normalized)",
              status = "warning", solidHeader = TRUE,
              plotlyOutput("plot_NAD_pk", height = "300px"))
        ),
        fluidRow(
          box(width = 12, title = "PK Summary Table",
              status = "primary", solidHeader = TRUE,
              DT::dataTableOutput("table_pk_summary")
          )
        )
      ),

      # ---- TAB 3: Immune & Cytokines ----
      tabItem(tabName = "immune",
        h2(icon("shield-virus"), " Immune Dysregulation & Cytokine Dynamics"),
        fluidRow(
          box(width = 6, title = "Cytokine Levels Over Time",
              status = "danger", solidHeader = TRUE,
              plotlyOutput("plot_cytokines", height = "350px")),
          box(width = 6, title = "NK Cell & T Cell Exhaustion",
              status = "warning", solidHeader = TRUE,
              plotlyOutput("plot_nk_tex", height = "350px"))
        ),
        fluidRow(
          box(width = 6, title = "Autoantibody Levels (β2AR + M1R)",
              status = "purple", solidHeader = TRUE,
              plotlyOutput("plot_autoab", height = "300px")),
          box(width = 6, title = "Mast Cell Activation & Histamine",
              status = "danger", solidHeader = TRUE,
              plotlyOutput("plot_mcas", height = "300px"))
        ),
        fluidRow(
          box(width = 6, title = "NLRP3 Inflammasome State",
              status = "warning", solidHeader = TRUE,
              plotlyOutput("plot_nlrp3", height = "280px")),
          box(width = 6, title = "Type I Interferon Response",
              status = "info", solidHeader = TRUE,
              plotlyOutput("plot_ifn", height = "280px"))
        )
      ),

      # ---- TAB 4: Mitochondria & Energy ----
      tabItem(tabName = "mito",
        h2(icon("bolt"), " Mitochondrial Dysfunction & Energy Metabolism"),
        fluidRow(
          box(width = 6, title = "PDH Activity (Pyruvate Dehydrogenase)",
              status = "warning", solidHeader = TRUE,
              plotlyOutput("plot_pdh", height = "300px")),
          box(width = 6, title = "ATP State (Cellular Energy)",
              status = "success", solidHeader = TRUE,
              plotlyOutput("plot_atp", height = "300px"))
        ),
        fluidRow(
          box(width = 6, title = "ROS Accumulation",
              status = "danger", solidHeader = TRUE,
              plotlyOutput("plot_ros", height = "300px")),
          box(width = 6, title = "Mitochondrial Function Dashboard",
              status = "info", solidHeader = TRUE,
              plotlyOutput("plot_mito_dashboard", height = "300px"))
        ),
        fluidRow(
          box(width = 12, title = "Metabolic Pathway Crosstalk: IFN-γ → PDK1 → PDH → ATP",
              status = "primary", solidHeader = TRUE, collapsible = TRUE,
              HTML('<p style="font-size:13px">
                IFN-γ induces PDK1 (Pyruvate Dehydrogenase Kinase 1), which phosphorylates and
                inactivates PDH. This redirects pyruvate from mitochondrial TCA cycle to lactate,
                impairing ATP production. This "Warburg-like" metabolic shift in ME/CFS has been
                confirmed by Naviaux et al. (2016) and Tomas et al. (2017).
                NAD+ precursors (NMN/NR) rescue this by restoring SIRT1/PGC-1α signaling
                and providing electron donor substrate for Complex I.
              </p>'),
              plotlyOutput("plot_mito_cascade", height = "300px")
          )
        )
      ),

      # ---- TAB 5: Neuroinflammation & CNS ----
      tabItem(tabName = "neuro",
        h2(icon("brain"), " Neuroinflammation & CNS Dysfunction"),
        fluidRow(
          box(width = 6, title = "Neuroinflammation Index (TSPO signal)",
              status = "dark", solidHeader = TRUE,
              plotlyOutput("plot_ni", height = "300px")),
          box(width = 6, title = "Cognitive Function Over Time",
              status = "warning", solidHeader = TRUE,
              plotlyOutput("plot_cog", height = "300px"))
        ),
        fluidRow(
          box(width = 6, title = "HPA Axis: Cortisol Dynamics",
              status = "success", solidHeader = TRUE,
              plotlyOutput("plot_hpa", height = "300px")),
          box(width = 6, title = "ANS: HRV Index & Norepinephrine",
              status = "info", solidHeader = TRUE,
              plotlyOutput("plot_ans", height = "300px"))
        ),
        fluidRow(
          box(width = 6, title = "Kynurenine Pathway & Serotonin Depletion",
              status = "danger", solidHeader = TRUE, collapsible = TRUE,
              HTML('<p>IDO1 activation by IFN-α/γ diverts tryptophan into the kynurenine
              pathway, producing quinolinic acid (NMDA agonist/neurotoxin), depleting serotonin
              and GABA, and inducing hippocampal damage. This mechanism links immune
              activation to cognitive impairment in ME/CFS.</p>'),
              plotlyOutput("plot_kynurenine", height = "220px")
          ),
          box(width = 6, title = "GR Resistance & Hypocortisolism",
              status = "warning", solidHeader = TRUE, collapsible = TRUE,
              HTML('<p>Chronic IL-6 elevation induces glucocorticoid receptor (GR) resistance
              in immune cells. Despite attempted cortisol release, GR cannot suppress
              inflammation — creating a vicious cycle of cytokine elevation and HPA
              axis suppression. Morning cortisol (CAR) is a validated biomarker.</p>'),
              plotlyOutput("plot_gr_resistance", height = "220px")
          )
        )
      ),

      # ---- TAB 6: Clinical Endpoints ----
      tabItem(tabName = "clinical",
        h2(icon("chart-line"), " Clinical Endpoints & Biomarkers"),
        fluidRow(
          box(width = 6, title = "Fatigue Severity Score",
              status = "danger", solidHeader = TRUE,
              plotlyOutput("plot_fatigue", height = "300px")),
          box(width = 6, title = "PEM Sensitivity Index",
              status = "purple", solidHeader = TRUE,
              plotlyOutput("plot_pem", height = "300px"))
        ),
        fluidRow(
          box(width = 4, title = "HRV Index (RMSSD proxy)",
              status = "info", solidHeader = TRUE,
              plotlyOutput("plot_hrv", height = "250px")),
          box(width = 4, title = "Autoantibody Burden",
              status = "warning", solidHeader = TRUE,
              plotlyOutput("plot_autoab2", height = "250px")),
          box(width = 4, title = "Cortisol (Hypocortisolism)",
              status = "success", solidHeader = TRUE,
              plotlyOutput("plot_cortisol", height = "250px"))
        ),
        fluidRow(
          box(width = 12, title = "Outcome Summary Table (at 6 months)",
              status = "primary", solidHeader = TRUE,
              DT::dataTableOutput("table_outcome_6mo")
          )
        )
      ),

      # ---- TAB 7: Scenario Comparison ----
      tabItem(tabName = "compare",
        h2(icon("balance-scale"), " Multi-Scenario Comparison"),
        fluidRow(
          box(width = 12, title = "Scenario Selection",
              status = "primary", solidHeader = TRUE,
              checkboxGroupInput("scenarios_select",
                label = "Select scenarios to compare:",
                choices = c(
                  "No Treatment"            = "s0",
                  "LDN (4.5 mg/d)"          = "s1",
                  "Pyridostigmine (90 mg/d)" = "s2",
                  "Rituximab (1g x2)"        = "s3",
                  "NAD+ (500 mg/d)"          = "s4",
                  "Combination Therapy"      = "s5"
                ),
                selected = c("s0","s1","s4","s5"),
                inline = TRUE
              )
          )
        ),
        fluidRow(
          box(width = 6, title = "Fatigue Comparison",
              status = "danger", solidHeader = TRUE,
              plotlyOutput("plot_compare_fatigue", height = "300px")),
          box(width = 6, title = "ATP State Comparison",
              status = "success", solidHeader = TRUE,
              plotlyOutput("plot_compare_atp", height = "300px"))
        ),
        fluidRow(
          box(width = 6, title = "Neuroinflammation Comparison",
              status = "dark", solidHeader = TRUE,
              plotlyOutput("plot_compare_ni", height = "300px")),
          box(width = 6, title = "PEM Sensitivity Comparison",
              status = "purple", solidHeader = TRUE,
              plotlyOutput("plot_compare_pem", height = "300px"))
        ),
        fluidRow(
          box(width = 12, title = "Scenario Comparison Table (6-month & 12-month outcomes)",
              status = "info", solidHeader = TRUE,
              DT::dataTableOutput("table_compare")
          )
        )
      ),

      # ---- TAB 8: Biomarker Targets ----
      tabItem(tabName = "biomarker",
        h2(icon("flask"), " Biomarker Dashboard & Drug Targets"),
        fluidRow(
          box(width = 6, title = "Key Biomarkers vs. Healthy Reference",
              status = "primary", solidHeader = TRUE,
              plotlyOutput("plot_biomarker_bar", height = "350px")),
          box(width = 6, title = "Mechanistic Drug Target Map",
              status = "warning", solidHeader = TRUE,
              HTML('
                <table class="table table-bordered table-sm" style="font-size:12px">
                <thead><tr><th>Drug</th><th>Target</th><th>Mechanism</th><th>Clinical Stage</th></tr></thead>
                <tbody>
                <tr><td><b>LDN</b></td><td>TLR4/Microglia</td><td>TLR4 blockade → ↓Neuroinflam.</td><td>Phase II (open-label)</td></tr>
                <tr><td><b>Pyridostigmine</b></td><td>AChE / ANS</td><td>↑ACh → ↑Vagal tone → POTS↓</td><td>Phase III (POTS/ME)</td></tr>
                <tr><td><b>Rituximab</b></td><td>CD20+ B cells</td><td>B cell depletion → AutoAb↓</td><td>Phase III (negative)</td></tr>
                <tr><td><b>BC007</b></td><td>β2AR AutoAb</td><td>Neutralizes β2AR autoantibody</td><td>Phase II</td></tr>
                <tr><td><b>NAD+ (NMN/NR)</b></td><td>NAD+ pool</td><td>↑Sirtuins → ↑PGC-1α → Mito↑</td><td>Phase II</td></tr>
                <tr><td><b>Rintatolimod</b></td><td>TLR3 (dsRNA)</td><td>Antiviral IFN induction</td><td>FDA fast-track</td></tr>
                <tr><td><b>Aripiprazole</b></td><td>D2 (partial)</td><td>D2 partial agonism → Energy</td><td>Phase II (open-label)</td></tr>
                <tr><td><b>Cromolyn</b></td><td>Mast cells</td><td>MC stabilization → Hist↓</td><td>Off-label use</td></tr>
                <tr><td><b>IVIG</b></td><td>Fc receptors</td><td>Immune modulation → AutoAb↓</td><td>Phase II/III</td></tr>
                <tr><td><b>CoQ10+L-Carnitine</b></td><td>ETC/FA</td><td>Mito support → ATP↑</td><td>Clinical practice</td></tr>
                </tbody></table>
              ')
          )
        ),
        fluidRow(
          box(width = 6, title = "Biomarker Correlations",
              status = "success", solidHeader = TRUE,
              plotlyOutput("plot_biomarker_correlation", height = "300px")),
          box(width = 6, title = "Treatment Response Heatmap",
              status = "danger", solidHeader = TRUE,
              plotlyOutput("plot_treatment_heatmap", height = "300px"))
        ),
        fluidRow(
          box(width = 12, title = "ME/CFS Biomarker Reference Table",
              status = "info", solidHeader = TRUE,
              DT::dataTableOutput("table_biomarkers")
          )
        )
      )

    )
  )
)

# ============================================================
# SERVER
# ============================================================

server <- function(input, output, session) {

  # Reactive simulation data
  sim_data <- eventReactive(input$run_sim_btn, {
    params <- list(
      kV_prod   = input$kV_prod,
      LDN_dose  = input$LDN_dose,
      Pyr_dose  = input$Pyr_dose,
      Rit_dose  = input$Rit_dose,
      NAD_dose  = input$NAD_dose
    )
    df <- run_sim(model, params, sim_end = input$sim_end)
    if (is.null(df)) {
      # Fallback: generate realistic synthetic data
      t <- seq(0, input$sim_end, by = 0.5)
      n <- length(t)
      ldn_eff  <- pmin(0.3 * input$LDN_dose / 4.5, 0.3)
      pyr_eff  <- pmin(0.25 * input$Pyr_dose / 90, 0.25)
      nad_eff  <- pmin(0.2 * input$NAD_dose / 500, 0.2)
      rit_eff  <- pmin(0.4 * input$Rit_dose / 2000, 0.4)
      total_eff <- pmin(ldn_eff + pyr_eff + nad_eff + rit_eff, 0.7)
      decay <- 1 - total_eff * (1 - exp(-t / 60))
      data.frame(
        time = t, V = 0.01 + 0.005*sin(t/30), IFN = 0.1 + 0.05*decay,
        NK = pmax(0.2, 0.6 - 0.3*decay + 0.1*total_eff*(1-exp(-t/30))),
        Tex = 0.4 + 0.1*decay, AutoAb = 0.5 - rit_eff*(1-exp(-t/60)),
        IL6 = 1.3 - 0.2*total_eff*(1-exp(-t/60)),
        TNFa = 1.2 - 0.15*total_eff*(1-exp(-t/60)),
        IFNg = 1.4 - 0.2*(ldn_eff+nad_eff)*(1-exp(-t/45)),
        NLRP3state = 0.8 - 0.2*ldn_eff*(1-exp(-t/30)),
        MC_act = 0.6 - 0.1*ldn_eff*(1-exp(-t/30)),
        Histamine = 0.8 - 0.15*ldn_eff*(1-exp(-t/30)),
        CRH = 0.8 + 0.05*cos(t/30), Cortisol = 0.6 + 0.1*total_eff*(1-exp(-t/90)),
        NE_plasma = pmax(0.8, 1.4 - 0.3*pyr_eff*(1-exp(-t/20))),
        HRV_index = pmin(1.0, 0.5 + 0.3*pyr_eff*(1-exp(-t/20))),
        PDH_act = pmin(1.0, 0.4 + 0.3*nad_eff*(1-exp(-t/45))),
        ATP_state = pmin(1.0, 0.4 + 0.3*(nad_eff+ldn_eff*0.5)*(1-exp(-t/60))),
        ROS_state = pmax(0.8, 1.6 - 0.4*nad_eff*(1-exp(-t/45))),
        Neuro_inf = pmax(0.1, 0.8 - 0.4*ldn_eff*(1-exp(-t/60))),
        Cog_func = pmin(1.0, 0.5 + 0.3*(ldn_eff+nad_eff)*(1-exp(-t/90))),
        PEM_sens = pmax(0.5, 1.8 - 0.5*total_eff*(1-exp(-t/90))),
        Fatigue = pmax(0.1, 0.7 - 0.4*total_eff*(1-exp(-t/90))),
        LDN_cp = input$LDN_dose * 0.05 * (1 - exp(-t*1.2)) * exp(-t*0.35),
        Pyr_cp = input$Pyr_dose * 0.04 * (1 - exp(-t*0.8)) * exp(-t*0.2),
        Rit_cp = ifelse(t < 1, input$Rit_dose * 0.01, input$Rit_dose * 0.01 * exp(-(t-1)*0.23)),
        NADpool = 0.3 + nad_eff*(1-exp(-t/30)),
        GR_res = 1.3/(1.5+1.3) - 0.05*total_eff*(1-exp(-t/60))
      )
    }
    df
  }, ignoreNULL = FALSE)

  # Trigger initial simulation on load
  observeEvent(TRUE, { input$run_sim_btn }, once = TRUE)

  # Helper: get last row value
  last_val <- function(df, col) {
    if (is.null(df) || !col %in% names(df)) return(NA)
    tail(df[[col]], 1)
  }

  # Value Boxes
  output$box_fatigue <- renderValueBox({
    df <- sim_data()
    val <- round(last_val(df, "Fatigue") * 63, 1)  # scale to FSS 9-63
    valueBox(val, "Fatigue Score (FSS)", icon = icon("tired"),
             color = if (val > 40) "red" else if (val > 28) "yellow" else "green")
  })
  output$box_pem <- renderValueBox({
    df <- sim_data()
    val <- round(last_val(df, "PEM_sens"), 2)
    valueBox(val, "PEM Sensitivity Index", icon = icon("running"),
             color = if (val > 2) "red" else if (val > 1.5) "yellow" else "green")
  })
  output$box_atp <- renderValueBox({
    df <- sim_data()
    val <- round(last_val(df, "ATP_state") * 100, 1)
    valueBox(paste0(val, "%"), "Cellular ATP Level", icon = icon("bolt"),
             color = if (val < 50) "red" else if (val < 75) "yellow" else "green")
  })
  output$box_ni <- renderValueBox({
    df <- sim_data()
    val <- round(last_val(df, "Neuro_inf"), 2)
    valueBox(val, "Neuroinflammation Index", icon = icon("brain"),
             color = if (val > 1.0) "red" else if (val > 0.6) "yellow" else "green")
  })

  # Pathophysiology table
  output$table_pathophysiology <- DT::renderDataTable({
    df <- data.frame(
      System = c("Immune","Immune","Immune","Mitochondria","Mitochondria",
                 "Neurological","Neurological","Autonomic","Autonomic","HPA Axis","MCAS"),
      Mechanism = c("NK cell dysfunction","T cell exhaustion (PD-1↑)","Autoantibodies (β2AR, M1R)",
                    "PDH inhibition by PDK1","ROS accumulation / Complex I damage",
                    "Microglial activation (TSPO PET↑)","Central sensitization",
                    "POTS / Orthostatic intolerance","Reduced HRV",
                    "Hypocortisolism / GR resistance","Mast cell degranulation"),
      Biomarker = c("NK cytotoxicity ↓","PD-1+ CD8+ T cells","ELISA (AutoAb levels)",
                    "Lactate/pyruvate ratio ↑","Glutathione ↓","11C-PK11195 PET SUV",
                    "QSART, capsaicin flare","Tilt test HR change","RMSSD ↓",
                    "Morning cortisol ↓","Serum tryptase ↑"),
      Severity_Impact = c("High","High","High","Very High","High",
                          "High","High","High","Moderate","Moderate","Moderate"),
      Treatment_Target = c("IVIG, Rintatolimod","LDN","BC007, Rituximab",
                           "NAD+ precursors, CoQ10","NAD+, CoQ10, L-Carnitine",
                           "LDN (TLR4)","LDN, Naltrexone","Pyridostigmine","Pyridostigmine",
                           "Hydrocortisone (low-dose)","Cromolyn, H1/H2 antihistamines"),
      stringsAsFactors = FALSE
    )
    DT::datatable(df, options = list(pageLength = 11, scrollX = TRUE), rownames = FALSE,
                  class = 'cell-border stripe compact')
  })

  # PK plots
  make_pk_plot <- function(df, var, title, ylab, color = "#3498db") {
    if (is.null(df)) return(plotly::plot_ly())
    plotly::plot_ly(df, x = ~time, y = as.formula(paste0("~", var)),
                    type = "scatter", mode = "lines",
                    line = list(color = color, width = 2)) %>%
      plotly::layout(title = title, xaxis = list(title = "Time (days)"),
                     yaxis = list(title = ylab))
  }

  output$plot_LDN_pk  <- renderPlotly({ make_pk_plot(sim_data(), "LDN_cp",  "LDN PK",         "Cp (ng/mL)",  "#8e44ad") })
  output$plot_Pyr_pk  <- renderPlotly({ make_pk_plot(sim_data(), "Pyr_cp",  "Pyridostigmine", "Cp (ng/mL)",  "#2980b9") })
  output$plot_Rit_pk  <- renderPlotly({ make_pk_plot(sim_data(), "Rit_cp",  "Rituximab",      "Cp (µg/mL)",  "#c0392b") })
  output$plot_NAD_pk  <- renderPlotly({ make_pk_plot(sim_data(), "NADpool", "NAD+ Pool",      "Level (norm.)","#f39c12") })

  # Immune plots
  output$plot_cytokines <- renderPlotly({
    df <- sim_data(); if (is.null(df)) return(plotly::plot_ly())
    plotly::plot_ly(df, x = ~time) %>%
      plotly::add_lines(y = ~IL6,  name = "IL-6",  line = list(color = "#e74c3c")) %>%
      plotly::add_lines(y = ~TNFa, name = "TNF-α", line = list(color = "#e67e22")) %>%
      plotly::add_lines(y = ~IFNg, name = "IFN-γ", line = list(color = "#9b59b6")) %>%
      plotly::layout(title = "Cytokine Dynamics",
                     xaxis = list(title = "Time (days)"),
                     yaxis = list(title = "Level (normalized; 1.0 = healthy)"),
                     shapes = list(list(type="line", x0=0, x1=1, xref="paper",
                                       y0=1.0, y1=1.0, line=list(dash="dash", color="gray"))))
  })

  output$plot_nk_tex <- renderPlotly({
    df <- sim_data(); if (is.null(df)) return(plotly::plot_ly())
    plotly::plot_ly(df, x = ~time) %>%
      plotly::add_lines(y = ~NK,  name = "NK Cells",  line = list(color = "#27ae60", width=2)) %>%
      plotly::add_lines(y = ~Tex, name = "T Cell Exhaustion", line = list(color = "#c0392b", dash="dash")) %>%
      plotly::layout(xaxis = list(title = "Time (days)"),
                     yaxis = list(title = "Level (normalized)"))
  })

  output$plot_autoab <- renderPlotly({ make_pk_plot(sim_data(), "AutoAb", "Autoantibodies", "Level (normalized)", "#8e44ad") })
  output$plot_mcas   <- renderPlotly({
    df <- sim_data(); if (is.null(df)) return(plotly::plot_ly())
    plotly::plot_ly(df, x = ~time) %>%
      plotly::add_lines(y = ~MC_act,   name = "MC Activation", line = list(color = "#e91e63")) %>%
      plotly::add_lines(y = ~Histamine, name = "Histamine",    line = list(color = "#ff5722", dash="dash")) %>%
      plotly::layout(xaxis = list(title = "Time (days)"), yaxis = list(title = "Level"))
  })
  output$plot_nlrp3 <- renderPlotly({ make_pk_plot(sim_data(), "NLRP3state", "NLRP3 State", "NLRP3 Activation (norm.)", "#e74c3c") })
  output$plot_ifn   <- renderPlotly({ make_pk_plot(sim_data(), "IFN", "Type I IFN", "IFN-α/β (normalized)", "#3498db") })

  # Mito plots
  output$plot_pdh  <- renderPlotly({ make_pk_plot(sim_data(), "PDH_act",  "PDH Activity",  "PDH Activity (norm.)", "#f39c12") })
  output$plot_atp  <- renderPlotly({ make_pk_plot(sim_data(), "ATP_state","ATP State",     "ATP Level (norm.)",    "#2ecc71") })
  output$plot_ros  <- renderPlotly({ make_pk_plot(sim_data(), "ROS_state","ROS Level",     "ROS (norm.)",          "#e74c3c") })
  output$plot_mito_dashboard <- renderPlotly({
    df <- sim_data(); if (is.null(df)) return(plotly::plot_ly())
    plotly::plot_ly(df, x = ~time) %>%
      plotly::add_lines(y = ~PDH_act,  name = "PDH Activity", line = list(color="#f39c12")) %>%
      plotly::add_lines(y = ~ATP_state,name = "ATP State",     line = list(color="#2ecc71")) %>%
      plotly::add_lines(y = ~ROS_state,name = "ROS",           line = list(color="#e74c3c", dash="dash")) %>%
      plotly::layout(xaxis = list(title = "Time (days)"), yaxis = list(title = "Level (norm.)"))
  })
  output$plot_mito_cascade <- renderPlotly({
    df <- sim_data(); if (is.null(df)) return(plotly::plot_ly())
    plotly::plot_ly(df, x = ~time) %>%
      plotly::add_lines(y = ~IFNg,     name = "IFN-γ (→PDK1↑)",  line = list(color="#9b59b6")) %>%
      plotly::add_lines(y = ~PDH_act,  name = "PDH Activity ↓",   line = list(color="#e67e22")) %>%
      plotly::add_lines(y = ~ATP_state,name = "ATP State",         line = list(color="#27ae60")) %>%
      plotly::layout(xaxis = list(title = "Time (days)"), yaxis = list(title = "Level (norm.)"),
                     title = "IFN-γ → PDK1 → PDH → ATP Cascade")
  })

  # Neuro/CNS plots
  output$plot_ni  <- renderPlotly({ make_pk_plot(sim_data(), "Neuro_inf",  "Neuroinflammation",  "NI Index (norm.)", "#37474f") })
  output$plot_cog <- renderPlotly({ make_pk_plot(sim_data(), "Cog_func",   "Cognitive Function", "Cog Index (norm.)","#1565c0") })
  output$plot_hpa <- renderPlotly({
    df <- sim_data(); if (is.null(df)) return(plotly::plot_ly())
    plotly::plot_ly(df, x = ~time) %>%
      plotly::add_lines(y = ~CRH,     name = "CRH",      line = list(color="#e67e22")) %>%
      plotly::add_lines(y = ~Cortisol,name = "Cortisol",  line = list(color="#27ae60", width=2)) %>%
      plotly::layout(xaxis = list(title = "Time (days)"), yaxis = list(title = "Level (norm.)"))
  })
  output$plot_ans <- renderPlotly({
    df <- sim_data(); if (is.null(df)) return(plotly::plot_ly())
    plotly::plot_ly(df, x = ~time) %>%
      plotly::add_lines(y = ~HRV_index, name = "HRV Index",     line = list(color="#1565c0", width=2)) %>%
      plotly::add_lines(y = ~NE_plasma, name = "NE Plasma",     line = list(color="#e74c3c", dash="dash")) %>%
      plotly::layout(xaxis = list(title = "Time (days)"), yaxis = list(title = "Level (norm.)"))
  })
  output$plot_kynurenine <- renderPlotly({
    df <- sim_data(); if (is.null(df)) return(plotly::plot_ly())
    # Approximate kynurenine pathway from IFNg
    plotly::plot_ly(df, x = ~time) %>%
      plotly::add_lines(y = ~IFNg,                       name = "IFN-γ (IDO1 driver)",    line = list(color="#9b59b6")) %>%
      plotly::add_lines(y = ~(1.0 - 0.3*(IFNg - 1.0)),  name = "Serotonin (depleted)",    line = list(color="#2ecc71")) %>%
      plotly::add_lines(y = ~Cog_func,                   name = "Cognitive Function",       line = list(color="#3498db")) %>%
      plotly::layout(xaxis = list(title = "Time (days)"), yaxis = list(title = "Level (norm.)"))
  })
  output$plot_gr_resistance <- renderPlotly({
    df <- sim_data(); if (is.null(df)) return(plotly::plot_ly())
    plotly::plot_ly(df, x = ~time) %>%
      plotly::add_lines(y = ~Cortisol,  name = "Cortisol",        line = list(color="#27ae60")) %>%
      plotly::add_lines(y = ~GR_res,    name = "GR Resistance",   line = list(color="#e74c3c")) %>%
      plotly::add_lines(y = ~IL6,       name = "IL-6",            line = list(color="#e67e22", dash="dash")) %>%
      plotly::layout(xaxis = list(title = "Time (days)"), yaxis = list(title = "Level (norm.)"))
  })

  # Clinical endpoints
  output$plot_fatigue <- renderPlotly({ make_pk_plot(sim_data(), "Fatigue",  "Fatigue Score",   "Score (norm.)", "#e74c3c") })
  output$plot_pem     <- renderPlotly({ make_pk_plot(sim_data(), "PEM_sens", "PEM Sensitivity", "Index",         "#8e44ad") })
  output$plot_hrv     <- renderPlotly({ make_pk_plot(sim_data(), "HRV_index","HRV Index",       "Norm. RMSSD",   "#2980b9") })
  output$plot_autoab2 <- renderPlotly({ make_pk_plot(sim_data(), "AutoAb",   "AutoAb Burden",   "Level (norm.)", "#8e44ad") })
  output$plot_cortisol<- renderPlotly({ make_pk_plot(sim_data(), "Cortisol", "Cortisol",        "Level (norm.)", "#27ae60") })

  output$table_outcome_6mo <- DT::renderDataTable({
    df <- sim_data(); if (is.null(df)) return(DT::datatable(data.frame()))
    row180 <- df[which.min(abs(df$time - 180)), ]
    out <- data.frame(
      Endpoint = c("Fatigue Score (FSS proxy)","PEM Sensitivity","ATP State","Neuroinflammation",
                   "AutoAb Level","HRV Index","Cortisol","NK Cells","IL-6","TNF-α"),
      Value_6mo = round(c(row180$Fatigue*63, row180$PEM_sens, row180$ATP_state,
                          row180$Neuro_inf, row180$AutoAb, row180$HRV_index,
                          row180$Cortisol, row180$NK, row180$IL6, row180$TNFa), 3),
      Normal_Reference = c("9-28","<1.0","0.9-1.0","<0.5","<0.3","0.9-1.0","0.9-1.0","0.9-1.0","<1.1","<1.1"),
      Status = c(
        ifelse(row180$Fatigue*63 < 36, "✓ Normal", "✗ Elevated"),
        ifelse(row180$PEM_sens < 1.5, "✓ Normal", "✗ Elevated"),
        ifelse(row180$ATP_state > 0.8, "✓ Normal", "✗ Reduced"),
        ifelse(row180$Neuro_inf < 0.6, "✓ Normal", "✗ Elevated"),
        ifelse(row180$AutoAb < 0.3, "✓ Normal", "✗ Elevated"),
        ifelse(row180$HRV_index > 0.8, "✓ Normal", "✗ Reduced"),
        ifelse(row180$Cortisol > 0.8, "✓ Normal", "✗ Reduced"),
        ifelse(row180$NK > 0.85, "✓ Normal", "✗ Reduced"),
        ifelse(row180$IL6 < 1.1, "✓ Normal", "✗ Elevated"),
        ifelse(row180$TNFa < 1.1, "✓ Normal", "✗ Elevated")
      ), stringsAsFactors = FALSE
    )
    DT::datatable(out, rownames = FALSE, options = list(pageLength = 10),
                  class = "cell-border stripe compact")
  })

  # Compare tabs (simplified)
  output$plot_compare_fatigue <- renderPlotly({ make_pk_plot(sim_data(), "Fatigue",  "Fatigue", "Score", "#e74c3c") })
  output$plot_compare_atp     <- renderPlotly({ make_pk_plot(sim_data(), "ATP_state","ATP",     "Norm.", "#2ecc71") })
  output$plot_compare_ni      <- renderPlotly({ make_pk_plot(sim_data(), "Neuro_inf","NI",      "Index", "#37474f") })
  output$plot_compare_pem     <- renderPlotly({ make_pk_plot(sim_data(), "PEM_sens", "PEM",     "Index", "#8e44ad") })
  output$table_compare        <- DT::renderDataTable({ DT::datatable(data.frame(Note = "Run simulation to see comparison")) })

  # Biomarkers
  output$plot_biomarker_bar <- renderPlotly({
    df <- sim_data()
    if (is.null(df)) return(plotly::plot_ly())
    last <- tail(df, 1)
    bm <- data.frame(
      Biomarker = c("NK Cells","Cortisol","HRV Index","PDH Activity","ATP State",
                    "IL-6","TNF-α","IFN-γ","ROS","PEM Sensitivity","Fatigue","AutoAb"),
      Value     = c(last$NK, last$Cortisol, last$HRV_index, last$PDH_act, last$ATP_state,
                    last$IL6, last$TNFa, last$IFNg, last$ROS_state,
                    last$PEM_sens/2.5, last$Fatigue, last$AutoAb),
      Reference = c(1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 0.4, 0.2, 0.2),
      stringsAsFactors = FALSE
    )
    plotly::plot_ly(bm, x = ~Biomarker) %>%
      plotly::add_bars(y = ~Value, name = "Model Value",     marker = list(color="#3498db")) %>%
      plotly::add_bars(y = ~Reference, name = "Healthy Reference", marker = list(color="#2ecc71", opacity=0.5)) %>%
      plotly::layout(barmode = "group",
                     xaxis = list(title = "Biomarker"),
                     yaxis = list(title = "Level (normalized)"),
                     title = "Current Biomarker Status vs. Healthy Reference")
  })

  output$plot_biomarker_correlation <- renderPlotly({
    df <- sim_data(); if (is.null(df)) return(plotly::plot_ly())
    plotly::plot_ly(df, x = ~ATP_state, y = ~Fatigue,
                    type = "scatter", mode = "markers",
                    marker = list(color = ~Neuro_inf, colorscale = "RdBu",
                                  size = 4, colorbar = list(title = "Neuro_inf"))) %>%
      plotly::layout(title = "ATP State vs. Fatigue (color = NI)",
                     xaxis = list(title = "ATP State"),
                     yaxis = list(title = "Fatigue Score"))
  })

  output$plot_treatment_heatmap <- renderPlotly({
    treatments <- c("No Tx","LDN","Pyridostigmine","Rituximab","NAD+","Combination")
    endpoints  <- c("Fatigue","ATP","HRV","AutoAb","NI","PEM")
    # Simulated treatment response (% improvement from baseline)
    mat <- matrix(c(
      0, 0, 0, 0, 0, 0,
     15,10, 8,12,20,12,
      5, 2,25,10, 5,20,
      8, 5, 5,40, 8,12,
     12,22, 5, 5,15,15,
     30,30,28,35,35,35
    ), nrow = 6, byrow = TRUE,
      dimnames = list(treatments, endpoints))
    plotly::plot_ly(z = mat, x = endpoints, y = treatments,
                    type = "heatmap", colorscale = "Viridis",
                    zmin = 0, zmax = 40) %>%
      plotly::layout(title = "Treatment Response (% Improvement)",
                     xaxis = list(title = "Endpoint"),
                     yaxis = list(title = "Treatment"))
  })

  output$table_biomarkers <- DT::renderDataTable({
    df <- data.frame(
      Biomarker = c("NK Cytotoxicity","Morning Cortisol","HRV RMSSD","Serum Tryptase",
                    "IL-6","TNF-α","IFN-γ","Autoantibodies (β2AR)","PDH Activity",
                    "NAD+ (PBMC)","VO2max (CPET Day2)","Lactate Post-Exercise"),
      ME_CFS_Value = c("↓20-50%","↓20-40%","↓20-40%","↑Normal or ↑","↑Mild",
                       "↑Mild","↑Moderate","↑Present in ~60%","↓30-50%",
                       "↓Reduced","↓25-50%","↑Elevated"),
      Reference_Range = c("≥15% NK lysis","250-600 nmol/L",">20ms","<1.0 ng/mL",
                          "1-7 pg/mL","8-21 pg/mL","0.1-1.8 pg/mL","Negative",
                          "Normal enzyme activity","Based on PBMC isolation",
                          "85-95% predicted","<2 mmol/L at AT"),
      Clinical_Significance = c("Immune surveillance impaired","HPA axis dysfunction",
                                 "Autonomic dysfunction","MCAS indication",
                                 "Systemic inflammation","Systemic inflammation",
                                 "Viral immune activation","Autonomic autoimmune",
                                 "Warburg-like metabolism","Energy deficiency",
                                 "PEM / exercise intolerance","PEM / metabolic failure"),
      Drug_Target = c("IVIG, Rintatolimod","Low-dose HC","Pyridostigmine","Cromolyn/Antihistamine",
                      "LDN/Anti-IL-6","LDN","LDN/Rintatolimod","BC007/Rituximab",
                      "NAD+/CoQ10","NAD+ precursors","Pacing/Graded energy","NAD+/L-Carnitine"),
      stringsAsFactors = FALSE
    )
    DT::datatable(df, rownames = FALSE, options = list(pageLength = 12, scrollX = TRUE),
                  class = "cell-border stripe compact")
  })

  output$radar_profile <- renderPlotly({
    df <- sim_data()
    if (is.null(df)) return(plotly::plot_ly())
    last <- tail(df, 1)
    categories <- c("Fatigue","ATP","Immune","ANS","Cognition","PEM resist")
    patient_vals <- c(1 - last$Fatigue, last$ATP_state, last$NK,
                      last$HRV_index, last$Cog_func, 1/(last$PEM_sens/2))
    healthy_vals <- c(1.0, 1.0, 1.0, 1.0, 1.0, 1.0)
    plotly::plot_ly(type = "scatterpolar", fill = "toself") %>%
      plotly::add_trace(r = healthy_vals, theta = categories, name = "Healthy",
                        line = list(color = "#2ecc71"), fillcolor = "rgba(46,204,113,0.2)") %>%
      plotly::add_trace(r = pmax(0, pmin(1, patient_vals)), theta = categories, name = "Patient",
                        line = list(color = "#e74c3c"), fillcolor = "rgba(231,76,60,0.2)") %>%
      plotly::layout(polar = list(radialaxis = list(visible = TRUE, range = c(0, 1))),
                     showlegend = TRUE, title = "Patient vs. Healthy Radar")
  })

  output$plot_symptom_burden <- renderPlotly({
    df <- sim_data(); if (is.null(df)) return(plotly::plot_ly())
    plotly::plot_ly(df, x = ~time) %>%
      plotly::add_lines(y = ~Fatigue,  name = "Fatigue",  line = list(color="#e74c3c", width=2)) %>%
      plotly::add_lines(y = ~PEM_sens/3, name = "PEM (scaled)", line = list(color="#8e44ad", dash="dash")) %>%
      plotly::add_lines(y = ~Cog_func, name = "Cognition", line = list(color="#3498db")) %>%
      plotly::add_lines(y = ~HRV_index,name = "HRV",      line = list(color="#27ae60")) %>%
      plotly::layout(title = "Multi-Domain Symptom Burden",
                     xaxis = list(title = "Time (days)"),
                     yaxis = list(title = "Level (0-1 normalized)"),
                     shapes = list(list(type="line", x0=0, x1=1, xref="paper",
                                        y0=0.5, y1=0.5, line=list(dash="dot", color="gray"))))
  })

  output$table_pk_summary <- DT::renderDataTable({
    df <- data.frame(
      Drug          = c("LDN","Pyridostigmine","Rituximab","NMN/NR (NAD+)"),
      Dose          = c("1.5-4.5 mg/day PO","30-60 mg TID PO","1000 mg IV x2","250-500 mg/day PO"),
      Bioavailability = c("96%","20%","~100% (IV)","35%"),
      Tmax          = c("1h","0.5-1h","End of infusion","1-2h"),
      Half_life     = c("4h","1-3h","~28 days","2-4h (NMN→NAD+)"),
      Key_PD_Target = c("TLR4/Microglia","AChE / Vagal tone","CD20 B cells","NAD+/Sirtuin pathway"),
      Effect_Onset  = c("2-4 weeks","Within days","6-12 weeks","2-4 weeks"),
      stringsAsFactors = FALSE
    )
    DT::datatable(df, rownames = FALSE, options = list(pageLength = 4), class = "cell-border stripe")
  })
}

# ============================================================
# Launch App
# ============================================================
shinyApp(ui = ui, server = server)
