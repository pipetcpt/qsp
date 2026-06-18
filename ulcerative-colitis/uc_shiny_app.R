# ============================================================
# Ulcerative Colitis QSP - Interactive Shiny Dashboard
# ============================================================
# Tabs:
#   1. Patient Profile
#   2. Drug PK
#   3. Cytokine & Immune Dynamics
#   4. Disease Activity
#   5. Mucosal Healing
#   6. Scenario Comparison
#   7. Biomarkers & Safety
# ============================================================

suppressPackageStartupMessages({
  library(shiny)
  library(shinydashboard)
  library(mrgsolve)
  library(dplyr)
  library(ggplot2)
  library(tidyr)
  library(DT)
})

# ============================================================
# MRGSOLVE MODEL CODE (inline)
# ============================================================
uc_model_code <- '
$PROB Ulcerative Colitis QSP Model v1.0 (Shiny)

$PARAM
CL_IFX=0.407, V1_IFX=3.28, V2_IFX=3.57, Q_IFX=0.484
kon_IFX=0.097, koff_IFX=0.001, kdeg_IFX=0.0693, Rbase_TNF=0.5
CL_VDZ=0.271, V1_VDZ=5.24
kon_VDZ=0.217, koff_VDZ=0.033, kdeg_VDZ=0.0528, Rbase_a4b7=2.0
ka_TOF=83.28, F_TOF=0.74
ka_OZA=28.8, CL_OZA=1056.0, V_OZA=3760.0, k_OZA_M1=0.432, CL_M1=252.0, V_M1=4850.0
ka_UST=0.44, CL_UST=0.192, V1_UST=4.62, F_UST_SC=0.615
kin_TNF=0.020, kout_TNF=0.040
kin_IL17=0.015, kout_IL17=0.030
kin_IL13=0.010, kout_IL13=0.025
kin_IL10=0.012, kout_IL10=0.030
kin_Th2=10.0, kout_Th2=0.030
kin_Th17=5.0, kout_Th17=0.030
kin_Treg=8.0, kout_Treg=0.025
kin_Neut=50.0, kout_Neut=0.15
kdam=0.002, krep=0.008, kCRP=0.50, kFC=2.00, Hill_n=2.0
EC50_TNF=0.30, EC50_IL13=0.20
Emax_IFX=0.92, EC50_IFX=1.5
Emax_VDZ=0.85, EC50_VDZ=10.0
Emax_TOF=0.88, EC50_TOF=50.0
Emax_OZA=0.75, EC50_OZA=0.5
Emax_UST=0.80, EC50_UST=1.0
MAYO0=9.0

$CMT
IFX_C1 IFX_C2 IFX_RC
VDZ_C1 VDZ_RC
TOF_GI TOF_C1
OZA_GI OZA_C1 OZA_M1
UST_DEPOT UST_C1
TNFa IL17 IL13 IL10
Th2 Th17 Treg Neutrophil
MayoScore MucosalHealing CRP FC

$INIT
IFX_C1=0, IFX_C2=0, IFX_RC=0
VDZ_C1=0, VDZ_RC=0
TOF_GI=0, TOF_C1=0
OZA_GI=0, OZA_C1=0, OZA_M1=0
UST_DEPOT=0, UST_C1=0
TNFa=0.5, IL17=0.5, IL13=0.4, IL10=0.4
Th2=400, Th17=200, Treg=300, Neutrophil=50
MayoScore=9.0, MucosalHealing=0.1, CRP=4.0, FC=10.0

$ODE
double IFX_conc_ugmL = (IFX_C1/V1_IFX)*1000.0;
double E_IFX = Emax_IFX*pow(IFX_conc_ugmL,Hill_n)/(pow(EC50_IFX,Hill_n)+pow(IFX_conc_ugmL,Hill_n));
dxdt_IFX_C1 = -(CL_IFX+Q_IFX)/V1_IFX*IFX_C1 + Q_IFX/V2_IFX*IFX_C2
              - kon_IFX*(IFX_C1/V1_IFX)*TNFa + koff_IFX*IFX_RC;
dxdt_IFX_C2 = Q_IFX/V1_IFX*IFX_C1 - Q_IFX/V2_IFX*IFX_C2;
dxdt_IFX_RC = kon_IFX*(IFX_C1/V1_IFX)*TNFa - koff_IFX*IFX_RC - kdeg_IFX*IFX_RC;

double VDZ_conc_ugmL = (VDZ_C1/V1_VDZ)*1000.0;
double E_VDZ = Emax_VDZ*pow(VDZ_conc_ugmL,Hill_n)/(pow(EC50_VDZ,Hill_n)+pow(VDZ_conc_ugmL,Hill_n));
dxdt_VDZ_C1 = -CL_VDZ/V1_VDZ*VDZ_C1 - kon_VDZ*(VDZ_C1/V1_VDZ)*Rbase_a4b7 + koff_VDZ*VDZ_RC;
dxdt_VDZ_RC = kon_VDZ*(VDZ_C1/V1_VDZ)*Rbase_a4b7 - koff_VDZ*VDZ_RC - kdeg_VDZ*VDZ_RC;

double TOF_conc_ngmL = (TOF_C1/87.0)*1e6;
double E_TOF = Emax_TOF*pow(TOF_conc_ngmL,Hill_n)/(pow(EC50_TOF,Hill_n)+pow(TOF_conc_ngmL,Hill_n));
dxdt_TOF_GI = -ka_TOF*TOF_GI;
dxdt_TOF_C1 = ka_TOF*F_TOF*TOF_GI - (1027.2/87.0)*TOF_C1;

double OZA_M1_conc_ngmL = (OZA_M1/V_M1)*1e6;
double E_OZA = Emax_OZA*pow(OZA_M1_conc_ngmL,Hill_n)/(pow(EC50_OZA,Hill_n)+pow(OZA_M1_conc_ngmL,Hill_n));
dxdt_OZA_GI = -ka_OZA*OZA_GI;
dxdt_OZA_C1 = ka_OZA*OZA_GI - (CL_OZA/V_OZA + k_OZA_M1)*OZA_C1;
dxdt_OZA_M1 = k_OZA_M1*(V_OZA/V_M1)*OZA_C1 - (CL_M1/V_M1)*OZA_M1;

double UST_conc_ugmL = (UST_C1/V1_UST)*1000.0;
double E_UST = Emax_UST*pow(UST_conc_ugmL,Hill_n)/(pow(EC50_UST,Hill_n)+pow(UST_conc_ugmL,Hill_n));
dxdt_UST_DEPOT = -ka_UST*UST_DEPOT;
dxdt_UST_C1 = ka_UST*F_UST_SC*UST_DEPOT - CL_UST/V1_UST*UST_C1;

double inh_TNF=E_IFX, inh_homing=E_VDZ, inh_JAK=E_TOF, inh_S1P=E_OZA, inh_IL12_23=E_UST;

dxdt_TNFa = kin_TNF*(1.0+2.0*MayoScore/12.0)*(1.0-inh_TNF) - kout_TNF*TNFa;
dxdt_IL17 = kin_IL17*(Th17/200.0)*(1.0-inh_JAK*0.6)*(1.0-inh_IL12_23*0.4) - kout_IL17*IL17;
dxdt_IL13 = kin_IL13*(Th2/400.0)*(1.0-inh_JAK*0.7)*(1.0-inh_homing*0.5) - kout_IL13*IL13;
dxdt_IL10 = kin_IL10*(Treg/300.0)*(1.0+inh_homing*0.2) - kout_IL10*IL10;

dxdt_Th2 = kin_Th2*(1.0-inh_homing)*(1.0-inh_S1P*0.6) - kout_Th2*Th2;
dxdt_Th17 = kin_Th17*(1.0-inh_JAK*0.5)*(1.0-inh_S1P*0.5)*(1.0-inh_IL12_23*0.3) - kout_Th17*Th17;
dxdt_Treg = kin_Treg*(1.0+inh_homing*0.3) - kout_Treg*Treg;
dxdt_Neutrophil = kin_Neut*(1.0+IL17/0.5)*(1.0-inh_JAK*0.4) - kout_Neut*Neutrophil;

double IL10_ratio=IL10/0.4, TNF_ratio=TNFa/Rbase_TNF, IL17_ratio=IL17/0.5, IL13_ratio=IL13/0.4;
double damage_driver = TNF_ratio*(1.0-IL10_ratio*0.3) + 0.5*IL13_ratio + 0.3*IL17_ratio;
dxdt_MayoScore = kdam*damage_driver*(12.0-MayoScore) - krep*IL10_ratio*MayoScore;
dxdt_MucosalHealing = krep*IL10_ratio*(1.0-MucosalHealing) - kdam*damage_driver*MucosalHealing;
dxdt_CRP = kCRP*(TNF_ratio+0.5*IL17_ratio-0.5*IL10_ratio-1.0) - 0.15*CRP;
dxdt_FC = kFC*(Neutrophil/50.0-1.0+0.5*MayoScore/12.0) - 0.10*FC;

$TABLE
double IFX_conc   = (IFX_C1/V1_IFX)*1000.0;
double VDZ_conc   = (VDZ_C1/V1_VDZ)*1000.0;
double TOF_conc_t = (TOF_C1/87.0)*1e6;
double UST_conc   = (UST_C1/V1_UST)*1000.0;
double OZA_M1_conc= (OZA_M1/V_M1)*1e6;
double Mayo_total = MayoScore;
double MH_index   = MucosalHealing;
double CRP_val    = CRP+1.0;
double FC_val     = FC+50.0;
double clin_remission = (MayoScore<=2.0) ? 1.0 : 0.0;
double MH_resp    = (MucosalHealing>=0.7) ? 1.0 : 0.0;
double deep_remission = ((MayoScore<=2.0)&&(MucosalHealing>=0.7)) ? 1.0 : 0.0;

$CAPTURE
IFX_conc VDZ_conc TOF_conc_t UST_conc OZA_M1_conc
Mayo_total MH_index CRP_val FC_val
TNFa IL17 IL13 IL10 Th2 Th17 Treg Neutrophil
clin_remission MH_resp deep_remission
'

# Pre-compile model at startup
UC_MOD <- mcode("UC_QSP_shiny", uc_model_code, quiet = TRUE)

# ============================================================
# HELPER: BUILD DOSING EVENTS
# ============================================================
build_events <- function(drug, bw = 70) {
  if (drug == "Placebo") {
    return(ev(time = 0, amt = 0, cmt = "IFX_C1"))
  } else if (drug == "Infliximab") {
    dose <- 5 * bw
    return(ev(time = c(0,14,42,98,154,210,266,322), amt = dose, cmt = "IFX_C1", rate = -2))
  } else if (drug == "Vedolizumab") {
    return(ev(time = c(0,14,42,98,154,210,266,322), amt = 300, cmt = "VDZ_C1", rate = -2))
  } else if (drug == "Tofacitinib") {
    t_ind  <- seq(0,   55.5, by = 0.5)
    t_main <- seq(56, 364.5, by = 0.5)
    return(ev(time = c(t_ind, t_main),
              amt  = c(rep(10, length(t_ind)), rep(5, length(t_main))),
              cmt  = "TOF_GI"))
  } else if (drug == "Ustekinumab") {
    ev_iv <- ev(time = 0,   amt = 520, cmt = "UST_C1")
    ev_sc <- ev(time = c(56,112,168,224,280,336), amt = 90, cmt = "UST_DEPOT")
    return(ev_iv + ev_sc)
  } else if (drug == "Ozanimod") {
    return(ev(time = seq(0, 364, by = 1), amt = 0.92, cmt = "OZA_GI"))
  }
  ev(time = 0, amt = 0, cmt = "IFX_C1")
}

# ============================================================
# HELPER: RUN ONE SIMULATION
# ============================================================
run_one <- function(mod, drug, mayo0 = 9, tend = 365) {
  events <- build_events(drug)
  mod2 <- param(mod, MAYO0 = mayo0)
  iv   <- init(mod2)
  iv["MayoScore"] <- mayo0
  mod2 <- init(mod2, iv)
  tryCatch({
    out <- mod2 %>% ev(events) %>%
      mrgsim(end = tend, delta = 1) %>%
      as.data.frame()
    out$drug <- drug
    out
  }, error = function(e) {
    message("Simulation error for ", drug, ": ", e$message)
    NULL
  })
}

# ============================================================
# CONSTANTS
# ============================================================
SCEN_COLORS <- c(
  Placebo      = "#999999",
  Infliximab   = "#E41A1C",
  Vedolizumab  = "#377EB8",
  Tofacitinib  = "#4DAF4A",
  Ustekinumab  = "#984EA3",
  Ozanimod     = "#FF7F00"
)
DRUGS <- names(SCEN_COLORS)

PK_PARAMS <- data.frame(
  Drug        = c("Infliximab","Vedolizumab","Tofacitinib","Ustekinumab","Ozanimod"),
  Route       = c("IV","IV","Oral","IV/SC","Oral"),
  Dose        = c("5 mg/kg wk0,2,6 Q8W","300 mg wk0,2,6 Q8W",
                  "10 mg BID (ind.), 5 mg BID","~520 mg IV; 90 mg SC Q8W","0.92 mg QD"),
  CL          = c("0.407 L/day","0.271 L/day","1027 L/day","0.192 L/day","1056 L/day"),
  Vd          = c("3.28 L","5.24 L","87 L","4.62 L","3760 L"),
  t_half      = c("~9-12 days","~25 days","~3 h","~20 days","~21 h"),
  Mechanism   = c("Anti-TNFa (TMDD)","Anti-a4b7 integrin","JAK1/3 inhibitor",
                  "Anti-IL-12/23 p40","S1PR1 modulator"),
  Trial       = c("ACT1/ACT2","GEMINI 1","OCTAVE","UNIFI","TRUE NORTH"),
  stringsAsFactors = FALSE
)

THEME_UC <- theme_bw() +
  theme(legend.position = "bottom", legend.title = element_blank(),
        strip.background = element_rect(fill = "#D6EAF8"),
        strip.text       = element_text(face = "bold"),
        plot.title       = element_text(face = "bold", hjust = 0.5))

# ============================================================
# UI
# ============================================================
ui <- dashboardPage(
  skin = "blue",

  dashboardHeader(title = "UC QSP Dashboard"),

  dashboardSidebar(
    sidebarMenu(
      menuItem("Patient Profile",     tabName = "profile",    icon = icon("user-md")),
      menuItem("Drug PK",             tabName = "pk",         icon = icon("chart-line")),
      menuItem("Cytokine & Immune",   tabName = "cytokines",  icon = icon("bacteria")),
      menuItem("Disease Activity",    tabName = "disease",    icon = icon("stethoscope")),
      menuItem("Mucosal Healing",     tabName = "mucosal",    icon = icon("heartbeat")),
      menuItem("Scenario Comparison", tabName = "comparison", icon = icon("balance-scale")),
      menuItem("Biomarkers & Safety", tabName = "safety",     icon = icon("shield-alt"))
    )
  ),

  dashboardBody(
    tags$head(tags$style(HTML(
      ".content-wrapper{background-color:#f5f5f5;}
       .box{border-radius:6px;}
       .skin-blue .main-header .logo{background-color:#1a5276;}
       .skin-blue .main-header .navbar{background-color:#1f618d;}"
    ))),

    tabItems(

      # ====================================================
      # TAB 1: PATIENT PROFILE
      # ====================================================
      tabItem(tabName = "profile",
        fluidRow(
          box(title = "Baseline Patient Parameters", width = 4, status = "primary", solidHeader = TRUE,
            numericInput("mayo0",  "Baseline Mayo Score (0-12):", value = 9,  min = 3,  max = 12,   step = 1),
            numericInput("crp0",   "Baseline CRP (mg/L):",        value = 25, min = 1,  max = 200,  step = 1),
            numericInput("fc0",    "Fecal Calprotectin (ug/g):",  value = 800,min = 50, max = 5000, step = 50),
            numericInput("dur_yr", "Disease Duration (years):",   value = 3,  min = 0,  max = 30,   step = 0.5),
            numericInput("bw",     "Body Weight (kg):",           value = 70, min = 40, max = 130,  step = 5),
            hr(),
            h5(strong("Prior Biologics:")),
            checkboxGroupInput("prior_bio", NULL,
              choices  = c("Anti-TNF","Vedolizumab","Ustekinumab","JAK inhibitor"),
              selected = NULL),
            hr(),
            h5(strong("Comorbidities:")),
            checkboxGroupInput("comorbid", NULL,
              choices  = c("PSC","Arthropathy","Uveitis","Anemia","Osteoporosis"),
              selected = NULL),
            hr(),
            selectInput("drug_profile", "Select Treatment:",
                        choices = DRUGS, selected = "Infliximab"),
            actionButton("run_sim", "Run Simulation",
                         icon = icon("play"), class = "btn-success btn-lg btn-block")
          ),

          box(title = "Current Disease Status", width = 4, status = "warning", solidHeader = TRUE,
            valueBoxOutput("vbox_mayo", width = 12),
            valueBoxOutput("vbox_crp",  width = 12),
            valueBoxOutput("vbox_fc",   width = 12)
          ),

          box(title = "Patient Summary", width = 4, status = "info", solidHeader = TRUE,
            tableOutput("patient_summary")
          )
        ),
        fluidRow(
          box(title = "Mayo Score Classification", width = 6, status = "info",
            tableOutput("mayo_class_table")
          ),
          box(title = "Treatment Selection Guide", width = 6, status = "warning",
            HTML('<b>Step-up therapy in moderate-severe UC:</b>
                  <ol>
                  <li><b>5-ASA (mesalamine)</b>: mild disease</li>
                  <li><b>Corticosteroids</b>: moderate flares, short-term</li>
                  <li><b>Immunomodulators (AZA/6-MP)</b>: steroid-dependent</li>
                  <li><b>Biologics / small molecules</b> (Mayo &gt;= 6):</li>
                  </ol>
                  <ul>
                  <li><b>Infliximab</b>: Anti-TNFa, rapid onset, IV</li>
                  <li><b>Vedolizumab</b>: Gut-selective alpha4beta7, IV</li>
                  <li><b>Ustekinumab</b>: Anti-IL-12/23, favourable safety</li>
                  <li><b>Tofacitinib</b>: Oral JAK1/3 inhibitor</li>
                  <li><b>Ozanimod</b>: Oral S1PR1 modulator</li>
                  </ul>')
          )
        )
      ),

      # ====================================================
      # TAB 2: DRUG PK
      # ====================================================
      tabItem(tabName = "pk",
        fluidRow(
          box(title = "PK Parameters & Trial References", width = 12, status = "primary", solidHeader = TRUE,
            DTOutput("pk_table")
          )
        ),
        fluidRow(
          box(title = "Options", width = 3, status = "info",
            selectInput("pk_drug", "Drug:",
                        choices = DRUGS[-1], selected = "Infliximab"),
            sliderInput("pk_weeks", "Duration (weeks):", min = 8, max = 52, value = 52, step = 4)
          ),
          box(title = "Concentration-Time Profile", width = 9, status = "primary",
            plotOutput("pk_plot", height = "340px")
          )
        ),
        fluidRow(
          box(title = "Trough Concentrations", width = 6, status = "warning",
            plotOutput("trough_plot", height = "300px")
          ),
          box(title = "Emax PD Relationship", width = 6, status = "info",
            plotOutput("emax_plot", height = "300px")
          )
        )
      ),

      # ====================================================
      # TAB 3: CYTOKINE & IMMUNE DYNAMICS
      # ====================================================
      tabItem(tabName = "cytokines",
        fluidRow(
          box(title = "Options", width = 3, status = "info",
            selectInput("cyto_drug",  "Treatment:", choices = DRUGS, selected = "Infliximab"),
            sliderInput("cyto_weeks", "Duration (weeks):", min = 8, max = 52, value = 52, step = 4),
            sliderInput("cyto_mayo0", "Baseline Mayo:", min = 3, max = 12, value = 9, step = 1)
          ),
          box(title = "Cytokine Time-Courses (nM)", width = 9, status = "primary",
            plotOutput("cyto_plot", height = "360px")
          )
        ),
        fluidRow(
          box(title = "Immune Cell Population Dynamics", width = 12, status = "info",
            plotOutput("immune_plot", height = "360px")
          )
        )
      ),

      # ====================================================
      # TAB 4: DISEASE ACTIVITY
      # ====================================================
      tabItem(tabName = "disease",
        fluidRow(
          box(title = "Options", width = 3, status = "info",
            selectInput("da_drug",   "Treatment:", choices = DRUGS, selected = "Infliximab"),
            sliderInput("da_mayo0",  "Baseline Mayo:", min = 3, max = 12, value = 9, step = 1),
            checkboxInput("da_all",  "Overlay all treatments", value = FALSE)
          ),
          box(title = "Mayo Total Score Over Time", width = 9, status = "danger",
            plotOutput("mayo_plot", height = "300px")
          )
        ),
        fluidRow(
          box(title = "Partial Mayo Score", width = 6, status = "warning",
            plotOutput("pmayo_plot", height = "260px")
          ),
          box(title = "CRP & Fecal Calprotectin", width = 6, status = "warning",
            plotOutput("biom_da_plot", height = "260px")
          )
        ),
        fluidRow(
          box(title = "Response Rates at Key Timepoints", width = 12, status = "success",
            tableOutput("response_table")
          )
        )
      ),

      # ====================================================
      # TAB 5: MUCOSAL HEALING
      # ====================================================
      tabItem(tabName = "mucosal",
        fluidRow(
          box(title = "Options", width = 3, status = "info",
            selectInput("mh_drug",  "Treatment:", choices = DRUGS, selected = "Infliximab"),
            sliderInput("mh_mayo0", "Baseline Mayo:", min = 3, max = 12, value = 9, step = 1)
          ),
          box(title = "Mucosal Healing Index Over Time", width = 9, status = "success",
            plotOutput("mh_plot", height = "300px")
          )
        ),
        fluidRow(
          box(title = "Endoscopic Activity Proxy", width = 6, status = "warning",
            plotOutput("endo_plot", height = "270px")
          ),
          box(title = "Histological Activity Proxy (Geboes)", width = 6, status = "warning",
            plotOutput("histo_plot", height = "270px")
          )
        ),
        fluidRow(
          box(title = "Deep Remission Over Time (Mayo <=2 + MH >=0.7)", width = 12, status = "success",
            plotOutput("deep_rem_plot", height = "280px")
          )
        )
      ),

      # ====================================================
      # TAB 6: SCENARIO COMPARISON
      # ====================================================
      tabItem(tabName = "comparison",
        fluidRow(
          box(title = "Options", width = 3, status = "primary",
            checkboxGroupInput("comp_drugs", "Treatments:",
              choices  = DRUGS,
              selected = DRUGS),
            selectInput("comp_ep", "Endpoint:",
              choices = c(
                "Mayo Score"           = "Mayo_total",
                "Mucosal Healing Index"= "MH_index",
                "CRP (relative)"       = "CRP_val",
                "Fecal Calprotectin"   = "FC_val",
                "TNF-alpha (nM)"       = "TNFa",
                "IL-13 (nM)"           = "IL13",
                "IL-17A (nM)"          = "IL17"
              ),
              selected = "Mayo_total"),
            sliderInput("comp_weeks", "Duration (weeks):", min = 8, max = 52, value = 52, step = 4),
            sliderInput("comp_mayo0", "Baseline Mayo:", min = 3, max = 12, value = 9, step = 1),
            actionButton("run_comp", "Update", icon = icon("sync"), class = "btn-primary btn-block")
          ),
          box(title = "Comparative Endpoint Time-Course", width = 9, status = "primary",
            plotOutput("comp_plot", height = "360px")
          )
        ),
        fluidRow(
          box(title = "Response at Week 8", width = 6, status = "info",
            DTOutput("comp_wk8")
          ),
          box(title = "Response at Week 52", width = 6, status = "success",
            DTOutput("comp_wk52")
          )
        )
      ),

      # ====================================================
      # TAB 7: BIOMARKERS & SAFETY
      # ====================================================
      tabItem(tabName = "safety",
        fluidRow(
          box(title = "Options", width = 3, status = "info",
            selectInput("saf_drug",  "Treatment:", choices = DRUGS, selected = "Tofacitinib"),
            sliderInput("saf_mayo0", "Baseline Mayo:", min = 3, max = 12, value = 9, step = 1)
          ),
          box(title = "CRP & Fecal Calprotectin Kinetics", width = 9, status = "primary",
            plotOutput("saf_biom", height = "300px")
          )
        ),
        fluidRow(
          box(title = "Neutrophil Change vs Placebo (Infection Risk Proxy)", width = 6, status = "warning",
            plotOutput("neut_plot", height = "270px"),
            p(em("JAK inhibitors reduce neutrophil activity; monitor for serious infections."))
          ),
          box(title = "Anemia Risk Index", width = 6, status = "warning",
            plotOutput("anemia_plot", height = "270px"),
            p(em("Driven by mucosal bleeding and hepcidin-mediated iron restriction."))
          )
        ),
        fluidRow(
          box(title = "Safety Considerations by Drug Class", width = 12, status = "danger",
            DTOutput("safety_table")
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

  # ---- Single-drug simulation (Profile tab) ----
  sim_single <- eventReactive(input$run_sim, {
    withProgress(message = paste("Simulating", input$drug_profile, "..."), {
      run_one(UC_MOD, input$drug_profile, input$mayo0)
    })
  }, ignoreNULL = FALSE)

  # ---- All-drugs simulation (Comparison tab) ----
  sim_comp <- reactiveVal(NULL)

  observe({
    if (is.null(sim_comp())) {
      withProgress(message = "Initialising model...", {
        res <- lapply(DRUGS, function(d) run_one(UC_MOD, d, 9))
        sim_comp(bind_rows(Filter(Negate(is.null), res)))
      })
    }
  })

  observeEvent(input$run_comp, {
    withProgress(message = "Running comparison...", value = 0, {
      drugs <- input$comp_drugs
      res <- lapply(seq_along(drugs), function(i) {
        incProgress(1/length(drugs), detail = drugs[i])
        run_one(UC_MOD, drugs[i], input$comp_mayo0, tend = input$comp_weeks * 7)
      })
      sim_comp(bind_rows(Filter(Negate(is.null), res)))
    })
  })

  # ============================================================
  # TAB 1 outputs
  # ============================================================
  output$vbox_mayo <- renderValueBox({
    col <- if (input$mayo0 >= 9) "red" else if (input$mayo0 >= 6) "orange" else "yellow"
    valueBox(input$mayo0, "Mayo Score",
             icon = icon("thermometer-three-quarters"), color = col)
  })
  output$vbox_crp <- renderValueBox({
    valueBox(paste0(input$crp0," mg/L"), "CRP",
             icon = icon("vial"),
             color = if (input$crp0 > 10) "red" else "green")
  })
  output$vbox_fc <- renderValueBox({
    valueBox(paste0(input$fc0," ug/g"), "Fecal Calprotectin",
             icon = icon("flask"),
             color = if (input$fc0 > 250) "red" else "green")
  })
  output$patient_summary <- renderTable({
    sev <- if (input$mayo0 >= 9) "Severe" else if (input$mayo0 >= 6) "Moderate" else "Mild"
    data.frame(
      Parameter = c("Mayo Score","CRP","Fecal Calprotectin","Disease Duration","Body Weight","Severity"),
      Value     = c(input$mayo0, paste0(input$crp0," mg/L"),
                    paste0(input$fc0," ug/g"), paste0(input$dur_yr," yr"),
                    paste0(input$bw," kg"), sev)
    )
  }, striped = TRUE)

  output$mayo_class_table <- renderTable({
    data.frame(
      Category   = c("Remission","Mild","Moderate","Severe"),
      Mayo       = c("0-2","3-5","6-10","11-12"),
      Description= c("Asymptomatic","Mild symptoms","Moderate disease","Fulminant")
    )
  }, striped = TRUE)

  # ============================================================
  # TAB 2 outputs
  # ============================================================
  output$pk_table <- renderDT({
    datatable(PK_PARAMS, rownames = FALSE,
              options = list(dom = 't', pageLength = 10, scrollX = TRUE),
              caption = "Population PK parameters and clinical trial sources")
  })

  pk_sim <- reactive({
    run_one(UC_MOD, input$pk_drug, tend = input$pk_weeks * 7)
  })

  get_conc_col <- function(drug) {
    switch(drug,
      Infliximab  = list(col = "IFX_conc",    unit = "ug/mL"),
      Vedolizumab = list(col = "VDZ_conc",    unit = "ug/mL"),
      Tofacitinib = list(col = "TOF_conc_t",  unit = "ng/mL"),
      Ustekinumab = list(col = "UST_conc",    unit = "ug/mL"),
      Ozanimod    = list(col = "OZA_M1_conc", unit = "ng/mL (CC112273)")
    )
  }

  output$pk_plot <- renderPlot({
    sim  <- pk_sim(); if (is.null(sim)) return(NULL)
    cc   <- get_conc_col(input$pk_drug)
    ggplot(sim, aes(x = time/7, y = .data[[cc$col]])) +
      geom_line(color = SCEN_COLORS[input$pk_drug], linewidth = 1.2) +
      scale_x_continuous(breaks = seq(0, input$pk_weeks, 8), name = "Time (weeks)") +
      labs(title = paste(input$pk_drug, "Concentration-Time Profile"),
           y = paste0("Concentration (", cc$unit, ")")) +
      THEME_UC
  })

  output$trough_plot <- renderPlot({
    sim <- pk_sim(); if (is.null(sim)) return(NULL)
    cc  <- get_conc_col(input$pk_drug)
    trough_times <- switch(input$pk_drug,
      Infliximab  = c(13,27,41,97,153,209,265,321),
      Vedolizumab = c(13,27,41,97,153,209,265,321),
      Tofacitinib = seq(12, input$pk_weeks*7, by = 14),
      Ustekinumab = c(55,111,167,223,279,335),
      Ozanimod    = seq(6, input$pk_weeks*7, by = 7))
    troughs <- sim[sim$time %in% trough_times, ]
    ggplot(sim, aes(x = time/7, y = .data[[cc$col]])) +
      geom_line(color = SCEN_COLORS[input$pk_drug], linewidth = 0.9, alpha = 0.7) +
      geom_point(data = troughs, aes(x = time/7, y = .data[[cc$col]]),
                 color = "red", size = 2.5) +
      scale_x_continuous(breaks = seq(0, input$pk_weeks, 8), name = "Time (weeks)") +
      labs(title = "Trough Concentrations (red dots)",
           y = paste0("Conc (", cc$unit, ")")) +
      THEME_UC
  })

  output$emax_plot <- renderPlot({
    emax <- switch(input$pk_drug, Infliximab=0.92, Vedolizumab=0.85,
                   Tofacitinib=0.88, Ustekinumab=0.80, Ozanimod=0.75)
    ec50 <- switch(input$pk_drug, Infliximab=1.5,  Vedolizumab=10.0,
                   Tofacitinib=50.0, Ustekinumab=1.0, Ozanimod=0.5)
    xmax <- switch(input$pk_drug, Infliximab=20, Vedolizumab=60,
                   Tofacitinib=200, Ustekinumab=12, Ozanimod=5)
    unit <- switch(input$pk_drug, Infliximab="ug/mL", Vedolizumab="ug/mL",
                   Tofacitinib="ng/mL", Ustekinumab="ug/mL", Ozanimod="ng/mL")
    cx <- seq(0, xmax, length.out = 300)
    df <- data.frame(conc = cx, effect = emax * cx^2 / (ec50^2 + cx^2))
    ggplot(df, aes(x = conc, y = effect)) +
      geom_line(color = SCEN_COLORS[input$pk_drug], linewidth = 1.4) +
      geom_vline(xintercept = ec50, linetype = "dashed", color = "red") +
      geom_hline(yintercept = emax/2, linetype = "dashed", color = "red") +
      annotate("text", x = ec50 * 1.1, y = 0.05,
               label = paste0("EC50 = ", ec50, " ", unit),
               color = "red", hjust = 0, size = 3.5) +
      labs(title = paste(input$pk_drug, "Emax Model"),
           x = paste0("Concentration (", unit, ")"),
           y = "Fractional Inhibition") +
      THEME_UC
  })

  # ============================================================
  # TAB 3 outputs
  # ============================================================
  cyto_sim <- reactive({
    run_one(UC_MOD, input$cyto_drug, input$cyto_mayo0, tend = input$cyto_weeks * 7)
  })

  output$cyto_plot <- renderPlot({
    sim <- cyto_sim(); if (is.null(sim)) return(NULL)
    long <- sim %>%
      select(time, TNFa, IL17, IL13, IL10) %>%
      pivot_longer(-time, names_to = "cytokine", values_to = "nM") %>%
      mutate(cytokine = factor(cytokine,
        levels = c("TNFa","IL13","IL17","IL10"),
        labels = c("TNF-alpha","IL-13","IL-17A","IL-10")))
    ggplot(long, aes(x = time/7, y = nM, color = cytokine)) +
      geom_line(linewidth = 1.0) +
      facet_wrap(~ cytokine, scales = "free_y", nrow = 2) +
      scale_x_continuous(breaks = seq(0, input$cyto_weeks, 8), name = "Time (weeks)") +
      scale_color_manual(values = c("TNF-alpha"="#E41A1C","IL-13"="#FF7F00",
                                    "IL-17A"="#4DAF4A","IL-10"="#377EB8")) +
      labs(title = paste("Cytokine Dynamics -", input$cyto_drug), y = "Concentration (nM)") +
      THEME_UC
  })

  output$immune_plot <- renderPlot({
    sim <- cyto_sim(); if (is.null(sim)) return(NULL)
    long <- sim %>%
      select(time, Th2, Th17, Treg, Neutrophil) %>%
      pivot_longer(-time, names_to = "cell", values_to = "count") %>%
      mutate(cell = factor(cell, levels = c("Th2","Th17","Treg","Neutrophil")))
    ggplot(long, aes(x = time/7, y = count, color = cell)) +
      geom_line(linewidth = 1.0) +
      facet_wrap(~ cell, scales = "free_y", nrow = 2) +
      scale_x_continuous(breaks = seq(0, input$cyto_weeks, 8), name = "Time (weeks)") +
      scale_color_manual(values = c(Th2="#E41A1C", Th17="#FF7F00",
                                    Treg="#4DAF4A", Neutrophil="#984EA3")) +
      labs(title = paste("Immune Cell Dynamics -", input$cyto_drug), y = "Count / Index") +
      THEME_UC
  })

  # ============================================================
  # TAB 4 outputs
  # ============================================================
  da_sim <- reactive({
    if (input$da_all)
      bind_rows(lapply(DRUGS, function(d) run_one(UC_MOD, d, input$da_mayo0)))
    else
      run_one(UC_MOD, input$da_drug, input$da_mayo0)
  })

  output$mayo_plot <- renderPlot({
    sim <- da_sim(); if (is.null(sim)) return(NULL)
    use_color <- input$da_all && "drug" %in% names(sim)
    p <- ggplot(sim, aes(x = time/7, y = Mayo_total,
                          color   = if (use_color) drug else NULL,
                          group   = if (use_color) drug else 1)) +
      geom_line(linewidth = 1.1) +
      geom_hline(yintercept = 2, linetype = "dashed") +
      annotate("text", x = max(sim$time)/7*0.92, y = 2.4,
               label = "Remission (Mayo <= 2)", size = 3) +
      scale_x_continuous(breaks = seq(0, max(sim$time)/7, 8), name = "Time (weeks)") +
      scale_y_continuous(limits = c(0, 12), name = "Mayo Score") +
      labs(title = if (input$da_all) "All Treatments" else paste("Mayo Score -", input$da_drug)) +
      THEME_UC
    if (use_color) p <- p + scale_color_manual(values = SCEN_COLORS)
    p
  })

  output$pmayo_plot <- renderPlot({
    sim <- da_sim(); if (is.null(sim)) return(NULL)
    use_color <- input$da_all && "drug" %in% names(sim)
    p <- ggplot(sim, aes(x = time/7, y = Mayo_total * 0.75,
                          color = if (use_color) drug else NULL,
                          group = if (use_color) drug else 1)) +
      geom_line(linewidth = 0.9) +
      geom_hline(yintercept = 2, linetype = "dashed") +
      scale_x_continuous(breaks = seq(0, max(sim$time)/7, 8), name = "Time (weeks)") +
      labs(title = "Partial Mayo Score", y = "Partial Mayo (0-9)") + THEME_UC
    if (use_color) p <- p + scale_color_manual(values = SCEN_COLORS)
    p
  })

  output$biom_da_plot <- renderPlot({
    sim <- da_sim(); if (is.null(sim)) return(NULL)
    use_color <- input$da_all && "drug" %in% names(sim)
    long <- sim %>%
      select(time, any_of("drug"), CRP_val, FC_val) %>%
      pivot_longer(c(CRP_val, FC_val), names_to = "bm", values_to = "val") %>%
      mutate(bm = factor(bm, labels = c("CRP (relative)","FC (ug/g)")))
    p <- ggplot(long, aes(x = time/7, y = val,
                           color    = if (use_color && "drug" %in% names(long)) drug else bm,
                           linetype = bm,
                           group    = if (use_color && "drug" %in% names(long))
                                        interaction(drug, bm) else bm)) +
      geom_line(linewidth = 0.9) +
      facet_wrap(~ bm, scales = "free_y") +
      scale_x_continuous(breaks = seq(0, max(sim$time)/7, 8), name = "Time (weeks)") +
      labs(title = "Biomarkers", y = "Value") + THEME_UC
    if (use_color && "drug" %in% names(long))
      p <- p + scale_color_manual(values = SCEN_COLORS)
    p
  })

  output$response_table <- renderTable({
    sim <- da_sim(); if (is.null(sim)) return(NULL)
    drug_col <- if ("drug" %in% names(sim)) "drug" else NULL
    sim %>%
      filter(time %in% c(0, 56, 112, 365)) %>%
      mutate(
        Week = case_when(time==0~"Baseline", time==56~"Week 8",
                         time==112~"Week 16", time==365~"Week 52"),
        Drug = if (!is.null(drug_col)) .data[[drug_col]] else input$da_drug
      ) %>%
      group_by(Drug, Week) %>%
      summarise(
        `Mayo Score`           = round(mean(Mayo_total), 1),
        `Clin. Remission (%)`  = round(mean(clin_remission)*100, 1),
        `Mucosal Healing (%)` = round(mean(MH_resp)*100, 1),
        `Deep Remission (%)`   = round(mean(deep_remission)*100, 1),
        `CRP (rel.)`           = round(mean(CRP_val), 2),
        .groups = "drop"
      )
  }, striped = TRUE)

  # ============================================================
  # TAB 5 outputs
  # ============================================================
  mh_sim <- reactive({
    run_one(UC_MOD, input$mh_drug, input$mh_mayo0)
  })

  output$mh_plot <- renderPlot({
    sim <- mh_sim(); if (is.null(sim)) return(NULL)
    ggplot(sim, aes(x = time/7, y = MH_index)) +
      geom_line(color = SCEN_COLORS[input$mh_drug], linewidth = 1.2) +
      geom_hline(yintercept = 0.7, linetype = "dashed", color = "darkgreen") +
      annotate("text", x = max(sim$time)/7*0.9, y = 0.73,
               label = "MH threshold (>=0.7)", size = 3, color = "darkgreen") +
      scale_x_continuous(breaks = seq(0, 52, 8), name = "Time (weeks)") +
      scale_y_continuous(limits = c(0, 1), name = "Mucosal Healing Index") +
      labs(title = paste("Mucosal Healing -", input$mh_drug)) + THEME_UC
  })

  output$endo_plot <- renderPlot({
    sim <- mh_sim(); if (is.null(sim)) return(NULL)
    ggplot(sim, aes(x = time/7, y = 1 - MH_index)) +
      geom_line(color = "#E41A1C", linewidth = 1.0) +
      geom_hline(yintercept = 0.3, linetype = "dashed", color = "darkgreen") +
      scale_x_continuous(breaks = seq(0, 52, 8), name = "Time (weeks)") +
      labs(title = "Endoscopic Activity Proxy (1 - MH Index)", y = "Activity") + THEME_UC
  })

  output$histo_plot <- renderPlot({
    sim <- mh_sim(); if (is.null(sim)) return(NULL)
    hprox <- pmax(0, (1 - sim$MH_index) * (sim$Neutrophil / 50))
    df <- data.frame(time = sim$time, hprox = hprox)
    ggplot(df, aes(x = time/7, y = hprox)) +
      geom_line(color = "#984EA3", linewidth = 1.0) +
      geom_hline(yintercept = 0.2, linetype = "dashed", color = "darkgreen") +
      scale_x_continuous(breaks = seq(0, 52, 8), name = "Time (weeks)") +
      labs(title = "Histological Activity Proxy", y = "Geboes Proxy Score") + THEME_UC
  })

  output$deep_rem_plot <- renderPlot({
    sim <- mh_sim(); if (is.null(sim)) return(NULL)
    ggplot(sim, aes(x = time/7, y = deep_remission)) +
      geom_line(color = SCEN_COLORS[input$mh_drug], linewidth = 1.2) +
      geom_hline(yintercept = 0.5, linetype = "dashed") +
      scale_x_continuous(breaks = seq(0, 52, 8), name = "Time (weeks)") +
      scale_y_continuous(limits = c(0, 1), name = "Deep Remission (0/1)") +
      labs(title = paste("Deep Remission -", input$mh_drug)) + THEME_UC
  })

  # ============================================================
  # TAB 6 outputs
  # ============================================================
  output$comp_plot <- renderPlot({
    sim <- sim_comp(); if (is.null(sim)) return(NULL)
    ep    <- input$comp_ep
    drugs <- input$comp_drugs
    sim_f <- sim[sim$drug %in% drugs, ]
    ggplot(sim_f, aes(x = time/7, y = .data[[ep]], color = drug)) +
      geom_line(linewidth = 1.1) +
      scale_x_continuous(breaks = seq(0, max(sim_f$time)/7, 8), name = "Time (weeks)") +
      scale_color_manual(values = SCEN_COLORS) +
      labs(title = paste("Scenario Comparison:", ep), y = ep) + THEME_UC
  })

  make_resp_df <- function(sim, day) {
    if (is.null(sim)) return(NULL)
    sim[sim$time == day, ] %>%
      group_by(drug) %>%
      summarise(
        `Mayo Score`          = round(mean(Mayo_total), 1),
        `Clin. Remission (%)` = round(mean(clin_remission)*100, 1),
        `Mucosal Healing (%)`= round(mean(MH_resp)*100, 1),
        `Deep Remission (%)`  = round(mean(deep_remission)*100, 1),
        .groups = "drop"
      )
  }
  output$comp_wk8  <- renderDT(datatable(make_resp_df(sim_comp(), 56),
    rownames=FALSE, options=list(dom='t', pageLength=10)))
  output$comp_wk52 <- renderDT(datatable(make_resp_df(sim_comp(), 365),
    rownames=FALSE, options=list(dom='t', pageLength=10)))

  # ============================================================
  # TAB 7 outputs
  # ============================================================
  saf_sim <- reactive({
    run_one(UC_MOD, input$saf_drug, input$saf_mayo0)
  })

  output$saf_biom <- renderPlot({
    sim <- saf_sim(); if (is.null(sim)) return(NULL)
    long <- sim %>%
      select(time, CRP_val, FC_val) %>%
      pivot_longer(-time, names_to = "bm", values_to = "val") %>%
      mutate(bm = factor(bm, labels = c("CRP (relative)","FC (ug/g)")))
    ggplot(long, aes(x = time/7, y = val, color = bm)) +
      geom_line(linewidth = 1.1) +
      facet_wrap(~ bm, scales = "free_y") +
      scale_x_continuous(breaks = seq(0, 52, 8), name = "Time (weeks)") +
      scale_color_manual(values = c("#E41A1C","#377EB8")) +
      labs(title = paste("Biomarker Kinetics -", input$saf_drug), y = "Value") + THEME_UC
  })

  output$neut_plot <- renderPlot({
    sim_d <- saf_sim()
    sim_p <- run_one(UC_MOD, "Placebo", input$saf_mayo0)
    if (is.null(sim_d) || is.null(sim_p)) return(NULL)
    df <- data.frame(
      time       = sim_d$time,
      pct_change = (sim_d$Neutrophil - sim_p$Neutrophil) / pmax(sim_p$Neutrophil, 0.1) * 100
    )
    ggplot(df, aes(x = time/7, y = pct_change)) +
      geom_line(color = SCEN_COLORS[input$saf_drug], linewidth = 1.1) +
      geom_hline(yintercept = 0, linetype = "dashed") +
      geom_hline(yintercept = -30, linetype = "dotted", color = "red") +
      annotate("text", x = 50, y = -28, label = "Clinically significant",
               size = 3, color = "red", hjust = 1) +
      scale_x_continuous(breaks = seq(0, 52, 8), name = "Time (weeks)") +
      labs(title = paste("Neutrophil Change vs Placebo -", input$saf_drug),
           y = "% Change") + THEME_UC
  })

  output$anemia_plot <- renderPlot({
    sim <- saf_sim(); if (is.null(sim)) return(NULL)
    df  <- data.frame(
      time  = sim$time,
      arisk = pmax(0, sim$Mayo_total/12*2 - sim$MH_index)
    )
    ggplot(df, aes(x = time/7, y = arisk)) +
      geom_area(fill = SCEN_COLORS[input$saf_drug], alpha = 0.3) +
      geom_line(color = SCEN_COLORS[input$saf_drug], linewidth = 1.0) +
      scale_x_continuous(breaks = seq(0, 52, 8), name = "Time (weeks)") +
      labs(title = paste("Anemia Risk Index -", input$saf_drug),
           y = "Risk Index (higher = worse)") + THEME_UC
  })

  output$safety_table <- renderDT({
    df <- data.frame(
      Drug              = c("Infliximab","Vedolizumab","Tofacitinib","Ustekinumab","Ozanimod"),
      `Infection Risk`  = c("Moderate (TB, fungal)","Low (gut-selective)",
                             "Moderate-High (HZV, URI)","Low-Moderate","Low"),
      Malignancy        = c("Lymphoma monitoring","Low","NMSC reported","Low","Low"),
      Cardiovascular    = c("Low","Low","MACE, DVT/PE (>50yr)","Low","Bradycardia (1st dose)"),
      Other             = c("Infusion rxn, demyelination","Nasopharyngitis",
                             "Lipid elevation, anaemia","Injection site","Macular oedema"),
      Monitoring        = c("TB, LFT, CBC","CBC, LFT","CBC, lipids, hep B","CBC, LFT","Ophthalmology, ECG"),
      check.names = FALSE, stringsAsFactors = FALSE
    )
    datatable(df, rownames = FALSE,
              options = list(dom = 't', scrollX = TRUE, pageLength = 10))
  })
}

# ============================================================
# LAUNCH
# ============================================================
shinyApp(ui = ui, server = server)
