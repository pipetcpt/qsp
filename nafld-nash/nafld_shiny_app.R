## ============================================================
## NAFLD/NASH QSP Interactive Dashboard — Shiny App
## Non-Alcoholic Fatty Liver Disease / Steatohepatitis
##
## Tabs:
##   1. Patient Profile & Disease Baseline
##   2. Drug PK Profiles
##   3. Hepatic Endpoints (Liver fat, NAS, ALT)
##   4. Fibrosis & Inflammation
##   5. Metabolic Biomarkers
##   6. Scenario Comparison
##   7. Biomarker Tracker
## ============================================================

library(shiny)
library(bslib)
library(mrgsolve)
library(dplyr)
library(ggplot2)
library(tidyr)
library(DT)
library(plotly)

## ── mrgsolve model (embedded) ─────────────────────────────
nafld_code <- '
$PARAM
KA_RSM=0.8, CL_RSM=4.5, V1_RSM=25, Q_RSM=8.0, V2_RSM=60,
FTUP_RSM=0.70, EC50_RSM=0.08, EMAX_RSM=0.85,
KA_OCA=1.2, CL_OCA=15, V1_OCA=20, EC50_OCA=0.05, EMAX_OCA=0.75,
KA_SEM=0.0045, CL_SEM=0.034, V1_SEM=3.2, EC50_SEM=0.025, EMAX_SEM=0.80,
KA_EMP=1.5, CL_EMP=10, V1_EMP=70, EC50_EMP=0.15, EMAX_EMP=0.90,
KLIN_LF=0.003, KOUT_LF=0.003, LF0=0.20,
DNL0=1.0, KDNL_IR=0.4,
IR0=2.5, KOUT_IR=0.02, KFFA_IR=0.15, KTNF_IR=0.25,
KUP0=0.5, KOUT_KUP=0.05, KLIP_KUP=0.2,
KOUT_TNF=0.12, KPROD_TNF=0.06, KKUP_TNF=0.5, TNF0=0.5,
KOUT_IL6=0.15, KPROD_IL6=0.05, KKUP_IL6=0.3, IL60=0.33,
KOUT_TGF=0.08, KPROD_TGF=0.012, KKUP_TGF=0.15, KLIP_TGF=0.10, TGF0=0.15,
KOUT_HSC=0.005, KTGF_HSC=0.15, HSC0=0.10,
KOUT_COL=0.003, KHSC_COL=0.05, COL0=0.15, KFIBREG=0.8,
KOUT_ALT=0.030, KREL_ALT=1.5, ALT0=45, KAPOP_ALT=8.0,
ADIPON0=6.0, KOUT_ADI=0.05, FXR0=0.5, LIPOTOX0=0.3,
WT0=95, KOUT_WT=0.0015,
DOSE_RSM=0, DOSE_OCA=0, DOSE_SEM=0, DOSE_EMP=0

$CMT
RSM_GUT RSM_CENT RSM_PERI RSM_LIVER
OCA_GUT OCA_CENT
SEM_SC SEM_CENT
EMP_GUT EMP_CENT
LIVER_FAT INS_RES KUPFFER TNFA IL6C TGFB HSC COLLAGEN ALT_CMT ADIPONECTIN BODY_WT

$MAIN
LIVER_FAT_0=0.20; INS_RES_0=IR0; KUPFFER_0=KUP0;
TNFA_0=TNF0; IL6C_0=IL60; TGFB_0=TGF0;
HSC_0=HSC0; COLLAGEN_0=COL0; ALT_CMT_0=ALT0;
ADIPONECTIN_0=ADIPON0; BODY_WT_0=WT0;

$ODE
double Cp_RSM=RSM_CENT/V1_RSM;
double E_RSM=EMAX_RSM*FTUP_RSM*Cp_RSM/(EC50_RSM+FTUP_RSM*Cp_RSM);
double Cp_OCA=OCA_CENT/V1_OCA;
double E_OCA=EMAX_OCA*Cp_OCA/(EC50_OCA+Cp_OCA);
double Cp_SEM=SEM_CENT/V1_SEM;
double E_SEM=EMAX_SEM*Cp_SEM/(EC50_SEM+Cp_SEM);
double Cp_EMP=EMP_CENT/V1_EMP;
double E_EMP=EMAX_EMP*Cp_EMP/(EC50_EMP+Cp_EMP);

dxdt_RSM_GUT=-KA_RSM*RSM_GUT;
dxdt_RSM_CENT=KA_RSM*RSM_GUT-(CL_RSM+Q_RSM)/V1_RSM*RSM_CENT+Q_RSM/V2_RSM*RSM_PERI;
dxdt_RSM_PERI=Q_RSM/V1_RSM*RSM_CENT-Q_RSM/V2_RSM*RSM_PERI;
dxdt_RSM_LIVER=FTUP_RSM*CL_RSM/V1_RSM*RSM_CENT-0.5*RSM_LIVER;

dxdt_OCA_GUT=-KA_OCA*OCA_GUT;
dxdt_OCA_CENT=KA_OCA*OCA_GUT-CL_OCA/V1_OCA*OCA_CENT;

dxdt_SEM_SC=-KA_SEM*SEM_SC;
dxdt_SEM_CENT=KA_SEM*SEM_SC-CL_SEM/V1_SEM*SEM_CENT;

dxdt_EMP_GUT=-KA_EMP*EMP_GUT;
dxdt_EMP_CENT=KA_EMP*EMP_GUT-CL_EMP/V1_EMP*EMP_CENT;

double WT_target=WT0*(1-0.15*E_SEM-0.05*E_EMP);
dxdt_BODY_WT=KOUT_WT*(WT_target-BODY_WT);

double ADIPON_ss=ADIPON0*(WT0/BODY_WT)*(1+0.3*E_SEM);
dxdt_ADIPONECTIN=KOUT_ADI*(ADIPON_ss-ADIPONECTIN);

double IR_ss=IR0*(1+KFFA_IR*(LIVER_FAT/LF0-1)+KTNF_IR*(TNFA/TNF0-1)-0.2*(ADIPONECTIN/ADIPON0-1))*(1-0.3*E_SEM-0.2*E_EMP-0.1*E_RSM);
IR_ss=(IR_ss<0.5)?0.5:IR_ss;
dxdt_INS_RES=KOUT_IR*(IR_ss-INS_RES);

double DNL=DNL0*(1+KDNL_IR*(INS_RES/IR0-1))*(1-0.4*E_RSM-0.25*E_OCA);
double LF_influx=KLIN_LF*DNL*(BODY_WT/WT0);
double LF_efflux=KOUT_LF*LIVER_FAT*(1+0.6*E_RSM)*(1+0.15*E_SEM);
dxdt_LIVER_FAT=LF_influx-LF_efflux;

double LIPOTOX=LIPOTOX0*(LIVER_FAT/LF0);
dxdt_KUPFFER=(KOUT_KUP*KUP0+KLIP_KUP*LIPOTOX/LIPOTOX0)/(1+0.3*E_OCA+0.2*E_SEM)-KOUT_KUP*KUPFFER;
dxdt_TNFA=KPROD_TNF+KKUP_TNF*KUPFFER-KOUT_TNF*TNFA;
dxdt_IL6C=KPROD_IL6+KKUP_IL6*KUPFFER-KOUT_IL6*IL6C;
dxdt_TGFB=(KPROD_TGF+KKUP_TGF*KUPFFER+KLIP_TGF*LIPOTOX/LIPOTOX0)/(1+0.4*E_OCA+0.2*E_RSM)-KOUT_TGF*TGFB;
dxdt_HSC=(KOUT_HSC*HSC0+KTGF_HSC*TGFB)/(1+0.4*E_OCA+0.15*E_RSM+0.1*E_SEM)-KOUT_HSC*HSC;
dxdt_COLLAGEN=(KOUT_COL*COL0+KHSC_COL*HSC)/(1+0.35*E_OCA+0.15*E_RSM)-KOUT_COL*COLLAGEN;
dxdt_ALT_CMT=KREL_ALT+KAPOP_ALT*TNFA*LIPOTOX-KOUT_ALT*ALT_CMT;

$TABLE
double LF_PCT=LIVER_FAT*100;
double FIB_SCORE=KFIBREG*COLLAGEN;
FIB_SCORE=(FIB_SCORE>4)?4:FIB_SCORE;
double NAS=0;
NAS+=(LF_PCT>=5)?((LF_PCT<33)?1:(LF_PCT<66)?2:3):0;
NAS+=(ALT_CMT>40)?1:0;
NAS+=(TNFA>TNF0*1.5)?1:0;
NAS=(NAS>8)?8:NAS;
double TG_SERUM=150*(LIVER_FAT/LF0)*(INS_RES/IR0)*(1-0.5*EMAX_RSM*FTUP_RSM*(RSM_CENT/V1_RSM)/(EC50_RSM+FTUP_RSM*(RSM_CENT/V1_RSM)));
double LDL_C=120*(1+0.1*(LIVER_FAT/LF0-1))*(1-0.25*E_RSM);
double PDFF=LF_PCT*0.85;
double FIB4=1.8*(COLLAGEN/COL0)*(ALT_CMT/ALT0);
double Cp_RSM_out=RSM_CENT/V1_RSM;
double Cp_OCA_out=OCA_CENT/V1_OCA;
double Cp_SEM_out=SEM_CENT/V1_SEM;
double Cp_EMP_out=EMP_CENT/V1_EMP;
double Week=time/168;

$CAPTURE
LF_PCT PDFF FIB_SCORE NAS ALT_CMT TNFA IL6C TGFB HSC COLLAGEN
KUPFFER INS_RES BODY_WT ADIPONECTIN TG_SERUM LDL_C FIB4
Cp_RSM_out Cp_OCA_out Cp_SEM_out Cp_EMP_out Week
'

## Compile once at startup
mod_global <- mcode("NAFLD_Shiny", nafld_code)

## ── UI ─────────────────────────────────────────────────────
ui <- page_navbar(
  title = div(
    span("NAFLD/NASH", style = "font-weight:bold; color:#880E4F;"),
    span(" QSP Dashboard", style = "color:#333;")
  ),
  theme = bs_theme(
    bootswatch  = "flatly",
    primary     = "#880E4F",
    success     = "#388E3C",
    info        = "#1565C0"
  ),
  bg = "#F8F9FA",

  ## ── Sidebar (shared across tabs) ─────────────────────────
  nav_panel(
    "Patient Profile",
    icon = icon("user-md"),
    layout_sidebar(
      sidebar = sidebar(
        width = 300,
        h5("Patient Characteristics", style = "color:#880E4F; font-weight:bold;"),
        sliderInput("bmi",   "BMI (kg/m²)",    min = 22, max = 50, value = 34, step = 0.5),
        sliderInput("age",   "Age (years)",     min = 18, max = 75, value = 48, step = 1),
        selectInput("sex",   "Sex",             choices = c("Male", "Female"), selected = "Male"),
        selectInput("t2dm",  "T2DM Comorbidity",choices = c("Yes", "No"), selected = "Yes"),
        selectInput("fib_stage", "Fibrosis Stage (baseline)",
                    choices = c("F0", "F1", "F2", "F3"), selected = "F2"),
        hr(),
        h5("Disease Activity", style = "color:#1565C0; font-weight:bold;"),
        numericInput("alt_base",  "Baseline ALT (U/L)",   value = 65, min = 10, max = 300),
        numericInput("pdff_base", "Baseline MRI-PDFF (%)",value = 18, min = 5,  max = 60),
        numericInput("homa_base", "Baseline HOMA-IR",     value = 3.2, min = 1,  max = 15, step = 0.1),
        hr(),
        actionButton("update_profile", "Update Patient Profile",
                     class = "btn-primary btn-sm", width = "100%")
      ),
      fluidRow(
        column(4,
          valueBox_shiny("Disease Activity", "NASH Confirmed",  "#E53935", icon("fire")),
          br(),
          valueBox_shiny("Fibrosis Risk",    "High (F2)",        "#F57F17", icon("chart-line"))
        ),
        column(8,
          h4("NAFLD/NASH Disease Overview", style = "margin-top:0;"),
          p("Non-Alcoholic Steatohepatitis (NASH) is characterized by hepatic steatosis,
            lobular inflammation, and hepatocellular ballooning. Advanced fibrosis (F3–F4)
            is the major determinant of long-term liver-related morbidity and mortality."),
          tags$hr(),
          h5("Key Pathophysiological Drivers"),
          tags$ul(
            tags$li(strong("Insulin Resistance:"), " Central driver of hepatic fat accumulation via de novo lipogenesis and adipose FFA flux"),
            tags$li(strong("Lipotoxicity:"), " Ceramide and DAG mediate ER stress, mitochondrial dysfunction, and hepatocyte apoptosis"),
            tags$li(strong("Gut-Liver Axis:"), " Dysbiosis → LPS → TLR4 → NF-κB → Kupffer cell activation → TNF-α/IL-1β"),
            tags$li(strong("TGF-β1 / HSC:"), " Master fibrogenic pathway; activated HSCs produce collagen I/III → ECM accumulation"),
            tags$li(strong("NASH Resolution:"), " ≥2-point NAS reduction without fibrosis worsening (FDA primary endpoint)")
          ),
          tags$hr(),
          h5("Treatment Landscape (2024–2025)"),
          tableOutput("treatment_landscape")
        )
      )
    )
  ),

  ## ── Tab 2: Drug PK ─────────────────────────────────────
  nav_panel(
    "Drug PK Profiles",
    icon = icon("pills"),
    layout_sidebar(
      sidebar = sidebar(
        width = 280,
        h5("Drug Selection & Dosing"),
        checkboxGroupInput("drugs", "Active Drugs",
                           choices  = c("Resmetirom 100mg QD"  = "rsm",
                                        "OCA 25mg QD"          = "oca",
                                        "Semaglutide 2.4mg QW" = "sem",
                                        "Empagliflozin 10mg QD"= "emp"),
                           selected = "rsm"),
        hr(),
        h6("Resmetirom PK Parameters"),
        sliderInput("ka_rsm", "ka (1/h)",  min = 0.2, max = 3,   value = 0.8, step = 0.1),
        sliderInput("cl_rsm", "CL (L/h)",  min = 1,   max = 15,  value = 4.5, step = 0.5),
        sliderInput("v1_rsm", "V1 (L)",    min = 5,   max = 80,  value = 25,  step = 5),
        actionButton("run_pk", "Simulate PK", class = "btn-success btn-sm", width = "100%")
      ),
      fluidRow(
        column(6, plotlyOutput("pk_rsm_plot", height = "300px")),
        column(6, plotlyOutput("pk_oca_plot", height = "300px"))
      ),
      fluidRow(
        column(6, plotlyOutput("pk_sem_plot", height = "300px")),
        column(6, plotlyOutput("pk_emp_plot", height = "300px"))
      ),
      fluidRow(
        column(12, DTOutput("pk_summary_table"))
      )
    )
  ),

  ## ── Tab 3: Hepatic Endpoints ───────────────────────────
  nav_panel(
    "Hepatic Endpoints",
    icon = icon("heartbeat"),
    layout_sidebar(
      sidebar = sidebar(
        width = 280,
        h5("Simulation Settings"),
        sliderInput("sim_weeks", "Duration (weeks)", min = 12, max = 104, value = 72, step = 4),
        selectInput("treatment_arm", "Treatment Arm",
                    choices = c("Placebo",
                                "Resmetirom 100mg QD",
                                "OCA 25mg QD",
                                "Semaglutide 2.4mg QW",
                                "Resmetirom + Semaglutide",
                                "Triple Combination")),
        hr(),
        h6("Endpoint Thresholds"),
        numericInput("nash_thresh", "NASH Resolution NAS target (≤)", value = 2, min = 0, max = 4),
        numericInput("pdff_thresh", "Responder PDFF threshold (%)", value = 10, min = 1, max = 30),
        actionButton("run_liver", "Run Simulation", class = "btn-success btn-sm", width = "100%")
      ),
      fluidRow(
        column(6, plotlyOutput("lf_plot",  height = "280px")),
        column(6, plotlyOutput("pdff_plot",height = "280px"))
      ),
      fluidRow(
        column(6, plotlyOutput("nas_plot", height = "280px")),
        column(6, plotlyOutput("alt_plot", height = "280px"))
      ),
      fluidRow(
        column(12,
          h5("Hepatic Endpoint Summary", style = "margin-top:15px;"),
          DTOutput("liver_table")
        )
      )
    )
  ),

  ## ── Tab 4: Fibrosis & Inflammation ────────────────────
  nav_panel(
    "Fibrosis & Inflammation",
    icon = icon("bacteria"),
    layout_sidebar(
      sidebar = sidebar(
        width = 280,
        h5("Fibrosis Parameters"),
        sliderInput("ktgf_hsc", "TGF-β→HSC rate", min = 0.05, max = 0.5, value = 0.15, step = 0.01),
        sliderInput("khsc_col", "HSC→Collagen rate",min = 0.01, max = 0.15, value = 0.05, step = 0.005),
        sliderInput("kout_col", "Collagen turnover (1/h)", min = 0.001, max = 0.01, value = 0.003, step = 0.001),
        hr(),
        h5("Inflammation Parameters"),
        sliderInput("klip_kup", "Lipotox→Kupffer", min = 0.05, max = 0.6, value = 0.2, step = 0.05),
        sliderInput("kkup_tgf", "Kupffer→TGF-β",  min = 0.05, max = 0.5, value = 0.15, step = 0.05),
        actionButton("run_fib", "Recalculate", class = "btn-warning btn-sm", width = "100%")
      ),
      fluidRow(
        column(6, plotlyOutput("fib_plot",  height = "280px")),
        column(6, plotlyOutput("hsc_plot",  height = "280px"))
      ),
      fluidRow(
        column(6, plotlyOutput("tgfb_plot", height = "280px")),
        column(6, plotlyOutput("tnfa_plot", height = "280px"))
      ),
      fluidRow(
        column(6, plotlyOutput("il6_plot",  height = "280px")),
        column(6, plotlyOutput("kup_plot",  height = "280px"))
      )
    )
  ),

  ## ── Tab 5: Metabolic Biomarkers ───────────────────────
  nav_panel(
    "Metabolic Biomarkers",
    icon = icon("flask"),
    layout_sidebar(
      sidebar = sidebar(
        width = 280,
        h5("Metabolic Profile"),
        sliderInput("ir0",   "Baseline HOMA-IR",   min = 1,  max = 12, value = 3.2, step = 0.1),
        sliderInput("wt0",   "Baseline Weight (kg)",min = 60, max = 150, value = 95, step = 1),
        sliderInput("adipon0","Baseline Adiponectin (µg/mL)", min = 2, max = 15, value = 6, step = 0.5),
        hr(),
        h5("Lipid Parameters"),
        numericInput("ldl_base",  "Baseline LDL-C (mg/dL)",  value = 120, min = 60, max = 250),
        numericInput("tg_base",   "Baseline TG (mg/dL)",     value = 185, min = 50, max = 500),
        actionButton("run_meta",  "Update", class = "btn-info btn-sm", width = "100%")
      ),
      fluidRow(
        column(6, plotlyOutput("homa_plot", height = "280px")),
        column(6, plotlyOutput("wt_plot",   height = "280px"))
      ),
      fluidRow(
        column(6, plotlyOutput("tg_plot",   height = "280px")),
        column(6, plotlyOutput("ldl_plot",  height = "280px"))
      ),
      fluidRow(
        column(6, plotlyOutput("adipon_plot",height = "280px")),
        column(6, plotlyOutput("fib4_plot", height = "280px"))
      )
    )
  ),

  ## ── Tab 6: Scenario Comparison ────────────────────────
  nav_panel(
    "Scenario Comparison",
    icon = icon("chart-bar"),
    layout_sidebar(
      sidebar = sidebar(
        width = 280,
        h5("Select Scenarios to Compare"),
        checkboxGroupInput("compare_arms", "Treatment Arms",
                           choices  = c("Placebo",
                                        "Resmetirom 100mg QD",
                                        "OCA 25mg QD",
                                        "Semaglutide 2.4mg QW",
                                        "Resmetirom + Semaglutide",
                                        "Triple Combination"),
                           selected = c("Placebo", "Resmetirom 100mg QD",
                                        "Semaglutide 2.4mg QW",
                                        "Resmetirom + Semaglutide")),
        hr(),
        sliderInput("compare_weeks", "Follow-up (weeks)", min = 24, max = 104, value = 72, step = 4),
        h6("Display Metric"),
        radioButtons("metric_select", NULL,
                     choices = c("Liver Fat (%)" = "LF_PCT",
                                 "Fibrosis Score"= "FIB_SCORE",
                                 "ALT (U/L)"     = "ALT_CMT",
                                 "NAS Score"     = "NAS",
                                 "HOMA-IR"       = "INS_RES",
                                 "TG (mg/dL)"    = "TG_SERUM"),
                     selected = "LF_PCT"),
        actionButton("run_compare", "Compare All", class = "btn-danger btn-sm", width = "100%")
      ),
      fluidRow(
        column(12, plotlyOutput("compare_plot", height = "400px"))
      ),
      fluidRow(
        column(12,
          h5("Week-72 Endpoint Summary Table", style = "margin-top:15px;"),
          DTOutput("compare_table")
        )
      )
    )
  ),

  ## ── Tab 7: Biomarker Tracker ──────────────────────────
  nav_panel(
    "Biomarker Tracker",
    icon = icon("vials"),
    layout_sidebar(
      sidebar = sidebar(
        width = 280,
        h5("Biomarker Panel Selection"),
        checkboxGroupInput("bm_panel", "Biomarkers",
                           choices = c(
                             "MRI-PDFF (%)"           = "PDFF_pct",
                             "Liver Stiffness (kPa)"  = "LSM_kPa",
                             "ALT (U/L)"              = "ALT_CMT",
                             "FIB-4 Index"            = "FIB4",
                             "HOMA-IR"                = "INS_RES",
                             "TG (mg/dL)"             = "TG_SERUM",
                             "LDL-C (mg/dL)"          = "LDL_C",
                             "Adiponectin (µg/mL)"    = "ADIPONECTIN",
                             "TNF-α (rel)"            = "TNFA",
                             "Body Weight (kg)"       = "BODY_WT"
                           ),
                           selected = c("PDFF_pct", "ALT_CMT", "FIB4", "INS_RES", "BODY_WT")),
        hr(),
        selectInput("bm_drug", "Drug for Biomarker Tracking",
                    choices = c("Resmetirom 100mg QD",
                                "OCA 25mg QD",
                                "Semaglutide 2.4mg QW",
                                "Resmetirom + Semaglutide")),
        actionButton("run_bm", "Track Biomarkers", class = "btn-primary btn-sm", width = "100%")
      ),
      fluidRow(
        column(12, plotlyOutput("bm_spider_plot", height = "400px"))
      ),
      fluidRow(
        column(12,
          h5("Change from Baseline at Key Timepoints"),
          DTOutput("bm_change_table")
        )
      )
    )
  )
)

## ── Helper functions ───────────────────────────────────────
valueBox_shiny <- function(title, value, color, icon_tag) {
  div(
    class = "card mb-3",
    style = paste0("background:", color, "; color:white; border-radius:8px; padding:12px;"),
    div(style = "font-size:13px; opacity:0.85;", title),
    div(style = "font-size:20px; font-weight:bold;", value),
    div(style = "font-size:28px; opacity:0.5; float:right; margin-top:-40px;", icon_tag)
  )
}

## Build event schedule for a scenario
get_events <- function(arm, weeks) {
  hrs <- weeks * 7 * 24
  ev_list <- list()
  addl_qd  <- ceiling(hrs / 24) - 1
  addl_qw  <- ceiling(hrs / 168) - 1

  if (arm %in% c("Resmetirom 100mg QD",    "Resmetirom + Semaglutide", "Triple Combination"))
    ev_list[["rsm"]] <- ev(amt = 100, cmt = "RSM_GUT", ii = 24,  addl = addl_qd)
  if (arm %in% c("OCA 25mg QD",            "Triple Combination"))
    ev_list[["oca"]] <- ev(amt = 25,  cmt = "OCA_GUT", ii = 24,  addl = addl_qd)
  if (arm %in% c("Semaglutide 2.4mg QW",   "Resmetirom + Semaglutide", "Triple Combination"))
    ev_list[["sem"]] <- ev(amt = 2.4, cmt = "SEM_SC",  ii = 168, addl = addl_qw)

  if (length(ev_list) == 0) return(ev(time = 0, amt = 0, cmt = "RSM_GUT"))
  Reduce("+", ev_list)
}

## Get param flags for a scenario
get_params <- function(arm) {
  list(
    DOSE_RSM = as.numeric(arm %in% c("Resmetirom 100mg QD",    "Resmetirom + Semaglutide", "Triple Combination")),
    DOSE_OCA = as.numeric(arm %in% c("OCA 25mg QD",            "Triple Combination")),
    DOSE_SEM = as.numeric(arm %in% c("Semaglutide 2.4mg QW",   "Resmetirom + Semaglutide", "Triple Combination")),
    DOSE_EMP = 0
  )
}

## Run a single scenario
sim_arm <- function(mod, arm, weeks = 72) {
  tg <- tgrid(0, weeks * 168, delta = 168)
  ev_arm <- get_events(arm, weeks)
  p_arm  <- get_params(arm)
  out <- mrgsim(do.call(param, c(list(mod), p_arm)),
                ev = ev_arm, tgrid = tg, obsonly = TRUE)
  df <- as.data.frame(out)
  df$Scenario <- arm
  df$PDFF_pct <- df$LF_PCT * 0.85
  df$LSM_kPa  <- 5 + 2 * df$FIB_SCORE   # proxy: 5 kPa baseline + fibrosis
  df
}

SCENARIO_COLORS <- c(
  "Placebo"                  = "#999999",
  "Resmetirom 100mg QD"     = "#E91E8C",
  "OCA 25mg QD"             = "#1565C0",
  "Semaglutide 2.4mg QW"   = "#388E3C",
  "Resmetirom + Semaglutide"= "#F57F17",
  "Triple Combination"      = "#B71C1C"
)

## ── Server ─────────────────────────────────────────────────
server <- function(input, output, session) {

  ## ── Patient Profile Tab ────────────────────────────────
  output$treatment_landscape <- renderTable({
    data.frame(
      Drug          = c("Resmetirom (Rezdiffra)", "OCA (Ocaliva)", "Semaglutide (Ozempic/Wegovy)",
                        "Empagliflozin", "Lanifibranor"),
      Target        = c("THRβ", "FXR", "GLP-1R", "SGLT2", "Pan-PPAR"),
      `FDA Status`  = c("Approved (Mar 2024)", "Phase 3", "Phase 3", "Phase 2b", "Phase 3"),
      `NASH Resolution %` = c("25.9%", "23%", "37%", "~28% (est)", "TBD"),
      stringsAsFactors = FALSE,
      check.names = FALSE
    )
  }, striped = TRUE, hover = TRUE, bordered = TRUE)

  ## ── Reactive: Run all-scenarios simulation ─────────────
  all_sims <- reactive({
    input$run_compare  # trigger
    req(input$compare_arms)
    isolate({
      weeks <- input$compare_weeks
      arms  <- c("Placebo", input$compare_arms)
      arms  <- unique(arms)
      bind_rows(lapply(arms, function(a) sim_arm(mod_global, a, weeks)))
    })
  }) %>% bindEvent(input$run_compare, ignoreNULL = FALSE)

  ## ── Reactive: single arm for liver tab ─────────────────
  liver_sims <- reactive({
    input$run_liver
    isolate({
      weeks <- input$sim_weeks
      arms  <- c("Placebo", input$treatment_arm)
      bind_rows(lapply(unique(arms), function(a) sim_arm(mod_global, a, weeks)))
    })
  }) %>% bindEvent(input$run_liver, ignoreNULL = FALSE)

  ## ── PK plots ───────────────────────────────────────────
  pk_data <- reactive({
    input$run_pk
    isolate({
      ev_rsm <- ev(amt = 100, cmt = "RSM_GUT", ii = 24, addl = 13)
      tg_pk  <- tgrid(0, 360, delta = 0.5)
      sim_rsm <- as.data.frame(
        mrgsim(param(mod_global, DOSE_RSM = 1,
                     KA_RSM = input$ka_rsm,
                     CL_RSM = input$cl_rsm,
                     V1_RSM = input$v1_rsm),
               ev = ev_rsm, tgrid = tg_pk, obsonly = TRUE)
      )
      sim_rsm$Drug <- "Resmetirom"
      sim_rsm
    })
  }) %>% bindEvent(input$run_pk, ignoreNULL = FALSE)

  make_pk_plot <- function(df, col, title, color) {
    p <- plot_ly(df, x = ~time, y = as.formula(paste0("~", col)),
                 type = "scatter", mode = "lines",
                 line = list(color = color, width = 2.5)) %>%
      layout(title = list(text = title, font = list(size = 13)),
             xaxis = list(title = "Time (hours)"),
             yaxis = list(title = "Concentration (µg/L)"),
             showlegend = FALSE)
    p
  }

  output$pk_rsm_plot <- renderPlotly({
    df <- pk_data()
    make_pk_plot(df, "Cp_RSM_out", "Resmetirom PK (100 mg QD)", "#E91E8C")
  })

  output$pk_oca_plot <- renderPlotly({
    ev_oca <- ev(amt = 25, cmt = "OCA_GUT", ii = 24, addl = 6)
    tg_pk  <- tgrid(0, 180, delta = 0.5)
    df_oca <- as.data.frame(mrgsim(param(mod_global, DOSE_OCA = 1),
                                   ev = ev_oca, tgrid = tg_pk, obsonly = TRUE))
    make_pk_plot(df_oca, "Cp_OCA_out", "OCA PK (25 mg QD)", "#1565C0")
  })

  output$pk_sem_plot <- renderPlotly({
    ev_sem <- ev(amt = 2.4, cmt = "SEM_SC", ii = 168, addl = 7)
    tg_pk  <- tgrid(0, 1344, delta = 6)
    df_sem <- as.data.frame(mrgsim(param(mod_global, DOSE_SEM = 1),
                                   ev = ev_sem, tgrid = tg_pk, obsonly = TRUE))
    make_pk_plot(df_sem, "Cp_SEM_out", "Semaglutide PK (2.4 mg QW)", "#388E3C")
  })

  output$pk_emp_plot <- renderPlotly({
    ev_emp <- ev(amt = 10, cmt = "EMP_GUT", ii = 24, addl = 6)
    tg_pk  <- tgrid(0, 180, delta = 0.5)
    df_emp <- as.data.frame(mrgsim(param(mod_global, DOSE_EMP = 1),
                                   ev = ev_emp, tgrid = tg_pk, obsonly = TRUE))
    make_pk_plot(df_emp, "Cp_EMP_out", "Empagliflozin PK (10 mg QD)", "#7B1FA2")
  })

  output$pk_summary_table <- renderDT({
    df <- pk_data()
    ss_df <- df %>% filter(time > 312)
    tibble(
      Parameter = c("Cmax (µg/L)", "Tmax (h)", "AUC₀₋₂₄ (µg·h/L)",
                    "Cmin (µg/L)", "t½ (h, effective)", "Hepatic Uptake (%)"),
      Value     = c(
        round(max(ss_df$Cp_RSM_out, na.rm=TRUE), 2),
        round(ss_df$time[which.max(ss_df$Cp_RSM_out)] %% 24, 1),
        round(sum(diff(ss_df$time[1:49]) * ss_df$Cp_RSM_out[1:48], na.rm=TRUE), 1),
        round(min(ss_df$Cp_RSM_out, na.rm=TRUE), 3),
        round(log(2) / (input$cl_rsm / input$v1_rsm), 1),
        "70 (OATP1B1)"
      )
    )
  }, options = list(dom = "t", pageLength = 10))

  ## ── Liver endpoint plots ───────────────────────────────
  make_plotly <- function(df, y_col, title, ylab, color_map,
                          hline = NULL, hline_label = NULL) {
    p <- plot_ly()
    arms <- unique(df$Scenario)
    for (a in arms) {
      sub <- df %>% filter(Scenario == a)
      col <- color_map[a]
      if (is.na(col)) col <- "#333333"
      p <- add_trace(p, data = sub, x = ~Week, y = as.formula(paste0("~", y_col)),
                     type = "scatter", mode = "lines",
                     name = a, line = list(color = col, width = 2))
    }
    if (!is.null(hline)) {
      p <- add_segments(p, x = 0, xend = max(df$Week), y = hline, yend = hline,
                        line = list(dash = "dot", color = "gray", width = 1),
                        showlegend = FALSE)
    }
    p %>% layout(
      title = list(text = title, font = list(size = 13)),
      xaxis = list(title = "Week"),
      yaxis = list(title = ylab),
      legend = list(orientation = "h", y = -0.25, font = list(size = 10))
    )
  }

  output$lf_plot  <- renderPlotly({
    df <- liver_sims()
    make_plotly(df, "LF_PCT",    "Hepatic Fat Fraction", "Liver Fat (%)", SCENARIO_COLORS, hline = 5)
  })
  output$pdff_plot<- renderPlotly({
    df <- liver_sims()
    make_plotly(df, "PDFF_pct",  "MRI-PDFF",              "PDFF (%)",      SCENARIO_COLORS, hline = 8)
  })
  output$nas_plot <- renderPlotly({
    df <- liver_sims()
    make_plotly(df, "NAS",       "NAFLD Activity Score",  "NAS (0-8)",     SCENARIO_COLORS, hline = 3)
  })
  output$alt_plot <- renderPlotly({
    df <- liver_sims()
    make_plotly(df, "ALT_CMT",   "Serum ALT",             "ALT (U/L)",     SCENARIO_COLORS, hline = 35)
  })

  output$liver_table <- renderDT({
    df <- liver_sims()
    tab <- df %>%
      filter(abs(Week - max(Week)) < 0.1) %>%
      select(Scenario, LF_PCT, PDFF_pct, FIB_SCORE, NAS, ALT_CMT) %>%
      mutate(across(where(is.numeric), ~ round(.x, 1)))
    names(tab) <- c("Scenario", "Liver Fat (%)", "PDFF (%)", "Fibrosis", "NAS", "ALT (U/L)")
    tab
  }, rownames = FALSE, options = list(dom = "t", pageLength = 10))

  ## ── Fibrosis / inflammation plots ─────────────────────
  fib_sims <- reactive({
    input$run_fib
    isolate({
      arms <- c("Placebo", "Resmetirom 100mg QD", "OCA 25mg QD",
                "Semaglutide 2.4mg QW", "Resmetirom + Semaglutide")
      mod_f <- param(mod_global,
                     KTGF_HSC = input$ktgf_hsc,
                     KHSC_COL = input$khsc_col,
                     KOUT_COL = input$kout_col,
                     KLIP_KUP = input$klip_kup,
                     KKUP_TGF = input$kkup_tgf)
      bind_rows(lapply(arms, function(a) {
        df <- sim_arm(mod_f, a, 72)
        df
      }))
    })
  }) %>% bindEvent(input$run_fib, ignoreNULL = FALSE)

  output$fib_plot  <- renderPlotly({
    df <- fib_sims()
    make_plotly(df, "FIB_SCORE", "Fibrosis Score",     "Stage (0-4)",   SCENARIO_COLORS)
  })
  output$hsc_plot  <- renderPlotly({
    df <- fib_sims()
    make_plotly(df, "HSC",       "HSC Activation",     "HSC (rel)",     SCENARIO_COLORS)
  })
  output$tgfb_plot <- renderPlotly({
    df <- fib_sims()
    make_plotly(df, "TGFB",      "TGF-β1",             "TGF-β1 (rel)",  SCENARIO_COLORS)
  })
  output$tnfa_plot <- renderPlotly({
    df <- fib_sims()
    make_plotly(df, "TNFA",      "TNF-α",              "TNF-α (rel)",   SCENARIO_COLORS)
  })
  output$il6_plot  <- renderPlotly({
    df <- fib_sims()
    make_plotly(df, "IL6C",      "IL-6",               "IL-6 (rel)",    SCENARIO_COLORS)
  })
  output$kup_plot  <- renderPlotly({
    df <- fib_sims()
    make_plotly(df, "KUPFFER",   "Kupffer Activation", "Kupffer (0-1)", SCENARIO_COLORS)
  })

  ## ── Metabolic biomarker plots ──────────────────────────
  meta_sims <- reactive({
    input$run_meta
    isolate({
      arms <- c("Placebo", "Resmetirom 100mg QD", "OCA 25mg QD",
                "Semaglutide 2.4mg QW", "Resmetirom + Semaglutide")
      mod_m <- param(mod_global,
                     IR0    = input$ir0,
                     WT0    = input$wt0,
                     ADIPON0= input$adipon0)
      bind_rows(lapply(arms, function(a) sim_arm(mod_m, a, 72)))
    })
  }) %>% bindEvent(input$run_meta, ignoreNULL = FALSE)

  output$homa_plot  <- renderPlotly({
    df <- meta_sims()
    make_plotly(df, "INS_RES",    "HOMA-IR",          "HOMA-IR",       SCENARIO_COLORS)
  })
  output$wt_plot    <- renderPlotly({
    df <- meta_sims()
    make_plotly(df, "BODY_WT",    "Body Weight",      "Weight (kg)",   SCENARIO_COLORS)
  })
  output$tg_plot    <- renderPlotly({
    df <- meta_sims()
    make_plotly(df, "TG_SERUM",   "Serum TG",         "TG (mg/dL)",    SCENARIO_COLORS, hline = 150)
  })
  output$ldl_plot   <- renderPlotly({
    df <- meta_sims()
    make_plotly(df, "LDL_C",      "LDL-C",            "LDL-C (mg/dL)", SCENARIO_COLORS, hline = 100)
  })
  output$adipon_plot<- renderPlotly({
    df <- meta_sims()
    make_plotly(df, "ADIPONECTIN","Adiponectin",       "µg/mL",         SCENARIO_COLORS)
  })
  output$fib4_plot  <- renderPlotly({
    df <- meta_sims()
    make_plotly(df, "FIB4",       "FIB-4 Index",      "FIB-4",         SCENARIO_COLORS, hline = 1.3)
  })

  ## ── Scenario comparison ────────────────────────────────
  output$compare_plot <- renderPlotly({
    df  <- all_sims()
    req(nrow(df) > 0)
    yv  <- input$metric_select
    make_plotly(df, yv,
                paste("Scenario Comparison:", yv),
                yv, SCENARIO_COLORS)
  })

  output$compare_table <- renderDT({
    df <- all_sims()
    req(nrow(df) > 0)
    tab <- df %>%
      filter(abs(Week - max(Week)) < 0.1) %>%
      select(Scenario, LF_PCT, FIB_SCORE, NAS, ALT_CMT,
             INS_RES, TG_SERUM, LDL_C, BODY_WT) %>%
      mutate(across(where(is.numeric), ~ round(.x, 1)))
    names(tab) <- c("Scenario", "Liver Fat (%)", "Fibrosis",
                    "NAS", "ALT (U/L)", "HOMA-IR",
                    "TG (mg/dL)", "LDL-C (mg/dL)", "Weight (kg)")
    tab
  }, rownames = FALSE,
     options = list(dom = "t", scrollX = TRUE, pageLength = 10))

  ## ── Biomarker tracker ─────────────────────────────────
  bm_sims <- reactive({
    input$run_bm
    isolate({
      arms <- c("Placebo", input$bm_drug)
      bind_rows(lapply(unique(arms), function(a) sim_arm(mod_global, a, 72)))
    })
  }) %>% bindEvent(input$run_bm, ignoreNULL = FALSE)

  output$bm_spider_plot <- renderPlotly({
    df <- bm_sims()
    bm <- intersect(input$bm_panel,
                    c("PDFF_pct","LSM_kPa","ALT_CMT","FIB4",
                      "INS_RES","TG_SERUM","LDL_C","ADIPONECTIN","TNFA","BODY_WT"))
    req(length(bm) >= 3)

    # Normalize to baseline (Week 0)
    base_df <- df %>% filter(Week < 0.1)
    end_df  <- df %>% filter(abs(Week - 72) < 0.5)

    plot_data <- lapply(unique(df$Scenario), function(sc) {
      b <- base_df %>% filter(Scenario == sc)
      e <- end_df  %>% filter(Scenario == sc)
      if (nrow(b) == 0 | nrow(e) == 0) return(NULL)
      vals <- sapply(bm, function(col) {
        bv <- mean(b[[col]], na.rm = TRUE)
        ev <- mean(e[[col]], na.rm = TRUE)
        if (bv == 0) return(1)
        ev / bv   # ratio to baseline
      })
      list(Scenario = sc, Biomarker = bm, Ratio = vals)
    })
    plot_data <- Filter(Negate(is.null), plot_data)

    p <- plot_ly(type = "scatterpolar", mode = "lines+markers")
    for (pd in plot_data) {
      col <- SCENARIO_COLORS[pd$Scenario]
      if (is.na(col)) col <- "#333333"
      p <- add_trace(p,
                     r    = c(pd$Ratio, pd$Ratio[1]),
                     theta= c(pd$Biomarker, pd$Biomarker[1]),
                     name = pd$Scenario,
                     line = list(color = col))
    }
    p %>% layout(
      polar  = list(radialaxis = list(visible = TRUE, range = c(0, 1.5))),
      title  = "Biomarker Spider Plot (Week 72 vs Baseline)",
      legend = list(orientation = "h", y = -0.15)
    )
  })

  output$bm_change_table <- renderDT({
    df <- bm_sims()
    bm_cols <- intersect(
      c("PDFF_pct","ALT_CMT","FIB4","INS_RES","TG_SERUM","LDL_C","ADIPONECTIN","BODY_WT"),
      names(df))

    timepoints <- c(12, 24, 48, 72)
    rows <- lapply(unique(df$Scenario), function(sc) {
      base_val <- df %>% filter(Scenario == sc, Week < 0.1) %>%
                  summarise(across(all_of(bm_cols), mean, na.rm = TRUE))
      lapply(timepoints, function(wk) {
        wk_val <- df %>% filter(Scenario == sc, abs(Week - wk) < 0.5) %>%
                  summarise(across(all_of(bm_cols), mean, na.rm = TRUE))
        if (nrow(wk_val) == 0) return(NULL)
        pct_change <- round((wk_val - base_val) / base_val * 100, 1)
        cbind(Scenario = sc, Week = wk, pct_change)
      })
    })
    bind_rows(unlist(rows, recursive = FALSE))
  }, rownames = FALSE, options = list(dom = "t", scrollX = TRUE, pageLength = 20))
}

## ── Launch ─────────────────────────────────────────────────
shinyApp(ui, server)
