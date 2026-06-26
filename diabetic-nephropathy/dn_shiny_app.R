## =============================================================================
## Diabetic Nephropathy QSP — Interactive Shiny Dashboard
## =============================================================================
## Tabs:
##  1. Patient Profile       — baseline demographics & CKD staging
##  2. Drug PK               — plasma concentration-time profiles
##  3. PD / Biomarkers       — TGF-β, ROS, ECM, Podocyte, AngII
##  4. Clinical Endpoints    — eGFR, UACR, SBP, CKD stage
##  5. Scenario Comparison   — side-by-side multi-treatment forest plots
##  6. GFR Slope Analysis    — annualized eGFR decline & ESKD risk
## =============================================================================

library(shiny)
library(shinydashboard)
library(mrgsolve)
library(dplyr)
library(ggplot2)
library(tidyr)
library(DT)
library(plotly)
library(scales)

## ---------------------------------------------------------------------------
## Embed the mrgsolve model
## ---------------------------------------------------------------------------
dn_code <- '
$PARAM
ka_acei=1.2  Vd_acei=25  CL_acei=5.5
ka_arb=0.9   Vd_arb=32   CL_arb=7.8
ka_sglt2=1.5 Vd_sglt2=74 CL_sglt2=9.2
ka_fine=1.1  Vd_fine=52  CL_fine=6.3
Emax_acei=0.90 EC50_acei=2.5
Emax_arb=0.85  EC50_arb=80
Emax_sglt2=0.85 EC50_sglt2=15
Emax_fine=0.88  EC50_fine=120
BG_base=8.5  BG_min=4.5  kBG_in=0.05  kBG_out=0.008
HbA1c_base=8.2
kAGE_in=0.01 kAGE_out=0.002
kAngII_in=0.15 kAngII_out=0.12
Pglo_base=48  Pglo_max=65  kPglo_AngII=8.0
kTGF_in=0.08  kTGF_out=0.06
TGF_AngII=0.35 TGF_AGE=0.20 TGF_ROS=0.25
kECM_in=0.04   kECM_out=0.01
ROS_base=1.0  kROS_in=0.12  kROS_out=0.10
ROS_BG=0.15   ROS_AngII=0.20  ROS_TGF=0.10
kPod_loss=0.003 Pod_TGF=0.8 Pod_AngII=0.6 Pod_ROS=0.5 Pod_Pglo=0.004 Pod_min=0.1
UACR_base=300 kUACR_Pod=200 kUACR_Pglo=3.0 kUACR_out=0.05
kTub_loss=0.002 Tub_UACR=0.001 Tub_hypoxia=0.3 kTub_out=0.005
kFib_in=0.012 kFib_out=0.002 Fib_TGF=0.6 Fib_Tub=0.4 Fib_base=0.1
GFR_0=55 GFR_min=5 kGFR_Fib=0.06 kGFR_ECM=0.04 kGFR_Pglo=0.003
SGLT2_GFR_dip=3.0
SBP_base=145 DBP_base=88 kSBP_AngII=10 kSBP_Na=5

$CMT
GI_acei CENT_acei GI_arb CENT_arb GI_sglt2 CENT_sglt2 GI_fine CENT_fine
BG AGE_cmpt AngII_cmpt TGF_cmpt ROS_cmpt ECM_cmpt Pod_cmpt UACR_cmpt Tub_cmpt Fib_cmpt GFR_cmpt

$INIT
GI_acei=0 CENT_acei=0 GI_arb=0 CENT_arb=0 GI_sglt2=0 CENT_sglt2=0 GI_fine=0 CENT_fine=0
BG=8.5 AGE_cmpt=1.0 AngII_cmpt=1.0 TGF_cmpt=1.0 ROS_cmpt=1.0 ECM_cmpt=1.0
Pod_cmpt=1.0 UACR_cmpt=300 Tub_cmpt=1.0 Fib_cmpt=0.1 GFR_cmpt=55

$ODE
dxdt_GI_acei   = -ka_acei * GI_acei;
dxdt_CENT_acei = ka_acei*GI_acei - (CL_acei/Vd_acei)*CENT_acei;
dxdt_GI_arb    = -ka_arb * GI_arb;
dxdt_CENT_arb  = ka_arb*GI_arb   - (CL_arb/Vd_arb)*CENT_arb;
dxdt_GI_sglt2  = -ka_sglt2*GI_sglt2;
dxdt_CENT_sglt2= ka_sglt2*GI_sglt2 - (CL_sglt2/Vd_sglt2)*CENT_sglt2;
dxdt_GI_fine   = -ka_fine*GI_fine;
dxdt_CENT_fine = ka_fine*GI_fine  - (CL_fine/Vd_fine)*CENT_fine;
double Cp_acei   = CENT_acei/Vd_acei*1000;
double E_acei    = Emax_acei*Cp_acei/(EC50_acei+Cp_acei);
double Cp_arb    = CENT_arb/Vd_arb*1000;
double E_arb     = Emax_arb*Cp_arb/(EC50_arb+Cp_arb);
double Cp_sglt2  = CENT_sglt2/Vd_sglt2*1000;
double E_sglt2   = Emax_sglt2*Cp_sglt2/(EC50_sglt2+Cp_sglt2);
double Cp_fine   = CENT_fine/Vd_fine*1000;
double E_fine    = Emax_fine*Cp_fine/(EC50_fine+Cp_fine);
double RAAS_block = 1-(1-E_acei)*(1-E_arb);
double BG_excess  = (BG>BG_min)?(BG-BG_min):0;
dxdt_BG = kBG_in*BG_base - kBG_out*BG - E_sglt2*0.015*BG;
dxdt_AGE_cmpt = kAGE_in*(BG_excess/(BG_base-BG_min)) - kAGE_out*AGE_cmpt;
dxdt_AngII_cmpt = kAngII_in*(1-RAAS_block) - kAngII_out*AngII_cmpt;
double TGF_dA = TGF_AngII*(AngII_cmpt-1)*(1-RAAS_block)*(1-E_fine);
double TGF_dG = TGF_AGE*(AGE_cmpt-1);
double TGF_dR = TGF_ROS*(ROS_cmpt-1);
double TGF_sglt2_red = E_sglt2*0.3*(TGF_cmpt-1);
dxdt_TGF_cmpt = kTGF_in + (TGF_dA>0?TGF_dA:0) + (TGF_dG>0?TGF_dG:0) + (TGF_dR>0?TGF_dR:0)
                - kTGF_out*TGF_cmpt - TGF_sglt2_red;
double ROS_in = kROS_in + ROS_BG*BG_excess/3.0 + ROS_AngII*(AngII_cmpt-1) + ROS_TGF*(TGF_cmpt-1);
dxdt_ROS_cmpt = ROS_in - kROS_out*ROS_cmpt;
dxdt_ECM_cmpt = kECM_in*TGF_cmpt*AngII_cmpt*(1-E_fine*0.5) - kECM_out*ECM_cmpt;
double Pod_eff = (Pod_cmpt>Pod_min)?Pod_cmpt:Pod_min;
double Pglo_curr = Pglo_base + kPglo_AngII*(AngII_cmpt-1)*(1-RAAS_block);
double Pglo_exc = (Pglo_curr>Pglo_base)?(Pglo_curr-Pglo_base):0;
double Pod_loss_rate = kPod_loss + Pod_TGF*(TGF_cmpt-1)*0.01
  + Pod_AngII*(AngII_cmpt-1)*0.005*(1-RAAS_block) + Pod_ROS*(ROS_cmpt-1)*0.004
  + Pod_Pglo*Pglo_exc*0.0005;
dxdt_Pod_cmpt = -Pod_loss_rate*Pod_eff;
double Pod_damage = 1.0-Pod_cmpt;
double UACR_ss = UACR_base*(1+kUACR_Pod*Pod_damage/UACR_base+kUACR_Pglo*Pglo_exc/UACR_base);
dxdt_UACR_cmpt = kUACR_out*(UACR_ss-UACR_cmpt);
double Tub_eff = (Tub_cmpt>0.05)?Tub_cmpt:0.05;
double Tub_loss_rate = kTub_loss + Tub_UACR*(UACR_cmpt-30)/1000 - Tub_hypoxia*E_sglt2*0.15
  + Tub_hypoxia*(BG_excess/5.0)*0.002;
dxdt_Tub_cmpt = -Tub_loss_rate*Tub_eff + kTub_out*(1.0-Tub_cmpt);
double Fib_in = kFib_in*(Fib_TGF*TGF_cmpt + Fib_Tub*(1.0-Tub_cmpt))*(1-E_fine*0.55);
dxdt_Fib_cmpt = Fib_in - kFib_out*Fib_cmpt;
double GFR_decline = kGFR_Fib*(Fib_cmpt-Fib_base) + kGFR_ECM*(ECM_cmpt-1.0) + kGFR_Pglo*Pglo_exc;
double GFR_floor = (GFR_cmpt>GFR_min)?1.0:0.0;
dxdt_GFR_cmpt = -(GFR_decline + E_sglt2*SGLT2_GFR_dip*0.05)*GFR_floor;

$TABLE
double Cp_acei_ng  = CENT_acei/Vd_acei*1000;
double Cp_arb_ng   = CENT_arb/Vd_arb*1000;
double Cp_sglt2_ng = CENT_sglt2/Vd_sglt2*1000;
double Cp_fine_ng  = CENT_fine/Vd_fine*1000;
double E_acei_t    = Emax_acei*Cp_acei_ng/(EC50_acei+Cp_acei_ng);
double E_arb_t     = Emax_arb*Cp_arb_ng/(EC50_arb+Cp_arb_ng);
double E_sglt2_t   = Emax_sglt2*Cp_sglt2_ng/(EC50_sglt2+Cp_sglt2_ng);
double E_fine_t    = Emax_fine*Cp_fine_ng/(EC50_fine+Cp_fine_ng);
double RAAS_block_t = 1-(1-E_acei_t)*(1-E_arb_t);
double SBP = SBP_base - 10*(E_acei_t + E_arb_t*0.8) - 5*E_sglt2_t;
double DBP = DBP_base - 5*(E_acei_t + E_arb_t*0.8);
double HbA1c = HbA1c_base - 0.5*E_sglt2_t;
double CKD_Stage;
if(GFR_cmpt>=90) CKD_Stage=1;
else if(GFR_cmpt>=60) CKD_Stage=2;
else if(GFR_cmpt>=45) CKD_Stage=3.1;
else if(GFR_cmpt>=30) CKD_Stage=3.2;
else if(GFR_cmpt>=15) CKD_Stage=4;
else CKD_Stage=5;

$CAPTURE
Cp_acei_ng Cp_arb_ng Cp_sglt2_ng Cp_fine_ng
E_acei_t E_arb_t E_sglt2_t E_fine_t RAAS_block_t
GFR_cmpt UACR_cmpt SBP DBP HbA1c
TGF_cmpt ROS_cmpt ECM_cmpt Pod_cmpt Fib_cmpt AngII_cmpt AGE_cmpt BG
CKD_Stage
'
mod <- mcode("DN_Shiny", dn_code, quiet=TRUE)

## ---------------------------------------------------------------------------
## Utility: run simulation
## ---------------------------------------------------------------------------
run_sim <- function(scenario_list, patient, end_days = 730) {
  # Initialize compartments from patient params
  # patient: list(GFR_0, UACR_0, BG_0, SBP_0, age, sex, dm_dur)
  init_vals <- list(
    BG = patient$BG_0,
    GFR_cmpt = patient$GFR_0,
    UACR_cmpt = patient$UACR_0
  )

  build_ev <- function(sc) {
    has_acei  <- sc$use_acei  && sc$dose_acei  > 0
    has_arb   <- sc$use_arb   && sc$dose_arb   > 0
    has_sglt2 <- sc$use_sglt2 && sc$dose_sglt2 > 0
    has_fine  <- sc$use_fine  && sc$dose_fine  > 0

    evl <- list()
    if (has_acei)  evl[[length(evl)+1]] <- ev(amt=sc$dose_acei,  cmt=1, ii=12, addl=end_days*2-1, time=0)
    if (has_arb)   evl[[length(evl)+1]] <- ev(amt=sc$dose_arb,   cmt=3, ii=24, addl=end_days-1,   time=0)
    if (has_sglt2) evl[[length(evl)+1]] <- ev(amt=sc$dose_sglt2, cmt=5, ii=24, addl=end_days-1,   time=0)
    if (has_fine)  evl[[length(evl)+1]] <- ev(amt=sc$dose_fine,  cmt=7, ii=24, addl=end_days-1,   time=0)
    if (length(evl) == 0) return(ev(time=0,amt=0,cmt=1))
    do.call(ev_c, evl)
  }

  results <- list()
  for (nm in names(scenario_list)) {
    sc   <- scenario_list[[nm]]
    evts <- build_ev(sc)
    out  <- mrgsim(mod, init=init_vals, ev=evts, end=end_days, delta=7, obsonly=TRUE)
    df   <- as.data.frame(out)
    df$Scenario  <- nm
    df$time_yr   <- df$time / 365
    results[[nm]] <- df
  }
  bind_rows(results)
}

## ---------------------------------------------------------------------------
## UI
## ---------------------------------------------------------------------------
ui <- dashboardPage(
  skin = "blue",

  dashboardHeader(
    title = tags$span(
      tags$img(src = "", height = 30),
      "Diabetic Nephropathy QSP"
    ),
    titleWidth = 300
  ),

  dashboardSidebar(
    width = 260,
    sidebarMenu(
      menuItem("Patient Profile",    tabName = "tab_patient",   icon = icon("user")),
      menuItem("Drug PK",            tabName = "tab_pk",        icon = icon("pills")),
      menuItem("PD / Biomarkers",    tabName = "tab_pd",        icon = icon("dna")),
      menuItem("Clinical Endpoints", tabName = "tab_endpoints", icon = icon("chart-line")),
      menuItem("Scenario Comparison",tabName = "tab_compare",   icon = icon("balance-scale")),
      menuItem("GFR Slope & ESKD",   tabName = "tab_slope",     icon = icon("tachometer-alt"))
    ),
    hr(),
    h5("Simulation Settings", style="padding-left:15px;color:#aaa;"),
    sliderInput("sim_years", "Horizon (years)", 1, 5, 3, step=0.5),
    actionButton("run_btn", "Run Simulation", class="btn-primary btn-block",
                 style="margin:10px 15px; width:230px;"),
    hr(),
    h5("Scenarios", style="padding-left:15px;color:#aaa;"),
    checkboxInput("sc_S0", "S0: No Treatment",    value=TRUE),
    checkboxInput("sc_S1", "S1: ACEi",            value=TRUE),
    checkboxInput("sc_S2", "S2: ARB",             value=FALSE),
    checkboxInput("sc_S3", "S3: SGLT2i",          value=TRUE),
    checkboxInput("sc_S4", "S4: ACEi+SGLT2i",     value=TRUE),
    checkboxInput("sc_S5", "S5: SGLT2i+Fine",     value=TRUE),
    checkboxInput("sc_S6", "S6: ACEi+SGLT2i+Fine",value=TRUE)
  ),

  dashboardBody(
    tags$head(tags$style(HTML("
      .content-wrapper { background: #f4f6f9; }
      .box { border-radius:8px; box-shadow: 0 2px 8px rgba(0,0,0,0.1); }
      .value-box { border-radius:8px; }
      h3 { font-weight: 600; }
    "))),

    tabItems(

      ## ================================================================
      ## TAB 1: Patient Profile
      ## ================================================================
      tabItem("tab_patient",
        fluidRow(
          box(title="Patient Baseline Parameters", width=4, status="primary",
            numericInput("pt_age",   "Age (years)",       value=60, min=30, max=80),
            selectInput( "pt_sex",   "Sex",               choices=c("Male","Female")),
            numericInput("pt_dm_dur","DM Duration (yrs)", value=12, min=1, max=40),
            numericInput("pt_GFR",   "Baseline eGFR",     value=55, min=5, max=120),
            numericInput("pt_UACR",  "Baseline UACR (mg/g)", value=300, min=3, max=3500),
            numericInput("pt_BG",    "Fasting Glucose (mmol/L)", value=8.5, min=4, max=20),
            numericInput("pt_SBP",   "Baseline SBP (mmHg)",     value=145, min=90, max=200),
            numericInput("pt_HbA1c", "HbA1c (%)",              value=8.2, min=5, max=15)
          ),
          box(title="Drug Doses", width=4, status="warning",
            h4("ACE Inhibitor (Enalapril)"),
            numericInput("dose_acei",  "ACEi dose (mg BID)",   value=10, min=0, max=40),
            h4("ARB (Losartan)"),
            numericInput("dose_arb",   "ARB dose (mg QD)",     value=100, min=0, max=200),
            h4("SGLT2i (Empagliflozin)"),
            numericInput("dose_sglt2", "SGLT2i dose (mg QD)",  value=25, min=0, max=50),
            h4("Finerenone"),
            numericInput("dose_fine",  "Finerenone dose (mg QD)", value=20, min=0, max=40)
          ),
          box(title="CKD Staging (Baseline)", width=4, status="info",
            valueBoxOutput("vb_gfr",  width=12),
            valueBoxOutput("vb_uacr", width=12),
            valueBoxOutput("vb_ckd",  width=12),
            plotOutput("plot_staging", height="250px")
          )
        ),
        fluidRow(
          box(title="Pathophysiology Overview", width=12, status="success",
            HTML("
              <div style='padding:10px;'>
                <h4>Diabetic Nephropathy — Key Mechanisms</h4>
                <table class='table table-bordered table-sm' style='font-size:13px;'>
                  <thead><tr><th>Pathway</th><th>Key Mediators</th><th>Clinical Impact</th></tr></thead>
                  <tbody>
                    <tr><td><b>Hemodynamic</b></td><td>AngII, PGLO↑, RAAS activation</td><td>Glomerular hypertension → GBM damage</td></tr>
                    <tr><td><b>Metabolic</b></td><td>AGE, PKC, Hexosamine, Polyol</td><td>Oxidative stress, endothelial dysfunction</td></tr>
                    <tr><td><b>TGF-β/Fibrosis</b></td><td>TGF-β1, Smad2/3, CTGF, ECM</td><td>Glomerulosclerosis, interstitial fibrosis</td></tr>
                    <tr><td><b>Inflammation</b></td><td>NF-κB, TNF-α, IL-1β, MCP-1</td><td>Macrophage infiltration, cytokine storm</td></tr>
                    <tr><td><b>Podocyte</b></td><td>Nephrin, Podocin, foot process</td><td>Slit diaphragm disruption → albuminuria</td></tr>
                  </tbody>
                </table>
              </div>
            ")
          )
        )
      ),

      ## ================================================================
      ## TAB 2: Drug PK
      ## ================================================================
      tabItem("tab_pk",
        fluidRow(
          box(title="PK Parameters", width=3, status="primary",
            selectInput("pk_drug", "Select Drug",
                        choices=c("ACEi (Enalapril)"="acei","ARB (Losartan)"="arb",
                                  "SGLT2i (Empa)"="sglt2","Finerenone"="fine")),
            numericInput("pk_dose", "Dose (mg)", value=25, min=1, max=200),
            numericInput("pk_days", "Duration (days)", value=3, min=1, max=7),
            numericInput("pk_freq", "Dosing interval (h)", value=24, min=8, max=48),
            actionButton("run_pk", "Plot PK", class="btn-info btn-block")
          ),
          box(title="Plasma Concentration-Time Profile", width=9, status="info",
            plotlyOutput("plot_pk_ct", height="380px")
          )
        ),
        fluidRow(
          box(title="PK Summary (Steady State ~Day 5)", width=6, status="success",
            DTOutput("table_pk_ss")
          ),
          box(title="Drug Effect (% Inhibition)", width=6, status="warning",
            plotlyOutput("plot_pk_effect", height="280px")
          )
        )
      ),

      ## ================================================================
      ## TAB 3: PD / Biomarkers
      ## ================================================================
      tabItem("tab_pd",
        fluidRow(
          valueBoxOutput("vb_tgf",  width=3),
          valueBoxOutput("vb_ros",  width=3),
          valueBoxOutput("vb_pod",  width=3),
          valueBoxOutput("vb_fib",  width=3)
        ),
        fluidRow(
          box(title="TGF-β Trajectory", width=6, status="danger",
            plotlyOutput("plot_tgf", height="280px")
          ),
          box(title="Oxidative Stress (ROS)", width=6, status="warning",
            plotlyOutput("plot_ros", height="280px")
          )
        ),
        fluidRow(
          box(title="ECM / Glomerulosclerosis", width=6, status="primary",
            plotlyOutput("plot_ecm", height="280px")
          ),
          box(title="Podocyte Integrity", width=6, status="success",
            plotlyOutput("plot_pod", height="280px")
          )
        ),
        fluidRow(
          box(title="Interstitial Fibrosis", width=6, status="info",
            plotlyOutput("plot_fib", height="280px")
          ),
          box(title="AngII / RAAS Activity", width=6, status="warning",
            plotlyOutput("plot_ang2", height="280px")
          )
        )
      ),

      ## ================================================================
      ## TAB 4: Clinical Endpoints
      ## ================================================================
      tabItem("tab_endpoints",
        fluidRow(
          box(title="eGFR Trajectory", width=6, status="primary",
            plotlyOutput("plot_efgr", height="350px")
          ),
          box(title="UACR (Proteinuria)", width=6, status="warning",
            plotlyOutput("plot_uacr", height="350px")
          )
        ),
        fluidRow(
          box(title="Systolic Blood Pressure", width=6, status="danger",
            plotlyOutput("plot_sbp", height="280px")
          ),
          box(title="HbA1c (Glycemic Control)", width=6, status="success",
            plotlyOutput("plot_hba1c", height="280px")
          )
        ),
        fluidRow(
          box(title="CKD Stage Progression", width=12, status="info",
            plotlyOutput("plot_ckd_stage", height="300px")
          )
        )
      ),

      ## ================================================================
      ## TAB 5: Scenario Comparison
      ## ================================================================
      tabItem("tab_compare",
        fluidRow(
          box(title="Year-End Outcome Summary", width=12, status="primary",
            DTOutput("table_compare"),
            br(),
            downloadButton("dl_compare", "Download CSV")
          )
        ),
        fluidRow(
          box(title="eGFR at Endpoint — Forest Plot", width=6, status="success",
            plotlyOutput("plot_forest_gfr", height="350px")
          ),
          box(title="UACR % Change — Forest Plot", width=6, status="warning",
            plotlyOutput("plot_forest_uacr", height="350px")
          )
        )
      ),

      ## ================================================================
      ## TAB 6: GFR Slope & ESKD
      ## ================================================================
      tabItem("tab_slope",
        fluidRow(
          box(title="Annualized eGFR Slope (mL/min/yr)", width=6, status="primary",
            plotlyOutput("plot_slope", height="350px")
          ),
          box(title="ESKD Risk Timeline", width=6, status="danger",
            plotlyOutput("plot_eskd", height="350px")
          )
        ),
        fluidRow(
          box(title="GFR Slope Details by Scenario", width=12, status="info",
            DTOutput("table_slope")
          )
        )
      )
    )
  )
)

## ---------------------------------------------------------------------------
## SERVER
## ---------------------------------------------------------------------------
server <- function(input, output, session) {

  ## Reactive: simulation results
  sim_data <- eventReactive(input$run_btn, {
    patient <- list(
      GFR_0  = input$pt_GFR,
      UACR_0 = input$pt_UACR,
      BG_0   = input$pt_BG,
      SBP_0  = input$pt_SBP
    )
    scenario_list <- list()
    if (input$sc_S0) scenario_list[["S0: No Treatment"]]   <- list(use_acei=F,use_arb=F,use_sglt2=F,use_fine=F,dose_acei=0,dose_arb=0,dose_sglt2=0,dose_fine=0)
    if (input$sc_S1) scenario_list[["S1: ACEi"]]           <- list(use_acei=T,use_arb=F,use_sglt2=F,use_fine=F,dose_acei=input$dose_acei,dose_arb=0,dose_sglt2=0,dose_fine=0)
    if (input$sc_S2) scenario_list[["S2: ARB"]]            <- list(use_acei=F,use_arb=T,use_sglt2=F,use_fine=F,dose_acei=0,dose_arb=input$dose_arb,dose_sglt2=0,dose_fine=0)
    if (input$sc_S3) scenario_list[["S3: SGLT2i"]]         <- list(use_acei=F,use_arb=F,use_sglt2=T,use_fine=F,dose_acei=0,dose_arb=0,dose_sglt2=input$dose_sglt2,dose_fine=0)
    if (input$sc_S4) scenario_list[["S4: ACEi+SGLT2i"]]   <- list(use_acei=T,use_arb=F,use_sglt2=T,use_fine=F,dose_acei=input$dose_acei,dose_arb=0,dose_sglt2=input$dose_sglt2,dose_fine=0)
    if (input$sc_S5) scenario_list[["S5: SGLT2i+Fine"]]   <- list(use_acei=F,use_arb=F,use_sglt2=T,use_fine=T,dose_acei=0,dose_arb=0,dose_sglt2=input$dose_sglt2,dose_fine=input$dose_fine)
    if (input$sc_S6) scenario_list[["S6: ACEi+SGLT2i+Fine"]] <- list(use_acei=T,use_arb=F,use_sglt2=T,use_fine=T,dose_acei=input$dose_acei,dose_arb=0,dose_sglt2=input$dose_sglt2,dose_fine=input$dose_fine)
    if (length(scenario_list) == 0) return(NULL)
    run_sim(scenario_list, patient, end_days = round(input$sim_years * 365))
  }, ignoreNULL=FALSE)

  ## Initialize with default run
  observe({ if (is.null(sim_data())) shinyjs::click("run_btn") })

  ## Tab 1 ValueBoxes
  output$vb_gfr <- renderValueBox({
    g <- input$pt_GFR
    stage <- if(g>=90)"G1" else if(g>=60)"G2" else if(g>=45)"G3a" else if(g>=30)"G3b" else if(g>=15)"G4" else "G5"
    col   <- if(g>=60)"green" else if(g>=30)"yellow" else "red"
    valueBox(paste0(g," mL/min"), paste("eGFR — CKD",stage), icon=icon("kidney"), color=col)
  })
  output$vb_uacr <- renderValueBox({
    u <- input$pt_UACR
    cat2 <- if(u<30)"A1 (normal)" else if(u<300)"A2 (micro)" else "A3 (macro)"
    col  <- if(u<30)"green" else if(u<300)"yellow" else "red"
    valueBox(paste0(u," mg/g"), paste("UACR —",cat2), icon=icon("tint"), color=col)
  })
  output$vb_ckd <- renderValueBox({
    g <- input$pt_GFR; u <- input$pt_UACR
    risk <- if(g>=60&&u<30)"Low" else if(g>=45||u<300)"Moderate" else if(g>=30||u<3500)"High" else "Very High"
    col  <- if(risk=="Low")"green" else if(risk=="Moderate")"yellow" else "red"
    valueBox(risk, "KDIGO Cardiovascular Risk", icon=icon("heartbeat"), color=col)
  })
  output$plot_staging <- renderPlot({
    df <- data.frame(
      Category = c("eGFR","UACR","HbA1c","SBP"),
      Value    = c(input$pt_GFR/120, input$pt_UACR/3500, input$pt_HbA1c/15, input$pt_SBP/200),
      Label    = c(paste0(input$pt_GFR," mL/min"),
                   paste0(input$pt_UACR," mg/g"),
                   paste0(input$pt_HbA1c,"%"),
                   paste0(input$pt_SBP," mmHg"))
    )
    ggplot(df, aes(Category, Value, fill=Category)) +
      geom_col() +
      geom_text(aes(label=Label), vjust=-0.3, size=3.5) +
      scale_y_continuous(limits=c(0,1.1)) +
      scale_fill_manual(values=c("eGFR"="#1E88E5","UACR"="#E53935","HbA1c"="#FB8C00","SBP"="#8E24AA")) +
      labs(y="Normalized value (0-1)", title="Patient Risk Profile") +
      theme_minimal(base_size=11) + theme(legend.position="none")
  })

  ## Helper: color palette
  get_colors <- function(scs) {
    pal <- c("S0: No Treatment"="#E53935","S1: ACEi"="#1E88E5","S2: ARB"="#43A047",
             "S3: SGLT2i"="#FB8C00","S4: ACEi+SGLT2i"="#8E24AA",
             "S5: SGLT2i+Fine"="#00897B","S6: ACEi+SGLT2i+Fine"="#F4511E")
    pal[scs]
  }

  ## Helper: generic plotly line
  make_plotly <- function(df, yvar, ytitle, title, log_y=FALSE) {
    if (is.null(df)) return(plotly_empty())
    df2 <- df %>% select(time_yr, Scenario, val=all_of(yvar))
    cols <- get_colors(unique(df2$Scenario))
    p <- plot_ly(df2, x=~time_yr, y=~val, color=~Scenario, colors=cols, type="scatter", mode="lines")
    p <- p %>% layout(title=title, xaxis=list(title="Time (years)"), yaxis=list(title=ytitle, type=if(log_y)"log" else "linear"))
    p
  }

  output$plot_efgr  <- renderPlotly({ make_plotly(sim_data(),"GFR_cmpt","eGFR (mL/min/1.73m²)","eGFR Trajectory") })
  output$plot_uacr  <- renderPlotly({ make_plotly(sim_data(),"UACR_cmpt","UACR (mg/g)","UACR (log scale)", log_y=TRUE) })
  output$plot_tgf   <- renderPlotly({ make_plotly(sim_data(),"TGF_cmpt","TGF-β (AU)","TGF-β") })
  output$plot_ros   <- renderPlotly({ make_plotly(sim_data(),"ROS_cmpt","ROS (AU)","Oxidative Stress") })
  output$plot_ecm   <- renderPlotly({ make_plotly(sim_data(),"ECM_cmpt","ECM (AU)","ECM Accumulation") })
  output$plot_pod   <- renderPlotly({ make_plotly(sim_data(),"Pod_cmpt","Podocyte (AU)","Podocyte Integrity") })
  output$plot_fib   <- renderPlotly({ make_plotly(sim_data(),"Fib_cmpt","Fibrosis (AU)","Interstitial Fibrosis") })
  output$plot_ang2  <- renderPlotly({ make_plotly(sim_data(),"AngII_cmpt","AngII (AU)","AngII / RAAS Activity") })
  output$plot_sbp   <- renderPlotly({ make_plotly(sim_data(),"SBP","SBP (mmHg)","Systolic BP") })
  output$plot_hba1c <- renderPlotly({ make_plotly(sim_data(),"HbA1c","HbA1c (%)","HbA1c") })

  output$plot_ckd_stage <- renderPlotly({
    df <- sim_data(); if (is.null(df)) return(plotly_empty())
    cols <- get_colors(unique(df$Scenario))
    plot_ly(df, x=~time_yr, y=~CKD_Stage, color=~Scenario, colors=cols, type="scatter", mode="lines") %>%
      layout(title="CKD Stage Over Time", xaxis=list(title="Time (years)"),
             yaxis=list(title="CKD Stage", tickvals=c(1,2,3.1,3.2,4,5),
                        ticktext=c("G1","G2","G3a","G3b","G4","G5")))
  })

  ## PD value boxes (at last time point, first scenario)
  get_last <- function(var) {
    df <- sim_data()
    if (is.null(df)) return(NA)
    df %>% filter(Scenario == df$Scenario[1]) %>%
      slice_max(time_yr, n=1) %>% pull({{var}})
  }
  output$vb_tgf <- renderValueBox({
    v <- get_last("TGF_cmpt")
    valueBox(round(v,2),"TGF-β (AU)",icon=icon("arrow-trend-up"), color=if(!is.na(v)&&v>1.5)"red" else "yellow")
  })
  output$vb_ros <- renderValueBox({
    v <- get_last("ROS_cmpt")
    valueBox(round(v,2),"Oxidative Stress",icon=icon("bolt"), color=if(!is.na(v)&&v>1.3)"orange" else "green")
  })
  output$vb_pod <- renderValueBox({
    v <- get_last("Pod_cmpt")
    valueBox(round(v,3),"Podocyte Integrity",icon=icon("shield-halved"),color=if(!is.na(v)&&v<0.7)"red" else "green")
  })
  output$vb_fib <- renderValueBox({
    v <- get_last("Fib_cmpt")
    valueBox(round(v,3),"Interstitial Fibrosis",icon=icon("layers"),color=if(!is.na(v)&&v>0.4)"red" else "yellow")
  })

  ## PK Tab
  pk_data <- eventReactive(input$run_pk, {
    drug_map <- list(
      acei  = list(cmt=1, ka=1.2, vd=25, cl=5.5, col="Cp_acei_ng"),
      arb   = list(cmt=3, ka=0.9, vd=32, cl=7.8, col="Cp_arb_ng"),
      sglt2 = list(cmt=5, ka=1.5, vd=74, cl=9.2, col="Cp_sglt2_ng"),
      fine  = list(cmt=7, ka=1.1, vd=52, cl=6.3, col="Cp_fine_ng")
    )
    dm <- drug_map[[input$pk_drug]]
    evts <- ev(amt=input$pk_dose, cmt=dm$cmt, ii=input$pk_freq, addl=input$pk_days*(24/input$pk_freq)-1, time=0)
    out  <- mrgsim(mod, ev=evts, end=input$pk_days*24, delta=0.5, obsonly=TRUE)
    df   <- as.data.frame(out)
    df$Cp <- df[[dm$col]]
    df$time_h <- df$time
    df
  })

  output$plot_pk_ct <- renderPlotly({
    df <- pk_data(); if (is.null(df)) return(plotly_empty())
    plot_ly(df, x=~time_h, y=~Cp, type="scatter", mode="lines",
            line=list(color="#1E88E5", width=2.5)) %>%
      layout(title=paste("PK Profile —", input$pk_drug, input$pk_dose,"mg"),
             xaxis=list(title="Time (h)"), yaxis=list(title="Plasma Conc (ng/mL)"))
  })

  output$table_pk_ss <- renderDT({
    df <- pk_data(); if (is.null(df)) return(datatable(data.frame()))
    ss <- df %>% filter(time_h >= max(time_h)-input$pk_freq) %>%
      summarise(Cmax=round(max(Cp),2), Cmin=round(min(Cp),2),
                Cavg=round(mean(Cp),2), AUC=round(sum(Cp)*0.5,1))
    datatable(ss, options=list(dom="t"))
  })

  output$plot_pk_effect <- renderPlotly({
    df <- pk_data(); if(is.null(df)) return(plotly_empty())
    eff_col <- switch(input$pk_drug, acei="E_acei_t", arb="E_arb_t", sglt2="E_sglt2_t", fine="E_fine_t")
    df$Effect_pct <- df[[eff_col]] * 100
    plot_ly(df, x=~time_h, y=~Effect_pct, type="scatter", mode="lines",
            line=list(color="#FB8C00", width=2)) %>%
      layout(title="Drug Effect (%)", xaxis=list(title="Time (h)"),
             yaxis=list(title="% Inhibition", range=c(0,100)))
  })

  ## Scenario Comparison Table
  compare_tab <- reactive({
    df <- sim_data(); if(is.null(df)) return(NULL)
    end_t <- max(df$time_yr, na.rm=TRUE)
    df %>%
      filter(abs(time_yr - end_t) < 0.1) %>%
      group_by(Scenario) %>%
      summarise(
        `eGFR (mL/min)` = round(mean(GFR_cmpt),1),
        `ΔGFR`          = round(mean(GFR_cmpt)-input$pt_GFR,1),
        `UACR (mg/g)`   = round(mean(UACR_cmpt),0),
        `ΔUACR%`        = round((mean(UACR_cmpt)-input$pt_UACR)/input$pt_UACR*100,1),
        `SBP (mmHg)`    = round(mean(SBP),1),
        `TGF-β`         = round(mean(TGF_cmpt),2),
        `Fibrosis`      = round(mean(Fib_cmpt),3),
        `Podocyte`      = round(mean(Pod_cmpt),3),
        `CKD Stage`     = round(mean(CKD_Stage),1),
        .groups="drop"
      ) %>% arrange(desc(`eGFR (mL/min)`))
  })

  output$table_compare <- renderDT({
    ct <- compare_tab(); if(is.null(ct)) return(datatable(data.frame()))
    datatable(ct, rownames=FALSE, options=list(scrollX=TRUE, pageLength=10)) %>%
      formatStyle("ΔGFR", color=styleInterval(c(-5,0), c("red","orange","green"))) %>%
      formatStyle("ΔUACR%", color=styleInterval(c(-20,0), c("green","orange","red")))
  })

  output$dl_compare <- downloadHandler(
    filename = function() paste0("DN_QSP_comparison_", Sys.Date(), ".csv"),
    content  = function(f) write.csv(compare_tab(), f, row.names=FALSE)
  )

  ## Forest plots
  output$plot_forest_gfr <- renderPlotly({
    ct <- compare_tab(); if(is.null(ct)) return(plotly_empty())
    ct <- ct %>% mutate(Scenario=factor(Scenario, levels=rev(Scenario)))
    cols <- get_colors(as.character(ct$Scenario))
    plot_ly(ct, x=~`eGFR (mL/min)`, y=~Scenario, type="bar", orientation="h",
            marker=list(color=unname(cols[as.character(ct$Scenario)]))) %>%
      add_segments(x=input$pt_GFR, xend=input$pt_GFR, y=0.5, yend=nrow(ct)+0.5,
                   line=list(dash="dash",color="red"), showlegend=FALSE) %>%
      layout(title="eGFR at Endpoint", xaxis=list(title="eGFR (mL/min/1.73m²)"))
  })

  output$plot_forest_uacr <- renderPlotly({
    ct <- compare_tab(); if(is.null(ct)) return(plotly_empty())
    ct <- ct %>% mutate(Scenario=factor(Scenario, levels=rev(Scenario)))
    cols <- get_colors(as.character(ct$Scenario))
    plot_ly(ct, x=~`ΔUACR%`, y=~Scenario, type="bar", orientation="h",
            marker=list(color=unname(cols[as.character(ct$Scenario)]))) %>%
      layout(title="UACR % Change from Baseline", xaxis=list(title="ΔUACR (%)"))
  })

  ## GFR Slope Tab
  slope_tab <- reactive({
    df <- sim_data(); if(is.null(df)) return(NULL)
    df %>% group_by(Scenario) %>%
      summarise(
        `GFR Slope (mL/yr)` = round(coef(lm(GFR_cmpt ~ time_yr))[2], 2),
        `Final eGFR`        = round(last(GFR_cmpt), 1),
        `ESKD (GFR<15)%`   = round(mean(GFR_cmpt < 15)*100, 1),
        .groups="drop"
      )
  })

  output$plot_slope <- renderPlotly({
    st <- slope_tab(); if(is.null(st)) return(plotly_empty())
    cols <- get_colors(st$Scenario)
    plot_ly(st, x=~Scenario, y=~`GFR Slope (mL/yr)`, type="bar",
            marker=list(color=unname(cols))) %>%
      layout(title="Annualized eGFR Slope",
             yaxis=list(title="mL/min/year"), xaxis=list(title=""))
  })

  output$plot_eskd <- renderPlotly({
    df <- sim_data(); if(is.null(df)) return(plotly_empty())
    cols <- get_colors(unique(df$Scenario))
    df2 <- df %>% mutate(below15 = GFR_cmpt < 15)
    plot_ly(df2, x=~time_yr, y=~GFR_cmpt, color=~Scenario, colors=cols,
            type="scatter", mode="lines") %>%
      add_segments(x=0, xend=max(df$time_yr), y=15, yend=15,
                   line=list(dash="dash",color="red",width=1.5),
                   showlegend=FALSE) %>%
      layout(title="ESKD Threshold (eGFR < 15)",
             xaxis=list(title="Time (years)"), yaxis=list(title="eGFR"))
  })

  output$table_slope <- renderDT({
    st <- slope_tab(); if(is.null(st)) return(datatable(data.frame()))
    datatable(st, rownames=FALSE, options=list(dom="t"))
  })
}

## ---------------------------------------------------------------------------
shinyApp(ui, server)
