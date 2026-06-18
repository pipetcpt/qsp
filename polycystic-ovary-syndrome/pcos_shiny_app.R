## ============================================================
## PCOS QSP Shiny Dashboard
## Polycystic Ovary Syndrome — Interactive Simulator
## Six tabs: Patient Profile · Hormone PK · Metabolic Markers ·
##            Clinical Endpoints · Scenario Comparison · Biomarkers
## ============================================================

library(shiny)
library(shinydashboard)
library(mrgsolve)
library(dplyr)
library(ggplot2)
library(tidyr)
library(DT)
library(plotly)

## ---- Embed model code (same as pcos_mrgsolve_model.R) -----
pcos_code <- '
$PARAM @annotated
k_GnRH_base:1.20:GnRH drive baseline
k_GnRH_E2neg:0.006:E2 neg feedback (per pg/mL)
k_GnRH_P4neg:0.30:P4 neg feedback (per ng/mL)
k_LH_stim:2.80:GnRH → LH rate (1/day)
k_LH_deg:0.85:LH degradation (1/day)
k_FSH_stim:1.10:GnRH → FSH rate (1/day)
k_FSH_deg:0.38:FSH degradation (1/day)
k_FSH_inhB:0.012:Inhibin B feedback on FSH
k_T_LH:2.00:LH → T production
k_T_IR:0.45:IR potentiation of CYP17A1
k_T_deg:0.38:T clearance (1/day)
k_E2_FSH:3.80:FSH → E2 production
k_E2_arT:0.09:T → E2 aromatization
k_E2_deg:0.65:E2 clearance (1/day)
k_P4_CL_max:14.0:Max P4 from CL (ng/mL/day)
k_P4_deg:1.60:P4 clearance (1/day)
k_AMH_base:9.00:AMH production rate
k_AMH_deg:0.045:AMH clearance (1/day)
AMH_PCOS_mult:2.80:PCOS AMH multiplier
AFC_ss_PCOS:19.0:AFC steady-state in PCOS
k_AFC_reg:0.04:AFC regression rate
k_DF_FSH:0.30:FSH → dominant follicle drive
k_DF_AMH_inh:0.28:AMH inhibits DF selection
FSH_thresh_DF:4.80:FSH threshold for DF
k_DF_decay:0.10:DF spontaneous decay
k_Ins_base:20.0:Insulin baseline parameter
k_Ins_prod:1.30:Glucose-stimulated insulin secretion
k_Ins_deg:0.42:Insulin clearance (1/day)
k_IR_PCOS:0.72:Intrinsic insulin resistance
k_Gluc_HGP:1.05:Hepatic glucose production
k_Gluc_util:0.013:Peripheral glucose utilization
IGFBP1_base:16.0:IGFBP-1 baseline (ng/mL)
k_IGFBP1_ins:0.85:Insulin suppression of IGFBP-1
SHBG_base:22.0:SHBG production set-point
k_SHBG_ins:0.55:Insulin suppression of SHBG
k_SHBG_EE:0.40:EE potentiation of SHBG
k_SHBG_eq:0.06:SHBG equilibration rate
CRP_base:3.20:hsCRP basal set-point
k_CRP_T:0.09:Testosterone pro-inflammation
k_CRP_BMI:0.14:BMI pro-inflammation
k_CRP_eq:0.10:CRP equilibration rate
BMI_ss_PCOS:28.8:BMI steady-state PCOS
k_BMI_drift:0.001:BMI background increase
HirsScore_base:12.0:Baseline FG hirsutism score
k_Hirs_FT:0.048:Free T driver of FG score
k_Hirs_decay:0.012:FG score decay rate
MET_ka:1.20:Metformin ka (1/h)
MET_F:0.50:Metformin bioavailability
MET_CL:45.0:Metformin CL (L/h)
MET_V1:80.0:Metformin V1 (L)
MET_V2:155.0:Metformin V2 (L)
MET_Q:12.0:Metformin Q (L/h)
MET_Emax_IR:0.35:Metformin max IR reduction
MET_EC50_IR:1.50:Metformin EC50 IR (mg/L)
MET_Emax_CYP17:0.26:Metformin max CYP17A1 inhibition
MET_EC50_CYP17:2.00:Metformin EC50 CYP17A1 (mg/L)
MET_Emax_BMI:0.004:Metformin BMI reduction
LET_ka:0.70:Letrozole ka (1/h)
LET_F:0.99:Letrozole F
LET_CL:2.10:Letrozole CL (L/h)
LET_V1:145.0:Letrozole V1 (L)
LET_Emax_AI:0.97:Letrozole max aromatase inhibition
LET_EC50_AI:0.012:Letrozole EC50 AI (mg/L)
EE_ka:1.50:EE ka (1/h)
EE_F:0.45:EE bioavailability
EE_CL:38.0:EE CL (L/h)
EE_V1:250.0:EE V1 (L)
EE_Emax_SHBG:3.50:EE max SHBG fold increase
EE_EC50_SHBG:0.06:EE EC50 SHBG (ng/mL)
EE_Emax_LH:0.92:EE max LH suppression
EE_EC50_LH:0.08:EE EC50 LH (ng/mL)
SPR_ka:1.80:Spiro ka (1/h)
SPR_F:0.90:Spiro F
SPR_CL:40.0:Spiro CL (L/h)
SPR_V1:65.0:Spiro V1 (L)
SPR_Ki_AR:0.18:Spiro AR Ki (mg/L)
SPR_Emax_T:0.32:Spiro max T production inhibition
CC_ka:0.80:Clomiphene ka (1/h)
CC_F:1.00:Clomiphene F
CC_CL:1.00:Clomiphene CL (L/h)
CC_V1:105.0:Clomiphene V1 (L)
CC_Emax_GnRH:0.70:Clomiphene GnRH increase
CC_EC50_GnRH:0.10:Clomiphene EC50 GnRH (mg/L)

$CMT GnRH_drive LH FSH T_total E2 P4 AMH AFC_state DF_state
     Insulin Glucose SHBG FreeT IGFBP1 CRP BMI HirsScore
     MET_gut MET_central MET_periph LET_plasma EE_plasma SPR_plasma CC_plasma

$INIT GnRH_drive=1.25 LH=12.5 FSH=5.0 T_total=68.0 E2=52.0 P4=0.4
      AMH=10.0 AFC_state=19.0 DF_state=0.04 Insulin=21.0 Glucose=101.0
      SHBG=21.0 FreeT=13.0 IGFBP1=15.0 CRP=3.6 BMI=29.0 HirsScore=12.0
      MET_gut=0 MET_central=0 MET_periph=0 LET_plasma=0 EE_plasma=0
      SPR_plasma=0 CC_plasma=0

$ODE
double Cp_MET=MET_central/MET_V1;double Cp_LET=LET_plasma/LET_V1;
double Cp_EE=EE_plasma/EE_V1;double Cp_SPR=SPR_plasma/SPR_V1;
double Cp_CC=CC_plasma/CC_V1;
double E_MET_IR=MET_Emax_IR*Cp_MET/(MET_EC50_IR+Cp_MET+1e-9);
double E_MET_CYP=MET_Emax_CYP17*Cp_MET/(MET_EC50_CYP17+Cp_MET+1e-9);
double E_MET_BMI=MET_Emax_BMI*Cp_MET/(MET_EC50_IR+Cp_MET+1e-9);
double E_LET_AI=LET_Emax_AI*Cp_LET/(LET_EC50_AI+Cp_LET+1e-9);
double EE_SHBG_fold=1.0+EE_Emax_SHBG*Cp_EE/(EE_EC50_SHBG+Cp_EE+1e-9);
double EE_LH_frac=1.0-EE_Emax_LH*Cp_EE/(EE_EC50_LH+Cp_EE+1e-9);
double SPR_AR_occ=Cp_SPR/(SPR_Ki_AR+Cp_SPR+1e-9);
double SPR_T_inh=SPR_Emax_T*Cp_SPR/(SPR_Ki_AR+Cp_SPR+1e-9);
double CC_GnRH_frac=1.0+CC_Emax_GnRH*Cp_CC/(CC_EC50_GnRH+Cp_CC+1e-9);
double E2_fb=k_GnRH_E2neg*E2;double P4_fb=k_GnRH_P4neg*P4;
double GnRH_ss=(k_GnRH_base/(1.0+E2_fb+P4_fb))*CC_GnRH_frac*EE_LH_frac;
dxdt_GnRH_drive=0.6*(GnRH_ss-GnRH_drive);
dxdt_LH=k_LH_stim*GnRH_drive-k_LH_deg*LH;
double InhB_p=AMH*0.55+1.0;double FSH_inhib=k_FSH_inhB*InhB_p;
dxdt_FSH=k_FSH_stim*GnRH_drive-(k_FSH_deg+FSH_inhib)*FSH;
double CYP17_act=(1.0+k_T_IR*Insulin/k_Ins_base)*(1.0-E_MET_CYP);
double T_prod=k_T_LH*LH*CYP17_act*(1.0-SPR_T_inh);
dxdt_T_total=T_prod-k_T_deg*T_total;
double aromatase_eff=1.0-E_LET_AI;
double E2_prod=(k_E2_FSH*FSH+k_E2_arT*T_total)*aromatase_eff;
dxdt_E2=E2_prod-k_E2_deg*E2;
double CL_rate=k_P4_CL_max*DF_state*DF_state;
dxdt_P4=CL_rate-k_P4_deg*P4;
double AMH_prod=k_AMH_base*AMH_PCOS_mult;
dxdt_AMH=AMH_prod-k_AMH_deg*AMH*AMH;
dxdt_AFC_state=k_AFC_reg*(AFC_ss_PCOS-AFC_state);
double FSH_eff_DF=(FSH>FSH_thresh_DF)?(FSH-FSH_thresh_DF)*k_DF_FSH:0.0;
double AMH_DF_inh=k_DF_AMH_inh*AMH;
double DF_drive=FSH_eff_DF/(1.0+AMH_DF_inh);
dxdt_DF_state=DF_drive*(1.0-DF_state)-k_DF_decay*DF_state;
double eff_IR=k_IR_PCOS*(1.0-E_MET_IR);
double Ins_prod=k_Ins_prod*Glucose/100.0;
dxdt_Insulin=Ins_prod-k_Ins_deg*Insulin;
double GU=k_Gluc_util*Glucose*Insulin/(21.0*(1.0+eff_IR));
double HGP=k_Gluc_HGP*(1.0-0.85*E_MET_IR);
dxdt_Glucose=HGP-GU;
double SHBG_ins_fac=1.0/(1.0+k_SHBG_ins*Insulin/k_Ins_base);
double SHBG_ss=SHBG_base*SHBG_ins_fac*EE_SHBG_fold;
dxdt_SHBG=k_SHBG_eq*(SHBG_ss-SHBG);
double SHBG_rel=22.0/(SHBG+1e-3);
double FreeT_ss=T_total*SHBG_rel*(1.0-SPR_AR_occ*0.5)*0.19;
dxdt_FreeT=0.5*(FreeT_ss-FreeT);
double IGFBP1_ss=IGFBP1_base/(1.0+k_IGFBP1_ins*Insulin/k_Ins_base);
dxdt_IGFBP1=0.12*(IGFBP1_ss-IGFBP1);
double CRP_ss=CRP_base*(1.0+k_CRP_T*T_total/65.0+k_CRP_BMI*(BMI-25.0)/25.0);
dxdt_CRP=k_CRP_eq*(CRP_ss-CRP);
dxdt_BMI=k_BMI_drift-E_MET_BMI;
double Hirs_drive=k_Hirs_FT*FreeT*(1.0-SPR_AR_occ);
dxdt_HirsScore=Hirs_drive-k_Hirs_decay*HirsScore;
dxdt_MET_gut=-MET_ka*MET_gut;
dxdt_MET_central=MET_F*MET_ka*MET_gut-(MET_CL+MET_Q)/MET_V1*MET_central+MET_Q/MET_V2*MET_periph;
dxdt_MET_periph=MET_Q/MET_V1*MET_central-MET_Q/MET_V2*MET_periph;
dxdt_LET_plasma=-LET_CL/LET_V1*LET_plasma;
dxdt_EE_plasma=-EE_CL/EE_V1*EE_plasma;
dxdt_SPR_plasma=-SPR_CL/SPR_V1*SPR_plasma;
dxdt_CC_plasma=-CC_CL/CC_V1*CC_plasma;

$TABLE
capture HOMA_IR=(Insulin*Glucose)/405.0;
capture FAI=T_total*0.0347/(SHBG+1e-3)*100.0;
capture LH_FSH=LH/(FSH+1e-3);
capture OvulProb=DF_state;
capture Cp_MET_out=MET_central/MET_V1;
capture Cp_LET_out=LET_plasma/LET_V1;
capture Cp_EE_out=EE_plasma/EE_V1;
capture Cp_SPR_out=SPR_plasma/SPR_V1;

$CAPTURE HOMA_IR FAI LH_FSH OvulProb Cp_MET_out Cp_LET_out Cp_EE_out Cp_SPR_out
'

## Compile model (cached on load)
pcos_mod <- mcode("PCOS_QSP_shiny", pcos_code)

## Dosing helper
make_dose_ev <- function(cmt, amount, freq_h = 8, days = 180) {
  times_h <- seq(0, days * 24, by = freq_h)
  ev(time = times_h / 24, cmt = cmt, amt = amount, evid = 1)
}

## Scenario builder
build_events <- function(use_met, met_dose, use_let, let_dose,
                          use_ocp, use_spr, spr_dose, use_cc) {
  evs <- list()
  if (use_met) evs[["met"]] <- make_dose_ev("MET_gut", met_dose / 3)
  if (use_let) {
    let_days <- unlist(lapply(0:5, function(i) (3 + i*28):(7 + i*28)))
    evs[["let"]] <- ev(time = let_days, cmt = "LET_plasma",
                       amt = let_dose * 0.99, evid = 1)
  }
  if (use_ocp) {
    ocp_days <- unlist(lapply(0:7, function(i) (i*28):(i*28+20)))
    ocp_days <- ocp_days[ocp_days <= 180]
    evs[["ocp"]] <- ev(time = ocp_days, cmt = "EE_plasma",
                       amt = 0.030 * 0.45, evid = 1)
  }
  if (use_spr) evs[["spr"]] <- make_dose_ev("SPR_plasma", spr_dose * 0.90, freq_h = 24)
  if (use_cc) {
    cc_days <- unlist(lapply(0:5, function(i) (3 + i*28):(7 + i*28)))
    evs[["cc"]] <- ev(time = cc_days, cmt = "CC_plasma", amt = 50 * 1.0, evid = 1)
  }
  if (length(evs) == 0) return(NULL)
  do.call(c, evs)
}

run_sim <- function(events, end_day = 180, mod = pcos_mod, ...) {
  extra_params <- list(...)
  if (length(extra_params) > 0) mod <- param(mod, .arg = extra_params)
  if (is.null(events))
    out <- mrgsim(mod, end = end_day, delta = 0.5) %>% as_tibble()
  else
    out <- mrgsim(mod, events = events, end = end_day, delta = 0.5) %>% as_tibble()
  out
}

## ---- UI ----
ui <- dashboardPage(
  skin = "purple",
  dashboardHeader(title = "PCOS QSP Simulator"),

  dashboardSidebar(
    sidebarMenu(
      menuItem("Patient Profile",      tabName = "profile",   icon = icon("user-md")),
      menuItem("Hormone Kinetics",     tabName = "hormones",  icon = icon("chart-line")),
      menuItem("Metabolic Endpoints",  tabName = "metabolic", icon = icon("heartbeat")),
      menuItem("Clinical Endpoints",   tabName = "clinical",  icon = icon("stethoscope")),
      menuItem("Scenario Comparison",  tabName = "scenarios", icon = icon("flask")),
      menuItem("Biomarker Dashboard",  tabName = "biomarkers",icon = icon("dna"))
    ),

    ## --- Drug selections ---
    hr(),
    h4("Treatment", style = "color:white; padding-left:10px"),
    checkboxInput("use_met",  "Metformin",     value = FALSE),
    conditionalPanel("input.use_met",
      sliderInput("met_dose", "Daily dose (mg)", 500, 2000, 1500, step = 250)
    ),
    checkboxInput("use_let",  "Letrozole",     value = FALSE),
    conditionalPanel("input.use_let",
      sliderInput("let_dose", "Dose (mg/day)", 2.5, 7.5, 2.5, step = 2.5)
    ),
    checkboxInput("use_ocp",  "Combined OCP",  value = FALSE),
    checkboxInput("use_spr",  "Spironolactone",value = FALSE),
    conditionalPanel("input.use_spr",
      sliderInput("spr_dose", "Daily dose (mg)", 50, 200, 100, step = 50)
    ),
    checkboxInput("use_cc",   "Clomiphene CC", value = FALSE),
    hr(),
    sliderInput("sim_days", "Simulation (days)", 30, 365, 180, step = 30)
  ),

  dashboardBody(
    tabItems(

      ## =========================================================
      ## TAB 1: Patient Profile
      ## =========================================================
      tabItem("profile",
        fluidRow(
          box(title = "Patient Parameters", width = 6, status = "primary",
            sliderInput("init_BMI",    "BMI (kg/m²)",        22, 42, 29, step = 0.5),
            sliderInput("init_IR",     "Insulin Resistance (k_IR)",  0, 1, 0.72, step = 0.05),
            sliderInput("init_AMH",    "AMH PCOS Multiplier",  1, 4, 2.8, step = 0.1),
            sliderInput("init_AFC",    "AFC (follicle count)", 5, 40, 19, step = 1)
          ),
          box(title = "Baseline Biomarkers", width = 6, status = "warning",
            tableOutput("baseline_table")
          )
        ),
        fluidRow(
          box(title = "PCOS Phenotype Classification (Rotterdam)", width = 12,
            status = "danger",
            HTML("
              <table class='table table-bordered'>
                <thead><tr><th>Phenotype</th><th>OA</th><th>HA</th><th>PCOM</th><th>Prevalence</th></tr></thead>
                <tbody>
                  <tr><td><b>A (Classic)</b></td><td>✓</td><td>✓</td><td>✓</td><td>~50%</td></tr>
                  <tr><td><b>B (Classic)</b></td><td>✓</td><td>✓</td><td>✗</td><td>~30%</td></tr>
                  <tr><td><b>C (Ovulatory)</b></td><td>✗</td><td>✓</td><td>✓</td><td>~12%</td></tr>
                  <tr><td><b>D (Normoandrogenic)</b></td><td>✓</td><td>✗</td><td>✓</td><td>~8%</td></tr>
                </tbody>
              </table>
              <p><i>OA = Oligo/Anovulation; HA = Hyperandrogenism; PCOM = Polycystic Ovary Morphology</i></p>
            ")
          )
        )
      ),

      ## =========================================================
      ## TAB 2: Hormone Kinetics
      ## =========================================================
      tabItem("hormones",
        fluidRow(
          box(title = "LH & FSH Dynamics", width = 6, status = "primary",
            plotlyOutput("plot_LH_FSH", height = 280)
          ),
          box(title = "Testosterone & Estradiol", width = 6, status = "warning",
            plotlyOutput("plot_T_E2", height = 280)
          )
        ),
        fluidRow(
          box(title = "Progesterone (P4)", width = 4, status = "success",
            plotlyOutput("plot_P4", height = 250)
          ),
          box(title = "AMH Level", width = 4, status = "danger",
            plotlyOutput("plot_AMH", height = 250)
          ),
          box(title = "SHBG & Free Testosterone", width = 4, status = "info",
            plotlyOutput("plot_SHBG", height = 250)
          )
        )
      ),

      ## =========================================================
      ## TAB 3: Metabolic Endpoints
      ## =========================================================
      tabItem("metabolic",
        fluidRow(
          box(title = "HOMA-IR (Insulin Resistance)", width = 6, status = "danger",
            plotlyOutput("plot_HOMA", height = 280)
          ),
          box(title = "Fasting Glucose & Insulin", width = 6, status = "warning",
            plotlyOutput("plot_GI", height = 280)
          )
        ),
        fluidRow(
          box(title = "BMI Trajectory", width = 4, status = "info",
            plotlyOutput("plot_BMI", height = 250)
          ),
          box(title = "IGFBP-1", width = 4, status = "primary",
            plotlyOutput("plot_IGFBP1", height = 250)
          ),
          box(title = "hsCRP (Inflammation)", width = 4, status = "warning",
            plotlyOutput("plot_CRP", height = 250)
          )
        )
      ),

      ## =========================================================
      ## TAB 4: Clinical Endpoints
      ## =========================================================
      tabItem("clinical",
        fluidRow(
          box(title = "Dominant Follicle / Ovulation Probability", width = 6,
            status = "success",
            plotlyOutput("plot_DF", height = 280)
          ),
          box(title = "Antral Follicle Count (AFC)", width = 6, status = "warning",
            plotlyOutput("plot_AFC", height = 280)
          )
        ),
        fluidRow(
          box(title = "Hirsutism (Ferriman-Gallwey Score)", width = 6, status = "danger",
            plotlyOutput("plot_hirs", height = 280)
          ),
          box(title = "Free Androgen Index (FAI)", width = 6, status = "warning",
            plotlyOutput("plot_FAI", height = 280)
          )
        )
      ),

      ## =========================================================
      ## TAB 5: Scenario Comparison
      ## =========================================================
      tabItem("scenarios",
        fluidRow(
          box(title = "Six-arm Treatment Comparison (Day 180)", width = 12,
            status = "primary",
            DT::dataTableOutput("scenario_table")
          )
        ),
        fluidRow(
          box(title = "Testosterone — All Scenarios", width = 6,
            plotlyOutput("plot_scen_T", height = 300)
          ),
          box(title = "HOMA-IR — All Scenarios", width = 6,
            plotlyOutput("plot_scen_HOMA", height = 300)
          )
        ),
        fluidRow(
          box(title = "Ovulation Probability — All Scenarios", width = 6,
            plotlyOutput("plot_scen_DF", height = 300)
          ),
          box(title = "FG Hirsutism Score — All Scenarios", width = 6,
            plotlyOutput("plot_scen_hirs", height = 300)
          )
        )
      ),

      ## =========================================================
      ## TAB 6: Biomarker Dashboard
      ## =========================================================
      tabItem("biomarkers",
        fluidRow(
          valueBoxOutput("vbox_LH",  width = 3),
          valueBoxOutput("vbox_T",   width = 3),
          valueBoxOutput("vbox_HOMA",width = 3),
          valueBoxOutput("vbox_AMH", width = 3)
        ),
        fluidRow(
          valueBoxOutput("vbox_FAI",  width = 3),
          valueBoxOutput("vbox_SHBG", width = 3),
          valueBoxOutput("vbox_CRP",  width = 3),
          valueBoxOutput("vbox_FG",   width = 3)
        ),
        fluidRow(
          box(title = "LH:FSH Ratio Over Time", width = 6,
            plotlyOutput("plot_LH_FSH_ratio", height = 260)
          ),
          box(title = "Risk Score Radar (Day 180)", width = 6,
            plotlyOutput("plot_radar", height = 260)
          )
        ),
        fluidRow(
          box(title = "Drug Plasma Concentrations", width = 12,
            status = "info",
            plotlyOutput("plot_drug_pk", height = 260)
          )
        )
      )
    )
  )
)

## ---- SERVER ----
server <- function(input, output, session) {

  ## Reactive: build events + run simulation
  sim_out <- reactive({
    events <- build_events(
      use_met = input$use_met, met_dose = input$met_dose,
      use_let = input$use_let, let_dose = input$let_dose,
      use_ocp = input$use_ocp,
      use_spr = input$use_spr, spr_dose = input$spr_dose,
      use_cc  = input$use_cc
    )
    run_sim(events, end_day = input$sim_days,
            AMH_PCOS_mult = input$init_AMH,
            k_IR_PCOS     = input$init_IR,
            BMI_ss_PCOS   = input$init_BMI,
            AFC_ss_PCOS   = input$init_AFC)
  })

  ## All-scenario comparison (fixed, run once)
  all_scenarios <- reactive({
    mod <- param(pcos_mod,
                 AMH_PCOS_mult = input$init_AMH,
                 k_IR_PCOS     = input$init_IR)

    run_one <- function(evs, label) {
      out <- if (is.null(evs)) mrgsim(mod, end = 180, delta = 1) else
                               mrgsim(mod, events = evs, end = 180, delta = 1)
      as_tibble(out) %>% mutate(Scenario = label)
    }

    ev_met <- make_dose_ev("MET_gut", 500)
    let_days <- unlist(lapply(0:5, function(i) (3+i*28):(7+i*28)))
    ev_let <- ev(time = let_days, cmt = "LET_plasma", amt = 2.5*0.99, evid = 1)
    ocp_days <- unlist(lapply(0:7, function(i) (i*28):(i*28+20))); ocp_days <- ocp_days[ocp_days<=180]
    ev_ocp <- ev(time = ocp_days, cmt = "EE_plasma", amt = 0.030*0.45, evid = 1)
    ev_spr <- make_dose_ev("SPR_plasma", 100*0.90, freq_h = 24)
    ev_met_let <- c(ev_met, ev_let)

    bind_rows(
      run_one(NULL,       "S1: Untreated"),
      run_one(ev_met,     "S2: Metformin"),
      run_one(ev_let,     "S3: Letrozole"),
      run_one(ev_ocp,     "S4: OCP"),
      run_one(ev_met_let, "S5: Met+Letrozole"),
      run_one(ev_spr,     "S6: Spironolactone")
    )
  })

  scen_colors <- c(
    "S1: Untreated"     = "#E74C3C",
    "S2: Metformin"     = "#27AE60",
    "S3: Letrozole"     = "#2980B9",
    "S4: OCP"           = "#8E44AD",
    "S5: Met+Letrozole" = "#E67E22",
    "S6: Spironolactone"= "#17A589"
  )

  ## Helper to make plotly line chart
  make_pl <- function(df, y_col, y_label, ref_line = NULL, ref_col = "darkgray") {
    p <- ggplot(df, aes(time, .data[[y_col]])) +
      geom_line(size = 0.9, color = "#2E4057") +
      labs(x = "Day", y = y_label) +
      theme_minimal(base_size = 11)
    if (!is.null(ref_line))
      p <- p + geom_hline(yintercept = ref_line, linetype = "dashed", color = ref_col)
    ggplotly(p)
  }

  ## --- TAB 1: Baseline table ---
  output$baseline_table <- renderTable({
    d <- sim_out() %>% filter(time < 0.1)
    data.frame(
      Marker  = c("LH (mIU/mL)","FSH (mIU/mL)","LH:FSH","T (ng/dL)","E2 (pg/mL)",
                  "SHBG (nmol/L)","FAI (%)","AMH (ng/mL)","Insulin (μIU/mL)",
                  "HOMA-IR","Glucose (mg/dL)","hsCRP (mg/L)","BMI","FG Score"),
      Value   = round(c(d$LH[1], d$FSH[1], d$LH_FSH[1], d$T_total[1], d$E2[1],
                        d$SHBG[1], d$FAI[1], d$AMH[1], d$Insulin[1],
                        d$HOMA_IR[1], d$Glucose[1], d$CRP[1], d$BMI[1], d$HirsScore[1]), 2),
      Normal  = c("<7","3-10","<2","<50","25-75",">30-120","<4.5","0.9-3.5","<15",
                  "<2.5","70-99","<1.0","<25","<8")
    )
  }, striped = TRUE, hover = TRUE, bordered = TRUE)

  ## --- TAB 2: Hormone Kinetics ---
  output$plot_LH_FSH <- renderPlotly({
    d <- sim_out() %>% select(time, LH, FSH) %>% pivot_longer(-time)
    p <- ggplot(d, aes(time, value, color = name)) +
      geom_line(size = 0.9) +
      scale_color_manual(values = c(LH = "#E74C3C", FSH = "#2980B9")) +
      geom_hline(yintercept = 10, linetype = "dashed", color = "#E74C3C", alpha = 0.5) +
      labs(x = "Day", y = "mIU/mL", color = "") + theme_minimal(base_size = 10)
    ggplotly(p)
  })

  output$plot_T_E2 <- renderPlotly({
    d <- sim_out() %>% select(time, T_total, E2) %>% pivot_longer(-time)
    p <- ggplot(d, aes(time, value, color = name)) +
      geom_line(size = 0.9) +
      scale_color_manual(values = c(T_total = "#E67E22", E2 = "#9B59B6")) +
      labs(x = "Day", y = "T: ng/dL | E2: pg/mL", color = "") +
      theme_minimal(base_size = 10)
    ggplotly(p)
  })

  output$plot_P4    <- renderPlotly(make_pl(sim_out(), "P4",      "P4 (ng/mL)",   1.5))
  output$plot_AMH   <- renderPlotly(make_pl(sim_out(), "AMH",     "AMH (ng/mL)",  3.5))
  output$plot_SHBG  <- renderPlotly({
    d <- sim_out() %>% select(time, SHBG, FreeT) %>% pivot_longer(-time)
    p <- ggplot(d, aes(time, value, color = name)) +
      geom_line(size=0.9) +
      scale_color_manual(values=c(SHBG="#1ABC9C", FreeT="#E74C3C")) +
      labs(x="Day", y="SHBG: nmol/L | FreeT: pg/mL", color="") +
      theme_minimal(base_size=10)
    ggplotly(p)
  })

  ## --- TAB 3: Metabolic ---
  output$plot_HOMA   <- renderPlotly(make_pl(sim_out(), "HOMA_IR", "HOMA-IR", 2.5))
  output$plot_GI     <- renderPlotly({
    d <- sim_out() %>% select(time, Insulin, Glucose) %>% pivot_longer(-time)
    p <- ggplot(d, aes(time, value, color=name)) + geom_line(size=0.9) +
      scale_color_manual(values=c(Insulin="#F39C12", Glucose="#16A085")) +
      labs(x="Day", y="Insulin: μIU/mL | Glucose: mg/dL", color="") +
      theme_minimal(base_size=10)
    ggplotly(p)
  })
  output$plot_BMI    <- renderPlotly(make_pl(sim_out(), "BMI",    "BMI (kg/m²)", 25))
  output$plot_IGFBP1 <- renderPlotly(make_pl(sim_out(), "IGFBP1","IGFBP-1 (ng/mL)"))
  output$plot_CRP    <- renderPlotly(make_pl(sim_out(), "CRP",   "hsCRP (mg/L)", 1.0))

  ## --- TAB 4: Clinical Endpoints ---
  output$plot_DF   <- renderPlotly(make_pl(sim_out(), "OvulProb", "DF / Ovulation Probability (0-1)"))
  output$plot_AFC  <- renderPlotly(make_pl(sim_out(), "AFC_state","Antral Follicle Count", 12))
  output$plot_hirs <- renderPlotly(make_pl(sim_out(), "HirsScore","Ferriman-Gallwey Score", 8))
  output$plot_FAI  <- renderPlotly(make_pl(sim_out(), "FAI",     "Free Androgen Index (%)", 4.5))

  ## --- TAB 5: Scenario Comparison ---
  output$scenario_table <- DT::renderDataTable({
    s <- all_scenarios() %>% filter(time == 180) %>%
      group_by(Scenario) %>%
      summarise(across(c(LH, FSH, LH_FSH, T_total, E2, P4, HOMA_IR,
                         SHBG, FAI, AMH, HirsScore, BMI, CRP, OvulProb),
                       ~ round(mean(.x), 2)), .groups = "drop") %>%
      rename(`LH (mIU/mL)` = LH, `FSH (mIU/mL)` = FSH, `LH:FSH` = LH_FSH,
             `T (ng/dL)` = T_total, `E2 (pg/mL)` = E2, `P4 (ng/mL)` = P4,
             `HOMA-IR` = HOMA_IR, `SHBG (nmol/L)` = SHBG, `FAI (%)` = FAI,
             `AMH (ng/mL)` = AMH, `FG Score` = HirsScore, `BMI` = BMI,
             `hsCRP` = CRP, `Ovul. Prob.` = OvulProb)
    DT::datatable(s, options = list(scrollX = TRUE, dom = "t"),
                  rownames = FALSE) %>%
      DT::formatStyle("HOMA-IR",
                      backgroundColor = DT::styleInterval(2.5, c("lightgreen", "#FFDDC1"))) %>%
      DT::formatStyle("T (ng/dL)",
                      backgroundColor = DT::styleInterval(50, c("lightgreen", "#FFDDC1")))
  })

  make_scen_plot <- function(y_col, y_label, ref = NULL) {
    d <- all_scenarios()
    p <- ggplot(d, aes(time, .data[[y_col]], color = Scenario)) +
      geom_line(size = 0.8, alpha = 0.9) +
      scale_color_manual(values = scen_colors) +
      labs(x = "Day", y = y_label, color = "") +
      theme_minimal(base_size = 10) + theme(legend.position = "bottom")
    if (!is.null(ref))
      p <- p + geom_hline(yintercept = ref, linetype = "dashed", color = "gray40")
    ggplotly(p)
  }

  output$plot_scen_T    <- renderPlotly(make_scen_plot("T_total",  "Testosterone (ng/dL)", 50))
  output$plot_scen_HOMA <- renderPlotly(make_scen_plot("HOMA_IR",  "HOMA-IR", 2.5))
  output$plot_scen_DF   <- renderPlotly(make_scen_plot("OvulProb", "Ovulation Probability"))
  output$plot_scen_hirs <- renderPlotly(make_scen_plot("HirsScore","FG Hirsutism Score", 8))

  ## --- TAB 6: Biomarker Dashboard ---
  last <- reactive(sim_out() %>% filter(time == max(time)) %>% slice(1))

  mk_vbox <- function(val, label, subtitle, icon_name, color) {
    valueBox(val, label, subtitle = subtitle, icon = icon(icon_name), color = color)
  }

  output$vbox_LH   <- renderValueBox(valueBox(round(last()$LH,1),    "LH (mIU/mL)",   icon=icon("arrow-up"),       color=if(last()$LH>10)"red" else "green"))
  output$vbox_T    <- renderValueBox(valueBox(round(last()$T_total,0),"T (ng/dL)",     icon=icon("venus-mars"),     color=if(last()$T_total>50)"red" else "green"))
  output$vbox_HOMA <- renderValueBox(valueBox(round(last()$HOMA_IR,2),"HOMA-IR",       icon=icon("tint"),           color=if(last()$HOMA_IR>2.5)"red" else "green"))
  output$vbox_AMH  <- renderValueBox(valueBox(round(last()$AMH,1),    "AMH (ng/mL)",   icon=icon("egg"),            color=if(last()$AMH>3.5)"yellow" else "green"))
  output$vbox_FAI  <- renderValueBox(valueBox(round(last()$FAI,2),    "FAI (%)",        icon=icon("thermometer"),   color=if(last()$FAI>4.5)"red" else "green"))
  output$vbox_SHBG <- renderValueBox(valueBox(round(last()$SHBG,1),   "SHBG (nmol/L)", icon=icon("shield-alt"),    color=if(last()$SHBG<30)"red" else "green"))
  output$vbox_CRP  <- renderValueBox(valueBox(round(last()$CRP,2),    "hsCRP (mg/L)",  icon=icon("fire"),          color=if(last()$CRP>3)"red" else "green"))
  output$vbox_FG   <- renderValueBox(valueBox(round(last()$HirsScore,1),"FG Score",    icon=icon("user"),          color=if(last()$HirsScore>=8)"red" else "green"))

  output$plot_LH_FSH_ratio <- renderPlotly(make_pl(sim_out(), "LH_FSH", "LH:FSH Ratio", 2.0))

  output$plot_radar <- renderPlotly({
    d <- last()
    ## Normalize to 0-1 scale (0=normal, 1=max abnormality)
    categories <- c("LH:FSH","T","HOMA-IR","FAI","hsCRP","FG Score")
    vals <- c(
      min((d$LH_FSH) / 4, 1),
      min((d$T_total) / 100, 1),
      min((d$HOMA_IR) / 5, 1),
      min((d$FAI) / 10, 1),
      min((d$CRP) / 8, 1),
      min((d$HirsScore) / 20, 1)
    )
    plot_ly(type = "scatterpolar", r = c(vals, vals[1]),
            theta = c(categories, categories[1]),
            fill = "toself", fillcolor = "rgba(231,76,60,0.3)",
            line = list(color = "#E74C3C")) %>%
      layout(polar = list(radialaxis = list(visible = TRUE, range = c(0,1))),
             showlegend = FALSE, title = "Risk Profile (normalized)")
  })

  output$plot_drug_pk <- renderPlotly({
    d <- sim_out() %>%
      select(time, Cp_MET_out, Cp_LET_out, Cp_EE_out, Cp_SPR_out) %>%
      pivot_longer(-time, names_to = "Drug", values_to = "Concentration")
    p <- ggplot(d, aes(time, Concentration, color = Drug)) +
      geom_line(size = 0.8, alpha = 0.9) +
      scale_color_manual(values = c(Cp_MET_out="#27AE60", Cp_LET_out="#E74C3C",
                                     Cp_EE_out="#9B59B6", Cp_SPR_out="#17A589")) +
      labs(x="Day", y="Concentration", color="Drug") +
      theme_minimal(base_size=10)
    ggplotly(p)
  })
}

## ---- LAUNCH ----
shinyApp(ui, server)
