## ============================================================
## CTEPH QSP Interactive Shiny Dashboard
## Chronic Thromboembolic Pulmonary Hypertension
## 7 Tabs: Patient Profile | PK | PD Signals | Hemodynamics |
##         Clinical Endpoints | Scenario Comparison | Biomarkers
## ============================================================

library(shiny)
library(shinydashboard)
library(mrgsolve)
library(dplyr)
library(tidyr)
library(ggplot2)
library(plotly)
library(DT)

## ---- mrgsolve model (embedded) -------------------------------
cteph_code <- '
$PARAM
CL_RIO=2.4, V1_RIO=30, V2_RIO=20, Q_RIO=1.8, ka_RIO=0.9, F_RIO=0.94, fm_RIO=0.15,
CL_M1=3.6, V_M1=25,
CL_MAC=1.1, V1_MAC=50, V2_MAC=40, Q_MAC=0.8, ka_MAC=0.35, F_MAC=0.75, fm_MAC=0.70,
CL_ACT=0.55, V_ACT=60,
CL_TREP=4.0, V_TREP=14, ka_TREP=2.0, F_TREP=0.79,
Emax_RIO=0.70, EC50_RIO=2.5, kin_cGMP=0.15, kout_cGMP=0.10, cGMP0=1.5,
Emax_ERA=0.55, EC50_ERA=0.8,
Emax_TREP=0.45, EC50_TREP=3.0, kin_cAMP=0.12, kout_cAMP=0.08, cAMP0=1.5,
TB0=0.75, kdeg_TB=0.002, kform_TB=0.003,
PVRf0=550, PVRv0=350, kin_PVRv=0.010, kout_PVRv=0.008,
ET1_0=3.5, kin_ET1=0.50, kout_ET1=0.14,
mPAP0=46, CO0=3.6, RVwork0=18, SaO2_0=90, BNP0=180, sixMWD0=342,
k6MWD_CO=45, k6MWD_SaO2=2.5,
PEA_done=0, BPA_done=0, BPA_sessions=0, AC_effect=1

$CMT
DEPOT_RIO C1_RIO C2_RIO MET_RIO
DEPOT_MAC C1_MAC C2_MAC MET_MAC
DEPOT_TREP C1_TREP
TB PVR_fixed PVR_var ET1 cGMP cAMP
RV_work mPAP CO SaO2 BNP sixMWD

$INIT
DEPOT_RIO=0,C1_RIO=0,C2_RIO=0,MET_RIO=0,
DEPOT_MAC=0,C1_MAC=0,C2_MAC=0,MET_MAC=0,
DEPOT_TREP=0,C1_TREP=0,
TB=0.75,PVR_fixed=550,PVR_var=350,ET1=3.5,cGMP=1.5,cAMP=1.5,
RV_work=18,mPAP=46,CO=3.6,SaO2=90,BNP=180,sixMWD=342

$ODE
double k10_RIO=CL_RIO/V1_RIO, k12_RIO=Q_RIO/V1_RIO, k21_RIO=Q_RIO/V2_RIO;
dxdt_DEPOT_RIO=-ka_RIO*DEPOT_RIO;
dxdt_C1_RIO=ka_RIO*F_RIO*DEPOT_RIO-(k10_RIO+k12_RIO)*C1_RIO+k21_RIO*C2_RIO;
dxdt_C2_RIO=k12_RIO*C1_RIO-k21_RIO*C2_RIO;
dxdt_MET_RIO=(fm_RIO*k10_RIO)*C1_RIO*V1_RIO/V_M1-(CL_M1/V_M1)*MET_RIO;
double Cp_RIO=C1_RIO/V1_RIO, Cm_RIO=MET_RIO/V_M1;

double k10_MAC=CL_MAC/V1_MAC, k12_MAC=Q_MAC/V1_MAC, k21_MAC=Q_MAC/V2_MAC;
dxdt_DEPOT_MAC=-ka_MAC*DEPOT_MAC;
dxdt_C1_MAC=ka_MAC*F_MAC*DEPOT_MAC-(k10_MAC+k12_MAC)*C1_MAC+k21_MAC*C2_MAC;
dxdt_C2_MAC=k12_MAC*C1_MAC-k21_MAC*C2_MAC;
dxdt_MET_MAC=(fm_MAC*k10_MAC)*C1_MAC*V1_MAC/V_ACT-(CL_ACT/V_ACT)*MET_MAC;
double Cp_MAC=C1_MAC/V1_MAC, Cm_MAC=MET_MAC/V_ACT;

dxdt_DEPOT_TREP=-ka_TREP*DEPOT_TREP;
dxdt_C1_TREP=ka_TREP*F_TREP*DEPOT_TREP-(CL_TREP/V_TREP)*C1_TREP;
double Cp_TREP=C1_TREP/V_TREP;

double E_RIO=(Emax_RIO*Cp_RIO)/(EC50_RIO+Cp_RIO);
double E_M1=(Emax_RIO*0.4*Cm_RIO)/(EC50_RIO*0.8+Cm_RIO);
double E_total_RIO=1.0+E_RIO+E_M1;
dxdt_cGMP=kin_cGMP*E_total_RIO-kout_cGMP*cGMP;

double E_TREP=(Emax_TREP*Cp_TREP)/(EC50_TREP+Cp_TREP);
double E_total_TREP=1.0+E_TREP;
dxdt_cAMP=kin_cAMP*E_total_TREP-kout_cAMP*cAMP;

double ERA_total=Cp_MAC+0.5*Cm_MAC;
double ET1_feedback=1.0+0.3*ERA_total/(0.5+ERA_total);
dxdt_ET1=kin_ET1*ET1_feedback-kout_ET1*ET1;

double AC_factor=AC_effect*0.6;
double kform_eff=kform_TB*(1.0-AC_factor);
double kdeg_eff=kdeg_TB*(1.0+0.5*E_TREP);
double TB_c=(TB<0)?0:(TB>1)?1:TB;
dxdt_TB=kform_eff*(1.0-TB_c)-kdeg_eff*TB_c;

double PEA_r=PEA_done*0.70;
double BPA_r=BPA_done*(1.0-pow(0.88,BPA_sessions));
double PVRf_t=PVR_fixed*(1.0-PEA_r)*(1.0-BPA_r);
if(PVRf_t<80) PVRf_t=80;
dxdt_PVR_fixed=-0.0015*(PVR_fixed-PVRf_t);

double cGMP_r=cGMP/cGMP0, cAMP_r=cAMP/cAMP0;
double ERA_eff=(Emax_ERA*ERA_total)/(EC50_ERA+ERA_total);
double ET1_r=ET1/ET1_0;
dxdt_PVR_var=kin_PVRv*ET1_r-kout_PVRv*(cGMP_r+cAMP_r+ERA_eff)*PVR_var;
if(PVR_var<100 && dxdt_PVR_var<0) dxdt_PVR_var=0;

double PVR_t=PVR_fixed+PVR_var;
double mPAP_t=(CO*PVR_t/80.0)+10.0;
dxdt_mPAP=0.02*(mPAP_t-mPAP);

double PVR_n=PVR_t/900.0; if(PVR_n>1) PVR_n=1;
double CO_t=5.2*(1.0-0.45*PVR_n);
double RV_adj=(RV_work<15)?(0.7+0.02*RV_work):1.0;
CO_t=CO_t*RV_adj; if(CO_t<1.5) CO_t=1.5;
dxdt_CO=0.015*(CO_t-CO);

double SV_e=(CO*1000.0)/80.0;
double mRAP_e=8.0+0.15*(mPAP-20.0);
double RVw_t=(mPAP-mRAP_e)*SV_e*0.0136;
dxdt_RV_work=0.01*(RVw_t-RV_work);

double SaO2_t=99.0-0.012*PVR_fixed-0.008*(CO0-CO)*10.0;
if(SaO2_t<70) SaO2_t=70; if(SaO2_t>99) SaO2_t=99;
dxdt_SaO2=0.05*(SaO2_t-SaO2);

double BNP_t=20.0*exp(0.08*(mPAP-20.0))*(5.0/(CO+0.01));
if(BNP_t>2000) BNP_t=2000;
dxdt_BNP=0.02*(BNP_t-BNP);

double s6_t=sixMWD0+k6MWD_CO*(CO-CO0)+k6MWD_SaO2*(SaO2-SaO2_0)-0.05*(mPAP-mPAP0);
if(s6_t<50) s6_t=50; if(s6_t>600) s6_t=600;
dxdt_sixMWD=0.008*(s6_t-sixMWD);

$TABLE
double Cp_RIO_out=C1_RIO/V1_RIO, Cm_RIO_out=MET_RIO/V_M1;
double Cp_MAC_out=C1_MAC/V1_MAC, Cm_MAC_out=MET_MAC/V_ACT;
double Cp_TREP_out=C1_TREP/V_TREP;
double PVR_total=PVR_fixed+PVR_var;
double WHO_FC=(sixMWD>440)?1.5:(sixMWD>300)?2.0:(sixMWD>150)?3.0:4.0;

$CAPTURE
Cp_RIO_out Cm_RIO_out Cp_MAC_out Cm_MAC_out Cp_TREP_out
ET1 cGMP cAMP TB PVR_fixed PVR_var PVR_total
mPAP CO SaO2 RV_work BNP sixMWD WHO_FC
'

mod_global <- mcode("cteph_shiny", cteph_code, quiet = TRUE)

## ---- Helper: run simulation ----------------------------------
run_sim <- function(
  use_rio, rio_dose, rio_freq,
  use_mac, mac_dose,
  use_trep, trep_dose,
  pvr_fixed_init, pvr_var_init, mpap_init, co_init,
  sao2_init, bnp_init, sixmwd_init, tb_init,
  pea_done, bpa_done, bpa_sessions, ac_effect,
  tend_wk = 52, delta_h = 24
) {
  tend  <- tend_wk * 7 * 24
  delta <- delta_h

  inits <- c(
    TB = tb_init, PVR_fixed = pvr_fixed_init, PVR_var = pvr_var_init,
    mPAP = mpap_init, CO = co_init, SaO2 = sao2_init,
    BNP = bnp_init, sixMWD = sixmwd_init,
    ET1 = 3.5, cGMP = 1.5, cAMP = 1.5, RV_work = 18,
    DEPOT_RIO=0, C1_RIO=0, C2_RIO=0, MET_RIO=0,
    DEPOT_MAC=0, C1_MAC=0, C2_MAC=0, MET_MAC=0,
    DEPOT_TREP=0, C1_TREP=0
  )

  params_ov <- list(
    PEA_done = as.numeric(pea_done), BPA_done = as.numeric(bpa_done),
    BPA_sessions = bpa_sessions, AC_effect = as.numeric(ac_effect),
    mPAP0 = mpap_init, CO0 = co_init, SaO2_0 = sao2_init,
    sixMWD0 = sixmwd_init
  )

  mod2 <- param(mod_global, params_ov)
  mod2 <- init(mod2, inits)

  events <- NULL
  if (use_rio && rio_dose > 0) {
    ii_h <- 24 / rio_freq
    ev_r  <- ev(amt = rio_dose, cmt = "DEPOT_RIO", ii = ii_h, addl = round(tend/ii_h))
    events <- if (is.null(events)) ev_r else events + ev_r
  }
  if (use_mac && mac_dose > 0) {
    ev_m  <- ev(amt = mac_dose, cmt = "DEPOT_MAC", ii = 24, addl = round(tend/24))
    events <- if (is.null(events)) ev_m else events + ev_m
  }
  if (use_trep && trep_dose > 0) {
    ev_t  <- ev(amt = trep_dose, cmt = "DEPOT_TREP", ii = 6, addl = round(tend/6))
    events <- if (is.null(events)) ev_t else events + ev_t
  }

  out <- if (is.null(events)) {
    mrgsim(mod2, end = tend, delta = delta)
  } else {
    mrgsim(mod2, events = events, end = tend, delta = delta)
  }
  df <- as.data.frame(out)
  df$time_weeks <- df$time / (7 * 24)
  df
}

## ============================================================
## SHINY UI
## ============================================================
ui <- dashboardPage(
  skin = "blue",
  dashboardHeader(title = "CTEPH QSP Dashboard"),

  dashboardSidebar(
    sidebarMenu(
      menuItem("Patient Profile",    tabName = "tab_patient",   icon = icon("user-md")),
      menuItem("Pharmacokinetics",   tabName = "tab_pk",        icon = icon("pills")),
      menuItem("PD Signals",         tabName = "tab_pd",        icon = icon("flask")),
      menuItem("Hemodynamics",       tabName = "tab_hemo",      icon = icon("heartbeat")),
      menuItem("Clinical Endpoints", tabName = "tab_endpoints", icon = icon("chart-line")),
      menuItem("Scenario Comparison",tabName = "tab_compare",   icon = icon("layer-group")),
      menuItem("Biomarkers",         tabName = "tab_biomarkers",icon = icon("microscope"))
    )
  ),

  dashboardBody(
    tags$head(tags$style(HTML("
      .content-wrapper, .right-side { background-color: #f4f6f9; }
      .box-header { background-color: #2C3E50 !important; color: white !important; }
      .nav-tabs-custom > .nav-tabs > li.active { border-top-color: #2980B9; }
      .skin-blue .main-header .logo { background-color: #2C3E50; }
      .skin-blue .main-header .navbar { background-color: #34495E; }
    "))),

    tabItems(

      ## ---- TAB 1: Patient Profile --------------------------------
      tabItem(tabName = "tab_patient",
        fluidRow(
          box(title = "Patient Characteristics", width = 6, status = "primary",
            h4("Baseline Hemodynamics"),
            sliderInput("pvr_fixed_init", "Fixed PVR (dyn·s/cm⁵):", 200, 1200, 550, 50),
            sliderInput("pvr_var_init",   "Variable PVR (dyn·s/cm⁵):", 50, 600, 350, 25),
            sliderInput("mpap_init",      "Baseline mPAP (mmHg):", 20, 80, 46, 1),
            sliderInput("co_init",        "Baseline CO (L/min):", 1.5, 6.0, 3.6, 0.1),
            sliderInput("sao2_init",      "Baseline SaO₂ (%):", 70, 99, 90, 1),
            sliderInput("tb_init",        "Thrombotic Burden (0-1):", 0.1, 1.0, 0.75, 0.05)
          ),
          box(title = "Biomarker Baseline", width = 6, status = "primary",
            sliderInput("bnp_init",    "Baseline BNP (pg/mL):", 20, 1000, 180, 10),
            sliderInput("sixmwd_init", "Baseline 6MWD (m):", 50, 600, 342, 10),
            hr(),
            h4("Procedural Intervention"),
            checkboxInput("pea_done",  "Post-PEA (Pulmonary Endarterectomy)", FALSE),
            checkboxInput("bpa_done",  "BPA Procedure Performed", FALSE),
            conditionalPanel("input.bpa_done == true",
              sliderInput("bpa_sessions", "Number of BPA Sessions:", 1, 10, 5, 1)
            ),
            checkboxInput("ac_effect", "On Anticoagulation", TRUE),
            hr(),
            h4("Simulation Duration"),
            sliderInput("tend_wk", "Simulation Duration (weeks):", 12, 104, 52, 4)
          )
        ),
        fluidRow(
          box(title = "Disease Overview", width = 12, status = "info",
            h4("Chronic Thromboembolic Pulmonary Hypertension (CTEPH)"),
            p("CTEPH is a distinct form of pulmonary hypertension caused by incomplete resolution
              of pulmonary emboli, leading to organized thrombi and secondary vascular remodeling."),
            p(strong("Diagnostic Criteria (ESC/ERS 2022):")),
            tags$ul(
              tags$li("mPAP > 20 mmHg at rest"),
              tags$li("PCWP ≤ 15 mmHg (pre-capillary)"),
              tags$li("Perfusion defects on V/Q scan after ≥3 months anticoagulation"),
              tags$li("CT-PA or conventional pulmonary angiography confirming chronic thrombi")
            ),
            tags$table(class = "table table-bordered",
              tags$thead(tags$tr(
                tags$th("Parameter"), tags$th("Normal"), tags$th("Mild CTEPH"),
                tags$th("Moderate CTEPH"), tags$th("Severe CTEPH")
              )),
              tags$tbody(
                tags$tr(tags$td("mPAP (mmHg)"), tags$td("<20"), tags$td("20-35"), tags$td("35-55"), tags$td(">55")),
                tags$tr(tags$td("PVR (dyn·s/cm⁵)"), tags$td("<240"), tags$td("240-400"), tags$td("400-800"), tags$td(">800")),
                tags$tr(tags$td("CO (L/min)"), tags$td(">4.5"), tags$td("3.5-4.5"), tags$td("2.5-3.5"), tags$td("<2.5")),
                tags$tr(tags$td("6MWD (m)"), tags$td(">500"), tags$td("400-500"), tags$td("250-400"), tags$td("<250"))
              )
            )
          )
        )
      ),

      ## ---- TAB 2: Pharmacokinetics -------------------------------
      tabItem(tabName = "tab_pk",
        fluidRow(
          box(title = "Drug Dosing", width = 4, status = "primary",
            h4("Riociguat (sGC Stimulator)"),
            checkboxInput("use_rio", "Enable Riociguat", TRUE),
            conditionalPanel("input.use_rio == true",
              sliderInput("rio_dose", "Dose (mg):", 0.5, 2.5, 2.5, 0.5),
              sliderInput("rio_freq", "Frequency (doses/day):", 1, 3, 3, 1)
            ),
            hr(),
            h4("Macitentan (ERA)"),
            checkboxInput("use_mac", "Enable Macitentan", FALSE),
            conditionalPanel("input.use_mac == true",
              sliderInput("mac_dose", "Dose (mg):", 3, 10, 10, 1)
            ),
            hr(),
            h4("Treprostinil (Prostacyclin)"),
            checkboxInput("use_trep", "Enable Treprostinil", FALSE),
            conditionalPanel("input.use_trep == true",
              sliderInput("trep_dose", "Dose (mg equiv.):", 0.1, 1.0, 0.5, 0.1)
            ),
            actionButton("run_sim", "Run Simulation", class = "btn-primary btn-lg",
                         style = "width:100%; margin-top:10px;")
          ),
          box(title = "PK Profiles Over Time", width = 8, status = "info",
            plotlyOutput("pk_plot", height = "450px")
          )
        ),
        fluidRow(
          box(title = "Riociguat PK Parameters", width = 4, status = "warning",
            tableOutput("pk_params_rio")
          ),
          box(title = "Macitentan PK Parameters", width = 4, status = "warning",
            tableOutput("pk_params_mac")
          ),
          box(title = "Treprostinil PK Parameters", width = 4, status = "warning",
            tableOutput("pk_params_trep")
          )
        )
      ),

      ## ---- TAB 3: PD Signals ------------------------------------
      tabItem(tabName = "tab_pd",
        fluidRow(
          box(title = "Second Messenger Dynamics", width = 6, status = "primary",
            plotlyOutput("cgmp_plot", height = "350px")
          ),
          box(title = "ET-1 & Vascular Tone", width = 6, status = "primary",
            plotlyOutput("et1_plot", height = "350px")
          )
        ),
        fluidRow(
          box(title = "PVR Components (Fixed vs Variable)", width = 6, status = "info",
            plotlyOutput("pvr_decomp_plot", height = "350px")
          ),
          box(title = "Thrombotic Burden Over Time", width = 6, status = "warning",
            plotlyOutput("tb_plot", height = "350px")
          )
        )
      ),

      ## ---- TAB 4: Hemodynamics ----------------------------------
      tabItem(tabName = "tab_hemo",
        fluidRow(
          box(title = "Mean Pulmonary Artery Pressure", width = 6, status = "danger",
            plotlyOutput("mpap_plot", height = "350px")
          ),
          box(title = "Total PVR", width = 6, status = "danger",
            plotlyOutput("pvr_total_plot", height = "350px")
          )
        ),
        fluidRow(
          box(title = "Cardiac Output", width = 6, status = "primary",
            plotlyOutput("co_plot", height = "350px")
          ),
          box(title = "Arterial O₂ Saturation (SaO₂)", width = 6, status = "success",
            plotlyOutput("sao2_plot", height = "350px")
          )
        ),
        fluidRow(
          box(title = "RV Stroke Work Index", width = 6, status = "warning",
            plotlyOutput("rvwork_plot", height = "350px")
          ),
          box(title = "Hemodynamic Summary at Selected Week", width = 6, status = "info",
            sliderInput("summary_week", "Select Week:", 1, 52, 16, 1),
            tableOutput("hemo_summary_table")
          )
        )
      ),

      ## ---- TAB 5: Clinical Endpoints ----------------------------
      tabItem(tabName = "tab_endpoints",
        fluidRow(
          box(title = "6-Minute Walk Distance (6MWD)", width = 6, status = "success",
            plotlyOutput("sixmwd_plot", height = "350px")
          ),
          box(title = "BNP (RV Stress Biomarker)", width = 6, status = "warning",
            plotlyOutput("bnp_plot", height = "350px")
          )
        ),
        fluidRow(
          box(title = "WHO Functional Class", width = 6, status = "primary",
            plotlyOutput("whofc_plot", height = "350px")
          ),
          box(title = "Clinical Response Waterfall", width = 6, status = "info",
            plotlyOutput("waterfall_plot", height = "350px")
          )
        ),
        fluidRow(
          box(title = "Key Efficacy Endpoints Summary (Week 16 & 52)", width = 12, status = "primary",
            DTOutput("endpoints_table")
          )
        )
      ),

      ## ---- TAB 6: Scenario Comparison ---------------------------
      tabItem(tabName = "tab_compare",
        fluidRow(
          box(title = "Predefined Clinical Scenarios", width = 3, status = "primary",
            checkboxGroupInput("scenarios_sel",
              "Select Scenarios:",
              choices = c(
                "No Treatment" = "s1",
                "Anticoagulation Only" = "s2",
                "Riociguat 2.5mg TID" = "s3",
                "Macitentan 10mg QD" = "s4",
                "Riociguat + Macitentan" = "s5",
                "BPA (5x) + Riociguat" = "s6",
                "Post-PEA + Combination" = "s7"
              ),
              selected = c("s1","s3","s5","s7")
            ),
            selectInput("compare_var",
              "Endpoint to Compare:",
              choices = c("mPAP","PVR_total","CO","SaO2","BNP","sixMWD","WHO_FC"),
              selected = "mPAP"
            ),
            actionButton("run_compare", "Run Comparison",
                         class = "btn-warning btn-lg", style = "width:100%;")
          ),
          box(title = "Scenario Comparison Plot", width = 9, status = "info",
            plotlyOutput("compare_plot", height = "500px")
          )
        ),
        fluidRow(
          box(title = "Scenario Comparison Table (Week 52)", width = 12, status = "primary",
            DTOutput("compare_table")
          )
        )
      ),

      ## ---- TAB 7: Biomarkers ------------------------------------
      tabItem(tabName = "tab_biomarkers",
        fluidRow(
          box(title = "Biomarker Trajectories", width = 8, status = "primary",
            plotlyOutput("biomarker_plot", height = "450px")
          ),
          box(title = "Risk Stratification", width = 4, status = "warning",
            h4("CTEPH Risk Markers"),
            p("ESC/ERS 2022 risk categories:"),
            tags$table(class = "table table-sm table-bordered",
              tags$thead(tags$tr(
                tags$th("Marker"), tags$th("Low"), tags$th("High")
              )),
              tags$tbody(
                tags$tr(tags$td("WHO FC"), tags$td("I-II"), tags$td("III-IV")),
                tags$tr(tags$td("6MWD (m)"), tags$td(">440"), tags$td("<165")),
                tags$tr(tags$td("BNP (pg/mL)"), tags$td("<50"), tags$td(">300")),
                tags$tr(tags$td("mPAP (mmHg)"), tags$td("<38"), tags$td(">54")),
                tags$tr(tags$td("PVR (WU)"), tags$td("<4"), tags$td(">9")),
                tags$tr(tags$td("CO (L/min)"), tags$td(">2.5"), tags$td("<2.0")),
                tags$tr(tags$td("TAPSE (mm)"), tags$td(">20"), tags$td("<12"))
              )
            ),
            hr(),
            h4("Current Patient Risk Profile"),
            uiOutput("risk_profile_ui")
          )
        ),
        fluidRow(
          box(title = "Drug Efficacy on Biomarkers at Week 16", width = 12, status = "info",
            plotlyOutput("biomarker_week16_plot", height = "350px")
          )
        )
      )
    ) # tabItems
  ) # dashboardBody
) # dashboardPage

## ============================================================
## SHINY SERVER
## ============================================================
server <- function(input, output, session) {

  ## ---- Reactive simulation ----------------------------------
  sim_data <- eventReactive(input$run_sim, {
    withProgress(message = "Running CTEPH simulation...", {
      run_sim(
        use_rio   = input$use_rio,
        rio_dose  = input$rio_dose,
        rio_freq  = input$rio_freq,
        use_mac   = input$use_mac,
        mac_dose  = input$mac_dose,
        use_trep  = input$use_trep,
        trep_dose = input$trep_dose,
        pvr_fixed_init = input$pvr_fixed_init,
        pvr_var_init   = input$pvr_var_init,
        mpap_init  = input$mpap_init,
        co_init    = input$co_init,
        sao2_init  = input$sao2_init,
        bnp_init   = input$bnp_init,
        sixmwd_init = input$sixmwd_init,
        tb_init    = input$tb_init,
        pea_done   = input$pea_done,
        bpa_done   = input$bpa_done,
        bpa_sessions = if(input$bpa_done) input$bpa_sessions else 0,
        ac_effect  = input$ac_effect,
        tend_wk    = input$tend_wk
      )
    })
  }, ignoreNULL = FALSE)

  ## ---- Default data at startup -----
  default_data <- reactive({
    run_sim(
      use_rio=TRUE, rio_dose=2.5, rio_freq=3,
      use_mac=FALSE, mac_dose=10,
      use_trep=FALSE, trep_dose=0.5,
      pvr_fixed_init=550, pvr_var_init=350,
      mpap_init=46, co_init=3.6, sao2_init=90,
      bnp_init=180, sixmwd_init=342, tb_init=0.75,
      pea_done=FALSE, bpa_done=FALSE, bpa_sessions=0,
      ac_effect=TRUE, tend_wk=52
    )
  })

  get_data <- reactive({
    if (input$run_sim == 0) default_data() else sim_data()
  })

  ## ---- TAB 2: PK plots --------------------------------------
  output$pk_plot <- renderPlotly({
    df <- get_data()
    pk_long <- df %>%
      select(time_weeks, Cp_RIO_out, Cm_RIO_out, Cp_MAC_out, Cm_MAC_out, Cp_TREP_out) %>%
      filter(time_weeks <= 2) %>%
      pivot_longer(-time_weeks, names_to="drug", values_to="conc") %>%
      filter(conc > 0.001)
    pk_long$drug <- factor(pk_long$drug,
      levels=c("Cp_RIO_out","Cm_RIO_out","Cp_MAC_out","Cm_MAC_out","Cp_TREP_out"),
      labels=c("Riociguat","M1 Metabolite","Macitentan","ACT-132577","Treprostinil"))
    p <- ggplot(pk_long, aes(x=time_weeks*168, y=conc, color=drug)) +
      geom_line(linewidth=1) +
      labs(title="PK Profiles (First 2 Weeks)", x="Time (hours)", y="Concentration (ng/mL)", color="") +
      theme_minimal()
    ggplotly(p)
  })

  output$pk_params_rio <- renderTable({
    data.frame(
      Parameter = c("CL (L/h)","Vc (L)","Vp (L)","Q (L/h)","ka (1/h)","F (%)","Tmax (h)"),
      Value     = c("2.4","30.0","20.0","1.8","0.9","94","1.5")
    )
  })
  output$pk_params_mac <- renderTable({
    data.frame(
      Parameter = c("CL (L/h)","Vc (L)","Vp (L)","Q (L/h)","ka (1/h)","F (%)","Tmax (h)"),
      Value     = c("1.1","50.0","40.0","0.8","0.35","75","8")
    )
  })
  output$pk_params_trep <- renderTable({
    data.frame(
      Parameter = c("CL (L/h)","Vc (L)","ka (1/h)","F_SC (%)","t1/2 (h)"),
      Value     = c("4.0","14.0","2.0","79","2.4")
    )
  })

  ## ---- TAB 3: PD plots --------------------------------------
  output$cgmp_plot <- renderPlotly({
    df <- get_data()
    p <- ggplot(df, aes(x=time_weeks)) +
      geom_line(aes(y=cGMP, color="cGMP"), linewidth=1.1) +
      geom_line(aes(y=cAMP, color="cAMP"), linewidth=1.1) +
      geom_hline(yintercept=1.5, linetype="dashed", color="gray") +
      scale_color_manual(values=c("cGMP"="#2980B9","cAMP"="#148F77")) +
      labs(title="Second Messengers", x="Weeks", y="pmol/mL", color="") +
      theme_minimal()
    ggplotly(p)
  })
  output$et1_plot <- renderPlotly({
    df <- get_data()
    p <- ggplot(df, aes(x=time_weeks, y=ET1)) +
      geom_line(color="#8E44AD", linewidth=1.2) +
      geom_hline(yintercept=3.5, linetype="dashed", color="gray") +
      labs(title="Plasma ET-1", x="Weeks", y="ET-1 (pg/mL)") + theme_minimal()
    ggplotly(p)
  })
  output$pvr_decomp_plot <- renderPlotly({
    df <- get_data() %>% select(time_weeks, PVR_fixed, PVR_var) %>%
      pivot_longer(-time_weeks, names_to="component", values_to="value")
    p <- ggplot(df, aes(x=time_weeks, y=value, fill=component)) +
      geom_area(alpha=0.7, position="stack") +
      scale_fill_manual(values=c("PVR_fixed"="#C0392B","PVR_var"="#E74C3C")) +
      labs(title="PVR: Fixed vs Variable", x="Weeks", y="dyn·s/cm⁵", fill="") +
      theme_minimal()
    ggplotly(p)
  })
  output$tb_plot <- renderPlotly({
    df <- get_data()
    p <- ggplot(df, aes(x=time_weeks, y=TB)) +
      geom_line(color="#E67E22", linewidth=1.2) +
      ylim(0,1) +
      labs(title="Thrombotic Burden (0=none, 1=max)", x="Weeks", y="TB (normalized)") +
      theme_minimal()
    ggplotly(p)
  })

  ## ---- TAB 4: Hemodynamics ----------------------------------
  output$mpap_plot <- renderPlotly({
    df <- get_data()
    p <- ggplot(df, aes(x=time_weeks, y=mPAP)) +
      geom_line(color="#E74C3C", linewidth=1.2) +
      geom_hline(yintercept=20, linetype="dashed", color="#27AE60", alpha=0.7) +
      labs(title="Mean PAP", x="Weeks", y="mPAP (mmHg)") + theme_minimal()
    ggplotly(p)
  })
  output$pvr_total_plot <- renderPlotly({
    df <- get_data() %>% mutate(PVR_total = PVR_fixed + PVR_var)
    p <- ggplot(df, aes(x=time_weeks, y=PVR_total)) +
      geom_line(color="#C0392B", linewidth=1.2) +
      labs(title="Total PVR", x="Weeks", y="PVR (dyn·s/cm⁵)") + theme_minimal()
    ggplotly(p)
  })
  output$co_plot <- renderPlotly({
    df <- get_data()
    p <- ggplot(df, aes(x=time_weeks, y=CO)) +
      geom_line(color="#2980B9", linewidth=1.2) +
      labs(title="Cardiac Output", x="Weeks", y="CO (L/min)") + theme_minimal()
    ggplotly(p)
  })
  output$sao2_plot <- renderPlotly({
    df <- get_data()
    p <- ggplot(df, aes(x=time_weeks, y=SaO2)) +
      geom_line(color="#27AE60", linewidth=1.2) +
      ylim(70,100) +
      labs(title="Arterial O₂ Saturation", x="Weeks", y="SaO₂ (%)") + theme_minimal()
    ggplotly(p)
  })
  output$rvwork_plot <- renderPlotly({
    df <- get_data()
    p <- ggplot(df, aes(x=time_weeks, y=RV_work)) +
      geom_line(color="#F39C12", linewidth=1.2) +
      labs(title="RV Stroke Work Index", x="Weeks", y="RVSWI (g·m/m²)") + theme_minimal()
    ggplotly(p)
  })
  output$hemo_summary_table <- renderTable({
    df <- get_data()
    wk <- input$summary_week
    row <- df %>% filter(abs(time_weeks - wk) == min(abs(time_weeks - wk))) %>% slice(1)
    data.frame(
      Parameter   = c("mPAP (mmHg)","PVR_fixed (dyn·s/cm⁵)","PVR_var (dyn·s/cm⁵)",
                      "Total PVR (dyn·s/cm⁵)","CO (L/min)","SaO₂ (%)","RV Work"),
      Value       = round(c(row$mPAP, row$PVR_fixed, row$PVR_var,
                            row$PVR_fixed+row$PVR_var,
                            row$CO, row$SaO2, row$RV_work), 1)
    )
  })

  ## ---- TAB 5: Clinical Endpoints ----------------------------
  output$sixmwd_plot <- renderPlotly({
    df <- get_data()
    p <- ggplot(df, aes(x=time_weeks, y=sixMWD)) +
      geom_line(color="#27AE60", linewidth=1.2) +
      labs(title="6MWD", x="Weeks", y="6MWD (m)") + theme_minimal()
    ggplotly(p)
  })
  output$bnp_plot <- renderPlotly({
    df <- get_data()
    p <- ggplot(df, aes(x=time_weeks, y=BNP)) +
      geom_line(color="#F39C12", linewidth=1.2) +
      labs(title="BNP", x="Weeks", y="BNP (pg/mL)") + theme_minimal()
    ggplotly(p)
  })
  output$whofc_plot <- renderPlotly({
    df <- get_data()
    p <- ggplot(df, aes(x=time_weeks, y=WHO_FC)) +
      geom_step(color="#2C3E50", linewidth=1.2) +
      scale_y_continuous(breaks=1:4, limits=c(1,4.5)) +
      labs(title="WHO Functional Class", x="Weeks", y="WHO FC") + theme_minimal()
    ggplotly(p)
  })
  output$waterfall_plot <- renderPlotly({
    df <- get_data()
    week16 <- df %>% filter(abs(time_weeks - 16) == min(abs(time_weeks - 16))) %>% slice(1)
    endpoints <- data.frame(
      Endpoint = c("ΔmPAP","ΔPVR","ΔCO","ΔSaO₂","ΔBNP","Δ6MWD"),
      Change   = c(
        round(week16$mPAP - input$mpap_init, 1),
        round((week16$PVR_fixed + week16$PVR_var) - (input$pvr_fixed_init + input$pvr_var_init), 0),
        round(week16$CO - input$co_init, 2),
        round(week16$SaO2 - input$sao2_init, 1),
        round(week16$BNP - input$bnp_init, 0),
        round(week16$sixMWD - input$sixmwd_init, 0)
      )
    )
    p <- ggplot(endpoints, aes(x=Endpoint, y=Change, fill=Change>0)) +
      geom_col() +
      scale_fill_manual(values=c("TRUE"="#27AE60","FALSE"="#E74C3C"),
                        labels=c("Worsening","Improvement")) +
      labs(title="Waterfall: Change from Baseline (Week 16)", x="", y="Change", fill="") +
      theme_minimal()
    ggplotly(p)
  })
  output$endpoints_table <- renderDT({
    df <- get_data()
    wk16 <- df %>% filter(abs(time_weeks-16)==min(abs(time_weeks-16))) %>% slice(1)
    wk52 <- df %>% filter(abs(time_weeks-52)==min(abs(time_weeks-52))) %>% slice(1)
    tbl <- data.frame(
      Endpoint     = c("mPAP (mmHg)","Total PVR (dyn·s/cm⁵)","CO (L/min)",
                       "SaO₂ (%)","BNP (pg/mL)","6MWD (m)","WHO FC"),
      Baseline     = c(input$mpap_init, input$pvr_fixed_init+input$pvr_var_init,
                       input$co_init, input$sao2_init, input$bnp_init,
                       input$sixmwd_init, 3.0),
      Wk16         = round(c(wk16$mPAP, wk16$PVR_fixed+wk16$PVR_var,
                             wk16$CO, wk16$SaO2, wk16$BNP, wk16$sixMWD, wk16$WHO_FC),1),
      Wk52         = round(c(wk52$mPAP, wk52$PVR_fixed+wk52$PVR_var,
                             wk52$CO, wk52$SaO2, wk52$BNP, wk52$sixMWD, wk52$WHO_FC),1)
    )
    tbl$`Δ Wk16` <- round(tbl$Wk16 - tbl$Baseline, 1)
    tbl$`Δ Wk52` <- round(tbl$Wk52 - tbl$Baseline, 1)
    datatable(tbl, options=list(pageLength=10, dom='t'), rownames=FALSE)
  })

  ## ---- TAB 6: Scenario Comparison ---------------------------
  compare_data <- eventReactive(input$run_compare, {
    withProgress(message="Running all scenarios...", {
      scen_defs <- list(
        s1 = list(rio=FALSE, mac=FALSE, trep=FALSE, pea=FALSE, bpa=FALSE, bpa_n=0, ac=FALSE),
        s2 = list(rio=FALSE, mac=FALSE, trep=FALSE, pea=FALSE, bpa=FALSE, bpa_n=0, ac=TRUE),
        s3 = list(rio=TRUE,  mac=FALSE, trep=FALSE, pea=FALSE, bpa=FALSE, bpa_n=0, ac=TRUE),
        s4 = list(rio=FALSE, mac=TRUE,  trep=FALSE, pea=FALSE, bpa=FALSE, bpa_n=0, ac=TRUE),
        s5 = list(rio=TRUE,  mac=TRUE,  trep=FALSE, pea=FALSE, bpa=FALSE, bpa_n=0, ac=TRUE),
        s6 = list(rio=TRUE,  mac=FALSE, trep=FALSE, pea=FALSE, bpa=TRUE,  bpa_n=5, ac=TRUE),
        s7 = list(rio=TRUE,  mac=TRUE,  trep=FALSE, pea=TRUE,  bpa=FALSE, bpa_n=0, ac=TRUE)
      )
      scen_labels <- c(
        s1="1. No Treatment", s2="2. Anticoagulation Only",
        s3="3. Riociguat 2.5mg TID", s4="4. Macitentan 10mg QD",
        s5="5. Riociguat + Macitentan", s6="6. BPA (5x) + Riociguat",
        s7="7. Post-PEA + Combination"
      )
      pea_inits <- list(
        s7 = list(pvr_fixed_init=165, mpap_init=30, co_init=4.5, sao2_init=95,
                  bnp_init=60, sixmwd_init=440)
      )
      sel <- input$scenarios_sel
      all_data <- lapply(sel, function(s) {
        d  <- scen_defs[[s]]
        pi <- if (!is.null(pea_inits[[s]])) pea_inits[[s]] else list(
          pvr_fixed_init=550, mpap_init=46, co_init=3.6,
          sao2_init=90, bnp_init=180, sixmwd_init=342
        )
        df <- run_sim(
          use_rio=d$rio, rio_dose=2.5, rio_freq=3,
          use_mac=d$mac, mac_dose=10,
          use_trep=d$trep, trep_dose=0.5,
          pvr_fixed_init=pi$pvr_fixed_init, pvr_var_init=350,
          mpap_init=pi$mpap_init, co_init=pi$co_init,
          sao2_init=pi$sao2_init, bnp_init=pi$bnp_init,
          sixmwd_init=pi$sixmwd_init, tb_init=0.75,
          pea_done=d$pea, bpa_done=d$bpa, bpa_sessions=d$bpa_n,
          ac_effect=d$ac, tend_wk=52
        )
        df$scenario <- scen_labels[s]
        df$PVR_total <- df$PVR_fixed + df$PVR_var
        df
      })
      bind_rows(all_data)
    })
  })

  output$compare_plot <- renderPlotly({
    df <- compare_data()
    req(nrow(df) > 0)
    var <- input$compare_var
    p <- ggplot(df, aes_string(x="time_weeks", y=var, color="scenario")) +
      geom_line(linewidth=1.1) +
      labs(title=paste("Scenario Comparison:", var), x="Weeks", y=var, color="Scenario") +
      theme_minimal() +
      theme(legend.position="bottom")
    ggplotly(p)
  })
  output$compare_table <- renderDT({
    df <- compare_data()
    req(nrow(df) > 0)
    df %>%
      filter(abs(time_weeks-52)==min(abs(time_weeks-52))) %>%
      group_by(scenario) %>% slice(1) %>%
      mutate(PVR_total=PVR_fixed+PVR_var) %>%
      select(scenario, mPAP, PVR_total, CO, SaO2, BNP, sixMWD, WHO_FC) %>%
      ungroup() %>%
      mutate_if(is.numeric, round, 1) %>%
      datatable(options=list(pageLength=10, dom='t'), rownames=FALSE)
  })

  ## ---- TAB 7: Biomarkers ------------------------------------
  output$biomarker_plot <- renderPlotly({
    df <- get_data()
    bm_long <- df %>%
      select(time_weeks, BNP, ET1, TB) %>%
      pivot_longer(-time_weeks, names_to="biomarker", values_to="value")
    p <- ggplot(bm_long, aes(x=time_weeks, y=value, color=biomarker)) +
      geom_line(linewidth=1.1) +
      facet_wrap(~biomarker, scales="free_y") +
      scale_color_manual(values=c("BNP"="#F39C12","ET1"="#8E44AD","TB"="#E74C3C")) +
      labs(title="Biomarker Trajectories", x="Weeks", y="Value", color="") +
      theme_minimal()
    ggplotly(p)
  })
  output$risk_profile_ui <- renderUI({
    df <- get_data()
    wk16 <- df %>% filter(abs(time_weeks-16)==min(abs(time_weeks-16))) %>% slice(1)
    risk_color <- function(val, low_thresh, high_thresh, low_good=FALSE) {
      if (low_good) {
        if (val <= low_thresh) return("green") else if (val >= high_thresh) return("red") else return("orange")
      } else {
        if (val >= low_thresh) return("green") else if (val <= high_thresh) return("red") else return("orange")
      }
    }
    fc_color <- risk_color(wk16$WHO_FC, 2, 3, TRUE)
    sixmwd_color <- risk_color(wk16$sixMWD, 440, 165)
    bnp_color <- risk_color(wk16$BNP, 50, 300, TRUE)
    tags$div(
      tags$p(style=paste0("color:",fc_color,"; font-weight:bold;"),
             paste("WHO FC:", round(wk16$WHO_FC, 1))),
      tags$p(style=paste0("color:",sixmwd_color,"; font-weight:bold;"),
             paste("6MWD:", round(wk16$sixMWD, 0), "m")),
      tags$p(style=paste0("color:",bnp_color,"; font-weight:bold;"),
             paste("BNP:", round(wk16$BNP, 0), "pg/mL")),
      tags$p(paste("mPAP:", round(wk16$mPAP, 1), "mmHg")),
      tags$p(paste("CO:", round(wk16$CO, 2), "L/min"))
    )
  })
  output$biomarker_week16_plot <- renderPlotly({
    df <- get_data()
    wk0  <- df %>% slice(1)
    wk16 <- df %>% filter(abs(time_weeks-16)==min(abs(time_weeks-16))) %>% slice(1)
    bm_df <- data.frame(
      Biomarker = c("BNP","ET-1","TB","mPAP"),
      Baseline  = c(wk0$BNP, wk0$ET1, wk0$TB*100, wk0$mPAP),
      Week16    = c(wk16$BNP, wk16$ET1, wk16$TB*100, wk16$mPAP)
    ) %>%
      pivot_longer(-Biomarker, names_to="Timepoint", values_to="Value")
    p <- ggplot(bm_df, aes(x=Biomarker, y=Value, fill=Timepoint)) +
      geom_col(position="dodge") +
      scale_fill_manual(values=c("Baseline"="#95A5A6","Week16"="#2980B9")) +
      labs(title="Biomarkers: Baseline vs Week 16", x="", y="Value", fill="") +
      theme_minimal()
    ggplotly(p)
  })
}

## ---- Launch app -----------------------------------------------
if (interactive()) {
  shinyApp(ui = ui, server = server)
}
