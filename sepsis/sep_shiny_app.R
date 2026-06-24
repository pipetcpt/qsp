## =============================================================================
## Sepsis & Septic Shock — Shiny Interactive Dashboard
## =============================================================================
## 8 Tabs: Patient Profile · Antibiotic PK · Cytokine/Immune PD ·
##         Hemodynamics & SOFA · Organ Function · Scenario Comparison ·
##         Biomarker Explorer · About
## =============================================================================

library(shiny)
library(shinydashboard)
library(ggplot2)
library(dplyr)
library(tidyr)
library(mrgsolve)

## ─────────────────────────────────────────────────────────────────────────────
## EMBEDDED mrgsolve MODEL (minimal inline)
## ─────────────────────────────────────────────────────────────────────────────
sepsis_code <- '
$PARAM
kgrow=1.2 kdeath0=0.05 Bmax=1e8
CL_abx=10 V1_abx=15 V2_abx=20 Q_abx=8 MIC=0.5 kmax_abx=5 EC50_abx=2
kprod_TNF=0.8 kdeg_TNF=0.6 kAmp_TNF=0.3 TNF0=0.02
kprod_IL6=1.2 kdeg_IL6=0.15 EC50_IL6=50
kprod_IL10=0.5 kdeg_IL10=0.2 EC50_IL10=5
kprod_IL1b=0.6 kdeg_IL1b=0.8
Neut_blood0=5000 kprod_N=200 kmarg_N=0.4 krestore_N=0.3 kdeath_N=0.08
Mac0=100 kact_Mac=0.5 kdact_Mac=0.1 kprod_Mac=10
kprod_C5a=0.3 kdeg_C5a=1.2 C5a_0=0.1
kprod_Thr=0.4 kdeg_Thr=0.5 kprod_Fib=0.2 kdeg_Fib=0.1 kprod_PAI=0.3 kdeg_PAI=0.2
kdam_End=0.15 krep_End=0.05
PF_base=400 kdam_lung=0.02 krep_lung=0.01
Cr_base=0.9 kprod_Cr=0.05 kclr_Cr=0.04
Bil_base=0.5 kprod_Bil=0.03 kclr_Bil=0.04
kprod_Lac=0.08 kclr_Lac=0.12
MAP_base=90 kdam_MAP=0.25 krep_MAP=0.15
Plt_base=250 kcons_Plt=0.04 kregen_Plt=0.02
CL_NE=120 V1_NE=18 Emax_NE=30 EC50_NE=0.05
CL_HC=25 V1_HC=40 Emax_HC=0.65 EC50_HC=200
CL_Toci=0.28 V1_Toci=4.2 Emax_Toci=0.90 EC50_Toci=1.5
useAbx=0 useNE=0 useHC=0 useToci=0 useFluid=0 FluidBoost=0
B0=1e5

$CMT BACT ABX1 ABX2 TNF IL6 IL10 IL1B NEUT_B NEUT_T MACS C5A
     THROMBIN FIBRIN PAI1 ENDOT PF_RATIO CREATININE BILIRUBIN LACTATE MAP_val PLT_COUNT
     NE_C HC_C TOCI_C

$INIT
BACT=1e5 ABX1=0 ABX2=0 TNF=0.02 IL6=0.01 IL10=0.02 IL1B=0.01
NEUT_B=5000 NEUT_T=200 MACS=10 C5A=0.1
THROMBIN=1.0 FIBRIN=1.0 PAI1=1.0 ENDOT=0.0
PF_RATIO=400 CREATININE=0.9 BILIRUBIN=0.5 LACTATE=1.0 MAP_val=90.0 PLT_COUNT=250.0
NE_C=0 HC_C=0 TOCI_C=0

$MAIN
double Cp_abx = ABX1/V1_abx;
double dC_abx = (Cp_abx > MIC) ? (Cp_abx-MIC) : 0.0;
double Ekill   = kmax_abx * pow(dC_abx,2.0)/(pow(EC50_abx,2.0)+pow(dC_abx,2.0));
double Inh_IL10 = 1.0/(1.0+IL10/5.0);
double E_HC   = Emax_HC*HC_C/(EC50_HC+HC_C);
double E_Toci = Emax_Toci*TOCI_C/(EC50_Toci+TOCI_C);
double E_NE   = Emax_NE*NE_C/(EC50_NE+NE_C);
double Inh_comb = (1.0-E_HC)*Inh_IL10;
double NFkB = BACT/(BACT+1e5);
double CI = (TNF/10.0+IL6/200.0+IL1B/5.0)/3.0;
double CI_norm = CI/(1.0+CI);

$ODE
double BactPh = 0.002*NEUT_T*BACT/(BACT+1e4);
double BactAK = useAbx*Ekill*BACT;
dxdt_BACT = kgrow*BACT*(1.0-BACT/Bmax) - kdeath0*BACT - BactPh - BactAK;
if(BACT<1.0) dxdt_BACT = -BACT;

dxdt_ABX1 = -CL_abx*Cp_abx - Q_abx*(ABX1/V1_abx-ABX2/V2_abx);
dxdt_ABX2 =  Q_abx*(ABX1/V1_abx-ABX2/V2_abx);

dxdt_TNF  = kprod_TNF*NFkB*Inh_comb + kAmp_TNF*TNF*Inh_comb - kdeg_TNF*(TNF-TNF0);
if(TNF<0) dxdt_TNF=0;
dxdt_IL6  = kprod_IL6*TNF/(EC50_IL6+TNF)*Inh_comb*(1.0-E_Toci) - kdeg_IL6*(IL6-0.01);
if(IL6<0) dxdt_IL6=0;
dxdt_IL10 = kprod_IL10*TNF/(EC50_IL10+TNF) - kdeg_IL10*(IL10-0.02);
if(IL10<0) dxdt_IL10=0;
dxdt_IL1B = kprod_IL1b*NFkB*Inh_comb - kdeg_IL1b*(IL1B-0.01);
if(IL1B<0) dxdt_IL1B=0;

double NR = kmarg_N*NEUT_B*IL6/(50.0+IL6);
dxdt_NEUT_B = kprod_N*(1.0+2.0*NFkB) - kdeath_N*NEUT_B - NR;
dxdt_NEUT_T = NR - kdeath_N*NEUT_T*(1.0+CI_norm);
dxdt_MACS   = kact_Mac*NFkB*(1.0-E_HC)*Mac0 + kprod_Mac - kdact_Mac*MACS;
dxdt_C5A    = kprod_C5a*NFkB - kdeg_C5a*(C5A-C5a_0);

dxdt_THROMBIN = kprod_Thr*CI_norm*(1.0+ENDOT) - kdeg_Thr*(THROMBIN-1.0);
dxdt_FIBRIN   = kprod_Fib*THROMBIN - kdeg_Fib*(1.0/PAI1)*FIBRIN;
dxdt_PAI1     = kprod_PAI*TNF/(5.0+TNF) - kdeg_PAI*(PAI1-1.0);

double EndD = kdam_End*CI_norm*(1.0-ENDOT);
double EndR = krep_End*ENDOT*(1.0-CI_norm);
dxdt_ENDOT  = EndD - EndR;
if(ENDOT>1.0) dxdt_ENDOT=0;
if(ENDOT<0.0) dxdt_ENDOT=0;

dxdt_PF_RATIO  = krep_lung*(PF_base-PF_RATIO)*(CI_norm<0.2?1.0:0.0)
                 - kdam_lung*CI_norm*ENDOT*PF_RATIO;
if(PF_RATIO<50) dxdt_PF_RATIO=0;

dxdt_CREATININE = kprod_Cr*(1.0+3.0*CI_norm*ENDOT) - kclr_Cr*CREATININE;
dxdt_BILIRUBIN  = kprod_Bil*(1.0+2.0*CI_norm) - kclr_Bil*BILIRUBIN;
dxdt_LACTATE    = kprod_Lac*CI_norm*(1.0+ENDOT) - kclr_Lac*LACTATE;

double NE_eff = useNE*E_NE;
double HC_map = useHC*0.3*E_HC;
double Fl_eff = useFluid*FluidBoost*2.0;
dxdt_MAP_val = NE_eff + HC_map + Fl_eff
               + krep_MAP*(MAP_base-MAP_val)*0.1
               - kdam_MAP*CI_norm*ENDOT;
if(MAP_val<30) dxdt_MAP_val=0;

dxdt_PLT_COUNT = kregen_Plt*(Plt_base-PLT_COUNT)*0.1 - kcons_Plt*FIBRIN*PLT_COUNT/200.0;
if(PLT_COUNT<10) dxdt_PLT_COUNT=0;

dxdt_NE_C  = -CL_NE/V1_NE*NE_C;
dxdt_HC_C  = -CL_HC/V1_HC*HC_C;
dxdt_TOCI_C= -CL_Toci/V1_Toci*TOCI_C;

$TABLE
int sl=0,sr=0,sliv=0,sc=0,sco=0,scns=0;
if(PF_RATIO>=400) sl=0; else if(PF_RATIO>=300) sl=1;
  else if(PF_RATIO>=200) sl=2; else if(PF_RATIO>=100) sl=3; else sl=4;
if(CREATININE<1.2) sr=0; else if(CREATININE<2.0) sr=1;
  else if(CREATININE<3.5) sr=2; else if(CREATININE<5.0) sr=3; else sr=4;
if(BILIRUBIN<1.2) sliv=0; else if(BILIRUBIN<2.0) sliv=1;
  else if(BILIRUBIN<6.0) sliv=2; else if(BILIRUBIN<12.0) sliv=3; else sliv=4;
if(MAP_val>=70) sc=0; else if(MAP_val>=65) sc=1; else if(useNE==1&&NE_C>0.005) sc=3; else sc=2;
if(PLT_COUNT>=150) sco=0; else if(PLT_COUNT>=100) sco=1;
  else if(PLT_COUNT>=50) sco=2; else if(PLT_COUNT>=20) sco=3; else sco=4;
double em=LACTATE+(TNF/10.0+IL6/200.0+IL1B/5.0)/3.0/(1.0+(TNF/10.0+IL6/200.0+IL1B/5.0)/3.0)*3.0;
if(em<1.5) scns=0; else if(em<2.5) scns=1;
  else if(em<4.0) scns=2; else if(em<6.0) scns=3; else scns=4;
int SOFA = sl+sr+sliv+sc+sco+scns;
double MortProb = 1.0/(1.0+exp(-(-6.5+0.45*SOFA)));
double logBACT = (BACT>0)?log10(BACT):0.0;
double shock_flag = ((MAP_val<65.0)&&(LACTATE>2.0))?1.0:0.0;
double IL6_pg = IL6*1000.0;
double CpAbx = ABX1/V1_abx;

$CAPTURE BACT ABX1 ABX2 TNF IL6 IL10 IL1B NEUT_B NEUT_T MACS C5A
         THROMBIN FIBRIN PAI1 ENDOT PF_RATIO CREATININE BILIRUBIN LACTATE MAP_val PLT_COUNT
         NE_C HC_C TOCI_C SOFA MortProb logBACT shock_flag IL6_pg CpAbx
         sl sr sliv sc sco scns
'

suppressMessages({
  base_mod <- mcode("sepsis_shiny", sepsis_code)
})

## ─────────────────────────────────────────────────────────────────────────────
## SIMULATION FUNCTION
## ─────────────────────────────────────────────────────────────────────────────
run_sim <- function(
    weight = 70, severity = "moderate",
    use_abx = TRUE, abx_delay = 0, abx_dose_mg = 1000, abx_ii = 8,
    use_ne  = FALSE, ne_rate = 0.1,
    use_hc  = FALSE, hc_dose = 200,
    use_toci= FALSE, toci_dose_per_kg = 8,
    use_fluid = TRUE, fluid_vol = 30,
    sim_hours = 168
) {
  # Severity → initial bacterial load + inflammatory profile
  B_init <- switch(severity,
    "mild"     = 1e4,
    "moderate" = 1e5,
    "severe"   = 5e5,
    "septic_shock" = 2e6, 1e5)

  kgrow_adj <- if (severity == "septic_shock") 1.4 else 1.2
  toci_dose  <- toci_dose_per_kg * weight

  p_list <- list(
    useAbx = as.integer(use_abx), useNE = as.integer(use_ne),
    useHC  = as.integer(use_hc),  useToci = as.integer(use_toci),
    useFluid = as.integer(use_fluid),
    FluidBoost = fluid_vol / 10,
    B0 = B_init, kgrow = kgrow_adj
  )

  mod_run <- base_mod %>%
    param(p_list) %>%
    init(BACT = B_init)

  ev_list <- ev()
  if (use_abx)  ev_list <- ev_list + ev(amt = abx_dose_mg, cmt = "ABX1", ii = abx_ii,
                                         addl = floor(sim_hours / abx_ii) - 1,
                                         time = abx_delay)
  if (use_ne)   ev_list <- ev_list + ev(amt = ne_rate * 18, cmt = "NE_C",
                                         rate = ne_rate * 18 / 4, time = 2, to = 72)
  if (use_hc)   ev_list <- ev_list + ev(amt = hc_dose / 4, cmt = "HC_C",
                                         ii = 6, addl = 27, time = 2)
  if (use_toci) ev_list <- ev_list + ev(amt = toci_dose, cmt = "TOCI_C", time = 6)

  out <- mrgsim(mod_run, events = ev_list,
                tgrid = seq(0, sim_hours, 0.5), obsonly = TRUE)
  as.data.frame(out)
}

## ─────────────────────────────────────────────────────────────────────────────
## UI
## ─────────────────────────────────────────────────────────────────────────────
ui <- dashboardPage(
  skin = "red",
  dashboardHeader(title = "Sepsis QSP Dashboard"),

  dashboardSidebar(
    sidebarMenu(
      menuItem("Patient Profile",       tabName = "tab_patient",  icon = icon("user-md")),
      menuItem("Antibiotic PK",         tabName = "tab_abx",      icon = icon("pills")),
      menuItem("Cytokine / Immune PD",  tabName = "tab_cytokine", icon = icon("dna")),
      menuItem("Hemodynamics & SOFA",   tabName = "tab_sofa",     icon = icon("heartbeat")),
      menuItem("Organ Function",        tabName = "tab_organ",    icon = icon("lungs")),
      menuItem("Scenario Comparison",   tabName = "tab_scenario", icon = icon("chart-bar")),
      menuItem("Biomarker Explorer",    tabName = "tab_biomarker",icon = icon("vial")),
      menuItem("About",                 tabName = "tab_about",    icon = icon("info-circle"))
    )
  ),

  dashboardBody(
    tags$head(tags$style(HTML("
      .box { border-top: 3px solid #d9534f; }
      .skin-red .main-header .navbar { background-color: #c0392b; }
    "))),
    tabItems(

      ## ── TAB 1: Patient Profile ─────────────────────────────────────────────
      tabItem("tab_patient",
        fluidRow(
          box(title = "Patient Parameters", status = "danger", solidHeader = TRUE, width = 4,
            numericInput("weight", "Body Weight (kg)", 70, 40, 150, 5),
            selectInput("severity", "Sepsis Severity",
                        choices = c("Mild" = "mild", "Moderate" = "moderate",
                                    "Severe" = "severe", "Septic Shock" = "septic_shock"),
                        selected = "severe"),
            numericInput("sim_hours", "Simulation Duration (h)", 168, 48, 336, 24),
            hr(),
            h4("Treatment Bundle"),
            checkboxInput("use_abx",   "Antibiotics (Meropenem)", TRUE),
            numericInput("abx_delay", "Antibiotic Delay (h from onset)", 1, 0, 24, 0.5),
            numericInput("abx_dose",  "Antibiotic Dose (mg/dose)", 1000, 250, 2000, 250),
            checkboxInput("use_fluid", "IV Fluid Resuscitation (30 mL/kg)", TRUE),
            checkboxInput("use_ne",    "Norepinephrine (vasopressor)", FALSE),
            numericInput("ne_rate",   "NE Infusion Rate (mcg/kg/min)", 0.1, 0.01, 1.0, 0.05),
            checkboxInput("use_hc",    "Hydrocortisone 200 mg/day", FALSE),
            checkboxInput("use_toci",  "Tocilizumab (8 mg/kg IV)", FALSE),
            actionButton("run_sim", "Run Simulation", class = "btn-danger btn-block",
                         icon = icon("play"))
          ),
          box(title = "SOFA at 24h — Summary Card", status = "danger",
              solidHeader = TRUE, width = 4,
            valueBoxOutput("vbox_sofa",  width = 12),
            valueBoxOutput("vbox_mort",  width = 12),
            valueBoxOutput("vbox_shock", width = 12),
            valueBoxOutput("vbox_bact",  width = 12)
          ),
          box(title = "Clinical Context — Sepsis-3 Criteria", status = "warning",
              solidHeader = TRUE, width = 4,
            tags$table(class = "table table-sm table-bordered",
              tags$tr(tags$th("Parameter"), tags$th("Threshold"), tags$th("Significance")),
              tags$tr(tags$td("SOFA Δ≥2"), tags$td("≥2 points"), tags$td("Sepsis (organ dysfunction)")),
              tags$tr(tags$td("MAP"), tags$td("<65 mmHg"), tags$td("Vasopressor required")),
              tags$tr(tags$td("Lactate"), tags$td(">2 mmol/L"), tags$td("Septic shock")),
              tags$tr(tags$td("PCT"), tags$td(">0.5 ng/mL"), tags$td("Bacterial infection")),
              tags$tr(tags$td("SOFA ≥11"), tags$td(">11"), tags$td("~40% 28-day mortality"))
            ),
            hr(),
            tags$p(tags$b("Antibiotic timing effect on mortality (Kumar 2006):")),
            tags$ul(
              tags$li("0–1h delay: survival ~80%"),
              tags$li("Each additional hour: ~7% decrease in survival"),
              tags$li("6h+ delay: survival <50%")
            )
          )
        ),
        fluidRow(
          box(title = "Time to Key Events", status = "info", solidHeader = TRUE, width = 12,
            plotOutput("plot_overview", height = "350px")
          )
        )
      ),

      ## ── TAB 2: Antibiotic PK ──────────────────────────────────────────────
      tabItem("tab_abx",
        fluidRow(
          box(title = "Antibiotic PK Parameters", status = "warning",
              solidHeader = TRUE, width = 3,
            numericInput("abx_dose2", "Meropenem Dose (mg)", 1000, 250, 2000, 250),
            numericInput("abx_ii2", "Dosing Interval (h)", 8, 4, 12, 4),
            numericInput("abx_mic", "MIC (mg/L)", 0.5, 0.01, 8, 0.25),
            numericInput("abx_delay2", "Delay from Onset (h)", 1, 0, 24, 0.5),
            numericInput("abx_cl", "CL (L/h)", 10, 2, 30, 1),
            numericInput("abx_v1", "V1 (L)", 15, 5, 40, 1),
            hr(),
            h5("PK/PD Target Attainment:"),
            tags$p("fT>MIC target: ≥40% of interval")
          ),
          box(title = "Plasma Concentration Profile", status = "warning",
              solidHeader = TRUE, width = 9,
            plotOutput("plot_abx_pk", height = "300px"),
            hr(),
            plotOutput("plot_abx_kill", height = "250px")
          )
        ),
        fluidRow(
          box(title = "Bacterial Burden — Effect of Antibiotic Timing", status = "danger",
              solidHeader = TRUE, width = 12,
            plotOutput("plot_bact_timing", height = "350px")
          )
        )
      ),

      ## ── TAB 3: Cytokine / Immune PD ──────────────────────────────────────
      tabItem("tab_cytokine",
        fluidRow(
          box(title = "Cytokine Storm Kinetics", status = "danger",
              solidHeader = TRUE, width = 6,
            plotOutput("plot_cytokines", height = "350px")
          ),
          box(title = "Innate Immune Response", status = "success",
              solidHeader = TRUE, width = 6,
            plotOutput("plot_immune", height = "350px")
          )
        ),
        fluidRow(
          box(title = "Coagulation & Endothelial Damage", status = "warning",
              solidHeader = TRUE, width = 6,
            plotOutput("plot_coag", height = "300px")
          ),
          box(title = "Complement & Macrophage Activation", status = "info",
              solidHeader = TRUE, width = 6,
            plotOutput("plot_complement", height = "300px")
          )
        )
      ),

      ## ── TAB 4: Hemodynamics & SOFA ────────────────────────────────────────
      tabItem("tab_sofa",
        fluidRow(
          box(title = "Mean Arterial Pressure & Shock", status = "danger",
              solidHeader = TRUE, width = 6,
            plotOutput("plot_map", height = "300px")
          ),
          box(title = "Blood Lactate — Tissue Hypoperfusion", status = "warning",
              solidHeader = TRUE, width = 6,
            plotOutput("plot_lactate", height = "300px")
          )
        ),
        fluidRow(
          box(title = "SOFA Score Components Over Time", status = "danger",
              solidHeader = TRUE, width = 7,
            plotOutput("plot_sofa_components", height = "350px")
          ),
          box(title = "28-day Mortality Probability", status = "danger",
              solidHeader = TRUE, width = 5,
            plotOutput("plot_mortality", height = "350px")
          )
        )
      ),

      ## ── TAB 5: Organ Function ─────────────────────────────────────────────
      tabItem("tab_organ",
        fluidRow(
          box(title = "Respiratory — PaO2/FiO2 Ratio (ARDS)", status = "primary",
              solidHeader = TRUE, width = 6,
            plotOutput("plot_pf", height = "250px")
          ),
          box(title = "Renal — Serum Creatinine (AKI)", status = "info",
              solidHeader = TRUE, width = 6,
            plotOutput("plot_creatinine", height = "250px")
          )
        ),
        fluidRow(
          box(title = "Hepatic — Bilirubin", status = "warning",
              solidHeader = TRUE, width = 6,
            plotOutput("plot_bili", height = "250px")
          ),
          box(title = "Coagulation — Platelets & Fibrin (DIC)", status = "danger",
              solidHeader = TRUE, width = 6,
            plotOutput("plot_platelets", height = "250px")
          )
        )
      ),

      ## ── TAB 6: Scenario Comparison ────────────────────────────────────────
      tabItem("tab_scenario",
        fluidRow(
          box(title = "Select Endpoint for Comparison", status = "primary",
              solidHeader = TRUE, width = 3,
            selectInput("sc_endpoint", "Endpoint",
                        choices = c("SOFA Score" = "SOFA",
                                    "MAP (mmHg)" = "MAP_val",
                                    "Lactate (mmol/L)" = "LACTATE",
                                    "Bacterial Load (log)" = "logBACT",
                                    "TNFα (ng/mL)" = "TNF",
                                    "IL-6 (pg/mL)" = "IL6_pg",
                                    "PaO2/FiO2" = "PF_RATIO",
                                    "Creatinine" = "CREATININE",
                                    "Mortality Risk" = "MortProb")),
            selectInput("sc_severity", "Patient Severity",
                        choices = c("Moderate" = "moderate", "Severe" = "severe",
                                    "Septic Shock" = "septic_shock"),
                        selected = "severe"),
            numericInput("sc_weight", "Weight (kg)", 70, 40, 150, 5),
            actionButton("run_scenarios", "Compare All Scenarios",
                         class = "btn-primary btn-block")
          ),
          box(title = "7-Scenario Endpoint Comparison", status = "primary",
              solidHeader = TRUE, width = 9,
            plotOutput("plot_scenarios", height = "450px")
          )
        ),
        fluidRow(
          box(title = "72-hour Outcome Summary Table", status = "success",
              solidHeader = TRUE, width = 12,
            tableOutput("table_scenarios")
          )
        )
      ),

      ## ── TAB 7: Biomarker Explorer ─────────────────────────────────────────
      tabItem("tab_biomarker",
        fluidRow(
          box(title = "Biomarker Selection", status = "primary",
              solidHeader = TRUE, width = 3,
            checkboxGroupInput("bm_vars", "Biomarkers to Display",
                               choices = c("TNFα" = "TNF", "IL-6 (pg/mL)" = "IL6_pg",
                                           "IL-10" = "IL10", "IL-1β" = "IL1B",
                                           "Neutrophils (blood)" = "NEUT_B",
                                           "Tissue Neutrophils" = "NEUT_T",
                                           "Macrophages" = "MACS",
                                           "C5a (Complement)" = "C5A",
                                           "Thrombin" = "THROMBIN",
                                           "Fibrin" = "FIBRIN",
                                           "PAI-1" = "PAI1",
                                           "Endothelial Damage" = "ENDOT"),
                               selected = c("TNF", "IL6_pg", "IL10", "NEUT_B")),
            hr(),
            selectInput("bm_scale", "Y-axis Scale",
                        choices = c("Linear" = "identity", "Log10" = "log10")),
            numericInput("bm_tmax", "Max Time (h)", 72, 12, 168, 12)
          ),
          box(title = "Biomarker Trajectories", status = "primary",
              solidHeader = TRUE, width = 9,
            plotOutput("plot_biomarkers", height = "500px")
          )
        )
      ),

      ## ── TAB 8: About ──────────────────────────────────────────────────────
      tabItem("tab_about",
        fluidRow(
          box(title = "About This QSP Model", status = "danger",
              solidHeader = TRUE, width = 12,
            h3("Sepsis & Septic Shock — QSP Model"),
            tags$p("This dashboard implements a 24-compartment Quantitative Systems Pharmacology (QSP)
                    model for sepsis and septic shock, simulating the multi-system pathophysiology
                    from initial bacterial inoculation through cytokine storm, organ failure, and
                    therapeutic intervention."),
            h4("Model Architecture"),
            tags$table(class = "table table-bordered",
              tags$tr(tags$th("Module"), tags$th("Components"), tags$th("Variables")),
              tags$tr(tags$td("Bacterial Dynamics"), tags$td("Growth/kill ODE"), tags$td("BACT")),
              tags$tr(tags$td("Antibiotic PK"), tags$td("2-cmpt meropenem"), tags$td("ABX1, ABX2")),
              tags$tr(tags$td("Cytokines"), tags$td("TNFα, IL-6, IL-10, IL-1β"), tags$td("TNF, IL6, IL10, IL1B")),
              tags$tr(tags$td("Innate Immunity"), tags$td("Neutrophils (B/T), Macrophages"), tags$td("NEUT_B, NEUT_T, MACS")),
              tags$tr(tags$td("Complement"), tags$td("C5a effector"), tags$td("C5A")),
              tags$tr(tags$td("Coagulation"), tags$td("Thrombin, Fibrin, PAI-1"), tags$td("THROMBIN, FIBRIN, PAI1")),
              tags$tr(tags$td("Endothelium"), tags$td("Damage index"), tags$td("ENDOT")),
              tags$tr(tags$td("Organ Function"), tags$td("PF ratio, Cr, Bili, Lac, MAP, Plt"), tags$td("6 variables")),
              tags$tr(tags$td("Drug PK"), tags$td("NE, Hydrocortisone, Tocilizumab"), tags$td("NE_C, HC_C, TOCI_C"))
            ),
            h4("Treatment Scenarios"),
            tags$ol(
              tags$li("S1: No treatment (natural history)"),
              tags$li("S2: Antibiotics alone"),
              tags$li("S3: Antibiotics + Norepinephrine"),
              tags$li("S4: Antibiotics + NE + Fluid resuscitation (Sepsis Bundle)"),
              tags$li("S5: Bundle + Hydrocortisone (ADRENAL/APROCCHSS)"),
              tags$li("S6: Bundle + HC + Tocilizumab (REMAP-CAP)"),
              tags$li("S7: Immunocompromised patient (high inoculum)")
            ),
            h4("Key References"),
            tags$ul(
              tags$li("Singer M et al. JAMA 2016;315:801 — Sepsis-3 Definition"),
              tags$li("Rivers E et al. NEJM 2001;345:1368 — EGDT"),
              tags$li("ADRENAL Trial. NEJM 2018;378:797 — Corticosteroids in septic shock"),
              tags$li("REMAP-CAP. NEJM 2021;384:1491 — Tocilizumab in sepsis"),
              tags$li("Kumar A et al. Crit Care Med 2006;34:1589 — Antibiotic timing")
            ),
            hr(),
            tags$p(tags$em("Generated by Claude Code Routine (CCR) — 2026-06-24"),
                   tags$br(), "QSP model for educational and research purposes only.")
          )
        )
      )
    )
  )
)

## ─────────────────────────────────────────────────────────────────────────────
## SERVER
## ─────────────────────────────────────────────────────────────────────────────
server <- function(input, output, session) {

  sim_data <- eventReactive(input$run_sim, {
    withProgress(message = "Running simulation...", value = 0, {
      run_sim(
        weight    = input$weight,
        severity  = input$severity,
        use_abx   = input$use_abx,
        abx_delay = input$abx_delay,
        abx_dose_mg = input$abx_dose,
        use_ne    = input$use_ne,
        ne_rate   = input$ne_rate,
        use_hc    = input$use_hc,
        use_toci  = input$use_toci,
        use_fluid = input$use_fluid,
        sim_hours = input$sim_hours
      )
    })
  }, ignoreNULL = FALSE)

  # Default simulation on load
  default_data <- reactive({
    run_sim(severity = "severe", use_abx = TRUE, use_ne = TRUE,
            use_fluid = TRUE, sim_hours = 168)
  })

  get_data <- reactive({
    if (input$run_sim == 0) default_data() else sim_data()
  })

  # Value boxes
  output$vbox_sofa <- renderValueBox({
    d24 <- get_data() %>% filter(abs(time - 24) < 0.6) %>% slice(1)
    valueBox(d24$SOFA, "SOFA at 24h", icon = icon("chart-line"),
             color = if(d24$SOFA < 6) "green" else if(d24$SOFA < 11) "yellow" else "red")
  })

  output$vbox_mort <- renderValueBox({
    d24 <- get_data() %>% filter(abs(time - 24) < 0.6) %>% slice(1)
    valueBox(paste0(round(d24$MortProb * 100, 1), "%"), "28d Mortality Risk",
             icon = icon("percent"), color = if(d24$MortProb < 0.2) "green" else "red")
  })

  output$vbox_shock <- renderValueBox({
    d24 <- get_data() %>% filter(abs(time - 24) < 0.6) %>% slice(1)
    valueBox(if(d24$shock_flag > 0.5) "YES" else "NO", "Septic Shock at 24h",
             icon = icon("exclamation-triangle"),
             color = if(d24$shock_flag < 0.5) "green" else "red")
  })

  output$vbox_bact <- renderValueBox({
    d48 <- get_data() %>% filter(abs(time - 48) < 0.6) %>% slice(1)
    valueBox(paste0(round(d48$logBACT, 1), " log CFU/mL"), "Bacteremia at 48h",
             icon = icon("bacterium"),
             color = if(d48$logBACT < 3) "green" else "orange")
  })

  # Overview plot
  output$plot_overview <- renderPlot({
    d <- get_data()
    d_long <- d %>%
      select(time, SOFA, MAP_val, LACTATE) %>%
      pivot_longer(-time)
    ggplot(d_long, aes(x = time / 24, y = value)) +
      geom_line(aes(color = name), linewidth = 1.2) +
      geom_hline(data = data.frame(name = c("MAP_val", "LACTATE"),
                                    yint = c(65, 2)),
                 aes(yintercept = yint), linetype = "dashed", color = "firebrick") +
      facet_wrap(~name, scales = "free_y") +
      labs(x = "Day", y = "Value") +
      theme_bw(12) + theme(legend.position = "none")
  })

  # Cytokines
  output$plot_cytokines <- renderPlot({
    d <- get_data() %>% select(time, TNF, IL6, IL10, IL1B) %>%
      pivot_longer(-time)
    ggplot(d, aes(x = time, y = value, color = name)) +
      geom_line(linewidth = 1.1) +
      scale_y_log10() +
      scale_color_manual(values = c(TNF = "#E74C3C", IL6 = "#E67E22",
                                     IL10 = "#3498DB", IL1B = "#8E44AD")) +
      labs(x = "Time (h)", y = "Concentration (ng/mL, log scale)",
           title = "Cytokine Dynamics", color = "Cytokine") +
      theme_bw(12)
  })

  # Immune
  output$plot_immune <- renderPlot({
    d <- get_data() %>% select(time, NEUT_B, NEUT_T, MACS) %>%
      pivot_longer(-time)
    ggplot(d, aes(x = time, y = value, color = name)) +
      geom_line(linewidth = 1.1) +
      scale_color_manual(values = c(NEUT_B = "#27AE60", NEUT_T = "#F39C12", MACS = "#8E44AD")) +
      labs(x = "Time (h)", y = "Cell count (relative units)",
           title = "Innate Immune Cell Dynamics", color = "Cell type") +
      theme_bw(12)
  })

  # Coagulation
  output$plot_coag <- renderPlot({
    d <- get_data() %>% select(time, THROMBIN, FIBRIN, PAI1, ENDOT) %>%
      pivot_longer(-time)
    ggplot(d, aes(x = time, y = value, color = name)) +
      geom_line(linewidth = 1.1) +
      labs(x = "Time (h)", y = "Relative Units / Index",
           title = "Coagulation & Endothelial Damage", color = "Marker") +
      theme_bw(12)
  })

  # Complement
  output$plot_complement <- renderPlot({
    d <- get_data() %>% select(time, C5A, MACS) %>% pivot_longer(-time)
    ggplot(d, aes(x = time, y = value, color = name)) +
      geom_line(linewidth = 1.1) +
      labs(x = "Time (h)", y = "Level (relative)", title = "Complement / Macrophages") +
      theme_bw(12)
  })

  # MAP
  output$plot_map <- renderPlot({
    d <- get_data()
    ggplot(d, aes(x = time, y = MAP_val)) +
      geom_line(color = "#E74C3C", linewidth = 1.3) +
      geom_hline(yintercept = 65, linetype = "dashed", color = "darkred") +
      annotate("text", x = 10, y = 67, label = "Shock threshold: 65 mmHg",
               color = "darkred", size = 3.5) +
      labs(x = "Time (h)", y = "MAP (mmHg)", title = "Mean Arterial Pressure") +
      theme_bw(12)
  })

  # Lactate
  output$plot_lactate <- renderPlot({
    d <- get_data()
    ggplot(d, aes(x = time, y = LACTATE)) +
      geom_line(color = "#F39C12", linewidth = 1.3) +
      geom_hline(yintercept = 2, linetype = "dashed", color = "darkred") +
      annotate("text", x = 10, y = 2.2, label = "Septic shock: >2 mmol/L",
               color = "darkred", size = 3.5) +
      labs(x = "Time (h)", y = "Lactate (mmol/L)", title = "Blood Lactate") +
      theme_bw(12)
  })

  # SOFA components
  output$plot_sofa_components <- renderPlot({
    d <- get_data() %>%
      select(time, sl, sr, sliv, sc, sco, scns, SOFA) %>%
      rename(Lung = sl, Renal = sr, Liver = sliv,
             Cardiovascular = sc, Coagulation = sco, CNS = scns) %>%
      pivot_longer(-c(time, SOFA))
    ggplot(d, aes(x = time, y = value, fill = name)) +
      geom_area(alpha = 0.75) +
      geom_line(data = get_data(), aes(x = time, y = SOFA),
                inherit.aes = FALSE, color = "black", linewidth = 1) +
      scale_fill_brewer(palette = "Set1") +
      labs(x = "Time (h)", y = "SOFA Sub-score", fill = "Organ",
           title = "SOFA Components (Stacked)") +
      theme_bw(12)
  })

  # Mortality
  output$plot_mortality <- renderPlot({
    d <- get_data()
    ggplot(d, aes(x = time, y = MortProb * 100)) +
      geom_line(color = "#C0392B", linewidth = 1.5) +
      geom_hline(yintercept = 40, linetype = "dashed", color = "grey40") +
      annotate("text", x = 10, y = 42, label = "~40% at SOFA=11", size = 3.5) +
      labs(x = "Time (h)", y = "28-day Mortality (%)",
           title = "Estimated Mortality Probability") +
      scale_y_continuous(limits = c(0, 100)) +
      theme_bw(12)
  })

  # Organ plots
  output$plot_pf <- renderPlot({
    d <- get_data()
    ggplot(d, aes(x = time, y = PF_RATIO)) +
      geom_line(color = "#2980B9", linewidth = 1.2) +
      geom_hline(yintercept = c(100, 200, 300), linetype = "dashed",
                 color = c("red", "orange", "gold"), alpha = 0.7) +
      annotate("text", x = 5, y = 85, label = "Severe ARDS", color = "red", size = 3) +
      labs(x = "Time (h)", y = "PaO2/FiO2 (mmHg)", title = "Respiratory Function") +
      theme_bw(12)
  })

  output$plot_creatinine <- renderPlot({
    d <- get_data()
    ggplot(d, aes(x = time, y = CREATININE)) +
      geom_line(color = "#16A085", linewidth = 1.2) +
      geom_hline(yintercept = 3.5, linetype = "dashed", color = "orange") +
      labs(x = "Time (h)", y = "Creatinine (mg/dL)", title = "Renal Function (AKI)") +
      theme_bw(12)
  })

  output$plot_bili <- renderPlot({
    d <- get_data()
    ggplot(d, aes(x = time, y = BILIRUBIN)) +
      geom_line(color = "#D4AC0D", linewidth = 1.2) +
      geom_hline(yintercept = 6.0, linetype = "dashed", color = "orange") +
      labs(x = "Time (h)", y = "Bilirubin (mg/dL)", title = "Hepatic Function") +
      theme_bw(12)
  })

  output$plot_platelets <- renderPlot({
    d <- get_data() %>% select(time, PLT_COUNT, FIBRIN) %>% pivot_longer(-time)
    ggplot(d, aes(x = time, y = value, color = name)) +
      geom_line(linewidth = 1.2) +
      scale_color_manual(values = c(PLT_COUNT = "#8E44AD", FIBRIN = "#E74C3C"),
                         labels = c("Fibrin Index", "Platelets (×10⁹/L)")) +
      labs(x = "Time (h)", y = "Level", color = "Marker",
           title = "Coagulation — DIC Markers") +
      theme_bw(12)
  })

  # Antibiotic PK tab
  output$plot_abx_pk <- renderPlot({
    ev_abx <- ev(amt = input$abx_dose2, cmt = "ABX1",
                 ii = input$abx_ii2, addl = 5, time = input$abx_delay2)
    pmod <- base_mod %>% param(list(useAbx = 1, MIC = input$abx_mic,
                                     CL_abx = input$abx_cl, V1_abx = input$abx_v1))
    d_abx <- mrgsim(pmod, events = ev_abx, tgrid = seq(0, 48, 0.1),
                    obsonly = TRUE) %>% as.data.frame()

    ggplot(d_abx, aes(x = time, y = CpAbx)) +
      geom_line(color = "#E67E22", linewidth = 1.2) +
      geom_hline(yintercept = input$abx_mic, linetype = "dashed", color = "red") +
      annotate("text", x = 2, y = input$abx_mic * 1.2,
               label = paste("MIC =", input$abx_mic, "mg/L"), color = "red", size = 3.5) +
      labs(x = "Time (h)", y = "Plasma Concentration (mg/L)",
           title = "Antibiotic PK — Meropenem") +
      theme_bw(12)
  })

  output$plot_abx_kill <- renderPlot({
    ev_abx <- ev(amt = input$abx_dose2, cmt = "ABX1",
                 ii = input$abx_ii2, addl = 5, time = input$abx_delay2)
    pmod <- base_mod %>% param(list(useAbx = 1, MIC = input$abx_mic,
                                     CL_abx = input$abx_cl, V1_abx = input$abx_v1,
                                     B0 = 1e5)) %>% init(BACT = 1e5)
    d_abx <- mrgsim(pmod, events = ev_abx, tgrid = seq(0, 48, 0.1),
                    obsonly = TRUE) %>% as.data.frame()

    ggplot(d_abx, aes(x = time, y = logBACT)) +
      geom_line(color = "#C0392B", linewidth = 1.2) +
      geom_hline(yintercept = 2, linetype = "dashed", color = "#27AE60") +
      annotate("text", x = 2, y = 2.3, label = "Clearance (<100 CFU/mL)",
               color = "#27AE60", size = 3.5) +
      labs(x = "Time (h)", y = "Log₁₀ Bacterial Load",
           title = "Bacterial Kill — Antibiotic Effect") +
      theme_bw(12)
  })

  output$plot_bact_timing <- renderPlot({
    delays <- c(0, 1, 3, 6, 12)
    colors <- c("#27AE60", "#2ECC71", "#F39C12", "#E74C3C", "#922B21")
    result <- lapply(seq_along(delays), function(i) {
      ev_abx <- ev(amt = 1000, cmt = "ABX1", ii = 8, addl = 5, time = delays[i])
      d <- mrgsim(base_mod %>% param(list(useAbx = 1)) %>% init(BACT = 1e5),
                  events = ev_abx, tgrid = seq(0, 72, 0.5), obsonly = TRUE) %>%
        as.data.frame() %>%
        mutate(delay_label = paste0(delays[i], "h delay"), color = colors[i])
      d
    })
    do.call(rbind, result) %>%
      ggplot(aes(x = time, y = logBACT, color = delay_label)) +
      geom_line(linewidth = 1.1) +
      scale_color_manual(values = setNames(colors, paste0(delays, "h delay"))) +
      geom_hline(yintercept = 2, linetype = "dashed") +
      labs(x = "Time (h)", y = "Log₁₀ Bacterial Load",
           title = "Impact of Antibiotic Timing on Bacterial Clearance",
           color = "Timing") +
      theme_bw(12)
  })

  # Scenario comparison
  scenario_data <- eventReactive(input$run_scenarios, {
    withProgress(message = "Running 7 scenarios...", value = 0, {
      scen_list <- list(
        list(n="S1: No Treatment",           abx=F, ne=F, hc=F, toci=F, fl=F),
        list(n="S2: Antibiotics Only",        abx=T, ne=F, hc=F, toci=F, fl=F),
        list(n="S3: Abx + NE",               abx=T, ne=T, hc=F, toci=F, fl=F),
        list(n="S4: Bundle (Abx+NE+Fluid)",  abx=T, ne=T, hc=F, toci=F, fl=T),
        list(n="S5: Bundle + HC",             abx=T, ne=T, hc=T, toci=F, fl=T),
        list(n="S6: Bundle + HC + Toci",     abx=T, ne=T, hc=T, toci=T, fl=T),
        list(n="S7: Immunocompromised",       abx=T, ne=T, hc=F, toci=F, fl=T)
      )
      colors <- c("#C0392B","#E67E22","#F1C40F","#27AE60","#2980B9","#8E44AD","#7F8C8D")
      lapply(seq_along(scen_list), function(i) {
        s <- scen_list[[i]]
        b0 <- if(i == 7) 1e6 else 1e5
        run_sim(weight = input$sc_weight, severity = input$sc_severity,
                use_abx = s$abx, use_ne = s$ne, use_hc = s$hc,
                use_toci = s$toci, use_fluid = s$fl) %>%
          mutate(scenario = s$n, color = colors[i])
      }) %>% bind_rows()
    })
  })

  output$plot_scenarios <- renderPlot({
    req(scenario_data())
    d <- scenario_data()
    colors_map <- d %>% distinct(scenario, color) %>% deframe()
    ggplot(d, aes(x = time, y = .data[[input$sc_endpoint]],
                   color = scenario)) +
      geom_line(linewidth = 1.1) +
      scale_color_manual(values = colors_map, name = "Scenario") +
      scale_x_continuous(breaks = seq(0, 168, 24),
                         labels = paste0("Day ", 0:7)) +
      labs(x = "Time", y = input$sc_endpoint,
           title = paste("Scenario Comparison —", input$sc_endpoint)) +
      theme_bw(12) +
      theme(legend.position = "bottom", legend.text = element_text(size = 8))
  })

  output$table_scenarios <- renderTable({
    req(scenario_data())
    scenario_data() %>%
      filter(abs(time - 72) < 0.6) %>%
      group_by(scenario) %>% slice(1) %>%
      select(scenario, SOFA, MAP_val, LACTATE, logBACT,
             TNF, IL6_pg, PF_RATIO, CREATININE, MortProb) %>%
      mutate(MortProb = paste0(round(MortProb * 100, 1), "%")) %>%
      rename(Scenario = scenario, `SOFA (72h)` = SOFA,
             `MAP (mmHg)` = MAP_val, `Lactate` = LACTATE,
             `Log Bacteria` = logBACT, `TNFα (ng/mL)` = TNF,
             `IL-6 (pg/mL)` = IL6_pg, `PaO2/FiO2` = PF_RATIO,
             `Creatinine` = CREATININE, `Mortality Risk` = MortProb) %>%
      ungroup()
  }, striped = TRUE, bordered = TRUE)

  # Biomarker Explorer
  output$plot_biomarkers <- renderPlot({
    req(length(input$bm_vars) > 0)
    d <- get_data() %>%
      filter(time <= input$bm_tmax) %>%
      select(time, all_of(input$bm_vars)) %>%
      pivot_longer(-time)
    ggplot(d, aes(x = time, y = value, color = name)) +
      geom_line(linewidth = 1.1) +
      scale_y_continuous(trans = input$bm_scale) +
      facet_wrap(~name, scales = "free_y") +
      labs(x = "Time (h)", y = "Value", title = "Biomarker Trajectories") +
      theme_bw(12) + theme(legend.position = "none")
  })
}

shinyApp(ui, server)
