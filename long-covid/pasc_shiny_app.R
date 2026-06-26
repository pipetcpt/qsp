################################################################################
## Long COVID (PASC) — Shiny Dashboard
## Post-Acute Sequelae of SARS-CoV-2 Infection
##
## Tabs:
##   1. Patient Profile & Disease Severity
##   2. Pharmacokinetics
##   3. PD: Viral & Immune Biomarkers
##   4. PD: Neurological & Autonomic
##   5. Clinical Endpoints & QoL
##   6. Scenario Comparison
##   7. Virtual Population Analysis
##   8. Biomarker Deep Dive
################################################################################

library(shiny)
library(shinydashboard)
library(mrgsolve)
library(dplyr)
library(ggplot2)
library(plotly)
library(DT)
library(tidyr)
library(purrr)
library(RColorBrewer)

# ============================================================
# MODEL CODE (inline)
# ============================================================
pasc_model_code <- '
$PARAM @annotated
kViral:0.80:/day:viral replication rate
kClear:1.20:/day:viral clearance rate
kReservoir:0.05:/day:seeding rate to reservoir
kActiv:0.02:/day:reservoir reactivation
kAntigen:0.30:/day:antigen clearance
IC50_nirm:0.003:ug/mL:nirmatrelvir IC50
kIFN:0.80:/day:IFN-I induction
kdIFN:0.25:/day:IFN-I decay
kCD8exh:0.10:/day:CD8 exhaustion
kdCD8:0.05:/day:CD8 recovery
kAutoAb:0.04:/day:autoAb production
kdAutoAb:0.008:/day:autoAb decay
kIL6:0.50:/day:IL-6 production
kdIL6:0.35:/day:IL-6 degradation
kTNF:0.30:/day:TNF production
kdTNF:0.45:/day:TNF degradation
IL6_basal:0.10::IL-6 baseline
EC50_met_IL6:500:ng/mL:metformin IL-6 EC50
kFibrin:0.15:/day:fibrin formation
kdFibrin:0.08:/day:fibrin lysis
kDdimer:0.20:/day:D-dimer generation
kdDdimer:0.12:/day:D-dimer clearance
kAnticoag:2.00:/day:anticoagulant fibrin effect
kBBB:0.20:/day:BBB disruption
kdBBB:0.05:/day:BBB restoration
kMicroglia:0.40:/day:microglial activation
kdMicroglia:0.10:/day:microglial deactivation
k5HT:0.15:/day:serotonin depletion
kd5HT:0.08:/day:serotonin recovery
kSSRI_5HT:0.25:/day:SSRI 5-HT restoration
EC50_LDN:2.00:ng/mL:LDN EC50
kAutoNom:0.10:/day:autonomic dysfunction onset
kdAutoNom:0.03:/day:autonomic recovery
kPOTS:0.25:/day:POTS severity
kROS:0.30:/day:ROS accumulation
kdROS:0.20:/day:ROS scavenging
kMitoDmg:0.12:/day:mitochondrial damage
kdMitoDmg:0.04:/day:mitochondrial recovery
kLactate:0.20:/day:lactate accumulation
kdLactate:0.30:/day:lactate clearance
kMet_mito:0.60:/day:metformin mito protection
F_nirm:0.74::nirmatrelvir bioavailability
ka_nirm:1.50:/day:nirmatrelvir ka
CL_nirm:8.00:L/h:nirmatrelvir CL
Vd_nirm:105.0:L:nirmatrelvir Vd
F_met:0.55::metformin F
ka_met:0.70:/day:metformin ka
CL_met:30.0:L/h:metformin CL
Vd_met:654.0:L:metformin Vd
F_sert:0.44::sertraline F
ka_sert:0.50:/day:sertraline ka
CL_sert:2.14:L/h:sertraline CL
Vd_sert:2052.0:L:sertraline Vd
F_LDN:0.96::LDN F
ka_LDN:2.40:/day:LDN ka
CL_LDN:95.0:L/h:LDN CL
Vd_LDN:1340.0:L:LDN Vd
use_nirm:0::nirmatrelvir switch
use_met:0::metformin switch
use_sert:0::sertraline switch
use_LDN:0::LDN switch
use_anticoag:0::anticoag switch

$CMT V_GUT V_PLASMA V_RES V_AG IFN CD8_exh Auto_Ab IL6 TNF
    Fibrin Ddimer BBB Microglia Serotonin AutNom ROS MitoDmg Lactate
    A_nirm C_nirm A_met C_met A_sert C_sert A_LDN C_LDN

$MAIN
double ke_nirm = (CL_nirm*1000)/(24.0*Vd_nirm);
double ke_met  = (CL_met*1000)/(24.0*Vd_met);
double ke_sert = (CL_sert*1000)/(24.0*Vd_sert);
double ke_LDN  = (CL_LDN*1000)/(24.0*Vd_LDN);

$ODE
dxdt_A_nirm = -ka_nirm*A_nirm;
dxdt_C_nirm = use_nirm*(ka_nirm*A_nirm*F_nirm/Vd_nirm) - ke_nirm*C_nirm;
dxdt_A_met  = -ka_met*A_met;
dxdt_C_met  = use_met*(ka_met*A_met*F_met/Vd_met) - ke_met*C_met;
dxdt_A_sert = -ka_sert*A_sert;
dxdt_C_sert = use_sert*(ka_sert*A_sert*F_sert/Vd_sert) - ke_sert*C_sert;
dxdt_A_LDN  = -ka_LDN*A_LDN;
dxdt_C_LDN  = use_LDN*(ka_LDN*A_LDN*F_LDN/Vd_LDN) - ke_LDN*C_LDN;
double nirm_eff = C_nirm/(IC50_nirm+C_nirm);
double met_IL6  = 1.0 - C_met/(EC50_met_IL6+C_met);
double met_mito = C_met/(EC50_met_IL6+C_met)*kMet_mito;
double sert_5HT = C_sert/(50.0+C_sert);
double LDN_eff  = C_LDN/(EC50_LDN+C_LDN);
dxdt_V_PLASMA = kViral*(1-nirm_eff)*V_PLASMA*(1-V_PLASMA/100) - kClear*IFN*V_PLASMA - kReservoir*V_PLASMA;
dxdt_V_RES    = kReservoir*V_PLASMA - kActiv*V_RES;
dxdt_V_AG     = kActiv*V_RES - kAntigen*V_AG;
dxdt_IFN      = kIFN*V_AG*(1-IFN) - kdIFN*IFN;
dxdt_CD8_exh  = kCD8exh*IFN*(1-CD8_exh) - kdCD8*(1-IFN)*CD8_exh;
dxdt_Auto_Ab  = kAutoAb*V_AG*(1-Auto_Ab) - kdAutoAb*Auto_Ab;
dxdt_IL6      = (kIL6*(V_AG+0.5*Auto_Ab)*met_IL6) - kdIL6*(IL6-IL6_basal);
if (IL6 < IL6_basal) dxdt_IL6 = 0;
dxdt_TNF      = kTNF*(V_AG+0.3*IL6)*(1-TNF) - kdTNF*TNF;
dxdt_Fibrin   = kFibrin*(IL6*0.5+Auto_Ab*0.5)*(1-Fibrin) - kdFibrin*Fibrin - use_anticoag*kAnticoag*Fibrin;
dxdt_Ddimer   = kDdimer*Fibrin - kdDdimer*Ddimer;
dxdt_BBB      = kBBB*(IL6+TNF+Fibrin)/3.0*(1-BBB) - kdBBB*BBB;
dxdt_Microglia= kMicroglia*BBB*(1-Microglia)*(1-LDN_eff) - kdMicroglia*(1-BBB)*Microglia;
dxdt_Serotonin= kd5HT*(1-Serotonin)+kSSRI_5HT*sert_5HT*(1-Serotonin) - k5HT*Microglia*Serotonin;
dxdt_AutNom   = kAutoNom*(Auto_Ab+0.5*CD8_exh)*(1-AutNom) - kdAutoNom*AutNom;
dxdt_ROS      = kROS*(Microglia+BBB+MitoDmg)/3.0 - kdROS*ROS;
dxdt_MitoDmg  = kMitoDmg*ROS*(1-MitoDmg) - kdMitoDmg*MitoDmg - met_mito*MitoDmg*0.1;
dxdt_Lactate  = kLactate*MitoDmg - kdLactate*Lactate;

$TABLE
double FSS = 1.0+6.0*(0.4*MitoDmg+0.3*Lactate+0.2*AutNom+0.1*Microglia);
if (FSS>7.0) FSS=7.0;
double VO2max = 100.0*(1.0-0.5*MitoDmg-0.3*Fibrin-0.2*Lactate);
if (VO2max<10.0) VO2max=10.0;
double MoCA = 30.0*(1.0-0.5*(1-Serotonin)-0.3*Microglia-0.2*BBB);
if (MoCA<0) MoCA=0;
double POTS_HR = 30.0*AutNom+15.0*(1.0-Serotonin);
if (POTS_HR<0) POTS_HR=0;
double SF36_PCS = 100.0-30.0*FSS/7.0-20.0*MitoDmg-15.0*AutNom-15.0*(1-VO2max/100.0)-10.0*Fibrin-10.0*Ddimer;
if (SF36_PCS<0) SF36_PCS=0;
double NfL = 5.0+45.0*(0.5*BBB+0.3*Microglia+0.2*(1-Serotonin));
capture FSS=FSS; capture VO2max=VO2max; capture MoCA=MoCA;
capture POTS_HR=POTS_HR; capture SF36_PCS=SF36_PCS; capture NfL_pg=NfL;
capture CRP_proxy=IL6*5.0; capture Ddimer_out=Ddimer;

$CAPTURE C_nirm C_met C_sert C_LDN V_PLASMA V_AG V_RES
         IL6 TNF CD8_exh Auto_Ab Fibrin BBB Microglia Serotonin AutNom MitoDmg Lactate

$INIT V_GUT=0 V_PLASMA=0.5 V_RES=3.0 V_AG=0.8
     IFN=0.4 CD8_exh=0.55 Auto_Ab=0.35 IL6=0.5 TNF=0.3
     Fibrin=0.45 Ddimer=0.4 BBB=0.4 Microglia=0.5 Serotonin=0.55
     AutNom=0.45 ROS=0.4 MitoDmg=0.45 Lactate=0.35
     A_nirm=0 C_nirm=0 A_met=0 C_met=0 A_sert=0 C_sert=0 A_LDN=0 C_LDN=0
'

mod_global <- mcode("pasc_shiny", pasc_model_code, quiet=TRUE)

# ============================================================
# UI
# ============================================================
ui <- dashboardPage(
  skin = "blue",
  dashboardHeader(title = "Long COVID (PASC) QSP Dashboard", titleWidth=350),

  dashboardSidebar(
    width = 300,
    sidebarMenu(
      menuItem("Patient Profile",    tabName="profile",    icon=icon("user")),
      menuItem("Pharmacokinetics",   tabName="pk",         icon=icon("pills")),
      menuItem("Viral & Immune",     tabName="viral",      icon=icon("virus")),
      menuItem("Neuro & Autonomic",  tabName="neuro",      icon=icon("brain")),
      menuItem("Clinical Endpoints", tabName="clinical",   icon=icon("chart-line")),
      menuItem("Scenario Comparison",tabName="compare",    icon=icon("balance-scale")),
      menuItem("Virtual Population", tabName="vpop",       icon=icon("users")),
      menuItem("Biomarker Panel",    tabName="biomarker",  icon=icon("flask"))
    ),
    hr(),
    h4("  Treatment Options", style="color:white; padding-left:10px"),
    checkboxInput("use_nirm",    "Nirmatrelvir 300mg BID ×15d", FALSE),
    checkboxInput("use_met",     "Metformin 500mg BID",          FALSE),
    checkboxInput("use_sert",    "Sertraline 50mg QD",            FALSE),
    checkboxInput("use_LDN",     "LDN 4.5mg QD",                  FALSE),
    checkboxInput("use_anticoag","Anticoagulation (Aspirin)",      FALSE),
    hr(),
    h4("  Simulation Period", style="color:white; padding-left:10px"),
    sliderInput("sim_days", "Days:", min=30, max=730, value=365, step=30),
    hr(),
    h4("  Patient Severity", style="color:white; padding-left:10px"),
    sliderInput("init_viral_res", "Viral Reservoir:", 0.5, 8, 3, 0.5),
    sliderInput("init_mito",      "Mitochondrial Damage:", 0.1, 0.9, 0.45, 0.05),
    sliderInput("init_autnom",    "Autonomic Dysfunction:", 0.0, 0.9, 0.45, 0.05),
    sliderInput("init_autoab",    "Autoantibodies:", 0.0, 0.9, 0.35, 0.05),
    actionButton("run_sim", "Run Simulation", class="btn-success btn-block",
                 style="margin-top:10px")
  ),

  dashboardBody(
    tags$head(tags$style(HTML("
      .content-wrapper { background-color: #f5f5f5; }
      .box-header { background-color: #2c3e50 !important; color: white !important; }
      .nav-tabs-custom > .nav-tabs > li.active { border-top-color: #3498db; }
    "))),
    tabItems(

      # ---- TAB 1: Patient Profile ----
      tabItem(tabName="profile",
        fluidRow(
          valueBoxOutput("vbox_fss", width=3),
          valueBoxOutput("vbox_vo2", width=3),
          valueBoxOutput("vbox_moca", width=3),
          valueBoxOutput("vbox_pots", width=3)
        ),
        fluidRow(
          box(title="Disease Mechanism Summary", width=8, status="primary",
            plotlyOutput("radar_plot", height="400px")),
          box(title="PASC Symptom Domains", width=4, status="info",
            h4("Current Patient Status"),
            tableOutput("status_table")
          )
        ),
        fluidRow(
          box(title="About Long COVID / PASC", width=12, status="warning",
            p("Post-Acute Sequelae of SARS-CoV-2 (PASC), commonly known as Long COVID,
               affects 10-30% of COVID-19 survivors. Symptoms persist >12 weeks post-infection
               and span multiple organ systems."),
            p(strong("Key Pathological Mechanisms:")),
            tags$ul(
              tags$li("Viral persistence in tissue reservoirs (GI, lymph nodes, CNS)"),
              tags$li("Immune dysregulation: T/B cell exhaustion, autoantibody generation"),
              tags$li("Endothelial injury with fibrin microthrombi (microclots)"),
              tags$li("Neuroinflammation and BBB disruption causing brain fog"),
              tags$li("Autonomic dysfunction / POTS (postural orthostatic tachycardia)"),
              tags$li("Mitochondrial dysfunction with impaired ATP production (PEM)"),
              tags$li("Gut dysbiosis amplifying systemic inflammation")
            )
          )
        )
      ),

      # ---- TAB 2: Pharmacokinetics ----
      tabItem(tabName="pk",
        fluidRow(
          box(title="Drug Concentration-Time Profiles", width=12, status="primary",
            plotlyOutput("pk_plot", height="450px"))
        ),
        fluidRow(
          box(title="PK Parameters Summary", width=6, status="info",
            DTOutput("pk_table")),
          box(title="Drug Mechanism of Action", width=6, status="success",
            h4("Nirmatrelvir (Paxlovid)"),
            p("3CL protease (Mpro) inhibitor. IC50 = 0.003 µg/mL. Reduces viral replication,
               decreasing antigen persistence and downstream immune activation.
               Extended 15-day course targets viral reservoirs."),
            h4("Metformin"),
            p("AMPK activator. Reduces mitochondrial Complex I-driven ROS, suppresses
               NF-κB (anti-IL-6, anti-TNF), and activates mTOR inhibition.
               COVID-OUT RCT: 41% reduction in long COVID incidence (HR=0.59)."),
            h4("Sertraline (SSRI)"),
            p("SERT inhibitor + σ1R agonist. Restores serotonin signaling, reduces platelet
               hyperactivation, and has direct anti-inflammatory effects via σ1R pathway."),
            h4("Low-Dose Naltrexone (LDN)"),
            p("At 1.5–4.5 mg/day: TLR4 antagonism → microglial suppression.
               Opioid receptor transient blockade → endorphin/enkephalin rebound.
               Reduces neuroinflammation and CNS fatigue signaling.")
          )
        )
      ),

      # ---- TAB 3: Viral & Immune ----
      tabItem(tabName="viral",
        fluidRow(
          box(title="Viral Kinetics", width=6, status="danger",
            plotlyOutput("viral_plot", height="350px")),
          box(title="Immune Biomarkers", width=6, status="warning",
            plotlyOutput("immune_plot", height="350px"))
        ),
        fluidRow(
          box(title="Cytokine Dynamics", width=6, status="warning",
            plotlyOutput("cytokine_plot", height="350px")),
          box(title="Coagulation Markers", width=6, status="danger",
            plotlyOutput("coag_plot", height="350px"))
        )
      ),

      # ---- TAB 4: Neuro & Autonomic ----
      tabItem(tabName="neuro",
        fluidRow(
          box(title="Neuroinflammation Cascade", width=6, status="primary",
            plotlyOutput("neuro_plot", height="350px")),
          box(title="Autonomic Dysfunction (POTS)", width=6, status="info",
            plotlyOutput("auto_plot", height="350px"))
        ),
        fluidRow(
          box(title="Mitochondrial & Energy Metabolism", width=6, status="warning",
            plotlyOutput("mito_plot", height="350px")),
          box(title="Serotonin & Neurotransmitter Recovery", width=6, status="success",
            plotlyOutput("sero_plot", height="350px"))
        )
      ),

      # ---- TAB 5: Clinical Endpoints ----
      tabItem(tabName="clinical",
        fluidRow(
          box(title="Fatigue Severity Score (FSS)", width=6, status="danger",
            plotlyOutput("fss_plot", height="300px")),
          box(title="VO2max & Exercise Capacity", width=6, status="warning",
            plotlyOutput("vo2_plot", height="300px"))
        ),
        fluidRow(
          box(title="Cognitive Function (MoCA)", width=6, status="info",
            plotlyOutput("moca_plot", height="300px")),
          box(title="SF-36 Physical Component Score", width=6, status="success",
            plotlyOutput("sf36_plot", height="300px"))
        ),
        fluidRow(
          box(title="Week 52 Endpoint Summary", width=12,
            DTOutput("endpoint_table"))
        )
      ),

      # ---- TAB 6: Scenario Comparison ----
      tabItem(tabName="compare",
        fluidRow(
          box(title="All Scenarios — FSS Over Time", width=12, status="primary",
            plotlyOutput("scenario_fss", height="400px"))
        ),
        fluidRow(
          box(title="Scenario Comparison: Key Endpoints at Week 52", width=8,
            plotlyOutput("scenario_bar", height="400px")),
          box(title="Scenario Summary Table", width=4,
            DTOutput("scenario_table"))
        )
      ),

      # ---- TAB 7: Virtual Population ----
      tabItem(tabName="vpop",
        fluidRow(
          box(title="Virtual Population Settings", width=4, status="info",
            numericInput("n_vp", "Number of Patients:", 100, min=20, max=500, step=20),
            actionButton("run_vpop", "Run VP Simulation", class="btn-warning btn-block")
          ),
          box(title="VP: FSS Distribution at Week 52", width=8,
            plotlyOutput("vpop_fss", height="350px"))
        ),
        fluidRow(
          box(title="VP: VO2max vs FSS (Week 52)", width=6,
            plotlyOutput("vpop_scatter", height="350px")),
          box(title="VP: Response Rate Analysis", width=6,
            plotlyOutput("vpop_response", height="350px"))
        )
      ),

      # ---- TAB 8: Biomarker Panel ----
      tabItem(tabName="biomarker",
        fluidRow(
          box(title="Neurofilament Light Chain (NfL)", width=6, status="primary",
            plotlyOutput("nfl_plot", height="300px")),
          box(title="CRP (proxy from IL-6)", width=6, status="warning",
            plotlyOutput("crp_plot", height="300px"))
        ),
        fluidRow(
          box(title="D-Dimer Trajectory", width=6, status="danger",
            plotlyOutput("ddimer_plot", height="300px")),
          box(title="Biomarker Reference Ranges", width=6, status="info",
            DTOutput("biomarker_ref_table"))
        )
      )
    )
  )
)

# ============================================================
# SERVER
# ============================================================
server <- function(input, output, session) {

  # Reactive simulation
  sim_data <- eventReactive(input$run_sim, {
    p <- list(
      use_nirm    = as.numeric(input$use_nirm),
      use_met     = as.numeric(input$use_met),
      use_sert    = as.numeric(input$use_sert),
      use_LDN     = as.numeric(input$use_LDN),
      use_anticoag= as.numeric(input$use_anticoag)
    )
    init_vals <- list(
      V_RES   = input$init_viral_res,
      MitoDmg = input$init_mito,
      AutNom  = input$init_autnom,
      Auto_Ab = input$init_autoab
    )
    events <- NULL
    if (input$use_nirm)     events <- c(events, ev(cmt="A_nirm", amt=300000, ii=0.5, addl=29, time=0))
    if (input$use_met)      events <- c(events, ev(cmt="A_met",  amt=500000, ii=0.5, addl=input$sim_days*2, time=0))
    if (input$use_sert)     events <- c(events, ev(cmt="A_sert", amt=50000, ii=1, addl=input$sim_days, time=0))
    if (input$use_LDN)      events <- c(events, ev(cmt="A_LDN",  amt=4500,  ii=1, addl=input$sim_days, time=0))

    if (is.null(events)) {
      out <- mod_global %>% param(p) %>% init(init_vals) %>%
        mrgsim(end=input$sim_days, delta=1) %>% as.data.frame()
    } else {
      out <- mod_global %>% param(p) %>% init(init_vals) %>%
        mrgsim(events=events, end=input$sim_days, delta=1) %>% as.data.frame()
    }
    out
  }, ignoreNULL=FALSE)

  # All scenarios data
  all_scen_data <- reactive({
    run_scen <- function(name, pars, evts) {
      if (is.null(evts)) {
        out <- mod_global %>% param(pars) %>%
          mrgsim(end=365, delta=1) %>% as.data.frame()
      } else {
        out <- mod_global %>% param(pars) %>%
          mrgsim(events=evts, end=365, delta=1) %>% as.data.frame()
      }
      out$scenario <- name; out
    }
    p0 <- list(use_nirm=0,use_met=0,use_sert=0,use_LDN=0,use_anticoag=0)
    ev_n  <- ev(cmt="A_nirm",amt=300000,ii=0.5,addl=29,time=0)
    ev_m  <- ev(cmt="A_met", amt=500000,ii=0.5,addl=729,time=0)
    ev_s  <- ev(cmt="A_sert",amt=50000, ii=1,  addl=364,time=0)
    ev_l  <- ev(cmt="A_LDN", amt=4500,  ii=1,  addl=364,time=0)
    bind_rows(
      run_scen("S1: No Treatment",    modifyList(p0, list(use_nirm=0,use_met=0,use_sert=0,use_LDN=0)), NULL),
      run_scen("S2: Nirmatrelvir",    modifyList(p0, list(use_nirm=1)), ev_n),
      run_scen("S3: Metformin",       modifyList(p0, list(use_met=1)),  ev_m),
      run_scen("S4: LDN",             modifyList(p0, list(use_LDN=1)),  ev_l),
      run_scen("S5: Sertraline",      modifyList(p0, list(use_sert=1)), ev_s),
      run_scen("S6: Nirm+Met",        modifyList(p0, list(use_nirm=1,use_met=1)), c(ev_n,ev_m)),
      run_scen("S7: Full Combo",      modifyList(p0, list(use_nirm=1,use_met=1,use_sert=1,use_LDN=1)), c(ev_n,ev_m,ev_s,ev_l))
    )
  })

  # Value boxes
  output$vbox_fss <- renderValueBox({
    d <- sim_data(); last <- tail(d, 1)
    fss_val <- round(last$FSS, 1)
    color <- ifelse(fss_val <= 4, "green", ifelse(fss_val <= 5.5, "yellow", "red"))
    valueBox(fss_val, "FSS (1-7)", icon=icon("tired"), color=color)
  })
  output$vbox_vo2 <- renderValueBox({
    d <- sim_data(); last <- tail(d, 1)
    v <- round(last$VO2max, 0)
    color <- ifelse(v >= 80, "green", ifelse(v >= 60, "yellow", "red"))
    valueBox(paste0(v,"%"), "VO2max (% Pred)", icon=icon("heartbeat"), color=color)
  })
  output$vbox_moca <- renderValueBox({
    d <- sim_data(); last <- tail(d, 1)
    m <- round(last$MoCA, 0)
    color <- ifelse(m >= 26, "green", ifelse(m >= 22, "yellow", "red"))
    valueBox(m, "MoCA Score (/30)", icon=icon("brain"), color=color)
  })
  output$vbox_pots <- renderValueBox({
    d <- sim_data(); last <- tail(d, 1)
    p <- round(last$POTS_HR, 0)
    color <- ifelse(p < 30, "green", ifelse(p < 40, "yellow", "red"))
    valueBox(paste0("+",p," bpm"), "Orthostatic ΔHR", icon=icon("heartbeat"), color=color)
  })

  # Status table
  output$status_table <- renderTable({
    d <- sim_data(); last <- tail(d, 1)
    tibble(
      Domain = c("Fatigue","Cognition","Exercise","Autonomic","Neuroinflam","Fibrin","IL-6"),
      Status = c(round(last$FSS,2), round(last$MoCA,1), round(last$VO2max,0),
                 round(last$POTS_HR,0), round(last$Microglia,3), round(last$Fibrin,3), round(last$IL6,3)),
      Unit   = c("FSS 1-7","MoCA 0-30","% pred","bpm","0-1 index","0-1 index","norm")
    )
  })

  # Radar plot
  output$radar_plot <- renderPlotly({
    d <- sim_data(); last <- tail(d, 1)
    domains <- c("Viral Burden","Immune Dysreg","Vascular","Neuro","Autonomic","Mitochondria")
    vals <- c(
      min(1, (last$V_AG + last$V_RES/5)/2),
      min(1, (last$CD8_exh + last$Auto_Ab)/2),
      min(1, (last$Fibrin + last$Ddimer)/2),
      min(1, (last$BBB + last$Microglia + (1-last$Serotonin))/3),
      min(1, last$AutNom),
      min(1, (last$MitoDmg + last$Lactate)/2)
    )
    plot_ly(type='scatterpolar', r=c(vals, vals[1]), theta=c(domains, domains[1]),
            fill='toself', fillcolor='rgba(52,152,219,0.3)',
            line=list(color='#2980b9')) %>%
      layout(polar=list(radialaxis=list(range=c(0,1))),
             title="Disease Domain Activity (end of simulation)")
  })

  # PK plot
  output$pk_plot <- renderPlotly({
    d <- sim_data()
    p <- plot_ly(d, x=~time) %>%
      add_lines(y=~C_nirm, name="Nirmatrelvir (ng/mL)", line=list(color="#9B59B6")) %>%
      add_lines(y=~C_met/100, name="Metformin/100 (ng/mL)", line=list(color="#F1C40F", dash="dash")) %>%
      add_lines(y=~C_sert*5, name="Sertraline×5 (ng/mL)", line=list(color="#3498DB", dash="dot")) %>%
      add_lines(y=~C_LDN*50, name="LDN×50 (ng/mL)", line=list(color="#2ECC71", dash="dashdot")) %>%
      layout(xaxis=list(title="Day"), yaxis=list(title="Concentration (scaled)"),
             title="PK Profiles (scaled to common axis)", legend=list(orientation="h"))
    p
  })

  # PK table
  output$pk_table <- renderDT({
    tribble(
      ~Drug, ~Dose, ~F, ~`t½(h)`, ~`Cmax(ng/mL)`, ~Mechanism,
      "Nirmatrelvir",  "300mg BID ×15d","74%","6.9","3000","3CL protease inhibitor",
      "Metformin",     "500mg BID",     "55%","6.0","1500","AMPK activator, Complex I",
      "Sertraline",    "50mg QD",       "44%","26h","20","SERT inhibitor, σ1R agonist",
      "LDN",           "4.5mg QD",      "96%","4h","5","TLR4 antagonist, opioid rebound"
    ) %>% datatable(options=list(dom="t"), rownames=FALSE)
  })

  # Viral plot
  output$viral_plot <- renderPlotly({
    d <- sim_data()
    plot_ly(d, x=~time) %>%
      add_lines(y=~V_PLASMA, name="Viral Load (plasma)", line=list(color="#E74C3C")) %>%
      add_lines(y=~V_RES/5, name="Reservoir/5", line=list(color="#C0392B", dash="dash")) %>%
      add_lines(y=~V_AG, name="Persistent Antigen", line=list(color="#E67E22")) %>%
      layout(xaxis=list(title="Day"), yaxis=list(title="Normalized Level"),
             title="Viral Compartments", legend=list(orientation="h"))
  })

  # Immune plot
  output$immune_plot <- renderPlotly({
    d <- sim_data()
    plot_ly(d, x=~time) %>%
      add_lines(y=~IFN,     name="Type I IFN",    line=list(color="#9B59B6")) %>%
      add_lines(y=~CD8_exh, name="CD8 Exhaustion",line=list(color="#8E44AD")) %>%
      add_lines(y=~Auto_Ab, name="Autoantibodies",line=list(color="#D35400")) %>%
      layout(xaxis=list(title="Day"), yaxis=list(title="Normalized Index"),
             title="Immune Dysfunction Markers", legend=list(orientation="h"))
  })

  output$cytokine_plot <- renderPlotly({
    d <- sim_data()
    plot_ly(d, x=~time) %>%
      add_lines(y=~IL6, name="IL-6 (norm)", line=list(color="#E74C3C")) %>%
      add_lines(y=~TNF, name="TNF-α (norm)", line=list(color="#E67E22")) %>%
      layout(xaxis=list(title="Day"), yaxis=list(title="Normalized Level"),
             title="Cytokine Dynamics")
  })

  output$coag_plot <- renderPlotly({
    d <- sim_data()
    plot_ly(d, x=~time) %>%
      add_lines(y=~Fibrin,    name="Fibrin Microclots", line=list(color="#C0392B")) %>%
      add_lines(y=~Ddimer_out,name="D-Dimer",           line=list(color="#E74C3C", dash="dash")) %>%
      layout(xaxis=list(title="Day"), yaxis=list(title="Normalized Level"),
             title="Coagulation/Fibrin Markers")
  })

  output$neuro_plot <- renderPlotly({
    d <- sim_data()
    plot_ly(d, x=~time) %>%
      add_lines(y=~BBB,      name="BBB Disruption",    line=list(color="#2C3E50")) %>%
      add_lines(y=~Microglia,name="Microglial Activ.", line=list(color="#2980B9")) %>%
      layout(xaxis=list(title="Day"), yaxis=list(title="Index (0-1)"),
             title="Neuroinflammation")
  })

  output$auto_plot <- renderPlotly({
    d <- sim_data()
    plot_ly(d, x=~time) %>%
      add_lines(y=~AutNom,  name="Autonomic Dysfunction", line=list(color="#27AE60")) %>%
      add_lines(y=~POTS_HR, name="POTS ΔHR (bpm/30)",    line=list(color="#16A085")) %>%
      add_lines(y=~rep(1, nrow(d)), name="POTS threshold/30",
                line=list(color="red", dash="dash")) %>%
      layout(xaxis=list(title="Day"), yaxis=list(title="Normalized/bpm"),
             title="Autonomic Dysfunction & POTS")
  })

  output$mito_plot <- renderPlotly({
    d <- sim_data()
    plot_ly(d, x=~time) %>%
      add_lines(y=~MitoDmg, name="Mitochondrial Damage", line=list(color="#F39C12")) %>%
      add_lines(y=~Lactate,  name="Lactate (norm)",        line=list(color="#E67E22")) %>%
      add_lines(y=~ROS,      name="ROS",                   line=list(color="#D35400")) %>%
      layout(xaxis=list(title="Day"), yaxis=list(title="Index (0-1)"),
             title="Mitochondrial Dysfunction & Energy")
  })

  output$sero_plot <- renderPlotly({
    d <- sim_data()
    plot_ly(d, x=~time) %>%
      add_lines(y=~Serotonin, name="CNS Serotonin", line=list(color="#3498DB")) %>%
      add_lines(y=~rep(1, nrow(d)), name="Normal Level",
                line=list(color="green", dash="dash")) %>%
      layout(xaxis=list(title="Day"), yaxis=list(title="Normalized (0-1)"),
             title="Serotonin & Neurotransmitter Recovery")
  })

  output$fss_plot <- renderPlotly({
    d <- sim_data()
    plot_ly(d, x=~time, y=~FSS, type='scatter', mode='lines',
            line=list(color="#E74C3C", width=2)) %>%
      add_lines(y=rep(4, nrow(d)), name="Responder threshold",
                line=list(dash="dash", color="gray")) %>%
      layout(xaxis=list(title="Day"), yaxis=list(title="FSS (1-7)"),
             title="Fatigue Severity Score")
  })

  output$vo2_plot <- renderPlotly({
    d <- sim_data()
    plot_ly(d, x=~time, y=~VO2max, type='scatter', mode='lines',
            line=list(color="#F39C12", width=2)) %>%
      layout(xaxis=list(title="Day"), yaxis=list(title="VO2max (% predicted)"),
             title="Exercise Capacity (VO2max)")
  })

  output$moca_plot <- renderPlotly({
    d <- sim_data()
    plot_ly(d, x=~time, y=~MoCA, type='scatter', mode='lines',
            line=list(color="#3498DB", width=2)) %>%
      add_lines(y=rep(26, nrow(d)), name="Normal (≥26)",
                line=list(dash="dash", color="green")) %>%
      layout(xaxis=list(title="Day"), yaxis=list(title="MoCA (0-30)"),
             title="Cognitive Function (MoCA)")
  })

  output$sf36_plot <- renderPlotly({
    d <- sim_data()
    plot_ly(d, x=~time, y=~SF36_PCS, type='scatter', mode='lines',
            line=list(color="#27AE60", width=2)) %>%
      layout(xaxis=list(title="Day"), yaxis=list(title="SF-36 PCS"),
             title="Quality of Life (SF-36 Physical Component)")
  })

  output$endpoint_table <- renderDT({
    d <- sim_data() %>% filter(time %in% c(0, 90, 180, 365)) %>%
      select(time, FSS, VO2max, MoCA, POTS_HR, SF36_PCS, NfL_pg, CRP_proxy) %>%
      mutate(across(where(is.numeric), ~round(., 2))) %>%
      rename(Day=time, `NfL (pg/mL)`=NfL_pg, `CRP proxy`=CRP_proxy)
    datatable(d, options=list(dom="t"), rownames=FALSE)
  })

  # Scenario comparison
  output$scenario_fss <- renderPlotly({
    d <- all_scen_data()
    colors <- c("#E74C3C","#E67E22","#F1C40F","#2ECC71","#3498DB","#9B59B6","#1ABC9C")
    scenarios <- unique(d$scenario)
    p <- plot_ly(d, x=~time, y=~FSS, color=~scenario, type='scatter', mode='lines',
                 colors=colors) %>%
      add_lines(y=rep(4, nrow(filter(d, scenario=="S1: No Treatment"))), x=filter(d, scenario=="S1: No Treatment")$time,
                name="Responder", line=list(dash="dash",color="black"), showlegend=FALSE) %>%
      layout(xaxis=list(title="Day"), yaxis=list(title="FSS (1-7)"),
             title="Fatigue Severity: All Scenarios",
             legend=list(orientation="h", y=-0.2))
    p
  })

  output$scenario_bar <- renderPlotly({
    d <- all_scen_data() %>%
      filter(time == 364) %>%
      select(scenario, FSS, VO2max, MoCA, SF36_PCS) %>%
      pivot_longer(-scenario, names_to="Endpoint", values_to="Value")
    plot_ly(d, x=~scenario, y=~Value, color=~Endpoint, type='bar') %>%
      layout(barmode='group', xaxis=list(title="", tickangle=30),
             yaxis=list(title="Value"), title="Week 52 Endpoints by Scenario")
  })

  output$scenario_table <- renderDT({
    all_scen_data() %>%
      filter(time == 364) %>%
      select(scenario, FSS, VO2max, MoCA, POTS_HR) %>%
      mutate(across(where(is.numeric), ~round(., 1))) %>%
      datatable(options=list(dom="t"), rownames=FALSE)
  })

  # Virtual population
  vp_results <- eventReactive(input$run_vpop, {
    n <- input$n_vp
    set.seed(123)
    params <- tibble(ID=1:n,
      kViral=rlnorm(n,log(0.8),0.3), kMitoDmg=rlnorm(n,log(0.12),0.35),
      kIL6=rlnorm(n,log(0.5),0.3),   kAutoAb=rlnorm(n,log(0.04),0.4))
    inits  <- tibble(ID=1:n,
      V_RES=rlnorm(n,log(3),0.5), MitoDmg=pmin(0.95,rlnorm(n,log(0.45),0.3)),
      AutNom=pmin(0.95,rlnorm(n,log(0.45),0.35)))
    ev_all <- c(ev(cmt="A_nirm",amt=300000,ii=0.5,addl=29,time=0),
                ev(cmt="A_met", amt=500000,ii=0.5,addl=729,time=0),
                ev(cmt="A_sert",amt=50000, ii=1,  addl=364,time=0),
                ev(cmt="A_LDN", amt=4500,  ii=1,  addl=364,time=0))
    trts <- list(
      "No Treatment" = list(p=list(use_nirm=0,use_met=0,use_sert=0,use_LDN=0,use_anticoag=0), ev=NULL),
      "Full Combo"   = list(p=list(use_nirm=1,use_met=1,use_sert=1,use_LDN=1,use_anticoag=0), ev=ev_all)
    )
    map_dfr(names(trts), function(tname) {
      trt <- trts[[tname]]
      map_dfr(1:n, function(i) {
        p_i <- modifyList(trt$p, list(
          kViral=params$kViral[i], kMitoDmg=params$kMitoDmg[i],
          kIL6=params$kIL6[i], kAutoAb=params$kAutoAb[i]))
        ini <- list(V_RES=inits$V_RES[i], MitoDmg=inits$MitoDmg[i], AutNom=inits$AutNom[i])
        if (is.null(trt$ev)) {
          out <- mod_global %>% param(p_i) %>% init(ini) %>%
            mrgsim(end=365, delta=365) %>% as.data.frame() %>% tail(1)
        } else {
          out <- mod_global %>% param(p_i) %>% init(ini) %>%
            mrgsim(events=trt$ev, end=365, delta=365) %>% as.data.frame() %>% tail(1)
        }
        out$ID <- i; out$trt <- tname; out
      })
    })
  })

  output$vpop_fss <- renderPlotly({
    req(vp_results())
    d <- vp_results()
    plot_ly(d, x=~FSS, color=~trt, type='histogram', opacity=0.6,
            colors=c("No Treatment"="#E74C3C","Full Combo"="#1ABC9C")) %>%
      layout(barmode='overlay', xaxis=list(title="FSS"), title="VP: FSS Distribution (Week 52)")
  })

  output$vpop_scatter <- renderPlotly({
    req(vp_results())
    d <- vp_results()
    plot_ly(d, x=~VO2max, y=~FSS, color=~trt, type='scatter', mode='markers',
            colors=c("No Treatment"="#E74C3C","Full Combo"="#1ABC9C"),
            marker=list(size=5, opacity=0.7)) %>%
      layout(xaxis=list(title="VO2max (%)"), yaxis=list(title="FSS"),
             title="VP: VO2max vs FSS")
  })

  output$vpop_response <- renderPlotly({
    req(vp_results())
    d <- vp_results() %>%
      group_by(trt) %>%
      summarise(
        FSS_resp   = mean(FSS <= 4)*100,
        VO2_resp   = mean(VO2max >= 70)*100,
        MoCA_resp  = mean(MoCA >= 26)*100,
        .groups="drop"
      ) %>%
      pivot_longer(-trt, names_to="Outcome", values_to="ResponseRate")
    plot_ly(d, x=~Outcome, y=~ResponseRate, color=~trt, type='bar',
            colors=c("No Treatment"="#E74C3C","Full Combo"="#1ABC9C")) %>%
      layout(yaxis=list(title="Response Rate (%)"), title="Response Rates")
  })

  # Biomarker tab
  output$nfl_plot <- renderPlotly({
    d <- sim_data()
    plot_ly(d, x=~time, y=~NfL_pg, type='scatter', mode='lines',
            line=list(color="#8E44AD", width=2)) %>%
      add_lines(y=rep(10, nrow(d)), name="Normal <10 pg/mL",
                line=list(dash="dash", color="green")) %>%
      layout(xaxis=list(title="Day"), yaxis=list(title="NfL (pg/mL)"),
             title="Neurofilament Light Chain")
  })

  output$crp_plot <- renderPlotly({
    d <- sim_data()
    plot_ly(d, x=~time, y=~CRP_proxy, type='scatter', mode='lines',
            line=list(color="#E74C3C", width=2)) %>%
      add_lines(y=rep(0.5, nrow(d)), name="Normal <0.5 (norm)",
                line=list(dash="dash", color="green")) %>%
      layout(xaxis=list(title="Day"), yaxis=list(title="CRP (normalized)"),
             title="CRP (Proxy from IL-6)")
  })

  output$ddimer_plot <- renderPlotly({
    d <- sim_data()
    plot_ly(d, x=~time, y=~Ddimer_out, type='scatter', mode='lines',
            line=list(color="#C0392B", width=2)) %>%
      layout(xaxis=list(title="Day"), yaxis=list(title="D-Dimer (normalized)"),
             title="D-Dimer Trajectory")
  })

  output$biomarker_ref_table <- renderDT({
    tribble(
      ~Biomarker, ~`Normal Range`, ~`Long COVID Typical`, ~Significance,
      "NfL (plasma)", "<10 pg/mL",   "15-50 pg/mL",    "Neuroaxonal injury",
      "GFAP (plasma)","<100 pg/mL",  "150-400 pg/mL",  "Astrocyte damage",
      "IL-6",         "<2 pg/mL",    "5-20 pg/mL",     "Systemic inflammation",
      "D-Dimer",      "<0.5 µg/mL",  "0.5-2.0 µg/mL",  "Fibrin microclots",
      "TNF-α",        "<8.1 pg/mL",  "10-30 pg/mL",    "Inflammatory cytokine",
      "IFN-β",        "~undetectable","elevated",       "Type I IFN dysregulation",
      "Anti-ACE2 Ab", "negative",    "variable +",     "Autonomic autoimmunity"
    ) %>% datatable(options=list(dom="t"), rownames=FALSE)
  })

}

shinyApp(ui, server)
