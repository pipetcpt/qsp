## ============================================================
## Shiny App: Primary Open-Angle Glaucoma (POAG) QSP Dashboard
## Tabs: (1) Patient Profile | (2) Drug PK | (3) IOP Dynamics
##       (4) Clinical Endpoints | (5) Scenario Comparison
##       (6) Biomarkers & Neuroprotection | (7) Sensitivity Analysis
## ============================================================

library(shiny)
library(shinydashboard)
library(mrgsolve)
library(dplyr)
library(ggplot2)
library(plotly)
library(DT)
library(tidyr)

# ---- Inline mrgsolve model code (abbreviated) ----
model_code <- '
$PARAM
ke_PGA=0.693 ke_BB=0.231 ke_CAI=0.139 ke_A2A=0.462 ke_ROCK=0.347
F_prod_base=2.5 F_uv_base=0.40 C_tm_base=0.25 P_ep=8.0
tau_IOP=0.1
Emax_PGA=1.0 EC50_PGA=2.0 hill_PGA=1.5
Emax_BB=0.30 EC50_BB=150.0 hill_BB=1.2
Emax_CAI=0.25 EC50_CAI=500.0 hill_CAI=1.0
Emax_A2A=0.25 EC50_A2A=1.5 hill_A2A=1.2 Emax_A2A_uv=0.10
Emax_ROCK=0.35 EC50_ROCK=30.0 hill_ROCK=1.5 Emax_ROCK_ep=3.0
k_ECM_prod=0.001 k_ECM_clear=0.0005 k_TM_IOP=0.01 k_C_tm_decline=2.0e-4 C_tm_min=0.10
BDNF0=50.0 k_BDNF_prod=0.05 k_BDNF_deg=0.001 k_BDNF_IOP=0.002 IOP_thresh_BDNF=18.0 k_A2A_BDNF=0.20
RGC0=1.2 k_RGC_base=5.7e-6 EC50_Casp3=0.3 hill_Casp3=2.0
k_Casp3_act=0.005 k_Casp3_inact=0.02 IOP_Casp3_thresh=15.0 hill_IOP_Casp3=2.5
RNFL0=100.0 k_RNFL_RGC=1.0 k_VF_RGC=25.0 n_VF=2.5 VF_MD_floor=-30.0
diurnal_amp=0.10

$CMT C_PGA C_BB C_CAI C_A2A C_ROCK F_aq F_uv C_tm IOP ECM_TM BDNF Casp3 RGC RNFL VF_MD

$MAIN
F_aq_0=F_prod_base; F_uv_0=F_uv_base; C_tm_0=C_tm_base;
IOP_0=(F_prod_base-F_uv_base)/C_tm_base+P_ep;
ECM_TM_0=0.0; BDNF_0=BDNF0; Casp3_0=0.01; RGC_0=RGC0; RNFL_0=RNFL0; VF_MD_0=0.0;

$ODE
dxdt_C_PGA  = -ke_PGA  * C_PGA;
dxdt_C_BB   = -ke_BB   * C_BB;
dxdt_C_CAI  = -ke_CAI  * C_CAI;
dxdt_C_A2A  = -ke_A2A  * C_A2A;
dxdt_C_ROCK = -ke_ROCK * C_ROCK;

double E_PGA  = Emax_PGA  * pow(C_PGA,  hill_PGA)  / (pow(EC50_PGA,  hill_PGA)  + pow(C_PGA,  hill_PGA)  + 1e-12);
double E_BB   = Emax_BB   * pow(C_BB,   hill_BB)   / (pow(EC50_BB,   hill_BB)   + pow(C_BB,   hill_BB)   + 1e-12);
double E_CAI  = Emax_CAI  * pow(C_CAI,  hill_CAI)  / (pow(EC50_CAI,  hill_CAI)  + pow(C_CAI,  hill_CAI)  + 1e-12);
double E_A2A  = Emax_A2A  * pow(C_A2A,  hill_A2A)  / (pow(EC50_A2A,  hill_A2A)  + pow(C_A2A,  hill_A2A)  + 1e-12);
double E_A2A_uv = Emax_A2A_uv * pow(C_A2A, hill_A2A) / (pow(EC50_A2A, hill_A2A) + pow(C_A2A, hill_A2A) + 1e-12);
double E_ROCK = Emax_ROCK * pow(C_ROCK, hill_ROCK) / (pow(EC50_ROCK, hill_ROCK) + pow(C_ROCK, hill_ROCK) + 1e-12);
double E_ROCK_ep = Emax_ROCK_ep * pow(C_ROCK, hill_ROCK) / (pow(EC50_ROCK, hill_ROCK) + pow(C_ROCK, hill_ROCK) + 1e-12);
double diurnal_factor = 1.0 + diurnal_amp * cos(2*3.14159*(fmod(SOLVERTIME,24.0)-9.0)/24.0);

double F_aq_target = F_prod_base * diurnal_factor * (1.0 - E_BB - E_CAI - E_A2A);
double F_uv_target = F_uv_base * (1.0 + E_PGA + E_A2A_uv);
double C_tm_disease = C_tm_base * exp(-k_C_tm_decline * ECM_TM);
if(C_tm_disease < C_tm_min) C_tm_disease = C_tm_min;
double C_tm_target = C_tm_disease * (1.0 + E_ROCK);
double P_ep_eff = P_ep - E_ROCK_ep;

dxdt_F_aq = 5.0 * (F_aq_target - F_aq);
dxdt_F_uv = 2.0 * (F_uv_target - F_uv);
dxdt_C_tm = 0.5 * (C_tm_target - C_tm);

double IOP_eq = (F_aq - F_uv) / C_tm + P_ep_eff;
if(IOP_eq < 5.0) IOP_eq = 5.0;
dxdt_IOP = (1.0/tau_IOP) * (IOP_eq - IOP);

double IOP_stress = (IOP > 18.0) ? (IOP - 18.0) : 0.0;
dxdt_ECM_TM = k_ECM_prod + k_TM_IOP * IOP_stress - k_ECM_clear * ECM_TM - 0.5 * E_ROCK * ECM_TM;
if(ECM_TM < 0.0) ECM_TM = 0.0;

double IOP_above = (IOP > IOP_thresh_BDNF) ? (IOP - IOP_thresh_BDNF) : 0.0;
double A2A_BDNF_effect = 1.0 + k_A2A_BDNF * E_A2A;
dxdt_BDNF = k_BDNF_prod * A2A_BDNF_effect - (k_BDNF_deg + k_BDNF_IOP * IOP_above) * BDNF;

double IOP_Casp_stress = pow(IOP / IOP_Casp3_thresh, hill_IOP_Casp3);
double BDNF_surv = BDNF / (EC50_Casp3 * 100 + BDNF);
dxdt_Casp3 = k_Casp3_act * IOP_Casp_stress * (1.0 - BDNF_surv) - k_Casp3_inact * Casp3;
if(Casp3 < 0.0) Casp3 = 0.0; if(Casp3 > 1.0) Casp3 = 1.0;

double Casp3_eff = pow(Casp3, hill_Casp3) / (pow(EC50_Casp3, hill_Casp3) + pow(Casp3, hill_Casp3) + 1e-12);
double k_RGC_loss = k_RGC_base + k_Casp3_act * 2.0 * Casp3_eff;
dxdt_RGC = -k_RGC_loss * RGC;
if(RGC < 0.0) RGC = 0.0;

dxdt_RNFL = k_RNFL_RGC * (RGC - RNFL / k_RNFL_RGC) * 0.001;
double RGC_frac = RGC / RGC0; if(RGC_frac < 0.01) RGC_frac = 0.01;
double VF_MD_eq = -k_VF_RGC * pow(1.0 - RGC_frac, n_VF);
if(VF_MD_eq < VF_MD_floor) VF_MD_eq = VF_MD_floor;
dxdt_VF_MD = 0.05 * (VF_MD_eq - VF_MD);

$CAPTURE IOP F_aq F_uv C_tm ECM_TM BDNF Casp3 RGC RNFL VF_MD
C_PGA C_BB C_CAI C_A2A C_ROCK
E_PGA E_BB E_CAI E_A2A E_ROCK
'

mod <- mread_cache("poag_shiny", tempdir(), model_code)

# ---- Color palette ----
COLS8 <- c("#D32F2F","#1976D2","#388E3C","#F57C00","#7B1FA2","#0097A7","#5D4037","#455A64")

# ---- Helper: build dosing events ----
build_events <- function(use_PGA, use_BB, use_CAI, use_A2A, use_ROCK, sim_yrs, iop_init) {
  end_h <- sim_yrs * 365 * 24
  n_days <- floor(end_h / 24)
  ev_list <- list()
  if(use_PGA) {
    times <- seq(21, n_days*24+21, by=24)
    ev_list[[length(ev_list)+1]] <- ev(cmt="C_PGA", amt=8,  time=times, rate=-2)
  }
  if(use_BB) {
    times <- c(sapply(0:(n_days-1), function(d) c(d*24+8, d*24+20)))
    ev_list[[length(ev_list)+1]] <- ev(cmt="C_BB",  amt=200, time=times, rate=-2)
  }
  if(use_CAI) {
    times <- c(sapply(0:(n_days-1), function(d) c(d*24+8, d*24+14, d*24+20)))
    ev_list[[length(ev_list)+1]] <- ev(cmt="C_CAI", amt=1500, time=times, rate=-2)
  }
  if(use_A2A) {
    times <- c(sapply(0:(n_days-1), function(d) c(d*24+8, d*24+20)))
    ev_list[[length(ev_list)+1]] <- ev(cmt="C_A2A", amt=0.5, time=times, rate=-2)
  }
  if(use_ROCK) {
    times <- seq(21, n_days*24+21, by=24)
    ev_list[[length(ev_list)+1]] <- ev(cmt="C_ROCK", amt=50, time=times, rate=-2)
  }
  if(length(ev_list)==0) return(NULL)
  Reduce(c, ev_list)
}

# ---- Run single simulation ----
run_sim <- function(iop_init=24, rnfl_init=88, vf_init=-3.5, bdnf_init=35,
                    use_PGA=FALSE, use_BB=FALSE, use_CAI=FALSE, use_A2A=FALSE, use_ROCK=FALSE,
                    sim_yrs=10) {
  C_tm_init <- 0.25
  init_vals <- c(C_PGA=0, C_BB=0, C_CAI=0, C_A2A=0, C_ROCK=0,
                 F_aq=2.5, F_uv=0.40, C_tm=C_tm_init, IOP=iop_init,
                 ECM_TM=0.15, BDNF=bdnf_init, Casp3=0.05,
                 RGC=1.1, RNFL=rnfl_init, VF_MD=vf_init)
  events <- build_events(use_PGA, use_BB, use_CAI, use_A2A, use_ROCK, sim_yrs, iop_init)
  end_h <- sim_yrs * 365 * 24
  if(is.null(events)) {
    out <- mod %>% init(init_vals) %>% mrgsim(end=end_h, delta=24)
  } else {
    out <- mod %>% init(init_vals) %>% mrgsim(events=events, end=end_h, delta=24)
  }
  as.data.frame(out) %>% mutate(time_yr = time / (365*24))
}

# ============================================================
# UI
# ============================================================
ui <- dashboardPage(
  skin = "blue",
  dashboardHeader(title = "POAG QSP Dashboard"),

  dashboardSidebar(
    sidebarMenu(
      menuItem("Patient Profile",     tabName = "tab_patient",   icon = icon("user-md")),
      menuItem("Drug PK",             tabName = "tab_pk",        icon = icon("pills")),
      menuItem("IOP Dynamics",        tabName = "tab_iop",       icon = icon("eye")),
      menuItem("Clinical Endpoints",  tabName = "tab_endpoints", icon = icon("chart-line")),
      menuItem("Scenario Comparison", tabName = "tab_scenarios", icon = icon("balance-scale")),
      menuItem("Biomarkers & Neuroprot.", tabName = "tab_bio",   icon = icon("dna")),
      menuItem("Sensitivity Analysis", tabName = "tab_sa",       icon = icon("sliders-h"))
    )
  ),

  dashboardBody(
    tags$head(tags$style(HTML(".content-wrapper{background-color:#f4f6f9;}"))),
    tabItems(

      # ------------------------------------------------------------------
      # TAB 1: Patient Profile
      # ------------------------------------------------------------------
      tabItem(tabName = "tab_patient",
        fluidRow(
          box(title = "Patient Characteristics", width = 4, solidHeader = TRUE,
              status = "primary",
              numericInput("iop_init",  "Baseline IOP (mmHg)", value = 24, min = 10, max = 40),
              numericInput("rnfl_init", "Baseline RNFL (μm)",  value = 88, min = 50, max = 120),
              numericInput("vf_init",   "Baseline VF-MD (dB)", value = -3.5, min = -30, max = 0),
              numericInput("bdnf_init", "Baseline BDNF (pg/mL)", value = 35, min = 5, max = 80),
              numericInput("sim_yrs",   "Simulation Duration (years)", value = 10, min = 1, max = 20),
              hr(),
              h5("Glaucoma Stage (Hodapp-Parrish-Anderson):"),
              verbatimTextOutput("stage_out")
          ),
          box(title = "Diagnosis Summary", width = 8, solidHeader = TRUE, status = "info",
              fluidRow(
                valueBoxOutput("vb_iop",  width = 4),
                valueBoxOutput("vb_rnfl", width = 4),
                valueBoxOutput("vb_vf",   width = 4)
              ),
              hr(),
              h4("Primary Open-Angle Glaucoma — Pathophysiology Overview"),
              tags$ul(
                tags$li("Aqueous production: Ciliary body NPE → Na⁺/K⁺-ATPase + CA-II → HCO₃⁻"),
                tags$li("Conventional outflow: Trabecular meshwork → Schlemm's canal → Episcleral veins"),
                tags$li("Uveoscleral outflow: Ciliary body stroma → Suprachoroidal space"),
                tags$li("Goldmann equation: IOP = (F_prod − F_uv)/C_tm + P_ep"),
                tags$li("Optic nerve damage: TPG → LC deformation → Axon transport block → BDNF ↓ → RGC apoptosis"),
                tags$li("Treatment targets: ↓ F_prod (BB, CAI, A2A), ↑ F_uv (PGA, A2A), ↑ C_tm (ROCK-I, pilocarpine)")
              )
          )
        )
      ),

      # ------------------------------------------------------------------
      # TAB 2: Drug PK
      # ------------------------------------------------------------------
      tabItem(tabName = "tab_pk",
        fluidRow(
          box(title = "Drug Selection", width = 3, solidHeader = TRUE, status = "primary",
              checkboxInput("use_PGA",  "Latanoprost 0.005% QD (PGA)",      value = TRUE),
              checkboxInput("use_BB",   "Timolol 0.5% BID (Beta-Blocker)",  value = FALSE),
              checkboxInput("use_CAI",  "Dorzolamide 2% TID (CAI)",         value = FALSE),
              checkboxInput("use_A2A",  "Brimonidine 0.2% BID (A2A)",       value = FALSE),
              checkboxInput("use_ROCK", "Netarsudil 0.02% QD (ROCK-I)",     value = FALSE),
              hr(),
              sliderInput("pk_days", "PK Observation Window (Days)", 1, 14, 7),
              actionButton("run_pk", "Run PK Simulation", class="btn-primary btn-block")
          ),
          box(title = "Aqueous Drug Concentrations", width = 9, solidHeader = TRUE,
              plotlyOutput("pk_plot", height = "420px"))
        ),
        fluidRow(
          box(title = "Drug PK Parameters", width = 12,
              DTOutput("pk_table"))
        )
      ),

      # ------------------------------------------------------------------
      # TAB 3: IOP Dynamics
      # ------------------------------------------------------------------
      tabItem(tabName = "tab_iop",
        fluidRow(
          box(title = "Treatment Configuration", width = 3, solidHeader = TRUE, status = "warning",
              checkboxInput("iop_PGA",  "PGA (Latanoprost)",   value = TRUE),
              checkboxInput("iop_BB",   "BB (Timolol)",        value = FALSE),
              checkboxInput("iop_CAI",  "CAI (Dorzolamide)",   value = FALSE),
              checkboxInput("iop_A2A",  "A2A (Brimonidine)",   value = FALSE),
              checkboxInput("iop_ROCK", "ROCK-I (Netarsudil)", value = FALSE),
              hr(),
              sliderInput("iop_horizon", "View Horizon (Days)", 1, 365, 90),
              actionButton("run_iop", "Simulate IOP", class="btn-warning btn-block")
          ),
          box(title = "IOP Over Time", width = 9, solidHeader = TRUE,
              plotlyOutput("iop_plot", height = "380px"))
        ),
        fluidRow(
          box(title = "Aqueous Flow Dynamics", width = 6,
              plotlyOutput("flow_plot", height = "300px")),
          box(title = "TM Outflow Facility (C_tm)", width = 6,
              plotlyOutput("ctm_plot", height = "300px"))
        )
      ),

      # ------------------------------------------------------------------
      # TAB 4: Clinical Endpoints
      # ------------------------------------------------------------------
      tabItem(tabName = "tab_endpoints",
        fluidRow(
          box(title = "Simulation Controls", width = 3, solidHeader = TRUE, status = "success",
              checkboxInput("ep_PGA",  "PGA",    value = TRUE),
              checkboxInput("ep_BB",   "BB",     value = FALSE),
              checkboxInput("ep_CAI",  "CAI",    value = FALSE),
              checkboxInput("ep_A2A",  "A2A",    value = FALSE),
              checkboxInput("ep_ROCK", "ROCK-I", value = FALSE),
              actionButton("run_ep", "Run Simulation", class="btn-success btn-block"),
              hr(),
              h5("Year 10 Summary:"),
              verbatimTextOutput("ep_summary")
          ),
          box(title = "VF-MD Progression", width = 9, solidHeader = TRUE,
              plotlyOutput("vf_plot", height = "380px"))
        ),
        fluidRow(
          box(title = "RNFL Thickness", width = 6,
              plotlyOutput("rnfl_plot", height = "280px")),
          box(title = "RGC Count", width = 6,
              plotlyOutput("rgc_plot",  height = "280px"))
        )
      ),

      # ------------------------------------------------------------------
      # TAB 5: Scenario Comparison
      # ------------------------------------------------------------------
      tabItem(tabName = "tab_scenarios",
        fluidRow(
          box(title = "Compare Treatment Scenarios", width = 12, solidHeader = TRUE,
              status = "info",
              actionButton("run_compare", "Run All Scenarios", class="btn-info btn-block"),
              br(), br(),
              p("Scenarios: Untreated | PGA | BB | CAI | A2A | ROCK-I | PGA+BB | Triple (PGA+BB+CAI)")
          )
        ),
        fluidRow(
          box(title = "IOP Comparison", width = 6,
              plotlyOutput("cmp_iop",  height = "320px")),
          box(title = "VF-MD Comparison", width = 6,
              plotlyOutput("cmp_vf",   height = "320px"))
        ),
        fluidRow(
          box(title = "RNFL Comparison", width = 6,
              plotlyOutput("cmp_rnfl", height = "280px")),
          box(title = "10-Year Results Table", width = 6,
              DTOutput("cmp_table"))
        )
      ),

      # ------------------------------------------------------------------
      # TAB 6: Biomarkers & Neuroprotection
      # ------------------------------------------------------------------
      tabItem(tabName = "tab_bio",
        fluidRow(
          box(title = "Biomarker Analysis", width = 3, solidHeader = TRUE, status = "purple",
              checkboxInput("bio_PGA",  "PGA",    value = FALSE),
              checkboxInput("bio_BB",   "BB",     value = FALSE),
              checkboxInput("bio_A2A",  "A2A (+ neuroprot.)", value = TRUE),
              checkboxInput("bio_ROCK", "ROCK-I", value = FALSE),
              actionButton("run_bio", "Run Biomarker Sim", class="btn-block",
                           style="background:#7B1FA2;color:white;"),
              hr(),
              p("A2A agonism (brimonidine) provides neuroprotection via",
                "BDNF upregulation, independent of IOP lowering.")
          ),
          box(title = "BDNF in Optic Nerve Head", width = 9, solidHeader = TRUE,
              plotlyOutput("bdnf_plot", height = "360px"))
        ),
        fluidRow(
          box(title = "Caspase-3 Apoptotic Activity", width = 6,
              plotlyOutput("casp_plot", height = "280px")),
          box(title = "ECM Accumulation in TM", width = 6,
              plotlyOutput("ecm_plot",  height = "280px"))
        )
      ),

      # ------------------------------------------------------------------
      # TAB 7: Sensitivity Analysis
      # ------------------------------------------------------------------
      tabItem(tabName = "tab_sa",
        fluidRow(
          box(title = "Parameter Sensitivity", width = 4, solidHeader = TRUE,
              status = "danger",
              sliderInput("sa_iop",    "Initial IOP (mmHg)",    16, 36, 24),
              sliderInput("sa_ctm",    "C_tm Baseline (×normal)", 0.5, 1.0, 0.83),
              sliderInput("sa_bdnf",   "Initial BDNF (% normal)", 40, 100, 70),
              sliderInput("sa_Emax_PGA", "PGA Emax (F_uv increase)", 0.5, 2.0, 1.0),
              actionButton("run_sa", "Run Sensitivity", class="btn-danger btn-block")
          ),
          box(title = "Tornado: Year-10 VF-MD Sensitivity to Parameters", width = 8,
              plotlyOutput("tornado_plot", height = "420px"))
        ),
        fluidRow(
          box(title = "IOP Reduction vs VF Preservation (10-yr)", width = 12,
              plotlyOutput("iop_vf_scatter", height = "350px"))
        )
      )
    )
  )
)

# ============================================================
# SERVER
# ============================================================
server <- function(input, output, session) {

  # ---- Reactive: run single scenario with current patient profile ----
  sim_data <- eventReactive(
    list(input$run_iop, input$run_ep, input$run_bio, input$run_pk), {
      run_sim(
        iop_init  = input$iop_init,
        rnfl_init = input$rnfl_init,
        vf_init   = input$vf_init,
        bdnf_init = input$bdnf_init,
        use_PGA   = isolate(input$iop_PGA  | input$ep_PGA  | input$bio_PGA),
        use_BB    = isolate(input$iop_BB   | input$ep_BB   | input$bio_BB),
        use_CAI   = isolate(input$iop_CAI  | input$ep_CAI),
        use_A2A   = isolate(input$iop_A2A  | input$ep_A2A  | input$bio_A2A),
        use_ROCK  = isolate(input$iop_ROCK | input$ep_ROCK | input$bio_ROCK),
        sim_yrs   = input$sim_yrs
      )
    }, ignoreNULL = FALSE)

  # ---- TAB 1 outputs ----
  output$stage_out <- renderText({
    vf <- input$vf_init
    if(vf > -6)       "Early Glaucoma (MD > −6 dB)"
    else if(vf > -12) "Moderate Glaucoma (−6 to −12 dB)"
    else              "Advanced Glaucoma (MD < −12 dB)"
  })
  output$vb_iop  <- renderValueBox(valueBox(paste0(input$iop_init, " mmHg"), "Baseline IOP", icon=icon("eye"), color="red"))
  output$vb_rnfl <- renderValueBox(valueBox(paste0(input$rnfl_init, " μm"),  "RNFL Thickness", icon=icon("layer-group"), color="blue"))
  output$vb_vf   <- renderValueBox(valueBox(paste0(input$vf_init, " dB"),   "VF Mean Deviation", icon=icon("chart-bar"), color="orange"))

  # ---- TAB 2: PK ----
  pk_sim <- eventReactive(input$run_pk, {
    end_h <- input$pk_days * 24
    events <- build_events(input$use_PGA, input$use_BB, input$use_CAI, input$use_A2A, input$use_ROCK, sim_yrs=0.5, iop_init=24)
    init_vals <- c(C_PGA=0,C_BB=0,C_CAI=0,C_A2A=0,C_ROCK=0,
                   F_aq=2.5,F_uv=0.4,C_tm=0.25,IOP=24,ECM_TM=0.15,BDNF=35,Casp3=0.05,RGC=1.1,RNFL=88,VF_MD=-3.5)
    if(is.null(events)) {
      mod %>% init(init_vals) %>% mrgsim(end=end_h, delta=0.5)
    } else {
      mod %>% init(init_vals) %>% mrgsim(events=events, end=end_h, delta=0.5)
    }
  }, ignoreNULL = FALSE)

  output$pk_plot <- renderPlotly({
    df <- as.data.frame(pk_sim())
    df_long <- df %>% select(time, C_PGA, C_BB, C_CAI, C_A2A, C_ROCK) %>%
      pivot_longer(-time, names_to = "Drug", values_to = "Conc")
    p <- ggplot(df_long, aes(x=time, y=Conc+1e-6, color=Drug)) +
      geom_line(size=0.9) + scale_y_log10() +
      scale_color_manual(values=COLS8[1:5],
        labels=c("PGA (ng/mL)","BB (ng/mL)","CAI (ng/mL)","A2A (ng/mL)","ROCK-I (ng/mL)")) +
      labs(x="Time (h)", y="Aqueous Concentration (ng/mL, log)", color="Drug Class") +
      theme_bw()
    ggplotly(p)
  })

  output$pk_table <- renderDT({
    data.frame(
      Drug         = c("Latanoprost acid (PGA)","Timolol (BB)","Dorzolamide (CAI)","Brimonidine (A2A)","Netarsudil (ROCK-I)"),
      Dosing       = c("50 μg QD PM","0.5% BID","2% TID","0.2% BID","0.02% QD"),
      ke_per_h     = c(0.693,0.231,0.139,0.462,0.347),
      t_half_h     = c(1.0,3.0,5.0,1.5,2.0),
      EC50         = c("2 ng/mL","150 ng/mL","500 ng/mL","1.5 ng/mL","30 ng/mL"),
      Mechanism    = c("FP-R → ↑F_uv","β₂-AR blockade → ↓F_prod","CA-II/IV inhib → ↓F_prod","α₂-AR → ↓F_prod + neuroprot.","ROCK1/2 inhib → ↑C_tm")
    )
  }, options=list(pageLength=5, dom="t"))

  # ---- TAB 3: IOP ----
  iop_data <- eventReactive(input$run_iop, {
    run_sim(iop_init=input$iop_init, rnfl_init=input$rnfl_init, vf_init=input$vf_init,
            bdnf_init=input$bdnf_init,
            use_PGA=input$iop_PGA, use_BB=input$iop_BB, use_CAI=input$iop_CAI,
            use_A2A=input$iop_A2A, use_ROCK=input$iop_ROCK, sim_yrs=input$sim_yrs)
  }, ignoreNULL=FALSE)

  output$iop_plot <- renderPlotly({
    df <- iop_data() %>% filter(time <= input$iop_horizon * 24)
    p <- ggplot(df, aes(x=time/24, y=IOP)) + geom_line(color="#1976D2", size=1) +
      geom_hline(yintercept=21, linetype="dashed", color="red") +
      geom_hline(yintercept=18, linetype="dotted", color="navy") +
      labs(x="Time (Days)", y="IOP (mmHg)", title="IOP Response to Treatment") + theme_bw()
    ggplotly(p)
  })
  output$flow_plot <- renderPlotly({
    df <- iop_data() %>% filter(time <= input$iop_horizon * 24) %>%
      select(time, F_aq, F_uv) %>% pivot_longer(-time)
    p <- ggplot(df, aes(x=time/24, y=value, color=name)) + geom_line(size=0.8) +
      scale_color_manual(values=c("#D32F2F","#388E3C"), labels=c("F_aq (prod)","F_uv (uveosc)")) +
      labs(x="Days", y="Flow (μL/min)", color="") + theme_bw()
    ggplotly(p)
  })
  output$ctm_plot <- renderPlotly({
    df <- iop_data()
    p <- ggplot(df, aes(x=time_yr, y=C_tm)) + geom_line(color="#F57C00", size=1) +
      geom_hline(yintercept=0.30, linetype="dashed", color="grey50") +
      labs(x="Years", y="C_tm (μL/min/mmHg)", title="TM Outflow Facility") + theme_bw()
    ggplotly(p)
  })

  # ---- TAB 4: Endpoints ----
  ep_data <- eventReactive(input$run_ep, {
    run_sim(iop_init=input$iop_init, rnfl_init=input$rnfl_init, vf_init=input$vf_init,
            bdnf_init=input$bdnf_init,
            use_PGA=input$ep_PGA, use_BB=input$ep_BB, use_CAI=input$ep_CAI,
            use_A2A=input$ep_A2A, use_ROCK=input$ep_ROCK, sim_yrs=input$sim_yrs)
  }, ignoreNULL=FALSE)

  output$vf_plot <- renderPlotly({
    df <- ep_data()
    p <- ggplot(df, aes(x=time_yr, y=VF_MD)) + geom_line(color="#388E3C", size=1.2) +
      geom_hline(yintercept=-6, linetype="dashed", color="#FF8C00") +
      geom_hline(yintercept=-12, linetype="dashed", color="#D32F2F") +
      labs(x="Years", y="VF-MD (dB)", title="Visual Field Progression") + theme_bw()
    ggplotly(p)
  })
  output$rnfl_plot <- renderPlotly({
    df <- ep_data()
    p <- ggplot(df, aes(x=time_yr, y=RNFL)) + geom_line(color="#7B1FA2", size=1) +
      geom_hline(yintercept=80, linetype="dashed", color="red") +
      labs(x="Years", y="RNFL (μm)") + theme_bw()
    ggplotly(p)
  })
  output$rgc_plot <- renderPlotly({
    df <- ep_data()
    p <- ggplot(df, aes(x=time_yr, y=RGC)) + geom_line(color="#D32F2F", size=1) +
      labs(x="Years", y="RGC (millions)") + theme_bw()
    ggplotly(p)
  })
  output$ep_summary <- renderText({
    df <- ep_data() %>% filter(time_yr >= (input$sim_yrs - 0.1))
    paste0(
      "IOP:     ", round(mean(df$IOP),1), " mmHg\n",
      "VF-MD:   ", round(mean(df$VF_MD),2), " dB\n",
      "RNFL:    ", round(mean(df$RNFL),1), " μm\n",
      "RGC:     ", round(mean(df$RGC),3), " M"
    )
  })

  # ---- TAB 5: Scenario Comparison ----
  cmp_data <- eventReactive(input$run_compare, {
    scens <- list(
      list(nm="Untreated",   F=F, B=F, C=F, A=F, R=F),
      list(nm="PGA",         F=T, B=F, C=F, A=F, R=F),
      list(nm="BB",          F=F, B=T, C=F, A=F, R=F),
      list(nm="CAI",         F=F, B=F, C=T, A=F, R=F),
      list(nm="A2A",         F=F, B=F, C=F, A=T, R=F),
      list(nm="ROCK-I",      F=F, B=F, C=F, A=F, R=T),
      list(nm="PGA+BB",      F=T, B=T, C=F, A=F, R=F),
      list(nm="Triple",      F=T, B=T, C=T, A=F, R=F)
    )
    df_all <- lapply(scens, function(s) {
      d <- run_sim(iop_init=input$iop_init, rnfl_init=input$rnfl_init,
                   vf_init=input$vf_init, bdnf_init=input$bdnf_init,
                   use_PGA=s$F, use_BB=s$B, use_CAI=s$C, use_A2A=s$A, use_ROCK=s$R,
                   sim_yrs=input$sim_yrs)
      d$scenario <- s$nm; d
    })
    do.call(rbind, df_all) %>% mutate(scenario = factor(scenario, levels=sapply(scens, `[[`, "nm")))
  })

  output$cmp_iop <- renderPlotly({
    df <- cmp_data()
    p <- ggplot(df, aes(x=time_yr, y=IOP, color=scenario)) +
      geom_line(size=0.8) + scale_color_manual(values=COLS8) +
      geom_hline(yintercept=18, linetype="dashed", color="navy", alpha=0.5) +
      labs(x="Years", y="IOP (mmHg)", color="") + theme_bw(base_size=10)
    ggplotly(p) %>% layout(legend=list(orientation="h"))
  })
  output$cmp_vf <- renderPlotly({
    df <- cmp_data()
    p <- ggplot(df, aes(x=time_yr, y=VF_MD, color=scenario)) +
      geom_line(size=0.8) + scale_color_manual(values=COLS8) +
      labs(x="Years", y="VF-MD (dB)", color="") + theme_bw(base_size=10)
    ggplotly(p) %>% layout(legend=list(orientation="h"))
  })
  output$cmp_rnfl <- renderPlotly({
    df <- cmp_data()
    p <- ggplot(df, aes(x=time_yr, y=RNFL, color=scenario)) +
      geom_line(size=0.8) + scale_color_manual(values=COLS8) +
      labs(x="Years", y="RNFL (μm)", color="") + theme_bw(base_size=10)
    ggplotly(p) %>% layout(legend=list(orientation="h"))
  })
  output$cmp_table <- renderDT({
    df <- cmp_data() %>% filter(time_yr >= (input$sim_yrs - 0.1)) %>%
      group_by(scenario) %>%
      summarise(IOP=round(mean(IOP),1), VF_MD=round(mean(VF_MD),2),
                RNFL=round(mean(RNFL),1), RGC_M=round(mean(RGC),3), .groups="drop")
    datatable(df, options=list(pageLength=8, dom="t"))
  })

  # ---- TAB 6: Biomarkers ----
  bio_data <- eventReactive(input$run_bio, {
    run_sim(iop_init=input$iop_init, rnfl_init=input$rnfl_init, vf_init=input$vf_init,
            bdnf_init=input$bdnf_init,
            use_PGA=input$bio_PGA, use_BB=input$bio_BB, use_A2A=input$bio_A2A, use_ROCK=input$bio_ROCK,
            sim_yrs=input$sim_yrs)
  }, ignoreNULL=FALSE)

  output$bdnf_plot <- renderPlotly({
    df <- bio_data()
    p <- ggplot(df, aes(x=time_yr, y=BDNF)) + geom_line(color="#1565C0", size=1.2) +
      geom_hline(yintercept=50, linetype="dashed", color="grey50") +
      labs(x="Years", y="BDNF (pg/mL)", title="BDNF in Optic Nerve Head") + theme_bw()
    ggplotly(p)
  })
  output$casp_plot <- renderPlotly({
    df <- bio_data()
    p <- ggplot(df, aes(x=time_yr, y=Casp3)) + geom_line(color="#D32F2F", size=1) +
      labs(x="Years", y="Caspase-3 Index (0–1)") + theme_bw()
    ggplotly(p)
  })
  output$ecm_plot <- renderPlotly({
    df <- bio_data()
    p <- ggplot(df, aes(x=time_yr, y=ECM_TM)) + geom_line(color="#F57C00", size=1) +
      labs(x="Years", y="ECM Accumulation Index") + theme_bw()
    ggplotly(p)
  })

  # ---- TAB 7: Sensitivity ----
  sa_data <- eventReactive(input$run_sa, {
    params_vary <- list(
      list(par="iop_init",  lo=16, hi=36, base=input$sa_iop),
      list(par="bdnf_init", lo=15, hi=70, base=input$sa_bdnf * 0.7),
      list(par="sim_yrs",   lo=5,  hi=15, base=10)
    )
    baseline_run <- run_sim(iop_init=input$sa_iop, bdnf_init=input$sa_bdnf*0.5,
                            sim_yrs=10, use_PGA=TRUE)
    base_vf <- mean(baseline_run$VF_MD[baseline_run$time_yr >= 9.9], na.rm=TRUE)

    results <- lapply(params_vary, function(p) {
      run_lo <- if(p$par=="iop_init") {
        run_sim(iop_init=p$lo, bdnf_init=35, sim_yrs=10, use_PGA=TRUE)
      } else if(p$par=="bdnf_init") {
        run_sim(iop_init=input$sa_iop, bdnf_init=p$lo, sim_yrs=10, use_PGA=TRUE)
      } else {
        run_sim(iop_init=input$sa_iop, sim_yrs=p$lo, use_PGA=TRUE)
      }
      run_hi <- if(p$par=="iop_init") {
        run_sim(iop_init=p$hi, bdnf_init=35, sim_yrs=10, use_PGA=TRUE)
      } else if(p$par=="bdnf_init") {
        run_sim(iop_init=input$sa_iop, bdnf_init=p$hi, sim_yrs=10, use_PGA=TRUE)
      } else {
        run_sim(iop_init=input$sa_iop, sim_yrs=p$hi, use_PGA=TRUE)
      }
      vf_lo <- mean(run_lo$VF_MD[run_lo$time_yr >= min(run_lo$time_yr)*0.99+9.0*0.01], na.rm=TRUE)
      vf_hi <- mean(run_hi$VF_MD[run_hi$time_yr >= min(run_hi$time_yr)*0.99+9.0*0.01], na.rm=TRUE)
      data.frame(Parameter=p$par, Base_VF=round(base_vf,2),
                 VF_lo=round(vf_lo,2), VF_hi=round(vf_hi,2))
    })
    do.call(rbind, results)
  })

  output$tornado_plot <- renderPlotly({
    df <- sa_data()
    df$range <- abs(df$VF_hi - df$VF_lo)
    df <- df[order(df$range),]
    p <- ggplot(df, aes(y=reorder(Parameter, range))) +
      geom_segment(aes(x=VF_lo, xend=VF_hi, yend=Parameter), size=10, color="#1976D2") +
      geom_vline(aes(xintercept=Base_VF[1]), linetype="dashed", color="red") +
      labs(x="Year-10 VF-MD (dB)", y="Parameter", title="Tornado Plot — Sensitivity to Parameters") +
      theme_bw()
    ggplotly(p)
  })

  output$iop_vf_scatter <- renderPlotly({
    iop_range <- seq(12, 32, by=2)
    vf10_list <- sapply(iop_range, function(iop_val) {
      d <- run_sim(iop_init=iop_val, sim_yrs=10, use_PGA=(iop_val > 22))
      mean(d$VF_MD[d$time_yr >= 9.9], na.rm=TRUE)
    })
    df <- data.frame(IOP_init=iop_range, VF10=vf10_list)
    p <- ggplot(df, aes(x=IOP_init, y=VF10)) +
      geom_line(color="#7B1FA2", size=1.2) + geom_point(size=3, color="#7B1FA2") +
      labs(x="Initial IOP (mmHg)", y="Year-10 VF-MD (dB)",
           title="IOP Level vs Long-term Visual Field Outcome") + theme_bw()
    ggplotly(p)
  })
}

# ============================================================
# Launch
# ============================================================
shinyApp(ui = ui, server = server)
