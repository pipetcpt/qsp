## ============================================================
## EGPA QSP Shiny Dashboard
## Eosinophilic Granulomatosis with Polyangiitis
## ============================================================
## Tabs:
##  1. Patient Profile & Disease Background
##  2. Pharmacokinetics (PK)
##  3. Eosinophil & Cytokine Dynamics (Core PD)
##  4. Vasculitis & Organ Damage
##  5. Clinical Endpoints (BVAS, FEV1, LVEF, eGFR)
##  6. Scenario Comparison
##  7. Biomarker Trajectory & Remission
## ============================================================

library(shiny)
library(shinydashboard)
library(mrgsolve)
library(dplyr)
library(ggplot2)
library(plotly)
library(DT)

## ---- Inline model (same PD skeleton as egpa_mrgsolve_model.R) ----
egpa_code <- '
$PARAM
ka_mepo=0.0041, F_mepo=0.80, CL_mepo=0.022, Vc_mepo=4.6,
Vp_mepo=3.3,    Q_mepo=0.029, kon_mepo=0.0035, koff_mepo=0.001,
kdeg_IL5c=0.030, kdeg_cplx=0.005,
ka_benra=0.0033, F_benra=0.59, CL_benra=0.014, Vc_benra=4.1,
Vp_benra=3.7,    Q_benra=0.035, IC50_benra=0.12, Kmax_benra=0.90,
ka_pred=0.53,    F_pred=0.82,  CL_pred=0.32,  Vc_pred=9.5, Ke0_pred=0.12,
CL_cyclo=4.2,    Vc_cyclo=28.0, CLm_cyclo=2.8, Vm_cyclo=15.0, kel_4OH=0.38,
CL_ritu=0.013,   Vc_ritu=3.0,   Vp_ritu=4.0,   Q_ritu=0.030,
kon_ritu=0.055,  koff_ritu=0.003, kdeg_ritu=0.008,
kprod_Th2=0.015, kdeg_Th2=0.010, Th2_0=1.5,
kprod_IL5=0.080, kdeg_IL5=0.050, IL5_0=0.10,
kprod_Eos=40.0,  kdeg_Eos=0.010, kmig_Eos=0.008, Eos_B0=500,
kprod_EosT=0.005, kdeg_EosT=0.008, EosT_0=50.0,
kprod_IgE=0.020, kdeg_IgE=0.0012, IgE_0=180,
kprod_ANCA=0.005, kdeg_ANCA=0.006,
kact_Vasc=0.003, kres_Vasc=0.012,
kact_Card=0.002, kres_Card=0.005,
kact_Nerv=0.0015, kres_Nerv=0.004,
Emax_pred=0.85,  EC50_pred=120, n_pred=1.5,
Emax_cyclo=0.90, EC50_cyclo=2.0,
Emax_ritu=0.95,  EC50_ritu=15.0,
ANCA_pos=0

$INIT
MEPO_DEPOT=0, MEPO_C=0, MEPO_P=0, MEPO_CPX=0,
BENRA_DEPOT=0, BENRA_C=0, BENRA_P=0,
PRED_DEPOT=0, PRED_C=0, PRED_EFF=0,
CYCLO_C=0, CYCLO_M=0,
RITU_C=0, RITU_P=0, RITU_CPX=0,
TH2=1.5, IL5=0.25, EOS_B=5000, EOS_T=200, IGE=350,
ANCA=0, VASC=3.0, CARD=0.5, NERV=1.0

$ODE
double IL5_free = IL5;
dxdt_MEPO_DEPOT = -ka_mepo * MEPO_DEPOT;
dxdt_MEPO_C     = (ka_mepo*MEPO_DEPOT*F_mepo)/Vc_mepo
                  -(CL_mepo/Vc_mepo)*MEPO_C -(Q_mepo/Vc_mepo)*MEPO_C
                  +(Q_mepo/Vp_mepo)*MEPO_P
                  -kon_mepo*MEPO_C*IL5_free*Vc_mepo +koff_mepo*MEPO_CPX;
dxdt_MEPO_P     = (Q_mepo/Vc_mepo)*MEPO_C -(Q_mepo/Vp_mepo)*MEPO_P;
dxdt_MEPO_CPX   = kon_mepo*MEPO_C*IL5_free*Vc_mepo
                  -koff_mepo*MEPO_CPX -kdeg_cplx*MEPO_CPX;
dxdt_BENRA_DEPOT= -ka_benra*BENRA_DEPOT;
dxdt_BENRA_C    = (ka_benra*BENRA_DEPOT*F_benra)/Vc_benra
                  -(CL_benra/Vc_benra)*BENRA_C -(Q_benra/Vc_benra)*BENRA_C
                  +(Q_benra/Vp_benra)*BENRA_P;
dxdt_BENRA_P    = (Q_benra/Vc_benra)*BENRA_C -(Q_benra/Vp_benra)*BENRA_P;
dxdt_PRED_DEPOT = -ka_pred*PRED_DEPOT;
dxdt_PRED_C     = (ka_pred*PRED_DEPOT*F_pred)/Vc_pred -(CL_pred/Vc_pred)*PRED_C;
dxdt_PRED_EFF   = Ke0_pred*(PRED_C - PRED_EFF);
dxdt_CYCLO_C    = -(CL_cyclo+CLm_cyclo)/Vc_cyclo*CYCLO_C;
dxdt_CYCLO_M    = (CLm_cyclo/Vc_cyclo)*CYCLO_C -kel_4OH*CYCLO_M;
dxdt_RITU_C     = -(CL_ritu/Vc_ritu)*RITU_C -(Q_ritu/Vc_ritu)*RITU_C
                  +(Q_ritu/Vp_ritu)*RITU_P -kon_ritu*RITU_C +koff_ritu*RITU_CPX;
dxdt_RITU_P     = (Q_ritu/Vc_ritu)*RITU_C -(Q_ritu/Vp_ritu)*RITU_P;
dxdt_RITU_CPX   = kon_ritu*RITU_C -(koff_ritu+kdeg_ritu)*RITU_CPX;
double E_pred_Th2 = PRED_EFF>0 ? Emax_pred*pow(PRED_EFF,n_pred)/(pow(EC50_pred,n_pred)+pow(PRED_EFF,n_pred)) : 0;
double E_pred_Eos = PRED_EFF>0 ? Emax_pred*PRED_EFF/(EC50_pred*0.6+PRED_EFF) : 0;
double E_cyclo    = CYCLO_M>0  ? Emax_cyclo*CYCLO_M/(EC50_cyclo+CYCLO_M) : 0;
double E_ritu     = RITU_CPX>0 ? Emax_ritu*RITU_CPX/(EC50_ritu+RITU_CPX) : 0;
double frac_mepo  = (IL5+MEPO_CPX/Vc_mepo)>0 ? MEPO_CPX/(MEPO_CPX+IL5*Vc_mepo) : 0;
double E_benra    = BENRA_C>0  ? Kmax_benra*BENRA_C/(IC50_benra+BENRA_C) : 0;
dxdt_TH2 = kprod_Th2*(1-E_pred_Th2)*(1-0.7*E_cyclo)*(1-0.5*E_ritu) -kdeg_Th2*TH2;
double IL5_sink = (kon_mepo*MEPO_C*IL5*Vc_mepo) -(koff_mepo*MEPO_CPX);
dxdt_IL5 = kprod_IL5*(TH2/Th2_0) -kdeg_IL5*IL5 -IL5_sink;
if(IL5<0) dxdt_IL5=0;
double kprod_EosB = kprod_Eos*(IL5/IL5_0)*(1-frac_mepo)*(1-E_benra)*(1-E_pred_Eos*0.80);
dxdt_EOS_B = kprod_EosB -(kdeg_Eos+kmig_Eos)*EOS_B;
if(EOS_B<0) dxdt_EOS_B=0;
dxdt_EOS_T = kmig_Eos*EOS_B*0.004*(1-E_benra)*(1-E_pred_Eos*0.70) -kdeg_EosT*EOS_T;
if(EOS_T<0) dxdt_EOS_T=0;
dxdt_IGE  = kprod_IgE*(TH2/Th2_0)*(1-0.5*E_ritu)*(1-0.4*E_cyclo) -kdeg_IgE*IGE;
dxdt_ANCA = ANCA_pos*(kprod_ANCA*(TH2/Th2_0)*(1-E_ritu)*(1-E_cyclo)*(1-0.4*E_pred_Th2)) -kdeg_ANCA*ANCA;
if(ANCA<0) dxdt_ANCA=0;
double Vasc_drive = kact_Vasc*(EOS_T/EosT_0 + ANCA_pos*1.5*ANCA/(ANCA+0.1));
double Vasc_ther  = kres_Vasc*(E_pred_Th2*0.6+E_cyclo*0.3+E_ritu*0.1);
dxdt_VASC = Vasc_drive*(10-VASC)/10 -(kres_Vasc+Vasc_ther)*VASC;
if(VASC<0) dxdt_VASC=0; if(VASC>10) dxdt_VASC=-kres_Vasc*VASC;
dxdt_CARD = kact_Card*(EOS_T/EosT_0)*(10-CARD)/10 -(kres_Card+0.4*E_pred_Eos)*CARD;
if(CARD<0) dxdt_CARD=0;
dxdt_NERV = kact_Nerv*(VASC/5+EOS_T/(EosT_0*2))*(10-NERV)/10 -(kres_Nerv+0.3*E_pred_Eos)*NERV;
if(NERV<0) dxdt_NERV=0;

$TABLE
capture BloodEos = EOS_B;
capture TissueEos = EOS_T;
capture IL5_lev   = IL5;
capture IgE_lev   = IGE;
capture ANCA_lev  = ANCA;
capture Vasc_act  = VASC;
capture Card_dmg  = CARD;
capture Nerv_dmg  = NERV;
capture MEPO_serum = MEPO_C;
capture BENRA_serum= BENRA_C;
capture PRED_ng   = PRED_C * 360.44;
capture CYCLO_ug  = CYCLO_M;
capture RITU_ug   = RITU_C;
double BVAS = 0;
BVAS += (EOS_B>1500)?3:(EOS_B>500?1:0);
BVAS += VASC*2.5;
BVAS += (ANCA>0.5)?5:0;
BVAS += CARD*1.5;
BVAS += NERV*1.5;
if(BVAS>63) BVAS=63;
capture BVAS_score = BVAS;
double FEV1 = 60 + 40*exp(-0.5*EOS_T/EosT_0)*(1-0.3*(1-exp(-VASC/5)));
capture FEV1_pct = FEV1;
double LVEF = 65 - CARD*4.5; if(LVEF<10) LVEF=10;
capture LVEF_pct = LVEF;
double eGFR_val = 90 - VASC*4 - ANCA_pos*ANCA*3; if(eGFR_val<5) eGFR_val=5;
capture eGFR = eGFR_val;
capture Remission = (BVAS<=1) ? 1.0 : 0.0;
capture EosSuppPct = 100*(1 - EOS_B/5000);
'

mod_egpa <- mcode("EGPA_shiny", egpa_code)

## ---- Helper: run simulation ----
run_sim <- function(anca_pos, eos_b0, vasc0, card0, nerv0, anca0,
                    use_pred, pred_dose_mg,
                    use_mepo, use_benra, use_cyclo, use_ritu,
                    sim_weeks = 104) {

  end_h <- sim_weeks * 168

  # Build events
  ev_all <- NULL

  if (use_pred && pred_dose_mg > 0) {
    taper_times <- seq(8, 26, by = 4) * 168
    taper_doses <- seq(pred_dose_mg, 7.5, length.out = length(taper_times))
    ev_pred_ind  <- ev(cmt = "PRED_DEPOT", amt = pred_dose_mg * 2.77,
                       ii = 24, addl = 7*8-1, time = 0)
    ev_pred_maint <- ev(cmt = "PRED_DEPOT", amt = 7.5 * 2.77,
                        ii = 24, addl = 9999, time = 26*168)
    ev_all <- ev_pred_ind + ev_pred_maint
  }

  if (use_mepo) {
    ev_m <- ev(cmt = "MEPO_DEPOT", amt = 300,
               ii = 28*24, addl = floor(sim_weeks/4), time = 0)
    ev_all <- if (is.null(ev_all)) ev_m else ev_all + ev_m
  }

  if (use_benra) {
    ev_b1 <- ev(cmt = "BENRA_DEPOT", amt = 30, ii = 28*24, addl = 2, time = 0)
    ev_b2 <- ev(cmt = "BENRA_DEPOT", amt = 30, ii = 56*24,
                addl = floor(sim_weeks/8), time = 3*28*24)
    ev_all <- if (is.null(ev_all)) (ev_b1+ev_b2) else ev_all + ev_b1 + ev_b2
  }

  if (use_cyclo) {
    for (i in 0:5) {
      evc <- ev(cmt = "CYCLO_C", amt = 750, time = i * 28 * 24)
      ev_all <- if (is.null(ev_all)) evc else ev_all + evc
    }
  }

  if (use_ritu) {
    ev_r <- ev(cmt = "RITU_C", amt = 1000, time = c(0, 14*24, 26*7*24, (26*7+14)*24))
    ev_all <- if (is.null(ev_all)) ev_r else ev_all + ev_r
  }

  if (is.null(ev_all)) {
    ev_all <- ev(time = 0, amt = 0, cmt = 1)
  }

  mod_egpa %>%
    param(ANCA_pos = as.numeric(anca_pos)) %>%
    init(EOS_B = eos_b0, VASC = vasc0, CARD = card0,
         NERV = nerv0, ANCA = anca0 * as.numeric(anca_pos)) %>%
    mrgsim(events = ev_all, end = end_h, delta = 24) %>%
    as_tibble() %>%
    mutate(Week = time / 168)
}

## ============================================================
## UI
## ============================================================
ui <- dashboardPage(
  skin = "red",

  dashboardHeader(
    title = "EGPA QSP Dashboard",
    titleWidth = 280
  ),

  dashboardSidebar(
    width = 280,
    sidebarMenu(
      menuItem("1. Patient Profile",    tabName = "profile",   icon = icon("user")),
      menuItem("2. Pharmacokinetics",   tabName = "pk",        icon = icon("flask")),
      menuItem("3. Eosinophil / IL-5",  tabName = "eos",       icon = icon("dna")),
      menuItem("4. Vasculitis / Organs",tabName = "vasc",      icon = icon("heart")),
      menuItem("5. Clinical Endpoints", tabName = "endpoints", icon = icon("chart-line")),
      menuItem("6. Scenario Compare",   tabName = "scenario",  icon = icon("layer-group")),
      menuItem("7. Biomarkers",         tabName = "bio",       icon = icon("microscope"))
    ),
    hr(),
    ## Patient / Disease Parameters
    h4("Patient Profile", style = "color:white; padding-left:15px;"),
    checkboxInput("anca_pos",  "ANCA-positive (anti-MPO)",  value = FALSE),
    sliderInput("eos_b0",    "Initial Blood Eosinophils",
                min = 500, max = 15000, value = 5000, step = 500),
    sliderInput("vasc0",     "Initial Vasculitis Activity",
                min = 0, max = 10, value = 3.0, step = 0.5),
    sliderInput("card0",     "Initial Cardiac Damage",
                min = 0, max = 10, value = 0.5, step = 0.5),
    sliderInput("nerv0",     "Initial Nerve Damage",
                min = 0, max = 10, value = 1.0, step = 0.5),
    conditionalPanel(
      condition = "input.anca_pos == true",
      sliderInput("anca0", "Initial ANCA Level (U/mL)",
                  min = 0, max = 5, value = 1.0, step = 0.1)
    ),
    hr(),
    ## Treatment Parameters
    h4("Treatment", style = "color:white; padding-left:15px;"),
    checkboxInput("use_pred",  "Prednisolone",       value = TRUE),
    sliderInput("pred_dose", "Pred Initial Dose (mg)",
                min = 10, max = 80, value = 50, step = 5),
    checkboxInput("use_mepo",  "Mepolizumab 300mg SC", value = FALSE),
    checkboxInput("use_benra", "Benralizumab 30mg SC", value = FALSE),
    checkboxInput("use_cyclo", "Cyclophosphamide IV",  value = FALSE),
    checkboxInput("use_ritu",  "Rituximab IV",         value = FALSE),
    hr(),
    sliderInput("sim_weeks", "Simulation (weeks)",
                min = 12, max = 156, value = 104, step = 4),
    actionButton("run_btn", "Run Simulation",
                 icon = icon("play"),
                 style = "background-color:#dc143c; color:white; width:100%;")
  ),

  dashboardBody(
    tabItems(

      ## ---- Tab 1: Patient Profile ----
      tabItem(tabName = "profile",
        fluidRow(
          box(title = "EGPA Disease Overview", width = 8, solidHeader = TRUE, status = "danger",
            HTML("
            <h4>Eosinophilic Granulomatosis with Polyangiitis (EGPA)</h4>
            <p>Formerly Churg-Strauss syndrome — a rare ANCA-associated vasculitis (≈1-3/million/year)
            characterised by three phases:</p>
            <ul>
              <li><b>Prodromal:</b> Asthma (often severe) + allergic rhinosinusitis + nasal polyps</li>
              <li><b>Eosinophilic:</b> Blood hypereosinophilia (>10% or >1.5×10⁹/L) + tissue eosinophilia
                  (lungs, heart, GI)</li>
              <li><b>Vasculitic:</b> Necrotizing small-medium vessel vasculitis — mononeuritis multiplex,
                  cardiac, renal, skin</li>
            </ul>
            <p><b>ANCA phenotype (~40%):</b> Anti-MPO positive; higher renal/neural involvement.
            <b>ANCA-negative (~60%):</b> Higher eosinophilic (cardiac/pulmonary) damage.</p>
            <h4>ACR/EULAR 2022 Classification Criteria</h4>
            <ul>
              <li>Obstructive airway disease</li>
              <li>Nasal polyps</li>
              <li>Peripheral neuropathy attributed to vasculitis</li>
              <li>Blood eosinophilia >1×10⁹/L</li>
              <li>Extravascular eosinophil-predominant inflammation on biopsy</li>
              <li>Positive ANCA (anti-MPO)</li>
            </ul>
            <h4>Treatment Algorithm</h4>
            <ul>
              <li>All patients: <b>Corticosteroids</b> (prednisolone 0.5–1 mg/kg/day, tapered)</li>
              <li>Non-severe/relapsing: <b>Mepolizumab</b> 300mg SC q4w (FDA-approved 2017)</li>
              <li>Severe (FFS ≥2) / ANCA+: <b>Cyclophosphamide</b> or <b>Rituximab</b></li>
              <li>Refractory: Benralizumab, dupilumab (investigational)</li>
            </ul>
            "),
            tags$hr(),
            h4("Five-Factor Score (FFS) — Prognosis"),
            HTML("
            <table border='1' style='width:100%; border-collapse:collapse;'>
              <tr><th>Factor</th><th>Points</th></tr>
              <tr><td>Creatinine > 150 µmol/L</td><td>+1</td></tr>
              <tr><td>Proteinuria > 1g/day</td><td>+1</td></tr>
              <tr><td>GI involvement (bleeding/perforation)</td><td>+1</td></tr>
              <tr><td>Cardiomyopathy</td><td>+1</td></tr>
              <tr><td>CNS involvement</td><td>+1</td></tr>
            </table>
            <p>FFS 0: 5-year mortality ~11% | FFS ≥2: ~26%</p>
            ")
          ),
          box(title = "Current Patient Parameters", width = 4,
              solidHeader = TRUE, status = "warning",
              tableOutput("patient_summary_tbl")
          )
        )
      ),

      ## ---- Tab 2: PK ----
      tabItem(tabName = "pk",
        fluidRow(
          box(title = "Drug Serum Concentrations", width = 12,
              solidHeader = TRUE, status = "danger",
              plotlyOutput("pk_plot", height = "450px"))
        ),
        fluidRow(
          box(title = "Prednisolone Plasma (ng/mL)", width = 6,
              plotlyOutput("pk_pred_plot", height = "300px")),
          box(title = "PK Parameters Summary", width = 6,
              tableOutput("pk_params_tbl"))
        )
      ),

      ## ---- Tab 3: Eosinophil / IL-5 ----
      tabItem(tabName = "eos",
        fluidRow(
          box(title = "Blood Eosinophil Count (cells/µL)", width = 6,
              solidHeader = TRUE, status = "danger",
              plotlyOutput("eos_blood_plot", height = "350px")),
          box(title = "Tissue Eosinophil Burden (AU)", width = 6,
              solidHeader = TRUE, status = "warning",
              plotlyOutput("eos_tissue_plot", height = "350px"))
        ),
        fluidRow(
          box(title = "Free IL-5 Level (nM)", width = 6,
              plotlyOutput("il5_plot", height = "300px")),
          box(title = "Total IgE (kU/L)", width = 6,
              plotlyOutput("ige_plot", height = "300px"))
        )
      ),

      ## ---- Tab 4: Vasculitis & Organ Damage ----
      tabItem(tabName = "vasc",
        fluidRow(
          box(title = "Vasculitis Activity (0-10)", width = 6,
              solidHeader = TRUE, status = "danger",
              plotlyOutput("vasc_plot", height = "350px")),
          box(title = "ANCA Level (U/mL) — ANCA+ only", width = 6,
              solidHeader = TRUE, status = "warning",
              plotlyOutput("anca_plot", height = "350px"))
        ),
        fluidRow(
          box(title = "Cardiac Damage Score (0-10)", width = 6,
              plotlyOutput("card_plot", height = "300px")),
          box(title = "Peripheral Nerve Damage Score (0-10)", width = 6,
              plotlyOutput("nerv_plot", height = "300px"))
        )
      ),

      ## ---- Tab 5: Clinical Endpoints ----
      tabItem(tabName = "endpoints",
        fluidRow(
          box(title = "BVAS Score (Disease Activity)", width = 6,
              solidHeader = TRUE, status = "danger",
              plotlyOutput("bvas_plot", height = "350px")),
          box(title = "FEV1% Predicted (Pulmonary)", width = 6,
              solidHeader = TRUE, status = "primary",
              plotlyOutput("fev1_plot", height = "350px"))
        ),
        fluidRow(
          box(title = "LVEF% (Cardiac Function)", width = 4,
              plotlyOutput("lvef_plot", height = "300px")),
          box(title = "eGFR (mL/min/1.73m²)", width = 4,
              plotlyOutput("egfr_plot", height = "300px")),
          box(title = "Clinical Endpoints at Week 52", width = 4,
              tableOutput("endpoints_tbl"))
        )
      ),

      ## ---- Tab 6: Scenario Comparison ----
      tabItem(tabName = "scenario",
        fluidRow(
          box(title = "Multi-Scenario: Blood Eosinophils", width = 6,
              solidHeader = TRUE, status = "danger",
              plotlyOutput("sc_eos_plot", height = "350px")),
          box(title = "Multi-Scenario: BVAS Score", width = 6,
              solidHeader = TRUE, status = "primary",
              plotlyOutput("sc_bvas_plot", height = "350px"))
        ),
        fluidRow(
          box(title = "Multi-Scenario Summary at Week 52", width = 12,
              DTOutput("scenario_tbl"))
        )
      ),

      ## ---- Tab 7: Biomarkers ----
      tabItem(tabName = "bio",
        fluidRow(
          box(title = "Eosinophil Suppression from Baseline (%)", width = 6,
              solidHeader = TRUE, status = "danger",
              plotlyOutput("eos_sup_plot", height = "350px")),
          box(title = "Remission Status over Time", width = 6,
              solidHeader = TRUE, status = "success",
              plotlyOutput("remission_plot", height = "350px"))
        ),
        fluidRow(
          box(title = "Biomarker Heatmap at Selected Timepoints", width = 12,
              plotlyOutput("bio_heat_plot", height = "350px"))
        )
      )
    )
  )
)

## ============================================================
## SERVER
## ============================================================
server <- function(input, output, session) {

  ## Reactive: single scenario (patient parameters)
  sim_data <- eventReactive(input$run_btn, {
    withProgress(message = "Running simulation...", value = 0.5, {
      run_sim(
        anca_pos  = input$anca_pos,
        eos_b0    = input$eos_b0,
        vasc0     = input$vasc0,
        card0     = input$card0,
        nerv0     = input$nerv0,
        anca0     = if (input$anca_pos) input$anca0 else 0,
        use_pred  = input$use_pred,
        pred_dose_mg = input$pred_dose,
        use_mepo  = input$use_mepo,
        use_benra = input$use_benra,
        use_cyclo = input$use_cyclo,
        use_ritu  = input$use_ritu,
        sim_weeks = input$sim_weeks
      )
    })
  }, ignoreNULL = FALSE)

  ## Multi-scenario (Tab 6): run 4 standard scenarios
  multi_sim <- eventReactive(input$run_btn, {
    withProgress(message = "Running scenarios...", value = 0.5, {
      scen <- list(
        "Prednisolone Only" = run_sim(FALSE, input$eos_b0, input$vasc0, input$card0, input$nerv0, 0,
                                       TRUE, 50, FALSE, FALSE, FALSE, FALSE, input$sim_weeks),
        "Mepolizumab + Pred" = run_sim(FALSE, input$eos_b0, input$vasc0, input$card0, input$nerv0, 0,
                                        TRUE, 50, TRUE, FALSE, FALSE, FALSE, input$sim_weeks),
        "Benralizumab + Pred" = run_sim(FALSE, input$eos_b0, input$vasc0, input$card0, input$nerv0, 0,
                                         TRUE, 50, FALSE, TRUE, FALSE, FALSE, input$sim_weeks),
        "Cyclophosp. + Pred (ANCA+)" = run_sim(TRUE, input$eos_b0, max(input$vasc0, 5), input$card0, input$nerv0, 1.0,
                                                 TRUE, 60, FALSE, FALSE, TRUE, FALSE, input$sim_weeks)
      )
      bind_rows(lapply(names(scen), function(n) mutate(scen[[n]], Scenario = n)))
    })
  }, ignoreNULL = FALSE)

  ## ---- Patient summary table ----
  output$patient_summary_tbl <- renderTable({
    data.frame(
      Parameter = c("Blood Eos (initial)", "Vasculitis Score", "Cardiac Score",
                    "Nerve Score", "ANCA Status",
                    "Treatment", "Simulation"),
      Value = c(
        paste0(input$eos_b0, " cells/µL"),
        paste0(input$vasc0, "/10"),
        paste0(input$card0, "/10"),
        paste0(input$nerv0, "/10"),
        if (input$anca_pos) "Positive (anti-MPO)" else "Negative",
        paste0(c("Pred"[input$use_pred], "Mepo"[input$use_mepo],
                 "Benra"[input$use_benra], "Cyclo"[input$use_cyclo],
                 "Ritu"[input$use_ritu]), collapse = " + "),
        paste0(input$sim_weeks, " weeks")
      )
    )
  })

  ## ---- PK Plots ----
  output$pk_plot <- renderPlotly({
    d <- sim_data()
    gg <- ggplot(d, aes(x = Week)) +
      geom_line(aes(y = MEPO_serum,  color = "Mepolizumab (µg/mL)"),  size = 0.9) +
      geom_line(aes(y = BENRA_serum, color = "Benralizumab (µg/mL)"), size = 0.9) +
      geom_line(aes(y = RITU_ug,     color = "Rituximab (µg/mL)"),    size = 0.9) +
      labs(x = "Week", y = "Serum Concentration", color = "Drug") +
      scale_color_manual(values = c("Mepolizumab (µg/mL)" = "#E63946",
                                    "Benralizumab (µg/mL)" = "#2A9D8F",
                                    "Rituximab (µg/mL)"    = "#457B9D")) +
      theme_bw()
    ggplotly(gg)
  })

  output$pk_pred_plot <- renderPlotly({
    d <- sim_data()
    gg <- ggplot(d, aes(x = Week, y = PRED_ng)) +
      geom_line(color = "#E63946", size = 1) +
      geom_hline(yintercept = 120, linetype = "dashed", color = "gray") +
      labs(x = "Week", y = "Prednisolone (ng/mL)") +
      theme_bw()
    ggplotly(gg)
  })

  output$pk_params_tbl <- renderTable({
    data.frame(
      Drug = c("Mepolizumab", "Benralizumab", "Prednisolone", "Cyclophosphamide", "Rituximab"),
      `Route/Dosing` = c("SC 300mg q4w", "SC 30mg q4w→q8w", "Oral taper", "IV 750mg/m² ×6", "IV 1000mg ×2"),
      `Half-life` = c("~16-22 d", "~15 d", "~2-3 h", "~7 h", "~22 d"),
      `Key Target` = c("IL-5", "IL-5Rα (ADCC)", "GR (broad)", "DNA alkylation", "CD20 (ADCC/CDC)")
    )
  })

  ## ---- Eosinophil Plots ----
  eos_plot_fn <- function(data, yvar, title, ylab, hline = NULL, hline_lab = NULL) {
    gg <- ggplot(data, aes_string(x = "Week", y = yvar)) +
      geom_line(color = "#8B008B", size = 1.1)
    if (!is.null(hline)) {
      gg <- gg + geom_hline(yintercept = hline, linetype = "dashed", color = "gray40") +
        annotate("text", x = max(data$Week) * 0.8, y = hline * 1.05, label = hline_lab, size = 3)
    }
    gg + labs(title = title, x = "Week", y = ylab) + theme_bw()
  }

  output$eos_blood_plot <- renderPlotly({
    ggplotly(eos_plot_fn(sim_data(), "BloodEos", "",
                         "Blood Eosinophils (cells/µL)",
                         hline = 500, hline_lab = "Normal (<500)"))
  })
  output$eos_tissue_plot <- renderPlotly({
    ggplotly(eos_plot_fn(sim_data(), "TissueEos", "", "Tissue Eosinophils (AU)"))
  })
  output$il5_plot <- renderPlotly({
    ggplotly(eos_plot_fn(sim_data(), "IL5_lev", "", "Free IL-5 (nM)"))
  })
  output$ige_plot <- renderPlotly({
    ggplotly(eos_plot_fn(sim_data(), "IgE_lev", "", "Total IgE (kU/L)"))
  })

  ## ---- Vasculitis / Organ Plots ----
  output$vasc_plot <- renderPlotly({
    ggplotly(eos_plot_fn(sim_data(), "Vasc_act", "", "Vasculitis Activity (0-10)"))
  })
  output$anca_plot <- renderPlotly({
    ggplotly(eos_plot_fn(sim_data(), "ANCA_lev", "", "ANCA Level (U/mL)"))
  })
  output$card_plot <- renderPlotly({
    ggplotly(eos_plot_fn(sim_data(), "Card_dmg", "", "Cardiac Damage (0-10)"))
  })
  output$nerv_plot <- renderPlotly({
    ggplotly(eos_plot_fn(sim_data(), "Nerv_dmg", "", "Nerve Damage (0-10)"))
  })

  ## ---- Clinical Endpoints ----
  output$bvas_plot <- renderPlotly({
    d <- sim_data()
    gg <- ggplot(d, aes(x = Week, y = BVAS_score)) +
      geom_line(color = "#E63946", size = 1.1) +
      geom_ribbon(aes(ymin = 0, ymax = 1), alpha = 0.15, fill = "green") +
      annotate("text", x = max(d$Week)*0.7, y = 0.5, label = "Remission zone (≤1)", size = 3) +
      labs(x = "Week", y = "BVAS (0-63)") + theme_bw()
    ggplotly(gg)
  })
  output$fev1_plot <- renderPlotly({
    d <- sim_data()
    gg <- ggplot(d, aes(x = Week, y = FEV1_pct)) +
      geom_line(color = "#457B9D", size = 1.1) +
      geom_hline(yintercept = 80, linetype = "dashed") +
      annotate("text", x = max(d$Week)*0.7, y = 82, label = "Normal (≥80%)", size = 3) +
      labs(x = "Week", y = "FEV1 (%predicted)") + ylim(30, 100) + theme_bw()
    ggplotly(gg)
  })
  output$lvef_plot <- renderPlotly({
    d <- sim_data()
    gg <- ggplot(d, aes(x = Week, y = LVEF_pct)) +
      geom_line(color = "#E76F51", size = 1) +
      geom_hline(yintercept = 55, linetype = "dashed") +
      labs(x = "Week", y = "LVEF (%)") + ylim(0, 75) + theme_bw()
    ggplotly(gg)
  })
  output$egfr_plot <- renderPlotly({
    d <- sim_data()
    gg <- ggplot(d, aes(x = Week, y = eGFR)) +
      geom_line(color = "#264653", size = 1) +
      geom_hline(yintercept = 60, linetype = "dashed") +
      labs(x = "Week", y = "eGFR (mL/min/1.73m²)") + theme_bw()
    ggplotly(gg)
  })
  output$endpoints_tbl <- renderTable({
    d <- sim_data()
    w52 <- d[which.min(abs(d$Week - 52)), ]
    data.frame(
      Endpoint = c("Blood Eos", "BVAS", "FEV1%", "LVEF%", "eGFR", "Remission"),
      Value = c(round(w52$BloodEos), round(w52$BVAS_score, 1),
                round(w52$FEV1_pct, 1), round(w52$LVEF_pct, 1),
                round(w52$eGFR, 1),
                ifelse(w52$Remission == 1, "YES", "NO"))
    )
  })

  ## ---- Scenario Comparison ----
  output$sc_eos_plot <- renderPlotly({
    d <- multi_sim()
    gg <- ggplot(d, aes(x = Week, y = BloodEos, color = Scenario)) +
      geom_line(size = 1) +
      geom_hline(yintercept = 500, linetype = "dashed") +
      labs(x = "Week", y = "Blood Eosinophils (cells/µL)") +
      scale_color_brewer(palette = "Set2") + theme_bw() +
      theme(legend.position = "bottom")
    ggplotly(gg)
  })
  output$sc_bvas_plot <- renderPlotly({
    d <- multi_sim()
    gg <- ggplot(d, aes(x = Week, y = BVAS_score, color = Scenario)) +
      geom_line(size = 1) +
      labs(x = "Week", y = "BVAS Score (0-63)") +
      scale_color_brewer(palette = "Set2") + theme_bw() +
      theme(legend.position = "bottom")
    ggplotly(gg)
  })
  output$scenario_tbl <- renderDT({
    d <- multi_sim()
    w52 <- d %>%
      group_by(Scenario) %>%
      filter(abs(Week - 52) == min(abs(Week - 52))) %>%
      slice(1) %>%
      select(Scenario, BloodEos, BVAS_score, FEV1_pct, LVEF_pct, eGFR, Remission) %>%
      mutate(across(where(is.numeric), ~round(., 1)),
             Remission = ifelse(Remission == 1, "YES", "NO"))
    names(w52) <- c("Scenario", "Blood Eos", "BVAS", "FEV1%", "LVEF%", "eGFR", "Remission")
    datatable(w52, options = list(pageLength = 10, dom = 't'), rownames = FALSE)
  })

  ## ---- Biomarkers ----
  output$eos_sup_plot <- renderPlotly({
    d <- sim_data()
    gg <- ggplot(d, aes(x = Week, y = EosSuppPct)) +
      geom_line(color = "#8B008B", size = 1.1) +
      geom_hline(yintercept = 90, linetype = "dashed", color = "blue") +
      annotate("text", x = max(d$Week)*0.5, y = 92, label = "90% target", size = 3) +
      labs(x = "Week", y = "Eosinophil Suppression (%)") +
      ylim(-20, 100) + theme_bw()
    ggplotly(gg)
  })
  output$remission_plot <- renderPlotly({
    d <- sim_data()
    gg <- ggplot(d, aes(x = Week, y = Remission)) +
      geom_step(color = "#2A9D8F", size = 1.2) +
      geom_ribbon(aes(ymin = 0, ymax = Remission), fill = "#2A9D8F", alpha = 0.3) +
      scale_y_continuous(breaks = c(0, 1), labels = c("Active", "Remission")) +
      labs(x = "Week", y = "Disease Status") + theme_bw()
    ggplotly(gg)
  })
  output$bio_heat_plot <- renderPlotly({
    d <- sim_data()
    timepoints <- c(0, 4, 12, 26, 52, 78, 104)
    d_sub <- d %>%
      filter(sapply(Week, function(w) any(abs(w - timepoints) < 0.6))) %>%
      group_by(Week = round(Week)) %>%
      slice(1)
    # Normalize biomarkers to 0-1
    bio_vars <- c("BloodEos", "BVAS_score", "FEV1_pct", "LVEF_pct", "eGFR", "IgE_lev", "ANCA_lev")
    d_long <- d_sub %>%
      select(Week, all_of(intersect(bio_vars, names(d)))) %>%
      tidyr::pivot_longer(-Week, names_to = "Biomarker", values_to = "Value")
    d_long <- d_long %>%
      group_by(Biomarker) %>%
      mutate(Norm = (Value - min(Value)) / (max(Value) - min(Value) + 1e-10))
    gg <- ggplot(d_long, aes(x = factor(Week), y = Biomarker, fill = Norm)) +
      geom_tile(color = "white") +
      scale_fill_gradient2(low = "#2A9D8F", mid = "#F4A261", high = "#E63946",
                           midpoint = 0.5) +
      labs(x = "Week", fill = "Normalized\nValue") + theme_bw()
    ggplotly(gg)
  })
}

## ============================================================
## Run App
## ============================================================
shinyApp(ui = ui, server = server)
