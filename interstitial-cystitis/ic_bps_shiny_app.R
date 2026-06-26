## =============================================================================
## IC/BPS Shiny Dashboard
## Disease: Interstitial Cystitis / Bladder Pain Syndrome
## Version: 1.0 | Date: 2026-06-26 | CCR/Claude
##
## Tabs (8):
##   1. Overview & Disease Map
##   2. Patient Profile
##   3. PK — Drug Concentrations
##   4. PD — Biomarkers
##   5. Clinical Endpoints
##   6. Scenario Comparison
##   7. Subtype Explorer (Hunner vs Non-Hunner)
##   8. Sensitivity Analysis
## =============================================================================

library(shiny)
library(shinydashboard)
library(mrgsolve)
library(dplyr)
library(ggplot2)
library(tidyr)
library(DT)
library(plotly)

## ---------------------------------------------------------------------------
## Embed model code (stripped-down from mrgsolve model file)
## ---------------------------------------------------------------------------
ic_bps_code <- '
$PROB IC/BPS QSP Shiny Model

$PARAM @annotated
k_GAG_syn:0.05:GAG synthesis
k_GAG_deg:0.08:GAG degradation
k_PERM_up:0.12:Permeability increase
k_PERM_res:0.04:Permeability restore
k_MC_act:0.10:Mast cell activation
k_MC_res:0.06:Mast cell resolution
k_HIST_rel:0.20:Histamine release
k_HIST_clr:0.30:Histamine clearance
k_SP_prod:0.15:SP production
k_SP_clr:0.25:SP clearance
k_NGF_prod:0.08:NGF production
k_NGF_clr:0.12:NGF clearance
k_IL6_prod:0.10:IL-6 production
k_IL6_clr:0.15:IL-6 clearance
k_TNF_prod:0.10:TNF production
k_TNF_clr:0.18:TNF clearance
k_CFIB_up:0.08:C-fiber sensitization
k_CFIB_res:0.03:C-fiber resolution
k_SPIN_up:0.05:Spinal sensitization
k_SPIN_res:0.02:Spinal resolution
k_CENT_up:0.04:Central sensitization
k_CENT_res:0.015:Central resolution
k_CAP_loss:0.003:Bladder cap loss
k_CAP_res:0.002:Bladder cap recovery
F_PPS:0.06:PPS bioavailability
ka_PPS:1.2:PPS absorption
CL_PPS:3.5:PPS clearance
V_PPS:8.0:PPS volume
Emax_GAG_PPS:0.8:PPS max GAG effect
IC50_GAG_PPS:0.5:PPS IC50 GAG
F_HYD:0.80:HYD bioavailability
ka_HYD:2.5:HYD absorption
CL_HYD:35.0:HYD clearance
V_HYD:350.0:HYD volume
IC50_HIST_HYD:10.0:HYD IC50 histamine
F_CsA:0.35:CsA bioavailability
ka_CsA:1.8:CsA absorption
CL_CsA:25.0:CsA clearance
V_CsA:400.0:CsA volume
IC50_Tcell_CsA:150.0:CsA IC50 T cell
F_AMI:0.50:AMI bioavailability
ka_AMI:2.0:AMI absorption
CL_AMI:110.0:AMI clearance
V_AMI:1500.0:AMI volume
IC50_HIST_AMI:25.0:AMI H1 IC50
EC50_PAIN_AMI:30.0:AMI pain IC50
BoNTA_eff:0.0:BoNTA effect
DMSO_eff:0.0:DMSO effect
Hunner:0.0:Hunner subtype
OLS0:14.0:Baseline OLS
PAIN0:6.5:Baseline pain

$INIT
PPS_GUT:0 PPS_CENT:0 HYD_GUT:0 HYD_CENT:0
CSA_GUT:0 CSA_CENT:0 AMI_GUT:0 AMI_CENT:0
GAG:0.35 PERM:0.65 MC:0.70 HIST:0.55
SP:0.60 NGF:0.55 IL6:0.50 TNF:0.45
C_FIBER:0.65 SPINAL:0.55 CENTRAL:0.50
CAP:180 PAIN:6.5 OLS:14.0

$ODE
double PPS_conc = PPS_CENT/V_PPS;
double PPS_urine = PPS_conc * 0.25;
dxdt_PPS_GUT  = -ka_PPS*PPS_GUT;
dxdt_PPS_CENT = F_PPS*ka_PPS*PPS_GUT - CL_PPS*PPS_conc;
double HYD_conc = HYD_CENT/V_HYD;
dxdt_HYD_GUT  = -ka_HYD*HYD_GUT;
dxdt_HYD_CENT = F_HYD*ka_HYD*HYD_GUT - CL_HYD*HYD_conc;
double CSA_conc = CSA_CENT/V_CsA;
dxdt_CSA_GUT  = -ka_CsA*CSA_GUT;
dxdt_CSA_CENT = F_CsA*ka_CsA*CSA_GUT - CL_CsA*CSA_conc;
double AMI_conc = AMI_CENT/V_AMI;
dxdt_AMI_GUT  = -ka_AMI*AMI_GUT;
dxdt_AMI_CENT = F_AMI*ka_AMI*AMI_GUT - CL_AMI*AMI_conc;
double E_PPS = Emax_GAG_PPS*PPS_urine/(IC50_GAG_PPS+PPS_urine);
double E_HYD = HYD_conc/(IC50_HIST_HYD+HYD_conc);
double E_CsA = CSA_conc/(IC50_Tcell_CsA+CSA_conc);
double E_AMI_H1 = AMI_conc/(IC50_HIST_AMI+AMI_conc);
double E_AMI_pain = AMI_conc/(EC50_PAIN_AMI+AMI_conc);
double E_BoNTA = BoNTA_eff*0.7;
double E_DMSO_MC = DMSO_eff*0.6;
double E_DMSO_IL6 = DMSO_eff*0.4;
double GAG_syn = k_GAG_syn*(1-GAG)+k_GAG_syn*E_PPS;
double GAG_deg = k_GAG_deg*PERM*(1+MC*0.5);
dxdt_GAG = GAG_syn - GAG_deg;
double PERM_up = k_PERM_up*(1-GAG)*(1+TNF*0.5+MC*0.3);
double PERM_res = k_PERM_res*GAG*(1+E_PPS*0.5);
dxdt_PERM = PERM_up - PERM_res;
double MC_act = k_MC_act*PERM*(1+SP*0.3)*(1-E_HYD*0.5)*(1-E_DMSO_MC);
double MC_res = k_MC_res*(1+E_HYD+E_DMSO_MC);
dxdt_MC = MC_act - MC_res*MC;
double HIST_rel = k_HIST_rel*MC;
double HIST_clr = k_HIST_clr*(1+E_HYD+E_AMI_H1*0.5);
dxdt_HIST = HIST_rel - HIST_clr*HIST;
double SP_prod = k_SP_prod*C_FIBER*(1+MC*0.2)*(1-E_BoNTA);
dxdt_SP = SP_prod - k_SP_clr*SP;
double NGF_prod = k_NGF_prod*(IL6+TNF)*(1+PERM*0.3);
dxdt_NGF = NGF_prod - k_NGF_clr*NGF;
double IL6_prod = k_IL6_prod*(MC+PERM*0.5)*(1+Hunner*0.8)*(1-E_CsA*0.4)*(1-E_DMSO_IL6);
dxdt_IL6 = IL6_prod - k_IL6_clr*IL6;
double TNF_prod = k_TNF_prod*(MC+IL6*0.3)*(1+Hunner*0.5)*(1-E_CsA*0.3)*(1-E_DMSO_MC*0.3);
dxdt_TNF = TNF_prod - k_TNF_clr*TNF;
double CFIB_drive = k_CFIB_up*(SP*0.4+HIST*0.3+NGF*0.3)*(1-E_BoNTA*0.5);
dxdt_C_FIBER = CFIB_drive - k_CFIB_res*C_FIBER;
double SPIN_drive = k_SPIN_up*C_FIBER*(1+SPINAL*0.2);
dxdt_SPINAL = SPIN_drive - k_SPIN_res*(1+E_AMI_pain*0.5)*SPINAL;
double CENT_drive = k_CENT_up*SPINAL*(1+CENTRAL*0.1);
dxdt_CENTRAL = CENT_drive - k_CENT_res*(1+E_AMI_pain*0.3)*CENTRAL;
dxdt_CAP = k_CAP_res*300 - k_CAP_loss*(IL6*0.5+MC*0.3+SPINAL*0.2)*CAP;
double PAIN_target = 10*(0.4*CENTRAL+0.3*SPINAL+0.3*SP)/(1+0.4*CENTRAL+0.3*SPINAL+0.3*SP);
dxdt_PAIN = 0.15*(PAIN_target-PAIN) - E_AMI_pain*PAIN*0.1;
double OLS_target = 0.5*(1-CAP/400)*20 + 0.5*PAIN;
if(OLS_target>20) OLS_target=20; if(OLS_target<0) OLS_target=0;
dxdt_OLS = 0.10*(OLS_target-OLS);

$TABLE
double PPS_conc_out = PPS_CENT/V_PPS;
double HYD_conc_out = HYD_CENT/V_HYD;
double CSA_conc_out = CSA_CENT/V_CsA;
double AMI_conc_out = AMI_CENT/V_AMI;
double FREQ = 1440.0/CAP*10.0;
if(FREQ<6) FREQ=6; if(FREQ>40) FREQ=40;
double OLS_improve = 100*(OLS0-OLS)/OLS0;
double CAP_pct = 100*CAP/400;

$CAPTURE
PPS_conc_out HYD_conc_out CSA_conc_out AMI_conc_out
GAG PERM MC HIST SP NGF IL6 TNF C_FIBER SPINAL CENTRAL
CAP PAIN OLS FREQ OLS_improve CAP_pct
'

mod <- suppressMessages(mcode("ic_bps_shiny", ic_bps_code))

## ---------------------------------------------------------------------------
## Simulation helper
## ---------------------------------------------------------------------------
run_sim <- function(input_params, end_day = 365) {
  p <- input_params
  scenario <- p$scenario

  init_state <- c(
    GAG = p$gag_init / 100, PERM = 1 - p$gag_init / 100,
    MC = p$mc_init / 100, HIST = p$mc_init / 100 * 0.75,
    SP = p$pain_init / 10 * 0.8, NGF = p$pain_init / 10 * 0.75,
    IL6 = p$pain_init / 10 * 0.70, TNF = p$pain_init / 10 * 0.65,
    C_FIBER = p$pain_init / 10 * 0.85, SPINAL = p$pain_init / 10 * 0.75,
    CENTRAL = p$pain_init / 10 * 0.70,
    CAP = p$cap_init, PAIN = p$pain_init, OLS = p$ols_init,
    PPS_GUT = 0, PPS_CENT = 0, HYD_GUT = 0, HYD_CENT = 0,
    CSA_GUT = 0, CSA_CENT = 0, AMI_GUT = 0, AMI_CENT = 0
  )

  m <- init(mod, init_state)
  m <- param(m, list(OLS0 = p$ols_init, PAIN0 = p$pain_init,
                      Hunner = ifelse(p$hunner, 1.0, 0.0),
                      DMSO_eff = ifelse(scenario == "DMSO", 0.80, 0.0),
                      BoNTA_eff = ifelse(scenario == "BoNTA", 0.85, 0.0)))

  events <- switch(scenario,
    "PPS"  = ev(amt = 100, ii = 8/24, addl = end_day * 3, cmt = 1),
    "HYD"  = ev(amt = 25,  ii = 1,    addl = end_day,     cmt = 3),
    "CsA"  = ev(amt = 105, ii = 0.5,  addl = end_day * 2, cmt = 5),
    "AMI"  = ev(amt = 25,  ii = 1,    addl = end_day,     cmt = 7),
    "DMSO" = ev(amt = 50,  ii = 1,    addl = 60,          cmt = 1),
    "Triple" = ev(
      amt = c(100, 25, 25), ii = c(8/24, 1, 1),
      addl = c(end_day * 3, end_day, end_day), cmt = c(1, 3, 7)
    ),
    NULL
  )

  if (is.null(events)) {
    out <- mrgsim(m, end = end_day, delta = 1, obsonly = TRUE)
  } else {
    out <- mrgsim(m, events = events, end = end_day, delta = 1, obsonly = TRUE)
  }
  as.data.frame(out)
}

## ---------------------------------------------------------------------------
## UI
## ---------------------------------------------------------------------------
ui <- dashboardPage(
  skin = "blue",
  dashboardHeader(title = "IC/BPS QSP Dashboard", titleWidth = 280),

  dashboardSidebar(
    width = 260,
    sidebarMenu(
      menuItem("Overview",          tabName = "overview",   icon = icon("map")),
      menuItem("Patient Profile",   tabName = "patient",    icon = icon("user")),
      menuItem("PK — Drug Levels",  tabName = "pk",         icon = icon("pills")),
      menuItem("PD — Biomarkers",   tabName = "pd",         icon = icon("flask")),
      menuItem("Clinical Endpoints",tabName = "endpoints",  icon = icon("chart-line")),
      menuItem("Scenario Comparison",tabName = "scenarios", icon = icon("code-compare")),
      menuItem("Subtype Explorer",  tabName = "subtypes",   icon = icon("microscope")),
      menuItem("Sensitivity Analysis",tabName = "sensitivity",icon = icon("sliders"))
    ),

    hr(),
    h5("Patient Profile", style = "padding-left:15px; color:#ccc"),
    sliderInput("gag_init",   "GAG Layer Index (%)", 10, 90, 35, 5),
    sliderInput("mc_init",    "Mast Cell Activity (%)", 20, 100, 70, 5),
    sliderInput("cap_init",   "Bladder Capacity (mL)", 60, 350, 180, 10),
    sliderInput("pain_init",  "VAS Pain (0-10)", 1, 10, 6.5, 0.5),
    sliderInput("ols_init",   "OLS Score (0-20)", 2, 20, 14, 1),
    checkboxInput("hunner",   "Hunner Subtype", FALSE),

    hr(),
    h5("Treatment", style = "padding-left:15px; color:#ccc"),
    selectInput("scenario", "Treatment Scenario",
      choices = list(
        "No Treatment (Natural History)" = "None",
        "PPS (Elmiron) 100mg TID"        = "PPS",
        "Hydroxyzine 25mg QD"            = "HYD",
        "Cyclosporine A (Hunner)"        = "CsA",
        "Amitriptyline 25mg QD"          = "AMI",
        "DMSO Intravesical (6 sessions)" = "DMSO",
        "BoNTA 100U Intravesical"        = "BoNTA",
        "Triple Combo (PPS+HYD+AMI)"     = "Triple"
      ),
      selected = "None"
    ),
    sliderInput("sim_days", "Simulation Duration (days)", 30, 730, 365, 30),
    actionButton("run_btn", "Run Simulation", class = "btn-primary btn-block",
                  icon = icon("play")),

    hr(),
    div(style = "padding:10px; font-size:10px; color:#aaa",
        "IC/BPS QSP Model v1.0", br(),
        "CCR/Claude · 2026-06-26")
  ),

  dashboardBody(
    tabItems(

      ## ---- Tab 1: Overview ----
      tabItem("overview",
        fluidRow(
          box(width = 12, title = "Interstitial Cystitis / Bladder Pain Syndrome — QSP Overview",
              status = "primary", solidHeader = TRUE,
            fluidRow(
              column(6,
                h4("Disease Summary"),
                p(strong("IC/BPS"), "is a chronic, debilitating condition characterized by bladder pain,
                  urgency, and frequency, in the absence of infection or other identifiable causes.
                  It affects ~2-7% of women and ~0.5% of men, with significant quality-of-life impairment."),
                br(),
                h4("QSP Model Highlights"),
                tags$ul(
                  tags$li("22 ODE compartments: 8 PK + 14 PD"),
                  tags$li("7 treatment scenarios with mechanistic PK/PD"),
                  tags$li("Two subtypes: Hunner (inflammatory) vs Non-Hunner (barrier dysfunction)"),
                  tags$li("Endpoints: VAS pain, O'Leary-Sant score, bladder capacity, voiding frequency")
                ),
                br(),
                h4("Key Pathways"),
                tags$ul(
                  tags$li("GAG layer deficiency → urothelial permeability↑"),
                  tags$li("Mast cell activation → histamine, tryptase, PGE2"),
                  tags$li("C-fiber sensitization (TRPV1, P2X3, substance P)"),
                  tags$li("Spinal wind-up → central sensitization"),
                  tags$li("Neurogenic inflammation → fibrosis → bladder capacity↓")
                )
              ),
              column(6,
                h4("Mechanistic Map Preview"),
                tags$img(src = "https://raw.githubusercontent.com/pipetcpt/qsp/main/interstitial-cystitis/ic_bps_qsp_model.png",
                         style = "max-width:100%; border:1px solid #ddd; border-radius:4px",
                         alt = "IC/BPS Mechanistic Map"),
                br(), br(),
                h4("Subtype Comparison"),
                tableOutput("subtype_table")
              )
            )
          )
        ),
        fluidRow(
          valueBoxOutput("vbox_prev", width = 3),
          valueBoxOutput("vbox_gag",  width = 3),
          valueBoxOutput("vbox_mast", width = 3),
          valueBoxOutput("vbox_pain", width = 3)
        )
      ),

      ## ---- Tab 2: Patient Profile ----
      tabItem("patient",
        fluidRow(
          box(width = 6, title = "Patient Characteristics", status = "info", solidHeader = TRUE,
            h5("Baseline Biomarker Profile"),
            plotlyOutput("radar_plot", height = "350px"),
            hr(),
            tableOutput("patient_summary")
          ),
          box(width = 6, title = "UPOINT Domain Score", status = "warning", solidHeader = TRUE,
            h5("UPOINT Phenotyping System"),
            p("UPOINT (Urinary/Psychosocial/Organ-specific/Infection/Neurological/Tenderness)
               guides multimodal IC/BPS treatment."),
            plotlyOutput("upoint_bar", height = "300px"),
            hr(),
            h5("IC Subtype Classification"),
            verbatimTextOutput("subtype_class")
          )
        )
      ),

      ## ---- Tab 3: PK ----
      tabItem("pk",
        fluidRow(
          box(width = 12, title = "Drug Concentration-Time Profiles", status = "primary", solidHeader = TRUE,
            plotlyOutput("pk_plot", height = "500px")
          )
        ),
        fluidRow(
          box(width = 6, title = "PK Parameters Summary", status = "info",
            tableOutput("pk_params_table")
          ),
          box(width = 6, title = "Steady-State Drug Levels", status = "warning",
            tableOutput("pk_ss_table")
          )
        )
      ),

      ## ---- Tab 4: PD Biomarkers ----
      tabItem("pd",
        fluidRow(
          box(width = 12, title = "PD Biomarker Dynamics", status = "primary", solidHeader = TRUE,
            selectInput("pd_marker", "Select Biomarker",
              choices = c("GAG Layer" = "GAG", "Permeability" = "PERM",
                          "Mast Cell Activity" = "MC", "Histamine" = "HIST",
                          "Substance P" = "SP", "NGF" = "NGF",
                          "IL-6" = "IL6", "TNF-alpha" = "TNF",
                          "C-Fiber Sensitization" = "C_FIBER",
                          "Spinal Sensitization" = "SPINAL",
                          "Central Sensitization" = "CENTRAL"),
              selected = "GAG"
            ),
            plotlyOutput("pd_plot", height = "350px")
          )
        ),
        fluidRow(
          box(width = 12, title = "All Biomarkers Heatmap (Day 0 vs Day 365)", status = "info",
            plotlyOutput("biomarker_heatmap", height = "300px")
          )
        )
      ),

      ## ---- Tab 5: Clinical Endpoints ----
      tabItem("endpoints",
        fluidRow(
          box(width = 6, title = "VAS Pain Score", status = "danger", solidHeader = TRUE,
            plotlyOutput("pain_plot", height = "280px")
          ),
          box(width = 6, title = "O'Leary-Sant Score", status = "warning", solidHeader = TRUE,
            plotlyOutput("ols_plot", height = "280px")
          )
        ),
        fluidRow(
          box(width = 6, title = "Functional Bladder Capacity", status = "info", solidHeader = TRUE,
            plotlyOutput("cap_plot", height = "280px")
          ),
          box(width = 6, title = "Voiding Frequency", status = "success", solidHeader = TRUE,
            plotlyOutput("freq_plot", height = "280px")
          )
        ),
        fluidRow(
          box(width = 12, title = "Clinical Outcomes Summary Table", status = "primary",
            DTOutput("endpoints_table")
          )
        )
      ),

      ## ---- Tab 6: Scenario Comparison ----
      tabItem("scenarios",
        fluidRow(
          box(width = 12, title = "All Treatment Scenarios — 1-Year Simulation", status = "primary", solidHeader = TRUE,
            plotlyOutput("scenario_plot", height = "600px")
          )
        ),
        fluidRow(
          box(width = 12, title = "Scenario Comparison Table (Day 365)", status = "info",
            DTOutput("scenario_table")
          )
        )
      ),

      ## ---- Tab 7: Subtype Explorer ----
      tabItem("subtypes",
        fluidRow(
          box(width = 6, title = "Hunner vs Non-Hunner Comparison", status = "danger", solidHeader = TRUE,
            plotlyOutput("subtype_pain", height = "300px")
          ),
          box(width = 6, title = "Immunological Profile by Subtype", status = "warning", solidHeader = TRUE,
            plotlyOutput("subtype_immune", height = "300px")
          )
        ),
        fluidRow(
          box(width = 12, title = "CsA Response: Hunner Subtype vs Non-Hunner", status = "info",
            fluidRow(
              column(6, plotlyOutput("hunner_csa", height = "280px")),
              column(6,
                h5("Hunner Subtype Characteristics"),
                tags$ul(
                  tags$li("Prevalence: ~5-15% of all IC/BPS"),
                  tags$li("Cystoscopy: distinctive mucosal lesions"),
                  tags$li("Histology: dense lymphocyte/plasma cell infiltration"),
                  tags$li("IL-6 and IFN-γ significantly elevated"),
                  tags$li("Better response to CsA, corticosteroids"),
                  tags$li("Higher risk of bladder cancer progression")
                ),
                br(),
                h5("Non-Hunner (Majority)"),
                tags$ul(
                  tags$li("Prevalence: ~85-95% of IC/BPS"),
                  tags$li("Primary mechanism: GAG layer deficiency"),
                  tags$li("Prominent mast cell infiltration"),
                  tags$li("Better response to PPS, intravesical treatments"),
                  tags$li("TRPV1 upregulation drives urgency/pain")
                )
              )
            )
          )
        )
      ),

      ## ---- Tab 8: Sensitivity Analysis ----
      tabItem("sensitivity",
        fluidRow(
          box(width = 4, title = "Parameter Selection", status = "info",
            selectInput("sa_param", "Parameter",
              choices = c(
                "GAG Synthesis Rate" = "k_GAG_syn",
                "GAG Degradation Rate" = "k_GAG_deg",
                "Mast Cell Activation" = "k_MC_act",
                "C-fiber Sensitization" = "k_CFIB_up",
                "Central Sensitization" = "k_CENT_up",
                "PPS Bioavailability" = "F_PPS",
                "PPS Emax GAG" = "Emax_GAG_PPS"
              ), selected = "k_GAG_syn"
            ),
            sliderInput("sa_range", "Parameter Range (fold change)", 0.1, 5.0, c(0.5, 3.0), 0.1),
            selectInput("sa_endpoint", "Endpoint",
              choices = c("VAS Pain" = "PAIN", "OLS Score" = "OLS",
                          "GAG Index" = "GAG", "Bladder Cap (mL)" = "CAP"),
              selected = "PAIN"
            ),
            numericInput("sa_n", "Number of samples", 15, 5, 30, 1),
            actionButton("run_sa", "Run Sensitivity", class = "btn-warning btn-block")
          ),
          box(width = 8, title = "One-Way Sensitivity Analysis", status = "primary", solidHeader = TRUE,
            plotlyOutput("sa_plot", height = "400px")
          )
        ),
        fluidRow(
          box(width = 12, title = "Tornado Diagram — Day 365 OLS Score", status = "info",
            plotlyOutput("tornado_plot", height = "350px")
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

  ## Reactive simulation result
  sim_result <- eventReactive(input$run_btn, {
    withProgress(message = "Running IC/BPS simulation...", value = 0.5, {
      params <- list(
        gag_init  = input$gag_init,
        mc_init   = input$mc_init,
        cap_init  = input$cap_init,
        pain_init = input$pain_init,
        ols_init  = input$ols_init,
        hunner    = input$hunner,
        scenario  = input$scenario
      )
      run_sim(params, end_day = input$sim_days)
    })
  }, ignoreNULL = FALSE)

  ## ---- VALUE BOXES ----
  output$vbox_prev <- renderValueBox({
    valueBox("2-7%", "IC/BPS prevalence (women)", icon = icon("female"), color = "blue")
  })
  output$vbox_gag <- renderValueBox({
    valueBox(paste0(input$gag_init, "%"), "GAG Layer Index", icon = icon("shield-alt"), color = "green")
  })
  output$vbox_mast <- renderValueBox({
    valueBox(paste0(input$mc_init, "%"), "Mast Cell Activity", icon = icon("circle-exclamation"), color = "orange")
  })
  output$vbox_pain <- renderValueBox({
    valueBox(input$pain_init, "VAS Pain (0-10)", icon = icon("face-grimace"), color = "red")
  })

  ## ---- SUBTYPE TABLE ----
  output$subtype_table <- renderTable({
    data.frame(
      Feature = c("Prevalence", "Key finding", "IL-6/TNF", "Mast cells",
                  "GAG defect", "Best treatment"),
      Hunner  = c("5-15%", "Mucosal lesions", "Very high", "Moderate",
                  "Secondary", "CsA, steroids"),
      `Non-Hunner` = c("85-95%", "Petechiae/glomerulations", "Moderate", "High",
                        "Primary", "PPS, BoNTA, HYD"),
      check.names = FALSE
    )
  })

  ## ---- PATIENT SUMMARY ----
  output$patient_summary <- renderTable({
    data.frame(
      Parameter   = c("GAG Layer Index", "Mast Cell Activity", "Bladder Capacity",
                      "VAS Pain", "OLS Score", "Subtype"),
      Value       = c(paste0(input$gag_init, "%"), paste0(input$mc_init, "%"),
                      paste0(input$cap_init, " mL"),
                      input$pain_init, input$ols_init,
                      ifelse(input$hunner, "Hunner (inflammatory)", "Non-Hunner (barrier)")),
      Ref_Range   = c("70-100%", "0-30%", "300-500 mL", "0-10", "0-20", "—")
    )
  })

  ## ---- SUBTYPE CLASSIFICATION ----
  output$subtype_class <- renderText({
    lines <- c()
    if (input$hunner) {
      lines <- c(lines, "SUBTYPE: Hunner (Inflammatory)")
      lines <- c(lines, "  → Dense lymphocyte/plasma cell infiltration")
      lines <- c(lines, "  → Recommended: Cystoscopic fulguration + CsA or oral steroids")
    } else {
      lines <- c(lines, "SUBTYPE: Non-Hunner (Barrier Dysfunction)")
      if (input$mc_init > 70)
        lines <- c(lines, "  → High mast cell activity: Hydroxyzine, DMSO, intravesical heparin")
      if (input$gag_init < 40)
        lines <- c(lines, "  → Severe GAG deficiency: PPS (Elmiron), intravesical hyaluronic acid")
      if (input$pain_init > 7)
        lines <- c(lines, "  → Severe pain / central sensitization: Amitriptyline, SNM, CBT")
      if (input$cap_init < 150)
        lines <- c(lines, "  → Severely reduced capacity: Botulinum toxin A, hydrodistension")
    }
    paste(lines, collapse = "\n")
  })

  ## ---- PK PLOT ----
  output$pk_plot <- renderPlotly({
    df <- sim_result()
    req(df)

    pk_long <- df %>%
      select(time, PPS_conc_out, HYD_conc_out, CSA_conc_out, AMI_conc_out) %>%
      pivot_longer(-time, names_to = "Drug", values_to = "Conc") %>%
      mutate(Drug = recode(Drug,
        PPS_conc_out = "PPS (mg/L)",
        HYD_conc_out = "Hydroxyzine (ng/mL equiv)",
        CSA_conc_out = "Cyclosporine A (ng/mL)",
        AMI_conc_out = "Amitriptyline (ng/mL)"
      ))

    p <- ggplot(pk_long, aes(x = time, y = Conc, color = Drug)) +
      geom_line(linewidth = 0.8) +
      facet_wrap(~Drug, scales = "free_y", ncol = 2) +
      labs(x = "Time (days)", y = "Drug Concentration", title = "Drug Concentration-Time Profiles") +
      theme_bw(base_size = 11) + theme(legend.position = "none")
    ggplotly(p)
  })

  ## ---- PK PARAMS TABLE ----
  output$pk_params_table <- renderTable({
    data.frame(
      Drug = c("PPS (Elmiron)", "Hydroxyzine", "Cyclosporine A", "Amitriptyline"),
      F_pct = c("6%", "80%", "35%", "50%"),
      t_half = c("4h", "20h", "~24h", "~24h"),
      Vd = c("8 L", "350 L", "400 L", "1500 L"),
      CL = c("3.5 L/day", "35 L/day", "25 L/day", "110 L/day"),
      Target_conc = c("Urine GAG", "H1 block", "Trough 100-200 ng/mL", "Pain mod."),
      check.names = FALSE
    )
  })

  ## ---- PK SS TABLE ----
  output$pk_ss_table <- renderTable({
    df <- sim_result()
    req(df)
    last_row <- tail(df, 1)
    data.frame(
      Drug = c("PPS", "Hydroxyzine", "CsA", "Amitriptyline"),
      Css = c(
        round(last_row$PPS_conc_out, 3),
        round(last_row$HYD_conc_out, 2),
        round(last_row$CSA_conc_out, 1),
        round(last_row$AMI_conc_out, 2)
      ),
      Unit = c("mg/L", "ng/mL equiv", "ng/mL", "ng/mL"),
      check.names = FALSE
    )
  })

  ## ---- PD PLOT ----
  output$pd_plot <- renderPlotly({
    df <- sim_result()
    req(df)
    marker <- input$pd_marker
    ref_lines <- list(GAG = 0.8, PERM = 0.15, MC = 0.3, HIST = 0.1,
                       SP = 0.2, NGF = 0.2, IL6 = 0.15, TNF = 0.15,
                       C_FIBER = 0.3, SPINAL = 0.2, CENTRAL = 0.2)

    p <- ggplot(df, aes_string(x = "time", y = marker)) +
      geom_line(color = "#1976D2", linewidth = 1.2) +
      geom_hline(yintercept = ref_lines[[marker]], linetype = "dashed", color = "red", alpha = 0.6) +
      annotate("text", x = max(df$time) * 0.7, y = ref_lines[[marker]] * 1.05,
               label = "Reference (healthy)", color = "red", size = 3.5) +
      labs(x = "Time (days)", y = marker,
           title = paste("PD Biomarker:", marker, "over time")) +
      theme_bw(base_size = 12)
    ggplotly(p)
  })

  ## ---- BIOMARKER HEATMAP ----
  output$biomarker_heatmap <- renderPlotly({
    df <- sim_result()
    req(df)
    markers <- c("GAG", "PERM", "MC", "HIST", "SP", "NGF", "IL6", "TNF", "C_FIBER", "SPINAL", "CENTRAL")
    d0 <- df %>% filter(time == min(time)) %>% select(all_of(markers)) %>% as.numeric()
    d365 <- df %>% filter(time == max(time)) %>% select(all_of(markers)) %>% as.numeric()

    heat_df <- data.frame(
      Biomarker = rep(markers, 2),
      Timepoint = rep(c("Baseline", "End of Treatment"), each = length(markers)),
      Value = c(d0, d365)
    )

    p <- ggplot(heat_df, aes(x = Timepoint, y = Biomarker, fill = Value)) +
      geom_tile(color = "white") +
      scale_fill_gradient2(low = "#1565C0", mid = "#FFFDE7", high = "#B71C1C",
                            midpoint = 0.5, limit = c(0, 1.5)) +
      labs(title = "Biomarker Heatmap: Baseline vs End of Treatment", fill = "Level") +
      theme_minimal(base_size = 11)
    ggplotly(p)
  })

  ## ---- PAIN PLOT ----
  output$pain_plot <- renderPlotly({
    df <- sim_result(); req(df)
    p <- ggplot(df, aes(x = time, y = PAIN)) +
      geom_line(color = "#C62828", linewidth = 1.2) +
      geom_hline(yintercept = 3, linetype = "dashed", color = "green") +
      scale_y_continuous(limits = c(0, 10)) +
      labs(x = "Days", y = "VAS Pain (0-10)", title = "VAS Pain Score") +
      theme_bw(base_size = 11)
    ggplotly(p)
  })

  ## ---- OLS PLOT ----
  output$ols_plot <- renderPlotly({
    df <- sim_result(); req(df)
    p <- ggplot(df, aes(x = time, y = OLS)) +
      geom_line(color = "#E65100", linewidth = 1.2) +
      geom_hline(yintercept = 8, linetype = "dashed", color = "green") +
      scale_y_continuous(limits = c(0, 20)) +
      labs(x = "Days", y = "OLS Score (0-20)", title = "O'Leary-Sant Symptom Score") +
      theme_bw(base_size = 11)
    ggplotly(p)
  })

  ## ---- CAPACITY PLOT ----
  output$cap_plot <- renderPlotly({
    df <- sim_result(); req(df)
    p <- ggplot(df, aes(x = time, y = CAP)) +
      geom_line(color = "#1565C0", linewidth = 1.2) +
      geom_hline(yintercept = 350, linetype = "dashed", color = "green") +
      labs(x = "Days", y = "Capacity (mL)", title = "Functional Bladder Capacity") +
      theme_bw(base_size = 11)
    ggplotly(p)
  })

  ## ---- FREQUENCY PLOT ----
  output$freq_plot <- renderPlotly({
    df <- sim_result(); req(df)
    p <- ggplot(df, aes(x = time, y = FREQ)) +
      geom_line(color = "#2E7D32", linewidth = 1.2) +
      geom_hline(yintercept = 8, linetype = "dashed", color = "green") +
      labs(x = "Days", y = "Voids/24h", title = "Voiding Frequency") +
      theme_bw(base_size = 11)
    ggplotly(p)
  })

  ## ---- ENDPOINTS TABLE ----
  output$endpoints_table <- renderDT({
    df <- sim_result(); req(df)
    df %>%
      filter(time %in% c(0, 30, 90, 180, 270, 365)) %>%
      select(time, PAIN, OLS, CAP, FREQ, GAG, MC, C_FIBER, CENTRAL, OLS_improve) %>%
      mutate(across(where(is.numeric), round, 2)) %>%
      rename(Day = time, `VAS Pain` = PAIN, `OLS Score` = OLS,
             `Cap (mL)` = CAP, `Freq/24h` = FREQ, `GAG Index` = GAG,
             `Mast Cell` = MC, `C-Fiber` = C_FIBER, `Central Sens.` = CENTRAL,
             `OLS Improve %` = OLS_improve) %>%
      datatable(options = list(pageLength = 10, scrollX = TRUE),
                class = "compact stripe hover")
  })

  ## ---- ALL SCENARIOS ----
  all_scenarios <- reactive({
    withProgress(message = "Running all 7 scenarios...", value = 0, {
      base_init <- list(gag_init = input$gag_init, mc_init = input$mc_init,
                         cap_init = input$cap_init, pain_init = input$pain_init,
                         ols_init = input$ols_init, hunner = input$hunner)
      scenarios <- c("None", "PPS", "HYD", "CsA", "AMI", "DMSO", "Triple")
      labels    <- c("S1: No Tx", "S2: PPS", "S3: HYD", "S4: CsA", "S5: AMI",
                     "S6: DMSO", "S7: Triple")
      bind_rows(lapply(seq_along(scenarios), function(i) {
        setProgress(i / length(scenarios))
        p <- c(base_init, list(scenario = scenarios[i]))
        run_sim(p, end_day = input$sim_days) %>% mutate(Scenario = labels[i])
      }))
    })
  })

  output$scenario_plot <- renderPlotly({
    df <- all_scenarios()
    req(df)
    colors <- c("#E53935","#1E88E5","#43A047","#8E24AA","#FB8C00","#00897B","#D81B60")
    p <- df %>%
      select(time, PAIN, OLS, CAP, GAG, MC, Scenario) %>%
      pivot_longer(c(PAIN, OLS, CAP, GAG, MC), names_to = "Endpoint", values_to = "Value") %>%
      ggplot(aes(x = time, y = Value, color = Scenario)) +
      geom_line(linewidth = 0.7) +
      facet_wrap(~Endpoint, scales = "free_y", ncol = 3) +
      scale_color_manual(values = colors) +
      labs(x = "Days", y = "Value", title = "IC/BPS — All Treatment Scenarios") +
      theme_bw(base_size = 10) +
      theme(legend.position = "bottom", legend.text = element_text(size = 8))
    ggplotly(p) %>% layout(legend = list(orientation = "h", y = -0.15))
  })

  output$scenario_table <- renderDT({
    df <- all_scenarios()
    req(df)
    df %>%
      filter(time == max(time)) %>%
      select(Scenario, PAIN, OLS, CAP, GAG, MC, CENTRAL, OLS_improve) %>%
      mutate(across(where(is.numeric), round, 2)) %>%
      rename(`VAS Pain` = PAIN, `OLS Score` = OLS, `Cap (mL)` = CAP,
             `GAG Index` = GAG, `Mast Cell` = MC, `Central Sens.` = CENTRAL,
             `OLS Improve %` = OLS_improve) %>%
      datatable(options = list(pageLength = 10, scrollX = TRUE),
                class = "compact stripe hover")
  })

  ## ---- SUBTYPE COMPARISONS ----
  output$subtype_pain <- renderPlotly({
    df_nonH <- run_sim(list(gag_init = input$gag_init, mc_init = input$mc_init,
                             cap_init = input$cap_init, pain_init = input$pain_init,
                             ols_init = input$ols_init, hunner = FALSE,
                             scenario = "PPS"), end_day = 365) %>%
      mutate(Subtype = "Non-Hunner (PPS)")
    df_H <- run_sim(list(gag_init = input$gag_init, mc_init = input$mc_init,
                          cap_init = input$cap_init, pain_init = input$pain_init,
                          ols_init = input$ols_init, hunner = TRUE,
                          scenario = "CsA"), end_day = 365) %>%
      mutate(Subtype = "Hunner (CsA)")
    df <- bind_rows(df_nonH, df_H)
    p <- ggplot(df, aes(x = time, y = PAIN, color = Subtype)) +
      geom_line(linewidth = 1.1) +
      scale_color_manual(values = c("Non-Hunner (PPS)" = "#1976D2", "Hunner (CsA)" = "#C62828")) +
      labs(x = "Days", y = "VAS Pain", title = "Pain Response by Subtype + Matched Therapy") +
      theme_bw(base_size = 11)
    ggplotly(p)
  })

  output$subtype_immune <- renderPlotly({
    df_nonH <- run_sim(list(gag_init = 35, mc_init = 75, cap_init = 180,
                             pain_init = 6.5, ols_init = 14, hunner = FALSE,
                             scenario = "None"), end_day = 365) %>%
      mutate(Subtype = "Non-Hunner")
    df_H <- run_sim(list(gag_init = 35, mc_init = 75, cap_init = 180,
                          pain_init = 7.5, ols_init = 16, hunner = TRUE,
                          scenario = "None"), end_day = 365) %>%
      mutate(Subtype = "Hunner")
    df <- bind_rows(df_nonH, df_H)
    p <- ggplot(df, aes(x = time, y = IL6, color = Subtype)) +
      geom_line(linewidth = 1.1) +
      scale_color_manual(values = c("Non-Hunner" = "#1976D2", "Hunner" = "#C62828")) +
      labs(x = "Days", y = "IL-6 Index", title = "IL-6 Dynamics (No Treatment) by Subtype") +
      theme_bw(base_size = 11)
    ggplotly(p)
  })

  output$hunner_csa <- renderPlotly({
    df_H_CsA <- run_sim(list(gag_init = 30, mc_init = 60, cap_init = 160,
                              pain_init = 7.5, ols_init = 16, hunner = TRUE,
                              scenario = "CsA"), end_day = 365) %>%
      mutate(Group = "Hunner + CsA")
    df_H_none <- run_sim(list(gag_init = 30, mc_init = 60, cap_init = 160,
                               pain_init = 7.5, ols_init = 16, hunner = TRUE,
                               scenario = "None"), end_day = 365) %>%
      mutate(Group = "Hunner + No Tx")
    df_nH_CsA <- run_sim(list(gag_init = 30, mc_init = 60, cap_init = 160,
                               pain_init = 7.5, ols_init = 16, hunner = FALSE,
                               scenario = "CsA"), end_day = 365) %>%
      mutate(Group = "Non-Hunner + CsA")
    df <- bind_rows(df_H_CsA, df_H_none, df_nH_CsA)
    p <- ggplot(df, aes(x = time, y = OLS, color = Group)) +
      geom_line(linewidth = 1.1) +
      scale_color_manual(values = c("Hunner + CsA" = "#8E24AA",
                                     "Hunner + No Tx" = "#C62828",
                                     "Non-Hunner + CsA" = "#1976D2")) +
      labs(x = "Days", y = "OLS Score", title = "CsA Efficacy by Subtype") +
      theme_bw(base_size = 11)
    ggplotly(p)
  })

  ## ---- SENSITIVITY ----
  sa_result <- eventReactive(input$run_sa, {
    base_params <- list(gag_init = input$gag_init, mc_init = input$mc_init,
                         cap_init = input$cap_init, pain_init = input$pain_init,
                         ols_init = input$ols_init, hunner = input$hunner,
                         scenario = input$scenario)

    param_name <- input$sa_param
    baseline_val <- as.numeric(param(mod)[param_name])
    fold_range <- seq(input$sa_range[1], input$sa_range[2], length.out = input$sa_n)
    param_vals <- fold_range * baseline_val
    endpoint   <- input$sa_endpoint

    results_list <- lapply(param_vals, function(pv) {
      init_state <- c(GAG = base_params$gag_init / 100,
                       PERM = 1 - base_params$gag_init / 100,
                       MC = base_params$mc_init / 100,
                       HIST = base_params$mc_init / 100 * 0.75,
                       SP = base_params$pain_init / 10 * 0.8,
                       NGF = base_params$pain_init / 10 * 0.75,
                       IL6 = base_params$pain_init / 10 * 0.70,
                       TNF = base_params$pain_init / 10 * 0.65,
                       C_FIBER = base_params$pain_init / 10 * 0.85,
                       SPINAL = base_params$pain_init / 10 * 0.75,
                       CENTRAL = base_params$pain_init / 10 * 0.70,
                       CAP = base_params$cap_init, PAIN = base_params$pain_init,
                       OLS = base_params$ols_init,
                       PPS_GUT = 0, PPS_CENT = 0, HYD_GUT = 0, HYD_CENT = 0,
                       CSA_GUT = 0, CSA_CENT = 0, AMI_GUT = 0, AMI_CENT = 0)

      m_sa <- param(mod, setNames(list(pv), param_name))
      m_sa <- init(m_sa, init_state)
      m_sa <- param(m_sa, list(OLS0 = base_params$ols_init, PAIN0 = base_params$pain_init))

      tryCatch({
        out <- mrgsim(m_sa, end = 365, delta = 7, obsonly = TRUE) %>% as.data.frame()
        val <- out %>% filter(time == max(time)) %>% pull(endpoint) %>% mean()
        data.frame(param_val = pv, fold = pv / baseline_val, endpoint_val = val)
      }, error = function(e) data.frame(param_val = pv, fold = pv / baseline_val, endpoint_val = NA))
    })

    bind_rows(results_list)
  })

  output$sa_plot <- renderPlotly({
    df <- sa_result()
    req(df)
    p <- ggplot(df %>% filter(!is.na(endpoint_val)),
                aes(x = fold, y = endpoint_val)) +
      geom_line(color = "#1565C0", linewidth = 1.2) +
      geom_point(color = "#1565C0", size = 2.5) +
      geom_vline(xintercept = 1, linetype = "dashed", color = "red") +
      labs(x = paste("Fold Change:", input$sa_param),
           y = paste("Endpoint at Day 365:", input$sa_endpoint),
           title = "One-Way Sensitivity Analysis") +
      theme_bw(base_size = 12)
    ggplotly(p)
  })

  output$tornado_plot <- renderPlotly({
    # Pre-computed tornado for key parameters
    params_list <- c("k_GAG_syn", "k_GAG_deg", "k_MC_act", "k_CFIB_up",
                      "k_CENT_up", "F_PPS", "Emax_GAG_PPS")
    param_labels <- c("GAG Synthesis", "GAG Degradation", "Mast Cell Act.",
                       "C-fiber Sensitization", "Central Sensitization",
                       "PPS Bioavailability", "PPS Emax GAG")

    init_state <- c(GAG = 0.35, PERM = 0.65, MC = 0.70, HIST = 0.55,
                     SP = 0.60, NGF = 0.55, IL6 = 0.50, TNF = 0.45,
                     C_FIBER = 0.65, SPINAL = 0.55, CENTRAL = 0.50,
                     CAP = 180, PAIN = 6.5, OLS = 14,
                     PPS_GUT = 0, PPS_CENT = 0, HYD_GUT = 0, HYD_CENT = 0,
                     CSA_GUT = 0, CSA_CENT = 0, AMI_GUT = 0, AMI_CENT = 0)

    tornado_vals <- lapply(params_list, function(pname) {
      base_val <- as.numeric(param(mod)[pname])
      e_pps <- ev(amt = 100, ii = 8/24, addl = 365 * 3, cmt = 1)

      get_ols <- function(fold) {
        tryCatch({
          m2 <- param(mod, setNames(list(base_val * fold), pname))
          m2 <- init(m2, init_state)
          m2 <- param(m2, list(OLS0 = 14, PAIN0 = 6.5))
          out <- mrgsim(m2, events = e_pps, end = 365, delta = 30, obsonly = TRUE)
          tail(as.data.frame(out)$OLS, 1)
        }, error = function(e) 14)
      }

      data.frame(
        Param = pname, Label = param_labels[match(pname, params_list)],
        Low  = get_ols(0.5),
        High = get_ols(2.0)
      )
    })

    td <- bind_rows(tornado_vals)
    base_ols <- 14
    td <- td %>%
      mutate(delta_low = Low - base_ols, delta_high = High - base_ols,
             range = abs(delta_high - delta_low)) %>%
      arrange(range)

    p <- plot_ly() %>%
      add_bars(data = td, y = ~Label, x = ~(Low - base_ols), orientation = "h",
               name = "Low (×0.5)", marker = list(color = "#1976D2")) %>%
      add_bars(data = td, y = ~Label, x = ~(High - base_ols), orientation = "h",
               name = "High (×2.0)", marker = list(color = "#C62828")) %>%
      layout(
        title = "Tornado Diagram — OLS Score at Day 365 (PPS scenario, deviation from base=14)",
        xaxis = list(title = "ΔOLS Score vs Baseline"),
        yaxis = list(title = ""),
        barmode = "overlay",
        shapes = list(list(type = "line", x0 = 0, x1 = 0, y0 = -0.5, y1 = nrow(td) - 0.5,
                            line = list(color = "black", width = 2, dash = "dash")))
      )
    p
  })

  ## ---- RADAR PLOT (Patient Profile) ----
  output$radar_plot <- renderPlotly({
    categories <- c("GAG Layer", "Mast Cell Suppress.",
                     "Bladder Cap.", "Pain Control", "OLS Control")
    values_pts <- c(input$gag_init / 100,
                     1 - input$mc_init / 100,
                     input$cap_init / 400,
                     1 - input$pain_init / 10,
                     1 - input$ols_init / 20)
    values_norm <- rep(1, 5)

    p <- plot_ly(type = "scatterpolar", fill = "toself") %>%
      add_trace(r = values_pts, theta = categories, name = "Patient",
                line = list(color = "#C62828"), fillcolor = "rgba(198,40,40,0.2)") %>%
      add_trace(r = values_norm, theta = categories, name = "Normal",
                line = list(color = "#1976D2", dash = "dash"), fillcolor = "rgba(25,118,210,0.1)") %>%
      layout(polar = list(radialaxis = list(visible = TRUE, range = c(0, 1))),
             showlegend = TRUE, title = "Patient Baseline Profile (1 = Normal)")
    p
  })

  ## ---- UPOINT BAR ----
  output$upoint_bar <- renderPlotly({
    scores <- data.frame(
      Domain = c("U: Urinary", "P: Psychosocial", "O: Organ", "I: Infection",
                  "N: Neurological", "T: Tenderness (pelvic floor)"),
      Score = c(
        min(input$ols_init / 2, 3),
        min((10 - input$pain_init) * 0.3, 3),
        min(1 - input$gag_init / 100, 1) * 3,
        0.5,
        min(input$pain_init / 10 * 3, 3),
        min(1 - input$cap_init / 400, 1) * 3
      )
    )
    p <- plot_ly(scores, x = ~Score, y = ~Domain, type = "bar", orientation = "h",
                  marker = list(color = "#1976D2")) %>%
      layout(title = "UPOINT Domain Scores (estimated)", xaxis = list(range = c(0, 3)),
             yaxis = list(title = ""))
    p
  })
}

## ---------------------------------------------------------------------------
## Launch
## ---------------------------------------------------------------------------
shinyApp(ui, server)
