## =============================================================================
## Polycythemia Vera QSP — Interactive Shiny Dashboard
## =============================================================================
## Tabs: (1) Overview  (2) Patient Profile  (3) PK Parameters
##        (4) PD & Hematologic Response  (5) Clinical Endpoints
##        (6) Scenario Comparison  (7) Biomarkers & Disease Progression
## =============================================================================

library(shiny)
library(shinydashboard)
library(mrgsolve)
library(dplyr)
library(ggplot2)
library(plotly)
library(DT)
library(patchwork)
library(scales)

## ---------------------------------------------------------------------------
## Embed mrgsolve model code inline
## ---------------------------------------------------------------------------
pv_code <- '
$PROB Polycythemia Vera QSP Model (Shiny)

$PARAM @annotated
// Disease parameters
JAK2_AB0  : 50    : JAK2 V617F allele burden baseline (%)
Emax_JAK2 : 1.5   : JAK2-driven proliferative fold-increase
k_clone   : 0.003 : Clone expansion rate (/day)
k_clone_d : 0.001 : Clone regression rate (/day)
// Kinetics
k_BFUE_in : 0.03  : BFU-E input rate
k_BFUE_d  : 0.05  : BFU-E differentiation rate
k_CFUE_d  : 0.10  : CFU-E differentiation rate
k_RETIC_BM: 0.15  : BM retic maturation rate
k_RETIC_C : 0.50  : Circ retic maturation rate
k_RBC_d   : 0.0083: RBC elimination
k_PLT_in  : 0.12  : Platelet production
k_PLT_d   : 0.10  : Platelet elimination
k_WBC_in  : 0.33  : WBC production
k_WBC_d   : 0.33  : WBC elimination
PLT_0     : 400   : Platelet baseline
WBC_0     : 8     : WBC baseline
BFUE_0    : 100   : BFU-E baseline
SPL_0     : 450   : Spleen baseline (mL)
k_SPL_in  : 0.002 : Spleen expansion rate
k_SPL_d   : 0.005 : Spleen regression rate
SPL_max   : 3000  : Max spleen volume
k_FIB_in  : 0.0005: Fibrosis rate
k_FIB_d   : 0.0002: Fibrosis regression
// Ruxolitinib PK
KA_RUX    : 2.4   : ka (/h)
F_RUX     : 0.95  : Bioavailability
VC_RUX    : 72    : Central volume (L)
VP_RUX    : 28    : Peripheral volume (L)
CL_RUX    : 22    : Clearance (L/h)
QP_RUX    : 16    : Inter-compartmental CL (L/h)
IC50_RUX  : 0.86  : IC50 (ng/mL, JAK2)
Imax_RUX  : 1.0   : Max inhibition
HILL_RUX  : 1.0   : Hill coefficient
// HYU PK/PD
KA_HYU    : 1.5   : HYU ka (/h)
F_HYU     : 0.80  : HYU bioavailability
VC_HYU    : 20    : HYU central volume
CL_HYU    : 5.0   : HYU clearance
IC50_HYU  : 150   : HYU IC50 (μM)
Imax_HYU  : 0.85  : HYU max inhibition
// IFN PK/PD
KA_IFN    : 0.03  : IFN ka (/h)
F_IFN     : 0.84  : IFN SC bioavailability
VC_IFN    : 6.0   : IFN central volume
CL_IFN    : 0.07  : IFN clearance
EC50_IFN  : 50    : IFN EC50 (pg/mL)
Emax_IFN  : 0.70  : IFN max clone suppression
// Treatment switches
USE_RUX   : 0     : Ruxolitinib (0/1)
USE_HYU   : 0     : Hydroxyurea (0/1)
USE_IFN   : 0     : PEG-IFN-a2a (0/1)
USE_ASP   : 1     : Aspirin (0/1)

$CMT DEPOT_RUX CENT_RUX PERI_RUX CENT_HYU SC_IFN CENT_IFN
     BFUE CFUE RETIC_BM RETIC_C RBC PLT WBC SPL FIBRO ALLELE

$INIT DEPOT_RUX=0 CENT_RUX=0 PERI_RUX=0 CENT_HYU=0 SC_IFN=0 CENT_IFN=0
      BFUE=100 CFUE=100 RETIC_BM=100 RETIC_C=100
      RBC=130 PLT=600 WBC=14 SPL=800 FIBRO=0.5 ALLELE=50

$ODE
double CP_RUX   = CENT_RUX / VC_RUX * 1000.0;
double CP_HYU_uM = (CENT_HYU / VC_HYU * 1000.0) * 13.2;
double CP_IFN   = CENT_IFN / VC_IFN * 1e6;

double INH_RUX = USE_RUX * Imax_RUX * pow(CP_RUX, HILL_RUX) /
                 (pow(IC50_RUX, HILL_RUX) + pow(CP_RUX, HILL_RUX));
double INH_HYU = USE_HYU * Imax_HYU * CP_HYU_uM / (IC50_HYU + CP_HYU_uM);
double EFF_IFN = USE_IFN * Emax_IFN * CP_IFN / (EC50_IFN + CP_IFN);
double INH_PROLIF = 1.0 - fmax(INH_RUX, INH_HYU);

double JAK2_effect = 1.0 + (Emax_JAK2 - 1.0) * (ALLELE / 100.0);

dxdt_DEPOT_RUX = -KA_RUX * DEPOT_RUX;
dxdt_CENT_RUX  = KA_RUX * DEPOT_RUX * F_RUX
                - (CL_RUX + QP_RUX) * (CENT_RUX / VC_RUX)
                + QP_RUX * (PERI_RUX / VP_RUX);
dxdt_PERI_RUX  = QP_RUX * (CENT_RUX / VC_RUX - PERI_RUX / VP_RUX);
dxdt_CENT_HYU  = KA_HYU * F_HYU * 500.0 / 24.0 * USE_HYU
                - CL_HYU * (CENT_HYU / VC_HYU);
dxdt_SC_IFN    = -KA_IFN * SC_IFN;
dxdt_CENT_IFN  = KA_IFN * F_IFN * SC_IFN - CL_IFN * (CENT_IFN / VC_IFN);

dxdt_BFUE     = k_BFUE_in * JAK2_effect * INH_PROLIF * BFUE_0
              - k_BFUE_d * BFUE;
dxdt_CFUE     = k_BFUE_d * BFUE - k_CFUE_d * CFUE;
dxdt_RETIC_BM = k_CFUE_d * CFUE - k_RETIC_BM * RETIC_BM;
dxdt_RETIC_C  = k_RETIC_BM * RETIC_BM - k_RETIC_C * RETIC_C;
dxdt_RBC      = k_RETIC_C * RETIC_C - k_RBC_d * RBC;
dxdt_PLT      = k_PLT_in * JAK2_effect * INH_PROLIF * PLT_0 - k_PLT_d * PLT;
dxdt_WBC      = k_WBC_in * JAK2_effect * INH_PROLIF * WBC_0 - k_WBC_d * WBC;

double BM_overflow = fmax(0.0, (RBC - 110.0) / 110.0);
dxdt_SPL  = k_SPL_in * BM_overflow * (SPL_max - SPL) / SPL_max
           - k_SPL_d * SPL * (INH_RUX + INH_HYU + 0.001);

dxdt_FIBRO = k_FIB_in * (ALLELE / 50.0) * (SPL / 800.0)
            - k_FIB_d * FIBRO * (INH_RUX + EFF_IFN + 0.001);
if (FIBRO >= 3.0) dxdt_FIBRO = fmin(0.0, dxdt_FIBRO);
if (FIBRO <= 0.0) dxdt_FIBRO = fmax(0.0, dxdt_FIBRO);

double k_allele_expansion = k_clone * (1.0 - ALLELE / 100.0);
dxdt_ALLELE = k_allele_expansion * ALLELE - (k_clone_d + EFF_IFN * 0.01) * ALLELE;
if (ALLELE >= 100) dxdt_ALLELE = fmin(0.0, dxdt_ALLELE);
if (ALLELE <= 0)   dxdt_ALLELE = fmax(0.0, dxdt_ALLELE);

$TABLE
double Hct = RBC * 0.35;
double Hgb = RBC * 0.115;
double EPO = 6.0 * exp(-0.03 * (RBC - 100.0));
double CP_RUX_out = CENT_RUX / VC_RUX * 1000.0;
double CP_HYU_out = (CENT_HYU / VC_HYU * 1000.0) * 13.2;
double CP_IFN_out = CENT_IFN / VC_IFN * 1e6;
double pSTAT5 = USE_RUX * Imax_RUX * CP_RUX_out /
                (IC50_RUX + CP_RUX_out);
double SVR = (800.0 - SPL) / 800.0 * 100.0;
double CHR = (Hct < 45.0 && PLT < 400.0 && WBC < 10.0 && SPL < 450.0) ? 1.0 : 0.0;
double THROMB_RISK = (Hct > 45 ? 0.06 : 0.02) + (PLT > 1000 ? 0.04 : 0.0);
double MPN_SAF = fmax(0, 20.0 + (Hct-45.0)*2.0 + (PLT-400.0)/50.0 +
                      (SPL-450.0)/100.0 - pSTAT5*30.0);
double MF_risk = 0.005 + FIBRO * 0.01 + (ALLELE > 80 ? 0.01 : 0.0);

capture Hct; capture Hgb; capture EPO;
capture PLT_out = PLT; capture WBC_out = WBC;
capture SPL_vol = SPL; capture SVR; capture CHR;
capture ALLELE_out = ALLELE; capture FIBRO_out = FIBRO;
capture pSTAT5; capture THROMB_RISK; capture MPN_SAF; capture MF_risk;
capture CP_RUX_out; capture CP_HYU_out; capture CP_IFN_out;
'

mod <- mcode("pv_shiny", pv_code)

## ---------------------------------------------------------------------------
## Helper: run simulation for a scenario
## ---------------------------------------------------------------------------
run_sim <- function(mod, params, rux_dose=10, hyu_dose=500, ifn_dose=45,
                    sim_weeks=48, delta=1) {
  evs <- NULL
  if (params$USE_RUX == 1) {
    evs <- ev(cmt="DEPOT_RUX", amt=rux_dose, ii=12,
              addl=sim_weeks*14-1, time=0)
  }
  if (params$USE_IFN == 1) {
    e2 <- ev(cmt="SC_IFN", amt=ifn_dose, ii=168,
             addl=sim_weeks-1, time=0)
    evs <- if(is.null(evs)) e2 else ev_seq(evs, e2)
  }
  sim <- mod %>%
    param(params) %>%
    {if(!is.null(evs)) mrgsim(., ev=evs, end=sim_weeks*7, delta=delta)
     else mrgsim(., end=sim_weeks*7, delta=delta)} %>%
    as.data.frame() %>%
    mutate(time_wk = time / 7)
  return(sim)
}

## ---------------------------------------------------------------------------
## UI
## ---------------------------------------------------------------------------
ui <- dashboardPage(
  skin = "red",
  dashboardHeader(title = "PV QSP Dashboard",
                  titleWidth = 280),
  dashboardSidebar(
    width = 280,
    sidebarMenu(
      menuItem("Overview",          tabName="tab_overview",   icon=icon("home")),
      menuItem("Patient Profile",   tabName="tab_patient",    icon=icon("user-md")),
      menuItem("Pharmacokinetics",  tabName="tab_pk",         icon=icon("flask")),
      menuItem("PD & Hematology",   tabName="tab_pd",         icon=icon("tint")),
      menuItem("Clinical Endpoints",tabName="tab_clinical",   icon=icon("stethoscope")),
      menuItem("Scenario Comparison",tabName="tab_compare",   icon=icon("chart-bar")),
      menuItem("Biomarkers",        tabName="tab_biomarker",  icon=icon("dna"))
    )
  ),
  dashboardBody(
    tags$head(tags$style(HTML("
      .main-header .logo { font-weight: bold; font-size: 16px; }
      .box-title { font-weight: bold; }
      .info-box { min-height: 65px; }
    "))),
    tabItems(

      ## ====== TAB 1: OVERVIEW =============================================
      tabItem(tabName="tab_overview",
        fluidRow(
          box(width=12, title="Polycythemia Vera QSP Model",
              solidHeader=TRUE, status="danger",
              h4("진성 다혈증 (Polycythemia Vera) — Quantitative Systems Pharmacology Model"),
              p("This dashboard implements a mechanistic QSP model for Polycythemia Vera (PV),
                 a BCR-ABL-negative myeloproliferative neoplasm driven by the JAK2 V617F gain-of-function mutation."),
              hr(),
              tags$ul(
                tags$li(strong("Disease Driver:"), " JAK2 V617F mutation → constitutive JAK-STAT5 signaling → uncontrolled erythropoiesis"),
                tags$li(strong("Key Drugs:"), " Ruxolitinib (JAK1/2 inhibitor), Hydroxyurea (cytoreductive), PEG-IFN-α2a (clonal suppression)"),
                tags$li(strong("Clinical Endpoints:"), " Hematocrit <45%, Platelet normalization, Spleen volume reduction (SVR35), Allele burden"),
                tags$li(strong("Model Compartments:"), " 16 ODEs covering erythropoiesis, thrombopoiesis, WBC, spleen, fibrosis, allele burden"),
                tags$li(strong("Key Trials Calibrated:"), " RESPONSE (ruxolitinib), PROUD-PV (IFN), REVEAL (natural history)")
              ),
              hr(),
              fluidRow(
                valueBox(16, "ODE Compartments", icon=icon("cogs"), color="red", width=3),
                valueBox(4, "Drug Scenarios", icon=icon("pills"), color="orange", width=3),
                valueBox(6, "Dashboard Tabs", icon=icon("chart-line"), color="blue", width=3),
                valueBox("JAK2 V617F", "Driver Mutation", icon=icon("dna"), color="purple", width=3)
              )
          )
        ),
        fluidRow(
          box(width=6, title="Disease Overview",
              img(src="https://upload.wikimedia.org/wikipedia/commons/thumb/a/a5/Polycythemia_vera_-_high_mag.jpg/640px-Polycythemia_vera_-_high_mag.jpg",
                  style="max-width:100%"),
              p(em("Hypercellular bone marrow in PV (H&E stain)"))
          ),
          box(width=6, title="Key Pathophysiology",
              plotOutput("overview_path_plot", height="250px"))
        )
      ),

      ## ====== TAB 2: PATIENT PROFILE ======================================
      tabItem(tabName="tab_patient",
        fluidRow(
          box(width=4, title="Patient Characteristics", status="primary", solidHeader=TRUE,
            sliderInput("AGE",  "Age (years)", 30, 85, 65, step=1),
            radioButtons("SEX", "Sex", choices=list("Male"=1,"Female"=0), inline=TRUE),
            sliderInput("JAK2_AB0", "JAK2 V617F Allele Burden (%)", 10, 95, 50, step=5),
            sliderInput("Hct_base", "Baseline Hematocrit (%)", 45, 75, 55, step=1),
            sliderInput("PLT_base", "Baseline Platelets (×10⁹/L)", 200, 1500, 600, step=50),
            sliderInput("SPL_base", "Baseline Spleen Volume (mL)", 300, 2000, 800, step=50)
          ),
          box(width=4, title="Risk Stratification", status="warning", solidHeader=TRUE,
            radioButtons("HIST_THROMB", "History of Thrombosis",
                         choices=list("No"=0,"Yes"=1), inline=TRUE),
            radioButtons("WBC_risk", "WBC >11×10⁹/L at diagnosis",
                         choices=list("No"=0,"Yes"=1), inline=TRUE),
            radioButtons("JAK2_homo", "JAK2 Homozygous (>80%)",
                         choices=list("No"=0,"Yes"=1), inline=TRUE),
            br(),
            h4("Risk Category"),
            uiOutput("risk_category_ui"),
            hr(),
            h4("ELN Treatment Recommendation"),
            uiOutput("treatment_rec_ui")
          ),
          box(width=4, title="Patient Summary", status="success", solidHeader=TRUE,
            h4("Current Inputs"),
            tableOutput("patient_summary_table"),
            hr(),
            h4("Annual Thrombosis Risk Estimate"),
            plotOutput("thromb_gauge", height="150px")
          )
        )
      ),

      ## ====== TAB 3: PHARMACOKINETICS =====================================
      tabItem(tabName="tab_pk",
        fluidRow(
          box(width=3, title="PK Parameters", status="primary", solidHeader=TRUE,
            h5("Ruxolitinib"),
            sliderInput("RUX_DOSE", "Dose (mg BID)", 5, 25, 10, step=5),
            sliderInput("CL_RUX",  "Clearance (L/h)", 10, 40, 22, step=1),
            sliderInput("VC_RUX",  "Central Volume (L)", 40, 120, 72, step=5),
            hr(),
            h5("Hydroxyurea"),
            sliderInput("HYU_DOSE", "Daily Dose (mg)", 500, 2000, 500, step=250),
            hr(),
            h5("PEG-IFN-α2a"),
            sliderInput("IFN_DOSE", "Weekly Dose (μg SC)", 45, 180, 45, step=45),
            hr(),
            checkboxGroupInput("PK_DRUGS", "Show Drug PK",
                               choices=c("Ruxolitinib"="rux",
                                         "Hydroxyurea"="hyu",
                                         "PEG-IFN"="ifn"),
                               selected="rux")
          ),
          box(width=9, title="Plasma Concentration-Time Profiles",
              status="info", solidHeader=TRUE,
            plotlyOutput("pk_plot", height="450px"),
            hr(),
            fluidRow(
              valueBoxOutput("rux_cmax_box", width=4),
              valueBoxOutput("rux_tmax_box", width=4),
              valueBoxOutput("rux_auc_box",  width=4)
            )
          )
        ),
        fluidRow(
          box(width=6, title="Ruxolitinib Dose-Concentration Curve",
              plotlyOutput("rux_dose_conc", height="300px")),
          box(width=6, title="JAK2 Inhibition vs Concentration",
              plotlyOutput("jak2_inhib_conc", height="300px"))
        )
      ),

      ## ====== TAB 4: PD & HEMATOLOGY =====================================
      tabItem(tabName="tab_pd",
        fluidRow(
          box(width=3, title="Treatment Selection", status="danger", solidHeader=TRUE,
            h5("Treatment Regimen"),
            checkboxInput("USE_RUX_pd", "Ruxolitinib", FALSE),
            checkboxInput("USE_HYU_pd", "Hydroxyurea", FALSE),
            checkboxInput("USE_IFN_pd", "PEG-IFN-α2a", FALSE),
            checkboxInput("USE_ASP_pd", "Aspirin", TRUE),
            hr(),
            h5("Simulation Duration"),
            sliderInput("SIM_WEEKS", "Simulation (weeks)", 12, 104, 48, step=4),
            hr(),
            h5("Ruxolitinib Dose"),
            sliderInput("RUX_DOSE_pd", "mg BID", 5, 25, 10, step=5),
            h5("HYU Daily Dose"),
            sliderInput("HYU_DOSE_pd", "mg/day", 500, 2000, 500, step=250),
            actionButton("run_sim_pd", "Run Simulation",
                         class="btn-danger btn-block", icon=icon("play"))
          ),
          box(width=9, title="Hematologic Response Over Time",
              status="info", solidHeader=TRUE,
            tabsetPanel(
              tabPanel("Hematocrit",
                plotlyOutput("pd_hct_plot", height="350px")),
              tabPanel("Platelets & WBC",
                plotlyOutput("pd_plt_wbc_plot", height="350px")),
              tabPanel("BFU-E & CFU-E (Erythropoiesis)",
                plotlyOutput("pd_erythroid_plot", height="350px")),
              tabPanel("EPO & pSTAT5",
                plotlyOutput("pd_epo_stat5_plot", height="350px"))
            )
          )
        ),
        fluidRow(
          box(width=6, title="Complete Hematologic Response (CHR)",
              plotlyOutput("chr_timeline", height="280px")),
          box(width=6, title="Symptom Burden (MPN-SAF TSS)",
              plotlyOutput("mpnsaf_plot", height="280px"))
        )
      ),

      ## ====== TAB 5: CLINICAL ENDPOINTS ===================================
      tabItem(tabName="tab_clinical",
        fluidRow(
          box(width=3, title="Endpoint Settings", status="primary", solidHeader=TRUE,
            h5("Risk Model Inputs"),
            sliderInput("AGE_clin", "Patient Age (y)", 30, 85, 65, step=1),
            radioButtons("HIST_clin", "Prior Thrombosis",
                         choices=list("No"=0,"Yes"=1), inline=TRUE),
            hr(),
            sliderInput("SIM_WK_clin", "Follow-up (weeks)", 24, 208, 104, step=8),
            hr(),
            h5("Selected Treatments"),
            checkboxInput("RUX_clin", "Ruxolitinib 10mg BID", FALSE),
            checkboxInput("HYU_clin", "Hydroxyurea 500mg/d", FALSE),
            checkboxInput("IFN_clin", "PEG-IFN-α2a 45μg/wk", FALSE),
            actionButton("run_clin", "Update", class="btn-primary btn-block")
          ),
          box(width=9, title="Clinical Outcomes", status="warning", solidHeader=TRUE,
            tabsetPanel(
              tabPanel("Spleen Volume & SVR35",
                plotlyOutput("spl_svr_plot", height="350px")),
              tabPanel("Thrombosis Risk",
                plotlyOutput("thromb_risk_plot", height="350px")),
              tabPanel("BM Fibrosis",
                plotlyOutput("fibro_plot", height="350px")),
              tabPanel("MF/AML Transformation Risk",
                plotlyOutput("transform_plot", height="350px"))
            )
          )
        ),
        fluidRow(
          box(width=12, title="ELN/IWG-MRT 2013 Response Assessment",
              status="success", solidHeader=TRUE,
            dataTableOutput("response_table"))
        )
      ),

      ## ====== TAB 6: SCENARIO COMPARISON ==================================
      tabItem(tabName="tab_compare",
        fluidRow(
          box(width=3, title="Comparison Settings", status="danger", solidHeader=TRUE,
            h5("Simulation Duration"),
            sliderInput("CMP_WEEKS", "Weeks", 24, 104, 48, step=4),
            hr(),
            h5("Ruxolitinib Dose"),
            sliderInput("RUX_DOSE_cmp", "mg BID", 5, 25, 10, step=5),
            h5("HYU Dose"),
            sliderInput("HYU_DOSE_cmp", "mg/day", 500, 2000, 500, step=250),
            h5("IFN Dose"),
            sliderInput("IFN_DOSE_cmp", "μg/wk", 45, 180, 45, step=45),
            hr(),
            actionButton("run_cmp", "Compare All Scenarios",
                         class="btn-danger btn-block", icon=icon("chart-bar"))
          ),
          box(width=9, title="All Treatment Scenarios — Head-to-Head",
              status="info", solidHeader=TRUE,
            tabsetPanel(
              tabPanel("Hematocrit", plotlyOutput("cmp_hct", height="350px")),
              tabPanel("Platelet Count", plotlyOutput("cmp_plt", height="350px")),
              tabPanel("Spleen Volume", plotlyOutput("cmp_spl", height="350px")),
              tabPanel("JAK2 Allele Burden", plotlyOutput("cmp_allele", height="350px")),
              tabPanel("Symptom Score", plotlyOutput("cmp_mpnsaf", height="350px"))
            )
          )
        ),
        fluidRow(
          box(width=12, title="24-Week & 48-Week Summary Table",
              dataTableOutput("cmp_summary_table"))
        )
      ),

      ## ====== TAB 7: BIOMARKERS ===========================================
      tabItem(tabName="tab_biomarker",
        fluidRow(
          box(width=4, title="Disease Progression Markers", status="purple", solidHeader=TRUE,
            plotlyOutput("allele_plot", height="280px"),
            hr(),
            plotlyOutput("fibro_grade_plot", height="280px")
          ),
          box(width=4, title="Erythropoietic Biomarkers", status="info", solidHeader=TRUE,
            plotlyOutput("epo_plot", height="280px"),
            hr(),
            plotlyOutput("retic_plot", height="280px")
          ),
          box(width=4, title="Pharmacodynamic Biomarkers", status="success", solidHeader=TRUE,
            plotlyOutput("pstat5_plot", height="280px"),
            hr(),
            h4("Biomarker Summary Table"),
            tableOutput("biomarker_table")
          )
        ),
        fluidRow(
          box(width=12, title="Annual Thrombosis Risk Trend",
              plotlyOutput("thromb_annual_plot", height="250px"))
        )
      )

    ) # end tabItems
  ) # end dashboardBody
) # end dashboardPage

## ---------------------------------------------------------------------------
## SERVER
## ---------------------------------------------------------------------------
server <- function(input, output, session) {

  ## ---- Risk Category ----
  output$risk_category_ui <- renderUI({
    age <- input$AGE
    hist <- input$HIST_THROMB == 1
    risk_cat <- if (age > 60 || hist) "HIGH RISK" else "LOW RISK"
    color <- if (risk_cat == "HIGH RISK") "red" else "green"
    tags$div(
      tags$span(risk_cat, style=sprintf("color:%s; font-size:20px; font-weight:bold;", color)),
      br(),
      if (risk_cat == "HIGH RISK")
        tags$p("Criteria: Age >60 or prior thrombosis → Cytoreductive therapy indicated")
      else
        tags$p("Criteria: Age ≤60 and no thrombosis → Phlebotomy + Aspirin first-line")
    )
  })

  output$treatment_rec_ui <- renderUI({
    age <- input$AGE
    hist <- input$HIST_THROMB == 1
    if (age > 60 || hist) {
      tags$ul(
        tags$li("First-line: Hydroxyurea 500-1000 mg/d"),
        tags$li("Alternative: Ruxolitinib (HYU-resistant/intolerant)"),
        tags$li("Consider: PEG-IFN-α2a (younger patients, childbearing)"),
        tags$li("All: Low-dose Aspirin + Phlebotomy to Hct <45%")
      )
    } else {
      tags$ul(
        tags$li("First-line: Phlebotomy to Hct <45%"),
        tags$li("All: Low-dose Aspirin 81mg/d"),
        tags$li("Cytoreductive if poor phlebotomy tolerance or symptomatic")
      )
    }
  })

  ## ---- Patient Summary Table ----
  output$patient_summary_table <- renderTable({
    data.frame(
      Parameter = c("Age", "Sex", "JAK2 Allele Burden",
                    "Baseline Hct", "Baseline Platelets", "Spleen"),
      Value = c(paste(input$AGE, "y"),
                ifelse(input$SEX == 1, "Male", "Female"),
                paste0(input$JAK2_AB0, "%"),
                paste0(input$Hct_base, "%"),
                paste0(input$PLT_base, " ×10⁹/L"),
                paste0(input$SPL_base, " mL"))
    )
  })

  ## ---- Thrombosis gauge ----
  output$thromb_gauge <- renderPlot({
    age <- input$AGE
    hist <- as.numeric(input$HIST_THROMB)
    Hct <- input$Hct_base
    risk <- (Hct > 45 ? 0.06 : 0.02) +
            (age > 60 ? (age-60)*0.002 : 0) +
            hist * 0.15
    risk <- min(risk, 0.30)
    df <- data.frame(x=1, y=risk, label=sprintf("%.1f%%", risk*100))
    ggplot(df, aes(x=x, y=y, fill=y)) +
      geom_bar(stat="identity", width=0.5) +
      geom_text(aes(label=label), vjust=-0.5, size=5, fontface="bold") +
      scale_y_continuous(limits=c(0,0.3), labels=percent_format()) +
      scale_fill_gradient(low="green", high="red", limits=c(0,0.3)) +
      labs(title="Annual Thrombosis Risk", x="", y="") +
      theme_minimal(base_size=12) + theme(legend.position="none",
                                           axis.text.x=element_blank())
  })

  ## ---- Overview path plot ----
  output$overview_path_plot <- renderPlot({
    df <- data.frame(
      Step = c("JAK2 V617F\nMutation", "Constitutive\nJAK-STAT5", "BFU-E/CFU-E\nHyperproliferation",
               "RBC Mass\nExpansion", "Elevated\nHematocrit", "Thrombosis\nRisk"),
      x = 1:6,
      y = c(3,4,4,3.5,3,2)
    )
    ggplot(df, aes(x=x, y=y, label=Step)) +
      geom_path(arrow=arrow(type="closed", length=unit(0.25,"cm")), color="#E74C3C", linewidth=1.2) +
      geom_point(size=5, color="#E74C3C") +
      geom_label(fontface="bold", size=2.8, fill="#FCE4EC") +
      theme_void(base_size=10) +
      labs(title="PV Pathophysiology Cascade")
  })

  ## ---- PK Simulation ----
  pk_sim_rux <- reactive({
    req("rux" %in% input$PK_DRUGS)
    ev_r <- ev(cmt="DEPOT_RUX", amt=input$RUX_DOSE, ii=12, addl=27, time=0)
    mod %>%
      param(USE_RUX=1, USE_HYU=0, USE_IFN=0,
            CL_RUX=input$CL_RUX, VC_RUX=input$VC_RUX) %>%
      mrgsim(ev=ev_r, end=168, delta=0.25) %>%
      as.data.frame() %>%
      mutate(time_h=time*24)
  })

  output$pk_plot <- renderPlotly({
    p <- ggplot()
    if ("rux" %in% input$PK_DRUGS) {
      df <- pk_sim_rux()
      p <- p + geom_line(data=df, aes(x=time_h, y=CP_RUX_out, color="Ruxolitinib (ng/mL)"), linewidth=1)
    }
    p <- p +
      scale_color_manual(values=c("Ruxolitinib (ng/mL)"="#E74C3C",
                                   "Hydroxyurea (μM)"="#3498DB",
                                   "PEG-IFN (pg/mL)"="#9B59B6")) +
      labs(title="Plasma Concentration-Time", x="Time (hours)", y="Concentration",
           color="Drug") +
      theme_bw(base_size=12)
    ggplotly(p, tooltip=c("x","y","color"))
  })

  output$rux_cmax_box <- renderValueBox({
    df <- pk_sim_rux()
    valueBox(round(max(df$CP_RUX_out[df$time_h <= 24]),1),
             "Cmax (ng/mL)", icon=icon("arrow-up"), color="red")
  })
  output$rux_tmax_box <- renderValueBox({
    df <- pk_sim_rux()
    tmax <- df$time_h[which.max(df$CP_RUX_out[df$time_h <= 24])]
    valueBox(round(tmax,1), "Tmax (h)", icon=icon("clock"), color="orange")
  })
  output$rux_auc_box <- renderValueBox({
    df <- pk_sim_rux()
    d24 <- df[df$time_h <= 24,]
    auc <- round(sum(diff(d24$time_h) * (head(d24$CP_RUX_out,-1)+tail(d24$CP_RUX_out,-1))/2), 0)
    valueBox(auc, "AUC₀₋₂₄ (ng·h/mL)", icon=icon("chart-area"), color="blue")
  })

  output$rux_dose_conc <- renderPlotly({
    doses <- c(5,10,15,20,25)
    auc_vals <- sapply(doses, function(d) {
      ev_d <- ev(cmt="DEPOT_RUX", amt=d, ii=12, addl=1, time=0)
      s <- mod %>% param(USE_RUX=1,USE_HYU=0,USE_IFN=0) %>%
        mrgsim(ev=ev_d,end=24,delta=0.5) %>% as.data.frame()
      sum(diff(s$time * 24) * (head(s$CP_RUX_out,-1)+tail(s$CP_RUX_out,-1))/2)
    })
    df <- data.frame(Dose=doses, AUC=auc_vals)
    p <- ggplot(df, aes(x=Dose,y=AUC)) +
      geom_point(size=4,color="#E74C3C") + geom_line(color="#E74C3C",linewidth=1) +
      labs(title="Ruxolitinib Dose vs AUC₀₋₂₄",x="Dose (mg BID)",y="AUC (ng·h/mL)") +
      theme_bw(base_size=11)
    ggplotly(p)
  })

  output$jak2_inhib_conc <- renderPlotly({
    cp <- seq(0,50,0.5)
    inh <- 1.0 * cp^1 / (0.86^1 + cp^1)
    df <- data.frame(Cp=cp, Inhibition=inh*100)
    p <- ggplot(df, aes(x=Cp,y=Inhibition)) +
      geom_line(color="#2ECC71",linewidth=1.5) +
      geom_vline(xintercept=0.86, linetype="dashed", color="red") +
      annotate("text",x=2,y=30,label="IC50 = 0.86 ng/mL",color="red",size=3) +
      labs(title="JAK2 Inhibition vs Ruxolitinib Cp",
           x="Concentration (ng/mL)",y="% Inhibition") +
      theme_bw(base_size=11)
    ggplotly(p)
  })

  ## ---- PD Simulation ----
  pd_data <- eventReactive(input$run_sim_pd, {
    params <- list(
      USE_RUX = as.integer(input$USE_RUX_pd),
      USE_HYU = as.integer(input$USE_HYU_pd),
      USE_IFN = as.integer(input$USE_IFN_pd),
      USE_ASP = as.integer(input$USE_ASP_pd)
    )
    evs <- NULL
    if (input$USE_RUX_pd) {
      evs <- ev(cmt="DEPOT_RUX", amt=input$RUX_DOSE_pd, ii=12,
                addl=input$SIM_WEEKS*14-1, time=0)
    }
    if (input$USE_IFN_pd) {
      e2 <- ev(cmt="SC_IFN", amt=45, ii=168, addl=input$SIM_WEEKS-1, time=0)
      evs <- if(is.null(evs)) e2 else ev_seq(evs, e2)
    }
    sim <- mod %>% param(params) %>%
      {if(!is.null(evs)) mrgsim(., ev=evs, end=input$SIM_WEEKS*7, delta=1)
       else mrgsim(., end=input$SIM_WEEKS*7, delta=1)} %>%
      as.data.frame() %>% mutate(time_wk=time/7)
    sim
  }, ignoreNULL=FALSE)

  output$pd_hct_plot <- renderPlotly({
    df <- pd_data()
    p <- ggplot(df, aes(x=time_wk, y=Hct)) +
      geom_line(color="#E74C3C", linewidth=1.2) +
      geom_hline(yintercept=45, linetype="dashed", color="black") +
      annotate("text", x=max(df$time_wk)*0.8, y=46.5,
               label="Target Hct <45%", size=3) +
      labs(title="Hematocrit Over Time", x="Time (weeks)", y="Hct (%)") +
      theme_bw(base_size=12)
    ggplotly(p, tooltip=c("x","y"))
  })

  output$pd_plt_wbc_plot <- renderPlotly({
    df <- pd_data()
    df_long <- df %>%
      select(time_wk, PLT_out, WBC_out) %>%
      tidyr::pivot_longer(cols=c(PLT_out, WBC_out),
                          names_to="Variable", values_to="Value") %>%
      mutate(Variable = ifelse(Variable=="PLT_out","Platelets (×10⁹/L)","WBC (×10⁹/L)"))
    p <- ggplot(df_long, aes(x=time_wk, y=Value, color=Variable)) +
      geom_line(linewidth=1.2) +
      facet_wrap(~Variable, scales="free_y") +
      scale_color_manual(values=c("Platelets (×10⁹/L)"="#9B59B6","WBC (×10⁹/L)"="#3498DB")) +
      labs(title="Platelet & WBC Counts", x="Time (weeks)", y="Count") +
      theme_bw(base_size=11) + theme(legend.position="none")
    ggplotly(p)
  })

  output$pd_erythroid_plot <- renderPlotly({
    df <- pd_data()
    df_long <- df %>%
      select(time_wk, BFUE, CFUE) %>%
      tidyr::pivot_longer(cols=c(BFUE, CFUE)) %>%
      mutate(name = ifelse(name=="BFUE","BFU-E Progenitors","CFU-E Progenitors"))
    p <- ggplot(df_long, aes(x=time_wk, y=value, color=name)) +
      geom_line(linewidth=1.2) +
      scale_color_manual(values=c("BFU-E Progenitors"="#E67E22","CFU-E Progenitors"="#F39C12")) +
      labs(title="Erythroid Progenitor Pools", x="Time (weeks)", y="Progenitor Count (a.u.)") +
      theme_bw(base_size=12)
    ggplotly(p)
  })

  output$pd_epo_stat5_plot <- renderPlotly({
    df <- pd_data()
    p1 <- ggplot(df, aes(x=time_wk, y=EPO)) +
      geom_line(color="#8E44AD", linewidth=1.2) +
      geom_hline(yintercept=10, linetype="dashed") +
      labs(title="Serum EPO (mU/mL)", x="Weeks", y="EPO") +
      theme_bw(base_size=11)
    p2 <- ggplot(df, aes(x=time_wk, y=pSTAT5*100)) +
      geom_line(color="#1ABC9C", linewidth=1.2) +
      labs(title="pSTAT5 Inhibition (%)", x="Weeks", y="% Inhibition") +
      theme_bw(base_size=11)
    subplot(ggplotly(p1), ggplotly(p2), nrows=1, shareX=TRUE)
  })

  output$chr_timeline <- renderPlotly({
    df <- pd_data()
    p <- ggplot(df, aes(x=time_wk, y=CHR)) +
      geom_area(fill="#2ECC71", alpha=0.4) +
      geom_line(color="#27AE60", linewidth=1.2) +
      scale_y_continuous(limits=c(0,1.1), labels=c("Not CHR","CHR")) +
      labs(title="Complete Hematologic Response", x="Time (weeks)", y="CHR Status") +
      theme_bw(base_size=12)
    ggplotly(p)
  })

  output$mpnsaf_plot <- renderPlotly({
    df <- pd_data()
    p <- ggplot(df, aes(x=time_wk, y=MPN_SAF)) +
      geom_line(color="#E74C3C", linewidth=1.2) +
      geom_hline(yintercept=20, linetype="dashed", color="orange") +
      annotate("text",x=max(df$time_wk)*0.8,y=22,
               label="Response threshold (TSS<20)",size=3,color="orange") +
      labs(title="MPN-SAF Total Symptom Score", x="Time (weeks)", y="TSS Score") +
      theme_bw(base_size=12)
    ggplotly(p)
  })

  ## ---- Clinical Endpoints Simulation ----
  clin_data <- eventReactive(input$run_clin, {
    params <- list(
      USE_RUX = as.integer(input$RUX_clin),
      USE_HYU = as.integer(input$HYU_clin),
      USE_IFN = as.integer(input$IFN_clin),
      USE_ASP = 1
    )
    evs <- NULL
    if (input$RUX_clin) {
      evs <- ev(cmt="DEPOT_RUX", amt=10, ii=12, addl=input$SIM_WK_clin*14-1, time=0)
    }
    if (input$IFN_clin) {
      e2 <- ev(cmt="SC_IFN", amt=45, ii=168, addl=input$SIM_WK_clin-1, time=0)
      evs <- if(is.null(evs)) e2 else ev_seq(evs, e2)
    }
    mod %>% param(params) %>%
      {if(!is.null(evs)) mrgsim(., ev=evs, end=input$SIM_WK_clin*7, delta=1)
       else mrgsim(., end=input$SIM_WK_clin*7, delta=1)} %>%
      as.data.frame() %>% mutate(time_wk=time/7)
  }, ignoreNULL=FALSE)

  output$spl_svr_plot <- renderPlotly({
    df <- clin_data()
    p1 <- ggplot(df, aes(x=time_wk, y=SPL_vol)) +
      geom_line(color="#FF9800", linewidth=1.2) +
      geom_hline(yintercept=450, linetype="dashed") +
      labs(title="Spleen Volume (mL)",x="Weeks",y="mL") + theme_bw(base_size=11)
    p2 <- ggplot(df, aes(x=time_wk, y=SVR)) +
      geom_line(color="#4CAF50", linewidth=1.2) +
      geom_hline(yintercept=35, linetype="dashed", color="red") +
      annotate("text",x=max(df$time_wk)*0.7,y=37,label="SVR35 threshold",size=3,color="red") +
      labs(title="Spleen Volume Reduction (%)",x="Weeks",y="SVR (%)") + theme_bw(base_size=11)
    subplot(ggplotly(p1), ggplotly(p2), nrows=1)
  })

  output$thromb_risk_plot <- renderPlotly({
    df <- clin_data()
    p <- ggplot(df, aes(x=time_wk, y=THROMB_RISK*100)) +
      geom_line(color="#E74C3C", linewidth=1.2) +
      geom_ribbon(aes(ymin=0, ymax=THROMB_RISK*100), fill="#E74C3C", alpha=0.2) +
      labs(title="Annual Thrombosis Risk (%)",x="Time (weeks)",y="Annual Risk (%)") +
      theme_bw(base_size=12)
    ggplotly(p)
  })

  output$fibro_plot <- renderPlotly({
    df <- clin_data()
    cols_grade <- c("Grade 0\n(0-0.5)"="#4CAF50","Grade 1\n(0.5-1.5)"="#FF9800",
                    "Grade 2\n(1.5-2.5)"="#FF5722","Grade 3\n(2.5-3)"="#B71C1C")
    p <- ggplot(df, aes(x=time_wk, y=FIBRO_out)) +
      geom_line(color="#795548", linewidth=1.5) +
      geom_hline(yintercept=c(0.5,1.5,2.5), linetype="dashed",
                 color=c("#4CAF50","#FF9800","#FF5722")) +
      annotate("text",x=max(df$time_wk)*0.9,y=c(0.25,1.0,2.0,2.75),
               label=c("MF-0","MF-1","MF-2","MF-3"),size=3,color="#795548") +
      scale_y_continuous(limits=c(0,3.2)) +
      labs(title="BM Reticulin Fibrosis Grade",x="Time (weeks)",y="Fibrosis Score") +
      theme_bw(base_size=12)
    ggplotly(p)
  })

  output$transform_plot <- renderPlotly({
    df <- clin_data()
    df$cum_MF_risk <- cumsum(df$MF_risk / 365) * 100
    p <- ggplot(df, aes(x=time_wk, y=cum_MF_risk)) +
      geom_line(color="#9C27B0", linewidth=1.5) +
      labs(title="Cumulative MF/AML Transformation Risk (%)",
           x="Time (weeks)",y="Cumulative Risk (%)") +
      theme_bw(base_size=12)
    ggplotly(p)
  })

  output$response_table <- renderDataTable({
    df <- clin_data()
    wk24 <- df[which.min(abs(df$time_wk - 24)),]
    wk48 <- df[which.min(abs(df$time_wk - 48)),]
    data.frame(
      Timepoint = c("Week 24","Week 48"),
      Hct_pct = round(c(wk24$Hct, wk48$Hct), 1),
      PLT_109L = round(c(wk24$PLT_out, wk48$PLT_out)),
      WBC_109L = round(c(wk24$WBC_out, wk48$WBC_out), 1),
      SPL_mL = round(c(wk24$SPL_vol, wk48$SPL_vol)),
      SVR_pct = round(c(wk24$SVR, wk48$SVR), 1),
      SVR35 = c(ifelse(wk24$SVR >= 35,"Yes","No"), ifelse(wk48$SVR >= 35,"Yes","No")),
      CHR = c(ifelse(wk24$CHR==1,"Complete","Not CHR"), ifelse(wk48$CHR==1,"Complete","Not CHR")),
      Allele_pct = round(c(wk24$ALLELE_out, wk48$ALLELE_out), 1),
      Fibrosis = round(c(wk24$FIBRO_out, wk48$FIBRO_out), 2)
    )
  }, options=list(pageLength=5, dom="t"))

  ## ---- Scenario Comparison ----
  cmp_data <- eventReactive(input$run_cmp, {
    scenarios <- list(
      "Untreated" = list(USE_RUX=0, USE_HYU=0, USE_IFN=0, USE_ASP=0),
      "Aspirin Only" = list(USE_RUX=0, USE_HYU=0, USE_IFN=0, USE_ASP=1),
      "Hydroxyurea" = list(USE_RUX=0, USE_HYU=1, USE_IFN=0, USE_ASP=1),
      "Ruxolitinib" = list(USE_RUX=1, USE_HYU=0, USE_IFN=0, USE_ASP=1),
      "PEG-IFN-α2a" = list(USE_RUX=0, USE_HYU=0, USE_IFN=1, USE_ASP=1),
      "RUX + ASP" = list(USE_RUX=1, USE_HYU=0, USE_IFN=0, USE_ASP=1)
    )
    bind_rows(lapply(names(scenarios), function(nm) {
      p <- scenarios[[nm]]
      evs <- NULL
      if (p$USE_RUX == 1)
        evs <- ev(cmt="DEPOT_RUX", amt=input$RUX_DOSE_cmp, ii=12,
                  addl=input$CMP_WEEKS*14-1, time=0)
      if (p$USE_IFN == 1) {
        e2 <- ev(cmt="SC_IFN", amt=input$IFN_DOSE_cmp, ii=168,
                 addl=input$CMP_WEEKS-1, time=0)
        evs <- if(is.null(evs)) e2 else ev_seq(evs, e2)
      }
      mod %>% param(p) %>%
        {if(!is.null(evs)) mrgsim(., ev=evs, end=input$CMP_WEEKS*7, delta=2)
         else mrgsim(., end=input$CMP_WEEKS*7, delta=2)} %>%
        as.data.frame() %>%
        mutate(Scenario=nm, time_wk=time/7)
    }))
  }, ignoreNULL=FALSE)

  make_cmp_plot <- function(var, title, ylab, hline=NULL, hline_label=NULL) {
    renderPlotly({
      df <- cmp_data()
      p <- ggplot(df, aes_string(x="time_wk", y=var, color="Scenario")) +
        geom_line(linewidth=1.2)
      if (!is.null(hline))
        p <- p + geom_hline(yintercept=hline, linetype="dashed") +
                 annotate("text", x=max(df$time_wk)*0.8, y=hline*1.04,
                          label=hline_label, size=3)
      p <- p + scale_color_brewer(palette="Set1") +
        labs(title=title, x="Time (weeks)", y=ylab) + theme_bw(base_size=11)
      ggplotly(p)
    })
  }

  output$cmp_hct    <- make_cmp_plot("Hct","Hematocrit","%",45,"Target 45%")
  output$cmp_plt    <- make_cmp_plot("PLT_out","Platelet Count","×10⁹/L",400,"Normal 400")
  output$cmp_spl    <- make_cmp_plot("SPL_vol","Spleen Volume","mL",450,"Baseline 450")
  output$cmp_allele <- make_cmp_plot("ALLELE_out","JAK2 V617F Allele Burden","%")
  output$cmp_mpnsaf <- make_cmp_plot("MPN_SAF","MPN-SAF TSS","Score",20,"Response")

  output$cmp_summary_table <- renderDataTable({
    df <- cmp_data()
    wk24 <- df %>% group_by(Scenario) %>%
      filter(time_wk == time_wk[which.min(abs(time_wk - 24))]) %>%
      slice(1) %>%
      select(Scenario, Hct, PLT_out, WBC_out, SPL_vol, ALLELE_out, CHR, SVR)
    wk48 <- df %>% group_by(Scenario) %>%
      filter(time_wk == time_wk[which.min(abs(time_wk - 48))]) %>%
      slice(1) %>%
      select(Scenario, Hct, PLT_out, WBC_out, ALLELE_out, CHR)
    bind_cols(
      wk24 %>% rename_with(~paste0(.,"_wk24"), -Scenario),
      wk48 %>% ungroup() %>% select(-Scenario) %>% rename_with(~paste0(.,"_wk48"))
    ) %>%
      mutate(across(where(is.numeric), ~round(.,1)))
  }, options=list(scrollX=TRUE, pageLength=6))

  ## ---- Biomarker Tab ----
  bio_data <- reactive({
    # Run untreated + ruxolitinib for biomarker comparison
    scen <- list(
      "Untreated" = list(USE_RUX=0, USE_HYU=0, USE_IFN=0, USE_ASP=0),
      "Ruxolitinib 10mg BID" = list(USE_RUX=1, USE_HYU=0, USE_IFN=0, USE_ASP=1),
      "PEG-IFN-α2a" = list(USE_RUX=0, USE_HYU=0, USE_IFN=1, USE_ASP=1)
    )
    bind_rows(lapply(names(scen), function(nm) {
      p <- scen[[nm]]
      evs <- NULL
      if (p$USE_RUX == 1)
        evs <- ev(cmt="DEPOT_RUX", amt=10, ii=12, addl=72*14-1, time=0)
      if (p$USE_IFN == 1)
        evs <- ev(cmt="SC_IFN", amt=45, ii=168, addl=71, time=0)
      mod %>% param(p) %>%
        {if(!is.null(evs)) mrgsim(., ev=evs, end=72*7, delta=2)
         else mrgsim(., end=72*7, delta=2)} %>%
        as.data.frame() %>% mutate(Scenario=nm, time_wk=time/7)
    }))
  })

  output$allele_plot <- renderPlotly({
    df <- bio_data()
    p <- ggplot(df, aes(x=time_wk, y=ALLELE_out, color=Scenario)) +
      geom_line(linewidth=1.2) +
      scale_color_manual(values=c("Untreated"="#E74C3C",
                                   "Ruxolitinib 10mg BID"="#2ECC71",
                                   "PEG-IFN-α2a"="#9B59B6")) +
      labs(title="JAK2 V617F Allele Burden (%)",x="Weeks",y="Allele Burden (%)") +
      theme_bw(base_size=11)
    ggplotly(p)
  })

  output$fibro_grade_plot <- renderPlotly({
    df <- bio_data()
    p <- ggplot(df, aes(x=time_wk, y=FIBRO_out, color=Scenario)) +
      geom_line(linewidth=1.2) +
      scale_y_continuous(limits=c(0,3)) +
      scale_color_manual(values=c("Untreated"="#E74C3C",
                                   "Ruxolitinib 10mg BID"="#2ECC71",
                                   "PEG-IFN-α2a"="#9B59B6")) +
      labs(title="BM Fibrosis Score",x="Weeks",y="Fibrosis (0-3)") +
      theme_bw(base_size=11)
    ggplotly(p)
  })

  output$epo_plot <- renderPlotly({
    df <- bio_data()
    p <- ggplot(df, aes(x=time_wk, y=EPO, color=Scenario)) +
      geom_line(linewidth=1.2) +
      scale_color_manual(values=c("Untreated"="#E74C3C",
                                   "Ruxolitinib 10mg BID"="#2ECC71",
                                   "PEG-IFN-α2a"="#9B59B6")) +
      labs(title="Serum EPO (mU/mL)",x="Weeks",y="EPO") +
      theme_bw(base_size=11)
    ggplotly(p)
  })

  output$retic_plot <- renderPlotly({
    df <- bio_data()
    p <- ggplot(df, aes(x=time_wk, y=RETIC_C, color=Scenario)) +
      geom_line(linewidth=1.2) +
      scale_color_manual(values=c("Untreated"="#E74C3C",
                                   "Ruxolitinib 10mg BID"="#2ECC71",
                                   "PEG-IFN-α2a"="#9B59B6")) +
      labs(title="Circulating Reticulocytes",x="Weeks",y="Count (a.u.)") +
      theme_bw(base_size=11)
    ggplotly(p)
  })

  output$pstat5_plot <- renderPlotly({
    df <- bio_data()
    p <- ggplot(df, aes(x=time_wk, y=pSTAT5*100, color=Scenario)) +
      geom_line(linewidth=1.2) +
      scale_color_manual(values=c("Untreated"="#E74C3C",
                                   "Ruxolitinib 10mg BID"="#2ECC71",
                                   "PEG-IFN-α2a"="#9B59B6")) +
      labs(title="pSTAT5 Inhibition (%)",x="Weeks",y="% Inhibition") +
      theme_bw(base_size=11)
    ggplotly(p)
  })

  output$biomarker_table <- renderTable({
    df <- bio_data()
    wk48 <- df %>% filter(time_wk >= 47 & time_wk <= 49) %>%
      group_by(Scenario) %>% slice(1) %>%
      select(Scenario, ALLELE_out, FIBRO_out, EPO, pSTAT5)
    wk48 %>%
      mutate(ALLELE_out=round(ALLELE_out,1),
             FIBRO_out=round(FIBRO_out,2),
             EPO=round(EPO,1),
             pSTAT5=round(pSTAT5*100,1)) %>%
      rename(`JAK2 Allele (%)`=ALLELE_out,
             `BM Fibrosis`=FIBRO_out,
             `EPO (mU/mL)`=EPO,
             `pSTAT5 inhib (%)`=pSTAT5)
  })

  output$thromb_annual_plot <- renderPlotly({
    df <- bio_data()
    p <- ggplot(df, aes(x=time_wk, y=THROMB_RISK*100, color=Scenario)) +
      geom_line(linewidth=1.2) +
      scale_color_manual(values=c("Untreated"="#E74C3C",
                                   "Ruxolitinib 10mg BID"="#2ECC71",
                                   "PEG-IFN-α2a"="#9B59B6")) +
      labs(title="Annual Thrombosis Risk (%)",x="Time (weeks)",y="Annual Risk (%)") +
      theme_bw(base_size=11)
    ggplotly(p)
  })
}

## ---------------------------------------------------------------------------
## Run App
## ---------------------------------------------------------------------------
shinyApp(ui = ui, server = server)
