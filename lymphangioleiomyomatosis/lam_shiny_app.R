## ============================================================
## Lymphangioleiomyomatosis (LAM) — Interactive Shiny QSP App
## ============================================================
## Tabs:
##   1. Patient Profile & Disease Severity
##   2. Drug PK (Concentration-Time)
##   3. mTOR Pathway Activity (PD Biomarkers)
##   4. Clinical Endpoints (FEV1, DLCO, Cyst Volume, 6MWT)
##   5. Treatment Scenario Comparison
##   6. Biomarker Dashboard (VEGF-D, MMP, S6K1, AML)
## ============================================================

library(shiny)
library(shinydashboard)
library(mrgsolve)
library(dplyr)
library(ggplot2)
library(plotly)
library(DT)
library(scales)

## ====================================================================
## mrgsolve MODEL (inline)
## ====================================================================
lam_code <- '
$PARAM
ka_siro=0.5, F_siro=0.15, CL_siro=10.0, V1_siro=62.5, Q_siro=8.0, V2_siro=1400
ka_ever=1.8, F_ever=0.30, CL_ever=20.0, V1_ever=90.0, Q_ever=15.0, V2_ever=500
Imax_siro=0.90, IC50_siro=2.5, hill_siro=1.5
Imax_ever=0.88, IC50_ever=3.0, hill_ever=1.4
mTORC1_LAM=4.0, kmTOR_act=0.5, kmTOR_deact=0.3
kS6K_on=0.8, kS6K_off=0.2
k4EBP_on=1.0, k4EBP_off=0.3
kLAM_prolif=0.008, kLAM_death=0.006, E2_LAMstim=1.4
kVEGFD_prod=5.0, kVEGFD_clear=0.015
kMMP_prod=0.1, kMMP_clear=0.2
E2_basal=1.0, kE2_turnover=0.05
Rheb_TSC_loss=2.5, kRheb_on=0.1, kRheb_off=0.05
kCyst_form=0.0005
kFEV1_decline=0.0137, FEV1_min=20
kDLCO_decline=0.011, DLCO_min=30
kAML_grow=0.0003, kAML_shrink=0.002

$CMT
SIRO_GUT SIRO_C SIRO_P EVER_GUT EVER_C EVER_P
RHEB_GTP MTORC1 S6K1_P EBPP1
LAM_CELLS VEGFD MMP_ACT ESTROGEN
CYST_VOL FEV1_PCT DLCO_PCT AML_VOL

$INIT
SIRO_GUT=0, SIRO_C=0, SIRO_P=0
EVER_GUT=0, EVER_C=0, EVER_P=0
RHEB_GTP=2.5, MTORC1=3.5, S6K1_P=3.0, EBPP1=3.0
LAM_CELLS=1.0, VEGFD=1500, MMP_ACT=2.0, ESTROGEN=1.0
CYST_VOL=25.0, FEV1_PCT=72.0, DLCO_PCT=60.0, AML_VOL=120

$ODE
double Cp_siro = SIRO_C / V1_siro * 1000;
double Cp_ever = EVER_C / V1_ever * 1000;
dxdt_SIRO_GUT = -ka_siro * SIRO_GUT;
dxdt_SIRO_C = ka_siro*F_siro*SIRO_GUT - (CL_siro/V1_siro)*SIRO_C - (Q_siro/V1_siro)*SIRO_C + (Q_siro/V2_siro)*SIRO_P;
dxdt_SIRO_P = (Q_siro/V1_siro)*SIRO_C - (Q_siro/V2_siro)*SIRO_P;
dxdt_EVER_GUT = -ka_ever * EVER_GUT;
dxdt_EVER_C = ka_ever*F_ever*EVER_GUT - (CL_ever/V1_ever)*EVER_C - (Q_ever/V1_ever)*EVER_C + (Q_ever/V2_ever)*EVER_P;
dxdt_EVER_P = (Q_ever/V1_ever)*EVER_C - (Q_ever/V2_ever)*EVER_P;
double I_siro = Imax_siro*pow(Cp_siro,hill_siro)/(pow(IC50_siro,hill_siro)+pow(Cp_siro,hill_siro));
double I_ever = Imax_ever*pow(Cp_ever,hill_ever)/(pow(IC50_ever,hill_ever)+pow(Cp_ever,hill_ever));
double I_drug = 1-(1-I_siro)*(1-I_ever);
double Rheb_target = Rheb_TSC_loss*(1-0.1*I_drug);
dxdt_RHEB_GTP = kRheb_on*(Rheb_target-RHEB_GTP) - kRheb_off*RHEB_GTP;
double mTOR_target = mTORC1_LAM*(RHEB_GTP/Rheb_TSC_loss)*(1-I_drug);
dxdt_MTORC1 = kmTOR_act*(mTOR_target-MTORC1) - kmTOR_deact*MTORC1*0.1;
dxdt_S6K1_P = kS6K_on*(3.0*(MTORC1/mTORC1_LAM)-S6K1_P) - kS6K_off*S6K1_P*0.1;
dxdt_EBPP1 = k4EBP_on*(3.0*(MTORC1/mTORC1_LAM)-EBPP1) - k4EBP_off*EBPP1*0.1;
dxdt_ESTROGEN = kE2_turnover*(E2_basal - ESTROGEN);
double LAM_prolif_rate = kLAM_prolif*(MTORC1/mTORC1_LAM)*ESTROGEN*E2_LAMstim;
dxdt_LAM_CELLS = LAM_prolif_rate*LAM_CELLS - kLAM_death*LAM_CELLS;
dxdt_VEGFD = kVEGFD_prod*LAM_CELLS - kVEGFD_clear*VEGFD;
dxdt_MMP_ACT = kMMP_prod*LAM_CELLS - kMMP_clear*MMP_ACT;
dxdt_CYST_VOL = kCyst_form*MMP_ACT*LAM_CELLS;
double FEV1_decline = kFEV1_decline*(CYST_VOL/25.0);
double FEV1_recov = 0.002*I_drug*(FEV1_PCT<75?1:0);
dxdt_FEV1_PCT = -FEV1_decline + FEV1_recov;
if(FEV1_PCT<FEV1_min) dxdt_FEV1_PCT=0;
double DLCO_decline = kDLCO_decline*(CYST_VOL/25.0);
dxdt_DLCO_PCT = -DLCO_decline + 0.001*I_drug*(DLCO_PCT<65?1:0);
if(DLCO_PCT<DLCO_min) dxdt_DLCO_PCT=0;
dxdt_AML_VOL = kAML_grow*LAM_CELLS - kAML_shrink*I_drug*AML_VOL;
if(AML_VOL<0) dxdt_AML_VOL=0;

$TABLE
capture Cp_siro_ngml = SIRO_C/V1_siro*1000;
capture Cp_ever_ngml = EVER_C/V1_ever*1000;
capture mTORC1_inhib_pct = (1-(1-(Imax_siro*pow(SIRO_C/V1_siro*1000,hill_siro)/(pow(IC50_siro,hill_siro)+pow(SIRO_C/V1_siro*1000,hill_siro))))*(1-(Imax_ever*pow(EVER_C/V1_ever*1000,hill_ever)/(pow(IC50_ever,hill_ever)+pow(EVER_C/V1_ever*1000,hill_ever)))))*100;
capture S6K1_phos = S6K1_P;
capture EBP1_phos = EBPP1;
capture VEGFD_pgmL = VEGFD;
capture MMP_normalized = MMP_ACT;
capture FEV1_pct = FEV1_PCT;
capture DLCO_pct = DLCO_PCT;
capture CystVol = CYST_VOL;
capture AML_mL = AML_VOL;
capture LAM_burden = LAM_CELLS;
'

lam_mod <- mcode("LAM_QSP_shiny", lam_code, quiet=TRUE)

## ====================================================================
## SIMULATION HELPER
## ====================================================================
run_sim <- function(mod, drug="Sirolimus", dose=2, duration_wk=104,
                    E2_level="Premenopausal",
                    incl_gnrh=FALSE) {
  # E2 baseline
  E2_val <- switch(E2_level,
                   "Premenopausal"   = 1.0,
                   "Postmenopausal"  = 0.2,
                   "GnRH Suppressed" = 0.05,
                   1.0)
  mod2 <- mod %>% param(E2_basal = E2_val)

  end_h <- duration_wk * 168
  evs <- NULL
  if(drug == "Sirolimus") {
    evs <- ev(amt=dose, cmt=1, ii=24, addl=duration_wk*7-1, time=0)
  } else if(drug == "Everolimus") {
    evs <- ev(amt=dose, cmt=4, ii=24, addl=duration_wk*7-1, time=0)
  } else if(drug == "Sirolimus + Everolimus") {
    evs <- ev_seq(ev(amt=dose, cmt=1, ii=24, addl=duration_wk*7-1),
                  ev(amt=10, cmt=4, ii=24, addl=duration_wk*7-1))
  }
  # else: no treatment

  if(!is.null(evs)) {
    out <- mod2 %>% ev(evs) %>%
      mrgsim(end=end_h, delta=168, add=seq(0,336,by=4)) %>%
      as_tibble()
  } else {
    out <- mod2 %>%
      mrgsim(end=end_h, delta=168) %>%
      as_tibble()
  }
  out %>% mutate(time_weeks=time/168, drug_label=drug)
}

## ====================================================================
## UI
## ====================================================================
ui <- dashboardPage(
  skin="blue",
  dashboardHeader(title="LAM QSP Model", titleWidth=280),

  dashboardSidebar(
    width=280,
    sidebarMenu(
      menuItem("Patient Profile",    tabName="tab_profile",  icon=icon("user-circle")),
      menuItem("Drug PK",            tabName="tab_pk",       icon=icon("flask")),
      menuItem("mTOR Pathway (PD)",  tabName="tab_pd",       icon=icon("dna")),
      menuItem("Clinical Endpoints", tabName="tab_endpoints",icon=icon("heartbeat")),
      menuItem("Scenario Comparison",tabName="tab_compare",  icon=icon("chart-line")),
      menuItem("Biomarker Dashboard",tabName="tab_biomarker",icon=icon("microscope"))
    ),
    hr(),
    h5("Treatment Settings", style="padding-left:15px; color:#ECF0F1"),
    selectInput("drug", "Drug",
                choices=c("No Treatment","Sirolimus","Everolimus","Sirolimus + Everolimus"),
                selected="Sirolimus"),
    conditionalPanel("input.drug=='Sirolimus'",
      sliderInput("dose_siro","Sirolimus Dose (mg/day)",1,5,2,0.5)),
    conditionalPanel("input.drug=='Everolimus'",
      sliderInput("dose_ever","Everolimus Dose (mg/day)",5,15,10,1)),
    selectInput("e2_level","Estrogen Status",
                choices=c("Premenopausal","Postmenopausal","GnRH Suppressed"),
                selected="Premenopausal"),
    sliderInput("sim_wk","Simulation (weeks)",52,208,104,step=26),
    hr(),
    h5("Patient Parameters", style="padding-left:15px; color:#ECF0F1"),
    sliderInput("fev1_base","Baseline FEV1 (%pred)",40,90,72,1),
    sliderInput("dlco_base","Baseline DLCO (%pred)",30,80,60,1),
    selectInput("lam_type","LAM Type",
                choices=c("Sporadic LAM","TSC-LAM"),
                selected="Sporadic LAM"),
    actionButton("run_sim","Run Simulation",class="btn-success btn-block",
                 icon=icon("play"))
  ),

  dashboardBody(
    tabItems(

      ## ==============================================================
      ## TAB 1: Patient Profile
      ## ==============================================================
      tabItem(tabName="tab_profile",
        h2("Patient Profile & Disease Overview"),
        fluidRow(
          valueBoxOutput("vb_fev1"),
          valueBoxOutput("vb_dlco"),
          valueBoxOutput("vb_vegfd")
        ),
        fluidRow(
          box(title="Disease Summary", width=6, status="primary",
            htmlOutput("disease_summary_text")
          ),
          box(title="Severity Assessment", width=6, status="warning",
            htmlOutput("severity_text"),
            plotlyOutput("severity_gauge", height="200px")
          )
        ),
        fluidRow(
          box(title="LAM Pathophysiology", width=12, status="info",
            tags$div(style="font-size:14px; line-height:1.8",
              tags$p(tags$strong("Lymphangioleiomyomatosis (LAM)"),
                "is a rare, slowly progressive cystic lung disease that predominantly affects women of childbearing age. ",
                "It is caused by inactivating mutations in TSC1 or TSC2 genes, leading to ",
                tags$strong("mTORC1 hyperactivation"), " in LAM cells (smooth muscle-like neoplastic cells). ",
                "LAM cells proliferate, invade, and secrete proteases (MMP-2, MMP-9) that destroy lung parenchyma, ",
                "forming thin-walled cysts.",
                tags$br(), tags$br(),
                tags$strong("Key features:"),
                tags$ul(
                  tags$li("Bilateral pulmonary cysts (HRCT)"),
                  tags$li("Obstructive or mixed spirometry pattern"),
                  tags$li("Elevated serum VEGF-D (>800 pg/mL, diagnostic)"),
                  tags$li("Association with tuberous sclerosis complex (TSC-LAM, 40% of TSC patients)"),
                  tags$li("Renal angiomyolipoma (AML) in 50-60%"),
                  tags$li("Risk of pneumothorax and chylothorax")
                )
              )
            )
          )
        )
      ),

      ## ==============================================================
      ## TAB 2: Drug PK
      ## ==============================================================
      tabItem(tabName="tab_pk",
        h2("Drug Pharmacokinetics"),
        fluidRow(
          box(title="Concentration-Time Profile (First 2 Weeks)", width=8,
              status="primary", plotlyOutput("pk_plot", height="380px")),
          box(title="PK Parameters", width=4, status="info",
            tags$table(class="table table-condensed",
              tags$thead(tags$tr(tags$th("Parameter"),tags$th("Sirolimus"),tags$th("Everolimus"))),
              tags$tbody(
                tags$tr(tags$td("Dose"),tags$td("2 mg QD"),tags$td("10 mg QD")),
                tags$tr(tags$td("Bioavailability (F)"),tags$td("~15%"),tags$td("~30%")),
                tags$tr(tags$td("Vss"),tags$td("~12 L/kg"),tags$td("~5 L/kg")),
                tags$tr(tags$td("CL/F"),tags$td("~65 L/h"),tags$td("~20 L/h")),
                tags$tr(tags$td("t½"),tags$td("~62 h"),tags$td("~28 h")),
                tags$tr(tags$td("Target trough"),tags$td("5-15 ng/mL"),tags$td("5-10 ng/mL")),
                tags$tr(tags$td("Metabolism"),tags$td("CYP3A4"),tags$td("CYP3A4")),
                tags$tr(tags$td("Protein binding"),tags$td("92%"),tags$td("74%"))
              )
            ),
            tags$p(tags$em("Ref: Groth 2001; Kovarik 2001; MILES PK substudy"))
          )
        ),
        fluidRow(
          box(title="Steady-State PK Simulation", width=8, status="primary",
              plotlyOutput("pk_ss_plot", height="350px")),
          box(title="Clinical Notes", width=4, status="warning",
            tags$ul(
              tags$li("Sirolimus requires TDM due to narrow therapeutic index and high inter-patient variability"),
              tags$li("Sirolimus t½ ~62h → steady-state achieved in 2-3 weeks"),
              tags$li("Drug-drug interactions: CYP3A4 inhibitors (ketoconazole) dramatically increase levels"),
              tags$li("Food effect: high-fat meal ↑ sirolimus Cmax by ~35%"),
              tags$li("Whole blood monitoring (not plasma) due to high RBC binding")
            )
          )
        )
      ),

      ## ==============================================================
      ## TAB 3: mTOR Pathway (PD Biomarkers)
      ## ==============================================================
      tabItem(tabName="tab_pd",
        h2("mTOR Pathway Activity (Pharmacodynamics)"),
        fluidRow(
          box(title="mTORC1 Inhibition Over Time", width=6, status="primary",
              plotlyOutput("pd_mtor_plot", height="350px")),
          box(title="S6K1-pT389 and 4E-BP1 (PD Biomarkers)", width=6,
              status="info", plotlyOutput("pd_s6k_plot", height="350px"))
        ),
        fluidRow(
          box(title="Rheb-GTP Dynamics", width=6, status="warning",
              plotlyOutput("pd_rheb_plot", height="300px")),
          box(title="PD Pathway Summary", width=6, status="success",
            tags$div(style="font-size:13px",
              tags$strong("mTOR Signaling in LAM:"),
              tags$ul(
                tags$li(tags$strong("TSC2 loss → Rheb-GTP↑:"),
                  " TSC1-TSC2 complex is the GTPase-activating protein (GAP) for Rheb.",
                  " Its loss allows Rheb to remain GTP-loaded (active)"),
                tags$li(tags$strong("Rheb-GTP activates mTORC1:"),
                  " 4-fold higher mTORC1 activity in LAM cells vs normal lung"),
                tags$li(tags$strong("S6K1-pT389:"),
                  " Primary mTORC1 substrate; used as PD readout in clinical studies"),
                tags$li(tags$strong("4E-BP1:"),
                  " mTORC1 substrate; hyperphosphorylation releases eIF4E for cap-dependent translation"),
                tags$li(tags$strong("Drug effect:"),
                  " Sirolimus/everolimus: FKBP12-mediated allosteric mTORC1 inhibition.",
                  " Note: mTORC2 is largely spared, maintaining partial Akt-pS473 activity (escape mechanism)")
              )
            )
          )
        )
      ),

      ## ==============================================================
      ## TAB 4: Clinical Endpoints
      ## ==============================================================
      tabItem(tabName="tab_endpoints",
        h2("Clinical Endpoints"),
        fluidRow(
          box(title="FEV1 (% Predicted)", width=6, status="primary",
              plotlyOutput("ep_fev1_plot", height="350px")),
          box(title="DLCO (% Predicted)", width=6, status="info",
              plotlyOutput("ep_dlco_plot", height="350px"))
        ),
        fluidRow(
          box(title="Lung Cyst Volume (CT)", width=6, status="warning",
              plotlyOutput("ep_cyst_plot", height="300px")),
          box(title="6-Minute Walk Distance (Estimated)", width=6,
              status="success", plotlyOutput("ep_6mwt_plot", height="300px"))
        ),
        fluidRow(
          box(title="Clinical Outcome Milestones", width=12, status="danger",
            tableOutput("milestone_table")
          )
        )
      ),

      ## ==============================================================
      ## TAB 5: Scenario Comparison
      ## ==============================================================
      tabItem(tabName="tab_compare",
        h2("Treatment Scenario Comparison"),
        fluidRow(
          box(title="FEV1 by Scenario", width=8, status="primary",
              plotlyOutput("compare_fev1_plot", height="380px")),
          box(title="Scenarios Modeled", width=4, status="info",
            tags$ol(
              tags$li(tags$strong("Untreated:"), " Natural history (~120 mL/yr FEV1 decline)"),
              tags$li(tags$strong("Sirolimus 2mg/day:"), " Standard of care (ERS 2022)"),
              tags$li(tags$strong("Everolimus 10mg/day:"), " Alternative mTOR inhibitor"),
              tags$li(tags$strong("Sirolimus → Stop (12mo):"), " MILES off-Rx phase"),
              tags$li(tags$strong("Everolimus + GnRH:"), " Combined mTOR + estrogen suppression")
            )
          )
        ),
        fluidRow(
          box(title="VEGF-D by Scenario", width=6, status="warning",
              plotlyOutput("compare_vegfd_plot", height="300px")),
          box(title="mTORC1 Inhibition by Scenario", width=6, status="success",
              plotlyOutput("compare_mtor_plot", height="300px"))
        ),
        fluidRow(
          box(title="Outcomes Table at 12 and 24 Months", width=12,
              DTOutput("scenario_table"))
        )
      ),

      ## ==============================================================
      ## TAB 6: Biomarker Dashboard
      ## ==============================================================
      tabItem(tabName="tab_biomarker",
        h2("Biomarker Dashboard"),
        fluidRow(
          valueBoxOutput("vb_vegfd_curr"),
          valueBoxOutput("vb_s6k1"),
          valueBoxOutput("vb_mmp")
        ),
        fluidRow(
          box(title="Serum VEGF-D (Diagnostic Biomarker)", width=6, status="primary",
              plotlyOutput("bio_vegfd_plot", height="320px")),
          box(title="S6K1-pT389 (PD Biomarker)", width=6, status="info",
              plotlyOutput("bio_s6k_plot", height="320px"))
        ),
        fluidRow(
          box(title="MMP Activity (Protease / Cyst Formation)", width=6,
              status="warning", plotlyOutput("bio_mmp_plot", height="300px")),
          box(title="Renal AML Volume", width=6, status="success",
              plotlyOutput("bio_aml_plot", height="300px"))
        ),
        fluidRow(
          box(title="Biomarker Reference Values", width=12,
            tags$table(class="table table-bordered",
              tags$thead(tags$tr(
                tags$th("Biomarker"), tags$th("Normal"), tags$th("LAM"),
                tags$th("Diagnostic?"), tags$th("Ref")
              )),
              tags$tbody(
                tags$tr(tags$td("Serum VEGF-D"),tags$td("<500 pg/mL"),
                        tags$td(">800 pg/mL"),tags$td("Yes (90% specificity)"),
                        tags$td("Young 2011")),
                tags$tr(tags$td("S6K1-pT389"),tags$td("Low (PBMCs)"),
                        tags$td("Elevated"),tags$td("Research only"),
                        tags$td("Goncharova 2011")),
                tags$tr(tags$td("MMP-2/MMP-9"),tags$td("Low"),
                        tags$td("Elevated in BAL/serum"),tags$td("Research"),
                        tags$td("Seyama 2006")),
                tags$tr(tags$td("Urine LAM cell clusters"),
                        tags$td("Absent"),tags$td("Present (30-40%)"),
                        tags$td("Suggestive"),tags$td("Schiavina 2007")),
                tags$tr(tags$td("CT cyst score"),tags$td("0 (absent)"),
                        tags$td("Warwick 1-4"),tags$td("Diagnostic (with VEGF-D)"),
                        tags$td("Johnson 2010"))
              )
            )
          )
        )
      )

    ) # end tabItems
  )
)

## ====================================================================
## SERVER
## ====================================================================
server <- function(input, output, session) {

  ## Reactive simulation
  sim_data <- eventReactive(input$run_sim, {
    dose_val <- if(input$drug=="Sirolimus") input$dose_siro else input$dose_ever

    run_sim(
      mod       = lam_mod,
      drug      = input$drug,
      dose      = dose_val,
      duration_wk = input$sim_wk,
      E2_level  = input$e2_level
    )
  }, ignoreNULL=FALSE)

  ## All 5 scenarios (for comparison tab)
  all_scenarios <- reactive({
    wk <- input$sim_wk
    e2_map <- c("Premenopausal"=1.0,"Postmenopausal"=0.2,"GnRH Suppressed"=0.05)
    E2_val <- e2_map[input$e2_level]

    # Scenario 1: Untreated
    s1 <- lam_mod %>% param(E2_basal=E2_val) %>%
      mrgsim(end=wk*168, delta=168) %>% as_tibble() %>%
      mutate(scenario="1_Untreated", time_weeks=time/168)

    # Scenario 2: Sirolimus 2mg
    ev2 <- ev(amt=2,cmt=1,ii=24,addl=wk*7-1)
    s2 <- lam_mod %>% param(E2_basal=E2_val) %>% ev(ev2) %>%
      mrgsim(end=wk*168, delta=168) %>% as_tibble() %>%
      mutate(scenario="2_Sirolimus_2mgQD", time_weeks=time/168)

    # Scenario 3: Everolimus 10mg
    ev3 <- ev(amt=10,cmt=4,ii=24,addl=wk*7-1)
    s3 <- lam_mod %>% param(E2_basal=E2_val) %>% ev(ev3) %>%
      mrgsim(end=wk*168, delta=168) %>% as_tibble() %>%
      mutate(scenario="3_Everolimus_10mgQD", time_weeks=time/168)

    # Scenario 4: Sirolimus 12 months, then stop
    ev4 <- ev(amt=2,cmt=1,ii=24,addl=52*7-1)
    s4 <- lam_mod %>% param(E2_basal=E2_val) %>% ev(ev4) %>%
      mrgsim(end=wk*168, delta=168) %>% as_tibble() %>%
      mutate(scenario="4_Siro_12mo_Stop", time_weeks=time/168)

    # Scenario 5: Everolimus + GnRH
    ev5 <- ev(amt=10,cmt=4,ii=24,addl=wk*7-1)
    s5 <- lam_mod %>% param(E2_basal=0.05) %>% ev(ev5) %>%
      mrgsim(end=wk*168, delta=168) %>% as_tibble() %>%
      mutate(scenario="5_Everolimus_GnRH", time_weeks=time/168)

    bind_rows(s1,s2,s3,s4,s5)
  })

  ## Value boxes (Tab 1)
  output$vb_fev1 <- renderValueBox({
    d <- sim_data()
    val <- round(tail(d$FEV1_pct,1),1)
    col <- if(val<30)"red" else if(val<50)"orange" else "blue"
    valueBox(paste0(val,"%"), "FEV1 (% predicted)", icon=icon("lungs"), color=col)
  })
  output$vb_dlco <- renderValueBox({
    d <- sim_data(); val <- round(tail(d$DLCO_pct,1),1)
    col <- if(val<40)"red" else if(val<55)"orange" else "green"
    valueBox(paste0(val,"%"), "DLCO (% predicted)", icon=icon("wind"), color=col)
  })
  output$vb_vegfd <- renderValueBox({
    d <- sim_data(); val <- round(d$VEGFD_pgmL[1],0)
    col <- if(val>800)"red" else "green"
    valueBox(paste0(val," pg/mL"), "Serum VEGF-D", icon=icon("tint"), color=col)
  })

  output$disease_summary_text <- renderUI({
    tags$div(
      tags$p(tags$strong("Type:"), input$lam_type),
      tags$p(tags$strong("Estrogen Status:"), input$e2_level),
      tags$p(tags$strong("Simulation Duration:"), paste(input$sim_wk,"weeks")),
      tags$p(tags$strong("Treatment:"), input$drug),
      tags$p(tags$strong("LAM Cell Burden (final):"),
             round(tail(sim_data()$LAM_burden,1),2), "× baseline")
    )
  })

  output$severity_text <- renderUI({
    d <- sim_data()
    fev1 <- tail(d$FEV1_pct,1)
    sev <- if(fev1>70)"Mild" else if(fev1>50)"Moderate" else if(fev1>30)"Severe" else "Very Severe"
    tags$p(tags$strong("LAM Severity: "), tags$span(sev, style=paste0(
      "color:",if(sev=="Mild")"green" else if(sev=="Moderate")"orange" else "red")))
  })

  output$severity_gauge <- renderPlotly({
    d <- sim_data(); fev1 <- tail(d$FEV1_pct,1)
    plot_ly(type="indicator", mode="gauge+number",
            value=fev1,
            title=list(text="FEV1 % predicted"),
            gauge=list(
              axis=list(range=list(0,100)),
              bar=list(color="steelblue"),
              steps=list(
                list(range=c(0,30),  color="#FF4444"),
                list(range=c(30,50), color="#FF9944"),
                list(range=c(50,70), color="#FFDD44"),
                list(range=c(70,100),color="#44CC44")
              )
            )) %>% layout(margin=list(t=30,b=0))
  })

  ## PK plots (Tab 2)
  output$pk_plot <- renderPlotly({
    d_pk <- sim_data() %>% filter(time_weeks <= 2)
    p <- ggplot(d_pk, aes(x=time, y=Cp_siro_ngml+Cp_ever_ngml)) +
      geom_line(color="steelblue",size=1.2) +
      geom_hline(yintercept=5, linetype="dashed", color="darkgreen", alpha=0.8) +
      geom_hline(yintercept=15, linetype="dashed", color="darkred", alpha=0.8) +
      labs(x="Hours", y="Blood Concentration (ng/mL)",
           title="Concentration-Time Profile") +
      theme_bw()
    ggplotly(p)
  })

  output$pk_ss_plot <- renderPlotly({
    d_ss <- sim_data() %>% filter(time_weeks >= 4, time_weeks <= 8)
    p <- ggplot(d_ss, aes(x=time_weeks*168, y=Cp_siro_ngml+Cp_ever_ngml)) +
      geom_line(color="darkorange",size=1.2) +
      geom_ribbon(aes(ymin=5, ymax=15), fill="green", alpha=0.1) +
      labs(x="Hours", y="Concentration (ng/mL)", title="Steady-State PK (Weeks 4-8)") +
      theme_bw()
    ggplotly(p)
  })

  ## PD plots (Tab 3)
  output$pd_mtor_plot <- renderPlotly({
    d <- sim_data()
    p <- ggplot(d, aes(x=time_weeks, y=mTORC1_inhib_pct)) +
      geom_line(color="darkblue",size=1.2) +
      labs(x="Weeks", y="mTORC1 Inhibition (%)", title="mTORC1 Inhibition") +
      ylim(0,100) + theme_bw()
    ggplotly(p)
  })

  output$pd_s6k_plot <- renderPlotly({
    d <- sim_data() %>%
      select(time_weeks, S6K1_phos, EBP1_phos) %>%
      tidyr::pivot_longer(cols=c(S6K1_phos,EBP1_phos),
                          names_to="marker", values_to="value")
    p <- ggplot(d, aes(x=time_weeks, y=value, color=marker)) +
      geom_line(size=1.2) +
      scale_color_manual(values=c("S6K1_phos"="steelblue","EBP1_phos"="darkorange")) +
      labs(x="Weeks", y="Phosphorylation (normalized)", title="S6K1-pT389 & 4E-BP1") +
      theme_bw()
    ggplotly(p)
  })

  output$pd_rheb_plot <- renderPlotly({
    d <- sim_data()
    p <- ggplot(d, aes(x=time_weeks, y=RHEB_GTP)) +
      geom_line(color="purple",size=1.2) +
      geom_hline(yintercept=1, linetype="dashed", color="gray") +
      labs(x="Weeks", y="Rheb-GTP (normalized)", title="Rheb-GTP Dynamics") +
      theme_bw()
    ggplotly(p)
  })

  ## Clinical endpoint plots (Tab 4)
  output$ep_fev1_plot <- renderPlotly({
    d <- sim_data()
    p <- ggplot(d, aes(x=time_weeks, y=FEV1_pct)) +
      geom_line(color="steelblue",size=1.2) +
      geom_hline(yintercept=30, linetype="dashed", color="red") +
      geom_ribbon(aes(ymin=30, ymax=50), fill="red", alpha=0.05) +
      annotate("text",x=5,y=32,label="Transplant threshold",size=3,color="red") +
      labs(x="Weeks", y="FEV1 (% predicted)") + theme_bw()
    ggplotly(p)
  })

  output$ep_dlco_plot <- renderPlotly({
    d <- sim_data()
    p <- ggplot(d, aes(x=time_weeks, y=DLCO_pct)) +
      geom_line(color="darkorange",size=1.2) +
      labs(x="Weeks", y="DLCO (% predicted)") + theme_bw()
    ggplotly(p)
  })

  output$ep_cyst_plot <- renderPlotly({
    d <- sim_data()
    p <- ggplot(d, aes(x=time_weeks, y=CystVol)) +
      geom_line(color="darkred",size=1.2) +
      labs(x="Weeks", y="Cyst Volume (% lung)", title="CT Cyst Volume") +
      theme_bw()
    ggplotly(p)
  })

  output$ep_6mwt_plot <- renderPlotly({
    d <- sim_data() %>%
      mutate(six_mwt_est = pmax(100, 450 - (100 - FEV1_pct)*3.5 - (60-DLCO_pct)*2))
    p <- ggplot(d, aes(x=time_weeks, y=six_mwt_est)) +
      geom_line(color="green4",size=1.2) +
      labs(x="Weeks", y="Estimated 6MWD (m)", title="6-Min Walk Distance (Estimated)") +
      theme_bw()
    ggplotly(p)
  })

  output$milestone_table <- renderTable({
    d <- sim_data()
    t_fev50 <- d %>% filter(FEV1_pct <= 50) %>% slice(1) %>% pull(time_weeks)
    t_fev30 <- d %>% filter(FEV1_pct <= 30) %>% slice(1) %>% pull(time_weeks)
    t_dlco40 <- d %>% filter(DLCO_pct <= 40) %>% slice(1) %>% pull(time_weeks)
    data.frame(
      Milestone=c("FEV1 ≤ 50% predicted","FEV1 ≤ 30% (transplant)","DLCO ≤ 40%"),
      `Time (weeks)`=c(
        ifelse(length(t_fev50)>0,round(t_fev50,0),"Not reached"),
        ifelse(length(t_fev30)>0,round(t_fev30,0),"Not reached"),
        ifelse(length(t_dlco40)>0,round(t_dlco40,0),"Not reached")
      )
    )
  }, rownames=FALSE)

  ## Scenario comparison plots (Tab 5)
  sc_cols <- c("1_Untreated"="#555555","2_Sirolimus_2mgQD"="#1F77B4",
               "3_Everolimus_10mgQD"="#FF7F0E","4_Siro_12mo_Stop"="#9467BD",
               "5_Everolimus_GnRH"="#2CA02C")

  output$compare_fev1_plot <- renderPlotly({
    d <- all_scenarios()
    p <- ggplot(d, aes(x=time_weeks, y=FEV1_pct, color=scenario)) +
      geom_line(size=1.1) +
      geom_hline(yintercept=30, linetype="dashed", color="red", alpha=0.7) +
      scale_color_manual(values=sc_cols, name="Scenario") +
      labs(x="Weeks", y="FEV1 (% predicted)", title="FEV1 by Treatment Scenario") +
      theme_bw()
    ggplotly(p)
  })

  output$compare_vegfd_plot <- renderPlotly({
    d <- all_scenarios()
    p <- ggplot(d, aes(x=time_weeks, y=VEGFD_pgmL, color=scenario)) +
      geom_line(size=1.1) +
      geom_hline(yintercept=800, linetype="dashed", color="darkred") +
      scale_color_manual(values=sc_cols) +
      labs(x="Weeks", y="VEGF-D (pg/mL)") + theme_bw()
    ggplotly(p)
  })

  output$compare_mtor_plot <- renderPlotly({
    d <- all_scenarios()
    p <- ggplot(d, aes(x=time_weeks, y=mTORC1_inhib_pct, color=scenario)) +
      geom_line(size=1.1) +
      scale_color_manual(values=sc_cols) +
      labs(x="Weeks", y="mTORC1 Inhibition (%)") + theme_bw()
    ggplotly(p)
  })

  output$scenario_table <- renderDT({
    d <- all_scenarios()
    tbl <- d %>%
      filter(abs(time_weeks-52)<2 | abs(time_weeks-104)<2) %>%
      mutate(timepoint=ifelse(abs(time_weeks-52)<2,"12 months","24 months")) %>%
      group_by(scenario, timepoint) %>%
      summarise(
        `FEV1 (%pred)` = round(mean(FEV1_pct),1),
        `DLCO (%pred)` = round(mean(DLCO_pct),1),
        `VEGF-D (pg/mL)` = round(mean(VEGFD_pgmL),0),
        `mTORC1 Inh (%)` = round(mean(mTORC1_inhib_pct),1),
        `AML (mL)` = round(mean(AML_mL),0),
        .groups="drop"
      )
    datatable(tbl, rownames=FALSE, options=list(pageLength=10)) %>%
      formatStyle("FEV1 (%pred)", backgroundColor=styleInterval(c(30,50),
                  c("#FF6666","#FFAA66","#AAFFAA")))
  })

  ## Biomarker dashboard (Tab 6)
  output$vb_vegfd_curr <- renderValueBox({
    d <- sim_data(); val <- round(tail(d$VEGFD_pgmL,1),0)
    col <- if(val>800)"red" else "green"
    valueBox(paste0(val," pg/mL"),"VEGF-D (final)",icon=icon("tint"),color=col)
  })
  output$vb_s6k1 <- renderValueBox({
    d <- sim_data(); val <- round(tail(d$S6K1_phos,1),2)
    valueBox(val,"S6K1-pT389 (norm.)",icon=icon("atom"),color="blue")
  })
  output$vb_mmp <- renderValueBox({
    d <- sim_data(); val <- round(tail(d$MMP_normalized,1),2)
    valueBox(val,"MMP Activity (norm.)",icon=icon("cut"),color="orange")
  })

  output$bio_vegfd_plot <- renderPlotly({
    d <- sim_data()
    p <- ggplot(d, aes(x=time_weeks, y=VEGFD_pgmL)) +
      geom_line(color="steelblue",size=1.2) +
      geom_hline(yintercept=800, linetype="dashed", color="red") +
      annotate("text",x=5,y=850,label=">800 pg/mL = Diagnostic",size=3,color="red") +
      labs(x="Weeks",y="VEGF-D (pg/mL)",title="Serum VEGF-D") +
      theme_bw()
    ggplotly(p)
  })

  output$bio_s6k_plot <- renderPlotly({
    d <- sim_data()
    p <- ggplot(d, aes(x=time_weeks, y=S6K1_phos)) +
      geom_line(color="purple",size=1.2) +
      labs(x="Weeks",y="S6K1-pT389 (normalized)",title="S6K1 Phosphorylation (PD Biomarker)") +
      theme_bw()
    ggplotly(p)
  })

  output$bio_mmp_plot <- renderPlotly({
    d <- sim_data()
    p <- ggplot(d, aes(x=time_weeks, y=MMP_normalized)) +
      geom_line(color="darkorange",size=1.2) +
      labs(x="Weeks",y="MMP Activity (normalized)",title="Matrix Metalloproteinase Activity") +
      theme_bw()
    ggplotly(p)
  })

  output$bio_aml_plot <- renderPlotly({
    d <- sim_data()
    p <- ggplot(d, aes(x=time_weeks, y=AML_mL)) +
      geom_line(color="brown",size=1.2) +
      geom_hline(yintercept=30, linetype="dashed", color="red") +
      annotate("text",x=5,y=35,label="Embolization threshold (>30mm)",size=3) +
      labs(x="Weeks",y="AML Volume (mL)",title="Renal Angiomyolipoma Volume") +
      theme_bw()
    ggplotly(p)
  })
}

## ====================================================================
## LAUNCH APP
## ====================================================================
shinyApp(ui=ui, server=server)
