## =============================================================================
## Polymyositis QSP — Interactive Shiny Dashboard
## =============================================================================
## Disease: Polymyositis (CD8+ T-cell Mediated Inflammatory Myopathy)
## Tabs  : Patient Profile · PK · Immunology · Muscle Pathology ·
##         Clinical Endpoints · Scenario Comparison · Biomarker Trajectories
## =============================================================================

library(shiny)
library(shinydashboard)
library(mrgsolve)
library(dplyr)
library(ggplot2)
library(plotly)
library(DT)
library(patchwork)

## ── Inline mrgsolve model (abbreviated for Shiny embedding) ─────────────────
PM_CODE <- '
$PARAM
  PRED_KA=1.5, PRED_F=0.82, PRED_CL=15.0, PRED_V1=45.0, PRED_Q=8.0, PRED_V2=90.0, KHB=0.7,
  MTX_KA=0.9,  MTX_F=0.72,  MTX_CL=5.2,  MTX_V1=28.0, MTX_KPG=0.05, MTX_KPGOUT=0.003,
  AZA_KA=0.8,  AZA_F=0.60,  AZA_K12=0.25,AZA_6MP_CL=2.5,AZA_V2=40.0, TPMT_act=0.6,
  RTX_CL=0.008,RTX_V1=3.5,  RTX_Q=0.025, RTX_V2=5.0,  RTX_KON=0.60, RTX_KOFF=0.01,
  RTX_KINT=0.03,RTX_CD20_0=150.0,
  JAKI_KA=1.2, JAKI_F=0.79, JAKI_CL=5.8, JAKI_V1=76.0, JAKI_Q=3.5,  JAKI_V2=40.0,
  CD8_PROD=5.0,CD8_KDEATH=0.02,CD8_KACT=0.15,CD8_KEXP=0.8,CD8_KMIG=0.20,CD8_KDEATH_E=0.08,
  CD4_PROD=4.0,CD4_KACT=0.10,CD4_KDEATH=0.05,TH1_IFNg_PROD=0.8,
  B_PROD=2.0,B_KACT=0.05,B_KDEATH=0.02,PC_PROD_RATE=0.10,PC_KDEATH=0.005,
  AB_PROD=0.5,AB_CL=0.015,
  IFNG_PROD_0=0.08,IFNG_KD=0.30,TNF_PROD_0=0.06,TNF_KD=0.50,
  IL6_PROD_0=0.10,IL6_KD=0.60,
  MUS_INF_KMAX=0.20,MUS_INF_KD=0.05,MUS_INF_EC50=15.0,
  MHC_IFNG_EC50=2.0,MHCI_KMAX=5.0,
  CK_BASE=120.0,CK_KPROD=0.10,CK_KD=0.15,CK_MAX_FACTOR=400.0,
  MMT_MAX=80.0,MMT_KDECLINE=0.008,MMT_KREC=0.015,
  PRED_IC50=2.0,PRED_IMAX=0.90,PRED_HILL=1.5,
  MTX_IC50=0.05,MTX_IMAX=0.75,
  TGN_IC50=180.0,TGN_IMAX=0.70,
  RTX_DEPL_EC50=0.5,RTX_DEPL_IMAX=0.97,
  JAKI_IC50=15.0,JAKI_IMAX=0.85

$CMT
  PRED_GUT PRED_C1 PRED_C2
  MTX_GUT MTX_C1 MTX_PG
  AZA_GUT SixTGN
  RTX_C1 RTX_C2 RTX_CD20 RTX_BOUND
  JAKI_C1 JAKI_C2
  CD8N CD8E CD4TH1 BCELL PLASMA_C AUTOAB
  IFNG TNFa IL6C MHCI MUS_INF CK_S MMT8

$INIT
  PRED_GUT=0, PRED_C1=0, PRED_C2=0,
  MTX_GUT=0, MTX_C1=0, MTX_PG=0,
  AZA_GUT=0, SixTGN=0,
  RTX_C1=0, RTX_C2=0, RTX_CD20=150, RTX_BOUND=0,
  JAKI_C1=0, JAKI_C2=0,
  CD8N=300, CD8E=80, CD4TH1=250, BCELL=120, PLASMA_C=50, AUTOAB=45,
  IFNG=18, TNFa=8, IL6C=12, MHCI=6.5, MUS_INF=7.0, CK_S=4800, MMT8=38

$MAIN
  double PREDNIS = KHB * PRED_C1 / PRED_V1 * 1000;
  double MTX_PG_conc = MTX_PG;
  double TGN_conc = SixTGN * 10.0;
  double JAKI_Cp = JAKI_C1 / JAKI_V1 * 1e6;
  double E_PRED = PRED_IMAX * pow(PREDNIS,PRED_HILL)/(pow(PRED_IC50,PRED_HILL)+pow(PREDNIS,PRED_HILL));
  double E_MTX  = MTX_IMAX * MTX_PG_conc/(MTX_IC50 + MTX_PG_conc);
  double E_TGN  = TGN_IMAX * TGN_conc/(TGN_IC50 + TGN_conc);
  double E_JAKI = JAKI_IMAX * JAKI_Cp/(JAKI_IC50 + JAKI_Cp);
  double RTX_RC = RTX_BOUND;
  double E_RTX_B = RTX_DEPL_IMAX * RTX_RC/(RTX_DEPL_EC50 + RTX_RC);
  double E_COMBO_CD8 = 1.0-(1.0-E_PRED)*(1.0-E_MTX)*(1.0-E_TGN);
  double E_COMBO_CYT = 1.0-(1.0-E_PRED)*(1.0-E_JAKI);
  double E_COMBO_B   = 1.0-(1.0-E_RTX_B)*(1.0-E_PRED)*(1.0-E_MTX);
  double MHCI_eq = 1.0+(MHCI_KMAX-1.0)*IFNG/(MHC_IFNG_EC50+IFNG);
  double Drug_suppr = (E_PRED*0.5+E_MTX*0.15+E_TGN*0.10+E_RTX_B*0.15+E_JAKI*0.20);

$ODE
  dxdt_PRED_GUT=-PRED_KA*PRED_GUT;
  dxdt_PRED_C1= PRED_KA*PRED_GUT-(PRED_CL+PRED_Q)/PRED_V1*PRED_C1+PRED_Q/PRED_V2*PRED_C2;
  dxdt_PRED_C2= PRED_Q/PRED_V1*PRED_C1-PRED_Q/PRED_V2*PRED_C2;
  dxdt_MTX_GUT=-MTX_KA*MTX_GUT;
  dxdt_MTX_C1= MTX_KA*MTX_GUT-(MTX_CL/28.0+MTX_KPG)*MTX_C1;
  dxdt_MTX_PG= MTX_KPG*MTX_C1-MTX_KPGOUT*MTX_PG;
  dxdt_AZA_GUT=-AZA_KA*AZA_GUT;
  double AZA_to_6TGN=AZA_K12*AZA_GUT*(1.0-0.3*TPMT_act);
  dxdt_SixTGN= AZA_to_6TGN-(AZA_6MP_CL/AZA_V2)*SixTGN;
  dxdt_RTX_C1=-(RTX_CL+RTX_Q)/RTX_V1*RTX_C1+RTX_Q/RTX_V2*RTX_C2
              -RTX_KON*RTX_C1/RTX_V1*RTX_CD20*RTX_V1+RTX_KOFF*RTX_BOUND;
  dxdt_RTX_C2= RTX_Q/RTX_V1*RTX_C1-RTX_Q/RTX_V2*RTX_C2;
  dxdt_RTX_CD20=-RTX_KON*(RTX_C1/RTX_V1)*RTX_CD20+RTX_KOFF*RTX_BOUND
               +RTX_KINT*(RTX_CD20_0-RTX_CD20);
  dxdt_RTX_BOUND=RTX_KON*(RTX_C1/RTX_V1)*RTX_CD20-(RTX_KOFF+RTX_KINT)*RTX_BOUND;
  dxdt_JAKI_C1=-(JAKI_CL+JAKI_Q)/JAKI_V1*JAKI_C1+JAKI_Q/JAKI_V2*JAKI_C2;
  dxdt_JAKI_C2= JAKI_Q/JAKI_V1*JAKI_C1-JAKI_Q/JAKI_V2*JAKI_C2;
  double MHCI_stim=MHCI/(1.0+MHCI);
  dxdt_CD8N= CD8_PROD-CD8_KDEATH*CD8N-CD8_KACT*CD8N*MHCI_stim*(1.0-E_COMBO_CD8);
  dxdt_CD8E= CD8_KACT*CD8N*MHCI_stim*(1.0-E_COMBO_CD8)*CD8_KEXP-CD8_KMIG*CD8E-CD8_KDEATH_E*CD8E;
  dxdt_CD4TH1= CD4_PROD*(1.0+IFNG/(2.0+IFNG))-CD4_KDEATH*CD4TH1-CD4_KACT*E_COMBO_CD8*CD4TH1;
  dxdt_BCELL= B_PROD-B_KDEATH*BCELL-B_KACT*(1.0-E_COMBO_B)*BCELL-E_RTX_B*B_KACT*BCELL;
  dxdt_PLASMA_C= PC_PROD_RATE*BCELL*(1.0-E_COMBO_B)-PC_KDEATH*PLASMA_C-E_RTX_B*0.3*PLASMA_C;
  dxdt_AUTOAB= AB_PROD*PLASMA_C*(1.0-E_COMBO_B)-AB_CL*AUTOAB;
  dxdt_IFNG= IFNG_PROD_0+TH1_IFNg_PROD*CD4TH1/100.0+0.05*CD8E-IFNG_KD*IFNG
             -E_COMBO_CYT*IFNG_PROD_0*5.0;
  dxdt_TNFa= TNF_PROD_0*(1.0+MUS_INF)+0.01*CD4TH1-TNF_KD*TNFa-E_COMBO_CYT*TNF_PROD_0*4.0;
  dxdt_IL6C= IL6_PROD_0*(1.0+0.5*MUS_INF)+0.008*CD4TH1-IL6_KD*IL6C-E_COMBO_CYT*IL6_PROD_0*5.0;
  dxdt_MHCI= 0.5*(MHCI_eq-MHCI);
  double CD8E_stim=MUS_INF_KMAX*CD8E/(MUS_INF_EC50+CD8E);
  double AB_stim=0.02*AUTOAB/(20.0+AUTOAB);
  dxdt_MUS_INF= CD8E_stim+AB_stim-MUS_INF_KD*MUS_INF*(1.0+Drug_suppr);
  dxdt_CK_S= CK_KPROD*MUS_INF*CK_BASE-CK_KD*(CK_S-CK_BASE);
  dxdt_MMT8= -MMT_KDECLINE*MUS_INF*MMT8+MMT_KREC*(MMT_MAX-MMT8)*Drug_suppr;

$TABLE
  double Cpred=KHB*PRED_C1/PRED_V1*1000;
  double Cjaki=JAKI_C1/JAKI_V1*1e6;
  double CK_fold=CK_S/120.0;
  double MMT8_pct=MMT8/80.0*100;
  double TIS_approx=(1.0-MUS_INF/10.0)*100;
  int Remission=(CK_S/200.0<1.5 && MMT8>72)?1:0;

$CAPTURE
  Cpred Cjaki CD8E CD4TH1 BCELL PLASMA_C AUTOAB
  IFNG TNFa IL6C MHCI MUS_INF CK_S MMT8
  CK_fold MMT8_pct TIS_approx Remission
'

pm_mod <- mcode("pm_shiny", PM_CODE)

## ── Run scenario function ────────────────────────────────────────────────────
run_scenario <- function(mod, scenario, pred_dose, mtx_dose, aza_dose,
                          rtx_dose, jaki_dose, duration_d,
                          init_ck, init_mmt, init_ifng) {
  m <- mod %>%
    init(CK_S = init_ck, MMT8 = init_mmt, IFNG = init_ifng,
         MUS_INF = 7 * init_ck / 4800,
         CD8E    = 80 * init_ck / 4800,
         AUTOAB  = 45 * init_ck / 4800,
         MHCI    = 6.5 * init_ifng / 18)

  evs <- NULL
  if (pred_dose > 0 && "Prednisone" %in% scenario) {
    evs <- c(evs, list(ev(amt = pred_dose, ii = 24, addl = duration_d - 1, cmt = 1, evid = 1)))
  }
  if (mtx_dose > 0 && "MTX" %in% scenario) {
    evs <- c(evs, list(ev(amt = mtx_dose, ii = 168, addl = ceiling(duration_d/7)-1, cmt = 4, evid = 1)))
  }
  if (aza_dose > 0 && "AZA" %in% scenario) {
    evs <- c(evs, list(ev(amt = aza_dose, ii = 24, addl = duration_d - 1, cmt = 7, evid = 1)))
  }
  if (rtx_dose > 0 && "Rituximab" %in% scenario) {
    evs <- c(evs, list(
      ev(amt = rtx_dose, cmt = 9, evid = 1, time = 0, rate = -2),
      ev(amt = rtx_dose, cmt = 9, evid = 1, time = 14*24, rate = -2)
    ))
  }
  if (jaki_dose > 0 && "JAKi" %in% scenario) {
    evs <- c(evs, list(ev(amt = jaki_dose, ii = 24, addl = duration_d - 1, cmt = 13, evid = 1)))
  }

  if (!is.null(evs)) {
    e_combined <- do.call(bind, evs)
    out <- m %>% ev(e_combined) %>% mrgsim(end = duration_d*24, delta = 12) %>% as.data.frame()
  } else {
    out <- m %>% mrgsim(end = duration_d*24, delta = 12) %>% as.data.frame()
  }
  out$time_d <- out$time / 24
  out
}

## ════════════════════════════════════════════════════════════════════════════
## UI
## ════════════════════════════════════════════════════════════════════════════
ui <- dashboardPage(
  skin = "blue",
  dashboardHeader(title = "Polymyositis QSP Dashboard"),

  dashboardSidebar(
    sidebarMenu(
      menuItem("Patient Profile",     tabName = "tab_patient",  icon = icon("user-md")),
      menuItem("PK Profiles",         tabName = "tab_pk",       icon = icon("flask")),
      menuItem("Immunology",          tabName = "tab_immune",   icon = icon("dna")),
      menuItem("Muscle Pathology",    tabName = "tab_muscle",   icon = icon("running")),
      menuItem("Clinical Endpoints",  tabName = "tab_clinical", icon = icon("chart-line")),
      menuItem("Scenario Comparison", tabName = "tab_compare",  icon = icon("balance-scale")),
      menuItem("Biomarker Panel",     tabName = "tab_biomark",  icon = icon("vials"))
    )
  ),

  dashboardBody(
    tags$head(tags$style(HTML("
      .content-wrapper { background-color: #F8F9FA; }
      .box { border-radius: 8px; box-shadow: 0 2px 6px rgba(0,0,0,0.1); }
      .value-box { border-radius: 8px; }
    "))),

    tabItems(
      ## ── TAB 1: Patient Profile ──────────────────────────────────────────
      tabItem(tabName = "tab_patient",
        fluidRow(
          box(title = "Disease Characteristics", width = 4, status = "primary",
            selectInput("autoab_type", "Autoantibody Type:",
              choices = c("Anti-Jo-1 (Anti-synthetase)" = "jo1",
                          "Anti-SRP (Signal Recognition Particle)" = "srp",
                          "Anti-MDA5 (Clinically Amyopathic)" = "mda5",
                          "Anti-Mi-2 (Classic PM/DM)" = "mi2",
                          "Seronegative" = "sero_neg"),
              selected = "jo1"),
            sliderInput("init_ck", "Initial Serum CK (U/L):", 500, 30000, 4800, 100),
            sliderInput("init_mmt", "Initial MMT-8 Score:", 0, 80, 38, 1),
            sliderInput("init_ifng", "Initial IFN-γ (pg/mL):", 0.5, 50, 18, 0.5),
            numericInput("weight", "Body Weight (kg):", 70, 40, 120, 5),
            numericInput("age", "Age (years):", 45, 18, 80, 1)
          ),
          box(title = "Treatment Selection", width = 4, status = "warning",
            checkboxGroupInput("scenario", "Active Treatments:",
              choices = c("Prednisone", "MTX", "AZA", "Rituximab", "JAKi"),
              selected = c("Prednisone", "MTX")),
            sliderInput("pred_dose",  "Prednisone dose (mg/day):", 5, 80, 60, 5),
            sliderInput("mtx_dose",   "MTX dose (mg/week):", 5, 25, 15, 2.5),
            sliderInput("aza_dose",   "AZA dose (mg/day):", 50, 250, 150, 25),
            sliderInput("rtx_dose",   "Rituximab dose (mg):", 500, 1000, 1000, 250),
            sliderInput("jaki_dose",  "JAKi (Baricitinib) dose (mg/day):", 2, 8, 4, 2),
            sliderInput("duration",   "Simulation duration (days):", 90, 730, 365, 30),
            actionButton("run_sim", "Run Simulation", icon = icon("play"),
                         class = "btn-success btn-lg", width = "100%")
          ),
          box(title = "PM Subtype Guide", width = 4, status = "info",
            h4("Autoantibody Clinical Correlations"),
            tableOutput("autoab_table"),
            hr(),
            h4("Disease Activity Classification"),
            tableOutput("das_table")
          )
        ),
        fluidRow(
          valueBoxOutput("vbox_ck",   width = 3),
          valueBoxOutput("vbox_mmt",  width = 3),
          valueBoxOutput("vbox_ifng", width = 3),
          valueBoxOutput("vbox_remission", width = 3)
        )
      ),

      ## ── TAB 2: PK Profiles ──────────────────────────────────────────────
      tabItem(tabName = "tab_pk",
        fluidRow(
          box(title = "Prednisolone Concentration-Time", width = 6, status = "primary",
            plotlyOutput("plot_pk_pred", height = "350px")),
          box(title = "MTX Polyglutamate Accumulation", width = 6, status = "warning",
            plotlyOutput("plot_pk_mtx", height = "350px"))
        ),
        fluidRow(
          box(title = "JAK Inhibitor (Baricitinib) PK", width = 6, status = "info",
            plotlyOutput("plot_pk_jaki", height = "350px")),
          box(title = "Rituximab — TMDD Model", width = 6, status = "danger",
            plotlyOutput("plot_pk_rtx", height = "350px"))
        ),
        fluidRow(
          box(title = "PK Summary Table", width = 12, status = "success",
            DTOutput("pk_table"))
        )
      ),

      ## ── TAB 3: Immunology ───────────────────────────────────────────────
      tabItem(tabName = "tab_immune",
        fluidRow(
          box(title = "CD8+ Effector T Cells in Muscle", width = 6, status = "danger",
            plotlyOutput("plot_cd8", height = "320px"),
            p("CD8+ T cells drive myofiber necrosis via perforin/granzyme B and FasL.")),
          box(title = "Th1 CD4+ T Cells", width = 6, status = "warning",
            plotlyOutput("plot_th1", height = "320px"),
            p("Th1 cells produce IFN-γ, which upregulates MHC-I on myofibers."))
        ),
        fluidRow(
          box(title = "B Cells & Plasma Cells", width = 6, status = "info",
            plotlyOutput("plot_bcell", height = "320px"),
            p("B cell depletion by rituximab takes weeks to reflect in autoantibody levels.")),
          box(title = "Autoantibody Level", width = 6, status = "primary",
            plotlyOutput("plot_autoab", height = "320px"),
            p("Autoantibodies (anti-Jo-1, SRP, MDA5) contribute to immune complex deposition."))
        )
      ),

      ## ── TAB 4: Muscle Pathology ─────────────────────────────────────────
      tabItem(tabName = "tab_muscle",
        fluidRow(
          box(title = "Muscle Inflammation Index", width = 6, status = "danger",
            plotlyOutput("plot_inf", height = "320px")),
          box(title = "MHC-I Expression on Myofibers", width = 6, status = "warning",
            plotlyOutput("plot_mhci", height = "320px"))
        ),
        fluidRow(
          box(title = "IFN-γ (Key Driver of MHC-I Upregulation)", width = 6, status = "info",
            plotlyOutput("plot_ifng", height = "320px")),
          box(title = "IL-6 & TNF-α", width = 6, status = "primary",
            plotlyOutput("plot_cyto", height = "320px"))
        ),
        fluidRow(
          box(title = "Muscle Pathology Explanation", width = 12, status = "success",
            fluidRow(
              column(4, h4("Endomysial Inflammation"),
                p("CD8+ T cells and macrophages infiltrate the endomysium, surrounding individual muscle fibers. This is the hallmark of PM histology.")),
              column(4, h4("MHC-I Upregulation"),
                p("IFN-γ drives MHC-I expression on normally MHC-I-negative myofibers. This enables CD8+ T cells to recognize and attack muscle fibers.")),
              column(4, h4("Necrosis & Regeneration"),
                p("Perforin/granzyme B from CD8+ cells causes necrosis. Satellite cells attempt regeneration (basophilic fibers on biopsy). CK leaks into serum."))
            )
          )
        )
      ),

      ## ── TAB 5: Clinical Endpoints ───────────────────────────────────────
      tabItem(tabName = "tab_clinical",
        fluidRow(
          box(title = "Serum CK Over Time", width = 6, status = "danger",
            plotlyOutput("plot_ck", height = "350px")),
          box(title = "MMT-8 Score Over Time", width = 6, status = "success",
            plotlyOutput("plot_mmt", height = "350px"))
        ),
        fluidRow(
          box(title = "Total Improvement Score (TIS)", width = 6, status = "primary",
            plotlyOutput("plot_tis", height = "320px"),
            p("TIS ≥20: Minimal; ≥40: Moderate; ≥60: Major improvement")),
          box(title = "Remission Status Over Time", width = 6, status = "success",
            plotlyOutput("plot_remission", height = "320px"),
            p("Remission: CK <1.5× ULN AND MMT-8 >72"))
        ),
        fluidRow(
          box(title = "Clinical Endpoints Timeline", width = 12, status = "info",
            DTOutput("endpoint_table"))
        )
      ),

      ## ── TAB 6: Scenario Comparison ──────────────────────────────────────
      tabItem(tabName = "tab_compare",
        fluidRow(
          box(title = "Run All 6 Scenarios", width = 3, status = "primary",
            p("Comparing all standard PM treatment strategies."),
            actionButton("run_all", "Compare All Scenarios",
                         icon = icon("chart-bar"), class = "btn-primary btn-lg",
                         width = "100%"),
            hr(),
            selectInput("compare_endpoint", "Endpoint to Display:",
              choices = c("Serum CK" = "CK_S",
                          "MMT-8 Score" = "MMT8",
                          "IFN-γ" = "IFNG",
                          "Muscle Inflammation" = "MUS_INF",
                          "Effector CD8+ T cells" = "CD8E",
                          "B Cells" = "BCELL",
                          "Autoantibodies" = "AUTOAB"))
          ),
          box(title = "Scenario Comparison Plot", width = 9, status = "primary",
            plotlyOutput("plot_compare", height = "500px"))
        ),
        fluidRow(
          box(title = "Day-365 Outcome Summary", width = 12, status = "success",
            DTOutput("compare_table"))
        )
      ),

      ## ── TAB 7: Biomarker Panel ──────────────────────────────────────────
      tabItem(tabName = "tab_biomark",
        fluidRow(
          box(title = "Biomarker Dashboard", width = 12, status = "info",
            fluidRow(
              column(6,
                h4("Standard Biomarkers"),
                plotlyOutput("plot_bm_standard", height = "400px")
              ),
              column(6,
                h4("IFN Signature & Novel Biomarkers"),
                plotlyOutput("plot_bm_novel", height = "400px")
              )
            )
          )
        ),
        fluidRow(
          box(title = "Biomarker Reference Ranges", width = 6, status = "warning",
            tableOutput("biomarker_ref")),
          box(title = "Drug Effect on Each Biomarker (Day 90)", width = 6,
              status = "success",
            plotlyOutput("plot_bm_drug_effect", height = "300px"))
        )
      )
    )
  )
)

## ════════════════════════════════════════════════════════════════════════════
## SERVER
## ════════════════════════════════════════════════════════════════════════════
server <- function(input, output, session) {

  ## Reactive simulation
  sim_data <- eventReactive(input$run_sim, {
    withProgress(message = "Running PM QSP simulation...", value = 0.5, {
      run_scenario(pm_mod,
        scenario    = input$scenario,
        pred_dose   = input$pred_dose,
        mtx_dose    = input$mtx_dose,
        aza_dose    = input$aza_dose,
        rtx_dose    = input$rtx_dose,
        jaki_dose   = input$jaki_dose,
        duration_d  = input$duration,
        init_ck     = input$init_ck,
        init_mmt    = input$init_mmt,
        init_ifng   = input$init_ifng
      )
    })
  }, ignoreNULL = FALSE)

  ## Run all 6 scenarios
  all_scenarios <- eventReactive(input$run_all, {
    withProgress(message = "Simulating all scenarios...", value = 0, {
      scen_list <- list(
        list(name="Untreated",         scenario=character(0), pred=0,  mtx=0,  aza=0,   rtx=0,    jaki=0),
        list(name="Prednisone Mono",   scenario="Prednisone", pred=60, mtx=0,  aza=0,   rtx=0,    jaki=0),
        list(name="Pred + MTX",        scenario=c("Prednisone","MTX"), pred=60,mtx=15,aza=0,rtx=0,jaki=0),
        list(name="Pred + AZA",        scenario=c("Prednisone","AZA"), pred=60,mtx=0,aza=150,rtx=0,jaki=0),
        list(name="Rituximab + Pred",  scenario=c("Prednisone","Rituximab"),pred=40,mtx=0,aza=0,rtx=1000,jaki=0),
        list(name="JAKi + Pred",       scenario=c("Prednisone","JAKi"), pred=20,mtx=0,aza=0,rtx=0,jaki=4)
      )
      results <- lapply(seq_along(scen_list), function(i) {
        setProgress(i / length(scen_list), message = paste("Simulating:", scen_list[[i]]$name))
        out <- run_scenario(pm_mod,
          scenario   = scen_list[[i]]$scenario,
          pred_dose  = scen_list[[i]]$pred,
          mtx_dose   = scen_list[[i]]$mtx,
          aza_dose   = scen_list[[i]]$aza,
          rtx_dose   = scen_list[[i]]$rtx,
          jaki_dose  = scen_list[[i]]$jaki,
          duration_d = 365,
          init_ck    = input$init_ck,
          init_mmt   = input$init_mmt,
          init_ifng  = input$init_ifng
        )
        out$Scenario <- scen_list[[i]]$name
        out
      })
      bind_rows(results)
    })
  })

  ## Helper: ggplotly wrapper
  gg2ly <- function(p, ...) ggplotly(p, tooltip = c("x", "y", "colour"), ...) %>%
    layout(legend = list(orientation = "h", y = -0.2))

  qsp_theme <- theme_bw(base_size = 11) +
    theme(legend.position = "bottom", panel.grid.minor = element_blank())

  ## ── Value Boxes ────────────────────────────────────────────────────────
  output$vbox_ck <- renderValueBox({
    d <- tail(sim_data(), 1)
    ck <- round(d$CK_S)
    status <- if (ck < 200) "green" else if (ck < 1000) "yellow" else "red"
    valueBox(paste(ck, "U/L"), "Serum CK at End", icon = icon("vial"), color = status)
  })
  output$vbox_mmt <- renderValueBox({
    d <- tail(sim_data(), 1)
    mmt <- round(d$MMT8, 1)
    status <- if (mmt > 72) "green" else if (mmt > 50) "yellow" else "red"
    valueBox(paste(mmt, "/80"), "MMT-8 Score", icon = icon("dumbbell"), color = status)
  })
  output$vbox_ifng <- renderValueBox({
    d <- tail(sim_data(), 1)
    ifng <- round(d$IFNG, 1)
    status <- if (ifng < 2) "green" else if (ifng < 8) "yellow" else "red"
    valueBox(paste(ifng, "pg/mL"), "IFN-γ", icon = icon("atom"), color = status)
  })
  output$vbox_remission <- renderValueBox({
    d <- tail(sim_data(), 1)
    rem <- d$Remission
    status <- if (rem == 1) "green" else "red"
    label <- if (rem == 1) "IN REMISSION" else "ACTIVE DISEASE"
    valueBox(label, "Disease Status", icon = icon("heartbeat"), color = status)
  })

  ## ── Static tables ─────────────────────────────────────────────────────
  output$autoab_table <- renderTable({
    data.frame(
      Antibody = c("Anti-Jo-1", "Anti-SRP", "Anti-MDA5", "Anti-Mi-2"),
      Frequency= c("~20% PM", "~5% PM", "~10% CADM", "~10% DM"),
      `ILD Risk`= c("High", "Moderate", "Very High", "Low"),
      `Prognosis`= c("Moderate", "Severe", "Variable", "Good"),
      check.names = FALSE
    )
  })
  output$das_table <- renderTable({
    data.frame(
      Category = c("Normal", "Mild", "Moderate", "Severe", "Remission"),
      `CK (U/L)`= c("<200", "200–1000", "1000–5000", ">5000", "<300"),
      `MMT-8`  = c("76–80", "65–75", "45–64", "<45", ">72"),
      check.names = FALSE
    )
  })
  output$biomarker_ref <- renderTable({
    data.frame(
      Biomarker = c("CK", "Aldolase", "LDH", "IFN-γ", "IL-6", "TNF-α", "Anti-Jo-1"),
      `Normal Range` = c("<200 U/L", "<8 U/L", "<250 U/L", "<5 pg/mL", "<5 pg/mL", "<8 pg/mL", "Negative"),
      `PM Active` = c(">1000 U/L", ">15 U/L", "↑", ">15 pg/mL", ">10 pg/mL", ">10 pg/mL", "Positive"),
      check.names = FALSE
    )
  })

  ## ── PK Plots ──────────────────────────────────────────────────────────
  output$plot_pk_pred <- renderPlotly({
    d <- sim_data()
    p <- ggplot(d, aes(time_d, Cpred)) + geom_line(color = "#1E88E5", linewidth = 1.1) +
      labs(x = "Time (days)", y = "Prednisolone (ng/mL)") + qsp_theme
    gg2ly(p)
  })
  output$plot_pk_mtx <- renderPlotly({
    d <- sim_data()
    p <- ggplot(d, aes(time_d, MTX_PG)) + geom_line(color = "#FB8C00", linewidth = 1.1) +
      labs(x = "Time (days)", y = "MTX Polyglutamates (μg/L equiv.)") + qsp_theme
    gg2ly(p)
  })
  output$plot_pk_jaki <- renderPlotly({
    d <- sim_data()
    p <- ggplot(d, aes(time_d, Cjaki)) + geom_line(color = "#8E24AA", linewidth = 1.1) +
      labs(x = "Time (days)", y = "Baricitinib Cp (nM)") + qsp_theme
    gg2ly(p)
  })
  output$plot_pk_rtx <- renderPlotly({
    d <- sim_data()
    p <- ggplot(d, aes(time_d, RTX_C1)) + geom_line(color = "#E53935", linewidth = 1.1) +
      labs(x = "Time (days)", y = "Rituximab Central Cp (mg)") + qsp_theme
    gg2ly(p)
  })
  output$pk_table <- renderDT({
    d <- sim_data()
    tpts <- c(0.5, 1, 2, 3, 7, 14, 30, 90, 180, 365)
    d_sub <- d[round(d$time_d, 1) %in% tpts, c("time_d", "Cpred", "MTX_PG", "Cjaki")]
    d_sub <- d_sub[!duplicated(round(d_sub$time_d, 0)), ]
    colnames(d_sub) <- c("Day", "Prednisolone (ng/mL)", "MTX-PG (μg/L)", "Baricitinib (nM)")
    datatable(d_sub, options = list(pageLength = 10), rownames = FALSE)
  })

  ## ── Immunology Plots ──────────────────────────────────────────────────
  make_im_plot <- function(col, ylab, col_hex) {
    d <- sim_data()
    ggplot(d, aes_string("time_d", col)) + geom_line(color = col_hex, linewidth = 1.1) +
      labs(x = "Time (days)", y = ylab) + qsp_theme
  }
  output$plot_cd8   <- renderPlotly(gg2ly(make_im_plot("CD8E",    "Effector CD8+ (cells/μL)", "#E53935")))
  output$plot_th1   <- renderPlotly(gg2ly(make_im_plot("CD4TH1",  "Th1 Cells (cells/μL)",     "#FB8C00")))
  output$plot_bcell <- renderPlotly(gg2ly(make_im_plot("BCELL",   "B Cells (cells/μL)",       "#1E88E5")))
  output$plot_autoab<- renderPlotly(gg2ly(make_im_plot("AUTOAB",  "Autoantibodies (U/mL)",    "#43A047")))

  ## ── Muscle Pathology Plots ─────────────────────────────────────────────
  output$plot_inf  <- renderPlotly(gg2ly(make_im_plot("MUS_INF", "Inflammation Index (0–10)", "#E53935")))
  output$plot_mhci <- renderPlotly(gg2ly(make_im_plot("MHCI",    "MHC-I Fold-Change",         "#8E24AA")))
  output$plot_ifng <- renderPlotly(gg2ly(make_im_plot("IFNG",    "IFN-γ (pg/mL)",             "#00ACC1")))
  output$plot_cyto <- renderPlotly({
    d <- sim_data() %>% tidyr::pivot_longer(cols = c("TNFa", "IL6C"),
                                             names_to = "Cytokine", values_to = "Conc")
    p <- ggplot(d, aes(time_d, Conc, color = Cytokine)) + geom_line(linewidth = 1.1) +
      scale_color_manual(values = c("TNFa" = "#E53935", "IL6C" = "#1E88E5"),
                         labels = c("TNFa" = "TNF-α", "IL6C" = "IL-6")) +
      labs(x = "Time (days)", y = "Cytokine (pg/mL)", color = NULL) + qsp_theme
    gg2ly(p)
  })

  ## ── Clinical Endpoint Plots ────────────────────────────────────────────
  output$plot_ck <- renderPlotly({
    d <- sim_data()
    p <- ggplot(d, aes(time_d, CK_S)) +
      geom_line(color = "#E53935", linewidth = 1.2) +
      geom_hline(yintercept = 200, linetype = "dashed", color = "gray") +
      scale_y_log10(labels = scales::comma) +
      labs(x = "Time (days)", y = "CK (U/L)", title = "Serum CK") + qsp_theme
    gg2ly(p)
  })
  output$plot_mmt <- renderPlotly({
    d <- sim_data()
    p <- ggplot(d, aes(time_d, MMT8)) +
      geom_line(color = "#43A047", linewidth = 1.2) +
      geom_hline(yintercept = 72, linetype = "dashed", color = "gray") +
      ylim(0, 80) +
      labs(x = "Time (days)", y = "MMT-8 Score") + qsp_theme
    gg2ly(p)
  })
  output$plot_tis <- renderPlotly({
    d <- sim_data()
    p <- ggplot(d, aes(time_d, TIS_approx)) +
      geom_line(color = "#1E88E5", linewidth = 1.2) +
      geom_hline(yintercept = 20, linetype="dashed", color="orange") +
      geom_hline(yintercept = 40, linetype="dashed", color="green") +
      labs(x = "Time (days)", y = "TIS (approx, %)") + qsp_theme
    gg2ly(p)
  })
  output$plot_remission <- renderPlotly({
    d <- sim_data()
    p <- ggplot(d, aes(time_d, Remission)) +
      geom_step(color = "#43A047", linewidth = 1.2) +
      scale_y_continuous(breaks = c(0, 1), labels = c("Active", "Remission")) +
      labs(x = "Time (days)", y = "Disease Status") + qsp_theme
    gg2ly(p)
  })
  output$endpoint_table <- renderDT({
    d <- sim_data()
    tpts_d <- c(0, 7, 14, 30, 60, 90, 180, 365)
    d2 <- d[!duplicated(round(d$time_d, 0)), ]
    d3 <- d2[round(d2$time_d, 0) %in% tpts_d,
              c("time_d","CK_S","MMT8","IFNG","MUS_INF","Remission")]
    d3 <- d3[!duplicated(d3$time_d), ]
    d3$CK_ULN <- round(d3$CK_S / 200, 1)
    d3$MMT8   <- round(d3$MMT8, 1)
    d3$IFNG   <- round(d3$IFNG, 1)
    d3$Remission <- ifelse(d3$Remission == 1, "YES", "NO")
    colnames(d3) <- c("Day","CK (U/L)","MMT-8","IFN-γ (pg/mL)",
                       "Muscle Inf.","Remission","CK × ULN")
    datatable(d3, options = list(pageLength = 10), rownames = FALSE)
  })

  ## ── Scenario Comparison Plots ──────────────────────────────────────────
  output$plot_compare <- renderPlotly({
    req(all_scenarios())
    d <- all_scenarios()
    col_var <- input$compare_endpoint
    scen_cols <- c("#E53935","#1E88E5","#43A047","#FB8C00","#8E24AA","#00ACC1")
    scen_levs <- c("Untreated","Prednisone Mono","Pred + MTX",
                   "Pred + AZA","Rituximab + Pred","JAKi + Pred")
    d$Scenario <- factor(d$Scenario, levels = scen_levs)
    p <- ggplot(d, aes_string("time_d", col_var, color = "Scenario")) +
      geom_line(linewidth = 0.9) +
      scale_color_manual(values = setNames(scen_cols, scen_levs)) +
      labs(x = "Time (days)", y = col_var, color = NULL) + qsp_theme
    if (col_var == "CK_S") p <- p + scale_y_log10(labels = scales::comma)
    gg2ly(p)
  })

  output$compare_table <- renderDT({
    req(all_scenarios())
    d <- all_scenarios()
    d_end <- d[round(d$time_d, 0) == 365 | d$time_d == max(d$time_d), ]
    d_end <- d_end[!duplicated(d_end$Scenario), ]
    out <- data.frame(
      Scenario    = d_end$Scenario,
      `CK (U/L)`  = round(d_end$CK_S),
      `MMT-8`     = round(d_end$MMT8, 1),
      `IFN-γ`     = round(d_end$IFNG, 1),
      `MUS Inf`   = round(d_end$MUS_INF, 2),
      Remission   = ifelse(d_end$Remission == 1, "YES", "NO"),
      check.names = FALSE
    )
    datatable(out, options = list(pageLength = 10), rownames = FALSE)
  })

  ## ── Biomarker Panel ────────────────────────────────────────────────────
  output$plot_bm_standard <- renderPlotly({
    d <- sim_data() %>%
      tidyr::pivot_longer(cols = c("CK_S","AUTOAB"),
                           names_to = "Biomarker", values_to = "Value")
    p <- ggplot(d, aes(time_d, Value, color = Biomarker)) +
      geom_line(linewidth = 1.1) +
      facet_wrap(~Biomarker, scales = "free_y") +
      scale_color_manual(values = c("CK_S" = "#E53935", "AUTOAB" = "#1E88E5")) +
      labs(x = "Time (days)", y = "Value", color = NULL) + qsp_theme
    gg2ly(p)
  })
  output$plot_bm_novel <- renderPlotly({
    d <- sim_data() %>%
      tidyr::pivot_longer(cols = c("IFNG","IL6C","TNFa","MHCI"),
                           names_to = "Marker", values_to = "Value")
    p <- ggplot(d, aes(time_d, Value, color = Marker)) +
      geom_line(linewidth = 1.0) +
      facet_wrap(~Marker, scales = "free_y") +
      scale_color_manual(values = c("IFNG"="#E53935","IL6C"="#1E88E5",
                                     "TNFa"="#FB8C00","MHCI"="#8E24AA")) +
      labs(x = "Time (days)", y = "Value", color = NULL) + qsp_theme
    gg2ly(p)
  })
  output$plot_bm_drug_effect <- renderPlotly({
    d <- sim_data()
    d90 <- d[round(d$time_d, 0) == 90 | d$time_d == max(d$time_d[d$time_d <= 90]), ]
    if (nrow(d90) == 0) return(NULL)
    d90 <- tail(d90, 1)
    bars <- data.frame(
      Biomarker = c("CK (×ULN)", "MMT8 (%max)", "IFN-γ (×norm)",
                    "CD8E (×norm)", "B cells (×norm)"),
      Value     = c(d90$CK_S / 200, d90$MMT8 / 80,
                    d90$IFNG / 0.5, d90$CD8E / 5,
                    d90$BCELL / 90)
    )
    p <- ggplot(bars, aes(Biomarker, Value, fill = Biomarker)) +
      geom_bar(stat = "identity") +
      geom_hline(yintercept = 1, linetype = "dashed") +
      scale_fill_brewer(palette = "Set2") +
      labs(x = NULL, y = "Relative Value (normalized to normal)") +
      theme_bw() + theme(legend.position = "none", axis.text.x = element_text(angle = 25, hjust = 1))
    gg2ly(p)
  })
}

## ── Launch ──────────────────────────────────────────────────────────────────
shinyApp(ui = ui, server = server)
