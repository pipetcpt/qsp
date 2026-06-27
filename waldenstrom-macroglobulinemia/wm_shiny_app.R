## ============================================================
##  Waldenström's Macroglobulinemia (WM) — Interactive Shiny Dashboard
##  Filename : wm_shiny_app.R
##  Run with : shiny::runApp("waldenstrom-macroglobulinemia/wm_shiny_app.R")
## ============================================================

library(shiny)
library(shinydashboard)
library(mrgsolve)
library(dplyr)
library(ggplot2)
library(tidyr)
library(DT)
library(plotly)
library(patchwork)

# ── Embedded mrgsolve model ─────────────────────────────────
code <- '
$PARAM
  IBR_dose=0, ka_ibr=1.2, F_ibr=0.10, CL_ibr=73, V_ibr=820,
  ZAN_dose=0, ka_zan=1.5, F_zan=0.60, CL_zan=29, V_zan=520,
  VEN_dose=0, ka_ven=0.40, F_ven=0.72, CL_ven=12, V_ven=256,
  RTX_dose=0, CL_RTX=0.21, V_RTX=4.5,
  BOR_Cp=0, EC50_BOR=10,
  EC50_ibr_btk=0.5, EC50_zan_btk=0.30, kout_BTK=0.008,
  MW_ibr=440, MW_zan=471, VEN_MW=868,
  NFkB_base=1.0, MYD88_drive=0.60, kBTK_NFkB=0.35,
  EC50_BOR_NFkB=50, kout_NFkB=0.50,
  kel0_LPC=0.002, kprolif_LPC=0.0045, NFkB_kLPC=0.30,
  LPC0=50, kconv_LPC=0.00015, kel0_PC=0.001, PC0=10, KMAX_BM=200,
  ksec_IgM=0.06, kel_IgM=0.004, IgM0=25.0,
  Hgb0=9.5, kprod_Hgb=0.030, kel_Hgb=0.0028, BMInf_Hgb=0.60,
  Visc_base=1.5, kIgM_visc=0.080, Visc_exp=1.4,
  BCL2_base=1.0, kNFkB_BCL2=0.40, kout_BCL2=0.10, EC50_VEN=0.5,
  EC50_RTX_CD20=0.05, kout_CD20=0.004, kADCC=0.025,
  kout_Prot=0.050, kprot_kill=0.012,
  NK0=100, kprod_NK=0.008, kel_NK=0.006,
  BENDA_flag=0, BENDA_kLPC=0.010,
  // Patient characteristics
  MYD88_mut=1,    // 1=L265P, 0=WT (very rare)
  CXCR4_mut=0,    // 1=WHIM mutation (reduced ibrutinib efficacy)
  cxcr4_factor=1.0 // reduces BTK response if CXCR4 WHIM

$CMT IBR_gut IBR_C ZAN_gut ZAN_C RTX_C VEN_gut VEN_C
     BTK_occ NFkB LPC PC IgM Hgb Visc BMInf BCL2 CD20 Apop Protsm NK

$MAIN
  _F(IBR_gut)=F_ibr; _F(ZAN_gut)=F_zan; _F(VEN_gut)=F_ven;
  LPC_0=LPC0; PC_0=PC0; IgM_0=IgM0; Hgb_0=Hgb0;
  Visc_0=Visc_base+kIgM_visc*pow(IgM0,Visc_exp);
  BMInf_0=(LPC0+PC0)/KMAX_BM;
  BCL2_0=BCL2_base; CD20_0=1.0; Protsm_0=100.0; NK_0=NK0;
  BTK_occ_0=0.0; NFkB_0=NFkB_base+MYD88_drive*MYD88_mut;

$ODE
  double ke_ibr=CL_ibr/V_ibr;
  dxdt_IBR_gut=-ka_ibr*IBR_gut;
  dxdt_IBR_C=ka_ibr*IBR_gut-ke_ibr*IBR_C;
  double Cp_ibr_nM=(IBR_C/V_ibr*1000.0)/MW_ibr*1000.0;

  double ke_zan=CL_zan/V_zan;
  dxdt_ZAN_gut=-ka_zan*ZAN_gut;
  dxdt_ZAN_C=ka_zan*ZAN_gut-ke_zan*ZAN_C;
  double Cp_zan_nM=(ZAN_C/V_zan*1000.0)/MW_zan*1000.0;

  dxdt_RTX_C=-(CL_RTX/V_RTX)*RTX_C;
  double Cp_RTX=RTX_C/V_RTX;

  double ke_ven=CL_ven/V_ven;
  dxdt_VEN_gut=-ka_ven*VEN_gut;
  dxdt_VEN_C=ka_ven*VEN_gut-ke_ven*VEN_C;
  double Cp_ven_uM=(VEN_C/V_ven*1000.0)/VEN_MW*1000.0;

  double IBR_eff=Cp_ibr_nM/(Cp_ibr_nM+EC50_ibr_btk)*cxcr4_factor;
  double ZAN_eff=Cp_zan_nM/(Cp_zan_nM+EC50_zan_btk)*cxcr4_factor;
  double BTK_input=fmax(IBR_eff,ZAN_eff);
  dxdt_BTK_occ=BTK_input*(1.0-BTK_occ)-kout_BTK*BTK_occ;

  double BOR_inhib_NF=EC50_BOR_NFkB/(EC50_BOR_NFkB+(100.0-Protsm));
  double NFkB_in=NFkB_base+MYD88_drive*MYD88_mut+kBTK_NFkB*(1.0-BTK_occ);
  dxdt_NFkB=NFkB_in*BOR_inhib_NF-kout_NFkB*NFkB;

  double VEN_inh=Cp_ven_uM/(Cp_ven_uM+EC50_VEN);
  dxdt_BCL2=kNFkB_BCL2*NFkB*(1.0-VEN_inh)-kout_BCL2*BCL2;

  double apo=BTK_occ*0.4+VEN_inh*0.5+(100.0-Protsm)/100.0*kprot_kill/kel0_LPC;
  dxdt_Apop=apo-0.5*Apop;

  double RTX_inh=Cp_RTX/(Cp_RTX+EC50_RTX_CD20);
  dxdt_CD20=kout_CD20*(1.0-CD20)-RTX_inh*CD20;
  dxdt_NK=kprod_NK*NK0-kel_NK*NK;

  double BOR_inhib=BOR_Cp/(BOR_Cp+EC50_BOR);
  dxdt_Protsm=kout_Prot*(100.0-Protsm)-BOR_inhib*Protsm;

  double total_tumor=LPC+PC;
  BMInf=fmin(total_tumor/KMAX_BM,1.0);
  dxdt_BMInf=0;

  double LPC_prolif=kprolif_LPC*(1.0+NFkB_kLPC*NFkB)*LPC*(1.0-total_tumor/KMAX_BM);
  double LPC_apop=(kel0_LPC+apo)*LPC;
  double LPC_ADCC=kADCC*NK*CD20*RTX_inh*LPC/(LPC+1.0);
  double LPC_BENDA=BENDA_flag*BENDA_kLPC*LPC;
  double LPC_conv=kconv_LPC*LPC;
  dxdt_LPC=LPC_prolif-LPC_apop-LPC_ADCC-LPC_BENDA-LPC_conv;

  double PC_apop=(kel0_PC+apo*0.6)*PC;
  double PC_ADCC=kADCC*0.5*NK*CD20*RTX_inh*PC/(PC+1.0);
  dxdt_PC=LPC_conv-PC_apop-PC_ADCC-BENDA_flag*BENDA_kLPC*0.7*PC;

  dxdt_IgM=ksec_IgM*PC-kel_IgM*IgM;
  double suppression=1.0-BMInf_Hgb*BMInf;
  dxdt_Hgb=kprod_Hgb*Hgb0*suppression-kel_Hgb*Hgb;
  dxdt_Visc=(Visc_base+kIgM_visc*pow(fmax(IgM,0),Visc_exp)-Visc)*0.5;

$TABLE
  double Cp_IBR_ng=IBR_C/V_ibr*1000.0;
  double Cp_ZAN_ng=ZAN_C/V_zan*1000.0;
  double Cp_RTX_mgl=Cp_RTX;
  double Cp_VEN_uM2=Cp_ven_uM;
  double BTK_pct=BTK_occ*100.0;
  double IgM_pct_change=(IgM-IgM0)/IgM0*100.0;

$CAPTURE Cp_IBR_ng Cp_ZAN_ng Cp_RTX_mgl Cp_VEN_uM2
         BTK_pct NFkB BCL2 LPC PC IgM Hgb Visc
         BMInf NK Protsm CD20 Apop IgM_pct_change
'
mod <- tryCatch(mcode("WM_shiny", code), error = function(e) NULL)

# ── Colour palette ───────────────────────────────────────────
pal <- c(
  "Watch & Wait"          = "#7F8C8D",
  "Ibrutinib"             = "#2980B9",
  "Ibrutinib + Rituximab" = "#1A5276",
  "Zanubrutinib"          = "#8E44AD",
  "R-Bendamustine"        = "#27AE60",
  "BDR"                   = "#D35400",
  "Venetoclax"            = "#C0392B"
)

# ── UI ───────────────────────────────────────────────────────
ui <- dashboardPage(
  skin = "purple",
  dashboardHeader(title = "WM QSP Dashboard"),

  dashboardSidebar(
    sidebarMenu(
      menuItem("Patient Profile",     tabName = "patient",   icon = icon("user")),
      menuItem("Drug PK",             tabName = "pk",        icon = icon("pills")),
      menuItem("PD Key Endpoints",    tabName = "pd",        icon = icon("chart-line")),
      menuItem("Clinical Endpoints",  tabName = "clinical",  icon = icon("hospital")),
      menuItem("Scenario Comparison", tabName = "scenario",  icon = icon("balance-scale")),
      menuItem("Biomarker Panel",     tabName = "biomarker", icon = icon("dna")),
      menuItem("About / References",  tabName = "about",     icon = icon("book"))
    )
  ),

  dashboardBody(
    tags$head(tags$style(HTML("
      .content-wrapper { background-color: #F4F6F7; }
      .box-header { background-color: #6C3483; color: white; }
    "))),

    tabItems(

      ## TAB 1 — Patient Profile ───────────────────────────────
      tabItem("patient",
        fluidRow(
          box(title = "Patient Characteristics", status = "primary", solidHeader = TRUE, width = 4,
            numericInput("IgM0",  "Baseline IgM (g/L):",   25,  min = 1,  max = 80),
            numericInput("LPC0",  "Baseline LPC (×10⁹):", 50,  min = 5,  max = 200),
            numericInput("Hgb0",  "Baseline Hgb (g/dL):",  9.5, min = 5,  max = 15),
            numericInput("PC0",   "Baseline PC (×10⁹):",  10,  min = 1,  max = 100),
            selectInput("MYD88",  "MYD88 status:",
                        choices = c("L265P (mutated ~95%)" = 1, "Wild-type (~5%)" = 0)),
            selectInput("CXCR4",  "CXCR4 WHIM mutation:",
                        choices = c("Absent" = 0, "Present (~35%)" = 1)),
            numericInput("sim_days", "Simulation (days):", 730, min = 90, max = 1825)
          ),
          box(title = "Baseline Disease Severity", status = "warning", solidHeader = TRUE, width = 4,
            h4("International Prognostic Scoring for WM (IPSSWM)"),
            selectInput("age_group", "Age:", choices = c("≤65" = 0, ">65" = 1)),
            selectInput("hgb_risk", "Hgb ≤11.5 g/dL:", choices = c("No" = 0, "Yes" = 1)),
            selectInput("plt_risk", "Platelets ≤100 ×10⁹/L:", choices = c("No" = 0, "Yes" = 1)),
            selectInput("b2m_risk", "β₂M >3 mg/L:", choices = c("No" = 0, "Yes" = 1)),
            selectInput("igm_risk", "IgM >7 g/dL:", choices = c("No" = 0, "Yes" = 1)),
            htmlOutput("ipss_result")
          ),
          box(title = "WM Disease Overview", status = "info", solidHeader = TRUE, width = 4,
            h4("Waldenström's Macroglobulinemia"),
            p("A lymphoplasmacytic lymphoma characterised by:"),
            tags$ul(
              tags$li("MYD88 L265P mutation (~95%)"),
              tags$li("CXCR4 WHIM mutation (~35%)"),
              tags$li("BM infiltration ≥10% LPC/PC"),
              tags$li("IgM monoclonal paraprotein"),
              tags$li("Key symptoms: anemia, hyperviscosity, neuropathy")
            ),
            p("Diagnostic criteria (Owen 2003): BM infiltration + IgM"),
            hr(),
            h5("Key Drug Targets:"),
            tags$ul(
              tags$li("BTK (ibrutinib, zanubrutinib)"),
              tags$li("CD20 (rituximab, ofatumumab)"),
              tags$li("BCL-2 (venetoclax)"),
              tags$li("26S Proteasome (bortezomib)")
            )
          )
        ),
        fluidRow(
          box(title = "Viscosity Risk Chart", status = "danger", solidHeader = TRUE, width = 6,
            plotOutput("visc_risk_chart", height = "280px")
          ),
          box(title = "IgM → Response Criteria (Owen/IWWM)", status = "success", solidHeader = TRUE, width = 6,
            tableOutput("response_table")
          )
        )
      ),

      ## TAB 2 — Drug PK ────────────────────────────────────────
      tabItem("pk",
        fluidRow(
          box(title = "Treatment Selection", status = "primary", solidHeader = TRUE, width = 3,
            checkboxGroupInput("pk_drugs", "Drugs to simulate:",
              choices = c("Ibrutinib (420 mg/d)" = "IBR",
                          "Zanubrutinib (160 mg BID)" = "ZAN",
                          "Rituximab (375 mg/m²)" = "RTX",
                          "Venetoclax (ramp 20→400)" = "VEN"),
              selected = c("IBR")
            ),
            numericInput("pk_days", "PK simulation (days):", 14, min = 1, max = 30),
            hr(),
            h5("Reference Cmax values:"),
            p("Ibrutinib: ~118 ng/mL"),
            p("Zanubrutinib: ~270 ng/mL"),
            p("Rituximab: ~150 mg/L (after infusion)"),
            p("Venetoclax: ~0.5-2 µM (at 400 mg)")
          ),
          box(title = "Plasma Concentration–Time Profiles", status = "info", solidHeader = TRUE, width = 9,
            plotlyOutput("pk_plot", height = "500px")
          )
        ),
        fluidRow(
          box(title = "BTK Occupancy Kinetics", status = "warning", solidHeader = TRUE, width = 6,
            plotlyOutput("btk_occ_plot", height = "320px")
          ),
          box(title = "PK Parameters Summary", status = "info", solidHeader = TRUE, width = 6,
            tableOutput("pk_table")
          )
        )
      ),

      ## TAB 3 — PD Key Endpoints ───────────────────────────────
      tabItem("pd",
        fluidRow(
          box(title = "PD Controls", status = "primary", solidHeader = TRUE, width = 3,
            selectInput("pd_scenario", "Treatment scenario:",
              choices = c("Watch & Wait", "Ibrutinib",
                          "Ibrutinib + Rituximab", "Zanubrutinib",
                          "R-Bendamustine", "BDR", "Venetoclax"),
              selected = "Ibrutinib"),
            hr(),
            h5("Key PD Pathway Nodes:"),
            p("NFkB: NF-κB activity (AU)"),
            p("BCL2: BCL-2 anti-apoptotic (AU)"),
            p("LPC: Lymphoplasmacytic cells"),
            p("PC: IgM-secreting plasma cells"),
            p("IgM: Serum IgM (g/L)")
          ),
          box(title = "NF-κB Activity & Downstream", status = "danger", solidHeader = TRUE, width = 9,
            plotlyOutput("nfkb_plot", height = "480px")
          )
        ),
        fluidRow(
          box(title = "LPC & PC Dynamics", status = "warning", solidHeader = TRUE, width = 6,
            plotlyOutput("lpc_pc_plot", height = "340px")
          ),
          box(title = "IgM Production vs Elimination", status = "info", solidHeader = TRUE, width = 6,
            plotlyOutput("igm_prod_plot", height = "340px")
          )
        )
      ),

      ## TAB 4 — Clinical Endpoints ─────────────────────────────
      tabItem("clinical",
        fluidRow(
          box(title = "Clinical Controls", status = "primary", solidHeader = TRUE, width = 3,
            selectInput("clin_scenario", "Treatment scenario:",
              choices = c("Watch & Wait", "Ibrutinib",
                          "Ibrutinib + Rituximab", "Zanubrutinib",
                          "R-Bendamustine", "BDR", "Venetoclax"),
              selected = "Ibrutinib + Rituximab"),
            hr(),
            h5("Thresholds:"),
            p("HVS: viscosity >3.5 cP"),
            p("Anemia: Hgb <10 g/dL"),
            p("VGPR: IgM ↓>90%"),
            p("PR:   IgM ↓>50%")
          ),
          box(title = "IgM & Hemoglobin Response", status = "success", solidHeader = TRUE, width = 9,
            plotlyOutput("clin_main_plot", height = "480px")
          )
        ),
        fluidRow(
          box(title = "Serum Viscosity", status = "danger", solidHeader = TRUE, width = 6,
            plotlyOutput("visc_plot", height = "340px")
          ),
          box(title = "BM Infiltration (%)", status = "warning", solidHeader = TRUE, width = 6,
            plotlyOutput("bm_plot", height = "340px")
          )
        )
      ),

      ## TAB 5 — Scenario Comparison ────────────────────────────
      tabItem("scenario",
        fluidRow(
          box(title = "Select Scenarios", status = "primary", solidHeader = TRUE, width = 3,
            checkboxGroupInput("cmp_scenarios", "Scenarios:",
              choices = c("Watch & Wait", "Ibrutinib",
                          "Ibrutinib + Rituximab", "Zanubrutinib",
                          "R-Bendamustine", "BDR", "Venetoclax"),
              selected = c("Watch & Wait", "Ibrutinib",
                           "Ibrutinib + Rituximab", "Zanubrutinib")
            ),
            selectInput("cmp_endpoint", "Primary endpoint:",
              choices = c("IgM (g/L)" = "IgM",
                          "Hemoglobin (g/dL)" = "Hgb",
                          "Serum Viscosity (cP)" = "Visc",
                          "BM Infiltration (%)" = "BMInf",
                          "LPC cells (×10⁹)" = "LPC"),
              selected = "IgM"
            ),
            numericInput("cmp_days", "Horizon (days):", 730, min = 90, max = 1825)
          ),
          box(title = "Multi-Scenario Trajectory", status = "info", solidHeader = TRUE, width = 9,
            plotlyOutput("cmp_plot", height = "480px")
          )
        ),
        fluidRow(
          box(title = "Response at 12 months (360 d)", status = "success", solidHeader = TRUE, width = 12,
            DTOutput("response_dt")
          )
        )
      ),

      ## TAB 6 — Biomarker Panel ────────────────────────────────
      tabItem("biomarker",
        fluidRow(
          box(title = "Biomarker Controls", status = "primary", solidHeader = TRUE, width = 3,
            selectInput("bio_scenario", "Treatment:",
              choices = c("Watch & Wait", "Ibrutinib",
                          "Ibrutinib + Rituximab", "Zanubrutinib",
                          "R-Bendamustine", "BDR", "Venetoclax"),
              selected = "Ibrutinib + Rituximab"),
            checkboxGroupInput("bio_vars", "Biomarkers to show:",
              choices = c("BTK Occupancy" = "BTK_pct",
                          "NF-κB Activity" = "NFkB",
                          "BCL-2 Activity" = "BCL2",
                          "NK Cells" = "NK",
                          "Proteasome" = "Protsm",
                          "CD20 Surface" = "CD20"),
              selected = c("BTK_pct", "NFkB", "BCL2", "NK")
            )
          ),
          box(title = "Pharmacodynamic Biomarker Trajectories", status = "warning", solidHeader = TRUE, width = 9,
            plotlyOutput("bio_plot", height = "500px")
          )
        ),
        fluidRow(
          box(title = "Mechanistic Waterfall: IgM at 6 months", status = "danger", solidHeader = TRUE, width = 6,
            plotlyOutput("waterfall_plot", height = "340px")
          ),
          box(title = "Biomarker Summary Table (Day 90 vs 180 vs 360)", status = "info", solidHeader = TRUE, width = 6,
            DTOutput("bio_table")
          )
        )
      ),

      ## TAB 7 — About ──────────────────────────────────────────
      tabItem("about",
        fluidRow(
          box(title = "Model Overview", status = "primary", solidHeader = TRUE, width = 6,
            h4("Waldenström's Macroglobulinemia QSP Model"),
            p("This dashboard simulates the quantitative systems pharmacology of WM."),
            h5("Mechanistic Map: 12 clusters, 100+ nodes"),
            tags$ul(
              tags$li("Drug PK: Ibrutinib, Zanubrutinib, Rituximab, Venetoclax, Bortezomib"),
              tags$li("BTK/BCR Signaling pathway"),
              tags$li("MYD88 L265P / NF-κB constitutive activation"),
              tags$li("PI3K/AKT/mTOR axis"),
              tags$li("CXCR4/BM homing (WHIM mutation)"),
              tags$li("Lymphoplasmacytic cell biology"),
              tags$li("Tumor microenvironment (BMSC, NK, Treg)"),
              tags$li("Apoptosis pathways (BCL-2/BAX/caspase)"),
              tags$li("IgM-mediated disease (HVS, neuropathy, cryoglobulinemia)"),
              tags$li("Drug PD target engagement"),
              tags$li("Clinical endpoints (Owen/IWWM criteria)"),
              tags$li("Genomic landscape")
            ),
            h5("ODE compartments: 20"),
            h5("Treatment scenarios: 7")
          ),
          box(title = "Key References", status = "info", solidHeader = TRUE, width = 6,
            h5("Pivotal Clinical Trials:"),
            tags$ul(
              tags$li("Treon SP et al. NEJM 2015 (Ibrutinib)"),
              tags$li("Dimopoulos MA et al. NEJM 2018 (INNOVATE: iR)"),
              tags$li("Tam CS et al. JCO 2020 (ASPEN: Zanubrutinib)"),
              tags$li("Rummel MJ et al. Lancet 2013 (R-Benda)"),
              tags$li("Dimopoulos MA et al. JCO 2013 (BDR)"),
              tags$li("Castillo JJ et al. Blood 2018 (Venetoclax)")
            ),
            h5("Mechanistic:"),
            tags$ul(
              tags$li("Treon SP et al. NEJM 2012 (MYD88 L265P discovery)"),
              tags$li("Hunter ZR et al. Blood 2014 (CXCR4 WHIM)"),
              tags$li("Yang G et al. Blood 2013 (MYD88→BTK→NF-κB)")
            ),
            h5("Model: Claude Code Routine (CCR) — 2026-06-27")
          )
        )
      )
    )
  )
)

# ── Server ───────────────────────────────────────────────────
server <- function(input, output, session) {

  if (is.null(mod)) {
    showNotification("mrgsolve model failed to compile — demo mode", type = "error")
    return()
  }

  # ── Reactive: build patient params ──────────────────────────
  pt_params <- reactive({
    list(
      IgM0       = input$IgM0,
      LPC0       = input$LPC0,
      Hgb0       = input$Hgb0,
      PC0        = input$PC0,
      MYD88_mut  = as.numeric(input$MYD88),
      CXCR4_mut  = as.numeric(input$CXCR4),
      cxcr4_factor = ifelse(as.numeric(input$CXCR4) == 1, 0.65, 1.0)
    )
  })

  # ── Reactive: run one scenario ───────────────────────────────
  run_sim <- function(scenario_name, days) {
    end_h  <- days * 24
    ka_ibr_dose <- 420 * 0.10
    ka_zan_dose <- 160 * 0.60

    e <- switch(scenario_name,
      "Watch & Wait"          = ev(amt = 0, cmt = "IBR_gut", time = 0),
      "Ibrutinib"             = ev(cmt = "IBR_gut", amt = ka_ibr_dose, ii = 24, addl = days - 1),
      "Ibrutinib + Rituximab" = ev_c(
          ev(cmt = "IBR_gut", amt = ka_ibr_dose, ii = 24, addl = days - 1),
          ev(cmt = "RTX_C",   amt = 375 * 1.8,
             time = seq(0, by = 28 * 24, length.out = 6))),
      "Zanubrutinib"          = ev(cmt = "ZAN_gut", amt = ka_zan_dose, ii = 12, addl = days * 2 - 1),
      "R-Bendamustine"        = ev_c(
          ev(cmt = "RTX_C", amt = 375 * 1.8,
             time = seq(0, by = 28 * 24, length.out = 6))),
      "BDR"                   = ev(cmt = "RTX_C", amt = 375 * 1.8,
                                    time = seq(0, by = 21 * 24, length.out = 6)),
      "Venetoclax"            = ev(cmt = "VEN_gut",
                                    amt = c(rep(20 * 0.72, 7), rep(50 * 0.72, 7),
                                            rep(100 * 0.72, 7), rep(200 * 0.72, 7),
                                            rep(400 * 0.72, days - 28)),
                                    time = seq(0, by = 24, length.out = days))
    )

    extra_p <- if (scenario_name == "BDR") list(BOR_Cp = 200) else list()

    mrgsim(mod,
      events  = e,
      param   = c(pt_params(), extra_p),
      end     = end_h,
      delta   = 24,
      obsonly = TRUE
    ) %>% as_tibble() %>%
      mutate(scenario = scenario_name, day = time / 24)
  }

  # ── IPSSWM Score ─────────────────────────────────────────────
  output$ipss_result <- renderUI({
    score <- as.numeric(input$age_group) + as.numeric(input$hgb_risk) +
             as.numeric(input$plt_risk)  + as.numeric(input$b2m_risk) +
             as.numeric(input$igm_risk)
    risk <- if (score <= 1) "Low" else if (score <= 2) "Intermediate" else "High"
    color <- if (risk == "Low") "green" else if (risk == "Intermediate") "orange" else "red"
    HTML(paste0("<b>IPSSWM Score: ", score, "</b><br>",
                "<span style='color:", color, ";font-size:18px;'>▶ ", risk,
                " Risk</span>"))
  })

  # ── Viscosity risk chart ──────────────────────────────────────
  output$visc_risk_chart <- renderPlot({
    igm_seq <- seq(0, 60, by = 0.5)
    visc_v  <- 1.5 + 0.08 * igm_seq^1.4
    df <- data.frame(IgM = igm_seq, Viscosity = visc_v)
    ggplot(df, aes(IgM, Viscosity)) +
      geom_ribbon(aes(ymin = 0, ymax = Viscosity,
                      fill = cut(Viscosity, c(0, 2, 3.5, Inf))), alpha = 0.4) +
      geom_line(linewidth = 1.5) +
      scale_fill_manual(values = c("#82E0AA", "#F9E79F", "#F1948A"),
                        labels = c("Normal", "Borderline", "HVS risk"),
                        name = "Zone") +
      geom_hline(yintercept = 3.5, color = "red", linetype = "dashed") +
      labs(x = "Serum IgM (g/L)", y = "Viscosity (cP)",
           title = "IgM → Viscosity relationship") +
      theme_bw()
  })

  # ── Response criteria table ───────────────────────────────────
  output$response_table <- renderTable({
    data.frame(
      Category = c("CR", "VGPR", "PR", "MR", "SD", "PD"),
      IgM_criteria = c("Normal + BM normal", "↓ >90%", "↓ >50%",
                        "↓ 25-50%", "↓<25% or ↑<25%", "↑ >25%"),
      Clinical = c("Resolution of symptoms", "Near-complete", "Major",
                    "Minor", "Stable", "Progression")
    )
  }, striped = TRUE, hover = TRUE)

  # ── TAB 2: PK plot ───────────────────────────────────────────
  output$pk_plot <- renderPlotly({
    days <- input$pk_days
    plots_list <- list()

    if ("IBR" %in% input$pk_drugs) {
      e <- ev(cmt = "IBR_gut", amt = 420 * 0.10, ii = 24, addl = days - 1)
      d <- mrgsim(mod, events = e, end = days * 24, delta = 1, obsonly = TRUE) %>%
           as_tibble() %>% mutate(day = time / 24, Cp = Cp_IBR_ng, drug = "Ibrutinib (ng/mL)")
      plots_list[["IBR"]] <- d
    }
    if ("ZAN" %in% input$pk_drugs) {
      e <- ev(cmt = "ZAN_gut", amt = 160 * 0.60, ii = 12, addl = days * 2 - 1)
      d <- mrgsim(mod, events = e, end = days * 24, delta = 1, obsonly = TRUE) %>%
           as_tibble() %>% mutate(day = time / 24, Cp = Cp_ZAN_ng, drug = "Zanubrutinib (ng/mL)")
      plots_list[["ZAN"]] <- d
    }
    if ("RTX" %in% input$pk_drugs) {
      e <- ev(cmt = "RTX_C", amt = 375 * 1.8)
      d <- mrgsim(mod, events = e, end = days * 24, delta = 1, obsonly = TRUE) %>%
           as_tibble() %>% mutate(day = time / 24, Cp = Cp_RTX_mgl, drug = "Rituximab (mg/L)")
      plots_list[["RTX"]] <- d
    }

    if (length(plots_list) == 0) return(plotly_empty())

    df_all <- bind_rows(plots_list)
    p <- ggplot(df_all, aes(day, Cp, color = drug)) +
      geom_line(linewidth = 1.1) +
      labs(x = "Day", y = "Plasma Concentration",
           color = "Drug", title = "Drug PK Profiles") +
      theme_bw() + facet_wrap(~drug, scales = "free_y")
    ggplotly(p)
  })

  output$btk_occ_plot <- renderPlotly({
    days <- input$pk_days
    e_ibr <- ev(cmt = "IBR_gut", amt = 420 * 0.10, ii = 24, addl = days - 1)
    e_zan <- ev(cmt = "ZAN_gut", amt = 160 * 0.60, ii = 12, addl = days * 2 - 1)
    d_ibr <- mrgsim(mod, e_ibr, end = days * 24, delta = 2, obsonly = TRUE) %>%
             as_tibble() %>% mutate(day = time/24, drug = "Ibrutinib")
    d_zan <- mrgsim(mod, e_zan, end = days * 24, delta = 2, obsonly = TRUE) %>%
             as_tibble() %>% mutate(day = time/24, drug = "Zanubrutinib")
    df <- bind_rows(d_ibr, d_zan)
    p <- ggplot(df, aes(day, BTK_pct, color = drug)) +
      geom_line(linewidth = 1.2) +
      geom_hline(yintercept = 95, linetype = "dashed", color = "red") +
      labs(x = "Day", y = "BTK Occupancy (%)",
           title = "BTK Occupancy: Ibrutinib vs Zanubrutinib") +
      theme_bw() + scale_color_manual(values = c("#2980B9","#8E44AD"))
    ggplotly(p)
  })

  output$pk_table <- renderTable({
    data.frame(
      Parameter = c("Ibrutinib CL/F", "Ibrutinib V/F", "Ibrutinib T½",
                    "Zanubrutinib CL/F", "Zanubrutinib V/F",
                    "Rituximab CL", "Rituximab V",
                    "Venetoclax CL", "Venetoclax V"),
      Value = c("73 L/h", "820 L", "~4-6 h",
                "29 L/h", "520 L",
                "0.21 L/h", "4.5 L",
                "12 L/h", "256 L"),
      Source = c("Wang 2014", "", "",
                 "Tam 2020","",
                 "Tobinai 2011","",
                 "Roberts 2016","")
    )
  })

  # ── TAB 3: PD plots ──────────────────────────────────────────
  pd_data <- reactive({
    tryCatch(run_sim(input$pd_scenario, input$sim_days), error = function(e) NULL)
  })

  output$nfkb_plot <- renderPlotly({
    df <- pd_data(); if (is.null(df)) return(plotly_empty())
    df_long <- df %>%
      select(day, `NF-κB` = NFkB, `BCL-2` = BCL2, `BTK Occ (%)` = BTK_pct) %>%
      pivot_longer(-day)
    p <- ggplot(df_long, aes(day, value, color = name)) +
      geom_line(linewidth = 1.1) +
      facet_wrap(~name, scales = "free_y") +
      labs(x = "Day", y = "", title = paste("NF-κB/BCL-2/BTK Occupancy —", input$pd_scenario)) +
      theme_bw()
    ggplotly(p)
  })

  output$lpc_pc_plot <- renderPlotly({
    df <- pd_data(); if (is.null(df)) return(plotly_empty())
    df_long <- df %>% select(day, LPC, PC) %>% pivot_longer(-day)
    p <- ggplot(df_long, aes(day, value, color = name)) +
      geom_line(linewidth = 1.2) +
      labs(x = "Day", y = "Cells (×10⁹)", title = "LPC & PC Dynamics") +
      scale_color_manual(values = c("#8E44AD","#C0392B")) + theme_bw()
    ggplotly(p)
  })

  output$igm_prod_plot <- renderPlotly({
    df <- pd_data(); if (is.null(df)) return(plotly_empty())
    p <- ggplot(df, aes(day, IgM)) +
      geom_line(linewidth = 1.4, color = "#1A5276") +
      geom_hline(yintercept = input$IgM0 * 0.1, linetype = "dashed", color = "#27AE60") +
      geom_hline(yintercept = input$IgM0 * 0.5, linetype = "dashed", color = "#F39C12") +
      labs(x = "Day", y = "IgM (g/L)", title = "Serum IgM over Time") + theme_bw()
    ggplotly(p)
  })

  # ── TAB 4: Clinical ──────────────────────────────────────────
  clin_data <- reactive({
    tryCatch(run_sim(input$clin_scenario, input$sim_days), error = function(e) NULL)
  })

  output$clin_main_plot <- renderPlotly({
    df <- clin_data(); if (is.null(df)) return(plotly_empty())
    df_long <- df %>%
      select(day, `IgM (g/L)` = IgM, `Hemoglobin (g/dL)` = Hgb) %>%
      pivot_longer(-day)
    p <- ggplot(df_long, aes(day, value, color = name)) +
      geom_line(linewidth = 1.3) + facet_wrap(~name, scales = "free_y") +
      geom_hline(data = data.frame(name = "Hemoglobin (g/dL)", y = 10),
                 aes(yintercept = y), linetype = "dashed", color = "red") +
      labs(x = "Day", y = "", title = paste("Clinical response —", input$clin_scenario)) +
      theme_bw()
    ggplotly(p)
  })

  output$visc_plot <- renderPlotly({
    df <- clin_data(); if (is.null(df)) return(plotly_empty())
    p <- ggplot(df, aes(day, Visc)) +
      geom_line(linewidth = 1.3, color = "#C0392B") +
      geom_hline(yintercept = 3.5, linetype = "dashed", color = "red") +
      annotate("text", x = max(df$day) * 0.8, y = 3.7, label = "HVS threshold", color = "red") +
      labs(x = "Day", y = "Viscosity (cP)", title = "Serum Viscosity") + theme_bw()
    ggplotly(p)
  })

  output$bm_plot <- renderPlotly({
    df <- clin_data(); if (is.null(df)) return(plotly_empty())
    p <- ggplot(df, aes(day, BMInf * 100)) +
      geom_line(linewidth = 1.3, color = "#8E44AD") +
      geom_hline(yintercept = 10, linetype = "dashed", color = "gray40") +
      labs(x = "Day", y = "BM Infiltration (%)", title = "Bone Marrow Infiltration") +
      theme_bw()
    ggplotly(p)
  })

  # ── TAB 5: Scenario comparison ───────────────────────────────
  cmp_data <- reactive({
    scenarios <- input$cmp_scenarios
    if (length(scenarios) == 0) return(NULL)
    days <- input$cmp_days
    bind_rows(lapply(scenarios, function(s) {
      tryCatch(run_sim(s, days), error = function(e) NULL)
    }))
  })

  output$cmp_plot <- renderPlotly({
    df <- cmp_data(); if (is.null(df)) return(plotly_empty())
    yvar <- input$cmp_endpoint
    p <- ggplot(df, aes(day, .data[[yvar]], color = scenario, linetype = scenario)) +
      geom_line(linewidth = 1.1, alpha = 0.9) +
      scale_color_manual(values = pal) +
      labs(x = "Day", y = yvar,
           title = paste("Multi-scenario comparison:", yvar),
           color = "Scenario", linetype = "Scenario") +
      theme_bw() + theme(legend.position = "right")
    ggplotly(p)
  })

  output$response_dt <- renderDT({
    df <- cmp_data(); if (is.null(df)) return(NULL)
    df %>% filter(day >= 359, day <= 361) %>%
      group_by(scenario) %>% slice(1) %>% ungroup() %>%
      mutate(
        IgM_red_pct = round((input$IgM0 - IgM) / input$IgM0 * 100, 1),
        Response = case_when(
          IgM_red_pct >= 90 ~ "VGPR",
          IgM_red_pct >= 50 ~ "PR",
          IgM_red_pct >= 25 ~ "MR",
          IgM_red_pct >= 0  ~ "SD",
          TRUE ~ "PD"),
        IgM_gL2   = round(IgM, 1),
        Hgb_gL    = round(Hgb, 1),
        Visc_cP2  = round(Visc, 2),
        BMInf_pct2 = round(BMInf * 100, 1)
      ) %>%
      select(Scenario = scenario, `IgM (g/L)` = IgM_gL2,
             `IgM Reduction (%)` = IgM_red_pct,
             `Hgb (g/dL)` = Hgb_gL,
             `BM Inf (%)` = BMInf_pct2,
             `Viscosity (cP)` = Visc_cP2, Response)
  }, rownames = FALSE, options = list(pageLength = 8))

  # ── TAB 6: Biomarkers ────────────────────────────────────────
  bio_data <- reactive({
    tryCatch(run_sim(input$bio_scenario, input$sim_days), error = function(e) NULL)
  })

  output$bio_plot <- renderPlotly({
    df <- bio_data(); vars <- input$bio_vars
    if (is.null(df) || length(vars) == 0) return(plotly_empty())
    df_long <- df %>% select(day, all_of(vars)) %>% pivot_longer(-day)
    p <- ggplot(df_long, aes(day, value, color = name)) +
      geom_line(linewidth = 1.1) + facet_wrap(~name, scales = "free_y") +
      labs(x = "Day", y = "Value (AU/%)", color = "Biomarker",
           title = paste("Biomarker dynamics —", input$bio_scenario)) +
      theme_bw()
    ggplotly(p)
  })

  output$waterfall_plot <- renderPlotly({
    scen_list <- c("Watch & Wait", "Ibrutinib", "Ibrutinib + Rituximab",
                   "Zanubrutinib", "R-Bendamustine", "BDR", "Venetoclax")
    igm0 <- input$IgM0
    wf <- bind_rows(lapply(scen_list, function(s) {
      tryCatch({
        d <- run_sim(s, 180)
        igm_d180 <- d %>% filter(day >= 179) %>% pull(IgM) %>% last()
        tibble(Scenario = s, IgM_red = (igm0 - igm_d180) / igm0 * 100)
      }, error = function(e) tibble(Scenario = s, IgM_red = NA_real_))
    })) %>% filter(!is.na(IgM_red)) %>% arrange(desc(IgM_red)) %>%
      mutate(color_flag = ifelse(IgM_red > 0, "Response", "Progression"))

    p <- ggplot(wf, aes(x = reorder(Scenario, IgM_red), y = IgM_red, fill = color_flag)) +
      geom_col() + coord_flip() +
      geom_hline(yintercept = c(25, 50, 90), linetype = "dashed", color = c("#F39C12","#27AE60","#1A5276")) +
      scale_fill_manual(values = c("Response" = "#2980B9", "Progression" = "#E74C3C")) +
      labs(x = NULL, y = "IgM Reduction (%)", title = "Waterfall: IgM reduction at 6 months") +
      theme_bw()
    ggplotly(p)
  })

  output$bio_table <- renderDT({
    df <- bio_data(); if (is.null(df)) return(NULL)
    df %>% filter(day %in% c(0, 90, 180, 360)) %>%
      select(Day = day, IgM, Hgb, BTK_pct, NFkB, BCL2, NK, BMInf) %>%
      mutate(across(where(is.numeric), ~ round(.x, 2)),
             BMInf = round(BMInf * 100, 1))
  }, rownames = FALSE, options = list(pageLength = 4))
}

shinyApp(ui, server)
