## ============================================================
## Aplastic Anemia (AA) QSP — Interactive Shiny Dashboard
## ============================================================
## Author  : Claude Code Routine (CCR)
## Date    : 2026-06-25
##
## Tabs:
##  1. Patient Profile & Disease Severity
##  2. Drug PK (ATG / CsA / EPAG)
##  3. Hematopoiesis & HSC Dynamics
##  4. Clinical Endpoints (Hgb, ANC, PLT)
##  5. Treatment Scenario Comparison
##  6. Biomarkers & Clonal Evolution (IFN-γ, PNH clone)
## ============================================================

library(shiny)
library(shinydashboard)
library(mrgsolve)
library(dplyr)
library(ggplot2)
library(tidyr)
library(purrr)
library(plotly)

## ── Load / inline the mrgsolve model code ──────────────────
aa_model_code <- '
$PARAM
ATG_CL=0.85, ATG_Vc=5.8, ATG_Vp=12.0, ATG_Q=0.30
ATG_Emax=0.92, ATG_EC50=0.15, ATG_hill=1.5
CsA_ka=0.55, CsA_F=0.34, CsA_Vc=4.5, CsA_CL=0.28
CsA_EC50=165, CsA_Emax=0.85
EPAG_ka=0.90, EPAG_F=0.52, EPAG_Vc=38.0, EPAG_CL=1.30
EPAG_EC50=60.0, EPAG_Emax=0.78
Danazol_ka=0.40, Danazol_F=0.20, Danazol_Vc=200.0, Danazol_CL=18.0
kTeff_prod=0.12, kTeff_death=0.08, kTreg_prod=0.03, kTreg_death=0.06
kHSC_self=0.04, kHSC_diff=0.06, kHSC_death=0.002
kHSC_IFNg=0.18, kHSC_TNFa=0.08, kHSC_EPAG=0.10
kCFUE_prod=0.30, kCFUE_mat=0.18, kRetic_mat=0.25, kRBC_death=0.0083
EPO_basal=15.0, EPO_max=200.0, EPO_k50=9.5, Hgb_normal=14.0
kCFUG_prod=0.35, kCFUG_mat=0.20, kANC_death=3.0, ANC_normal=3.0
kMK_prod=0.15, kMK_mat=0.06, kPLT_prod=4.0, kPLT_death=0.10, PLT_normal=200.0
kBM_repair=0.02, kBM_damage=0.05
kIFNg_prod=2.5, kIFNg_CL=1.2
kTNFa_prod=1.0, kTNFa_CL=2.0
kIL2_prod=0.8, kIL2_CL=3.0
kPNH_growth=0.005, kPNH_max=0.60
use_ATG=0, use_CsA=0, use_EPAG=0, use_Danazol=0, use_HSCT=0
severity=1.0

$CMT ATG_C ATG_P CsA_C EPAG_C Danazol_C
     Teff Treg HSC CFU_E Retic RBC CFU_G ANC_pool MK PLT_pool
     BM_score IFNg_c TNFa_c IL2_c PNH_clone

$MAIN
double ATG_Eff = (use_ATG > 0.5) ?
  ATG_Emax * pow(ATG_C, ATG_hill) / (pow(ATG_EC50, ATG_hill) + pow(ATG_C, ATG_hill) + 1e-12) : 0.0;
double CsA_Eff = (use_CsA > 0.5) ?
  CsA_Emax * CsA_C / (CsA_EC50 + CsA_C + 1e-12) : 0.0;
double EPAG_Eff = (use_EPAG > 0.5) ?
  EPAG_Emax * EPAG_C / (EPAG_EC50 + EPAG_C + 1e-12) : 0.0;
double Danazol_Eff = (use_Danazol > 0.5) ?
  0.30 * Danazol_C / (400.0 + Danazol_C + 1e-12) : 0.0;
double HSCT_Eff = use_HSCT;
double Hgb_cur = (RBC / 100.0) * Hgb_normal;
double EPO_fold = 1.0 + (EPO_max/EPO_basal - 1.0) *
                  pow((1.0 - Hgb_cur/Hgb_normal), 2.0) /
                  (pow(EPO_k50/Hgb_normal, 2.0) + pow((1.0 - Hgb_cur/Hgb_normal), 2.0) + 1e-6);
double IFNg_unit = IFNg_c / 100.0;
double TNFa_unit = TNFa_c / 100.0;
double Treg_suppress = Treg / (Treg + 5.0);
double HSC_kill = (kHSC_IFNg * IFNg_unit + kHSC_TNFa * TNFa_unit) *
                  (HSC / 100.0) * severity * (1.0 - 0.5 * Treg_suppress);

$INIT
ATG_C=0, ATG_P=0, CsA_C=0, EPAG_C=0, Danazol_C=0
Teff=10.0, Treg=1.0, HSC=20.0
CFU_E=25.0, Retic=20.0, RBC=55.0
CFU_G=15.0, ANC_pool=10.0
MK=12.0, PLT_pool=8.0
BM_score=0.12, IFNg_c=45.0, TNFa_c=25.0, IL2_c=20.0, PNH_clone=0.02

$ODE
dxdt_ATG_C   = -ATG_CL/ATG_Vc*ATG_C - ATG_Q/ATG_Vc*ATG_C + ATG_Q/ATG_Vp*ATG_P;
dxdt_ATG_P   =  ATG_Q/ATG_Vc*ATG_C - ATG_Q/ATG_Vp*ATG_P;
dxdt_CsA_C   = -CsA_CL/CsA_Vc*CsA_C;
dxdt_EPAG_C  = -EPAG_CL/EPAG_Vc*EPAG_C;
dxdt_Danazol_C = -Danazol_CL/Danazol_Vc*Danazol_C;

double Teff_prolif = kTeff_prod*Teff*(1.0+(IL2_c/50.0))*(1.0-CsA_Eff)*(1.0-ATG_Eff);
double Teff_death  = (kTeff_death + ATG_Eff*0.5)*Teff;
double Treg_inh    = 0.3*(Treg/(Treg+2.0))*Teff;
dxdt_Teff = Teff_prolif - Teff_death - Treg_inh - HSCT_Eff*5.0*Teff;
dxdt_Treg = kTreg_prod*(1.0+0.5*ATG_Eff)*5.0 - kTreg_death*Treg - HSCT_Eff*0.5*Treg;

double HSC_selfR = kHSC_self*HSC*(1.0+EPAG_Eff*kHSC_EPAG/kHSC_self);
dxdt_HSC = (HSC<105.0) ? (HSC_selfR - kHSC_diff*HSC - kHSC_death*HSC - HSC_kill) : 0.0;

dxdt_CFU_E  = kCFUE_prod*(HSC/100.0)*EPO_fold - kCFUE_mat*CFU_E - kHSC_IFNg*0.5*IFNg_unit*CFU_E;
dxdt_Retic  = kCFUE_mat*CFU_E - kRetic_mat*Retic;
dxdt_RBC    = kRetic_mat*Retic*(1.0+Danazol_Eff) - kRBC_death*RBC;
dxdt_CFU_G  = kCFUG_prod*(HSC/100.0) - kCFUG_mat*CFU_G;
dxdt_ANC_pool = kCFUG_mat*CFU_G - kANC_death*ANC_pool;
dxdt_MK     = kMK_prod*(HSC/100.0)*(1.0+EPAG_Eff) - kMK_mat*MK;
dxdt_PLT_pool = kMK_mat*MK*kPLT_prod - kPLT_death*PLT_pool;
dxdt_BM_score = kBM_repair*(HSC-20.0)/100.0 - kBM_damage*(IFNg_unit+TNFa_unit*0.5)*BM_score;
dxdt_IFNg_c = kIFNg_prod*Teff - kIFNg_CL*IFNg_c;
dxdt_TNFa_c = kTNFa_prod*Teff - kTNFa_CL*TNFa_c;
dxdt_IL2_c  = kIL2_prod*Teff*(1.0-CsA_Eff) - kIL2_CL*IL2_c;
dxdt_PNH_clone = kPNH_growth*PNH_clone*(1.0-PNH_clone) - ATG_Eff*0.1*PNH_clone;

$TABLE
capture Hgb      = (RBC/100.0)*Hgb_normal;
capture ANC      = (ANC_pool/100.0)*ANC_normal;
capture PLT      = (PLT_pool/100.0)*PLT_normal;
capture ARC      = (Retic/100.0)*80.0;
capture BM_cell  = BM_score*100.0;
capture CR_flag  = (Hgb>11.0 && ANC>1.0 && PLT>100.0) ? 1.0 : 0.0;
capture PR_flag  = (Hgb>=8.0 && PLT>20.0 && ANC>0.5)  ? 1.0 : 0.0;
capture ATG_conc = ATG_C;
capture CsA_conc = CsA_C;
capture EPAG_conc= EPAG_C;
'

aa_mod <- mcode("aa_shiny", aa_model_code, quiet = TRUE)

## ── Helper: simulate one scenario ──────────────────────────
simulate_scenario <- function(params, duration_days = 365) {
  mod2 <- param(aa_mod,
    use_ATG     = as.numeric(params$use_ATG),
    use_CsA     = as.numeric(params$use_CsA),
    use_EPAG    = as.numeric(params$use_EPAG),
    use_Danazol = as.numeric(params$use_Danazol),
    severity    = params$severity,
    ATG_Emax    = params$ATG_Emax,
    EPAG_EC50   = params$EPAG_EC50,
    CsA_EC50    = params$CsA_EC50
  )
  mod2 <- init(mod2,
    HSC      = params$HSC_init,
    RBC      = params$HSC_init * 2.5,
    ANC_pool = params$HSC_init * 0.5,
    PLT_pool = params$HSC_init * 0.4
  )

  evlist <- list()
  if (params$use_ATG) {
    n_days <- if (params$ATG_type == "hATG") 4 else 5
    dose_per_day <- if (params$ATG_type == "hATG") 40 * 70 else 3.5 * 70
    evlist$atg <- ev(
      time = 0:(n_days - 1) * 24,
      amt  = dose_per_day,
      rate = dose_per_day / 12,
      cmt  = "ATG_C"
    )
  }
  if (params$use_CsA) {
    evlist$csa <- ev(
      time = 0:(duration_days - 1) * 24,
      amt  = 300,
      cmt  = "CsA_C",
      rate = -2
    )
  }
  if (params$use_EPAG) {
    start_h <- params$EPAG_start_day * 24
    epag_days <- seq(start_h, start_h + (params$EPAG_duration - 1) * 24, by = 24)
    evlist$epag <- ev(
      time = epag_days,
      amt  = params$EPAG_dose,
      cmt  = "EPAG_C",
      rate = -2
    )
  }
  if (params$use_Danazol) {
    evlist$danazol <- ev(
      time = 0:(duration_days - 1) * 24,
      amt  = 400,  # mg/d total
      cmt  = "Danazol_C",
      rate = -2
    )
  }

  all_ev <- if (length(evlist) > 0) do.call(c, evlist) else NULL

  out <- if (!is.null(all_ev)) {
    mrgsim(mod2, events = all_ev, end = duration_days * 24, delta = 24)
  } else {
    mrgsim(mod2, end = duration_days * 24, delta = 24)
  }
  as.data.frame(out) %>%
    mutate(time_days = time / 24)
}

## ============================================================
## UI
## ============================================================
ui <- dashboardPage(
  skin = "red",
  dashboardHeader(title = "Aplastic Anemia QSP"),

  dashboardSidebar(
    sidebarMenu(
      menuItem("Patient Profile",        tabName = "tab_patient",  icon = icon("user")),
      menuItem("Drug PK",                tabName = "tab_pk",       icon = icon("pills")),
      menuItem("Hematopoiesis",          tabName = "tab_hemato",   icon = icon("dna")),
      menuItem("Clinical Endpoints",     tabName = "tab_clinical", icon = icon("heartbeat")),
      menuItem("Scenario Comparison",    tabName = "tab_scenario", icon = icon("chart-bar")),
      menuItem("Biomarkers & Clones",    tabName = "tab_bio",      icon = icon("microscope"))
    ),
    hr(),
    h5("Disease Parameters", style = "margin-left:10px; color:white;"),
    sliderInput("severity", "Disease Severity",
                min = 0.5, max = 1.5, value = 1.0, step = 0.1),
    sliderInput("HSC_init", "Initial HSC Pool (%)",
                min = 5, max = 40, value = 20, step = 1),
    sliderInput("duration_days", "Simulation Duration (days)",
                min = 90, max = 730, value = 365, step = 30),
    hr(),
    h5("Treatment Selection", style = "margin-left:10px; color:white;"),
    checkboxInput("use_ATG",     "ATG (IST)",      value = FALSE),
    conditionalPanel("input.use_ATG",
      radioButtons("ATG_type", "ATG Type:",
                   choices = c("hATG (ATGAM)" = "hATG", "rATG (Thymoglobulin)" = "rATG"),
                   selected = "hATG", inline = TRUE)
    ),
    checkboxInput("use_CsA",     "Cyclosporine",   value = FALSE),
    checkboxInput("use_EPAG",    "Eltrombopag",    value = FALSE),
    conditionalPanel("input.use_EPAG",
      sliderInput("EPAG_dose",  "EPAG Dose (mg/d)", min = 75, max = 200, value = 150, step = 25),
      sliderInput("EPAG_start", "EPAG Start (day)", min = 1,  max = 60,  value = 14,  step = 1),
      sliderInput("EPAG_dur",   "EPAG Duration (d)",min = 30, max = 360, value = 180, step = 30)
    ),
    checkboxInput("use_Danazol", "Danazol (androgen)", value = FALSE),
    actionButton("run_sim", "Run Simulation", icon = icon("play"),
                 style = "color:#fff; background:#e74c3c; border:none; margin:10px; width:90%;")
  ),

  dashboardBody(
    tags$head(tags$style(HTML("
      .skin-red .main-header .logo { background:#922B21; }
      .skin-red .main-header .navbar { background:#C0392B; }
      .content-wrapper { background:#f4f6f8; }
      .box { border-radius:8px; }
    "))),

    tabItems(

      ## ── TAB 1: PATIENT PROFILE ──────────────────────────
      tabItem("tab_patient",
        fluidRow(
          valueBoxOutput("vbox_hgb",  width = 3),
          valueBoxOutput("vbox_anc",  width = 3),
          valueBoxOutput("vbox_plt",  width = 3),
          valueBoxOutput("vbox_bm",   width = 3)
        ),
        fluidRow(
          box(title = "Severity Classification", width = 4, status = "danger",
            tableOutput("tbl_severity")
          ),
          box(title = "Initial State Summary", width = 8, status = "primary",
            plotOutput("plot_radar", height = 300)
          )
        ),
        fluidRow(
          box(title = "Disease Mechanism Overview", width = 12, status = "info",
            column(6,
              h4("Pathophysiology"),
              tags$ul(
                tags$li("T-cell mediated destruction of HSCs via IFN-γ / FasL / Perforin"),
                tags$li("Treg deficiency allows autoreactive clone expansion"),
                tags$li("BM microenvironment injury: MSC dysfunction, adipose replacement"),
                tags$li("PNH clone may expand via immune escape (GPI-anchor deficiency)")
              ),
              h4("Treatment Approach"),
              tags$ol(
                tags$li("IST: ATG + CsA ± EPAG — for non-HSCT candidates"),
                tags$li("HSCT: curative for young patients with MSD donor"),
                tags$li("Supportive: transfusions, G-CSF, antimicrobial prophylaxis")
              )
            ),
            column(6,
              h4("Severity Criteria (Camitta, 1976; modified)"),
              tags$table(class = "table table-bordered table-sm",
                tags$thead(tags$tr(
                  tags$th("Category"), tags$th("ANC"), tags$th("PLT"), tags$th("ARC or BM")
                )),
                tags$tbody(
                  tags$tr(tags$td("nSAA"), tags$td("0.5–1.5"), tags$td(">20"), tags$td("any")),
                  tags$tr(tags$td("SAA"),  tags$td("<0.5"),    tags$td("<20 or Hgb<8"), tags$td("<40 K/μL")),
                  tags$tr(tags$td("VSAA"), tags$td("<0.2"),    tags$td("any"),          tags$td("any"))
                )
              )
            )
          )
        )
      ),

      ## ── TAB 2: DRUG PK ──────────────────────────────────
      tabItem("tab_pk",
        fluidRow(
          box(title = "ATG Concentration–Time Profile", width = 6, status = "danger",
            plotlyOutput("plot_pk_atg", height = 280)
          ),
          box(title = "Cyclosporine Trough Trend", width = 6, status = "warning",
            plotlyOutput("plot_pk_csa", height = 280)
          )
        ),
        fluidRow(
          box(title = "Eltrombopag Concentration", width = 6, status = "success",
            plotlyOutput("plot_pk_epag", height = 280)
          ),
          box(title = "PK Parameters Summary", width = 6, status = "info",
            tableOutput("tbl_pk_params")
          )
        ),
        fluidRow(
          box(title = "ATG PD Effect: T-Cell Depletion", width = 12, status = "primary",
            plotlyOutput("plot_pk_pd", height = 280)
          )
        )
      ),

      ## ── TAB 3: HEMATOPOIESIS ────────────────────────────
      tabItem("tab_hemato",
        fluidRow(
          box(title = "HSC Pool Recovery (% of Normal)", width = 6, status = "primary",
            plotlyOutput("plot_hsc", height = 280)
          ),
          box(title = "BM Cellularity (%)", width = 6, status = "info",
            plotlyOutput("plot_bm", height = 280)
          )
        ),
        fluidRow(
          box(title = "Erythroid Lineage (CFU-E → Retic → RBC)", width = 6, status = "danger",
            plotlyOutput("plot_erythro", height = 280)
          ),
          box(title = "Myeloid & Megakaryocyte Lineage", width = 6, status = "warning",
            plotlyOutput("plot_myeloid", height = 280)
          )
        ),
        fluidRow(
          box(title = "Effector T Cell vs Treg Dynamics", width = 12, status = "success",
            plotlyOutput("plot_tcell", height = 280)
          )
        )
      ),

      ## ── TAB 4: CLINICAL ENDPOINTS ───────────────────────
      tabItem("tab_clinical",
        fluidRow(
          box(title = "Hemoglobin (g/dL)", width = 4, status = "danger",
            plotlyOutput("plot_hgb", height = 250)
          ),
          box(title = "ANC (×10⁹/L)", width = 4, status = "warning",
            plotlyOutput("plot_anc", height = 250)
          ),
          box(title = "Platelet Count (×10⁹/L)", width = 4, status = "info",
            plotlyOutput("plot_plt", height = 250)
          )
        ),
        fluidRow(
          box(title = "Absolute Reticulocyte Count (×10⁹/L)", width = 6, status = "success",
            plotlyOutput("plot_arc", height = 250)
          ),
          box(title = "Response Classification Over Time", width = 6, status = "primary",
            plotlyOutput("plot_response", height = 250)
          )
        ),
        fluidRow(
          box(title = "Transfusion Need Probability", width = 12, status = "danger",
            plotlyOutput("plot_transfusion", height = 220)
          )
        )
      ),

      ## ── TAB 5: SCENARIO COMPARISON ──────────────────────
      tabItem("tab_scenario",
        fluidRow(
          box(title = "Scenario Definitions", width = 12, status = "primary",
            tableOutput("tbl_scenarios")
          )
        ),
        fluidRow(
          box(title = "Hemoglobin — All Scenarios", width = 6, status = "danger",
            plotlyOutput("plot_sc_hgb", height = 280)
          ),
          box(title = "ANC — All Scenarios", width = 6, status = "warning",
            plotlyOutput("plot_sc_anc", height = 280)
          )
        ),
        fluidRow(
          box(title = "Platelet — All Scenarios", width = 6, status = "info",
            plotlyOutput("plot_sc_plt", height = 280)
          ),
          box(title = "HSC Recovery — All Scenarios", width = 6, status = "success",
            plotlyOutput("plot_sc_hsc", height = 280)
          )
        ),
        fluidRow(
          box(title = "Response Rate at Day 180 (Model Prediction)", width = 12, status = "primary",
            plotlyOutput("plot_sc_response", height = 300)
          )
        )
      ),

      ## ── TAB 6: BIOMARKERS & CLONES ──────────────────────
      tabItem("tab_bio",
        fluidRow(
          box(title = "IFN-γ Concentration Over Time", width = 6, status = "warning",
            plotlyOutput("plot_ifng", height = 280)
          ),
          box(title = "TNF-α and IL-2", width = 6, status = "danger",
            plotlyOutput("plot_tnf_il2", height = 280)
          )
        ),
        fluidRow(
          box(title = "PNH Clone Size (Fraction)", width = 6, status = "info",
            plotlyOutput("plot_pnh", height = 280)
          ),
          box(title = "Biomarker Reference Ranges", width = 6, status = "primary",
            tableOutput("tbl_biomarker")
          )
        ),
        fluidRow(
          box(title = "Clonal Evolution Risk Summary", width = 12, status = "danger",
            column(6,
              h4("PNH Clone Dynamics"),
              p("PNH clones (GPI-anchor deficient cells) expand in AA due to relative",
                "immune escape. 50–60% of SAA patients have detectable PNH clones at",
                "diagnosis (>1% by flow cytometry). Clone size may expand during IST."),
              h4("Risk of Transformation"),
              tags$ul(
                tags$li("MDS/AML transformation: 10–15% at 10 years"),
                tags$li("Monosomy 7 most common cytogenetic event"),
                tags$li("Supported by EPAG (monitor PLT rise)"),
                tags$li("Annual cytogenetics recommended for high-risk patients")
              )
            ),
            column(6,
              plotlyOutput("plot_clone_risk", height = 250)
            )
          )
        )
      )
    )
  )
)

## ============================================================
## SERVER
## ============================================================
server <- function(input, output, session) {

  ## ── Reactive simulation ─────────────────────────────────
  sim_result <- eventReactive(input$run_sim, {
    params <- list(
      use_ATG      = input$use_ATG,
      use_CsA      = input$use_CsA,
      use_EPAG     = input$use_EPAG,
      use_Danazol  = input$use_Danazol,
      use_HSCT     = FALSE,
      severity     = input$severity,
      HSC_init     = input$HSC_init,
      ATG_type     = if (input$use_ATG) input$ATG_type else "hATG",
      ATG_Emax     = 0.92,
      EPAG_EC50    = 60.0,
      EPAG_dose    = if (input$use_EPAG) input$EPAG_dose else 150,
      EPAG_start_day = if (input$use_EPAG) input$EPAG_start else 14,
      EPAG_duration = if (input$use_EPAG) input$EPAG_dur else 180,
      CsA_EC50     = 165
    )
    simulate_scenario(params, duration_days = input$duration_days)
  }, ignoreNULL = FALSE)

  ## default on load
  default_result <- reactive({
    params <- list(
      use_ATG = FALSE, use_CsA = FALSE, use_EPAG = FALSE,
      use_Danazol = FALSE, use_HSCT = FALSE,
      severity = 1.0, HSC_init = 20.0, ATG_type = "hATG",
      ATG_Emax = 0.92, EPAG_EC50 = 60.0, EPAG_dose = 150,
      EPAG_start_day = 14, EPAG_duration = 180, CsA_EC50 = 165
    )
    simulate_scenario(params, duration_days = 365)
  })

  res <- reactive({
    if (input$run_sim == 0) default_result() else sim_result()
  })

  ## ── Scenario comparison data ─────────────────────────────
  scenario_data <- reactive({
    scen_defs <- list(
      list(label = "No Treatment",           use_ATG = FALSE, use_CsA = FALSE, use_EPAG = FALSE),
      list(label = "hATG + CsA",             use_ATG = TRUE,  use_CsA = TRUE,  use_EPAG = FALSE, ATG_type = "hATG"),
      list(label = "hATG + CsA + EPAG",      use_ATG = TRUE,  use_CsA = TRUE,  use_EPAG = TRUE,  ATG_type = "hATG"),
      list(label = "rATG + CsA + EPAG",      use_ATG = TRUE,  use_CsA = TRUE,  use_EPAG = TRUE,  ATG_type = "rATG"),
      list(label = "EPAG Monotherapy",        use_ATG = FALSE, use_CsA = FALSE, use_EPAG = TRUE,  ATG_type = "hATG")
    )
    map_dfr(scen_defs, function(s) {
      params <- list(
        use_ATG = s$use_ATG, use_CsA = s$use_CsA, use_EPAG = s$use_EPAG,
        use_Danazol = FALSE, use_HSCT = FALSE,
        severity = input$severity, HSC_init = input$HSC_init,
        ATG_type = if (!is.null(s$ATG_type)) s$ATG_type else "hATG",
        ATG_Emax = 0.92, EPAG_EC50 = 60.0, EPAG_dose = 150,
        EPAG_start_day = 14, EPAG_duration = 180, CsA_EC50 = 165
      )
      simulate_scenario(params, 365) %>% mutate(Scenario = s$label)
    })
  })

  ## ── Tab 1: Value Boxes ────────────────────────────────────
  final_vals <- reactive({
    r <- res(); r[nrow(r), ]
  })
  output$vbox_hgb <- renderValueBox({
    v <- round(final_vals()$Hgb, 1)
    col <- if (v >= 11) "green" else if (v >= 8) "yellow" else "red"
    valueBox(paste(v, "g/dL"), "Hemoglobin", icon = icon("tint"), color = col)
  })
  output$vbox_anc <- renderValueBox({
    v <- round(final_vals()$ANC, 2)
    col <- if (v >= 1.0) "green" else if (v >= 0.5) "yellow" else "red"
    valueBox(paste(v, "×10⁹/L"), "ANC", icon = icon("shield-alt"), color = col)
  })
  output$vbox_plt <- renderValueBox({
    v <- round(final_vals()$PLT, 0)
    col <- if (v >= 100) "green" else if (v >= 20) "yellow" else "red"
    valueBox(paste(v, "×10⁹/L"), "Platelet Count", icon = icon("circle"), color = col)
  })
  output$vbox_bm <- renderValueBox({
    v <- round(final_vals()$BM_cell, 0)
    col <- if (v >= 30) "green" else if (v >= 15) "yellow" else "red"
    valueBox(paste(v, "%"), "BM Cellularity", icon = icon("microscope"), color = col)
  })

  output$tbl_severity <- renderTable({
    data.frame(
      Category  = c("VSAA", "SAA", "nSAA"),
      ANC_cutoff = c("<0.2", "0.2–0.5", "0.5–1.5"),
      PLT_cutoff = c("<20 K", "<20 K", ">20 K"),
      Required   = c("Any 2", "Any 2 of 3", "Mild")
    )
  }, striped = TRUE, bordered = TRUE)

  output$plot_radar <- renderPlot({
    r0 <- res()[1, ]
    df <- data.frame(
      Variable = c("Hgb/14", "ANC/3", "PLT/200", "HSC%", "BM%"),
      Value    = c(r0$Hgb / 14, r0$ANC / 3, r0$PLT / 200,
                   r0$HSC / 100, r0$BM_cell / 100) * 100
    )
    ggplot(df, aes(x = Variable, y = Value)) +
      geom_col(fill = "#E74C3C", alpha = 0.7) +
      geom_hline(yintercept = 100, linetype = "dashed", color = "darkgreen") +
      ylim(0, 110) +
      labs(title = "Initial Values (% of Normal)", y = "% Normal", x = "") +
      theme_minimal(base_size = 13)
  })

  ## ── Tab 2: PK Plots ──────────────────────────────────────
  make_ply <- function(p) ggplotly(p, tooltip = c("x", "y"))

  output$plot_pk_atg <- renderPlotly({
    make_ply(
      ggplot(res(), aes(x = time_days, y = ATG_conc)) +
        geom_area(fill = "#E74C3C", alpha = 0.5) + geom_line(color = "#922B21") +
        labs(x = "Day", y = "ATG (mg/L)") + theme_bw()
    )
  })
  output$plot_pk_csa <- renderPlotly({
    make_ply(
      ggplot(res(), aes(x = time_days, y = CsA_conc)) +
        geom_line(color = "#F39C12", linewidth = 1.0) +
        geom_hline(yintercept = 165, linetype = "dashed") +
        labs(x = "Day", y = "CsA (ng/mL)") + theme_bw()
    )
  })
  output$plot_pk_epag <- renderPlotly({
    make_ply(
      ggplot(res(), aes(x = time_days, y = EPAG_conc)) +
        geom_area(fill = "#2ECC71", alpha = 0.4) + geom_line(color = "#1A8A4A") +
        geom_hline(yintercept = 60, linetype = "dashed") +
        labs(x = "Day", y = "EPAG (μg/mL)") + theme_bw()
    )
  })
  output$tbl_pk_params <- renderTable({
    data.frame(
      Drug    = c("hATG", "rATG", "CsA", "EPAG", "Danazol"),
      Dose    = c("40 mg/kg×4d", "3.5 mg/kg×5d", "5 mg/kg/d", "150 mg/d", "400-800 mg/d"),
      Route   = c("IV", "IV", "PO BID", "PO QD", "PO"),
      t_half  = c("~7 h", "~12 h", "~27 h", "~21 h", "~3 h"),
      Target  = c("T cells", "T cells", "150-250 ng/mL", "≥70 μg/mL", "Clinical")
    )
  }, striped = TRUE, bordered = TRUE)

  output$plot_pk_pd <- renderPlotly({
    r <- res() %>% mutate(
      ATG_eff = 0.92 * ATG_conc^1.5 / (0.15^1.5 + ATG_conc^1.5 + 1e-9)
    )
    make_ply(
      ggplot(r, aes(x = time_days, y = ATG_eff * 100)) +
        geom_line(color = "#8E44AD", linewidth = 1.1) +
        labs(x = "Day", y = "T-Cell Depletion Effect (%)") + theme_bw()
    )
  })

  ## ── Tab 3: Hematopoiesis ─────────────────────────────────
  output$plot_hsc <- renderPlotly({
    make_ply(
      ggplot(res(), aes(x = time_days, y = HSC)) +
        geom_area(fill = "#3498DB", alpha = 0.4) + geom_line(color = "#1A5276") +
        geom_hline(yintercept = 100, linetype = "dashed", color = "darkgreen") +
        labs(x = "Day", y = "HSC Pool (% Normal)") + theme_bw()
    )
  })
  output$plot_bm <- renderPlotly({
    make_ply(
      ggplot(res(), aes(x = time_days, y = BM_cell)) +
        geom_line(color = "#8E44AD", linewidth = 1.1) +
        geom_hline(yintercept = 25, linetype = "dashed") +
        labs(x = "Day", y = "BM Cellularity (%)") + theme_bw()
    )
  })
  output$plot_erythro <- renderPlotly({
    r <- res() %>% select(time_days, CFU_E, Retic, RBC) %>%
      pivot_longer(-time_days, names_to = "Compartment", values_to = "Value")
    make_ply(
      ggplot(r, aes(x = time_days, y = Value, color = Compartment)) +
        geom_line(linewidth = 1.0) +
        scale_color_manual(values = c(CFU_E = "#E74C3C", Retic = "#F39C12", RBC = "#C0392B")) +
        labs(x = "Day", y = "% of Normal") + theme_bw()
    )
  })
  output$plot_myeloid <- renderPlotly({
    r <- res() %>% select(time_days, CFU_G, ANC_pool, MK, PLT_pool) %>%
      pivot_longer(-time_days, names_to = "Compartment", values_to = "Value")
    make_ply(
      ggplot(r, aes(x = time_days, y = Value, color = Compartment)) +
        geom_line(linewidth = 1.0) +
        labs(x = "Day", y = "% of Normal") + theme_bw()
    )
  })
  output$plot_tcell <- renderPlotly({
    r <- res() %>% select(time_days, Teff, Treg) %>%
      pivot_longer(-time_days, names_to = "Cell", values_to = "Count")
    make_ply(
      ggplot(r, aes(x = time_days, y = Count, color = Cell)) +
        geom_line(linewidth = 1.2) +
        scale_color_manual(values = c(Teff = "#E74C3C", Treg = "#2ECC71")) +
        labs(x = "Day", y = "Cell Count (×10⁶/kg)") + theme_bw()
    )
  })

  ## ── Tab 4: Clinical Endpoints ────────────────────────────
  output$plot_hgb <- renderPlotly({
    make_ply(
      ggplot(res(), aes(x = time_days, y = Hgb)) +
        geom_area(fill = "#E74C3C", alpha = 0.4) + geom_line(color = "#922B21", linewidth = 1) +
        geom_hline(yintercept = c(8, 11), linetype = "dashed", color = c("red", "green")) +
        labs(x = "Day", y = "Hgb (g/dL)") + theme_bw()
    )
  })
  output$plot_anc <- renderPlotly({
    make_ply(
      ggplot(res(), aes(x = time_days, y = ANC)) +
        geom_area(fill = "#F39C12", alpha = 0.4) + geom_line(color = "#B7770D", linewidth = 1) +
        geom_hline(yintercept = c(0.2, 0.5, 1.0), linetype = "dashed") +
        labs(x = "Day", y = "ANC (×10⁹/L)") + theme_bw()
    )
  })
  output$plot_plt <- renderPlotly({
    make_ply(
      ggplot(res(), aes(x = time_days, y = PLT)) +
        geom_area(fill = "#3498DB", alpha = 0.4) + geom_line(color = "#1A5276", linewidth = 1) +
        geom_hline(yintercept = c(10, 20, 100), linetype = "dashed") +
        labs(x = "Day", y = "PLT (×10⁹/L)") + theme_bw()
    )
  })
  output$plot_arc <- renderPlotly({
    make_ply(
      ggplot(res(), aes(x = time_days, y = ARC)) +
        geom_line(color = "#2ECC71", linewidth = 1.1) +
        geom_hline(yintercept = 60, linetype = "dashed") +
        labs(x = "Day", y = "ARC (×10⁹/L)") + theme_bw()
    )
  })
  output$plot_response <- renderPlotly({
    r <- res() %>%
      mutate(Response = case_when(
        CR_flag == 1 ~ "CR",
        PR_flag == 1 ~ "PR",
        TRUE         ~ "NR"
      ))
    make_ply(
      ggplot(r, aes(x = time_days, fill = Response)) +
        geom_histogram(binwidth = 14, position = "stack") +
        scale_fill_manual(values = c(CR = "#2ECC71", PR = "#F39C12", NR = "#E74C3C")) +
        labs(x = "Day", y = "Count") + theme_bw()
    )
  })
  output$plot_transfusion <- renderPlotly({
    r <- res() %>% mutate(
      pRBC_need = as.numeric(Hgb < 8.0),
      PLT_need  = as.numeric(PLT < 10.0)
    ) %>% select(time_days, pRBC_need, PLT_need) %>%
      pivot_longer(-time_days, names_to = "Type", values_to = "Need")
    make_ply(
      ggplot(r, aes(x = time_days, y = Need, color = Type)) +
        geom_step(linewidth = 1.0) +
        scale_color_manual(values = c(pRBC_need = "#E74C3C", PLT_need = "#3498DB"),
                           labels = c("pRBC transfusion", "Platelet transfusion")) +
        labs(x = "Day", y = "Transfusion Needed (1=Yes)") + theme_bw()
    )
  })

  ## ── Tab 5: Scenario Comparison ──────────────────────────
  output$tbl_scenarios <- renderTable({
    data.frame(
      Scenario = c("No Treatment", "hATG+CsA", "hATG+CsA+EPAG", "rATG+CsA+EPAG", "EPAG mono"),
      ATG      = c("—", "hATG", "hATG", "rATG", "—"),
      CsA      = c("—", "Yes", "Yes", "Yes", "—"),
      EPAG     = c("—", "—", "Day 14", "Day 14", "Day 1"),
      `Expected CR (%)` = c("~10", "~50", "~68", "~74", "~35"),
      Reference = c("Natural history", "Scheinberg 2011", "Townsley 2017",
                    "Peffault 2021", "Olnes 2012")
    )
  }, striped = TRUE, bordered = TRUE)

  make_scen_ply <- function(var, ylab) {
    cols <- c(
      "No Treatment"     = "#E74C3C",
      "hATG + CsA"       = "#3498DB",
      "hATG + CsA + EPAG"= "#2ECC71",
      "rATG + CsA + EPAG"= "#F39C12",
      "EPAG Monotherapy" = "#9B59B6"
    )
    ggplotly(
      ggplot(scenario_data(), aes(x = time_days, y = .data[[var]], color = Scenario)) +
        geom_line(linewidth = 1.0) +
        scale_color_manual(values = cols) +
        labs(x = "Day", y = ylab) + theme_bw(base_size = 11),
      tooltip = c("x", "y", "colour")
    )
  }
  output$plot_sc_hgb <- renderPlotly(make_scen_ply("Hgb", "Hgb (g/dL)"))
  output$plot_sc_anc <- renderPlotly(make_scen_ply("ANC", "ANC (×10⁹/L)"))
  output$plot_sc_plt <- renderPlotly(make_scen_ply("PLT", "PLT (×10⁹/L)"))
  output$plot_sc_hsc <- renderPlotly(make_scen_ply("HSC", "HSC Pool (% Normal)"))

  output$plot_sc_response <- renderPlotly({
    d180 <- scenario_data() %>%
      filter(time_days >= 179, time_days <= 181) %>%
      group_by(Scenario) %>%
      summarize(
        CR = mean(CR_flag) * 100,
        PR = mean(PR_flag - CR_flag) * 100,
        NR = 100 - CR - PR,
        .groups = "drop"
      ) %>%
      pivot_longer(-Scenario, names_to = "Response", values_to = "Rate") %>%
      mutate(Response = factor(Response, levels = c("NR", "PR", "CR")))

    ggplotly(
      ggplot(d180, aes(x = Scenario, y = pmax(Rate, 0), fill = Response)) +
        geom_col(position = "stack") +
        scale_fill_manual(values = c(CR = "#2ECC71", PR = "#F39C12", NR = "#E74C3C")) +
        labs(x = "", y = "Response Rate (%)") +
        coord_flip() + theme_bw(base_size = 11)
    )
  })

  ## ── Tab 6: Biomarkers & Clones ───────────────────────────
  output$plot_ifng <- renderPlotly({
    make_ply(
      ggplot(res(), aes(x = time_days, y = IFNg_c)) +
        geom_area(fill = "#F39C12", alpha = 0.4) + geom_line(color = "#B7770D") +
        labs(x = "Day", y = "IFN-γ (pg/mL)") + theme_bw()
    )
  })
  output$plot_tnf_il2 <- renderPlotly({
    r <- res() %>% select(time_days, TNFa_c, IL2_c) %>%
      pivot_longer(-time_days, names_to = "Cytokine", values_to = "Conc")
    make_ply(
      ggplot(r, aes(x = time_days, y = Conc, color = Cytokine)) +
        geom_line(linewidth = 1.0) +
        scale_color_manual(values = c(TNFa_c = "#E74C3C", IL2_c = "#9B59B6")) +
        labs(x = "Day", y = "Concentration (pg/mL)") + theme_bw()
    )
  })
  output$plot_pnh <- renderPlotly({
    r <- res() %>% mutate(PNH_pct = PNH_clone * 100)
    make_ply(
      ggplot(r, aes(x = time_days, y = PNH_pct)) +
        geom_area(fill = "#8E44AD", alpha = 0.3) + geom_line(color = "#5B2C6F") +
        labs(x = "Day", y = "PNH Clone (%)") + theme_bw()
    )
  })
  output$tbl_biomarker <- renderTable({
    data.frame(
      Biomarker         = c("IFN-γ", "TNF-α", "IL-2", "IFN-γ/IL-4 ratio", "sFasL", "PNH clone"),
      `Normal Range`    = c("<10 pg/mL", "<15 pg/mL", "<10 pg/mL", "<2", "<300 pg/mL", "<1%"),
      `Active AA Range` = c("30–200 pg/mL", "20–80 pg/mL", "15–60 pg/mL", ">10", "elevated", "2–50%"),
      Significance      = c("Key driver", "Synergistic HSC kill", "T-cell amplifier",
                             "Th1 skew", "Apoptosis signal", "Immune escape marker")
    )
  }, striped = TRUE, bordered = TRUE)

  output$plot_clone_risk <- renderPlotly({
    r <- res() %>% mutate(
      MDS_risk = pmin(PNH_clone * 0.15 + (1 - HSC / 100) * 0.08, 0.25) * 100
    )
    make_ply(
      ggplot(r, aes(x = time_days, y = MDS_risk)) +
        geom_area(fill = "#E74C3C", alpha = 0.4) + geom_line(color = "#922B21") +
        labs(x = "Day", y = "Estimated MDS/Clonal Risk (%)") + theme_bw()
    )
  })
}

## ============================================================
## Run
## ============================================================
shinyApp(ui = ui, server = server)
