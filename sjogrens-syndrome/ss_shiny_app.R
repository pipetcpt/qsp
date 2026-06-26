## ============================================================
## Sjögren's Syndrome QSP – Interactive Shiny Dashboard
## ============================================================
## Tabs:
##   1. Patient Profile & Disease Activity
##   2. Drug PK Profiles
##   3. Disease Biomarkers (IFN-I, BAFF, B cells, Anti-SSA)
##   4. Glandular Function (UWSF, Schirmer, ocular)
##   5. Scenario Comparison (ESSDAI, biomarkers, endpoints)
##   6. Lymphoma Risk & FFS Score
## ============================================================

library(shiny)
library(shinydashboard)
library(mrgsolve)
library(dplyr)
library(tidyr)
library(ggplot2)
library(plotly)
library(DT)
library(shinyWidgets)

# ── mrgsolve model (inline) ──────────────────────────────────────────────────
SS_MODEL_CODE <- '
$PARAM @annotated
TVCL_HCQ: 18.0  : HCQ CL (L/h)
TVV1_HCQ: 800   : HCQ central Vd (L)
TVV2_HCQ: 25000 : HCQ peripheral Vd (L)
TVQ_HCQ : 25    : HCQ Q (L/h)
TVF_HCQ : 0.74  : HCQ bioavailability
TVKA_HCQ: 0.05  : HCQ Ka (1/h)
IC50_HCQ: 200   : HCQ IC50 TLR inh (ng/mL)
EMAX_HCQ: 0.80  : HCQ Emax TLR inh
TVCL_PIL: 7.5   : Pilocarpine CL (L/h)
TVV_PIL : 38    : Pilocarpine Vd (L)
TVF_PIL : 0.80  : Pilocarpine F (oral)
TVKA_PIL: 1.5   : Pilocarpine Ka (1/h)
EC50_PIL: 80    : Pilocarpine EC50 UWSF (ng/mL)
EMAX_PIL: 0.55  : Pilocarpine Emax UWSF
TVCL_RTX: 0.20  : RTX CL (L/h)
TVV1_RTX: 3.6   : RTX central Vd (L)
TVV2_RTX: 4.2   : RTX peripheral Vd (L)
TVQ_RTX : 0.08  : RTX Q (L/h)
EC50_RTX: 50    : RTX EC50 B-depl (mcg/mL)
EMAX_RTX: 0.90  : RTX Emax B-depl
HILL_RTX: 2.0   : RTX Hill coef
TVCL_IAN: 0.18  : Ianalumab CL (L/h)
TVV1_IAN: 3.8   : Ianalumab central Vd (L)
TVV2_IAN: 4.5   : Ianalumab peripheral Vd (L)
TVQ_IAN : 0.09  : Ianalumab Q (L/h)
EC50_IAN: 8.0   : Ianalumab EC50 B-depl (mcg/mL)
EMAX_IAN: 0.88  : Ianalumab Emax B-depl
KIN_IFN : 0.05  : IFN-I Kin
KOUT_IFN: 0.05  : IFN-I Kout
STIM_IFN: 1.5   : IFN-I disease stim factor
KIN_BAF : 0.04  : BAFF Kin
KOUT_BAF: 0.04  : BAFF Kout
STIM_BAF: 1.4   : BAFF STIM by IFN
KIN_BC  : 0.002 : B-cell Kin
KOUT_BC : 0.002 : B-cell Kout
STIM_BC : 1.8   : B-cell STIM by BAFF
KIN_PC  : 0.001 : Plasma cell Kin
KOUT_PC : 0.001 : Plasma cell Kout
KIN_AB  : 0.0003: Anti-SSA Kin
KOUT_AB : 0.0003: Anti-SSA Kout
STIM_AB : 2.0   : Anti-SSA stim by PC
KIN_SAL : 0.01  : Salivary Kin
KOUT_SAL: 0.01  : Salivary Kout
INH_SAL : 0.60  : Salivary inh by IFN-I
BASE_SAL: 0.45  : Salivary baseline fraction
KIN_LAC : 0.008 : Lacrimal Kin
KOUT_LAC: 0.008 : Lacrimal Kout
INH_LAC : 0.55  : Lacrimal inh by IFN-I
BASE_LAC: 0.42  : Lacrimal baseline fraction
KIN_ESD : 0.005 : ESSDAI Kin
KOUT_ESD: 0.005 : ESSDAI Kout
SCAL_ESD: 30    : ESSDAI baseline (median)
LR_RATE : 0.0001: Lymphoma risk rate
TVCL_PRD: 13    : Prednisolone CL (L/h)
TVV_PRD : 50    : Prednisolone Vd (L)
TVKA_PRD: 1.2   : Prednisolone Ka
TVF_PRD : 0.82  : Prednisolone F

$CMT @annotated
DEPOT_HCQ : HCQ depot
C1_HCQ    : HCQ central (mg)
C2_HCQ    : HCQ peripheral (mg)
DEPOT_PIL : Pilocarpine depot
C_PIL     : Pilocarpine central (mg)
C1_RTX    : RTX central (mg)
C2_RTX    : RTX peripheral (mg)
C1_IAN    : Ianalumab central (mg)
C2_IAN    : Ianalumab peripheral (mg)
DEPOT_PRD : Prednisolone depot
C_PRD     : Prednisolone central (mg)
IFN       : IFN-I index
BCELL     : B-cell pool
BAFF_pd   : BAFF serum
AB_SSA    : Anti-SSA titer
SAL       : Salivary function
LAC       : Lacrimal function
ESSDAI_pd : ESSDAI proxy
PLASMA_C  : Plasma cells
LYMPHOMA  : Lymphoma risk

$MAIN
double CL_HCQ = TVCL_HCQ; double V1_HCQ = TVV1_HCQ;
double V2_HCQ = TVV2_HCQ; double Q_HCQ  = TVQ_HCQ;
double KA_HCQ = TVKA_HCQ; double F1_HCQ = TVF_HCQ;
double CL_PIL = TVCL_PIL; double V_PIL  = TVV_PIL;
double KA_PIL = TVKA_PIL; double F1_PIL = TVF_PIL;
double CL_RTX = TVCL_RTX; double V1_RTX = TVV1_RTX;
double V2_RTX = TVV2_RTX; double Q_RTX  = TVQ_RTX;
double CL_IAN = TVCL_IAN; double V1_IAN = TVV1_IAN;
double V2_IAN = TVV2_IAN; double Q_IAN  = TVQ_IAN;
double CL_PRD = TVCL_PRD; double V_PRD  = TVV_PRD;
double KA_PRD = TVKA_PRD; double F1_PRD = TVF_PRD;
F_DEPOT_HCQ = F1_HCQ;
F_DEPOT_PIL = F1_PIL;
F_DEPOT_PRD = F1_PRD;

$ODE
double cp_HCQ = C1_HCQ/TVV1_HCQ;
double cp_PIL = C_PIL/TVV_PIL;
double cp_RTX = C1_RTX/TVV1_RTX;
double cp_IAN = C1_IAN/TVV1_IAN;
double cp_PRD = C_PRD/TVV_PRD;
double HCQ_inh = (EMAX_HCQ*cp_HCQ)/(IC50_HCQ+cp_HCQ);
double PIL_stim = (EMAX_PIL*cp_PIL)/(EC50_PIL+cp_PIL);
double RTX_dep = (EMAX_RTX*pow(cp_RTX,HILL_RTX))/(pow(EC50_RTX,HILL_RTX)+pow(cp_RTX,HILL_RTX));
double IAN_dep = (EMAX_IAN*cp_IAN)/(EC50_IAN+cp_IAN);
double Bcell_supp = 1.0 - fmax(RTX_dep, IAN_dep);
double PRD_inh = 0.5*(cp_PRD/(300.0+cp_PRD));
double IFN_eff  = fmax(IFN,0.001);
double BAFF_eff = fmax(BAFF_pd,0.001);
double BC_eff   = fmax(BCELL,0.001);
double PC_eff   = fmax(PLASMA_C,0.001);

dxdt_DEPOT_HCQ = -TVKA_HCQ*DEPOT_HCQ;
dxdt_C1_HCQ    = TVKA_HCQ*DEPOT_HCQ - (TVCL_HCQ/TVV1_HCQ)*C1_HCQ
                 - (TVQ_HCQ/TVV1_HCQ)*C1_HCQ + (TVQ_HCQ/TVV2_HCQ)*C2_HCQ;
dxdt_C2_HCQ    = (TVQ_HCQ/TVV1_HCQ)*C1_HCQ - (TVQ_HCQ/TVV2_HCQ)*C2_HCQ;

dxdt_DEPOT_PIL = -TVKA_PIL*DEPOT_PIL;
dxdt_C_PIL     = TVKA_PIL*DEPOT_PIL - (TVCL_PIL/TVV_PIL)*C_PIL;

dxdt_C1_RTX    = -(TVCL_RTX/TVV1_RTX)*C1_RTX - (TVQ_RTX/TVV1_RTX)*C1_RTX
                 + (TVQ_RTX/TVV2_RTX)*C2_RTX;
dxdt_C2_RTX    = (TVQ_RTX/TVV1_RTX)*C1_RTX - (TVQ_RTX/TVV2_RTX)*C2_RTX;

dxdt_C1_IAN    = -(TVCL_IAN/TVV1_IAN)*C1_IAN - (TVQ_IAN/TVV1_IAN)*C1_IAN
                 + (TVQ_IAN/TVV2_IAN)*C2_IAN;
dxdt_C2_IAN    = (TVQ_IAN/TVV1_IAN)*C1_IAN - (TVQ_IAN/TVV2_IAN)*C2_IAN;

dxdt_DEPOT_PRD = -TVKA_PRD*DEPOT_PRD;
dxdt_C_PRD     = TVKA_PRD*DEPOT_PRD - (TVCL_PRD/TVV_PRD)*C_PRD;

double IFN_kin  = KIN_IFN*STIM_IFN*(1-HCQ_inh)*(1-PRD_inh);
dxdt_IFN       = IFN_kin - KOUT_IFN*IFN_eff;
double BAFF_kin = KIN_BAF*(1+(STIM_BAF-1)*IFN_eff)*(1-PRD_inh);
dxdt_BAFF_pd   = BAFF_kin - KOUT_BAF*BAFF_eff;
double BC_kin   = KIN_BC*(1+(STIM_BC-1)*BAFF_eff)*Bcell_supp;
dxdt_BCELL     = BC_kin - KOUT_BC*BC_eff;
double PC_kin   = KIN_PC*BC_eff*(1-0.3*RTX_dep);
dxdt_PLASMA_C  = PC_kin - KOUT_PC*PC_eff;
double AB_kin   = KIN_AB*(1+(STIM_AB-1)*PC_eff);
dxdt_AB_SSA    = AB_kin - KOUT_AB*fmax(AB_SSA,0.001);
double SAL_tgt  = BASE_SAL/(1+INH_SAL*(IFN_eff-1));
double SAL_kin  = KIN_SAL*SAL_tgt*(1+PIL_stim);
dxdt_SAL       = SAL_kin - KOUT_SAL*SAL;
double LAC_tgt  = BASE_LAC/(1+INH_LAC*(IFN_eff-1));
double LAC_kin  = KIN_LAC*LAC_tgt*(1+0.4*PIL_stim);
dxdt_LAC       = LAC_kin - KOUT_LAC*LAC;
double ESD_drv  = 0.4*BC_eff + 0.3*IFN_eff + 0.3*fmax(AB_SSA,0.001);
double ESD_tgt  = (SCAL_ESD/123.0)*ESD_drv;
dxdt_ESSDAI_pd = KIN_ESD*ESD_tgt - KOUT_ESD*ESSDAI_pd;
dxdt_LYMPHOMA  = LR_RATE*BC_eff*BAFF_eff;

$CAPTURE
cp_HCQ cp_PIL cp_RTX cp_IAN cp_PRD
IFN BCELL BAFF_pd AB_SSA SAL LAC ESSDAI_pd PLASMA_C LYMPHOMA
RTX_dep IAN_dep HCQ_inh PIL_stim

$INIT
DEPOT_HCQ=0 C1_HCQ=0 C2_HCQ=0
DEPOT_PIL=0 C_PIL=0
C1_RTX=0   C2_RTX=0
C1_IAN=0   C2_IAN=0
DEPOT_PRD=0 C_PRD=0
IFN=1.5 BCELL=1.8 BAFF_pd=1.6 AB_SSA=2.5
SAL=0.42 LAC=0.40 ESSDAI_pd=0.24 PLASMA_C=1.6 LYMPHOMA=0
'

ss_mod <<- suppressMessages(mcode("ss_shiny", SS_MODEL_CODE))

# ── Simulation function ───────────────────────────────────────────────────────
run_sim <- function(dur_weeks = 52,
                    use_hcq   = TRUE,  hcq_dose   = 400,
                    use_pilo  = FALSE, pilo_dose  = 5,
                    use_rtx   = FALSE, rtx_scheme = "2x1g",
                    use_ian   = FALSE, ian_dose   = 300,
                    use_prd   = FALSE, prd_dose   = 10) {
  sim_h <- dur_weeks * 7 * 24
  ev_list <- list()

  if (use_hcq && hcq_dose > 0) {
    t_hcq <- seq(0, sim_h - 1, by = 12)
    ev_list$hcq <- ev(time = t_hcq, amt = hcq_dose / 2,
                      cmt = "DEPOT_HCQ", evid = 1)
  }
  if (use_pilo && pilo_dose > 0) {
    t_pil <- seq(0, sim_h - 1, by = 6)
    ev_list$pilo <- ev(time = t_pil, amt = pilo_dose,
                       cmt = "DEPOT_PIL", evid = 1)
  }
  if (use_rtx) {
    rtx_t <- if (rtx_scheme == "2x1g") c(0, 2*7*24) else
              if (rtx_scheme == "4x375") c(0, 7*24, 2*7*24, 3*7*24) else
              seq(0, sim_h - 1, by = 24 * 7 * 6)
    ev_list$rtx <- ev(time = rtx_t, amt = 1000,
                      cmt = "C1_RTX", evid = 1)
  }
  if (use_ian && ian_dose > 0) {
    t_ian <- seq(0, sim_h - 1, by = 4 * 7 * 24)
    ev_list$ian <- ev(time = t_ian, amt = ian_dose,
                      cmt = "C1_IAN", evid = 1)
  }
  if (use_prd && prd_dose > 0) {
    t_prd <- seq(0, sim_h - 1, by = 24)
    ev_list$prd <- ev(time = t_prd, amt = prd_dose,
                      cmt = "DEPOT_PRD", evid = 1)
  }

  if (length(ev_list) == 0) {
    dose_ev <- ev(time = 0, amt = 0, cmt = 1, evid = 0)
  } else {
    dose_ev <- do.call(c, ev_list)
  }

  ss_mod %>%
    ev(dose_ev) %>%
    mrgsim(end = sim_h, delta = 24) %>%
    as.data.frame() %>%
    mutate(
      time_weeks    = time / (7 * 24),
      ESSDAI_abs    = ESSDAI_pd * 123,
      UWSF_mL15     = SAL * 1.5,
      Schirmer_mm   = LAC * 15,
      AntiSSA_EU    = AB_SSA * 250,
      BAFF_pgmL     = BAFF_pd * 1800,
      IFN_score     = IFN,
      B_pct         = BCELL * 100,
      ESSPRI_dry    = pmax(0, 10 * (1 - SAL / 0.8)),
      ESSPRI_fat    = pmax(0, 10 * (IFN / 2.5)),
      ESSPRI_pain   = pmax(0, 10 * (AB_SSA / 4.0)),
      ESSPRI_total  = (ESSPRI_dry + ESSPRI_fat + ESSPRI_pain) / 3
    )
}

# ── Pre-compute all 5 scenarios ───────────────────────────────────────────────
precompute_scenarios <- function() {
  list(
    "No Treatment"              = run_sim(),
    "HCQ 400 mg/d"              = run_sim(use_hcq = TRUE, hcq_dose = 400),
    "HCQ + Pilocarpine 5mg QID" = run_sim(use_hcq = TRUE, use_pilo = TRUE),
    "HCQ + Rituximab 2×1g"      = run_sim(use_hcq = TRUE, use_rtx = TRUE),
    "HCQ + Ianalumab 300mg q4w" = run_sim(use_hcq = TRUE, use_ian = TRUE)
  )
}

# ─────────────────────────────────────────────────────────────────────────────
# UI
# ─────────────────────────────────────────────────────────────────────────────
ui <- dashboardPage(
  skin = "blue",
  dashboardHeader(
    title = tags$span(
      tags$img(src = "", height = "20px"),
      "Sjögren's Syndrome QSP Dashboard"
    ),
    titleWidth = 320
  ),
  dashboardSidebar(
    width = 270,
    sidebarMenu(
      menuItem("Patient Profile",    tabName = "profile",   icon = icon("user-md")),
      menuItem("Drug PK Profiles",   tabName = "pk",        icon = icon("pills")),
      menuItem("Disease Biomarkers", tabName = "biomarkers",icon = icon("dna")),
      menuItem("Glandular Function", tabName = "gland",     icon = icon("tint")),
      menuItem("Scenario Comparison",tabName = "scenario",  icon = icon("chart-bar")),
      menuItem("Lymphoma Risk",      tabName = "lymphoma",  icon = icon("exclamation-triangle"))
    ),
    hr(),
    h5("Patient Parameters", style = "padding-left:10px; color:#ccc"),
    sliderInput("age",   "Age (years)",    30, 80, 52, 1),
    sliderInput("wt",    "Weight (kg)",    40, 120, 65, 1),
    selectInput("sex", "Sex",
                choices = c("Female" = "F", "Male" = "M"), selected = "F"),
    sliderInput("dis_dur", "Disease Duration (years)", 0, 20, 4, 0.5),
    hr(),
    h5("Treatment Selection", style = "padding-left:10px; color:#ccc"),
    checkboxInput("use_hcq",  "HCQ (400 mg/d)",          TRUE),
    checkboxInput("use_pilo", "Pilocarpine (5 mg QID)",   FALSE),
    checkboxInput("use_rtx",  "Rituximab (2×1g)",         FALSE),
    checkboxInput("use_ian",  "Ianalumab (300 mg q4w)",   FALSE),
    checkboxInput("use_prd",  "Prednisolone",              FALSE),
    conditionalPanel("input.use_prd",
      sliderInput("prd_dose", "Prednisolone dose (mg/d)", 2.5, 30, 10, 2.5)
    ),
    sliderInput("dur_wks", "Simulation Duration (weeks)", 12, 104, 52, 4),
    actionButton("run_btn", "Run Simulation", class = "btn-primary",
                 style = "width:90%; margin:5px 5% 5px 5%;")
  ),
  dashboardBody(
    tags$head(tags$style(HTML(
      ".content-wrapper { background-color: #f4f6f9; }
       .info-box { min-height: 75px; }
       .info-box-icon { height: 75px; line-height: 75px; }
       .info-box-content { padding-top: 8px; }"
    ))),
    tabItems(
      # ── Tab 1: Patient Profile ──────────────────────────────────────────────
      tabItem(
        tabName = "profile",
        fluidRow(
          infoBoxOutput("box_essdai",   width = 3),
          infoBoxOutput("box_esspri",   width = 3),
          infoBoxOutput("box_uwsf",     width = 3),
          infoBoxOutput("box_schirmer", width = 3)
        ),
        fluidRow(
          box(title = "ESSDAI Over Time", width = 6, status = "primary",
              plotlyOutput("plot_essdai_profile", height = "300px")),
          box(title = "ESSPRI Over Time", width = 6, status = "info",
              plotlyOutput("plot_esspri_profile", height = "300px"))
        ),
        fluidRow(
          box(title = "Classification Criteria Summary", width = 12, status = "warning",
              collapsible = TRUE,
              DTOutput("table_criteria"))
        )
      ),
      # ── Tab 2: Drug PK ─────────────────────────────────────────────────────
      tabItem(
        tabName = "pk",
        fluidRow(
          box(title = "HCQ Plasma Concentration", width = 6, status = "primary",
              sliderInput("hcq_dose_pk", "HCQ dose (mg/day)", 200, 600, 400, 100),
              plotlyOutput("plot_hcq_pk", height = "280px")),
          box(title = "Pilocarpine PK", width = 6, status = "success",
              sliderInput("pilo_dose_pk", "Pilocarpine dose (mg)", 2.5, 10, 5, 2.5),
              plotlyOutput("plot_pilo_pk", height = "280px"))
        ),
        fluidRow(
          box(title = "Rituximab PK (2×1g, weeks 0 & 2)", width = 6, status = "warning",
              plotlyOutput("plot_rtx_pk", height = "280px")),
          box(title = "Ianalumab PK (300 mg SC q4w)", width = 6, status = "danger",
              sliderInput("ian_dose_pk", "Ianalumab dose (mg)", 50, 600, 300, 50),
              plotlyOutput("plot_ian_pk", height = "280px"))
        )
      ),
      # ── Tab 3: Biomarkers ──────────────────────────────────────────────────
      tabItem(
        tabName = "biomarkers",
        fluidRow(
          box(title = "Type I IFN Index", width = 6, status = "danger",
              plotlyOutput("plot_ifn", height = "260px")),
          box(title = "BAFF Serum Level (pg/mL)", width = 6, status = "warning",
              plotlyOutput("plot_baff", height = "260px"))
        ),
        fluidRow(
          box(title = "B-cell Pool (% of Normal)", width = 6, status = "primary",
              plotlyOutput("plot_bcell", height = "260px")),
          box(title = "Anti-SSA/Ro Titer (EU)", width = 6, status = "info",
              plotlyOutput("plot_antiSSA", height = "260px"))
        ),
        fluidRow(
          box(title = "Biomarker Summary Table (Wk 24 & 52)", width = 12, status = "success",
              DTOutput("table_biomarkers"))
        )
      ),
      # ── Tab 4: Glandular Function ──────────────────────────────────────────
      tabItem(
        tabName = "gland",
        fluidRow(
          box(title = "Salivary Function – UWSF (mL/15 min)", width = 6, status = "primary",
              plotlyOutput("plot_uwsf", height = "280px")),
          box(title = "Lacrimal Function – Schirmer (mm/5 min)", width = 6, status = "info",
              plotlyOutput("plot_schirmer", height = "280px"))
        ),
        fluidRow(
          box(title = "Pilocarpine Effect on UWSF", width = 6, status = "success",
              plotlyOutput("plot_pilo_uwsf", height = "280px")),
          box(title = "Glandular Function Summary", width = 6, status = "warning",
              tableOutput("table_gland"))
        )
      ),
      # ── Tab 5: Scenario Comparison ─────────────────────────────────────────
      tabItem(
        tabName = "scenario",
        fluidRow(
          box(title = "ESSDAI – All Scenarios", width = 6, status = "primary",
              plotlyOutput("plot_scen_essdai", height = "280px")),
          box(title = "B-cell Pool – All Scenarios", width = 6, status = "warning",
              plotlyOutput("plot_scen_bcell", height = "280px"))
        ),
        fluidRow(
          box(title = "Anti-SSA Titer – All Scenarios", width = 6, status = "info",
              plotlyOutput("plot_scen_ssa", height = "280px")),
          box(title = "UWSF – All Scenarios", width = 6, status = "success",
              plotlyOutput("plot_scen_uwsf", height = "280px"))
        ),
        fluidRow(
          box(title = "Endpoint Comparison Table (Week 24)", width = 12, status = "danger",
              DTOutput("table_scenarios"))
        )
      ),
      # ── Tab 6: Lymphoma Risk ───────────────────────────────────────────────
      tabItem(
        tabName = "lymphoma",
        fluidRow(
          valueBoxOutput("vbox_ffs",    width = 4),
          valueBoxOutput("vbox_malt",   width = 4),
          valueBoxOutput("vbox_dlbcl",  width = 4)
        ),
        fluidRow(
          box(title = "Lymphoma Risk Accumulator Over Time", width = 6, status = "danger",
              plotlyOutput("plot_lymphoma", height = "280px")),
          box(title = "CXCL13 & B-cell Dynamics (Ectopic GC markers)", width = 6,
              status = "warning",
              plotlyOutput("plot_ectgc", height = "280px"))
        ),
        fluidRow(
          box(title = "FFS Score Components & Risk Classification", width = 12,
              status = "primary",
              DTOutput("table_ffs"))
        )
      )
    )
  )
)

# ─────────────────────────────────────────────────────────────────────────────
# Server
# ─────────────────────────────────────────────────────────────────────────────
server <- function(input, output, session) {

  # Reactive simulation (custom treatment)
  sim_data <- eventReactive(input$run_btn, {
    run_sim(
      dur_weeks = input$dur_wks,
      use_hcq   = input$use_hcq,   hcq_dose  = 400,
      use_pilo  = input$use_pilo,  pilo_dose = 5,
      use_rtx   = input$use_rtx,
      use_ian   = input$use_ian,   ian_dose  = 300,
      use_prd   = input$use_prd,   prd_dose  = input$prd_dose
    )
  }, ignoreNULL = FALSE)

  # Pre-computed scenario comparison
  all_scen <- reactive({
    withProgress(message = "Running 5 scenarios...", {
      precompute_scenarios()
    })
  })

  combined <- reactive({
    sc <- all_scen()
    bind_rows(lapply(names(sc), function(nm) mutate(sc[[nm]], scenario = nm)))
  })

  # ── Tab 1: Info boxes ──────────────────────────────────────────────────────
  output$box_essdai <- renderInfoBox({
    d <- sim_data()
    wk24 <- d %>% filter(abs(time_weeks - 24) < 0.5)
    val  <- round(mean(wk24$ESSDAI_abs), 1)
    col  <- if (val < 13) "green" else if (val < 28) "yellow" else "red"
    infoBox("ESSDAI at Wk 24", paste0(val, " / 123"),
            subtitle = if (val < 13) "Low activity" else if (val < 28) "Moderate" else "High",
            icon = icon("stethoscope"), color = col)
  })

  output$box_esspri <- renderInfoBox({
    d <- sim_data()
    wk24 <- d %>% filter(abs(time_weeks - 24) < 0.5)
    val  <- round(mean(wk24$ESSPRI_total), 1)
    infoBox("ESSPRI at Wk 24", paste0(val, " / 10"),
            subtitle = if (val < 5) "Acceptable" else "Unacceptable",
            icon = icon("heart"), color = if (val < 5) "green" else "red")
  })

  output$box_uwsf <- renderInfoBox({
    d <- sim_data()
    wk24 <- d %>% filter(abs(time_weeks - 24) < 0.5)
    val  <- round(mean(wk24$UWSF_mL15), 2)
    infoBox("UWSF at Wk 24", paste0(val, " mL/15 min"),
            subtitle = if (val >= 1.5) "Normal" else if (val >= 0.5) "Reduced" else "Severe xerostomia",
            icon = icon("tint"), color = if (val >= 1.5) "green" else if (val >= 0.5) "yellow" else "red")
  })

  output$box_schirmer <- renderInfoBox({
    d <- sim_data()
    wk24 <- d %>% filter(abs(time_weeks - 24) < 0.5)
    val  <- round(mean(wk24$Schirmer_mm), 1)
    infoBox("Schirmer at Wk 24", paste0(val, " mm/5 min"),
            subtitle = if (val > 5) "Negative" else "Positive (≤5 mm)",
            icon = icon("eye"), color = if (val > 5) "green" else "red")
  })

  output$plot_essdai_profile <- renderPlotly({
    d <- sim_data()
    p <- ggplot(d, aes(time_weeks, ESSDAI_abs)) +
      geom_line(colour = "#1976D2", linewidth = 1) +
      geom_hline(yintercept = 27, linetype = "dashed", colour = "orange") +
      geom_hline(yintercept = 13, linetype = "dashed", colour = "green") +
      annotate("text", x = max(d$time_weeks)*0.02, y = 28.5,
               label = "ΔESSDAI ≥3 threshold", hjust = 0, size = 3) +
      labs(x = "Weeks", y = "ESSDAI") + theme_bw(base_size = 11)
    ggplotly(p)
  })

  output$plot_esspri_profile <- renderPlotly({
    d <- sim_data()
    p <- ggplot(d, aes(time_weeks)) +
      geom_line(aes(y = ESSPRI_total, colour = "Total"), linewidth = 1) +
      geom_line(aes(y = ESSPRI_dry,   colour = "Dryness"), linetype = "dashed") +
      geom_line(aes(y = ESSPRI_fat,   colour = "Fatigue"),  linetype = "dashed") +
      geom_line(aes(y = ESSPRI_pain,  colour = "Pain"),     linetype = "dashed") +
      geom_hline(yintercept = 5, linetype = "dotted", colour = "red") +
      scale_colour_manual(values = c(Total="#D32F2F", Dryness="#1976D2",
                                     Fatigue="#388E3C", Pain="#F57C00")) +
      labs(x = "Weeks", y = "ESSPRI (0-10)", colour = "Domain") +
      theme_bw(base_size = 11)
    ggplotly(p)
  })

  output$table_criteria <- renderDT({
    data.frame(
      Criterion = c("Labial gland biopsy (focus score ≥1)",
                    "Anti-SSA/Ro positivity",
                    "Ocular staining score ≥5",
                    "Schirmer's test ≤5 mm/5 min",
                    "UWSF ≤0.1 mL/min",
                    "Ocular symptoms (OSDI)"),
      Weight = c(3, 3, 1, 1, 1, 1),
      `2016 ACR/EULAR` = c("Yes","Yes","Yes","Yes","Yes","Yes")
    )
  }, options = list(pageLength = 6, dom = 't'))

  # ── Tab 2: PK ──────────────────────────────────────────────────────────────
  output$plot_hcq_pk <- renderPlotly({
    d_hcq <- run_sim(use_hcq = TRUE, hcq_dose = input$hcq_dose_pk, dur_weeks = 12)
    p <- ggplot(d_hcq, aes(time_weeks, cp_HCQ)) +
      geom_line(colour = "#1976D2", linewidth = 1) +
      geom_hline(yintercept = 200, linetype = "dashed", colour = "red") +
      annotate("text", x = 1, y = 210, label = "IC50 = 200 ng/mL", hjust = 0, size = 3) +
      labs(x = "Weeks", y = "HCQ (ng/mL)") + theme_bw(base_size = 10)
    ggplotly(p)
  })

  output$plot_pilo_pk <- renderPlotly({
    d_pil <- run_sim(use_pilo = TRUE, pilo_dose = input$pilo_dose_pk, dur_weeks = 2)
    p <- ggplot(d_pil, aes(time_weeks * 7 * 24, cp_PIL)) +
      geom_line(colour = "#388E3C", linewidth = 1) +
      labs(x = "Hours", y = "Pilocarpine (ng/mL)") + theme_bw(base_size = 10)
    ggplotly(p)
  })

  output$plot_rtx_pk <- renderPlotly({
    d_rtx <- run_sim(use_rtx = TRUE, dur_weeks = 26)
    p <- ggplot(d_rtx, aes(time_weeks, cp_RTX)) +
      geom_line(colour = "#7B1FA2", linewidth = 1) +
      geom_hline(yintercept = 50, linetype = "dashed", colour = "red") +
      labs(x = "Weeks", y = "Rituximab (mcg/mL)") + theme_bw(base_size = 10)
    ggplotly(p)
  })

  output$plot_ian_pk <- renderPlotly({
    ian_t <- seq(0, 51 * 7 * 24, by = 4 * 7 * 24)
    d_ian <- run_sim(use_ian = TRUE, ian_dose = input$ian_dose_pk, dur_weeks = 52)
    p <- ggplot(d_ian, aes(time_weeks, cp_IAN)) +
      geom_line(colour = "#F57C00", linewidth = 1) +
      geom_hline(yintercept = 8, linetype = "dashed", colour = "red") +
      labs(x = "Weeks", y = "Ianalumab (mcg/mL)") + theme_bw(base_size = 10)
    ggplotly(p)
  })

  # ── Tab 3: Biomarkers ──────────────────────────────────────────────────────
  mk_biom_plot <- function(df, yvar, ytitle, col, hline = NULL) {
    p <- ggplot(df, aes(time_weeks, .data[[yvar]])) +
      geom_line(colour = col, linewidth = 1) +
      labs(x = "Weeks", y = ytitle) +
      theme_bw(base_size = 10)
    if (!is.null(hline)) p <- p + geom_hline(yintercept = hline, linetype = "dashed")
    ggplotly(p)
  }

  output$plot_ifn      <- renderPlotly(mk_biom_plot(sim_data(), "IFN_score",  "IFN-I index", "#D32F2F", 1.0))
  output$plot_baff     <- renderPlotly(mk_biom_plot(sim_data(), "BAFF_pgmL",  "BAFF (pg/mL)", "#F57C00", 1800))
  output$plot_bcell    <- renderPlotly(mk_biom_plot(sim_data(), "B_pct",      "B cells (%)", "#1976D2", 100))
  output$plot_antiSSA  <- renderPlotly(mk_biom_plot(sim_data(), "AntiSSA_EU", "Anti-SSA (EU)", "#7B1FA2", 250))

  output$table_biomarkers <- renderDT({
    d <- sim_data()
    bind_rows(
      d %>% filter(abs(time_weeks - 24) < 0.5) %>% slice(1) %>%
        mutate(Timepoint = "Week 24"),
      d %>% filter(abs(time_weeks - 52) < 0.5) %>% slice(1) %>%
        mutate(Timepoint = "Week 52")
    ) %>%
      select(Timepoint, IFN_score, BAFF_pgmL, B_pct, AntiSSA_EU,
             ESSDAI_abs, UWSF_mL15, Schirmer_mm) %>%
      mutate(across(where(is.numeric), ~round(., 2)))
  }, options = list(dom = 't', pageLength = 2))

  # ── Tab 4: Glandular Function ──────────────────────────────────────────────
  output$plot_uwsf <- renderPlotly({
    d <- sim_data()
    p <- ggplot(d, aes(time_weeks, UWSF_mL15)) +
      geom_line(colour = "#1976D2", linewidth = 1) +
      geom_hline(yintercept = 1.5, linetype = "dashed", colour = "green") +
      geom_hline(yintercept = 0.1, linetype = "dotted", colour = "red") +
      labs(x = "Weeks", y = "UWSF (mL/15 min)") + theme_bw(base_size = 10)
    ggplotly(p)
  })

  output$plot_schirmer <- renderPlotly({
    d <- sim_data()
    p <- ggplot(d, aes(time_weeks, Schirmer_mm)) +
      geom_line(colour = "#006064", linewidth = 1) +
      geom_hline(yintercept = 5, linetype = "dashed", colour = "red") +
      annotate("text", x = 1, y = 5.5, label = "≤5 mm = KCS positive",
               hjust = 0, size = 3) +
      labs(x = "Weeks", y = "Schirmer (mm/5 min)") + theme_bw(base_size = 10)
    ggplotly(p)
  })

  output$plot_pilo_uwsf <- renderPlotly({
    d_no  <- run_sim(use_pilo = FALSE, dur_weeks = 12) %>% mutate(pilo = "Off")
    d_yes <- run_sim(use_pilo = TRUE,  dur_weeks = 12) %>% mutate(pilo = "On (5 mg QID)")
    d_comb <- bind_rows(d_no, d_yes)
    p <- ggplot(d_comb, aes(time_weeks, UWSF_mL15, colour = pilo)) +
      geom_line(linewidth = 1) +
      scale_colour_manual(values = c("Off" = "#D32F2F", "On (5 mg QID)" = "#388E3C")) +
      labs(x = "Weeks", y = "UWSF (mL/15 min)", colour = "Pilocarpine") +
      theme_bw(base_size = 10)
    ggplotly(p)
  })

  output$table_gland <- renderTable({
    d <- sim_data()
    wk24 <- d %>% filter(abs(time_weeks - 24) < 0.5) %>% slice(1)
    wk52 <- d %>% filter(abs(time_weeks - 52) < 0.5) %>% slice(1)
    data.frame(
      Endpoint       = c("UWSF (mL/15 min)", "Schirmer (mm/5 min)",
                          "UWSF norm threshold", "Schirmer pos. threshold"),
      `Week 24`      = c(round(wk24$UWSF_mL15, 2), round(wk24$Schirmer_mm, 1), "≥1.5", ">5"),
      `Week 52`      = c(round(wk52$UWSF_mL15, 2), round(wk52$Schirmer_mm, 1), "≥1.5", ">5"),
      Status_Wk24    = c(
        if (wk24$UWSF_mL15 >= 1.5) "Normal" else "Reduced",
        if (wk24$Schirmer_mm > 5)  "Negative" else "Positive (KCS)",
        "", ""
      )
    )
  })

  # ── Tab 5: Scenario Comparison ─────────────────────────────────────────────
  scen_colors <- c(
    "No Treatment"              = "#D32F2F",
    "HCQ 400 mg/d"              = "#1976D2",
    "HCQ + Pilocarpine 5mg QID" = "#388E3C",
    "HCQ + Rituximab 2×1g"      = "#7B1FA2",
    "HCQ + Ianalumab 300mg q4w" = "#F57C00"
  )

  mk_scen_plot <- function(df, yvar, ytitle) {
    p <- ggplot(df, aes(time_weeks, .data[[yvar]], colour = scenario)) +
      geom_line(linewidth = 0.9) +
      scale_colour_manual(values = scen_colors) +
      labs(x = "Weeks", y = ytitle, colour = "") +
      theme_bw(base_size = 10) +
      theme(legend.position = "bottom", legend.text = element_text(size = 8))
    ggplotly(p) %>% layout(legend = list(orientation = "h", y = -0.2))
  }

  output$plot_scen_essdai <- renderPlotly(mk_scen_plot(combined(), "ESSDAI_abs",  "ESSDAI (0-123)"))
  output$plot_scen_bcell  <- renderPlotly(mk_scen_plot(combined(), "B_pct",       "B cells (%)"))
  output$plot_scen_ssa    <- renderPlotly(mk_scen_plot(combined(), "AntiSSA_EU",  "Anti-SSA (EU)"))
  output$plot_scen_uwsf   <- renderPlotly(mk_scen_plot(combined(), "UWSF_mL15",  "UWSF (mL/15 min)"))

  output$table_scenarios <- renderDT({
    combined() %>%
      filter(abs(time_weeks - 24) < 0.5) %>%
      group_by(scenario) %>%
      summarise(
        ESSDAI          = round(mean(ESSDAI_abs), 1),
        `ΔESSDAI`       = round(mean(ESSDAI_abs) - 30, 1),
        `Response (%)`  = as.integer(mean(ESSDAI_abs) < 27) * 100,
        `UWSF (mL/15m)` = round(mean(UWSF_mL15),  2),
        `Schirmer (mm)` = round(mean(Schirmer_mm), 1),
        `Anti-SSA (EU)` = round(mean(AntiSSA_EU),  0),
        `B-cells (%)`   = round(mean(B_pct),       1),
        .groups = "drop"
      )
  }, options = list(dom = 't', pageLength = 5))

  # ── Tab 6: Lymphoma Risk ───────────────────────────────────────────────────
  ffs_vals <- reactive({
    d <- sim_data()
    wk52 <- d %>% filter(abs(time_weeks - 52) < 0.5) %>% slice(1)
    list(
      ffs = sum(c(
        wk52$BAFF_pgmL > 2200,       # C4 surrogate (low C4 → +1)
        wk52$B_pct > 150,             # β2-MG proxy
        wk52$AntiSSA_EU > 400,        # cryoglobulin proxy
        wk52$ESSDAI_abs > 40          # parotid enlargement
      )),
      risk = wk52$LYMPHOMA,
      essdai = wk52$ESSDAI_abs
    )
  })

  output$vbox_ffs <- renderValueBox({
    fv <- ffs_vals()
    col <- if (fv$ffs == 0) "green" else if (fv$ffs == 1) "yellow" else "red"
    lbl <- if (fv$ffs == 0) "Low (<1%/yr)" else if (fv$ffs == 1) "Intermediate" else "High (≥7%/yr)"
    valueBox(paste("FFS =", fv$ffs), lbl, icon = icon("shield-alt"), color = col)
  })

  output$vbox_malt <- renderValueBox({
    fv <- ffs_vals()
    ann_risk <- fv$risk * 8760  # hours to years
    valueBox(
      sprintf("%.4f%%", ann_risk * 100),
      "Modelled MALT lymphoma rate/yr",
      icon = icon("exclamation-triangle"),
      color = if (ann_risk < 0.01) "green" else if (ann_risk < 0.05) "yellow" else "red"
    )
  })

  output$vbox_dlbcl <- renderValueBox({
    fv <- ffs_vals()
    valueBox(
      if (fv$ffs < 2) "Low" else "Monitor",
      "DLBCL transformation risk",
      icon = icon("radiation-alt"), color = if (fv$ffs < 2) "green" else "red"
    )
  })

  output$plot_lymphoma <- renderPlotly({
    d <- sim_data()
    p <- ggplot(d, aes(time_weeks, LYMPHOMA)) +
      geom_area(fill = "#EF9A9A", alpha = 0.4) +
      geom_line(colour = "#D32F2F", linewidth = 1) +
      labs(x = "Weeks", y = "Cumulative lymphoma risk index") +
      theme_bw(base_size = 10)
    ggplotly(p)
  })

  output$plot_ectgc <- renderPlotly({
    d <- sim_data()
    p <- ggplot(d, aes(time_weeks)) +
      geom_line(aes(y = B_pct,       colour = "B-cell pool (%)"), linewidth = 1) +
      geom_line(aes(y = BAFF_pgmL/18, colour = "BAFF×0.056 (scaled)"), linewidth = 1) +
      scale_colour_manual(values = c("B-cell pool (%)" = "#1976D2",
                                     "BAFF×0.056 (scaled)" = "#F57C00")) +
      labs(x = "Weeks", y = "Scaled index", colour = "") +
      theme_bw(base_size = 10)
    ggplotly(p)
  })

  output$table_ffs <- renderDT({
    d <- sim_data()
    wk52 <- d %>% filter(abs(time_weeks - 52) < 0.5) %>% slice(1)
    data.frame(
      `FFS Component`   = c("Low C4 (cryoglobulin-related)",
                             "β₂-MG > 3 mg/L",
                             "Cryoglobulins present",
                             "Parotid enlargement"),
      `Model Surrogate` = c("BAFF > 2200 pg/mL",
                             "B-cell pool > 150%",
                             "Anti-SSA > 400 EU",
                             "ESSDAI > 40"),
      `Present?`        = c(
        if (wk52$BAFF_pgmL > 2200) "Yes (+1)" else "No",
        if (wk52$B_pct > 150)      "Yes (+1)" else "No",
        if (wk52$AntiSSA_EU > 400)  "Yes (+1)" else "No",
        if (wk52$ESSDAI_abs > 40)  "Yes (+1)" else "No"
      )
    )
  }, options = list(dom = 't', pageLength = 4))
}

shinyApp(ui = ui, server = server)
