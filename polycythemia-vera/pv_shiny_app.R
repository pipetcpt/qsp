## ============================================================
## Polycythemia Vera (PV) — Shiny Interactive Dashboard
## ============================================================
## Tabs:
##   1. Patient Profile & Risk Stratification
##   2. Drug Pharmacokinetics
##   3. Hematologic Response (HCT, RBC, PLT, WBC)
##   4. Spleen & Symptom Response
##   5. Disease Progression (Fibrosis, Allele Burden)
##   6. Clinical Endpoints & Scenario Comparison
##   7. Biomarker Panel
##   8. About / Model Description
## ============================================================

library(shiny)
library(shinydashboard)
library(mrgsolve)
library(dplyr)
library(ggplot2)
library(tidyr)
library(DT)
library(plotly)

## ============================================================
## Inline mrgsolve model (same as pv_mrgsolve_model.R but simplified)
## ============================================================
pv_model_code <- '
$PROB PV QSP Shiny Model (22-compartment)

$PARAM @annotated
k_clone_grow:0.012:Mutant HSC growth(/day)
k_clone_death:0.003:Mutant HSC death(/day)
k_wt_grow:0.008:WT HSC growth(/day)
k_wt_death:0.005:WT HSC death(/day)
tau_allele:30:Allele burden time constant(days)
kin_stat5:0.80:STAT5 activation(/day)
kout_stat5:0.75:STAT5 deactivation(/day)
STAT5_base:1.0:Baseline STAT5(normalized)
k_bfu_prod:18.0:Mutant BFU-E production
k_bfu_prod_wt:6.0:WT BFU-E production
k_bfu_mat:0.15:BFU-E maturation(/day)
k_rbc_prod:0.028:RBC production(/day)
k_rbc_dest:0.00833:RBC destruction(/day)
RBC_base:5.5:Baseline RBC(x10^12/L)
HCT_per_RBC:0.082:HCT fraction per RBC unit
phlebotomy_eff:0.045:RBC removed/phlebotomy
k_plt_prod:2.5:PLT production
k_plt_dest:0.10:PLT destruction(/day)
PLT_base:350:Baseline PLT(x10^9/L)
k_wbc_prod:0.30:WBC production
k_wbc_dest:0.14:WBC destruction(/day)
WBC_base:8.5:Baseline WBC(x10^9/L)
k_spl_grow:0.0015:Spleen growth rate
k_spl_base:0.0004:Spleen regulation(/day)
Spleen_norm:400:Normal spleen volume(mL)
k_tss_rise:0.04:TSS increase rate
k_tss_fall:0.035:TSS natural fall(/day)
TSS_base:8.0:Baseline TSS
k_mf:0.00008:MF progression rate
TGFb_driver:1.5:TGF-beta fibrosis driver
MF_max:3.0:Maximum fibrosis grade
k_thromb_base:0.005:Baseline thrombosis hazard(/day)
HCT_thr_threshold:45.0:HCT threshold for thrombosis(%)
k_thromb_HCT:0.0008:HCT excess thrombosis factor
k_thromb_PLT:0.0001:PLT excess thrombosis factor
ka_rux:2.16:Rux absorption(/h)
CL_rux:17.7:Rux clearance(L/h)
Vd1_rux:50.0:Rux central volume(L)
Vd2_rux:22.0:Rux peripheral volume(L)
CLD_rux:12.0:Rux inter-compartmental CL(L/h)
Emax_rux:0.92:Rux max JAK2 inhibition
IC50_rux:450:Rux IC50(ng/mL)
Hill_rux:1.2:Rux Hill coefficient
ka_hu:2.88:HU absorption(/h)
CL_hu:18.5:HU clearance(L/h)
Vd_hu:50.0:HU volume(L)
Emax_hu_rbc:0.65:HU max RBC suppression
Emax_hu_plt:0.70:HU max PLT suppression
Emax_hu_wbc:0.75:HU max WBC suppression
EC50_hu:8.5:HU EC50(mg/L)
Hill_hu:1.1:HU Hill coefficient
ka_ifn:0.15:IFN absorption(/day)
CL_ifn:1.2:IFN clearance(L/day)
Vd_ifn:12.0:IFN volume(L)
Emax_ifn_clone:0.90:IFN max clone suppression
EC50_ifn:250:IFN EC50(IU/mL)
Hill_ifn:1.4:IFN Hill coefficient
ka_asp:4.32:Aspirin absorption(/h)
CL_asp:25.0:Aspirin clearance(L/h)
Vd_asp:8.0:Aspirin volume(L)
Emax_asp_plt:0.85:Aspirin max PLT inhibition
EC50_asp:0.05:Aspirin EC50(mg/L)
ka_fed:0.48:Fedratinib absorption(/h)
CL_fed:4.2:Fedratinib clearance(L/h)
Vd_fed:212:Fedratinib volume(L)
Emax_fed:0.88:Fedratinib max JAK2 inhibition
IC50_fed:600:Fedratinib IC50(ng/mL)
Hill_fed:1.2:Fedratinib Hill coefficient
use_phlebotomy:0:Phlebotomy(0/1)
phlebotomy_freq:0:Phlebotomy sessions/year
use_HU:0:Hydroxyurea(0/1)
HU_daily_dose:0:HU daily dose(mg)
use_RUX:0:Ruxolitinib(0/1)
RUX_BID_dose:0:Ruxolitinib BID dose(mg)
use_IFN:0:Ropeg-IFN(0/1)
use_ASP:0:Aspirin(0/1)
ASP_daily_dose:0:Aspirin dose(mg/day)
use_FED:0:Fedratinib(0/1)
FED_daily_dose:0:Fedratinib daily dose(mg)

$INIT
mut_clone=0.30
wt_clone=0.70
STAT5=1.40
allele_burden=45.0
BFU_E_mut=80.0
BFU_E_wt=40.0
RBC=6.8
HCT=55.0
PLT=550
WBC=12.5
Spleen_vol=850
MPN_SAF=20.0
MF_score=0.30
Thromb_hazard=0.0
RUX_gut=0.0
RUX_cent=0.0
RUX_periph=0.0
HU_cent=0.0
IFN_sc=0.0
IFN_cent=0.0
ASP_cent=0.0
FED_cent=0.0

$MAIN
double RUX_C = RUX_cent;
double rux_inhib = use_RUX*Emax_rux*pow(RUX_C,Hill_rux)/(pow(IC50_rux,Hill_rux)+pow(RUX_C,Hill_rux));
double FED_C = FED_cent;
double fed_inhib = use_FED*Emax_fed*pow(FED_C,Hill_fed)/(pow(IC50_fed,Hill_fed)+pow(FED_C,Hill_fed));
double jak2_inhib = 1.0-(1.0-rux_inhib)*(1.0-fed_inhib);
double HU_C = HU_cent;
double hu_rbc = use_HU*Emax_hu_rbc*pow(HU_C,Hill_hu)/(pow(EC50_hu,Hill_hu)+pow(HU_C,Hill_hu));
double hu_plt = use_HU*Emax_hu_plt*pow(HU_C,Hill_hu)/(pow(EC50_hu,Hill_hu)+pow(HU_C,Hill_hu));
double hu_wbc = use_HU*Emax_hu_wbc*pow(HU_C,Hill_hu)/(pow(EC50_hu,Hill_hu)+pow(HU_C,Hill_hu));
double IFN_C = IFN_cent;
double ifn_eff = use_IFN*Emax_ifn_clone*pow(IFN_C,Hill_ifn)/(pow(EC50_ifn,Hill_ifn)+pow(IFN_C,Hill_ifn));
double asp_inh = use_ASP*Emax_asp_plt*ASP_cent/(EC50_asp+ASP_cent);
double phlebotomy_drain = use_phlebotomy*phlebotomy_freq*phlebotomy_eff/365.0;
double STAT5_eff = STAT5*(1.0-jak2_inhib);
double HCT_excess = (HCT>HCT_thr_threshold)?(HCT-HCT_thr_threshold):0.0;
double PLT_excess = (PLT>400.0)?(PLT-400.0)/100.0:0.0;

$ODE
dxdt_mut_clone = k_clone_grow*mut_clone*(1.0-mut_clone-wt_clone)-k_clone_death*mut_clone-ifn_eff*mut_clone;
dxdt_wt_clone = k_wt_grow*wt_clone*(1.0-mut_clone-wt_clone)-k_wt_death*wt_clone;
dxdt_allele_burden = (mut_clone/(mut_clone+wt_clone+1e-6)*100.0-allele_burden)/tau_allele;
dxdt_STAT5 = kin_stat5*(allele_burden/50.0)*(1.0-jak2_inhib)-kout_stat5*STAT5;
dxdt_BFU_E_mut = k_bfu_prod*STAT5_eff*(1.0-hu_rbc)-k_bfu_mat*BFU_E_mut-ifn_eff*BFU_E_mut;
dxdt_BFU_E_wt = k_bfu_prod_wt*(1.0-hu_rbc)-k_bfu_mat*BFU_E_wt;
dxdt_RBC = k_rbc_prod*(BFU_E_mut+BFU_E_wt)-k_rbc_dest*RBC-phlebotomy_drain;
dxdt_HCT = (RBC*HCT_per_RBC*100.0-HCT)/5.0;
dxdt_PLT = k_plt_prod*STAT5_eff*(1.0-hu_plt)-k_plt_dest*PLT;
dxdt_WBC = k_wbc_prod*STAT5_eff*(1.0-hu_wbc)-k_wbc_dest*WBC;
dxdt_Spleen_vol = k_spl_grow*(STAT5_eff-STAT5_base)*Spleen_vol*(1.0-jak2_inhib)-k_spl_base*(Spleen_vol-Spleen_norm);
dxdt_MPN_SAF = k_tss_rise*(STAT5_eff-STAT5_base)*20.0-k_tss_fall*MPN_SAF*(1.0+jak2_inhib*1.5);
dxdt_MF_score = k_mf*TGFb_driver*STAT5_eff*(MF_max-MF_score);
dxdt_Thromb_hazard = k_thromb_base+k_thromb_HCT*HCT_excess+k_thromb_PLT*PLT_excess-asp_inh*k_thromb_base*0.40;
dxdt_RUX_gut = -ka_rux*RUX_gut;
dxdt_RUX_cent = ka_rux*RUX_gut/Vd1_rux*1000.0-(CL_rux+CLD_rux)/Vd1_rux*RUX_cent+CLD_rux/Vd1_rux*RUX_periph;
dxdt_RUX_periph = CLD_rux/Vd1_rux*RUX_cent-CLD_rux/Vd2_rux*RUX_periph;
dxdt_HU_cent = use_HU*HU_daily_dose/24.0/Vd_hu-CL_hu/Vd_hu*HU_cent;
dxdt_IFN_sc = -ka_ifn*IFN_sc;
dxdt_IFN_cent = ka_ifn*IFN_sc/Vd_ifn-CL_ifn/Vd_ifn*IFN_cent;
dxdt_ASP_cent = use_ASP*ASP_daily_dose/24.0/Vd_asp-CL_asp/Vd_asp*ASP_cent;
dxdt_FED_cent = use_FED*FED_daily_dose/24.0/Vd_fed*1000.0-CL_fed/Vd_fed*FED_cent;

$TABLE
double tss_clamped = (MPN_SAF<0)?0:((MPN_SAF>100)?100:MPN_SAF);
double svr35 = (Spleen_vol <= 850*0.65)?1:0;
double hct_ctrl = (HCT<45.0)?1:0;
double tss50 = (tss_clamped <= 10.0)?1:0;

$CAPTURE HCT allele_burden PLT WBC Spleen_vol tss_clamped MF_score Thromb_hazard
RUX_cent HU_cent IFN_cent FED_cent RBC svr35 hct_ctrl tss50
'

pv_mod <- mcode("pv_shiny", pv_model_code, quiet = TRUE)

## ============================================================
## Helper: Run simulation for a given scenario
## ============================================================
run_pv_sim <- function(sim_years = 5,
                       init_hct = 55, init_allele = 45, init_spleen = 850,
                       init_plt = 550, init_wbc = 12.5,
                       use_phlebotomy = 0, phlebotomy_freq = 0,
                       use_HU = 0, HU_daily_dose = 0,
                       use_RUX = 0, RUX_BID_dose = 0,
                       use_IFN = 0, IFN_q2w_mcg = 0,
                       use_ASP = 0, ASP_daily_dose = 0,
                       use_FED = 0, FED_daily_dose = 0) {

  end_days <- sim_years * 365

  # Compute initial RBC from HCT
  init_rbc <- init_hct / (0.082 * 100)
  # Compute initial BFU-E proportional to allele burden
  init_bfu_mut <- 80.0 * (init_allele / 45)
  init_bfu_wt  <- 40.0 * (1 - (init_allele - 45) / 100)
  init_clone <- init_allele / 200  # simplified
  init_wt    <- 1 - init_clone - 0.1

  # Build IFN events
  ifn_events <- NULL
  if (use_IFN == 1 && IFN_q2w_mcg > 0) {
    dose_IU <- IFN_q2w_mcg * 3e4
    ifn_events <- ev(time = seq(1, end_days, by = 14),
                     cmt = "IFN_sc", amt = dose_IU, evid = 1)
  }

  # Build ruxolitinib BID events
  rux_events <- NULL
  if (use_RUX == 1 && RUX_BID_dose > 0) {
    times_bid <- c(seq(1, end_days, by = 1), seq(1.5, end_days + 0.5, by = 1))
    rux_events <- ev(time = sort(times_bid),
                     cmt = "RUX_gut", amt = RUX_BID_dose, evid = 1)
  }

  # Merge events
  all_events <- NULL
  if (!is.null(ifn_events) && !is.null(rux_events)) {
    all_events <- c(ifn_events, rux_events)
  } else if (!is.null(ifn_events)) {
    all_events <- ifn_events
  } else if (!is.null(rux_events)) {
    all_events <- rux_events
  }

  mod_run <- pv_mod %>%
    param(
      use_phlebotomy = use_phlebotomy, phlebotomy_freq = phlebotomy_freq,
      use_HU = use_HU, HU_daily_dose = HU_daily_dose,
      use_RUX = use_RUX, RUX_BID_dose = RUX_BID_dose,
      use_IFN = use_IFN,
      use_ASP = use_ASP, ASP_daily_dose = ASP_daily_dose,
      use_FED = use_FED, FED_daily_dose = FED_daily_dose
    ) %>%
    init(mut_clone = init_clone, wt_clone = max(0.1, init_wt),
         allele_burden = init_allele, RBC = init_rbc,
         HCT = init_hct, PLT = init_plt, WBC = init_wbc,
         Spleen_vol = init_spleen,
         BFU_E_mut = init_bfu_mut, BFU_E_wt = init_bfu_wt)

  if (!is.null(all_events)) {
    sim_result <- mod_run %>% mrgsim(events = all_events, end = end_days, delta = 1)
  } else {
    sim_result <- mod_run %>% mrgsim(end = end_days, delta = 1)
  }

  as.data.frame(sim_result) %>% mutate(time_yr = time / 365)
}

## ============================================================
## UI
## ============================================================
ui <- dashboardPage(
  skin = "red",

  dashboardHeader(
    title = tags$span(
      tags$img(src = "", width = 0),
      "PV QSP Model"
    ),
    titleWidth = 260
  ),

  dashboardSidebar(
    width = 260,
    sidebarMenu(
      menuItem("Patient Profile",      tabName = "tab_patient",  icon = icon("user-md")),
      menuItem("Drug PK",              tabName = "tab_pk",       icon = icon("flask")),
      menuItem("Hematologic Response", tabName = "tab_heme",     icon = icon("tint")),
      menuItem("Spleen & Symptoms",    tabName = "tab_spleen",   icon = icon("stethoscope")),
      menuItem("Disease Progression",  tabName = "tab_prog",     icon = icon("chart-line")),
      menuItem("Scenario Comparison",  tabName = "tab_scenario", icon = icon("table")),
      menuItem("Biomarker Panel",      tabName = "tab_bio",      icon = icon("vial")),
      menuItem("About / References",   tabName = "tab_about",    icon = icon("info-circle"))
    ),

    hr(),
    h5("  Patient Parameters", style = "color:#F5CBA7; margin-left:10px;"),
    sliderInput("init_hct",    "Baseline HCT (%)",          min = 42, max = 70, value = 55, step = 1),
    sliderInput("init_allele", "JAK2 Allele Burden (%)",    min = 10, max = 100, value = 45, step = 5),
    sliderInput("init_spleen", "Baseline Spleen (mL)",      min = 400, max = 3000, value = 850, step = 50),
    sliderInput("init_plt",    "Baseline PLT (×10⁹/L)",     min = 300, max = 2000, value = 550, step = 50),
    sliderInput("init_wbc",    "Baseline WBC (×10⁹/L)",     min = 7, max = 30, value = 12.5, step = 0.5),
    sliderInput("sim_years",   "Simulation Duration (yr)",  min = 1, max = 10, value = 5, step = 1)
  ),

  dashboardBody(
    tags$head(tags$style(HTML("
      .skin-red .main-header .logo { background-color: #922B21; }
      .skin-red .main-header .navbar { background-color: #C0392B; }
      .content-wrapper { background-color: #F9F9F9; }
      .box-header { background-color: #2C3E50 !important; color: white !important; }
      .nav-tabs-custom > .tab-content { padding: 10px; }
    "))),

    tabItems(

      ## =========================================================
      ## Tab 1: Patient Profile & Risk Stratification
      ## =========================================================
      tabItem(tabName = "tab_patient",
        fluidRow(
          box(title = "WHO 2022 PV Diagnostic Criteria", width = 6, solidHeader = TRUE,
              status = "danger",
              h4("Major Criteria"),
              p("1. Hemoglobin >16.5g/dL (M) / >16g/dL (F), or HCT >49% (M) / >48% (F)"),
              p("2. BM biopsy: hypercellularity (trilineage), pleomorphic MK"),
              p("3. JAK2 V617F or exon 12 mutation"),
              h4("Minor Criterion"),
              p("Subnormal EPO level (<3 IU/L)"),
              br(),
              p(strong("Diagnosis:"), " 3 major, OR first 2 major + minor criterion")
          ),
          box(title = "ELN Risk Stratification", width = 6, solidHeader = TRUE,
              status = "warning",
              DT::dataTableOutput("risk_table")
          )
        ),
        fluidRow(
          box(title = "Patient Summary (Current Parameters)", width = 12, solidHeader = TRUE,
              status = "primary",
              fluidRow(
                column(3, valueBoxOutput("vbox_hct", width = 12)),
                column(3, valueBoxOutput("vbox_allele", width = 12)),
                column(3, valueBoxOutput("vbox_spleen", width = 12)),
                column(3, valueBoxOutput("vbox_plt", width = 12))
              )
          )
        ),
        fluidRow(
          box(title = "PV — Key Pathophysiology", width = 12, solidHeader = TRUE,
              status = "info",
              p(strong("JAK2 V617F mutation (>95% of PV cases):"),
                " Gain-of-function mutation in the JH2 pseudokinase domain removes
                autoinhibition → constitutive JAK2 kinase activity → persistent STAT5/3 activation →
                EPO-independent erythropoiesis, thrombopoietin-independent megakaryopoiesis."),
              p(strong("Thrombosis mechanism:"),
                " Erythrocytosis ↑ blood viscosity; platelet-leukocyte aggregates via P-selectin;
                NETs from neutrophils; endothelial dysfunction (NO ↓); TXA2 from activated platelets.
                Annual risk 2–5% for major CV events."),
              p(strong("MF transformation:"),
                " Abnormal megakaryocytes release TGF-β1, PDGF, VEGF → reticulin fibrosis →
                post-PV myelofibrosis (10–25% at 15yr) → AML (5–10% lifetime).")
          )
        )
      ),

      ## =========================================================
      ## Tab 2: Drug Pharmacokinetics
      ## =========================================================
      tabItem(tabName = "tab_pk",
        fluidRow(
          box(title = "Treatment Selection", width = 3, solidHeader = TRUE, status = "primary",
              checkboxInput("use_phlebotomy",  "Phlebotomy",             FALSE),
              conditionalPanel("input.use_phlebotomy",
                sliderInput("phlebotomy_freq", "Phlebotomy (sessions/yr)", 1, 20, 8, 1)
              ),
              checkboxInput("use_HU",  "Hydroxyurea",     FALSE),
              conditionalPanel("input.use_HU",
                sliderInput("HU_daily_dose", "HU Dose (mg/day)", 500, 2500, 1500, 250)
              ),
              checkboxInput("use_RUX", "Ruxolitinib",     FALSE),
              conditionalPanel("input.use_RUX",
                sliderInput("RUX_BID_dose", "Ruxolitinib (mg BID)", 5, 25, 10, 5)
              ),
              checkboxInput("use_IFN", "Ropeg-IFN α-2b",  FALSE),
              conditionalPanel("input.use_IFN",
                sliderInput("IFN_q2w_mcg", "Ropeg-IFN (mcg q2w)", 50, 250, 100, 50)
              ),
              checkboxInput("use_ASP", "Aspirin 100mg/day", TRUE),
              checkboxInput("use_FED", "Fedratinib",       FALSE),
              conditionalPanel("input.use_FED",
                sliderInput("FED_daily_dose", "Fedratinib (mg/day)", 200, 400, 400, 100)
              ),
              actionButton("run_sim", "Run Simulation", class = "btn-danger btn-block",
                           icon = icon("play"))
          ),
          box(title = "Ruxolitinib PK — Multiple Dose (ng/mL)", width = 9, solidHeader = TRUE,
              status = "info",
              plotlyOutput("plot_pk_rux", height = "350px"),
              p(em("2-compartment model; ka=2.16/h, CL=17.7L/h, Vd=72L, t½~3h.
                   IC50(JAK2)~450ng/mL (dashed line). Calibrated: RESPONSE trial."))
          )
        ),
        fluidRow(
          box(title = "Hydroxyurea PK (mg/L)", width = 6, solidHeader = TRUE, status = "success",
              plotlyOutput("plot_pk_hu", height = "280px"),
              p(em("1-compartment; CL=18.5L/h. EC50~8.5mg/L. Half-life ~3.5h."))
          ),
          box(title = "Ropeginterferon-α2b PK (IU/mL, q2w SC)", width = 6, solidHeader = TRUE,
              status = "primary",
              plotlyOutput("plot_pk_ifn", height = "280px"),
              p(em("SC absorption; ka=0.15/day, t½~80-130h. EC50 clone suppression ~250 IU/mL."))
          )
        )
      ),

      ## =========================================================
      ## Tab 3: Hematologic Response
      ## =========================================================
      tabItem(tabName = "tab_heme",
        fluidRow(
          box(title = "Hematocrit (%)", width = 6, solidHeader = TRUE, status = "danger",
              plotlyOutput("plot_hct", height = "320px"),
              p(em("Target <45% — primary endpoint (CYTOREDUCE). Dashed: HCT=45%."))
          ),
          box(title = "RBC Count (×10¹²/L)", width = 6, solidHeader = TRUE, status = "warning",
              plotlyOutput("plot_rbc", height = "320px"),
              p(em("Normal range: 4.5–5.9 (M) / 4.1–5.1 (F) ×10¹²/L."))
          )
        ),
        fluidRow(
          box(title = "Platelet Count (×10⁹/L)", width = 6, solidHeader = TRUE, status = "success",
              plotlyOutput("plot_plt", height = "320px"),
              p(em("Normal: 150–400. PLT>1500 → acquired VWD / bleeding risk.
                   PLT>400 at year 5 predicts major thrombosis (ECLAP study)."))
          ),
          box(title = "WBC Count (×10⁹/L)", width = 6, solidHeader = TRUE, status = "info",
              plotlyOutput("plot_wbc", height = "320px"),
              p(em("Normal: 4–11. WBC>11 at diagnosis = independent thrombosis risk factor (ECLAP)."))
          )
        )
      ),

      ## =========================================================
      ## Tab 4: Spleen & Symptom Response
      ## =========================================================
      tabItem(tabName = "tab_spleen",
        fluidRow(
          box(title = "Spleen Volume (mL)", width = 6, solidHeader = TRUE, status = "warning",
              plotlyOutput("plot_spleen", height = "320px"),
              p(em("SVR35 = ≥35% reduction from baseline (primary endpoint RESPONSE/RESPONSE-2).
                   Dashed line: SVR35 threshold."))
          ),
          box(title = "MPN-SAF Total Symptom Score (0–100)", width = 6, solidHeader = TRUE,
              status = "primary",
              plotlyOutput("plot_tss", height = "320px"),
              p(em("MPN-SAF-TSS: fatigue, early satiety, abdominal discomfort, inactivity, night sweats,
                   pruritus, bone pain, fever, unintentional weight loss.
                   TSS50 (≥50% reduction) = secondary endpoint RESPONSE."))
          )
        ),
        fluidRow(
          box(title = "Spleen Response Waterfall (Week 32)", width = 6, solidHeader = TRUE,
              status = "info",
              plotlyOutput("plot_spleen_waterfall", height = "280px")
          ),
          box(title = "Key PV Symptoms Overview", width = 6, solidHeader = TRUE, status = "success",
              DT::dataTableOutput("tbl_symptoms")
          )
        )
      ),

      ## =========================================================
      ## Tab 5: Disease Progression
      ## =========================================================
      tabItem(tabName = "tab_prog",
        fluidRow(
          box(title = "JAK2 V617F Allele Burden (%)", width = 6, solidHeader = TRUE, status = "danger",
              plotlyOutput("plot_allele", height = "320px"),
              p(em("Target: <50% (response), <10% (molecular remission).
                   Only IFN reduces allele burden by eliminating the mutant clone.
                   JAK inhibitors suppress signaling without clone elimination."))
          ),
          box(title = "Bone Marrow Fibrosis Score (0–3)", width = 6, solidHeader = TRUE,
              status = "warning",
              plotlyOutput("plot_mf", height = "320px"),
              p(em("Grade ≥2 = post-PV MF. Progressive fibrosis driven by
                   TGF-β1 from abnormal megakaryocytes.
                   10–25% PV patients develop MF at 15 years."))
          )
        ),
        fluidRow(
          box(title = "Cumulative Thrombosis Hazard", width = 6, solidHeader = TRUE,
              status = "danger",
              plotlyOutput("plot_thromb", height = "280px"),
              p(em("Composite of HCT excess, PLT excess, and NO deficiency effects.
                   Annual thrombosis rate ~2.5% (low-risk), >5% (high-risk).
                   Aspirin reduces microvascular events by ~65%."))
          ),
          box(title = "PV → Post-PV MF → AML Pathway", width = 6, solidHeader = TRUE,
              status = "info",
              DT::dataTableOutput("tbl_transformation")
          )
        )
      ),

      ## =========================================================
      ## Tab 6: Scenario Comparison
      ## =========================================================
      tabItem(tabName = "tab_scenario",
        fluidRow(
          box(title = "Scenario Comparison — 6 Treatment Arms", width = 12,
              solidHeader = TRUE, status = "primary",
              p("All 6 scenarios run simultaneously with fixed patient parameters from sidebar."),
              actionButton("run_all_scenarios", "Run All 6 Scenarios",
                           class = "btn-primary", icon = icon("refresh")),
              br(), br()
          )
        ),
        fluidRow(
          box(title = "Hematocrit — All Scenarios", width = 6, solidHeader = TRUE, status = "danger",
              plotlyOutput("plot_comp_hct", height = "300px")
          ),
          box(title = "Allele Burden — All Scenarios", width = 6, solidHeader = TRUE, status = "warning",
              plotlyOutput("plot_comp_allele", height = "300px")
          )
        ),
        fluidRow(
          box(title = "Response at Week 32 (RESPONSE Design)", width = 12,
              solidHeader = TRUE, status = "success",
              DT::dataTableOutput("tbl_response_wk32")
          )
        )
      ),

      ## =========================================================
      ## Tab 7: Biomarker Panel
      ## =========================================================
      tabItem(tabName = "tab_bio",
        fluidRow(
          box(title = "Biomarker Timeline", width = 8, solidHeader = TRUE, status = "primary",
              selectInput("bio_var", "Select Biomarker:",
                          choices = c(
                            "HCT (%)" = "HCT",
                            "Allele Burden (%)" = "allele_burden",
                            "PLT (×10⁹/L)" = "PLT",
                            "WBC (×10⁹/L)" = "WBC",
                            "Spleen Volume (mL)" = "Spleen_vol",
                            "MPN-SAF TSS" = "tss_clamped",
                            "MF Score" = "MF_score",
                            "Ruxolitinib Cp (ng/mL)" = "RUX_cent",
                            "HU Cp (mg/L)" = "HU_cent",
                            "IFN-α Cp (IU/mL)" = "IFN_cent"
                          ),
                          selected = "allele_burden"),
              plotlyOutput("plot_biomarker", height = "350px")
          ),
          box(title = "Biomarker Reference Ranges", width = 4, solidHeader = TRUE, status = "info",
              DT::dataTableOutput("tbl_biomarker_ref")
          )
        ),
        fluidRow(
          box(title = "JAK2 Inhibition Profile (PD)", width = 12, solidHeader = TRUE, status = "warning",
              plotlyOutput("plot_jak2_inhib", height = "300px"),
              p(em("Effective JAK2 inhibition over time (combined ruxolitinib + fedratinib, if active).
                   Target: >80% inhibition for sustained HCT control."))
          )
        )
      ),

      ## =========================================================
      ## Tab 8: About / Model Description
      ## =========================================================
      tabItem(tabName = "tab_about",
        fluidRow(
          box(title = "About this QSP Model", width = 8, solidHeader = TRUE, status = "info",
              h4("Model Overview"),
              p("This Quantitative Systems Pharmacology (QSP) model simulates the
                 pathophysiology and treatment of Polycythemia Vera (PV), a JAK2V617F-driven
                 myeloproliferative neoplasm characterized by uncontrolled erythrocytosis,
                 thrombosis, and risk of transformation to myelofibrosis and AML."),
              h4("Model Structure (22 ODEs)"),
              tags$ul(
                tags$li("Clonal HSC compartments (mutant / wild-type)"),
                tags$li("JAK-STAT5 signaling (allele-burden driven)"),
                tags$li("Erythroid progenitors (BFU-E mutant / wt)"),
                tags$li("Mature blood cells (RBC, HCT, PLT, WBC)"),
                tags$li("Spleen volume (EMH-driven)"),
                tags$li("MPN-SAF symptom score, Fibrosis progression"),
                tags$li("Cumulative thrombosis hazard"),
                tags$li("Drug PK: Ruxolitinib (2-cmt), HU, Ropeg-IFN (SC), Aspirin, Fedratinib")
              ),
              h4("Clinical Calibration"),
              tags$ul(
                tags$li("RESPONSE (NEJM 2015): Ruxolitinib vs BAT — SVR35 38% vs 1%, HCT control 60% vs 20%"),
                tags$li("RESPONSE-2 (Blood 2017): Ruxolitinib in non-enlarged spleen"),
                tags$li("PROUD-PV / CONTINUATION-PV (Blood 2019; Leukemia 2020): Ropeg-IFN allele burden"),
                tags$li("MAJIC-PV (Lancet Haematol 2017): Ruxolitinib vs HU head-to-head"),
                tags$li("CYTOREDUCE (NEJM 2013): HCT <45% reduces thrombosis 4-fold"),
                tags$li("ECLAP (Lancet 2004): Aspirin reduces microvascular events")
              )
          ),
          box(title = "Key References", width = 4, solidHeader = TRUE, status = "success",
              p(tags$a("Vannucchi AM, NEJM 2015", href="#"), " — RESPONSE trial (ruxolitinib)"),
              p(tags$a("Kiladjian JJ, NEJM 2013", href="#"), " — CYTOREDUCE (HCT target)"),
              p(tags$a("Gisslinger H, Leukemia 2020", href="#"), " — PROUD-PV (ropeg-IFN)"),
              p(tags$a("Harrison CN, Lancet Haematol 2017", href="#"), " — MAJIC-PV"),
              p(tags$a("Finazzi G, Lancet 2004", href="#"), " — ECLAP aspirin trial"),
              p(tags$a("Barbui T, Leukemia 2018", href="#"), " — ELN recommendations 2018"),
              p(tags$a("Arber DA, Blood 2022", href="#"), " — WHO 2022 classification"),
              p(tags$a("Mughal TI, Leukemia 2018", href="#"), " — JAK2 V617F biology review"),
              hr(),
              p("Directory: ", code("polycythemia-vera/")),
              p("Author: Catholic Univ. Seoul — QSP CCR 2026-06-24"),
              p("Model version: 1.0 | Compartments: 22 | Scenarios: 6")
          )
        )
      )
    )
  )
)

## ============================================================
## Server
## ============================================================
server <- function(input, output, session) {

  ## ---- Reactive: Run user-selected simulation ----
  sim_data <- eventReactive(input$run_sim, {
    withProgress(message = "Running PV simulation...", {
      run_pv_sim(
        sim_years      = input$sim_years,
        init_hct       = input$init_hct,
        init_allele    = input$init_allele,
        init_spleen    = input$init_spleen,
        init_plt       = input$init_plt,
        init_wbc       = input$init_wbc,
        use_phlebotomy = as.integer(input$use_phlebotomy),
        phlebotomy_freq = if(input$use_phlebotomy) input$phlebotomy_freq else 0,
        use_HU         = as.integer(input$use_HU),
        HU_daily_dose  = if(input$use_HU) input$HU_daily_dose else 0,
        use_RUX        = as.integer(input$use_RUX),
        RUX_BID_dose   = if(input$use_RUX) input$RUX_BID_dose else 0,
        use_IFN        = as.integer(input$use_IFN),
        IFN_q2w_mcg    = if(input$use_IFN) input$IFN_q2w_mcg else 0,
        use_ASP        = as.integer(input$use_ASP),
        ASP_daily_dose = if(input$use_ASP) 100 else 0,
        use_FED        = as.integer(input$use_FED),
        FED_daily_dose = if(input$use_FED) input$FED_daily_dose else 0
      )
    })
  }, ignoreNULL = FALSE)

  ## ---- Reactive: Run all 6 pre-defined scenarios ----
  all_scenarios_data <- eventReactive(input$run_all_scenarios, {
    withProgress(message = "Running 6 PV scenarios...", {
      params <- list(
        init_hct = input$init_hct, init_allele = input$init_allele,
        init_spleen = input$init_spleen, init_plt = input$init_plt,
        init_wbc = input$init_wbc, sim_years = input$sim_years
      )

      s0 <- do.call(run_pv_sim, c(params)) %>% mutate(scenario = "S0: Natural History")
      s1 <- do.call(run_pv_sim, c(params, list(use_phlebotomy=1, phlebotomy_freq=8, use_ASP=1, ASP_daily_dose=100))) %>%
            mutate(scenario = "S1: Phlebotomy+ASA")
      s2 <- do.call(run_pv_sim, c(params, list(use_HU=1, HU_daily_dose=1500, use_phlebotomy=1, phlebotomy_freq=3, use_ASP=1, ASP_daily_dose=100))) %>%
            mutate(scenario = "S2: HU+Phlebotomy+ASA")
      s3 <- do.call(run_pv_sim, c(params, list(use_RUX=1, RUX_BID_dose=10, use_phlebotomy=1, phlebotomy_freq=2, use_ASP=1, ASP_daily_dose=100))) %>%
            mutate(scenario = "S3: Ruxolitinib 10mg BID")
      s4 <- do.call(run_pv_sim, c(params, list(use_IFN=1, IFN_q2w_mcg=100, use_phlebotomy=1, phlebotomy_freq=4, use_ASP=1, ASP_daily_dose=100))) %>%
            mutate(scenario = "S4: Ropeg-IFN 100mcg q2w")
      s5 <- do.call(run_pv_sim, c(params, list(use_FED=1, FED_daily_dose=400, use_phlebotomy=1, phlebotomy_freq=2, use_ASP=1, ASP_daily_dose=100))) %>%
            mutate(scenario = "S5: Fedratinib 400mg/day")

      bind_rows(s0, s1, s2, s3, s4, s5)
    })
  })

  ## ---- Value boxes ----
  output$vbox_hct    <- renderValueBox(valueBox(
    paste0(input$init_hct, "%"), "Baseline HCT",
    icon = icon("tint"), color = if(input$init_hct > 49) "red" else "yellow"))
  output$vbox_allele <- renderValueBox(valueBox(
    paste0(input$init_allele, "%"), "JAK2 Allele Burden",
    icon = icon("dna"), color = if(input$init_allele > 50) "red" else "orange"))
  output$vbox_spleen <- renderValueBox(valueBox(
    paste0(input$init_spleen, " mL"), "Spleen Volume",
    icon = icon("stethoscope"), color = if(input$init_spleen > 1000) "red" else "yellow"))
  output$vbox_plt    <- renderValueBox(valueBox(
    paste0(input$init_plt, "×10⁹/L"), "Platelet Count",
    icon = icon("circle"), color = if(input$init_plt > 600) "orange" else "green"))

  ## ---- Risk stratification table ----
  output$risk_table <- DT::renderDataTable({
    data.frame(
      Risk       = c("Low", "High"),
      Criteria   = c("Age <60yr AND no prior thrombosis AND PLT <1500",
                     "Age ≥60yr OR prior thrombosis"),
      Management = c("Phlebotomy + Low-dose aspirin (100mg/day)",
                     "Cytoreductive therapy (HU or Ruxolitinib) + Aspirin"),
      Thrombosis_Rate = c("~2.5%/yr", "~5%/yr")
    ) %>% setNames(c("Risk Category", "Criteria", "First-line Management", "Thrombosis Rate"))
  }, options = list(dom = 't', scrollX = TRUE), rownames = FALSE)

  ## ---- Hematology plots ----
  make_line_plot <- function(data, y_var, y_label, y_int = NULL, y_int_label = NULL,
                             title_str = "", color = "#2980B9") {
    p <- plot_ly(data, x = ~time_yr, y = as.formula(paste0("~", y_var)),
                 type = "scatter", mode = "lines",
                 line = list(color = color, width = 2),
                 name = y_label) %>%
      layout(title = title_str, xaxis = list(title = "Time (years)"),
             yaxis = list(title = y_label))
    if (!is.null(y_int)) {
      p <- p %>% add_lines(x = ~time_yr, y = y_int, line = list(dash = "dash", color = "red", width = 1.5),
                           name = y_int_label %||% "Threshold")
    }
    p
  }

  output$plot_hct <- renderPlotly({
    req(sim_data())
    df <- sim_data()
    plot_ly(df, x = ~time_yr, y = ~HCT, type = "scatter", mode = "lines",
            line = list(color = "#E74C3C", width = 2.5), name = "HCT (%)") %>%
      add_lines(x = ~time_yr, y = 45, line = list(dash = "dash", color = "black", width = 1.5),
                name = "Target <45%") %>%
      layout(xaxis = list(title = "Time (years)"), yaxis = list(title = "HCT (%)"),
             showlegend = TRUE)
  })

  output$plot_rbc <- renderPlotly({
    req(sim_data())
    df <- sim_data()
    plot_ly(df, x = ~time_yr, y = ~RBC, type = "scatter", mode = "lines",
            line = list(color = "#C0392B", width = 2.5), name = "RBC (×10¹²/L)") %>%
      layout(xaxis = list(title = "Time (years)"), yaxis = list(title = "RBC (×10¹²/L)"))
  })

  output$plot_plt <- renderPlotly({
    req(sim_data())
    df <- sim_data()
    plot_ly(df, x = ~time_yr, y = ~PLT, type = "scatter", mode = "lines",
            line = list(color = "#148F77", width = 2.5), name = "PLT") %>%
      add_lines(x = ~time_yr, y = 400, line = list(dash = "dash", color = "orange", width = 1.5),
                name = "Normal upper limit (400)") %>%
      add_lines(x = ~time_yr, y = 1500, line = list(dash = "dot", color = "red", width = 1.5),
                name = "Bleeding risk (1500)") %>%
      layout(xaxis = list(title = "Time (years)"), yaxis = list(title = "PLT (×10⁹/L)"))
  })

  output$plot_wbc <- renderPlotly({
    req(sim_data())
    df <- sim_data()
    plot_ly(df, x = ~time_yr, y = ~WBC, type = "scatter", mode = "lines",
            line = list(color = "#27AE60", width = 2.5), name = "WBC") %>%
      add_lines(x = ~time_yr, y = 11, line = list(dash = "dash", color = "orange", width = 1.5),
                name = "Upper normal (11)") %>%
      layout(xaxis = list(title = "Time (years)"), yaxis = list(title = "WBC (×10⁹/L)"))
  })

  ## ---- Spleen & Symptoms ----
  output$plot_spleen <- renderPlotly({
    req(sim_data())
    df <- sim_data()
    svr_thr <- input$init_spleen * 0.65
    plot_ly(df, x = ~time_yr, y = ~Spleen_vol, type = "scatter", mode = "lines",
            line = list(color = "#F39C12", width = 2.5), name = "Spleen Vol (mL)") %>%
      add_lines(x = ~time_yr, y = svr_thr, line = list(dash = "dash", color = "#2980B9", width = 1.5),
                name = "SVR35 threshold") %>%
      layout(xaxis = list(title = "Time (years)"), yaxis = list(title = "Spleen Volume (mL)"))
  })

  output$plot_tss <- renderPlotly({
    req(sim_data())
    df <- sim_data()
    tss_thr <- 20.0 * 0.5
    plot_ly(df, x = ~time_yr, y = ~tss_clamped, type = "scatter", mode = "lines",
            line = list(color = "#8E44AD", width = 2.5), name = "MPN-SAF TSS") %>%
      add_lines(x = ~time_yr, y = tss_thr, line = list(dash = "dash", color = "#E74C3C", width = 1.5),
                name = "TSS50 threshold") %>%
      layout(xaxis = list(title = "Time (years)"), yaxis = list(title = "MPN-SAF TSS (0–100)"))
  })

  output$plot_spleen_waterfall <- renderPlotly({
    req(sim_data())
    df <- sim_data()
    wk32 <- df %>% filter(abs(time - 224) < 1)
    if (nrow(wk32) == 0) return(NULL)
    pct_change <- (wk32$Spleen_vol[1] - input$init_spleen) / input$init_spleen * 100
    plot_ly(x = c("Baseline", "Week 32"),
            y = c(input$init_spleen, wk32$Spleen_vol[1]),
            type = "bar",
            marker = list(color = c("#F39C12", if(pct_change < -35) "#27AE60" else "#E74C3C"))) %>%
      layout(yaxis = list(title = "Spleen Volume (mL)"),
             title = paste0("Spleen Change at Wk32: ", round(pct_change, 1), "%"))
  })

  output$tbl_symptoms <- DT::renderDataTable({
    data.frame(
      Symptom = c("Fatigue", "Pruritus (aquagenic)", "Headache", "Plethora", "Splenomegaly",
                  "Erythromelalgia", "Night sweats", "Visual disturbances", "DVT/PE"),
      Frequency  = c("85%", "40–70%", "60%", "50%", "40%", "5–30%", "25%", "10%", "5–30%"),
      Mechanism  = c("Cytokines (IL-1β, TNF-α)", "Histamine + IL-13 (basophils)", "Hyperviscosity",
                     "Erythrocytosis", "EMH (JAK2-driven)", "TXA2/platelet occlusion",
                     "TNF-α", "Microvascular", "HCT/PLT/WBC ↑")
    )
  }, options = list(dom = 't', pageLength = 9, scrollX = TRUE), rownames = FALSE)

  ## ---- Disease Progression ----
  output$plot_allele <- renderPlotly({
    req(sim_data())
    df <- sim_data()
    plot_ly(df, x = ~time_yr, y = ~allele_burden, type = "scatter", mode = "lines",
            line = list(color = "#E74C3C", width = 2.5), name = "Allele Burden (%)") %>%
      add_lines(x = ~time_yr, y = 50, line = list(dash = "dot", color = "red", width = 1.5),
                name = "50% (high burden)") %>%
      add_lines(x = ~time_yr, y = 10, line = list(dash = "dash", color = "green", width = 1.5),
                name = "10% (mol. remission)") %>%
      layout(xaxis = list(title = "Time (years)"), yaxis = list(title = "JAK2 V617F Allele Burden (%)"))
  })

  output$plot_mf <- renderPlotly({
    req(sim_data())
    df <- sim_data()
    plot_ly(df, x = ~time_yr, y = ~MF_score, type = "scatter", mode = "lines",
            line = list(color = "#7D3C98", width = 2.5), name = "Fibrosis Grade") %>%
      add_lines(x = ~time_yr, y = 1.0, line = list(dash = "dot", color = "orange", width = 1.5),
                name = "MF-1 boundary") %>%
      add_lines(x = ~time_yr, y = 2.0, line = list(dash = "dash", color = "red", width = 1.5),
                name = "MF-2 (Post-PV MF)") %>%
      layout(xaxis = list(title = "Time (years)"), yaxis = list(title = "WHO Fibrosis Grade (0–3)"))
  })

  output$plot_thromb <- renderPlotly({
    req(sim_data())
    df <- sim_data()
    plot_ly(df, x = ~time_yr, y = ~Thromb_hazard, type = "scatter", mode = "lines",
            line = list(color = "#C0392B", width = 2.5), name = "Cumul. Thrombosis Hazard") %>%
      layout(xaxis = list(title = "Time (years)"), yaxis = list(title = "Cumulative Hazard"))
  })

  output$tbl_transformation <- DT::renderDataTable({
    data.frame(
      Stage = c("Polycythemia Vera", "Post-PV MF (MF-2/3)", "Blast Phase (AML)"),
      Prevalence = c("~44/100,000", "~10–25% at 15yr", "~5–10% lifetime"),
      Median_OS  = c("~14 years", "~5 years", "<1 year"),
      Key_Driver = c("JAK2V617F, HCT↑, PLT↑", "TGF-β1, reticulin fibrosis grade ≥2", "IDH1/2, ASXL1 mutations"),
      Treatment  = c("Phlebotomy/HU/RUX/IFN", "Ruxolitinib (MF dose), HSCT", "AML-directed therapy, HSCT")
    )
  }, options = list(dom = 't', scrollX = TRUE), rownames = FALSE)

  ## ---- PK Plots ----
  output$plot_pk_rux <- renderPlotly({
    # Simulate 4 days BID 10mg ruxolitinib
    ev_rux <- ev(time = c(seq(0, 96, 24), seq(12, 108, 24)),
                 cmt = "RUX_gut", amt = 10, evid = 1)
    sim_pk <- pv_mod %>%
      param(use_RUX = 1, RUX_BID_dose = 10) %>%
      mrgsim(events = ev_rux, end = 120, delta = 0.5) %>%
      as.data.frame()

    plot_ly(sim_pk, x = ~time, y = ~RUX_cent, type = "scatter", mode = "lines",
            line = list(color = "#2980B9", width = 2), name = "Ruxolitinib Cp") %>%
      add_lines(x = ~time, y = 450, line = list(dash = "dash", color = "red", width = 1.5),
                name = "IC50 JAK2 (~450 ng/mL)") %>%
      layout(xaxis = list(title = "Time (hours)"), yaxis = list(title = "Cp (ng/mL)"),
             title = "Ruxolitinib 10mg BID — Multiple Dose PK")
  })

  output$plot_pk_hu <- renderPlotly({
    sim_hu <- pv_mod %>%
      param(use_HU = 1, HU_daily_dose = 1500) %>%
      mrgsim(end = 72, delta = 0.25) %>%
      as.data.frame()

    plot_ly(sim_hu, x = ~time, y = ~HU_cent, type = "scatter", mode = "lines",
            line = list(color = "#27AE60", width = 2), name = "HU Cp (mg/L)") %>%
      add_lines(x = ~time, y = 8.5, line = list(dash = "dash", color = "orange", width = 1.5),
                name = "EC50 (~8.5 mg/L)") %>%
      layout(xaxis = list(title = "Time (hours)"), yaxis = list(title = "Cp (mg/L)"),
             title = "Hydroxyurea 1500mg/day — Steady-state PK")
  })

  output$plot_pk_ifn <- renderPlotly({
    ev_ifn <- ev(time = c(0, 14, 28), cmt = "IFN_sc", amt = 100 * 3e4, evid = 1)
    sim_ifn <- pv_mod %>%
      param(use_IFN = 1) %>%
      mrgsim(events = ev_ifn, end = 42, delta = 0.1) %>%
      as.data.frame()

    plot_ly(sim_ifn, x = ~time, y = ~IFN_cent, type = "scatter", mode = "lines",
            line = list(color = "#1B4F72", width = 2), name = "IFN-α Cp (IU/mL)") %>%
      add_lines(x = ~time, y = 250, line = list(dash = "dash", color = "purple", width = 1.5),
                name = "EC50 clone suppression (~250 IU/mL)") %>%
      layout(xaxis = list(title = "Time (days)"), yaxis = list(title = "IFN-α Cp (IU/mL)"),
             title = "Ropeginterferon-α2b 100mcg q2w SC — PK")
  })

  ## ---- Scenario Comparison ----
  output$plot_comp_hct <- renderPlotly({
    req(all_scenarios_data())
    df <- all_scenarios_data()
    plot_ly(df, x = ~time_yr, y = ~HCT, color = ~scenario,
            type = "scatter", mode = "lines") %>%
      add_lines(x = ~range(df$time_yr), y = c(45, 45),
                line = list(dash = "dash", color = "black", width = 1), name = "Target <45%",
                inherit = FALSE) %>%
      layout(xaxis = list(title = "Time (years)"), yaxis = list(title = "HCT (%)"))
  })

  output$plot_comp_allele <- renderPlotly({
    req(all_scenarios_data())
    df <- all_scenarios_data()
    plot_ly(df, x = ~time_yr, y = ~allele_burden, color = ~scenario,
            type = "scatter", mode = "lines") %>%
      layout(xaxis = list(title = "Time (years)"), yaxis = list(title = "JAK2 Allele Burden (%)"))
  })

  output$tbl_response_wk32 <- DT::renderDataTable({
    req(all_scenarios_data())
    df <- all_scenarios_data()
    df %>%
      filter(abs(time - 224) < 1) %>%
      mutate(
        HCT_ctrl  = ifelse(HCT < 45, "Yes ✓", "No ✗"),
        SVR35     = ifelse(Spleen_vol < input$init_spleen * 0.65, "Yes ✓", "No ✗"),
        TSS50     = ifelse(tss_clamped < 10.0, "Yes ✓", "No ✗"),
        CHR       = ifelse(HCT < 45 & PLT < 400 & WBC < 10, "Yes ✓", "No ✗")
      ) %>%
      select(scenario, HCT, PLT, WBC, Spleen_vol, tss_clamped, allele_burden,
             HCT_ctrl, SVR35, TSS50, CHR) %>%
      mutate(across(where(is.numeric), ~round(.x, 1))) %>%
      setNames(c("Scenario", "HCT(%)", "PLT(10⁹/L)", "WBC(10⁹/L)", "Spleen(mL)",
                 "MPN-SAF TSS", "Allele(%)", "HCT Control", "SVR35", "TSS50", "CHR"))
  }, options = list(scrollX = TRUE), rownames = FALSE)

  ## ---- Biomarkers ----
  output$plot_biomarker <- renderPlotly({
    req(sim_data())
    df <- sim_data()
    yvar <- input$bio_var
    if (!yvar %in% names(df)) return(NULL)
    plot_ly(df, x = ~time_yr, y = df[[yvar]], type = "scatter", mode = "lines",
            line = list(color = "#2C3E50", width = 2), name = yvar) %>%
      layout(xaxis = list(title = "Time (years)"), yaxis = list(title = yvar))
  })

  output$tbl_biomarker_ref <- DT::renderDataTable({
    data.frame(
      Biomarker = c("HCT", "Allele Burden", "PLT", "WBC", "Spleen", "TSS", "MF Score", "Rux Cp"),
      Normal_Range = c("<45% (M)", "<2%", "150–400", "4–11", "<500mL", "<10", "MF-0", "–"),
      PV_Range     = c("49–70%", "20–100%", "400–2000", "10–30", "500–3000mL", "10–50", "MF-0→3", "–"),
      Intervention_Target = c("<45%", "<50% response\n<10% mol. remission", "<400", "<10", "SVR≥35%", "TSS50", "MF-0/1", ">IC50 (450ng/mL)")
    )
  }, options = list(dom = 't', pageLength = 8), rownames = FALSE)

  output$plot_jak2_inhib <- renderPlotly({
    req(sim_data())
    df <- sim_data()
    if ("RUX_cent" %in% names(df)) {
      rux_inh <- 0.92 * df$RUX_cent^1.2 / (450^1.2 + df$RUX_cent^1.2)
      plot_ly(df, x = ~time_yr, y = rux_inh * 100, type = "scatter", mode = "lines",
              line = list(color = "#2980B9", width = 2), name = "JAK2 Inhibition (%)") %>%
        add_lines(x = ~time_yr, y = 80, line = list(dash = "dash", color = "red", width = 1.5),
                  name = "80% threshold") %>%
        layout(xaxis = list(title = "Time (years)"), yaxis = list(title = "JAK2 Inhibition (%)"),
               title = "Effective JAK2 Inhibition over Time")
    }
  })
}

## ============================================================
## Run the application
## ============================================================
shinyApp(ui = ui, server = server)
