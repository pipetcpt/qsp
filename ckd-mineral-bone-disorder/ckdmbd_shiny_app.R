## ============================================================
## CKD-MBD QSP Shiny Application
## FGF23 · Klotho · PTH · Vitamin D · Bone · Vascular Calcification
## 7 Tabs: Patient Profile · PK · PTH/Mineral · Bone · CV · Scenarios · Biomarkers
## ============================================================

library(shiny)
library(shinydashboard)
library(mrgsolve)
library(dplyr)
library(tidyr)
library(ggplot2)
library(plotly)
library(DT)

# ─── Inline mrgsolve model ────────────────────────────────────────────────────
ckdmbd_code <- '
$PARAM @annotated
GFR0         : 15    : Baseline GFR (mL/min/1.73m2)
Pi_in        : 1200  : Dietary phosphate intake (mg/day)
Ca_in        : 800   : Dietary calcium intake (mg/day)
kabs_Pi      : 0.60  : Fractional Pi absorption
kPi_urine    : 0.0008: Renal Pi clearance per GFR
kFGF23_syn   : 0.02  : FGF23 synthesis rate
kFGF23_deg   : 0.35  : FGF23 degradation
EC50_Pi_FGF23: 5.5   : Pi EC50 for FGF23 stimulation
Emax_Pi_FGF23: 8.0   : Max Pi-driven FGF23 fold
kKlotho_deg  : 0.12  : Klotho degradation
kKlotho_GFR  : 0.04  : GFR-dependent Klotho synthesis
kPTH_syn     : 0.50  : PTH synthesis rate
kPTH_deg     : 1.20  : PTH degradation
EC50_Ca_PTH  : 1.10  : Ca EC50 for PTH suppression
IC50_VitD_PTH: 30    : VitD IC50 for PTH
IC50_CaSR_cin: 15    : Cinacalcet IC50 for CaSR
Hill_PTH     : 3.5   : Hill coefficient
k25_syn      : 0.08  : 25-OH-D synthesis
k25_deg      : 0.02  : 25-OH-D catabolism
k125_syn     : 0.015 : Calcitriol synthesis
k125_deg     : 0.30  : Calcitriol catabolism
Imax_FGF23_CYP27B1 : 0.80 : FGF23 max CYP27B1 inhibition
IC50_FGF23_CYP27B1 : 200  : FGF23 IC50 CYP27B1
kabs_Ca      : 0.30  : Ca absorption fraction
kCa_urine    : 0.006 : Renal Ca excretion
kCa_bone     : 0.008 : Bone Ca exchange
kOB_syn      : 0.008 : Osteoblast synthesis
kOB_deg      : 0.02  : Osteoblast apoptosis
kOC_syn      : 0.006 : Osteoclast synthesis
kOC_deg      : 0.03  : Osteoclast apoptosis
kBMD_form    : 0.0003: BMD formation rate
kBMD_res     : 0.0004: BMD resorption rate
BMD_ss       : 1.0   : Normal BMD
kVASC_calc   : 0.0001: Vascular Ca-Pi deposition
VC_threshold : 55    : Ca×Pi calcification threshold
Sev_Emax     : 0.65  : Sevelamer Pi binding max
Sev_EC50     : 800   : Sevelamer ED50 (mg/d)
Cin_ka       : 1.2   : Cinacalcet absorption
Cin_F        : 0.21  : Cinacalcet bioavailability
Cin_V        : 1000  : Cinacalcet volume
Cin_CL       : 250   : Cinacalcet clearance
Par_V        : 34    : Paricalcitol volume
Par_CL       : 17    : Paricalcitol clearance
Par_Emax     : 0.90  : Paricalcitol VDR Emax
Par_EC50     : 0.2   : Paricalcitol EC50
Etel_V       : 11.5  : Etelcalcetide volume
Etel_CL      : 0.60  : Etelcalcetide clearance
Etel_Emax    : 0.85  : Etelcalcetide CaSR Emax
Etel_EC50    : 50    : Etelcalcetide EC50
Den_ka       : 0.006 : Denosumab absorption
Den_F        : 0.62  : Denosumab bioavailability
Den_V        : 2.8   : Denosumab volume
Den_CL       : 0.008 : Denosumab clearance
Den_Imax     : 0.95  : Denosumab RANKL inhibition
DOSE_Sev     : 0     : Sevelamer dose (mg/d)
DOSE_CaCO3   : 0     : CaCO3 elemental Ca (mg/d)

$CMT @annotated
Pi          : Serum phosphate (mg/dL)
FGF23       : Plasma FGF23 (pg/mL)
Klotho      : Soluble Klotho (rel.)
PTH         : Intact PTH (pg/mL)
VitD25      : 25-OH-D (nmol/L)
VitD_act    : 1,25-OH2D calcitriol (pg/mL)
Ca          : Serum calcium (mg/dL)
OB          : Osteoblast activity
OC          : Osteoclast activity
BMD         : Bone mineral density (rel.)
VC          : Vascular calcification (AU)
CIN_GUT     : Cinacalcet gut (mg)
CIN_PLASMA  : Cinacalcet plasma (ng/mL)
PAR_PLASMA  : Paricalcitol plasma (ng/mL)
ETEL_PLASMA : Etelcalcetide plasma (ng/mL)
DEN_DEPOT   : Denosumab depot (mg)
DEN_PLASMA  : Denosumab plasma (mg/L)

$MAIN
double Ca_ion   = Ca / 2.51;
double CaSR_act = pow(Ca_ion, Hill_PTH) / (pow(EC50_Ca_PTH, Hill_PTH) + pow(Ca_ion, Hill_PTH));
double Cin_CaSR  = CIN_PLASMA / (IC50_CaSR_cin + CIN_PLASMA);
double Etel_CaSR = ETEL_PLASMA / (Etel_EC50 + ETEL_PLASMA);
double CaSR_total = CaSR_act + (1.0 - CaSR_act) * (Cin_CaSR * 0.7 + Etel_CaSR * Etel_Emax);
double VitD_VDR  = VitD_act / (VitD_act + IC50_VitD_PTH);
double Par_VDR   = (PAR_PLASMA * Par_Emax) / (PAR_PLASMA + Par_EC50);
double VDR_act   = fmin(1.0, VitD_VDR + Par_VDR);
double Sev_bind  = (DOSE_Sev * Sev_Emax) / (DOSE_Sev + Sev_EC50);
double CaCO3_Ca  = DOSE_CaCO3 * 0.40 / 10.0;
double Pi_stim   = 1.0 + Emax_Pi_FGF23 * pow(fmax(Pi - 3.5, 0), 2) /
    (pow(EC50_Pi_FGF23 - 3.5, 2) + pow(fmax(Pi - 3.5, 0), 2));
double GFR_act   = GFR0 * (1.0 - 0.005 * (VC + 0.001));
double kPi_clear = kPi_urine * GFR_act;
double FGF23_inh = 1.0 - Imax_FGF23_CYP27B1 * FGF23 / (IC50_FGF23_CYP27B1 + FGF23);
double PTH_stim  = 1.0 + 0.5 * PTH / (PTH + 100);
double Den_RANKL = (DEN_PLASMA * Den_Imax) / (DEN_PLASMA + 0.03);
double RANKL_eff = kOC_syn * OB * (1.0 - 0.5 * VDR_act) * (1.0 - Den_RANKL);
double Scler_eff = fmin(1.0, 0.3 + 0.5 * (1.0 - GFR_act / 90.0));

$ODE
dxdt_Pi        = (kabs_Pi * (1.0 - Sev_bind) * Pi_in / 100.0) + 0.005 * BMD
               - kPi_clear * Pi - 0.002 * Pi;
dxdt_FGF23     = kFGF23_syn * Pi_stim * (1.0 + VitD_act / 50.0) - kFGF23_deg * FGF23;
dxdt_Klotho    = kKlotho_GFR * GFR_act - kKlotho_deg * Klotho;
dxdt_PTH       = kPTH_syn * (1.0 - CaSR_total) * (1.0 - VDR_act * 0.6)
               + 0.02 * Pi_stim - kPTH_deg * PTH;
dxdt_VitD25    = k25_syn - k25_deg * VitD25;
dxdt_VitD_act  = k125_syn * VitD25 * FGF23_inh * PTH_stim
               + PAR_PLASMA * Par_Emax * 0.1
               - k125_deg * VitD_act
               - 0.05 * VitD_act * VitD_act / (50.0 + VitD_act);
dxdt_Ca        = (kabs_Ca * (1.0 + 0.5 * VDR_act) * Ca_in / 200.0) + CaCO3_Ca
               + 0.003 * PTH * BMD - kCa_urine * GFR_act * Ca
               - kCa_bone * (OB - OC * 0.7);
dxdt_OB        = kOB_syn * (1.0 + 0.4 * VDR_act) * (1.0 - 0.3 * Scler_eff)
               - kOB_deg * OB;
dxdt_OC        = RANKL_eff * (1.0 + 0.6 * PTH / (PTH + 150)) - kOC_deg * OC;
dxdt_BMD       = kBMD_form * OB - kBMD_res * OC - 0.0001 * (PTH / 65 - 1) * BMD;
dxdt_VC        = kVASC_calc * fmax(0, Ca * Pi - VC_threshold) * 0.5
               + 0.00002 * Pi * (Pi - 3.5) - 0.0002 * VC;
dxdt_CIN_GUT    = -Cin_ka * CIN_GUT;
dxdt_CIN_PLASMA = Cin_ka * CIN_GUT * Cin_F / Cin_V - Cin_CL / Cin_V * CIN_PLASMA;
dxdt_PAR_PLASMA = -Par_CL / Par_V * PAR_PLASMA;
dxdt_ETEL_PLASMA = -Etel_CL / Etel_V * ETEL_PLASMA;
dxdt_DEN_DEPOT  = -Den_ka * DEN_DEPOT;
dxdt_DEN_PLASMA = Den_ka * DEN_DEPOT * Den_F / Den_V - Den_CL / Den_V * DEN_PLASMA;

$TABLE
double iPTH  = PTH;
double sPi   = Pi;
double sCa   = Ca;
double s25D  = VitD25;
double s125D = VitD_act;
double sFGF23 = FGF23;
double sKlotho = Klotho;
double sBMD  = BMD;
double sVC   = VC;
double CaP   = Ca * Pi;
double sGFR  = GFR_act;

$CAPTURE iPTH sPi sCa s25D s125D sFGF23 sKlotho sBMD sVC CaP sGFR
'

mod <- mcode("ckdmbd_shiny", ckdmbd_code, quiet = TRUE)

# ─── Default initial conditions (CKD G5) ─────────────────────────────────────
default_init <- list(
  Pi = 6.2, FGF23 = 800, Klotho = 0.35, PTH = 420,
  VitD25 = 18, VitD_act = 15, Ca = 8.8,
  OB = 0.7, OC = 1.4, BMD = 0.82, VC = 25,
  CIN_GUT = 0, CIN_PLASMA = 0, PAR_PLASMA = 0,
  ETEL_PLASMA = 0, DEN_DEPOT = 0, DEN_PLASMA = 0
)

run_sim <- function(params, dose_sev = 0, dose_cin = 0,
                    dose_par = 0, par_freq = 2,
                    dose_etel = 0, etel_freq = 2,
                    dose_den = 0, den_freq = 180,
                    gfr0 = 15, sim_days = 365) {

  tg <- tgrid(0, sim_days, 1)
  events <- NULL

  if (dose_cin > 0) {
    events <- c(events, list(ev(cmt = "CIN_GUT", amt = dose_cin,
                                ii = 24, addl = sim_days, time = 0)))
  }
  if (dose_par > 0) {
    events <- c(events, list(ev(cmt = "PAR_PLASMA", amt = dose_par,
                                ii = par_freq * 24,
                                addl = floor(sim_days / par_freq), time = 0)))
  }
  if (dose_etel > 0) {
    events <- c(events, list(ev(cmt = "ETEL_PLASMA", amt = dose_etel * 1000,
                                ii = etel_freq * 24,
                                addl = floor(sim_days / etel_freq), time = 0)))
  }
  if (dose_den > 0) {
    events <- c(events, list(ev(cmt = "DEN_DEPOT", amt = dose_den,
                                ii = den_freq * 24,
                                addl = floor(sim_days / den_freq), time = 0)))
  }

  all_ev <- if (!is.null(events)) Reduce(function(a, b) ev_seq(a, b), events) else NULL

  m <- mod %>%
    init(default_init) %>%
    param(GFR0 = gfr0, DOSE_Sev = dose_sev)

  if (!is.null(all_ev)) {
    m <- m %>% mrgsim(events = all_ev, tgrid = tg)
  } else {
    m <- m %>% mrgsim(tgrid = tg)
  }
  as_tibble(m)
}

# ─── KDIGO targets ────────────────────────────────────────────────────────────
kdigo_targets <- data.frame(
  Analyte    = c("iPTH (pg/mL)", "Pi (mg/dL)", "Ca (mg/dL)", "25-OH-D (nmol/L)", "Ca×Pi"),
  Lower      = c(150, NA, 8.4, 75, NA),
  Upper      = c(600, 5.5, 10.2, NA, 55),
  KDIGO_note = c("2–9× ULN for CKD G5D", "Normal range", "Normal range",
                 "Sufficiency", "Calcification risk"),
  stringsAsFactors = FALSE
)

# ═══════════════════════════════════════════════════════════════════════════════
# UI
# ═══════════════════════════════════════════════════════════════════════════════
ui <- dashboardPage(
  skin = "blue",
  dashboardHeader(title = "CKD-MBD QSP Model"),
  dashboardSidebar(
    sidebarMenu(
      menuItem("Patient Profile",    tabName = "tab_patient",   icon = icon("user")),
      menuItem("Drug PK",            tabName = "tab_pk",        icon = icon("pills")),
      menuItem("PTH & Minerals",     tabName = "tab_pthminerals", icon = icon("flask")),
      menuItem("Bone Disease",       tabName = "tab_bone",      icon = icon("bone")),
      menuItem("Cardiovascular",     tabName = "tab_cv",        icon = icon("heartbeat")),
      menuItem("Scenario Comparison",tabName = "tab_scenarios", icon = icon("chart-bar")),
      menuItem("Biomarkers",         tabName = "tab_biomarkers", icon = icon("vials"))
    )
  ),
  dashboardBody(
    tabItems(

      # ── Tab 1: Patient Profile ─────────────────────────────────────────────
      tabItem(tabName = "tab_patient",
        fluidRow(
          box(title = "CKD Stage & Patient Parameters", width = 4, status = "primary",
            sliderInput("gfr0", "Baseline GFR (mL/min/1.73m²)", 5, 90, 15, step = 5),
            sliderInput("sim_days", "Simulation Duration (days)", 90, 730, 365, step = 30),
            sliderInput("pi_in", "Dietary Pi Intake (mg/day)", 600, 2000, 1200, step = 100),
            sliderInput("ca_in", "Dietary Ca Intake (mg/day)", 400, 1500, 800, step = 100),
            hr(),
            actionButton("run_sim", "Run Simulation", class = "btn-primary btn-block")
          ),
          box(title = "CKD Stage Reference", width = 4, status = "info",
            tableOutput("ckd_table")
          ),
          box(title = "Baseline Lab Summary", width = 4, status = "warning",
            tableOutput("baseline_labs")
          )
        ),
        fluidRow(
          box(title = "Disease Progression Overview (iPTH & FGF23)", width = 12,
            plotlyOutput("patient_overview_plot", height = "320px")
          )
        )
      ),

      # ── Tab 2: Drug PK ─────────────────────────────────────────────────────
      tabItem(tabName = "tab_pk",
        fluidRow(
          box(title = "Phosphate Binders", width = 4, status = "primary",
            sliderInput("dose_sev", "Sevelamer (mg/day)", 0, 7200, 2400, step = 400),
            sliderInput("dose_caco3", "CaCO₃ Elemental Ca (mg/day)", 0, 2000, 0, step = 200),
            h5("Note: Sevelamer — non-absorbed; CaCO₃ — absorbed")
          ),
          box(title = "Calcimimetics", width = 4, status = "warning",
            selectInput("calci_type", "Calcimimetic Agent",
                        choices = c("None", "Cinacalcet (oral)", "Etelcalcetide (IV)")),
            conditionalPanel("input.calci_type == 'Cinacalcet (oral)'",
              sliderInput("dose_cin", "Cinacalcet Dose (mg/day)", 0, 180, 60, step = 30)
            ),
            conditionalPanel("input.calci_type == 'Etelcalcetide (IV)'",
              sliderInput("dose_etel", "Etelcalcetide Dose (mg/session)", 0, 15, 5, step = 1),
              sliderInput("etel_freq", "Dosing Interval (days)", 1, 7, 2, step = 1)
            )
          ),
          box(title = "Vitamin D & Other", width = 4, status = "success",
            checkboxInput("use_par", "Paricalcitol (IV)", FALSE),
            conditionalPanel("input.use_par",
              sliderInput("dose_par", "Paricalcitol (mcg/session)", 1, 8, 4, step = 1),
              sliderInput("par_freq", "Dosing Interval (days)", 1, 7, 2, step = 1)
            ),
            checkboxInput("use_den", "Denosumab (SC)", FALSE),
            conditionalPanel("input.use_den",
              sliderInput("dose_den", "Denosumab (mg)", 30, 120, 60, step = 30),
              sliderInput("den_freq", "Dosing Interval (days)", 90, 360, 180, step = 30)
            )
          )
        ),
        fluidRow(
          box(title = "Drug Plasma Concentration — Calcimimetic", width = 6,
            plotlyOutput("pk_cin_plot", height = "280px")),
          box(title = "Drug Plasma Concentration — Paricalcitol", width = 6,
            plotlyOutput("pk_par_plot", height = "280px"))
        ),
        fluidRow(
          box(title = "Drug Plasma Concentration — Denosumab", width = 6,
            plotlyOutput("pk_den_plot", height = "280px")),
          box(title = "PK Summary (Cmax, Tmax, AUC)", width = 6,
            tableOutput("pk_summary_table"))
        )
      ),

      # ── Tab 3: PTH & Mineral Metabolism ───────────────────────────────────
      tabItem(tabName = "tab_pthminerals",
        fluidRow(
          box(title = "Intact PTH (pg/mL)", width = 6,
            plotlyOutput("pth_plot", height = "280px")),
          box(title = "Serum Phosphate (mg/dL)", width = 6,
            plotlyOutput("pi_plot", height = "280px"))
        ),
        fluidRow(
          box(title = "Serum Calcium (mg/dL)", width = 6,
            plotlyOutput("ca_plot", height = "280px")),
          box(title = "Ca × Pi Product", width = 6,
            plotlyOutput("cap_plot", height = "280px"))
        ),
        fluidRow(
          box(title = "FGF23 (pg/mL)", width = 6,
            plotlyOutput("fgf23_plot", height = "280px")),
          box(title = "Vitamin D Status", width = 6,
            plotlyOutput("vitd_plot", height = "280px"))
        ),
        fluidRow(
          box(title = "KDIGO 2017 Target Attainment at Day 365", width = 12,
            tableOutput("kdigo_attainment_table"))
        )
      ),

      # ── Tab 4: Bone Disease ────────────────────────────────────────────────
      tabItem(tabName = "tab_bone",
        fluidRow(
          box(title = "Bone Mineral Density (rel. to normal)", width = 6,
            plotlyOutput("bmd_plot", height = "300px")),
          box(title = "Osteoblast vs Osteoclast Activity", width = 6,
            plotlyOutput("ob_oc_plot", height = "300px"))
        ),
        fluidRow(
          box(title = "Bone Remodeling Markers", width = 12,
            p(strong("P1NP (Formation): "), "∝ Osteoblast activity"),
            p(strong("CTX (Resorption): "), "∝ Osteoclast activity"),
            plotlyOutput("bone_markers_plot", height = "280px")
          )
        ),
        fluidRow(
          box(title = "Bone Disease Classification Guide", width = 6, status = "info",
            tableOutput("bone_class_table")),
          box(title = "Fracture Risk Estimate", width = 6,
            plotlyOutput("fracture_risk_plot", height = "250px"))
        )
      ),

      # ── Tab 5: Cardiovascular ─────────────────────────────────────────────
      tabItem(tabName = "tab_cv",
        fluidRow(
          box(title = "Vascular Calcification Score (AU)", width = 6,
            plotlyOutput("vc_plot", height = "300px")),
          box(title = "GFR Trajectory", width = 6,
            plotlyOutput("gfr_plot", height = "300px"))
        ),
        fluidRow(
          box(title = "Soluble Klotho (rel. units)", width = 6,
            plotlyOutput("klotho_plot", height = "280px")),
          box(title = "CV Risk Assessment", width = 6, status = "danger",
            tableOutput("cv_risk_table"))
        ),
        fluidRow(
          box(title = "Key CV Drivers — Correlation Heatmap", width = 12,
            plotlyOutput("cv_heatmap", height = "280px"))
        )
      ),

      # ── Tab 6: Scenario Comparison ────────────────────────────────────────
      tabItem(tabName = "tab_scenarios",
        fluidRow(
          box(title = "Predefined Treatment Scenarios", width = 12, status = "primary",
            p("Comparing 5 standard-of-care strategies at 365 days."),
            plotlyOutput("scenario_pth_plot", height = "280px")
          )
        ),
        fluidRow(
          box(title = "Scenario — Serum Phosphate", width = 6,
            plotlyOutput("scenario_pi_plot", height = "250px")),
          box(title = "Scenario — BMD", width = 6,
            plotlyOutput("scenario_bmd_plot", height = "250px"))
        ),
        fluidRow(
          box(title = "Scenario — Vascular Calcification", width = 6,
            plotlyOutput("scenario_vc_plot", height = "250px")),
          box(title = "Scenario Summary Table (Day 365)", width = 6,
            DTOutput("scenario_summary_dt"))
        )
      ),

      # ── Tab 7: Biomarkers ─────────────────────────────────────────────────
      tabItem(tabName = "tab_biomarkers",
        fluidRow(
          box(title = "Biomarker Dashboard — Key Lab Values Over Time", width = 12,
            plotlyOutput("biomarker_all_plot", height = "400px")
          )
        ),
        fluidRow(
          box(title = "Biomarker Radar Chart (Day 90 vs Day 365)", width = 6,
            plotlyOutput("biomarker_radar", height = "320px")),
          box(title = "Biomarker Reference Ranges", width = 6, status = "info",
            tableOutput("biomarker_ref_table"))
        ),
        fluidRow(
          box(title = "FGF23 & Klotho Longitudinal", width = 12,
            plotlyOutput("fgf23_klotho_plot", height = "280px"))
        )
      )
    ) # end tabItems
  ) # end dashboardBody
)

# ═══════════════════════════════════════════════════════════════════════════════
# SERVER
# ═══════════════════════════════════════════════════════════════════════════════
server <- function(input, output, session) {

  sim_data <- eventReactive(input$run_sim, {
    req(input$gfr0, input$sim_days)
    dose_cin  <- if (!is.null(input$calci_type) && input$calci_type == "Cinacalcet (oral)") input$dose_cin else 0
    dose_etel <- if (!is.null(input$calci_type) && input$calci_type == "Etelcalcetide (IV)") input$dose_etel else 0
    dose_par  <- if (isTRUE(input$use_par)) input$dose_par else 0
    dose_den  <- if (isTRUE(input$use_den)) input$dose_den else 0

    run_sim(params = list(),
            dose_sev  = input$dose_sev,
            dose_cin  = dose_cin,
            dose_par  = dose_par, par_freq = input$par_freq %||% 2,
            dose_etel = dose_etel, etel_freq = input$etel_freq %||% 2,
            dose_den  = dose_den, den_freq = input$den_freq %||% 180,
            gfr0 = input$gfr0, sim_days = input$sim_days)
  }, ignoreNULL = FALSE)

  # Pre-compute scenarios once
  scenarios_data <- reactive({
    scens <- list(
      "S1: No Treatment"           = run_sim(list(), gfr0 = 15, sim_days = 365),
      "S2: Sevelamer 2400mg/d"     = run_sim(list(), dose_sev = 2400, gfr0 = 15, sim_days = 365),
      "S3: Cinacalcet 60mg/d"      = run_sim(list(), dose_cin = 60, gfr0 = 15, sim_days = 365),
      "S4: Paricalcitol 4mcg 3x/wk"= run_sim(list(), dose_par = 4, par_freq = 2, gfr0 = 15, sim_days = 365),
      "S5: Sev + Cinacalcet"       = run_sim(list(), dose_sev = 2400, dose_cin = 60, gfr0 = 15, sim_days = 365)
    )
    bind_rows(mapply(function(d, n) mutate(d, Scenario = n), scens, names(scens), SIMPLIFY = FALSE))
  })

  # ── Patient Overview ──────────────────────────────────────────────────────
  output$ckd_table <- renderTable({
    data.frame(
      Stage = c("G1", "G2", "G3a", "G3b", "G4", "G5"),
      GFR   = c("≥90", "60–89", "45–59", "30–44", "15–29", "<15"),
      Risk  = c("Normal", "Mild↓", "Mild-Mod↓", "Mod-Sev↓", "Severe↓", "Kidney Failure")
    )
  }, striped = TRUE, bordered = TRUE, hover = TRUE)

  output$baseline_labs <- renderTable({
    data.frame(
      Lab       = c("iPTH", "Pi", "Ca", "25-OH-D", "FGF23", "Klotho"),
      Value     = c("420 pg/mL", "6.2 mg/dL", "8.8 mg/dL", "18 nmol/L", "800 pg/mL", "0.35 rel."),
      Target    = c("150–600", "<5.5", "8.4–10.2", ">75", "<200", ">0.5"),
      Status    = c("↑", "↑", "Low-normal", "↓", "↑↑", "↓")
    )
  }, striped = TRUE, bordered = TRUE)

  output$patient_overview_plot <- renderPlotly({
    d <- sim_data()
    p <- plot_ly(d, x = ~time) %>%
      add_lines(y = ~iPTH,   name = "iPTH (pg/mL)",   line = list(color = "red")) %>%
      add_lines(y = ~sFGF23, name = "FGF23 (pg/mL)",  line = list(color = "orange")) %>%
      layout(title = "iPTH & FGF23 Over Time",
             xaxis = list(title = "Day"),
             yaxis = list(title = "pg/mL"))
    p
  })

  # ── PK Plots ─────────────────────────────────────────────────────────────
  output$pk_cin_plot <- renderPlotly({
    d <- sim_data()
    plot_ly(d, x = ~time, y = ~CIN_PLASMA, type = "scatter", mode = "lines",
            line = list(color = "steelblue")) %>%
      layout(title = "Cinacalcet Plasma (ng/mL)", xaxis = list(title = "Day"),
             yaxis = list(title = "ng/mL"))
  })
  output$pk_par_plot <- renderPlotly({
    d <- sim_data()
    plot_ly(d, x = ~time, y = ~PAR_PLASMA, type = "scatter", mode = "lines",
            line = list(color = "darkorange")) %>%
      layout(title = "Paricalcitol Plasma (ng/mL)", xaxis = list(title = "Day"),
             yaxis = list(title = "ng/mL"))
  })
  output$pk_den_plot <- renderPlotly({
    d <- sim_data()
    plot_ly(d, x = ~time, y = ~DEN_PLASMA, type = "scatter", mode = "lines",
            line = list(color = "purple")) %>%
      layout(title = "Denosumab Plasma (mg/L)", xaxis = list(title = "Day"),
             yaxis = list(title = "mg/L"))
  })
  output$pk_summary_table <- renderTable({
    d <- sim_data()
    data.frame(
      Drug = c("Cinacalcet", "Paricalcitol", "Denosumab"),
      Cmax = c(round(max(d$CIN_PLASMA, na.rm=TRUE), 2),
               round(max(d$PAR_PLASMA, na.rm=TRUE), 2),
               round(max(d$DEN_PLASMA, na.rm=TRUE), 3)),
      Tmax_day = c(d$time[which.max(d$CIN_PLASMA)],
                   d$time[which.max(d$PAR_PLASMA)],
                   d$time[which.max(d$DEN_PLASMA)])
    )
  }, striped = TRUE, bordered = TRUE)

  # ── PTH & Minerals ────────────────────────────────────────────────────────
  mk_plotly_line <- function(d, y_var, title_str, ylab,
                             lo = NULL, hi = NULL, color = "steelblue") {
    p <- plot_ly(d, x = ~time, y = ~get(y_var), type = "scatter", mode = "lines",
                 name = title_str, line = list(color = color)) %>%
      layout(title = title_str, xaxis = list(title = "Day"),
             yaxis = list(title = ylab))
    if (!is.null(lo)) p <- p %>% add_lines(x = range(d$time), y = c(lo, lo),
                                           line = list(dash = "dash", color = "green4"),
                                           name = "Target Low", showlegend = FALSE)
    if (!is.null(hi)) p <- p %>% add_lines(x = range(d$time), y = c(hi, hi),
                                           line = list(dash = "dash", color = "red3"),
                                           name = "Target High", showlegend = FALSE)
    p
  }

  output$pth_plot  <- renderPlotly({ mk_plotly_line(sim_data(), "iPTH",  "iPTH",     "pg/mL",  150, 600, "red3") })
  output$pi_plot   <- renderPlotly({ mk_plotly_line(sim_data(), "sPi",   "Serum Pi", "mg/dL",  NA,  5.5, "darkorange") })
  output$ca_plot   <- renderPlotly({ mk_plotly_line(sim_data(), "sCa",   "Serum Ca", "mg/dL",  8.4, 10.2, "steelblue") })
  output$cap_plot  <- renderPlotly({ mk_plotly_line(sim_data(), "CaP",   "Ca×Pi",    "mg²/dL²", NA, 55, "brown") })
  output$fgf23_plot<- renderPlotly({ mk_plotly_line(sim_data(), "sFGF23","FGF23",    "pg/mL",  NA,  NA, "darkorange") })
  output$vitd_plot <- renderPlotly({
    d <- sim_data()
    plot_ly(d, x = ~time) %>%
      add_lines(y = ~s25D,  name = "25-OH-D (nmol/L)",     line = list(color = "gold")) %>%
      add_lines(y = ~s125D, name = "1,25-OH₂D (pg/mL)",    line = list(color = "orange")) %>%
      layout(title = "Vitamin D Status", xaxis = list(title = "Day"),
             yaxis = list(title = "Concentration"))
  })
  output$kdigo_attainment_table <- renderTable({
    d <- sim_data()
    d365 <- d[which.min(abs(d$time - 365)), ]
    data.frame(
      Parameter   = c("iPTH (pg/mL)", "Pi (mg/dL)", "Ca (mg/dL)", "Ca×Pi", "25-OH-D (nmol/L)"),
      Value_Day365 = round(c(d365$iPTH, d365$sPi, d365$sCa, d365$CaP, d365$s25D), 1),
      Target       = c("150–600", "<5.5", "8.4–10.2", "<55", ">75"),
      Met          = c(
        ifelse(d365$iPTH >= 150 & d365$iPTH <= 600, "YES", "NO"),
        ifelse(d365$sPi <= 5.5, "YES", "NO"),
        ifelse(d365$sCa >= 8.4 & d365$sCa <= 10.2, "YES", "NO"),
        ifelse(d365$CaP < 55, "YES", "NO"),
        ifelse(d365$s25D > 75, "YES", "NO")
      )
    )
  }, striped = TRUE, bordered = TRUE)

  # ── Bone ─────────────────────────────────────────────────────────────────
  output$bmd_plot <- renderPlotly({
    d <- sim_data()
    plot_ly(d, x = ~time, y = ~sBMD, type = "scatter", mode = "lines",
            line = list(color = "saddlebrown")) %>%
      add_lines(x = range(d$time), y = c(1, 1), line = list(dash = "dash", color = "green4"),
                showlegend = FALSE) %>%
      layout(title = "Bone Mineral Density", xaxis = list(title = "Day"),
             yaxis = list(title = "BMD (rel. to normal)"))
  })
  output$ob_oc_plot <- renderPlotly({
    d <- sim_data()
    plot_ly(d, x = ~time) %>%
      add_lines(y = ~OB_act, name = "Osteoblast", line = list(color = "blue")) %>%
      add_lines(y = ~OC_act, name = "Osteoclast", line = list(color = "red")) %>%
      layout(title = "Bone Cell Activity", xaxis = list(title = "Day"),
             yaxis = list(title = "Relative Activity"))
  })
  output$bone_markers_plot <- renderPlotly({
    d <- sim_data() %>% mutate(P1NP = OB_act * 80, CTX = OC_act * 0.6)
    plot_ly(d, x = ~time) %>%
      add_lines(y = ~P1NP, name = "P1NP (ng/mL)", line = list(color = "blue")) %>%
      add_lines(y = ~CTX,  name = "CTX (ng/mL)",  line = list(color = "red")) %>%
      layout(title = "Bone Turnover Markers", xaxis = list(title = "Day"),
             yaxis = list(title = "ng/mL"))
  })
  output$bone_class_table <- renderTable({
    data.frame(
      Type = c("Osteitis Fibrosa", "Adynamic Bone", "Mixed Uremic", "Osteomalacia"),
      Turnover = c("High", "Low", "Mixed", "Low"),
      iPTH     = c(">600", "<150", "150–600", "Variable"),
      Cause    = c("2°HPT", "Oversuppression", "Both", "VitD def.")
    )
  }, striped = TRUE, bordered = TRUE)
  output$fracture_risk_plot <- renderPlotly({
    d <- sim_data() %>% mutate(Fracture_risk = 100 * (1 - sBMD) * (1 + OC_act))
    plot_ly(d, x = ~time, y = ~Fracture_risk, type = "scatter", mode = "lines",
            line = list(color = "red3")) %>%
      layout(title = "Fracture Risk Index", xaxis = list(title = "Day"),
             yaxis = list(title = "Risk Index (AU)"))
  })

  # ── CV Tab ────────────────────────────────────────────────────────────────
  output$vc_plot <- renderPlotly({
    d <- sim_data()
    plot_ly(d, x = ~time, y = ~sVC, type = "scatter", mode = "lines",
            line = list(color = "darkred")) %>%
      layout(title = "Vascular Calcification", xaxis = list(title = "Day"),
             yaxis = list(title = "Calcification Score (AU)"))
  })
  output$gfr_plot <- renderPlotly({
    d <- sim_data()
    plot_ly(d, x = ~time, y = ~sGFR, type = "scatter", mode = "lines",
            line = list(color = "steelblue")) %>%
      layout(title = "GFR Trajectory", xaxis = list(title = "Day"),
             yaxis = list(title = "GFR (mL/min/1.73m²)"))
  })
  output$klotho_plot <- renderPlotly({
    d <- sim_data()
    plot_ly(d, x = ~time, y = ~sKlotho, type = "scatter", mode = "lines",
            line = list(color = "purple")) %>%
      layout(title = "Soluble Klotho", xaxis = list(title = "Day"),
             yaxis = list(title = "Klotho (rel. units)"))
  })
  output$cv_risk_table <- renderTable({
    d <- sim_data()
    d365 <- d[which.min(abs(d$time - 365)), ]
    data.frame(
      Factor      = c("Ca×Pi product", "FGF23", "Klotho deficiency", "Vascular Calc. Score"),
      Value       = round(c(d365$CaP, d365$sFGF23, 1 - d365$sKlotho, d365$sVC), 1),
      Risk_Cutoff = c(">55", ">800", ">0.65", ">50"),
      Risk_Level  = c(
        ifelse(d365$CaP > 55, "High", "Normal"),
        ifelse(d365$sFGF23 > 800, "High", "Moderate"),
        ifelse(1 - d365$sKlotho > 0.65, "High", "Moderate"),
        ifelse(d365$sVC > 50, "High", "Moderate")
      )
    )
  }, striped = TRUE, bordered = TRUE)
  output$cv_heatmap <- renderPlotly({
    d <- sim_data() %>% filter(time %% 30 == 0)
    plot_ly(d, x = ~time) %>%
      add_lines(y = ~sVC/max(sVC), name = "VC (norm)", line = list(color = "red")) %>%
      add_lines(y = ~CaP/80,       name = "Ca×Pi (norm)", line = list(color = "orange")) %>%
      add_lines(y = ~sFGF23/max(sFGF23), name = "FGF23 (norm)", line = list(color = "brown")) %>%
      add_lines(y = ~sKlotho,      name = "Klotho (rel)", line = list(color = "purple")) %>%
      layout(title = "CV Risk Drivers (Normalized)", xaxis = list(title = "Day"),
             yaxis = list(title = "Normalized Value (0–1)"))
  })

  # ── Scenario Comparison ───────────────────────────────────────────────────
  output$scenario_pth_plot <- renderPlotly({
    d <- scenarios_data()
    colors <- c("gray40", "steelblue", "red3", "darkorange", "green4")
    p <- plot_ly()
    for (i in seq_along(unique(d$Scenario))) {
      sc <- unique(d$Scenario)[i]
      sub <- filter(d, Scenario == sc)
      p <- p %>% add_lines(data = sub, x = ~time, y = ~iPTH, name = sc,
                            line = list(color = colors[i]))
    }
    p %>% layout(title = "iPTH by Scenario",
                 shapes = list(
                   list(type = "line", x0 = 0, x1 = 365, y0 = 150, y1 = 150,
                        line = list(dash = "dash", color = "green4")),
                   list(type = "line", x0 = 0, x1 = 365, y0 = 600, y1 = 600,
                        line = list(dash = "dash", color = "red3"))
                 ),
                 xaxis = list(title = "Day"),
                 yaxis = list(title = "iPTH (pg/mL)"))
  })
  output$scenario_pi_plot <- renderPlotly({
    d <- scenarios_data()
    plot_ly(d, x = ~time, y = ~sPi, color = ~Scenario, type = "scatter", mode = "lines") %>%
      layout(title = "Serum Pi", xaxis = list(title = "Day"),
             yaxis = list(title = "Pi (mg/dL)"))
  })
  output$scenario_bmd_plot <- renderPlotly({
    d <- scenarios_data()
    plot_ly(d, x = ~time, y = ~sBMD, color = ~Scenario, type = "scatter", mode = "lines") %>%
      layout(title = "BMD", xaxis = list(title = "Day"),
             yaxis = list(title = "BMD (rel.)"))
  })
  output$scenario_vc_plot <- renderPlotly({
    d <- scenarios_data()
    plot_ly(d, x = ~time, y = ~sVC, color = ~Scenario, type = "scatter", mode = "lines") %>%
      layout(title = "Vascular Calcification", xaxis = list(title = "Day"),
             yaxis = list(title = "VC Score (AU)"))
  })
  output$scenario_summary_dt <- renderDT({
    d <- scenarios_data() %>%
      filter(time == 365) %>%
      select(Scenario, iPTH, sPi, sCa, CaP, sBMD, sVC) %>%
      mutate(across(where(is.numeric), ~round(., 1)))
    datatable(d, options = list(pageLength = 5, dom = "t"),
              rownames = FALSE)
  })

  # ── Biomarkers ────────────────────────────────────────────────────────────
  output$biomarker_all_plot <- renderPlotly({
    d <- sim_data()
    plot_ly(d, x = ~time) %>%
      add_lines(y = ~iPTH / max(iPTH, na.rm=TRUE),   name = "iPTH (norm)") %>%
      add_lines(y = ~sPi  / 10,                       name = "Pi/10") %>%
      add_lines(y = ~sCa  / 10,                       name = "Ca/10") %>%
      add_lines(y = ~sFGF23 / max(sFGF23, na.rm=TRUE), name = "FGF23 (norm)") %>%
      add_lines(y = ~sKlotho,                          name = "Klotho (rel)") %>%
      add_lines(y = ~sBMD,                             name = "BMD (rel)") %>%
      layout(title = "All Biomarkers (Normalized)", xaxis = list(title = "Day"),
             yaxis = list(title = "Normalized Value"))
  })
  output$biomarker_radar <- renderPlotly({
    d <- sim_data()
    d90  <- d[which.min(abs(d$time - 90)), ]
    d365 <- d[which.min(abs(d$time - 365)), ]
    cats <- c("PTH ctrl", "Pi ctrl", "Ca ctrl", "FGF23 ctrl", "Klotho", "BMD", "VC ctrl")
    vals90 <- c(
      1 - pmin(d90$iPTH, 600) / 600,
      1 - d90$sPi / 8,
      pmin(d90$sCa / 10, 1),
      1 - pmin(d90$sFGF23, 1000) / 1000,
      d90$sKlotho,
      d90$sBMD,
      1 - d90$sVC / 100
    )
    vals365 <- c(
      1 - pmin(d365$iPTH, 600) / 600,
      1 - d365$sPi / 8,
      pmin(d365$sCa / 10, 1),
      1 - pmin(d365$sFGF23, 1000) / 1000,
      d365$sKlotho,
      d365$sBMD,
      1 - d365$sVC / 100
    )
    plot_ly(type = "scatterpolar", mode = "lines+markers", fill = "toself") %>%
      add_trace(r = vals90,  theta = cats, name = "Day 90") %>%
      add_trace(r = vals365, theta = cats, name = "Day 365") %>%
      layout(title = "Biomarker Spider (higher = better)", polar = list(radialaxis = list(visible = TRUE, range = c(0, 1))))
  })
  output$biomarker_ref_table <- renderTable({
    data.frame(
      Biomarker   = c("iPTH", "Phosphate", "Calcium", "25-OH-D", "FGF23", "Klotho"),
      Normal      = c("15–65 pg/mL", "2.5–4.5 mg/dL", "8.5–10.2 mg/dL", ">75 nmol/L", "<100 pg/mL", ">1.0 rel."),
      CKD_G5      = c("150–600 pg/mL", "<5.5 mg/dL", "8.4–10.2 mg/dL", ">50 nmol/L", "600–10000 pg/mL", "0.1–0.4 rel."),
      Source      = c("KDIGO 2017", "KDIGO 2017", "KDIGO 2017", "KDIGO 2017", "Ix 2011", "Hu 2011")
    )
  }, striped = TRUE, bordered = TRUE)
  output$fgf23_klotho_plot <- renderPlotly({
    d <- sim_data()
    plot_ly(d, x = ~time) %>%
      add_lines(y = ~sFGF23,  name = "FGF23 (pg/mL)",    line = list(color = "orange")) %>%
      add_lines(y = ~sKlotho * 200, name = "Klotho×200 (rel)", line = list(color = "purple")) %>%
      layout(title = "FGF23 & Klotho Axis", xaxis = list(title = "Day"),
             yaxis = list(title = "FGF23 (pg/mL) | Klotho×200"))
  })
}

# ─── Helpers ─────────────────────────────────────────────────────────────────
`%||%` <- function(a, b) if (!is.null(a)) a else b

shinyApp(ui, server)
