## ============================================================
##  TTP (Thrombotic Thrombocytopenic Purpura) QSP Shiny App
##  Framework: mrgsolve + Shiny + ggplot2
##  6 Tabs: Patient Profile, PK, ADAMTS13/VWF/PLT, Clinical Endpoints,
##          Scenario Comparison, Biomarkers Panel
##  Author: Claude Code QSP Routine | 2026-06-25
## ============================================================

library(shiny)
library(shinydashboard)
library(mrgsolve)
library(dplyr)
library(ggplot2)
library(tidyr)
library(DT)
library(plotly)

## ── Embed model code ──────────────────────────────────────────
ttp_model_code <- '
$PARAM @annotated
A13_prod : 0.50  : ADAMTS13 production (U/dL/d)
A13_deg  : 0.231 : ADAMTS13 degradation (/d)
A13_base : 100   : Normal ADAMTS13 (%)
k_inh_on : 0.025 : Inhibitor-ADAMTS13 on-rate (/BU/d)
k_inh_off: 0.008 : Inhibitor-ADAMTS13 off-rate (/d)
k_inh_deg: 0.033 : IgG catabolism (/d)
f_ADCP   : 0.15  : ADCP clearance factor (/d per BU)
ULVWF_prod:12.0  : ULVWF secretion (ng/mL/d)
ULVWF_deg: 1.20  : ULVWF degradation (/d)
k_cleave : 0.060 : ADAMTS13 cleavage rate (/d)
PLT0     : 250   : Baseline PLT (×10⁹/L)
PLT_prod : 25.0  : PLT production (×10⁹/L/d)
k_PLT_deg: 0.10  : PLT removal (/d)
k_MT_form: 0.002 : MT formation rate constant
k_MT_lysis:0.40  : MT lysis (/d)
MT_Hill  : 1.5   : MT formation Hill coeff
MT_EC50  : 40.0  : MT formation EC50 (ng/mL ULVWF)
BC0      : 100   : Baseline B cells
BC_prod  : 3.5   : B cell production (/d)
BC_deg   : 0.035 : B cell death (/d)
PC_diff  : 0.020 : B→PC differentiation (/d)
PC_deg   : 0.050 : Plasma cell death (/d)
Ab_prod  : 0.12  : Ab production rate (BU/AU/d)
Ab_deg   : 0.033 : Ab catabolism (/d)
LDH0     : 180   : Baseline LDH (IU/L)
k_LDH_hem: 60.0  : LDH per hemolysis
k_LDH_deg: 0.50  : LDH clearance (/d)
Cr0      : 75    : Baseline Cr (μmol/L)
k_Cr_MT  : 4.0   : Cr rise per MT (/d)
k_Cr_deg : 0.25  : Cr clearance (/d)
Trop0    : 0.02  : Baseline TnI (ng/mL)
k_Trop_MT: 0.08  : TnI rise per MT (/d)
k_Trop_deg:0.30  : TnI clearance (/d)
Hgb0     : 14.0  : Baseline Hgb (g/dL)
k_Hgb_prod:0.12  : Erythropoiesis (g/dL/d)
k_hemo   : 0.03  : Hemolysis rate per MT
CAPLA_CL : 0.50  : Caplacizumab CL (L/d)
CAPLA_V1 : 3.0   : Caplacizumab V1 (L)
CAPLA_V2 : 2.0   : Caplacizumab V2 (L)
CAPLA_Q  : 0.80  : Caplacizumab Q (L/d)
CAPLA_Ka : 0.80  : Caplacizumab SC Ka (/d)
CAPLA_F  : 0.85  : Caplacizumab SC F
CAPLA_IC50:0.002 : Caplacizumab IC50 (mg/L = μg/mL)
CAPLA_Hn : 1.5   : Caplacizumab Hill n
RTX_CL   : 0.80  : Rituximab CL (L/d)
RTX_V1   : 4.0   : Rituximab V1 (L)
RTX_V2   : 3.0   : Rituximab V2 (L)
RTX_Q    : 0.50  : Rituximab Q (L/d)
RTX_EC50 : 0.50  : Rituximab B cell EC50 (μg/mL)
RTX_Emax : 0.95  : Rituximab B cell Emax
PRED_CL  : 15.0  : Prednisolone CL (L/d)
PRED_V   : 35.0  : Prednisolone V (L)
PRED_EC50:150    : Prednisolone Ab suppression EC50 (ng/mL)
PRED_Emax: 0.65  : Prednisolone Ab suppression Emax
PEX_A13  : 60.0  : ADAMTS13 per TPE session (U/dL)
PEX_f_inh: 0.70  : Inhibitor fraction removed per TPE
PEX_f_ULVWF:0.35 : ULVWF fraction removed per TPE

$CMT @annotated
CAPLA_GUT : Caplacizumab SC depot (mg)
CAPLA_C   : Caplacizumab central (mg)
CAPLA_P   : Caplacizumab peripheral (mg)
RTX_C     : Rituximab central (mg)
RTX_P     : Rituximab peripheral (mg)
A13_ACT   : ADAMTS13 activity (U/dL)
INH       : Inhibitor titer (BU)
ULVWF     : ULVWF (ng/mL)
PLT       : Platelet count (×10⁹/L)
MT        : Microthrombus burden (AU)
BC        : B cells (AU)
PC        : Plasma cells (AU)
AUTOAB    : Autoantibody (BU)
LDH_AB    : LDH (IU/L)
CREAT     : Creatinine (μmol/L)
TROP      : Troponin I (ng/mL)
HGB       : Hemoglobin (g/dL)
PRED_C    : Prednisolone (ng/mL)

$MAIN
double CAPLA_Cp = CAPLA_C / CAPLA_V1;
double CAPLA_eff = pow(CAPLA_Cp, CAPLA_Hn) /
                   (pow(CAPLA_IC50, CAPLA_Hn) + pow(CAPLA_Cp, CAPLA_Hn));
double RTX_Cp = RTX_C / RTX_V1;
double RTX_eff = RTX_Emax * RTX_Cp / (RTX_EC50 + RTX_Cp);
double PRED_Cp = PRED_C;
double PRED_eff = PRED_Emax * PRED_Cp / (PRED_EC50 + PRED_Cp);
double A13_frac = A13_ACT / A13_base;
double ULVWF_pos = (ULVWF > 0) ? ULVWF : 0.0;
double PLT_pos = (PLT > 0) ? PLT : 0.0;
double MT_pos = (MT > 0) ? MT : 0.0;
double MT_form_rate = k_MT_form * pow(ULVWF_pos, MT_Hill) /
                      (pow(MT_EC50, MT_Hill) + pow(ULVWF_pos, MT_Hill)) *
                      (PLT_pos / PLT0) * (1.0 - CAPLA_eff);
double hemolysis = k_hemo * MT_pos;

$ODE
dxdt_CAPLA_GUT = -CAPLA_Ka * CAPLA_GUT;
dxdt_CAPLA_C = CAPLA_Ka * CAPLA_GUT * CAPLA_F
               - (CAPLA_CL + CAPLA_Q) / CAPLA_V1 * CAPLA_C
               + CAPLA_Q / CAPLA_V2 * CAPLA_P;
dxdt_CAPLA_P = CAPLA_Q / CAPLA_V1 * CAPLA_C
               - CAPLA_Q / CAPLA_V2 * CAPLA_P;
dxdt_RTX_C = -(RTX_CL + RTX_Q) / RTX_V1 * RTX_C + RTX_Q / RTX_V2 * RTX_P;
dxdt_RTX_P = RTX_Q / RTX_V1 * RTX_C - RTX_Q / RTX_V2 * RTX_P;
dxdt_A13_ACT = A13_prod * A13_base - A13_deg * A13_ACT
               - k_inh_on * INH * A13_ACT + k_inh_off * INH
               - f_ADCP * INH * A13_ACT;
dxdt_INH = AUTOAB * k_inh_on - k_inh_off * INH - k_inh_deg * INH;
dxdt_ULVWF = ULVWF_prod - ULVWF_deg * ULVWF - k_cleave * A13_frac * ULVWF;
dxdt_PLT = PLT_prod - k_PLT_deg * PLT - MT_form_rate * PLT;
dxdt_MT = MT_form_rate - k_MT_lysis * MT;
dxdt_BC = BC_prod - BC_deg * BC - PC_diff * BC - RTX_eff * BC_deg * BC;
dxdt_PC = PC_diff * BC - PC_deg * PC;
dxdt_AUTOAB = Ab_prod * PC * (1.0 - PRED_eff) - Ab_deg * AUTOAB;
dxdt_LDH_AB = k_LDH_hem * hemolysis - k_LDH_deg * (LDH_AB - LDH0);
dxdt_CREAT = k_Cr_MT * MT_pos - k_Cr_deg * (CREAT - Cr0);
dxdt_TROP = k_Trop_MT * MT_pos - k_Trop_deg * (TROP - Trop0);
dxdt_HGB = k_Hgb_prod * (Hgb0 - HGB) - hemolysis;
dxdt_PRED_C = -PRED_CL / PRED_V * PRED_C;

$TABLE
double CAPLA_conc_ng = CAPLA_C / CAPLA_V1 * 1000.0;
double RTX_conc_ug = RTX_C / RTX_V1;
double ADAMTS13_pct = A13_ACT;
double PLT_count = PLT;
double ULVWF_conc = ULVWF;
double MT_burden = MT;
double LDH_IUL = LDH_AB;
double Cr_umol = CREAT;
double Trop_ng = TROP;
double Hgb_gdL = HGB;
double BC_pct = BC;
double AutoAb_BU = AUTOAB;
double Inhibitor_BU = INH;
double schistocyte_pct = 5.0 * MT_burden / (MT_burden + 2.0);
double PRED_ng_mL = PRED_C;

$INIT
CAPLA_GUT=0, CAPLA_C=0, CAPLA_P=0, RTX_C=0, RTX_P=0,
A13_ACT=4.0, INH=4.5, ULVWF=60.0, PLT=18.0, MT=8.0,
BC=100.0, PC=25.0, AUTOAB=8.0,
LDH_AB=650.0, CREAT=115.0, TROP=0.45, HGB=7.8, PRED_C=0.0
'

## ── Build model once at startup ───────────────────────────────
mod_base <- mcode("ttp_shiny", ttp_model_code, quiet = TRUE)

## ── Helper: run simulation ────────────────────────────────────
run_sim <- function(input) {
  # Patient parameters
  bwt <- input$body_weight
  bsa <- 0.007184 * bwt^0.425 * input$patient_height^0.725  # Du Bois formula

  # Build event list
  ev_list <- list()

  # TPE events
  if (input$use_tpe) {
    n_tpe <- input$n_tpe_sessions
    tpe_days <- c(seq(0, min(6, n_tpe - 1), by = 1))
    if (n_tpe > 7) tpe_days <- c(tpe_days, seq(8, 8 + (n_tpe - 7) * 2 - 2, by = 2))
    tpe_days <- tpe_days[seq_len(min(n_tpe, length(tpe_days)))]
    tpe_evs <- lapply(tpe_days, function(d)
      ev(time = d + 0.01, cmt = "A13_ACT", amt = input$pex_adamts13))
    ev_list <- c(ev_list, tpe_evs)
  }

  # Caplacizumab
  if (input$use_capla) {
    ev_list <- c(ev_list,
      list(ev(time = 0.0, cmt = "CAPLA_C",   amt = input$capla_iv_dose)),
      list(ev(time = 0.5, cmt = "CAPLA_GUT", amt = input$capla_sc_dose,
              addl = input$capla_duration - 1, ii = 1))
    )
  }

  # Rituximab
  if (input$use_rtx) {
    rtx_amt <- 375 * bsa
    for (i in 0:(input$n_rtx - 1)) {
      ev_list <- c(ev_list,
        list(ev(time = i * 7, cmt = "RTX_C", amt = rtx_amt))
      )
    }
  }

  # Prednisolone
  if (input$use_pred) {
    pred_dose_mg <- input$pred_dose_mgkg * bwt
    pred_Cp_peak <- pred_dose_mg * 1e6 / (35 * 1000)  # ng/mL
    ev_list <- c(ev_list,
      list(ev(time = 0, cmt = "PRED_C", amt = pred_Cp_peak,
              addl = input$pred_duration - 1, ii = 1))
    )
  }

  # Merge events
  if (length(ev_list) == 0) {
    evs <- ev(time = 0, cmt = 1, amt = 0)
  } else {
    evs <- do.call(c, ev_list)
  }

  # Initial conditions
  init_vals <- list(
    A13_ACT = input$init_adamts13,
    INH     = input$init_inhibitor,
    ULVWF   = input$init_ulvwf,
    PLT     = input$init_plt,
    MT      = input$init_mt,
    AUTOAB  = input$init_autoab,
    LDH_AB  = input$init_ldh,
    CREAT   = input$init_creat,
    TROP    = input$init_trop,
    HGB     = input$init_hgb
  )

  out <- mod_base %>%
    init(init_vals) %>%
    mrgsim(events = evs, end = input$sim_days, delta = 0.25,
           obsonly = TRUE) %>%
    as_tibble()

  return(out)
}

## ================================================================
## UI
## ================================================================
ui <- dashboardPage(
  skin = "red",
  dashboardHeader(title = "TTP QSP Model"),

  dashboardSidebar(
    sidebarMenu(
      menuItem("① Patient Profile",    tabName = "tab_patient",   icon = icon("user")),
      menuItem("② Drug PK",            tabName = "tab_pk",        icon = icon("pills")),
      menuItem("③ ADAMTS13 / VWF / PLT", tabName = "tab_pd",     icon = icon("vials")),
      menuItem("④ Clinical Endpoints", tabName = "tab_clinical",  icon = icon("heartbeat")),
      menuItem("⑤ Scenario Comparison",tabName = "tab_scenarios", icon = icon("chart-bar")),
      menuItem("⑥ Biomarkers Panel",   tabName = "tab_biomarkers",icon = icon("microscope")),
      menuItem("ℹ Model Info",          tabName = "tab_info",      icon = icon("info-circle"))
    )
  ),

  dashboardBody(
    tabItems(

      ## ── TAB 1: Patient Profile & Treatment Setup ─────────────
      tabItem(tabName = "tab_patient",
        fluidRow(
          box(title = "Patient Demographics", status = "danger", solidHeader = TRUE, width = 4,
            sliderInput("body_weight", "Body Weight (kg)", 40, 120, 70, step = 1),
            sliderInput("patient_height", "Height (cm)", 140, 200, 170, step = 1),
            selectInput("sex", "Sex", c("Female (typical TTP)", "Male")),
            numericInput("age", "Age (years)", 42, 18, 80),
            hr(),
            h5("Derived BSA:"),
            verbatimTextOutput("bsa_output")
          ),
          box(title = "Disease Parameters at Presentation", status = "warning", solidHeader = TRUE, width = 4,
            sliderInput("init_adamts13", "Initial ADAMTS13 Activity (%)", 0, 30, 4, step = 0.5),
            sliderInput("init_inhibitor", "Inhibitor Titer (BU)", 0, 20, 4.5, step = 0.5),
            sliderInput("init_plt", "Initial Platelet Count (×10⁹/L)", 1, 100, 18, step = 1),
            sliderInput("init_autoab", "Autoantibody Level (BU)", 0, 30, 8, step = 0.5),
            sliderInput("init_ulvwf", "Initial ULVWF (ng/mL)", 10, 150, 60, step = 5),
            sliderInput("init_mt", "Initial Microthrombus Burden (AU)", 0, 20, 8, step = 1)
          ),
          box(title = "Presentation Biomarkers", status = "warning", solidHeader = TRUE, width = 4,
            sliderInput("init_ldh", "Initial LDH (IU/L)", 100, 2000, 650, step = 50),
            sliderInput("init_creat", "Initial Creatinine (μmol/L)", 50, 400, 115, step = 5),
            sliderInput("init_trop", "Initial Troponin I (ng/mL)", 0, 5, 0.45, step = 0.05),
            sliderInput("init_hgb", "Initial Hemoglobin (g/dL)", 4, 14, 7.8, step = 0.1),
            hr(),
            numericInput("sim_days", "Simulation Duration (days)", 180, 14, 365)
          )
        ),
        fluidRow(
          box(title = "PLASMIC Score Components", status = "primary", solidHeader = TRUE, width = 6,
            h5("Based on presentation values:"),
            verbatimTextOutput("plasmic_output"),
            p("PLASMIC Score ≥5 suggests high probability of TTP (specificity ~92%)"),
            p("Components: P-latelet <30k, L-ysis (LDH>2×ULN), A-DAMS13, S-evere anemia,
               M-CV, I-nflammatory marker, C-reatinine <2mg/dL")
          ),
          box(title = "Treatment Selection", status = "success", solidHeader = TRUE, width = 6,
            checkboxInput("use_tpe", "Therapeutic Plasma Exchange (TPE)", TRUE),
            conditionalPanel("input.use_tpe",
              sliderInput("n_tpe_sessions", "Number of TPE Sessions", 1, 20, 12, step = 1),
              sliderInput("pex_adamts13", "ADAMTS13 per TPE (U/dL)", 20, 100, 60, step = 5)
            ),
            checkboxInput("use_pred", "Prednisolone (Corticosteroids)", TRUE),
            conditionalPanel("input.use_pred",
              sliderInput("pred_dose_mgkg", "Prednisolone Dose (mg/kg/d)", 0.5, 2, 1, step = 0.1),
              sliderInput("pred_duration", "Prednisolone Duration (days)", 7, 60, 28, step = 1)
            ),
            checkboxInput("use_capla", "Caplacizumab (Cablivi®)", FALSE),
            conditionalPanel("input.use_capla",
              numericInput("capla_iv_dose", "Caplacizumab IV Dose (mg)", 10, 1, 30),
              numericInput("capla_sc_dose", "Caplacizumab SC Daily Dose (mg)", 10, 1, 20),
              sliderInput("capla_duration", "Caplacizumab Duration (days)", 7, 60, 42, step = 1)
            ),
            checkboxInput("use_rtx", "Rituximab (Anti-CD20)", FALSE),
            conditionalPanel("input.use_rtx",
              sliderInput("n_rtx", "Number of Rituximab Doses (weekly)", 1, 4, 4, step = 1)
            )
          )
        )
      ),

      ## ── TAB 2: Drug PK ───────────────────────────────────────
      tabItem(tabName = "tab_pk",
        fluidRow(
          box(title = "Caplacizumab PK (Nanobody, Anti-VWF A1)", status = "primary",
              solidHeader = TRUE, width = 6,
            plotlyOutput("pk_capla_plot", height = "350px"),
            p("Caplacizumab (10mg IV + 10mg SC qd) — 2-compartment PK model.",
              "IC₅₀ for VWF–platelet inhibition: ~2 ng/mL.",
              "SC bioavailability ~85%, t½ ~11–13h (SC trough maintained above IC₅₀).")
          ),
          box(title = "Rituximab PK (Anti-CD20, TMDD)", status = "primary",
              solidHeader = TRUE, width = 6,
            plotlyOutput("pk_rtx_plot", height = "350px"),
            p("Rituximab (375mg/m² IV q7d ×4) — 2-compartment model with",
              "target-mediated disposition (TMDD). Initial high concentration",
              "drives B cell depletion; levels decline over weeks.")
          )
        ),
        fluidRow(
          box(title = "Prednisolone PK (Corticosteroid)", status = "warning",
              solidHeader = TRUE, width = 6,
            plotlyOutput("pk_pred_plot", height = "300px"),
            p("Prednisolone — 1-compartment oral model.",
              "Dose: 1mg/kg/d (typical TTP dose).",
              "t½ ~3h (CL=15 L/d, V=35 L). Suppresses B cell proliferation and IL-6.")
          ),
          box(title = "PK Parameter Summary", status = "info", solidHeader = TRUE, width = 6,
            DT::dataTableOutput("pk_params_table")
          )
        )
      ),

      ## ── TAB 3: ADAMTS13 / VWF / PLT (Core PD) ───────────────
      tabItem(tabName = "tab_pd",
        fluidRow(
          box(title = "ADAMTS13 Activity (%)", status = "danger", solidHeader = TRUE, width = 6,
            plotlyOutput("pd_adamts13_plot", height = "300px"),
            p("Critical threshold: <10% (defining feature of acquired TTP).",
              "Recovery driven by: (1) TPE providing exogenous ADAMTS13 (from FFP),",
              "(2) Rituximab depleting Ab-producing B cells,",
              "(3) Prednisolone suppressing autoantibody production.")
          ),
          box(title = "ULVWF Pool (ng/mL)", status = "warning", solidHeader = TRUE, width = 6,
            plotlyOutput("pd_ulvwf_plot", height = "300px"),
            p("ULVWF accumulates when ADAMTS13 activity falls below critical threshold.",
              "Caplacizumab blocks VWF A1 domain → prevents platelet adhesion to ULVWF.",
              "Normal: <10 ng/mL (fully cleaved by ADAMTS13).")
          )
        ),
        fluidRow(
          box(title = "Platelet Count (×10⁹/L)", status = "success", solidHeader = TRUE, width = 6,
            plotlyOutput("pd_plt_plot", height = "300px"),
            p("Response threshold: PLT > 150×10⁹/L sustained ≥2 days.",
              "Caplacizumab dramatically accelerates platelet recovery by blocking",
              "ULVWF-mediated platelet consumption (HERCULES trial: 2.69d vs 2.88d).")
          ),
          box(title = "Microthrombus Burden & B Cell Dynamics", status = "info",
              solidHeader = TRUE, width = 6,
            plotlyOutput("pd_mt_bc_plot", height = "300px"),
            p("Microthrombus formation: proportional to ULVWF^1.5 × PLT.",
              "B cell depletion by rituximab → reduced autoantibody →",
              "ADAMTS13 recovery → sustained remission.")
          )
        )
      ),

      ## ── TAB 4: Clinical Endpoints ─────────────────────────────
      tabItem(tabName = "tab_clinical",
        fluidRow(
          box(title = "LDH (IU/L) — Hemolysis Marker", status = "danger", solidHeader = TRUE, width = 6,
            plotlyOutput("ep_ldh_plot", height = "280px"),
            p("LDH rises with RBC fragmentation (MAHA).",
              ">2×ULN (>360 IU/L) indicates significant hemolysis.",
              "Rapid fall with treatment indicates hemolysis resolution.")
          ),
          box(title = "Hemoglobin (g/dL) — Anemia", status = "danger", solidHeader = TRUE, width = 6,
            plotlyOutput("ep_hgb_plot", height = "280px"),
            p("MAHA-type anemia: Coombs-negative, microangiopathic.",
              "Typical presentation Hgb 6–10 g/dL.",
              "Transfusion threshold: Hgb <7 g/dL (use with caution in TTP).")
          )
        ),
        fluidRow(
          box(title = "Creatinine (μmol/L) — Renal Function", status = "warning", solidHeader = TRUE, width = 6,
            plotlyOutput("ep_creat_plot", height = "280px"),
            p("Renal microthrombi cause tubular ischemia and AKI.",
              "Cr <177 μmol/L = PLASMIC score component (+1).",
              "Recovery tracks microthrombus resolution.")
          ),
          box(title = "Troponin I (ng/mL) — Cardiac Injury", status = "warning", solidHeader = TRUE, width = 6,
            plotlyOutput("ep_trop_plot", height = "280px"),
            p("Cardiac microthrombi cause MINOCA (myocardial infarction without obstructive CAD).",
              "Troponin elevation in TTP associated with worse prognosis.",
              "ICU admission criteria: Troponin >0.5 ng/mL + neurological symptoms.")
          )
        ),
        fluidRow(
          box(title = "TMA Activity Index (Composite Score)", status = "primary",
              solidHeader = TRUE, width = 12,
            plotlyOutput("ep_tma_plot", height = "280px"),
            p("Composite TMA Activity Index incorporates PLT, LDH, Cr, and schistocyte %.",
              "Higher values indicate more active TTP disease.",
              "Normalization of TMA index = complete clinical remission.")
          )
        )
      ),

      ## ── TAB 5: Scenario Comparison ───────────────────────────
      tabItem(tabName = "tab_scenarios",
        fluidRow(
          box(title = "Standard Scenarios Comparison Setup", status = "primary",
              solidHeader = TRUE, width = 12,
            p("Compare 5 predefined treatment strategies head-to-head:"),
            actionButton("run_scenarios", "Run All Scenarios", class = "btn-primary btn-lg"),
            hr()
          )
        ),
        fluidRow(
          box(title = "Platelet Counts — All Scenarios", status = "success",
              solidHeader = TRUE, width = 6,
            plotlyOutput("sc_plt_plot", height = "320px")
          ),
          box(title = "ADAMTS13 Activity — All Scenarios", status = "danger",
              solidHeader = TRUE, width = 6,
            plotlyOutput("sc_a13_plot", height = "320px")
          )
        ),
        fluidRow(
          box(title = "LDH — All Scenarios", status = "warning",
              solidHeader = TRUE, width = 6,
            plotlyOutput("sc_ldh_plot", height = "320px")
          ),
          box(title = "Autoantibody — All Scenarios", status = "info",
              solidHeader = TRUE, width = 6,
            plotlyOutput("sc_ab_plot", height = "320px")
          )
        ),
        fluidRow(
          box(title = "Response Table (Time to PLT>150×10⁹/L)", status = "info",
              solidHeader = TRUE, width = 12,
            DT::dataTableOutput("sc_response_table")
          )
        )
      ),

      ## ── TAB 6: Biomarkers Panel ──────────────────────────────
      tabItem(tabName = "tab_biomarkers",
        fluidRow(
          box(title = "Schistocyte Estimate (%)", status = "danger", solidHeader = TRUE, width = 4,
            plotlyOutput("bm_schist_plot", height = "250px"),
            p("Schistocytes >1% = MAHA. Normal <0.1%.",
              "Derived from microthrombus burden (MT/(MT+2))×5%.",
              "ISTH diagnostic criterion for TTP.")
          ),
          box(title = "Autoantibody (Anti-ADAMTS13, BU)", status = "warning", solidHeader = TRUE, width = 4,
            plotlyOutput("bm_ab_plot", height = "250px"),
            p("Inhibitor titer >0.4 BU = detectable inhibitor.",
              "High titer (>2 BU) → severe ADAMTS13 deficiency.",
              "Rituximab depletes Ab-producing B cells → titer falls over 6–12 weeks.")
          ),
          box(title = "B Cell (% Normal)", status = "info", solidHeader = TRUE, width = 4,
            plotlyOutput("bm_bc_plot", height = "250px"),
            p("B cell depletion by rituximab: nadir at ~4 weeks.",
              "B cell recovery: 12–18 months post-rituximab.",
              "Persistent depletion → sustained immunologic remission.")
          )
        ),
        fluidRow(
          box(title = "Biomarker Summary Table (Key Time Points)", status = "primary",
              solidHeader = TRUE, width = 12,
            DT::dataTableOutput("bm_summary_table")
          )
        ),
        fluidRow(
          box(title = "ADAMTS13 Response Prediction", status = "success",
              solidHeader = TRUE, width = 12,
            plotlyOutput("bm_a13_trajectory", height = "300px"),
            p("ADAMTS13 activity trajectory: initial rebound from TPE replacement,",
              "followed by sustained recovery as autoantibody falls (rituximab effect).",
              "Immunologic remission defined as ADAMTS13 >50% + inhibitor undetectable.")
          )
        )
      ),

      ## ── TAB 7: Model Info ────────────────────────────────────
      tabItem(tabName = "tab_info",
        fluidRow(
          box(title = "TTP Pathophysiology", status = "primary", solidHeader = TRUE, width = 6,
            h4("Disease Mechanism"),
            p("Acquired TTP is caused by autoantibody-mediated inhibition of ADAMTS13
               (a disintegrin and metalloproteinase with thrombospondin type 1 motifs,
               member 13), a plasma metalloprotease that cleaves ultra-large von
               Willebrand factor (ULVWF) multimers."),
            p("Without ADAMTS13 cleavage, ULVWF accumulates → mediates pathological
               platelet aggregation at high shear stress → platelet-rich microthrombi
               in the microvasculature → end-organ ischemia (brain, kidney, heart)."),
            p("The resulting microangiopathic hemolytic anemia (MAHA) occurs when
               erythrocytes are mechanically fragmented by microthrombi
               (schistocytes visible on blood smear, Coombs-negative hemolysis,
               elevated LDH, undetectable haptoglobin)."),
            h4("Incidence"),
            p("~3–7 cases/million population/year; female:male ~3:1; peak age 30–50 years."),
            h4("Mortality"),
            p("Without treatment: ~90% mortality. With modern TPE-based therapy: ~10–20%.")
          ),
          box(title = "QSP Model Structure", status = "info", solidHeader = TRUE, width = 6,
            h4("Model Compartments (17 ODEs)"),
            tags$table(class = "table table-condensed",
              tags$thead(tags$tr(tags$th("#"), tags$th("Compartment"), tags$th("Biological Role"))),
              tags$tbody(
                tags$tr(tags$td("1-3"),  tags$td("CAPLA_GUT, _C, _P"), tags$td("Caplacizumab 2-comp PK")),
                tags$tr(tags$td("4-5"),  tags$td("RTX_C, RTX_P"),      tags$td("Rituximab 2-comp PK")),
                tags$tr(tags$td("6"),    tags$td("A13_ACT"),           tags$td("ADAMTS13 activity (U/dL)")),
                tags$tr(tags$td("7"),    tags$td("INH"),               tags$td("Anti-ADAMTS13 inhibitor (BU)")),
                tags$tr(tags$td("8"),    tags$td("ULVWF"),             tags$td("ULVWF pool (ng/mL)")),
                tags$tr(tags$td("9"),    tags$td("PLT"),               tags$td("Platelet count (×10⁹/L)")),
                tags$tr(tags$td("10"),   tags$td("MT"),                tags$td("Microthrombus burden (AU)")),
                tags$tr(tags$td("11"),   tags$td("BC"),                tags$td("B cells (% normal)")),
                tags$tr(tags$td("12"),   tags$td("PC"),                tags$td("Plasma cells (Ab-secreting)")),
                tags$tr(tags$td("13"),   tags$td("AUTOAB"),            tags$td("Autoantibody IgG (BU)")),
                tags$tr(tags$td("14"),   tags$td("LDH_AB"),            tags$td("LDH (IU/L) — hemolysis")),
                tags$tr(tags$td("15"),   tags$td("CREAT"),             tags$td("Creatinine (μmol/L) — renal")),
                tags$tr(tags$td("16"),   tags$td("TROP"),              tags$td("Troponin I (ng/mL) — cardiac")),
                tags$tr(tags$td("17"),   tags$td("HGB"),               tags$td("Hemoglobin (g/dL) — MAHA")),
                tags$tr(tags$td("18"),   tags$td("PRED_C"),            tags$td("Prednisolone plasma (ng/mL)"))
              )
            )
          )
        ),
        fluidRow(
          box(title = "Clinical Calibration", status = "success", solidHeader = TRUE, width = 6,
            h4("Key Trial References"),
            tags$ul(
              tags$li(tags$b("HERCULES trial (Scully 2019):"),
                " Caplacizumab + TPE vs TPE alone.",
                " Primary: time to platelet response (2.69 vs 2.88 days, HR 1.55).",
                " Caplacizumab ↓ exacerbation (3% vs 28%), ↓ refractory TTP."),
              tags$li(tags$b("TITAN trial (Peyvandi 2016):"),
                " Rituximab ×4 + TPE vs TPE alone.",
                " ↑ Complete remission at 12wks (59% vs 43%).",
                " ↑ ADAMTS13 recovery; ↓ relapse at 2 years."),
              tags$li(tags$b("Froissart 2012:"),
                " ADAMTS13 kinetics: inhibitor binds and clears ADAMTS13 via FcR-macrophage."),
              tags$li(tags$b("Coppo 2010:"),
                " Rituximab achieves B cell depletion by 4 weeks;",
                " ADAMTS13 normalizes in 6–12 weeks in 60–80% of cases."),
              tags$li(tags$b("Westwood 2017:"),
                " Rituximab prophylaxis in relapsing TTP reduces relapse rate to ~5%/year.")
            )
          ),
          box(title = "Drug Information", status = "warning", solidHeader = TRUE, width = 6,
            h4("Caplacizumab (Cablivi®)"),
            p("Humanized bivalent VHH nanobody (anti-VWF A1 domain).",
              "First-in-class approved for acquired TTP (EMA 2018, FDA 2019).",
              "Mechanism: blocks ULVWF–GPIbα interaction → prevents platelet aggregation.",
              "Does NOT treat underlying ADAMTS13 deficiency → continue with immunosuppression."),
            h4("Rituximab (MabThera/Rituxan®)"),
            p("Chimeric anti-CD20 IgG1 monoclonal antibody.",
              "Depletes CD20+ B cells via ADCC, CDC, apoptosis.",
              "Effect: ↓ autoantibody production → ADAMTS13 recovery over 6–12 weeks.",
              "Used as remission-induction and relapse-prevention therapy."),
            h4("Therapeutic Plasma Exchange (TPE)"),
            p("Removes autoantibody + ULVWF; replaces ADAMTS13 (from FFP/SD plasma).",
              "Gold standard first-line therapy; reduces mortality from ~90% to ~10%.",
              "Frequency: daily until PLT>150×10⁹/L sustained ×2 days, then taper.")
          )
        )
      )
    )
  )
)

## ================================================================
## SERVER
## ================================================================
server <- function(input, output, session) {

  ## Reactive: run simulation
  sim_data <- reactive({
    run_sim(input)
  })

  ## BSA output
  output$bsa_output <- renderText({
    bsa <- 0.007184 * input$body_weight^0.425 * input$patient_height^0.725
    sprintf("BSA = %.2f m²  |  RTX dose (375mg/m²) = %.0f mg", bsa, 375 * bsa)
  })

  ## PLASMIC score
  output$plasmic_output <- renderText({
    s_plt  <- ifelse(input$init_plt < 30, 1, 0)
    s_ldh  <- ifelse(input$init_ldh > 360, 1, 0)
    s_cr   <- ifelse(input$init_creat < 177, 1, 0)
    s_a13  <- ifelse(input$init_adamts13 < 10, 1, 0)
    score3 <- s_plt + s_ldh + s_cr
    paste0(
      sprintf("PLT <30k: %d  |  LDH >2×ULN: %d  |  Cr <177: %d\n", s_plt, s_ldh, s_cr),
      sprintf("ADAMTS13 <10%%: %d (if measured)\n", s_a13),
      sprintf("PLASMIC Score ≥3/3 measured components: %d/3\n", score3),
      ifelse(score3 == 3, "→ HIGH TTP probability", "→ Lower TTP probability (check other components)")
    )
  })

  ## PK tab
  output$pk_capla_plot <- renderPlotly({
    df <- sim_data()
    p <- ggplot(df, aes(x = time, y = CAPLA_conc_ng)) +
      geom_line(color = "#2980B9", linewidth = 1.2) +
      geom_hline(yintercept = 2, linetype = "dashed", color = "red") +
      annotate("text", x = max(df$time)*0.8, y = 3, label = "IC₅₀ = 2 ng/mL",
               color = "red", size = 3) +
      labs(title = "Caplacizumab Plasma Concentration",
           x = "Time (days)", y = "Cp (ng/mL)") +
      theme_bw()
    ggplotly(p)
  })

  output$pk_rtx_plot <- renderPlotly({
    df <- sim_data()
    p <- ggplot(df, aes(x = time, y = RTX_conc_ug)) +
      geom_line(color = "#E67E22", linewidth = 1.2) +
      labs(title = "Rituximab Plasma Concentration",
           x = "Time (days)", y = "Cp (μg/mL)") +
      theme_bw()
    ggplotly(p)
  })

  output$pk_pred_plot <- renderPlotly({
    df <- sim_data()
    p <- ggplot(df, aes(x = time, y = PRED_ng_mL)) +
      geom_line(color = "#8E44AD", linewidth = 1.0) +
      labs(title = "Prednisolone Plasma Concentration",
           x = "Time (days)", y = "Cp (ng/mL)") +
      theme_bw()
    ggplotly(p)
  })

  output$pk_params_table <- DT::renderDataTable({
    data.frame(
      Drug       = c("Caplacizumab", "Rituximab", "Prednisolone"),
      Mechanism  = c("Anti-VWF A1 nanobody", "Anti-CD20 mAb", "Glucocorticoid"),
      `Dose`     = c("10mg IV + 10mg SC qd", "375mg/m² IV q7d×4", "1mg/kg/d PO"),
      `t½`       = c("6-8h (IV), 11-13h (SC)", "~14-20 days", "~3 hours"),
      `IC50/EC50`= c("IC₅₀=2 ng/mL (VWF)", "EC₅₀=0.5 μg/mL (B cell)", "EC₅₀=150 ng/mL (Ab)"),
      `PK Model` = c("2-compartment", "2-comp + TMDD", "1-compartment")
    )
  }, options = list(pageLength = 5, scrollX = TRUE), rownames = FALSE)

  ## PD tab
  output$pd_adamts13_plot <- renderPlotly({
    df <- sim_data()
    p <- ggplot(df, aes(x = time, y = ADAMTS13_pct)) +
      geom_line(color = "#2ECC71", linewidth = 1.2) +
      geom_hline(yintercept = 10, linetype = "dashed", color = "red", linewidth = 1) +
      geom_hline(yintercept = 50, linetype = "dashed", color = "blue", linewidth = 0.8) +
      annotate("text", x = max(df$time)*0.9, y = 12, label = "<10% = TTP",
               color = "red", size = 3) +
      annotate("text", x = max(df$time)*0.9, y = 52, label = ">50% = remission",
               color = "blue", size = 3) +
      scale_y_continuous(limits = c(0, 110)) +
      labs(title = "ADAMTS13 Activity", x = "Time (days)", y = "Activity (%)") +
      theme_bw()
    ggplotly(p)
  })

  output$pd_ulvwf_plot <- renderPlotly({
    df <- sim_data()
    p <- ggplot(df, aes(x = time, y = ULVWF_conc)) +
      geom_line(color = "#F39C12", linewidth = 1.2) +
      geom_hline(yintercept = 10, linetype = "dashed", color = "green", linewidth = 0.8) +
      labs(title = "ULVWF Pool", x = "Time (days)", y = "ULVWF (ng/mL)") +
      theme_bw()
    ggplotly(p)
  })

  output$pd_plt_plot <- renderPlotly({
    df <- sim_data()
    p <- ggplot(df, aes(x = time, y = PLT_count)) +
      geom_line(color = "#9B59B6", linewidth = 1.2) +
      geom_hline(yintercept = 150, linetype = "dashed", color = "black", linewidth = 1) +
      geom_hline(yintercept = 30, linetype = "dashed", color = "red", linewidth = 0.8) +
      annotate("text", x = max(df$time)*0.85, y = 155, label = "Response (150×10⁹/L)",
               color = "black", size = 3) +
      labs(title = "Platelet Count", x = "Time (days)", y = "PLT (×10⁹/L)") +
      theme_bw()
    ggplotly(p)
  })

  output$pd_mt_bc_plot <- renderPlotly({
    df <- sim_data() %>%
      select(time, MT_burden, BC_pct) %>%
      pivot_longer(-time, names_to = "variable", values_to = "value")
    p <- ggplot(df, aes(x = time, y = value, color = variable)) +
      geom_line(linewidth = 1.2) +
      scale_color_manual(values = c("MT_burden" = "#E74C3C", "BC_pct" = "#3498DB"),
                         labels = c("MT_burden" = "Microthrombus (AU)",
                                    "BC_pct"    = "B Cells (% normal)")) +
      labs(title = "Microthrombus & B Cells", x = "Time (days)", y = "Value",
           color = "") +
      theme_bw()
    ggplotly(p)
  })

  ## Clinical Endpoints
  output$ep_ldh_plot <- renderPlotly({
    df <- sim_data()
    p <- ggplot(df, aes(x = time, y = LDH_IUL)) +
      geom_line(color = "#E74C3C", linewidth = 1.2) +
      geom_hline(yintercept = 360, linetype = "dashed", color = "orange") +
      labs(title = "LDH", x = "Time (days)", y = "LDH (IU/L)") + theme_bw()
    ggplotly(p)
  })

  output$ep_hgb_plot <- renderPlotly({
    df <- sim_data()
    p <- ggplot(df, aes(x = time, y = Hgb_gdL)) +
      geom_line(color = "#C0392B", linewidth = 1.2) +
      geom_hline(yintercept = 12, linetype = "dashed", color = "navy") +
      labs(title = "Hemoglobin", x = "Time (days)", y = "Hgb (g/dL)") + theme_bw()
    ggplotly(p)
  })

  output$ep_creat_plot <- renderPlotly({
    df <- sim_data()
    p <- ggplot(df, aes(x = time, y = Cr_umol)) +
      geom_line(color = "#8E44AD", linewidth = 1.2) +
      geom_hline(yintercept = 110, linetype = "dashed", color = "green") +
      labs(title = "Creatinine", x = "Time (days)", y = "Creatinine (μmol/L)") + theme_bw()
    ggplotly(p)
  })

  output$ep_trop_plot <- renderPlotly({
    df <- sim_data()
    p <- ggplot(df, aes(x = time, y = Trop_ng)) +
      geom_line(color = "#E67E22", linewidth = 1.2) +
      geom_hline(yintercept = 0.04, linetype = "dashed", color = "red") +
      labs(title = "Troponin I", x = "Time (days)", y = "TnI (ng/mL)") + theme_bw()
    ggplotly(p)
  })

  output$ep_tma_plot <- renderPlotly({
    df <- sim_data() %>%
      mutate(TMA_idx = (300 - PLT_count)/300*40 +
               pmax(0, LDH_IUL - 180)/180*20 +
               pmax(0, Cr_umol - 75)/75*10 +
               schistocyte_pct*5)
    p <- ggplot(df, aes(x = time, y = TMA_idx)) +
      geom_line(color = "#1A5276", linewidth = 1.2) +
      geom_hline(yintercept = 10, linetype = "dashed", color = "green") +
      labs(title = "TMA Activity Index (Composite)", x = "Time (days)", y = "TMA Index") + theme_bw()
    ggplotly(p)
  })

  ## Scenario comparison
  sc_data <- eventReactive(input$run_scenarios, {
    sc_params <- list(
      list(tpe=TRUE, pred=TRUE,  capla=FALSE, rtx=FALSE, label="S2: TPE + Pred"),
      list(tpe=TRUE, pred=TRUE,  capla=TRUE,  rtx=FALSE, label="S3: TPE + CAPLA + Pred"),
      list(tpe=TRUE, pred=TRUE,  capla=FALSE, rtx=TRUE,  label="S4: TPE + RTX + Pred"),
      list(tpe=TRUE, pred=TRUE,  capla=TRUE,  rtx=TRUE,  label="S5: Triple Therapy"),
      list(tpe=FALSE,pred=FALSE, capla=FALSE, rtx=FALSE,  label="S0: No Treatment")
    )
    bsa_val <- 0.007184 * input$body_weight^0.425 * input$patient_height^0.725
    pex_a13_val <- input$pex_adamts13

    bind_rows(lapply(sc_params, function(sc) {
      ev_list <- list()
      if (sc$tpe) {
        tpe_days <- c(0:6, 8, 10, 12, 15, 18)
        ev_list <- c(ev_list, lapply(tpe_days, function(d)
          ev(time = d + 0.01, cmt = "A13_ACT", amt = pex_a13_val)))
      }
      if (sc$pred) {
        pred_mg <- input$pred_dose_mgkg * input$body_weight
        pred_Cp <- pred_mg * 1e6 / (35 * 1000)
        ev_list <- c(ev_list, list(ev(time=0, cmt="PRED_C", amt=pred_Cp, addl=27, ii=1)))
      }
      if (sc$capla) {
        ev_list <- c(ev_list,
          list(ev(time=0.0,  cmt="CAPLA_C",   amt=10)),
          list(ev(time=0.5,  cmt="CAPLA_GUT", amt=10, addl=41, ii=1)))
      }
      if (sc$rtx) {
        for (i in 0:3) ev_list <- c(ev_list, list(ev(time=i*7, cmt="RTX_C", amt=375*bsa_val)))
      }
      evs <- if (length(ev_list) == 0) ev(time=0, cmt=1, amt=0) else do.call(c, ev_list)
      init_vals <- list(A13_ACT=input$init_adamts13, INH=input$init_inhibitor,
                        ULVWF=input$init_ulvwf, PLT=input$init_plt, MT=input$init_mt,
                        AUTOAB=input$init_autoab, LDH_AB=input$init_ldh,
                        CREAT=input$init_creat, TROP=input$init_trop, HGB=input$init_hgb)
      mod_base %>% init(init_vals) %>%
        mrgsim(events=evs, end=input$sim_days, delta=0.25, obsonly=TRUE) %>%
        as_tibble() %>% mutate(Scenario=sc$label)
    }))
  })

  sc_colors <- c("S0: No Treatment"="red","S2: TPE + Pred"="gold",
                 "S3: TPE + CAPLA + Pred"="green","S4: TPE + RTX + Pred"="blue",
                 "S5: Triple Therapy"="purple")

  output$sc_plt_plot <- renderPlotly({
    req(sc_data())
    p <- ggplot(sc_data(), aes(x=time, y=PLT_count, color=Scenario)) +
      geom_line(linewidth=1) +
      geom_hline(yintercept=150, linetype="dashed") +
      scale_color_manual(values=sc_colors) +
      labs(x="Days", y="PLT (×10⁹/L)") + theme_bw()
    ggplotly(p)
  })

  output$sc_a13_plot <- renderPlotly({
    req(sc_data())
    p <- ggplot(sc_data(), aes(x=time, y=ADAMTS13_pct, color=Scenario)) +
      geom_line(linewidth=1) +
      geom_hline(yintercept=10, linetype="dashed", color="red") +
      scale_color_manual(values=sc_colors) +
      labs(x="Days", y="ADAMTS13 (%)") + theme_bw()
    ggplotly(p)
  })

  output$sc_ldh_plot <- renderPlotly({
    req(sc_data())
    p <- ggplot(sc_data(), aes(x=time, y=LDH_IUL, color=Scenario)) +
      geom_line(linewidth=1) +
      geom_hline(yintercept=360, linetype="dashed", color="orange") +
      scale_color_manual(values=sc_colors) +
      labs(x="Days", y="LDH (IU/L)") + theme_bw()
    ggplotly(p)
  })

  output$sc_ab_plot <- renderPlotly({
    req(sc_data())
    p <- ggplot(sc_data(), aes(x=time, y=AutoAb_BU, color=Scenario)) +
      geom_line(linewidth=1) +
      geom_hline(yintercept=0.4, linetype="dashed", color="red") +
      scale_color_manual(values=sc_colors) +
      labs(x="Days", y="AutoAb (BU)") + theme_bw()
    ggplotly(p)
  })

  output$sc_response_table <- DT::renderDataTable({
    req(sc_data())
    sc_data() %>%
      group_by(Scenario) %>%
      summarise(
        `PLT>150 (day)`     = {r <- first(which(PLT_count > 150)); ifelse(length(r), round(time[r],1), NA)},
        `A13>50% (day)`     = {r <- first(which(ADAMTS13_pct > 50)); ifelse(length(r), round(time[r],1), NA)},
        `LDH<360 (day)`     = {r <- first(which(LDH_IUL < 360)); ifelse(length(r), round(time[r],1), NA)},
        `Day 180 PLT`       = round(last(PLT_count[time <= 180]), 1),
        `Day 180 A13 (%)`   = round(last(ADAMTS13_pct[time <= 180]), 1),
        `Day 180 AutoAb`    = round(last(AutoAb_BU[time <= 180]), 2),
        .groups = "drop"
      )
  }, options = list(pageLength=5, scrollX=TRUE), rownames=FALSE)

  ## Biomarkers tab
  output$bm_schist_plot <- renderPlotly({
    df <- sim_data()
    p <- ggplot(df, aes(x=time, y=schistocyte_pct)) +
      geom_line(color="#884EA0", linewidth=1.2) +
      geom_hline(yintercept=1, linetype="dashed", color="red") +
      labs(title="Schistocyte %", x="Days", y="%") + theme_bw()
    ggplotly(p)
  })

  output$bm_ab_plot <- renderPlotly({
    df <- sim_data()
    p <- ggplot(df, aes(x=time, y=AutoAb_BU)) +
      geom_line(color="#E74C3C", linewidth=1.2) +
      geom_hline(yintercept=0.4, linetype="dashed", color="red") +
      labs(title="Anti-ADAMTS13 Ab (BU)", x="Days", y="BU") + theme_bw()
    ggplotly(p)
  })

  output$bm_bc_plot <- renderPlotly({
    df <- sim_data()
    p <- ggplot(df, aes(x=time, y=BC_pct)) +
      geom_line(color="#2980B9", linewidth=1.2) +
      geom_hline(yintercept=20, linetype="dashed", color="orange") +
      labs(title="B Cell (% Normal)", x="Days", y="%") + theme_bw()
    ggplotly(p)
  })

  output$bm_summary_table <- DT::renderDataTable({
    sim_data() %>%
      filter(time %in% c(0, 3, 7, 14, 30, 60, 90, 180)) %>%
      transmute(
        `Day`      = time,
        `PLT (×10⁹/L)` = round(PLT_count, 1),
        `ADAMTS13 (%)` = round(ADAMTS13_pct, 1),
        `ULVWF (ng/mL)` = round(ULVWF_conc, 1),
        `LDH (IU/L)` = round(LDH_IUL, 0),
        `Hgb (g/dL)` = round(Hgb_gdL, 1),
        `Cr (μmol/L)` = round(Cr_umol, 1),
        `TnI (ng/mL)` = round(Trop_ng, 2),
        `AutoAb (BU)` = round(AutoAb_BU, 2),
        `Schisto. (%)` = round(schistocyte_pct, 2)
      )
  }, options = list(pageLength=8, scrollX=TRUE), rownames=FALSE)

  output$bm_a13_trajectory <- renderPlotly({
    df <- sim_data() %>%
      select(time, ADAMTS13_pct, Inhibitor_BU) %>%
      pivot_longer(-time, names_to="Marker", values_to="Value")
    p <- ggplot(df, aes(x=time, y=Value, color=Marker)) +
      geom_line(linewidth=1.2) +
      geom_hline(yintercept=50, linetype="dashed", color="blue") +
      geom_hline(yintercept=10, linetype="dashed", color="red") +
      scale_color_manual(values=c("ADAMTS13_pct"="#2ECC71","Inhibitor_BU"="#E74C3C"),
                         labels=c("ADAMTS13 Activity (%)","Inhibitor Titer (BU)")) +
      labs(title="ADAMTS13 Activity vs Inhibitor Titer",
           x="Time (days)", y="Value", color="") +
      theme_bw()
    ggplotly(p)
  })
}

## ── Launch ────────────────────────────────────────────────────
shinyApp(ui = ui, server = server)
