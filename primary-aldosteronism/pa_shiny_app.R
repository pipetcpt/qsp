# ============================================================
# Primary Aldosteronism (Conn's Syndrome) — QSP Shiny App
# ============================================================
# Tabs: 1) Patient Profile  2) RAAS / Drug PK  3) Aldosterone Panel
#       4) Ion Homeostasis  5) Cardiovascular & Organ Damage
#       6) Scenario Comparison  7) Biomarker Explorer
# ============================================================

library(shiny)
library(shinydashboard)
library(mrgsolve)
library(dplyr)
library(tidyr)
library(ggplot2)
library(plotly)
library(DT)

# ── inline mrgsolve model code ─────────────────────────────────────────────
pa_model_code <- '
$PARAM @annotated
// APA/BAH disease
APA_severity : 1.0 : APA severity (0=none, 1=typical APA)
BAH_severity : 0.0 : BAH severity (0=none, 1=bilateral hyperplasia)
surgery      : 0   : Adrenalectomy flag (0/1)
surgery_effect : 0.90 : Efficacy of adrenalectomy on APA

// Spironolactone PK
spiro_dose  : 0    : Spironolactone daily dose (mg)
ka_spiro    : 1.4  : Absorption rate constant (1/h)
Vd_spiro    : 80   : Volume of distribution (L)
CL_spiro    : 18   : Clearance (L/h)
ka_canr     : 0.06 : Canrenone formation (1/h)
CL_canr     : 2.8  : Canrenone clearance (L/h)
Vd_canr     : 60   : Canrenone Vd (L)

// Eplerenone PK
eple_dose   : 0
ka_eple     : 1.2
Vd_eple     : 130
CL_eple     : 21

// Finerenone PK
fine_dose   : 0
ka_fine     : 1.5
Vd_fine     : 52
CL_fine     : 9.6

// ACEi (Ramipril) PK
acei_dose   : 0
ka_acei     : 1.0
Vd_acei     : 500
CL_acei     : 14

// CCB (Amlodipine) PK
ccb_dose    : 0
ka_ccb      : 0.5
Vd_ccb      : 2100
CL_ccb      : 35

// RAAS physiology
Renin_base  : 1.0
Renin_kout  : 0.15
AngII_kout  : 0.8
Aldo_kout   : 0.12
Aldo_base   : 8.0

// Renal / ion physiology
ENaC_kout   : 0.5
Na_kout     : 0.05
K_kout      : 0.08
HCO3_kout   : 0.04
Vol_kout     : 0.03

// Cardiovascular
MAP_kout    : 0.02
TPR_kout    : 0.04
GFR_kout    : 0.01

// Target organ damage
CardFib_kin  : 0.0003
CardFib_kout : 0.0002
LVMi_kout    : 0.005

// PD parameters (MR antagonism)
IC50_spiro  : 1.2
IC50_canr   : 0.8
IC50_eple   : 2.5
IC50_fine   : 0.65
Emax_MRA    : 0.88

// ACEi effect on AngII
IC50_acei_angII : 0.3
Emax_acei       : 0.75

// CCB effect on TPR/MAP
IC50_ccb_tpr : 0.05
Emax_ccb     : 0.40

$INIT
C_spiro  = 0
C_canr   = 0
C_eple   = 0
C_fine   = 0
C_acei   = 0
C_ccb    = 0
Renin_c  = 1.0
AngII_c  = 1.0
Aldo_c   = 8.0
ENaC_act = 1.0
Na_c     = 140
K_c      = 4.0
HCO3_c   = 24.0
Vol_c    = 1.0
MAP_c    = 95.0
TPR_c    = 1.0
GFR_c    = 1.0
CardFib  = 0.0
LVMi_c   = 90.0
APA_act  = 1.0
CYP11B2_c = 1.0
ARR_c    = 30.0
HOMA_proxy = 1.0

$MAIN
double dose_interval = 24.0;
double spiro_Css  = (spiro_dose/1000.0 * ka_spiro) / (CL_spiro + 1e-9);
double eple_Css   = (eple_dose/1000.0  * ka_eple)  / (CL_eple  + 1e-9);
double fine_Css   = (fine_dose/1000.0  * ka_fine)  / (CL_fine  + 1e-9);
double acei_Css   = (acei_dose/1000.0  * ka_acei)  / (CL_acei  + 1e-9);
double ccb_Css    = (ccb_dose/1000.0   * ka_ccb)   / (CL_ccb   + 1e-9);

$ODE
// Drug PK (simplified 1-cmt)
dxdt_C_spiro = -CL_spiro/Vd_spiro * C_spiro + ka_spiro * spiro_dose/1000.0;
dxdt_C_canr  = ka_canr * C_spiro - CL_canr/Vd_canr * C_canr;
dxdt_C_eple  = -CL_eple/Vd_eple  * C_eple  + ka_eple  * eple_dose/1000.0;
dxdt_C_fine  = -CL_fine/Vd_fine  * C_fine  + ka_fine  * fine_dose/1000.0;
dxdt_C_acei  = -CL_acei/Vd_acei  * C_acei  + ka_acei  * acei_dose/1000.0;
dxdt_C_ccb   = -CL_ccb/Vd_ccb    * C_ccb   + ka_ccb   * ccb_dose/1000.0;

// MRA PD — composite inhibition of ENaC
double MRA_inh = Emax_MRA * (
  (C_spiro/(C_spiro + IC50_spiro) +
   C_canr /(C_canr  + IC50_canr)  +
   C_eple /(C_eple  + IC50_eple)  +
   C_fine /(C_fine  + IC50_fine)) / 4.0
);
MRA_inh = (MRA_inh > 0.95) ? 0.95 : MRA_inh;

// ACEi suppression of AngII
double ACEi_eff = Emax_acei * C_acei / (C_acei + IC50_acei_angII);

// CCB on TPR
double CCB_eff  = Emax_ccb * C_ccb / (C_ccb + IC50_ccb_tpr);

// Surgery effect
double surg_eff = surgery * surgery_effect;

// RAAS cascade
double APA_prod  = APA_severity * (1.0 - surg_eff) * APA_act;
double BAH_prod  = BAH_severity * 0.5;
double Renin_kin = Renin_base * (1.0 - 0.7 * (Aldo_c/Aldo_base > 2.0 ? 1.0 : Aldo_c/Aldo_base/2.0));
dxdt_Renin_c  = Renin_kin - Renin_kout * Renin_c;
double AngII_kin = 0.6 * Renin_c * (1.0 - ACEi_eff);
dxdt_AngII_c  = AngII_kin - AngII_kout * AngII_c;
double Aldo_kin  = Aldo_base * (0.3 * AngII_c + APA_prod * 3.5 + BAH_prod + 0.2);
dxdt_Aldo_c   = Aldo_kin - Aldo_kout * Aldo_c;

// APA activity (slow structural change; surgery zeroes it)
dxdt_APA_act  = -0.001 * surg_eff * APA_act;
dxdt_CYP11B2_c = 0.05 * (APA_prod - CYP11B2_c);

// ENaC activation by aldosterone (MRA blocks)
double ENaC_kin  = 0.8 * (Aldo_c / 8.0) * (1.0 - MRA_inh);
dxdt_ENaC_act = ENaC_kin - ENaC_kout * ENaC_act;

// Sodium (retention driven by ENaC)
double Na_kin   = 140.0 * 0.05 + 1.5 * (ENaC_act - 1.0);
dxdt_Na_c  = Na_kin - Na_kout * Na_c;

// Potassium (wasting driven by ENaC lumenal electronegativity)
double K_kin    = 4.0 * 0.08 - 0.6 * (ENaC_act - 1.0) * (ENaC_act > 1.0 ? 1.0 : 0.0);
dxdt_K_c   = K_kin - K_kout * K_c;

// Bicarbonate (alkalosis from H+ secretion)
double HCO3_kin = 24.0 * 0.04 + 0.15 * (ENaC_act - 1.0);
dxdt_HCO3_c = HCO3_kin - HCO3_kout * HCO3_c;

// Volume (Na retention → volume expansion)
double Vol_kin  = 1.0 * 0.03 + 0.04 * (Na_c - 140.0) / 140.0;
dxdt_Vol_c  = Vol_kin - Vol_kout * Vol_c;

// MAP and TPR
double MAP_kin  = 95.0 * 0.02 + 0.8 * (Vol_c - 1.0) + 5.0 * (TPR_c - 1.0);
dxdt_MAP_c = MAP_kin - MAP_kout * MAP_c;
double TPR_kin  = 1.0 * 0.04 + 0.3 * (Aldo_c / 8.0 - 1.0) * 0.04;
dxdt_TPR_c = TPR_kin * (1.0 - CCB_eff) - TPR_kout * TPR_c;

// GFR (may decline with high MAP and fibrosis)
double GFR_kin  = 1.0 * 0.01 - 0.001 * (MAP_c - 95.0) / 95.0;
dxdt_GFR_c = GFR_kin - GFR_kout * GFR_c;

// Target organ damage
double extra_fine = (C_fine > 0) ? 0.0001 : 0.0;
dxdt_CardFib = CardFib_kin * (Aldo_c / 8.0) * (MAP_c / 95.0) - (CardFib_kout + extra_fine) * CardFib;
dxdt_LVMi_c  = 0.005 * (MAP_c / 95.0 + CardFib * 2.0) - LVMi_kout * LVMi_c;

// Biomarkers
double ARR_num  = (Aldo_c * 18.0); // PAC ng/dL → pmol/L conversion proxy
double ARR_den  = (Renin_c * 0.5 + 0.01);
dxdt_ARR_c  = 0.1 * (ARR_num / ARR_den - ARR_c);
dxdt_HOMA_proxy = 0.01 * (Na_c / 140.0 * 1.2 - HOMA_proxy); // proxy IR from hypernatremia

$CAPTURE
C_spiro C_canr C_eple C_fine C_acei C_ccb
Renin_c AngII_c Aldo_c ENaC_act
Na_c K_c HCO3_c Vol_c
MAP_c TPR_c GFR_c
CardFib LVMi_c ARR_c HOMA_proxy
'

# Load model once at startup
pa_mod <- mread("pa", tempdir(), pa_model_code)

# Helper: run simulation
run_sim <- function(params, end_time = 365, delta = 1) {
  param_list <- as.list(params)
  out <- pa_mod %>%
    param(param_list) %>%
    init(Renin_c = if (!is.null(params$Renin_init)) params$Renin_init else 0.3,
         Aldo_c  = if (!is.null(params$Aldo_init))  params$Aldo_init  else 28.0,
         K_c     = if (!is.null(params$K_init))     params$K_init     else 3.2,
         MAP_c   = if (!is.null(params$MAP_init))   params$MAP_init   else 110.0,
         LVMi_c  = if (!is.null(params$LVMi_init))  params$LVMi_init  else 110.0,
         ARR_c   = if (!is.null(params$ARR_init))   params$ARR_init   else 100.0) %>%
    mrgsim(end = end_time, delta = delta) %>%
    as_tibble()
  out
}

# ── UI ────────────────────────────────────────────────────────────────────────
ui <- dashboardPage(
  skin = "red",
  dashboardHeader(title = "Primary Aldosteronism QSP"),
  dashboardSidebar(
    sidebarMenu(
      menuItem("Patient Profile",       tabName = "profile",   icon = icon("user")),
      menuItem("RAAS / Drug PK",        tabName = "raas",      icon = icon("pills")),
      menuItem("Aldosterone Panel",     tabName = "aldo",      icon = icon("vials")),
      menuItem("Ion Homeostasis",       tabName = "ions",      icon = icon("tint")),
      menuItem("Cardiovascular & TOD",  tabName = "cv",        icon = icon("heart")),
      menuItem("Scenario Comparison",   tabName = "scenarios", icon = icon("chart-bar")),
      menuItem("Biomarker Explorer",    tabName = "biomarker", icon = icon("microscope"))
    )
  ),
  dashboardBody(
    tags$head(tags$style(HTML("
      .content-wrapper { background: #f9f9f9; }
      .box { border-radius: 6px; }
      .ggplotly-output { width: 100% !important; }
    "))),
    tabItems(

      # ── Tab 1: Patient Profile ───────────────────────────────────────────
      tabItem(tabName = "profile",
        fluidRow(
          box(title = "Disease Subtype & Severity", width = 4, solidHeader = TRUE, status = "danger",
            sliderInput("APA_severity", "APA Severity (0=none, 1=typical)", 0, 1, 0.9, 0.05),
            sliderInput("BAH_severity", "BAH Severity (0=none, 1=bilateral)", 0, 1, 0, 0.05),
            checkboxInput("surgery", "Adrenalectomy performed", FALSE),
            sliderInput("surgery_effect", "Adrenalectomy Efficacy", 0.5, 1.0, 0.90, 0.05),
            hr(),
            h5("Baseline Laboratory Values"),
            sliderInput("Aldo_init", "Baseline PAC (ng/dL → scaled)", 5, 50, 28, 1),
            sliderInput("K_init",    "Baseline K⁺ (mmol/L)",         2.5, 4.5, 3.2, 0.1),
            sliderInput("MAP_init",  "Baseline MAP (mmHg)",           80, 140, 110, 1)
          ),
          box(title = "Simulation Settings", width = 4, solidHeader = TRUE, status = "warning",
            sliderInput("sim_duration", "Simulation Duration (days)", 30, 730, 365, 30),
            hr(),
            h5("Drug Doses (mg/day)"),
            sliderInput("spiro_dose", "Spironolactone",  0, 400, 0, 25),
            sliderInput("eple_dose",  "Eplerenone",      0, 200, 0, 25),
            sliderInput("fine_dose",  "Finerenone",      0, 40,  0, 5),
            sliderInput("acei_dose",  "Ramipril (ACEi)", 0, 20,  0, 1),
            sliderInput("ccb_dose",   "Amlodipine (CCB)",0, 20,  0, 1)
          ),
          box(title = "About This Model", width = 4, solidHeader = TRUE, status = "info",
            tags$p("This QSP model simulates Primary Aldosteronism (Conn's Syndrome) pathophysiology with:"),
            tags$ul(
              tags$li("23 ODE compartments"),
              tags$li("Full RAAS cascade (Renin → AngII → Aldosterone)"),
              tags$li("Adrenal adenoma (APA) & bilateral hyperplasia (BAH)"),
              tags$li("ENaC/ROMK renal ion handling"),
              tags$li("6 drug classes (MRA, ACEi, CCB, adrenalectomy)"),
              tags$li("Target organ damage (cardiac fibrosis, LVH)")
            ),
            hr(),
            tags$p(strong("Disease prevalence:"), " 5-10% of hypertensive patients"),
            tags$p(strong("Key biomarkers:"), " ARR > 30 (ng/dL)/(ng/mL/h), PAC > 15 ng/dL"),
            tags$p(strong("Mutation hotspots:"), " KCNJ5, CACNA1D, ATP1A1, ATP2B3")
          )
        ),
        fluidRow(
          box(title = "Model Diagram Summary", width = 12, solidHeader = TRUE, status = "primary",
            tags$div(style = "text-align:center;",
              tags$p("Adrenal Cortex → CYP11B2 → Aldosterone → MR → ENaC/SGK1 → Na⁺ retention / K⁺ wasting"),
              tags$p("↓"),
              tags$p("Volume expansion → Hypertension → LVH / Cardiac Fibrosis / Renal Damage"),
              tags$p("Renin suppression (ARR elevation) — hallmark of primary vs secondary aldosteronism")
            )
          )
        )
      ),

      # ── Tab 2: RAAS / Drug PK ────────────────────────────────────────────
      tabItem(tabName = "raas",
        fluidRow(
          box(title = "RAAS Cascade Over Time", width = 8, solidHeader = TRUE, status = "danger",
            plotlyOutput("plot_raas", height = "400px")
          ),
          box(title = "Interpretation", width = 4, solidHeader = TRUE, status = "info",
            tags$ul(
              tags$li("Renin is suppressed (< 1 ng/mL/h) in untreated APA due to volume feedback"),
              tags$li("AngII remains low-normal despite high aldosterone (autonomous APA production bypasses AngII stimulus)"),
              tags$li("Aldosterone rebounds post-adrenalectomy normalization"),
              tags$li("ACEi reduces AngII but has limited effect on APA-driven aldosterone")
            ),
            hr(),
            tags$p(strong("Normal ranges:")),
            tags$p("Renin: 0.5–4 ng/mL/h | AngII: 10–60 pg/mL | Aldosterone: 4–31 ng/dL")
          )
        ),
        fluidRow(
          box(title = "Drug Plasma Concentrations", width = 6, solidHeader = TRUE, status = "warning",
            plotlyOutput("plot_pk", height = "350px")
          ),
          box(title = "MR Occupancy & Inhibition", width = 6, solidHeader = TRUE, status = "primary",
            plotlyOutput("plot_mr_occ", height = "350px")
          )
        )
      ),

      # ── Tab 3: Aldosterone Panel ─────────────────────────────────────────
      tabItem(tabName = "aldo",
        fluidRow(
          box(title = "Aldosterone & ARR Trajectory", width = 8, solidHeader = TRUE, status = "danger",
            plotlyOutput("plot_aldo_arr", height = "400px")
          ),
          box(title = "ARR Diagnostic Thresholds", width = 4, solidHeader = TRUE, status = "info",
            tags$ul(
              tags$li("ARR ≥ 30 (ng/dL per ng/mL/h): screening positive"),
              tags$li("PAC ≥ 10 ng/dL required for confirmation"),
              tags$li("Confirmatory tests: salt load, fludrocortisone suppression, captopril challenge"),
              tags$li("AVS lateralization ratio ≥ 4:1 → APA (unilateral)"),
              tags$li("AVS lateralization ratio < 3:1 → BAH (bilateral)")
            ),
            hr(),
            tags$p("MRA therapy normalises ARR by suppressing the renin-suppression loop (renin rebounds).")
          )
        ),
        fluidRow(
          box(title = "CYP11B2 Activity", width = 6, solidHeader = TRUE, status = "warning",
            plotlyOutput("plot_cyp", height = "300px")
          ),
          box(title = "ENaC Activation", width = 6, solidHeader = TRUE, status = "primary",
            plotlyOutput("plot_enac", height = "300px")
          )
        )
      ),

      # ── Tab 4: Ion Homeostasis ────────────────────────────────────────────
      tabItem(tabName = "ions",
        fluidRow(
          box(title = "Electrolyte Dynamics", width = 8, solidHeader = TRUE, status = "danger",
            plotlyOutput("plot_ions", height = "420px")
          ),
          box(title = "Ion Physiology Notes", width = 4, solidHeader = TRUE, status = "info",
            tags$p(strong("Hypokalemia (K⁺ < 3.5 mEq/L):")),
            tags$ul(
              tags$li("ENaC activation creates lumenal electronegativity"),
              tags$li("ROMK-driven K⁺ secretion → urinary K⁺ wasting"),
              tags$li("Hypokalemia in only 9–37% of PA (often normokalemic)")
            ),
            hr(),
            tags$p(strong("Metabolic Alkalosis (HCO₃⁻ > 26 mEq/L):")),
            tags$ul(
              tags$li("H⁺/K⁺-ATPase H⁺ secretion in parallel with Na⁺ retention"),
              tags$li("K⁺ shifts into cells → further alkalosis")
            ),
            hr(),
            tags$p(strong("Volume expansion:")),
            tags$ul(
              tags$li("Na⁺ retention → ~10-15% plasma volume increase"),
              tags$li("Partial escape via ANP/BNP (escape phenomenon)")
            )
          )
        ),
        fluidRow(
          box(title = "Volume & Renal Function", width = 6, solidHeader = TRUE, status = "warning",
            plotlyOutput("plot_vol_gfr", height = "300px")
          ),
          box(title = "Ion Summary Table", width = 6, solidHeader = TRUE, status = "primary",
            DTOutput("ion_table")
          )
        )
      ),

      # ── Tab 5: Cardiovascular & TOD ──────────────────────────────────────
      tabItem(tabName = "cv",
        fluidRow(
          box(title = "Blood Pressure & Cardiac Fibrosis", width = 8, solidHeader = TRUE, status = "danger",
            plotlyOutput("plot_cv_main", height = "400px")
          ),
          box(title = "Cardiovascular Risk Notes", width = 4, solidHeader = TRUE, status = "info",
            tags$p("PA patients have ~4× higher cardiovascular event rate than essential HTN of equal BP."),
            tags$ul(
              tags$li("Aldosterone directly promotes cardiac fibrosis via MR on fibroblasts"),
              tags$li("LV mass index (LVMi) elevated beyond that explained by BP"),
              tags$li("Finerenone (non-steroidal MRA) shows superior anti-fibrotic effect"),
              tags$li("Endothelial dysfunction, oxidative stress mediated by MR")
            ),
            hr(),
            tags$p(strong("LVMi targets:")),
            tags$p("Normal: <95 g/m² (F), <115 g/m² (M)"),
            tags$p("Adrenalectomy reduces LVMi by ~20% over 6-12 months")
          )
        ),
        fluidRow(
          box(title = "LV Mass Index Trajectory", width = 6, solidHeader = TRUE, status = "warning",
            plotlyOutput("plot_lvmi", height = "300px")
          ),
          box(title = "GFR Over Time", width = 6, solidHeader = TRUE, status = "primary",
            plotlyOutput("plot_gfr", height = "300px")
          )
        )
      ),

      # ── Tab 6: Scenario Comparison ───────────────────────────────────────
      tabItem(tabName = "scenarios",
        fluidRow(
          box(title = "Predefined Scenarios", width = 12, solidHeader = TRUE, status = "danger",
            fluidRow(
              column(3, checkboxGroupInput("scen_select", "Select Scenarios:",
                choices = c(
                  "Untreated APA"         = "s1",
                  "Adrenalectomy"         = "s2",
                  "Spironolactone 100mg"  = "s3",
                  "Eplerenone 100mg"      = "s4",
                  "Finerenone 20mg"       = "s5",
                  "Spiro + Amlodipine"    = "s6",
                  "Normal Control"        = "s7"
                ),
                selected = c("s1","s2","s3","s5","s7")
              )),
              column(3,
                selectInput("scen_endpoint", "Endpoint:",
                  choices = c("MAP (mmHg)" = "MAP_c",
                              "K⁺ (mEq/L)" = "K_c",
                              "Aldo (ng/dL)" = "Aldo_c",
                              "ARR" = "ARR_c",
                              "CardFib" = "CardFib",
                              "LVMi" = "LVMi_c",
                              "GFR" = "GFR_c")),
                sliderInput("scen_duration", "Duration (days)", 90, 730, 365, 30)
              ),
              column(6,
                plotlyOutput("plot_scenarios", height = "380px")
              )
            )
          )
        ),
        fluidRow(
          box(title = "Scenario Summary at End-Point", width = 12, solidHeader = TRUE, status = "warning",
            DTOutput("scenario_table")
          )
        )
      ),

      # ── Tab 7: Biomarker Explorer ────────────────────────────────────────
      tabItem(tabName = "biomarker",
        fluidRow(
          box(title = "Biomarker Panel", width = 8, solidHeader = TRUE, status = "danger",
            plotlyOutput("plot_biomarkers", height = "450px")
          ),
          box(title = "Biomarker Interpretation", width = 4, solidHeader = TRUE, status = "info",
            tags$h5("ARR (Aldosterone-to-Renin Ratio)"),
            tags$ul(
              tags$li("Best screening test for PA"),
              tags$li("Threshold: > 30 (ng/dL)/(ng/mL/h)"),
              tags$li("Renin suppression is key driver (not just aldosterone elevation)")
            ),
            hr(),
            tags$h5("PAC (Plasma Aldosterone Concentration)"),
            tags$ul(
              tags$li("Elevated in APA: typically 15–100 ng/dL"),
              tags$li("Normal: 4–31 ng/dL"),
              tags$li("Post-saline infusion > 10 ng/dL confirms PA")
            ),
            hr(),
            tags$h5("Insulin Resistance Proxy"),
            tags$ul(
              tags$li("Hypernatremia and volume expansion worsen IR"),
              tags$li("Aldosterone activates IRS-1 serine phosphorylation"),
              tags$li("Treatment improves metabolic markers")
            )
          )
        ),
        fluidRow(
          box(title = "Dose-Response: ARR vs MRA Dose", width = 6, solidHeader = TRUE, status = "warning",
            plotlyOutput("plot_dose_arr", height = "300px")
          ),
          box(title = "K⁺ Recovery with MRA Treatment", width = 6, solidHeader = TRUE, status = "primary",
            plotlyOutput("plot_k_dose", height = "300px")
          )
        )
      )
    )
  )
)

# ── Server ────────────────────────────────────────────────────────────────────
server <- function(input, output, session) {

  # Reactive: collect current patient parameters
  cur_params <- reactive({
    list(
      APA_severity = input$APA_severity,
      BAH_severity = input$BAH_severity,
      surgery      = as.numeric(input$surgery),
      surgery_effect = input$surgery_effect,
      spiro_dose   = input$spiro_dose,
      eple_dose    = input$eple_dose,
      fine_dose    = input$fine_dose,
      acei_dose    = input$acei_dose,
      ccb_dose     = input$ccb_dose,
      Aldo_init    = input$Aldo_init,
      K_init       = input$K_init,
      MAP_init     = input$MAP_init,
      Renin_init   = 0.3,
      LVMi_init    = 110.0,
      ARR_init     = 100.0
    )
  })

  # Run simulation
  sim_data <- reactive({
    run_sim(cur_params(), end_time = input$sim_duration)
  })

  # ── Tab 2: RAAS plots ─────────────────────────────────────────────────────
  output$plot_raas <- renderPlotly({
    d <- sim_data()
    p <- ggplot(d, aes(x = time)) +
      geom_line(aes(y = Renin_c, color = "Renin (scaled)"), size = 0.9) +
      geom_line(aes(y = AngII_c, color = "AngII (scaled)"), size = 0.9) +
      geom_line(aes(y = Aldo_c / 10, color = "Aldosterone /10"), size = 0.9) +
      scale_color_manual(values = c("Renin (scaled)" = "#e74c3c",
                                    "AngII (scaled)" = "#f39c12",
                                    "Aldosterone /10" = "#8e44ad")) +
      labs(x = "Day", y = "Relative Units", color = "", title = "RAAS Cascade") +
      theme_bw() + theme(legend.position = "bottom")
    ggplotly(p)
  })

  output$plot_pk <- renderPlotly({
    d <- sim_data()
    d_long <- d %>% select(time, C_spiro, C_canr, C_eple, C_fine, C_acei, C_ccb) %>%
      pivot_longer(-time, names_to = "drug", values_to = "conc")
    p <- ggplot(d_long, aes(x = time, y = conc, color = drug)) +
      geom_line(size = 0.8) +
      labs(x = "Day", y = "Plasma Conc (μg/L)", title = "Drug PK", color = "Drug") +
      theme_bw() + theme(legend.position = "bottom")
    ggplotly(p)
  })

  output$plot_mr_occ <- renderPlotly({
    d <- sim_data()
    ic50s <- c(C_spiro = 1.2, C_canr = 0.8, C_eple = 2.5, C_fine = 0.65)
    d$MR_occ <- with(d, 0.88 * ((C_spiro/(C_spiro+1.2) + C_canr/(C_canr+0.8) +
                                  C_eple/(C_eple+2.5)  + C_fine/(C_fine+0.65)) / 4))
    p <- ggplot(d, aes(x = time, y = MR_occ * 100)) +
      geom_line(color = "#2980b9", size = 1) +
      geom_hline(yintercept = 50, linetype = "dashed", color = "gray50") +
      labs(x = "Day", y = "MR Occupancy (%)", title = "Mineralocorticoid Receptor Blockade") +
      ylim(0, 100) + theme_bw()
    ggplotly(p)
  })

  # ── Tab 3: Aldosterone ────────────────────────────────────────────────────
  output$plot_aldo_arr <- renderPlotly({
    d <- sim_data()
    p <- ggplot(d, aes(x = time)) +
      geom_line(aes(y = Aldo_c, color = "PAC (ng/dL)"), size = 1) +
      geom_line(aes(y = ARR_c / 5, color = "ARR /5"), size = 1) +
      geom_hline(yintercept = 10, linetype = "dashed", color = "#e74c3c") +
      annotate("text", x = max(d$time)*0.8, y = 10.5, label = "PAC threshold (10 ng/dL)", color = "#e74c3c", size = 3) +
      scale_color_manual(values = c("PAC (ng/dL)" = "#8e44ad", "ARR /5" = "#e67e22")) +
      labs(x = "Day", y = "Value (see scale)", title = "Aldosterone Panel", color = "") +
      theme_bw() + theme(legend.position = "bottom")
    ggplotly(p)
  })

  output$plot_cyp <- renderPlotly({
    d <- sim_data()
    p <- ggplot(d, aes(x = time, y = CYP11B2_c)) +
      geom_line(color = "#c0392b", size = 1) +
      labs(x = "Day", y = "CYP11B2 Activity (relative)", title = "Aldosterone Synthase (CYP11B2)") +
      theme_bw()
    ggplotly(p)
  })

  output$plot_enac <- renderPlotly({
    d <- sim_data()
    p <- ggplot(d, aes(x = time, y = ENaC_act)) +
      geom_line(color = "#1abc9c", size = 1) +
      geom_hline(yintercept = 1.0, linetype = "dashed", color = "gray60") +
      labs(x = "Day", y = "ENaC Activation (relative)", title = "ENaC / SGK1 Activity") +
      theme_bw()
    ggplotly(p)
  })

  # ── Tab 4: Ions ───────────────────────────────────────────────────────────
  output$plot_ions <- renderPlotly({
    d <- sim_data()
    p <- ggplot(d, aes(x = time)) +
      geom_line(aes(y = Na_c,   color = "Na⁺ (mEq/L)"),  size = 0.9) +
      geom_line(aes(y = K_c * 20, color = "K⁺ ×20"),     size = 0.9) +
      geom_line(aes(y = HCO3_c, color = "HCO₃⁻ (mEq/L)"), size = 0.9) +
      scale_color_manual(values = c("Na⁺ (mEq/L)" = "#2980b9",
                                    "K⁺ ×20" = "#e74c3c",
                                    "HCO₃⁻ (mEq/L)" = "#27ae60")) +
      labs(x = "Day", y = "Concentration", title = "Electrolyte Panel", color = "") +
      theme_bw() + theme(legend.position = "bottom")
    ggplotly(p)
  })

  output$plot_vol_gfr <- renderPlotly({
    d <- sim_data()
    p <- ggplot(d, aes(x = time)) +
      geom_line(aes(y = Vol_c, color = "Volume (relative)"), size = 1) +
      geom_line(aes(y = GFR_c, color = "GFR (relative)"), size = 1) +
      scale_color_manual(values = c("Volume (relative)" = "#3498db", "GFR (relative)" = "#e67e22")) +
      labs(x = "Day", y = "Relative Units", title = "Volume & GFR", color = "") +
      theme_bw() + theme(legend.position = "bottom")
    ggplotly(p)
  })

  output$ion_table <- renderDT({
    d <- tail(sim_data(), 1)
    df <- data.frame(
      Marker   = c("Na⁺", "K⁺", "HCO₃⁻", "Volume", "GFR"),
      Value    = round(c(d$Na_c, d$K_c, d$HCO3_c, d$Vol_c, d$GFR_c), 2),
      Unit     = c("mEq/L","mEq/L","mEq/L","relative","relative"),
      Normal   = c("136–145","3.5–5.0","22–26","1.0","1.0")
    )
    datatable(df, options = list(pageLength = 5, dom = "t"), rownames = FALSE)
  })

  # ── Tab 5: CV ─────────────────────────────────────────────────────────────
  output$plot_cv_main <- renderPlotly({
    d <- sim_data()
    p <- ggplot(d, aes(x = time)) +
      geom_line(aes(y = MAP_c, color = "MAP (mmHg)"), size = 1) +
      geom_line(aes(y = CardFib * 500 + 80, color = "Cardiac Fibrosis ×500+80"), size = 1) +
      geom_hline(yintercept = 100, linetype = "dashed", color = "gray50") +
      scale_color_manual(values = c("MAP (mmHg)" = "#e74c3c", "Cardiac Fibrosis ×500+80" = "#8e44ad")) +
      labs(x = "Day", y = "MAP (mmHg) / Scaled CardFib", title = "Cardiovascular Endpoints", color = "") +
      theme_bw() + theme(legend.position = "bottom")
    ggplotly(p)
  })

  output$plot_lvmi <- renderPlotly({
    d <- sim_data()
    p <- ggplot(d, aes(x = time, y = LVMi_c)) +
      geom_line(color = "#c0392b", size = 1) +
      geom_hline(yintercept = 115, linetype = "dashed", color = "#e74c3c") +
      annotate("text", x = max(d$time)*0.7, y = 117, label = "LVH threshold ♂ (115 g/m²)", color = "#e74c3c", size = 3) +
      labs(x = "Day", y = "LVMi (g/m²)", title = "Left Ventricular Mass Index") +
      theme_bw()
    ggplotly(p)
  })

  output$plot_gfr <- renderPlotly({
    d <- sim_data()
    p <- ggplot(d, aes(x = time, y = GFR_c)) +
      geom_line(color = "#27ae60", size = 1) +
      geom_hline(yintercept = 0.6, linetype = "dashed", color = "orange") +
      labs(x = "Day", y = "GFR (relative to normal)", title = "Glomerular Filtration Rate") +
      theme_bw()
    ggplotly(p)
  })

  # ── Tab 6: Scenarios ──────────────────────────────────────────────────────
  scenario_params <- list(
    s1 = list(APA_severity=1.0, BAH_severity=0, surgery=0, spiro_dose=0, eple_dose=0,
              fine_dose=0, acei_dose=0, ccb_dose=0, Aldo_init=28, K_init=3.2,
              MAP_init=110, Renin_init=0.3, LVMi_init=110, ARR_init=100),
    s2 = list(APA_severity=1.0, BAH_severity=0, surgery=1, surgery_effect=0.90, spiro_dose=0,
              eple_dose=0, fine_dose=0, acei_dose=0, ccb_dose=0, Aldo_init=28, K_init=3.2,
              MAP_init=110, Renin_init=0.3, LVMi_init=110, ARR_init=100),
    s3 = list(APA_severity=0.5, BAH_severity=0.5, surgery=0, spiro_dose=100, eple_dose=0,
              fine_dose=0, acei_dose=0, ccb_dose=0, Aldo_init=20, K_init=3.4,
              MAP_init=108, Renin_init=0.4, LVMi_init=108, ARR_init=80),
    s4 = list(APA_severity=0.5, BAH_severity=0.5, surgery=0, spiro_dose=0, eple_dose=100,
              fine_dose=0, acei_dose=0, ccb_dose=0, Aldo_init=20, K_init=3.4,
              MAP_init=108, Renin_init=0.4, LVMi_init=108, ARR_init=80),
    s5 = list(APA_severity=0.5, BAH_severity=0.5, surgery=0, spiro_dose=0, eple_dose=0,
              fine_dose=20, acei_dose=0, ccb_dose=0, Aldo_init=20, K_init=3.4,
              MAP_init=108, Renin_init=0.4, LVMi_init=108, ARR_init=80),
    s6 = list(APA_severity=1.0, BAH_severity=0, surgery=0, spiro_dose=100, eple_dose=0,
              fine_dose=0, acei_dose=0, ccb_dose=10, Aldo_init=28, K_init=3.2,
              MAP_init=110, Renin_init=0.3, LVMi_init=110, ARR_init=100),
    s7 = list(APA_severity=0, BAH_severity=0, surgery=0, spiro_dose=0, eple_dose=0,
              fine_dose=0, acei_dose=0, ccb_dose=0, Aldo_init=8, K_init=4.0,
              MAP_init=90, Renin_init=1.0, LVMi_init=90, ARR_init=15)
  )

  scenario_labels <- c(s1="Untreated APA", s2="Adrenalectomy", s3="Spiro 100mg",
                       s4="Eplerenone 100mg", s5="Finerenone 20mg",
                       s6="Spiro+Amlodipine", s7="Normal Control")
  scenario_colors <- c(s1="#e74c3c", s2="#27ae60", s3="#3498db",
                       s4="#f39c12", s5="#9b59b6", s6="#1abc9c", s7="#95a5a6")

  output$plot_scenarios <- renderPlotly({
    req(input$scen_select)
    endpoint <- input$scen_endpoint
    sims <- lapply(input$scen_select, function(s) {
      d <- run_sim(scenario_params[[s]], end_time = input$scen_duration)
      d$scenario <- scenario_labels[s]
      d$color    <- scenario_colors[s]
      d
    })
    combined <- bind_rows(sims)
    p <- ggplot(combined, aes_string(x = "time", y = endpoint, color = "scenario")) +
      geom_line(size = 1) +
      scale_color_manual(values = setNames(unique(combined$color), unique(combined$scenario))) +
      labs(x = "Day", y = endpoint, color = "Scenario", title = paste("Endpoint:", endpoint)) +
      theme_bw() + theme(legend.position = "bottom")
    ggplotly(p)
  })

  output$scenario_table <- renderDT({
    all_scens <- names(scenario_params)
    endpoint  <- input$scen_endpoint
    if (is.null(endpoint)) endpoint <- "MAP_c"
    rows <- lapply(all_scens, function(s) {
      d <- run_sim(scenario_params[[s]], end_time = 365)
      last <- tail(d, 1)
      data.frame(Scenario = scenario_labels[s],
                 MAP    = round(last$MAP_c, 1),
                 K      = round(last$K_c, 2),
                 Aldo   = round(last$Aldo_c, 1),
                 ARR    = round(last$ARR_c, 0),
                 CardFib= round(last$CardFib, 4),
                 LVMi   = round(last$LVMi_c, 1),
                 GFR    = round(last$GFR_c, 3),
                 stringsAsFactors = FALSE)
    })
    df <- bind_rows(rows)
    datatable(df, options = list(pageLength = 8, dom = "t"), rownames = FALSE)
  })

  # ── Tab 7: Biomarker ──────────────────────────────────────────────────────
  output$plot_biomarkers <- renderPlotly({
    d <- sim_data()
    p <- ggplot(d, aes(x = time)) +
      geom_line(aes(y = ARR_c,      color = "ARR"), size = 1) +
      geom_line(aes(y = Aldo_c * 3, color = "PAC ×3"), size = 1) +
      geom_line(aes(y = HOMA_proxy * 30, color = "HOMA proxy ×30"), size = 1) +
      geom_hline(yintercept = 30, linetype = "dashed", color = "#e74c3c") +
      annotate("text", x = max(d$time)*0.8, y = 32, label = "ARR threshold = 30", color = "#e74c3c", size = 3) +
      scale_color_manual(values = c("ARR" = "#e74c3c", "PAC ×3" = "#8e44ad",
                                    "HOMA proxy ×30" = "#f39c12")) +
      labs(x = "Day", y = "Value (scaled)", title = "Biomarker Panel", color = "") +
      theme_bw() + theme(legend.position = "bottom")
    ggplotly(p)
  })

  output$plot_dose_arr <- renderPlotly({
    doses <- seq(0, 200, by = 25)
    end_arr <- sapply(doses, function(d) {
      p <- list(APA_severity=0.8, BAH_severity=0.2, surgery=0, spiro_dose=d,
                eple_dose=0, fine_dose=0, acei_dose=0, ccb_dose=0,
                Aldo_init=25, K_init=3.3, MAP_init=108,
                Renin_init=0.35, LVMi_init=108, ARR_init=90)
      s <- run_sim(p, end_time = 90)
      tail(s$ARR_c, 1)
    })
    df <- data.frame(Dose = doses, ARR = end_arr)
    p <- ggplot(df, aes(x = Dose, y = ARR)) +
      geom_line(color = "#e74c3c", size = 1) +
      geom_point(color = "#c0392b", size = 2) +
      geom_hline(yintercept = 30, linetype = "dashed", color = "gray50") +
      labs(x = "Spironolactone Dose (mg/day)", y = "ARR at Day 90",
           title = "ARR vs Spironolactone Dose") +
      theme_bw()
    ggplotly(p)
  })

  output$plot_k_dose <- renderPlotly({
    doses <- seq(0, 200, by = 25)
    end_k <- sapply(doses, function(d) {
      p <- list(APA_severity=0.8, BAH_severity=0.2, surgery=0, spiro_dose=d,
                eple_dose=0, fine_dose=0, acei_dose=0, ccb_dose=0,
                Aldo_init=25, K_init=3.3, MAP_init=108,
                Renin_init=0.35, LVMi_init=108, ARR_init=90)
      s <- run_sim(p, end_time = 90)
      tail(s$K_c, 1)
    })
    df <- data.frame(Dose = doses, K = end_k)
    p <- ggplot(df, aes(x = Dose, y = K)) +
      geom_line(color = "#27ae60", size = 1) +
      geom_point(color = "#1e8449", size = 2) +
      geom_hline(yintercept = 3.5, linetype = "dashed", color = "orange") +
      annotate("text", x = 150, y = 3.52, label = "Hypokalemia threshold", color = "orange", size = 3) +
      labs(x = "Spironolactone Dose (mg/day)", y = "K⁺ at Day 90 (mEq/L)",
           title = "K⁺ Recovery vs Spironolactone Dose") +
      theme_bw()
    ggplotly(p)
  })
}

shinyApp(ui, server)
